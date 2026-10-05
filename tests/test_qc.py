import json
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

from variantrail.errors import NotFoundError
from variantrail.server import Handler
from variantrail.service import VariantRail

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
PREAMBLE = "##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1>"

VCF = "\n".join(
    [
        PREAMBLE,
        HEADER,
        "chr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42",
        "chr1\t11856378\trs2\tG\tT\t20\tPASS\tDP=42",
        "chr1\t11856379\t.\tA\tG\t.\tq10\tDP=5",
        "chr7\t55019017\t.\tG\tGA\t50\tPASS\tDP=31",
        "chr7\t55019017\t.\tG\tT\t50\tPASS\t.",
        "chr12\t25245350\t.\tC\tA,T\t80\tPASS\tDP=60",
        "chr17\t43093445\t.\tC\tT\t99\tPASS\tDP=55",
    ]
) + "\n"


def single(chrom="chr1", pos=11856378, ref="G", alt="A", qual="60", filter_value="PASS", info="DP=42"):
    return "\n".join([PREAMBLE, HEADER, f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter_value}\t{info}"]) + "\n"


CLASS_VCF = "\n".join(
    [
        PREAMBLE,
        HEADER,
        "chr1\t100\t.\tN\tA\t50\tPASS\tDP=10",
        "chr1\t200\t.\tAC\tGT\t50\tPASS\tDP=10",
        "chr1\t300\t.\tAC\tA\t50\tPASS\tDP=10",
        "chr1\t400\t.\tG\tGA\t50\tPASS\tDP=10",
        "chr1\t500\t.\tA\tN\t50\tPASS\tDP=10",
    ]
) + "\n"


class RunQcServiceTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.service = VariantRail(str(Path(self.directory.name) / "test.db"))

    def tearDown(self):
        self.directory.cleanup()

    def start_run(self, params=None, run_id="run-1"):
        self.service.create_sample({"id": "s1", "vcf": VCF}, "s1")
        body = {"id": run_id}
        if params is not None:
            body["params"] = params
        return self.service.create_run("s1", body, "r1")

    def test_qc_summarizes_only_retained_records_and_alleles(self):
        run = self.start_run()
        payload = self.service.run_qc("run-1")
        self.assertEqual(
            {"run_id", "sample_id", "sample_sha256", "provenance_head", "records", "alleles"},
            set(payload),
        )
        self.assertEqual("run-1", payload["run_id"])
        self.assertEqual("s1", payload["sample_id"])
        self.assertEqual(run["sample_sha256"], payload["sample_sha256"])
        self.assertEqual(run["provenance_head"], payload["provenance_head"])
        self.assertEqual(
            {"total", "pass", "missing_qual", "missing_dp", "multiallelic", "qual", "dp"},
            set(payload["records"]),
        )
        records = payload["records"]
        self.assertEqual(7, records["total"])
        self.assertEqual(6, records["pass"])  # every record except the q10 one
        self.assertEqual(1, records["missing_qual"])
        self.assertEqual(1, records["missing_dp"])  # the chr7 stop-gained record has no DP
        self.assertEqual(1, records["multiallelic"])  # only C>A,T carries two annotations
        self.assertEqual({"count": 6, "min": 20.0, "max": 99.0}, records["qual"])
        self.assertEqual({"count": 6, "min": 5, "max": 60}, records["dp"])

    def test_qc_classifies_alleles_and_computes_ti_tv(self):
        self.start_run()
        alleles = self.service.run_qc("run-1")["alleles"]
        self.assertEqual(
            {
                "total", "snv", "mnv", "insertion", "deletion", "non_acgt_snv",
                "transitions", "transversions", "ti_tv_ratio",
            },
            set(alleles),
        )
        # G>A, G>T, A>G, G>T, C>A, C>T, C>T are SNVs; G>GA is an insertion.
        self.assertEqual(8, alleles["total"])
        self.assertEqual(7, alleles["snv"])
        self.assertEqual(0, alleles["mnv"])
        self.assertEqual(1, alleles["insertion"])
        self.assertEqual(0, alleles["deletion"])
        self.assertEqual(0, alleles["non_acgt_snv"])
        # transitions: G>A, A>G, C>T, C>T ; transversions: G>T, G>T, C>A
        self.assertEqual(4, alleles["transitions"])
        self.assertEqual(3, alleles["transversions"])
        self.assertEqual(4 / 3, alleles["ti_tv_ratio"])

    def test_qc_reflects_record_and_allele_filters(self):
        self.start_run({"pass_only": True, "min_qual": 30, "min_dp": 10})
        payload = self.service.run_qc("run-1")
        records = payload["records"]
        self.assertEqual(4, records["total"])
        self.assertEqual(4, records["pass"])
        self.assertEqual(0, records["missing_qual"])
        self.assertEqual(0, records["missing_dp"])
        self.assertEqual(1, records["multiallelic"])
        self.assertEqual({"count": 4, "min": 50.0, "max": 99.0}, records["qual"])
        self.assertEqual({"count": 4, "min": 31, "max": 60}, records["dp"])
        alleles = payload["alleles"]
        self.assertEqual(5, alleles["total"])
        self.assertEqual(4, alleles["snv"])
        self.assertEqual(1, alleles["insertion"])
        self.assertEqual(3, alleles["transitions"])
        self.assertEqual(1, alleles["transversions"])
        self.assertEqual(3.0, alleles["ti_tv_ratio"])

    def test_qc_counts_annotations_after_allele_level_filtering(self):
        # The C>T,G record keeps only the tabled HIGH allele; the surviving
        # variant therefore has one annotation and is not multiallelic.
        vcf = "\n".join([PREAMBLE, HEADER, "chr17\t43093445\t.\tC\tT,G\t80\tPASS\tDP=60"]) + "\n"
        self.service.create_sample({"id": "s1", "vcf": vcf}, "s1")
        self.service.create_run("s1", {"id": "run-1", "params": {"impacts": ["HIGH"]}}, "r1")
        payload = self.service.run_qc("run-1")
        self.assertEqual(1, payload["records"]["total"])
        self.assertEqual(0, payload["records"]["multiallelic"])
        self.assertEqual(1, payload["alleles"]["total"])
        self.assertEqual(1, payload["alleles"]["snv"])
        self.assertEqual(1, payload["alleles"]["transitions"])
        self.assertEqual(0, payload["alleles"]["transversions"])

    def test_qc_classifies_mnv_indels_and_n_snvs_with_ratio_null(self):
        self.service.create_sample({"id": "s1", "vcf": CLASS_VCF}, "s1")
        self.service.create_run("s1", {"id": "run-1"}, "r1")
        alleles = self.service.run_qc("run-1")["alleles"]
        self.assertEqual(5, alleles["total"])
        self.assertEqual(0, alleles["snv"])
        self.assertEqual(1, alleles["mnv"])
        self.assertEqual(1, alleles["insertion"])
        self.assertEqual(1, alleles["deletion"])
        self.assertEqual(2, alleles["non_acgt_snv"])
        self.assertEqual(0, alleles["transitions"])
        self.assertEqual(0, alleles["transversions"])
        self.assertIsNone(alleles["ti_tv_ratio"])

    def test_qc_ratio_is_null_when_only_transitions_exist(self):
        self.service.create_sample({"id": "s1", "vcf": single(ref="G", alt="A")}, "s1")
        self.service.create_run("s1", {"id": "run-1"}, "r1")
        alleles = self.service.run_qc("run-1")["alleles"]
        self.assertEqual(1, alleles["snv"])
        self.assertEqual(1, alleles["transitions"])
        self.assertEqual(0, alleles["transversions"])
        self.assertIsNone(alleles["ti_tv_ratio"])

    def test_qc_empty_run_zero_fills_every_count(self):
        run = self.start_run({"genes": ["NO_SUCH_GENE"]})
        payload = self.service.run_qc("run-1")
        self.assertEqual("run-1", payload["run_id"])
        self.assertEqual("s1", payload["sample_id"])
        self.assertEqual(run["sample_sha256"], payload["sample_sha256"])
        self.assertEqual(run["provenance_head"], payload["provenance_head"])
        records = payload["records"]
        self.assertEqual(0, records["total"])
        self.assertEqual(0, records["pass"])
        self.assertEqual(0, records["missing_qual"])
        self.assertEqual(0, records["missing_dp"])
        self.assertEqual(0, records["multiallelic"])
        self.assertEqual({"count": 0, "min": None, "max": None}, records["qual"])
        self.assertEqual({"count": 0, "min": None, "max": None}, records["dp"])
        alleles = payload["alleles"]
        self.assertEqual(0, alleles["total"])
        self.assertEqual(0, alleles["snv"])
        self.assertEqual(0, alleles["mnv"])
        self.assertEqual(0, alleles["insertion"])
        self.assertEqual(0, alleles["deletion"])
        self.assertEqual(0, alleles["non_acgt_snv"])
        self.assertEqual(0, alleles["transitions"])
        self.assertEqual(0, alleles["transversions"])
        self.assertIsNone(alleles["ti_tv_ratio"])

    def test_qc_unknown_run_is_not_found(self):
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.run_qc("missing")

    def test_qc_is_deterministic_and_read_only(self):
        self.start_run()
        before = self.service._document("run-1")
        first = self.service.run_qc("run-1")
        second = self.service.run_qc("run-1")
        self.assertEqual(first, second)
        self.assertEqual(
            json.dumps(first, ensure_ascii=False, separators=(",", ":")),
            json.dumps(second, ensure_ascii=False, separators=(",", ":")),
        )
        for format in ("jsonl", "tsv", "vcf"):
            self.assertEqual(
                self.service.run_export("run-1", format),
                self.service.run_export("run-1", format),
            )
        # The stored run, variants and exports are untouched by QC requests.
        self.assertEqual(before, self.service._document("run-1"))
        self.assertTrue(self.service.run_provenance("run-1")["verified"])


class QcHttpTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        service = VariantRail(str(Path(self.directory.name) / "test.db"))
        service.create_sample({"id": "s1", "vcf": VCF}, "s1")
        service.create_run("s1", {"id": "run-1"}, "r1")
        service.create_run("s1", {"id": "run-2", "params": {"genes": ["NO_SUCH_GENE"]}}, "r2")
        Handler.service = service
        self.service = service
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.port = self.server.server_address[1]
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.directory.cleanup()

    def get(self, path):
        request = urllib.request.Request(f"http://127.0.0.1:{self.port}{path}")
        try:
            with urllib.request.urlopen(request) as response:
                return response.status, response.headers["Content-Type"], response.read()
        except urllib.error.HTTPError as error:
            body = error.read()
            content_type = error.headers["Content-Type"]
            status = error.code
            error.close()
            return status, content_type, body

    def test_qc_endpoint_returns_the_summary(self):
        status, content_type, body = self.get("/runs/run-1/qc")
        self.assertEqual(200, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        payload = json.loads(body)
        self.assertEqual("run-1", payload["run_id"])
        self.assertEqual("s1", payload["sample_id"])
        self.assertEqual(64, len(payload["sample_sha256"]))
        self.assertEqual(64, len(payload["provenance_head"]))
        records = payload["records"]
        self.assertEqual(7, records["total"])
        self.assertEqual(6, records["pass"])
        self.assertEqual(1, records["missing_qual"])
        self.assertEqual(1, records["missing_dp"])
        self.assertEqual(1, records["multiallelic"])
        self.assertEqual({"count": 6, "min": 20.0, "max": 99.0}, records["qual"])
        self.assertEqual({"count": 6, "min": 5, "max": 60}, records["dp"])
        alleles = payload["alleles"]
        self.assertEqual(
            {"total": 8, "snv": 7, "mnv": 0, "insertion": 1, "deletion": 0,
             "non_acgt_snv": 0, "transitions": 4, "transversions": 3},
            {key: alleles[key] for key in alleles if key != "ti_tv_ratio"},
        )
        self.assertEqual(4 / 3, alleles["ti_tv_ratio"])
        # Repeated requests are byte-for-byte identical.
        _, _, again = self.get("/runs/run-1/qc")
        self.assertEqual(body, again)

    def test_qc_endpoint_serves_empty_run(self):
        status, _, body = self.get("/runs/run-2/qc")
        self.assertEqual(200, status)
        payload = json.loads(body)
        self.assertEqual("run-2", payload["run_id"])
        self.assertEqual(0, payload["records"]["total"])
        self.assertEqual({"count": 0, "min": None, "max": None}, payload["records"]["qual"])
        self.assertEqual({"count": 0, "min": None, "max": None}, payload["records"]["dp"])
        self.assertEqual(0, payload["alleles"]["total"])
        self.assertIsNone(payload["alleles"]["ti_tv_ratio"])

    def test_qc_endpoint_errors_use_the_standard_error_object(self):
        status, content_type, body = self.get("/runs/missing/qc")
        self.assertEqual(404, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        error = json.loads(body)["error"]
        self.assertEqual("not_found", error["code"])
        self.assertEqual("run missing was not found", error["message"])

    def test_qc_endpoint_coexists_with_other_run_routes(self):
        status, _, body = self.get("/runs/run-1/qc")
        self.assertEqual(200, status)
        status, _, variants_body = self.get("/runs/run-1/variants")
        self.assertEqual(200, status)
        self.assertEqual(7, json.loads(variants_body)["count"])
        status, _, run_body = self.get("/runs/run-1")
        self.assertEqual(200, status)
        self.assertEqual("succeeded", json.loads(run_body)["status"])
        # Unknown sub-paths keep the existing route-not-found behaviour.
        status, _, body = self.get("/runs/run-1/qc/extra")
        self.assertEqual(404, status)
        self.assertEqual("not_found", json.loads(body)["error"]["code"])

    def test_qc_request_does_not_modify_the_run(self):
        before = self.service._document("run-1")
        self.get("/runs/run-1/qc")
        self.get("/runs/run-1/qc")
        self.assertEqual(before, self.service._document("run-1"))
        self.assertTrue(self.service.run_provenance("run-1")["verified"])


if __name__ == "__main__":
    unittest.main()

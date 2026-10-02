import json
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

from variantrail.errors import ConflictError, NotFoundError, ValidationError
from variantrail.export import TSV_HEADER
from variantrail.server import Handler
from variantrail.service import VariantRail

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
PREAMBLE = "##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1>"

# Lines 5-11 of VCF hold the data records used by most tests.
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


class VariantRailTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.service = VariantRail(str(Path(self.directory.name) / "test.db"))

    def tearDown(self):
        self.directory.cleanup()

    def sample(self, key="s1", vcf=VCF, sample_id="s1"):
        return self.service.create_sample({"id": sample_id, "vcf": vcf}, key)

    def start_run(self, params=None, key="r1", run_id="run-1", sample_id="s1"):
        body = {"id": run_id}
        if params is not None:
            body["params"] = params
        return self.service.create_run(sample_id, body, key)

    # --- ingestion ---------------------------------------------------------

    def test_sample_parses_headers_and_reports_shape(self):
        created = self.sample()
        self.assertEqual(
            {"chromosomes": ["chr1", "chr7", "chr12", "chr17"], "id": "s1", "meta_lines": 3, "records": 7},
            {key: value for key, value in created.items() if key != "sha256"},
        )
        self.assertEqual(64, len(created["sha256"]))

    def test_sample_carries_its_own_record_and_meta_lines(self):
        self.sample()
        document = self.service.get_sample("s1")
        self.assertEqual(3, len(document["meta"]))
        self.assertEqual("##fileformat=VCFv4.2", document["meta"][0]["raw"])
        self.assertEqual("VCFv4.2", document["meta"][0]["value"])
        self.assertEqual("GRCh38", document["meta"][1]["value"])
        self.assertEqual("chr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42", document["records"][0]["raw"])
        self.assertEqual(5, document["records"][0]["line"])

    def test_lowercase_sequence_fields_are_normalised(self):
        self.sample(vcf=single(ref="g", alt="a"))
        document = self.service.get_sample("s1")
        self.assertEqual("G", document["records"][0]["ref"])
        self.assertEqual(["A"], document["records"][0]["alts"])

    def test_crlf_line_endings_are_accepted(self):
        created = self.sample(vcf=VCF.replace("\n", "\r\n"))
        self.assertEqual(7, created["records"])

    def test_missing_trailing_newline_is_accepted(self):
        created = self.sample(vcf=VCF.rstrip("\n"))
        self.assertEqual(7, created["records"])

    def test_invalid_vcf_reports_the_line_number(self):
        header_only = PREAMBLE + "\n" + HEADER
        cases = [
            ("chr1\t1\t.\tA\tG\t.\tPASS\t.", "line 1"),
            ("##reference=GRCh38\n" + HEADER + "\nchr1\t1\t.\tA\tG\t.\tPASS\t.", "line 1"),
            ("##fileformat=VCFv3.0\n" + HEADER + "\nchr1\t1\t.\tA\tG\t.\tPASS\t.", "line 1"),
            (PREAMBLE + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\nchr1\t1\t.\tA\tG\t.\tPASS\t.", "line 4"),
            (header_only + "\nchr1\t1\t.\tA\tG\t.\tPASS", "line 5"),
            (header_only + "\nchr1\t0\t.\tA\tG\t.\tPASS\t.", "line 5"),
            (header_only + "\nchr1\t1\t.\tA\tG\tabc\tPASS\t.", "line 5"),
            (header_only + "\nchr1\t1\t.\tA\t<DEL>\t5\tPASS\t.", "line 5"),
            (header_only + "\nchr1\t1\t.\tA\tG\t5\tPASS\tDP=1;DP=2", "line 5"),
            (header_only + "\nchr1\t1\t.\tA\tG\t5\tPASS\tDP=x", "line 5"),
            (header_only + "\nchr1\t1\t.\tA\tG\t5\tPASS;q10\t.", "line 5"),
            (header_only + "\n\nchr1\t1\t.\tA\tG\t5\tPASS\t.", "line 5"),
            (header_only + "\n##extra=1\nchr1\t1\t.\tA\tG\t5\tPASS\t.", "line 5"),
            ("##fileformat=VCFv4.2\n" + HEADER + "\n##extra=1\nchr1\t1\t.\tA\tG\t5\tPASS\t.", "line 3"),
        ]
        for vcf, pattern in cases:
            with self.subTest(vcf=vcf), self.assertRaisesRegex(ValidationError, pattern):
                self.sample(vcf=vcf, key=f"k-{pattern}-{len(vcf)}")

    def test_structural_vcf_failures_are_rejected(self):
        for vcf, pattern in [
            (PREAMBLE, "no #CHROM header"),
            (PREAMBLE + "\n" + HEADER + "\n", "no data records"),
            ("   \n", "must not be empty"),
        ]:
            with self.subTest(vcf=vcf), self.assertRaisesRegex(ValidationError, pattern):
                self.sample(vcf=vcf)

    def test_unordered_and_duplicate_variants_are_rejected(self):
        with self.assertRaisesRegex(ValidationError, "sorted by POS"):
            self.sample(vcf="\n".join([PREAMBLE, HEADER, "chr1\t200\t.\tA\tG\t5\tPASS\t.", "chr1\t100\t.\tA\tG\t5\tPASS\t."]))
        with self.assertRaisesRegex(ValidationError, "duplicate variant"):
            self.sample(vcf="\n".join([PREAMBLE, HEADER, "chr1\t100\t.\tA\tG\t5\tPASS\t.", "chr1\t100\t.\tA\tG\t5\tPASS\t."]))
        with self.assertRaisesRegex(ValidationError, "more than one block"):
            self.sample(
                vcf="\n".join(
                    [
                        PREAMBLE,
                        HEADER,
                        "chr1\t100\t.\tA\tG\t5\tPASS\t.",
                        "chr2\t100\t.\tA\tG\t5\tPASS\t.",
                        "chr1\t200\t.\tA\tG\t5\tPASS\t.",
                    ]
                )
            )

    def test_sample_body_and_identifier_are_validated(self):
        with self.assertRaisesRegex(ValidationError, "exactly id and vcf"):
            self.service.create_sample({"id": "s1", "vcf": VCF, "extra": 1}, "k1")
        with self.assertRaisesRegex(ValidationError, "sample id"):
            self.service.create_sample({"id": "bad id", "vcf": VCF}, "k2")
        with self.assertRaisesRegex(ValidationError, "vcf must be a string"):
            self.service.create_sample({"id": "s1", "vcf": 5}, "k3")

    # --- runs --------------------------------------------------------------

    def test_default_run_keeps_every_record(self):
        self.sample()
        run = self.start_run()
        statistics = run["statistics"]
        self.assertEqual("succeeded", run["status"])
        self.assertEqual(7, statistics["records_total"])
        self.assertEqual(7, statistics["records_kept"])
        self.assertEqual(8, statistics["alleles_total"])
        self.assertEqual(8, statistics["alleles_kept"])
        self.assertEqual({"allele_filtered": 0, "dp_below_min": 0, "not_pass": 0, "qual_below_min": 0}, statistics["filter_counts"])

    def test_run_applies_qual_dp_and_filter_thresholds(self):
        self.sample()
        run = self.start_run({"pass_only": True, "min_qual": 30, "min_dp": 10})
        statistics = run["statistics"]
        self.assertEqual(4, statistics["records_kept"])
        self.assertEqual(3, statistics["records_filtered"])
        self.assertEqual({"allele_filtered": 0, "dp_below_min": 1, "not_pass": 1, "qual_below_min": 1}, statistics["filter_counts"])
        self.assertEqual(5, statistics["alleles_kept"])
        self.assertEqual({"frameshift_variant": 1, "missense_variant": 3, "stop_gained": 1}, statistics["consequence_counts"])
        self.assertEqual({"HIGH": 2, "MODERATE": 3}, statistics["impact_counts"])
        self.assertEqual({"BRCA1": 1, "EGFR": 1, "KRAS": 2, "MTHFR": 1}, statistics["gene_counts"])
        self.assertEqual(statistics["records_filtered"], sum(statistics["filter_counts"].values()))

    def test_variants_are_annotated_in_file_order(self):
        self.sample()
        run = self.start_run({"pass_only": True, "min_qual": 30, "min_dp": 10})
        payload = self.service.run_variants(run["id"])
        self.assertEqual(4, payload["count"])
        self.assertEqual([5, 8, 10, 11], [variant["line"] for variant in payload["variants"]])
        self.assertEqual(
            [
                {"allele": "A", "allele_index": 1, "consequence": "missense_variant", "gene": "MTHFR", "impact": "MODERATE"}
            ],
            payload["variants"][0]["annotations"],
        )
        self.assertEqual(
            [
                {"allele": "A", "allele_index": 1, "consequence": "missense_variant", "gene": "KRAS", "impact": "MODERATE"},
                {"allele": "T", "allele_index": 2, "consequence": "missense_variant", "gene": "KRAS", "impact": "MODERATE"},
            ],
            payload["variants"][2]["annotations"],
        )

    def test_original_line_is_preserved_on_every_variant(self):
        self.sample()
        run = self.start_run()
        for variant in self.service.run_variants(run["id"])["variants"]:
            self.assertEqual(variant["raw"].split("\t")[0], variant["chrom"])
            self.assertEqual(int(variant["raw"].split("\t")[1]), variant["pos"])

    def test_untabled_variants_fall_back_to_intergenic_annotation(self):
        self.sample(vcf=single(pos=11856379, ref="A", alt="G", qual=".", filter_value=".", info="."))
        run = self.start_run()
        annotations = self.service.run_variants(run["id"])["variants"][0]["annotations"]
        self.assertEqual(
            [{"allele": "G", "allele_index": 1, "consequence": "intergenic_variant", "gene": None, "impact": "MODIFIER"}],
            annotations,
        )
        self.assertEqual({"NA": 1}, run["statistics"]["gene_counts"])
        self.assertIsNone(self.service.run_variants(run["id"])["variants"][0]["qual"])
        self.assertEqual([], self.service.run_variants(run["id"])["variants"][0]["filter"])

    def test_gene_allowlist_filters_alleles_and_drops_empty_records(self):
        self.sample()
        run = self.start_run({"genes": ["EGFR"]})
        statistics = run["statistics"]
        self.assertEqual(2, statistics["records_kept"])
        self.assertEqual(5, statistics["records_filtered"])
        self.assertEqual(5, statistics["filter_counts"]["allele_filtered"])
        self.assertEqual(statistics["records_filtered"], sum(statistics["filter_counts"].values()))
        self.assertEqual({"EGFR": 2}, statistics["gene_counts"])

    def test_impact_allowlist_keeps_only_matching_alleles(self):
        self.sample()
        run = self.start_run({"impacts": ["HIGH"]})
        variants = self.service.run_variants(run["id"])["variants"]
        self.assertEqual(4, len(variants))
        self.assertEqual([6, 8, 9, 11], [variant["line"] for variant in variants])
        self.assertEqual({"HIGH": 4}, run["statistics"]["impact_counts"])
        self.assertEqual(4, run["statistics"]["alleles_kept"])
        self.assertEqual(3, run["statistics"]["records_filtered"])
        self.assertEqual(
            {"allele_filtered": 3, "dp_below_min": 0, "not_pass": 0, "qual_below_min": 0},
            run["statistics"]["filter_counts"],
        )

    def test_partial_allele_filter_keeps_the_record(self):
        vcf = "\n".join([PREAMBLE, HEADER, "chr17\t43093445\t.\tC\tT,G\t80\tPASS\tDP=60"]) + "\n"
        self.sample(vcf=vcf)
        run = self.start_run({"impacts": ["HIGH"]})
        statistics = run["statistics"]
        self.assertEqual(1, statistics["records_kept"])
        self.assertEqual(2, statistics["alleles_total"])
        self.assertEqual(1, statistics["alleles_kept"])
        self.assertEqual(0, statistics["filter_counts"]["allele_filtered"])
        variant = self.service.run_variants(run["id"])["variants"][0]
        self.assertEqual(["T", "G"], variant["alts"])
        self.assertEqual(["T"], [entry["allele"] for entry in variant["annotations"]])
        self.assertEqual({"BRCA1": 1}, statistics["gene_counts"])

    def test_params_validation(self):
        self.sample()
        for params, pattern in [
            ({"min_qual": -1}, "min_qual"),
            ({"min_qual": "30"}, "min_qual"),
            ({"min_qual": True}, "min_qual"),
            ({"min_dp": 1.5}, "min_dp"),
            ({"pass_only": "yes"}, "pass_only"),
            ({"genes": "EGFR"}, "params.genes"),
            ({"genes": [""]}, "params.genes"),
            ({"genes": ["EGFR", "EGFR"]}, "params.genes"),
            ({"impacts": ["high"]}, "params.impacts"),
            ({"depth": 3}, "unknown fields"),
        ]:
            with self.subTest(params=params), self.assertRaisesRegex(ValidationError, pattern):
                self.start_run(params=params, key=f"k-{sorted(params)[0]}")

    def test_run_body_must_be_exactly_id_and_params(self):
        self.sample()
        with self.assertRaisesRegex(ValidationError, "exactly id and optional params"):
            self.service.create_run("s1", {"id": "run-1", "min_qual": 5}, "k1")
        with self.assertRaisesRegex(ValidationError, "exactly id and optional params"):
            self.service.create_run("s1", {"params": {}}, "k2")

    def test_run_and_sample_lookups_are_not_found(self):
        with self.assertRaises(NotFoundError):
            self.service.create_run("missing", {"id": "run-1"}, "k1")
        with self.assertRaises(NotFoundError):
            self.service.get_run("missing")
        with self.assertRaises(NotFoundError):
            self.service.run_variants("missing")
        with self.assertRaises(NotFoundError):
            self.service.run_provenance("missing")

    def test_duplicate_run_id_conflicts(self):
        self.sample()
        self.start_run(key="k1")
        with self.assertRaisesRegex(ConflictError, "already exists"):
            self.start_run(params={"min_qual": 1}, key="k2")

    # --- idempotency -------------------------------------------------------

    def test_idempotency_key_is_required(self):
        with self.assertRaisesRegex(ValidationError, "Idempotency-Key"):
            self.service.create_sample({"id": "s1", "vcf": VCF}, None)
        self.sample()
        with self.assertRaisesRegex(ValidationError, "Idempotency-Key"):
            self.service.create_run("s1", {"id": "run-1"}, None)

    def test_repeated_key_returns_the_first_result(self):
        first = self.sample(key="same")
        second = self.sample(key="same")
        self.assertEqual(first, second)
        run_first = self.start_run(params={"min_qual": 30}, key="rk")
        run_second = self.start_run(params={"min_qual": 999}, key="rk")
        self.assertEqual(run_first, run_second)
        self.assertEqual(30.0, run_second["params"]["min_qual"])

    def test_key_cannot_be_reused_for_another_operation(self):
        self.sample(key="shared")
        with self.assertRaises(ConflictError):
            self.start_run(key="shared")
        with self.assertRaises(ConflictError):
            self.sample(key="shared", sample_id="s2")

    # --- provenance --------------------------------------------------------

    def test_provenance_has_four_chained_steps(self):
        self.sample()
        run = self.start_run({"pass_only": True, "min_qual": 30})
        payload = self.service.run_provenance(run["id"])
        self.assertTrue(payload["verified"])
        self.assertTrue(payload["chain_verified"])
        self.assertTrue(payload["reproduction_verified"])
        self.assertEqual(4, len(payload["steps"]))
        self.assertEqual(["ingest", "filter", "annotate", "summarize"], [entry["name"] for entry in payload["steps"]])
        self.assertIsNone(payload["steps"][0]["previous_sha256"])
        self.assertEqual(run["provenance_head"], payload["head_sha256"])
        self.assertEqual(run["provenance_length"], 4)
        for previous, current in zip(payload["steps"], payload["steps"][1:]):
            self.assertEqual(previous["hash"], current["previous_sha256"])
            self.assertEqual(previous["output_sha256"], current["input_sha256"])
        self.assertEqual(64, len(payload["head_sha256"]))

    def test_same_input_and_params_reproduce_identical_hashes(self):
        self.sample(key="k1", sample_id="s1")
        self.sample(key="k2", sample_id="s2")
        first = self.start_run(params={"min_dp": 10}, key="r1", run_id="run-1", sample_id="s1")
        second = self.start_run(params={"min_dp": 10}, key="r2", run_id="run-2", sample_id="s2")
        self.assertEqual(first["statistics"], second["statistics"])
        self.assertEqual(first["provenance_head"], second["provenance_head"])
        steps_first = self.service.run_provenance("run-1")["steps"]
        steps_second = self.service.run_provenance("run-2")["steps"]
        self.assertEqual(steps_first, steps_second)

    def test_ingest_hash_covers_the_original_text(self):
        self.sample()
        run = self.start_run()
        steps = self.service.run_provenance(run["id"])["steps"]
        self.assertEqual(self.service.get_sample("s1")["sha256"], steps[0]["input_sha256"])
        self.assertNotEqual(steps[0]["input_sha256"], steps[0]["output_sha256"])

    def test_annotation_step_binds_the_table_content(self):
        from variantrail.annotate import table_sha256

        self.sample()
        run = self.start_run()
        steps = self.service.run_provenance(run["id"])["steps"]
        self.assertEqual(table_sha256(), steps[2]["params"]["annotation_table_sha256"])

    def test_tampered_chain_is_detected(self):
        self.sample()
        self.start_run()
        document = self.service._document("run-1")
        document["provenance"][1]["params"]["min_qual"] = 999.0
        self.service.store.connection.execute(
            "UPDATE runs SET document = ? WHERE id = ?", (self.service.store.encode(document), "run-1")
        )
        payload = self.service.run_provenance("run-1")
        self.assertFalse(payload["chain_verified"])
        self.assertFalse(payload["verified"])

    def test_tampered_statistics_are_detected_by_reproduction(self):
        self.sample()
        self.start_run()
        document = self.service._document("run-1")
        document["statistics"]["records_kept"] = 999
        self.service.store.connection.execute(
            "UPDATE runs SET document = ? WHERE id = ?", (self.service.store.encode(document), "run-1")
        )
        payload = self.service.run_provenance("run-1")
        self.assertTrue(payload["chain_verified"])
        self.assertFalse(payload["reproduction_verified"])
        self.assertFalse(payload["verified"])

    def test_re_running_does_not_change_earlier_hashes(self):
        self.sample()
        first = self.start_run(params={"min_qual": 30}, key="r1", run_id="run-1")
        second = self.start_run(params={"min_qual": 30}, key="r2", run_id="run-2")
        self.assertEqual(first["provenance_head"], second["provenance_head"])
        self.assertEqual(
            self.service.run_provenance("run-1")["steps"],
            self.service.run_provenance("run-2")["steps"],
        )

    # --- comparison --------------------------------------------------------

    def test_compare_partitions_alleles_into_shared_and_left_only(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(params={"impacts": ["HIGH"]}, key="r2", run_id="run-2")
        payload = self.service.compare_runs("run-1", "run-2")
        self.assertEqual(
            {"left_only": 4, "right_only": 0, "shared": 4, "union": 8},
            payload["counts"],
        )
        self.assertEqual(8, payload["left"]["allele_count"])
        self.assertEqual(4, payload["right"]["allele_count"])
        self.assertEqual("run-1", payload["left"]["run_id"])
        self.assertEqual("run-2", payload["right"]["run_id"])
        self.assertEqual("s1", payload["left"]["sample_id"])
        self.assertEqual(
            self.service.get_run("run-1")["sample_sha256"], payload["left"]["sample_sha256"]
        )
        self.assertEqual(
            [
                ("chr1", 11856378, "G", "T"),
                ("chr17", 43093445, "C", "T"),
                ("chr7", 55019017, "G", "GA"),
                ("chr7", 55019017, "G", "T"),
            ],
            [(a["chrom"], a["pos"], a["ref"], a["alt"]) for a in payload["shared"]],
        )
        self.assertEqual(
            [
                ("chr1", 11856378, "G", "A"),
                ("chr1", 11856379, "A", "G"),
                ("chr12", 25245350, "C", "A"),
                ("chr12", 25245350, "C", "T"),
            ],
            [(a["chrom"], a["pos"], a["ref"], a["alt"]) for a in payload["left_only"]],
        )
        self.assertEqual([], payload["right_only"])
        entry = payload["shared"][0]
        self.assertEqual(
            {"alt", "chrom", "consequence", "gene", "impact", "pos", "ref"}, set(entry)
        )
        self.assertEqual("MTHFR", entry["gene"])
        self.assertEqual("stop_gained", entry["consequence"])
        self.assertEqual("HIGH", entry["impact"])
        self.assertIsNone(payload["left_only"][1]["gene"])

    def test_compare_summaries_group_by_gene_and_impact(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(params={"impacts": ["HIGH"]}, key="r2", run_id="run-2")
        payload = self.service.compare_runs("run-1", "run-2")
        self.assertEqual(
            {
                "BRCA1": {"left_only": 0, "right_only": 0, "shared": 1},
                "EGFR": {"left_only": 0, "right_only": 0, "shared": 2},
                "KRAS": {"left_only": 2, "right_only": 0, "shared": 0},
                "MTHFR": {"left_only": 1, "right_only": 0, "shared": 1},
                "NA": {"left_only": 1, "right_only": 0, "shared": 0},
            },
            payload["gene_summary"],
        )
        self.assertEqual(
            {
                "HIGH": {"left_only": 0, "right_only": 0, "shared": 4},
                "MODERATE": {"left_only": 3, "right_only": 0, "shared": 0},
                "MODIFIER": {"left_only": 1, "right_only": 0, "shared": 0},
            },
            payload["impact_summary"],
        )
        self.assertEqual(sorted(payload["gene_summary"]), list(payload["gene_summary"]))
        self.assertEqual(sorted(payload["impact_summary"]), list(payload["impact_summary"]))

    def test_compare_is_symmetric_and_byte_stable(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(params={"impacts": ["HIGH"]}, key="r2", run_id="run-2")
        forward = self.service.compare_runs("run-1", "run-2")
        reverse = self.service.compare_runs("run-2", "run-1")
        self.assertEqual(forward["shared"], reverse["shared"])
        self.assertEqual(forward["left_only"], reverse["right_only"])
        self.assertEqual(forward["right_only"], reverse["left_only"])
        self.assertEqual(forward["counts"]["union"], reverse["counts"]["union"])
        self.assertEqual(forward, self.service.compare_runs("run-1", "run-2"))

    def test_compare_same_run_reports_everything_shared(self):
        self.sample()
        self.start_run()
        payload = self.service.compare_runs("run-1", "run-1")
        self.assertEqual(
            {"left_only": 0, "right_only": 0, "shared": 8, "union": 8},
            payload["counts"],
        )
        self.assertEqual([], payload["left_only"])
        self.assertEqual([], payload["right_only"])
        self.assertEqual(8, len(payload["shared"]))

    def test_compare_works_across_samples(self):
        self.sample(key="s1", sample_id="s1")
        self.sample(key="s2", sample_id="s2", vcf=single(pos=5227002, chrom="chr11", ref="T", alt="C"))
        self.start_run(key="r1", run_id="run-1", sample_id="s1")
        self.start_run(key="r2", run_id="run-2", sample_id="s2")
        payload = self.service.compare_runs("run-1", "run-2")
        self.assertEqual("s1", payload["left"]["sample_id"])
        self.assertEqual("s2", payload["right"]["sample_id"])
        self.assertNotEqual(payload["left"]["sample_sha256"], payload["right"]["sample_sha256"])
        self.assertEqual({"left_only": 8, "right_only": 1, "shared": 0, "union": 9}, payload["counts"])
        self.assertEqual("HBB", payload["right_only"][0]["gene"])

    def test_compare_missing_run_is_not_found(self):
        self.sample()
        self.start_run()
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.compare_runs("run-1", "missing")
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.compare_runs("missing", "run-1")

    def test_compare_does_not_modify_runs_or_provenance(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(params={"impacts": ["HIGH"]}, key="r2", run_id="run-2")
        before_left = self.service._document("run-1")
        before_right = self.service._document("run-2")
        self.service.compare_runs("run-1", "run-2")
        self.assertEqual(before_left, self.service._document("run-1"))
        self.assertEqual(before_right, self.service._document("run-2"))
        self.assertTrue(self.service.run_provenance("run-1")["verified"])
        self.assertTrue(self.service.run_provenance("run-2")["verified"])

    # --- cohort comparison -------------------------------------------------

    def _three_runs(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(params={"impacts": ["HIGH"]}, key="r2", run_id="run-2")
        self.start_run(params={"genes": ["EGFR"]}, key="r3", run_id="run-3")

    def test_cohort_reports_runs_in_path_order(self):
        self._three_runs()
        payload = self.service.cohort_runs("run-1", "run-2,run-3")
        self.assertEqual(3, payload["run_count"])
        self.assertEqual(
            [
                {"run_id": "run-1", "sample_id": "s1", "sample_sha256": self.service.get_run("run-1")["sample_sha256"], "allele_count": 8},
                {"run_id": "run-2", "sample_id": "s1", "sample_sha256": self.service.get_run("run-2")["sample_sha256"], "allele_count": 4},
                {"run_id": "run-3", "sample_id": "s1", "sample_sha256": self.service.get_run("run-3")["sample_sha256"], "allele_count": 2},
            ],
            payload["runs"],
        )

    def test_cohort_partitions_alleles_by_presence(self):
        self._three_runs()
        payload = self.service.cohort_runs("run-1", "run-2,run-3")
        self.assertEqual(8, len(payload["alleles"]))
        self.assertEqual(
            [
                ("chr1", 11856378, "G", "A"),
                ("chr1", 11856378, "G", "T"),
                ("chr1", 11856379, "A", "G"),
                ("chr12", 25245350, "C", "A"),
                ("chr12", 25245350, "C", "T"),
                ("chr17", 43093445, "C", "T"),
                ("chr7", 55019017, "G", "GA"),
                ("chr7", 55019017, "G", "T"),
            ],
            [(a["chrom"], a["pos"], a["ref"], a["alt"]) for a in payload["alleles"]],
        )
        by_key = {(a["chrom"], a["pos"], a["ref"], a["alt"]): a for a in payload["alleles"]}
        self.assertEqual(1, by_key[("chr1", 11856378, "G", "A")]["run_count"])
        self.assertEqual(["run-1"], by_key[("chr1", 11856378, "G", "A")]["present_in"])
        self.assertEqual(2, by_key[("chr1", 11856378, "G", "T")]["run_count"])
        self.assertEqual(["run-1", "run-2"], by_key[("chr1", 11856378, "G", "T")]["present_in"])
        self.assertEqual(3, by_key[("chr7", 55019017, "G", "GA")]["run_count"])
        self.assertEqual(["run-1", "run-2", "run-3"], by_key[("chr7", 55019017, "G", "GA")]["present_in"])
        entry = by_key[("chr1", 11856378, "G", "T")]
        self.assertEqual(
            {"alt", "chrom", "consequence", "gene", "impact", "pos", "present_in", "ref", "run_count"},
            set(entry),
        )
        self.assertEqual("MTHFR", entry["gene"])
        self.assertEqual("stop_gained", entry["consequence"])
        self.assertEqual("HIGH", entry["impact"])
        self.assertIsNone(by_key[("chr1", 11856379, "A", "G")]["gene"])

    def test_cohort_counts_and_frequency_summary(self):
        self._three_runs()
        payload = self.service.cohort_runs("run-1", "run-2,run-3")
        self.assertEqual(
            {
                "unique_alleles": 8,
                "core_alleles": 2,
                "variable_alleles": 2,
                "private_alleles": {"run-1": 4, "run-2": 0, "run-3": 0},
            },
            payload["counts"],
        )
        self.assertEqual({"1": 4, "2": 2, "3": 2}, payload["frequency_summary"])
        self.assertEqual(["1", "2", "3"], list(payload["frequency_summary"]))

    def test_cohort_frequency_summary_fills_missing_buckets_with_zero(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        self.start_run(key="r2", run_id="run-2")
        self.start_run(key="r3", run_id="run-3")
        payload = self.service.cohort_runs("run-1", "run-2,run-3")
        self.assertEqual(8, payload["counts"]["core_alleles"])
        self.assertEqual(0, payload["counts"]["variable_alleles"])
        self.assertEqual({"1": 0, "2": 0, "3": 8}, payload["frequency_summary"])

    def test_cohort_summaries_group_by_gene_and_impact(self):
        self._three_runs()
        payload = self.service.cohort_runs("run-1", "run-2,run-3")
        self.assertEqual(
            {
                "BRCA1": {"unique_alleles": 1, "core_alleles": 0, "variable_alleles": 1},
                "EGFR": {"unique_alleles": 2, "core_alleles": 2, "variable_alleles": 0},
                "KRAS": {"unique_alleles": 2, "core_alleles": 0, "variable_alleles": 0},
                "MTHFR": {"unique_alleles": 2, "core_alleles": 0, "variable_alleles": 1},
                "NA": {"unique_alleles": 1, "core_alleles": 0, "variable_alleles": 0},
            },
            payload["gene_summary"],
        )
        self.assertEqual(
            {
                "HIGH": {"unique_alleles": 4, "core_alleles": 2, "variable_alleles": 2},
                "MODERATE": {"unique_alleles": 3, "core_alleles": 0, "variable_alleles": 0},
                "MODIFIER": {"unique_alleles": 1, "core_alleles": 0, "variable_alleles": 0},
            },
            payload["impact_summary"],
        )
        self.assertEqual(sorted(payload["gene_summary"]), list(payload["gene_summary"]))
        self.assertEqual(sorted(payload["impact_summary"]), list(payload["impact_summary"]))

    def test_cohort_reordering_only_changes_order_dependent_fields(self):
        self._three_runs()
        forward = self.service.cohort_runs("run-1", "run-2,run-3")
        reordered = self.service.cohort_runs("run-2", "run-3,run-1")

        def normalize(payload):
            return {
                key: (
                    [{**allele, "present_in": sorted(allele["present_in"])} for allele in payload[key]]
                    if key == "alleles"
                    else payload[key]
                )
                for key in payload
                if key != "runs"
            }

        self.assertEqual(normalize(forward), normalize(reordered))
        self.assertEqual(["run-2", "run-3", "run-1"], [run["run_id"] for run in reordered["runs"]])
        by_key = {(a["chrom"], a["pos"], a["ref"], a["alt"]): a for a in reordered["alleles"]}
        self.assertEqual(["run-2", "run-3", "run-1"], by_key[("chr7", 55019017, "G", "GA")]["present_in"])
        self.assertEqual(["run-1"], by_key[("chr12", 25245350, "C", "A")]["present_in"])
        self.assertEqual(["run-2", "run-1"], by_key[("chr17", 43093445, "C", "T")]["present_in"])

    def test_cohort_is_byte_stable_across_requests(self):
        self._three_runs()
        self.assertEqual(
            self.service.cohort_runs("run-1", "run-2,run-3"),
            self.service.cohort_runs("run-1", "run-2,run-3"),
        )

    def test_cohort_requires_two_runs(self):
        self.sample()
        self.start_run()
        with self.assertRaisesRegex(ValidationError, "cohort requires at least two distinct runs"):
            self.service.cohort_runs("run-1", None)

    def test_cohort_rejects_duplicate_run_ids(self):
        self._three_runs()
        with self.assertRaisesRegex(ConflictError, "cohort run ids must be distinct"):
            self.service.cohort_runs("run-1", "run-1")
        with self.assertRaisesRegex(ConflictError, "cohort run ids must be distinct"):
            self.service.cohort_runs("run-1", "run-2,run-2")
        with self.assertRaisesRegex(ConflictError, "cohort run ids must be distinct"):
            self.service.cohort_runs("run-2", "run-1,run-2")

    def test_cohort_missing_run_is_not_found(self):
        self._three_runs()
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.cohort_runs("run-1", "missing")
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.cohort_runs("missing", "run-1")
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.cohort_runs("run-1", "run-2,missing")

    def test_cohort_validates_shape_before_existence(self):
        self.sample()
        self.start_run(key="r1", run_id="run-1")
        with self.assertRaisesRegex(ValidationError, "cohort requires at least two distinct runs"):
            self.service.cohort_runs("missing", None)
        with self.assertRaisesRegex(ConflictError, "cohort run ids must be distinct"):
            self.service.cohort_runs("missing", "missing")

    def test_cohort_does_not_modify_runs_or_provenance(self):
        self._three_runs()
        before = {identifier: self.service._document(identifier) for identifier in ("run-1", "run-2", "run-3")}
        self.service.cohort_runs("run-3", "run-1,run-2")
        for identifier, document in before.items():
            self.assertEqual(document, self.service._document(identifier))
            self.assertTrue(self.service.run_provenance(identifier)["verified"])

    # --- exports -----------------------------------------------------------

    def test_jsonl_export_has_one_canonical_object_per_line_in_file_order(self):
        self.sample()
        self.start_run()
        variants = self.service.run_variants("run-1")["variants"]
        content_type, body = self.service.run_export("run-1", "jsonl")
        self.assertEqual("application/x-ndjson; charset=utf-8", content_type)
        text = body.decode("utf-8")
        lines = text.split("\n")
        self.assertEqual("", lines[-1])  # exactly one trailing LF, no blank line
        self.assertEqual(len(variants), len(lines) - 1)
        self.assertNotIn("\r", text)
        for line, variant in zip(lines[:-1], variants):
            self.assertEqual(json.dumps(variant, ensure_ascii=False, separators=(",", ":"), sort_keys=True), line)
            self.assertEqual(variant, json.loads(line))
        self.assertEqual(variants[0], json.loads(lines[0]))

    def test_tsv_export_has_fixed_header_and_one_row_per_annotation(self):
        self.sample()
        self.start_run()
        content_type, body = self.service.run_export("run-1", "tsv")
        self.assertEqual("text/tab-separated-values; charset=utf-8", content_type)
        lines = body.decode("utf-8").split("\n")
        self.assertEqual(TSV_HEADER, lines[0])
        self.assertEqual("", lines[-1])
        rows = [line.split("\t") for line in lines[1:-1]]
        self.assertEqual(8, len(rows))  # seven records, one with two ALTs
        self.assertEqual(
            ["chr1", "11856378", "rs1", "G", "A", "1", "MTHFR", "missense_variant", "MODERATE", "60.0", "42", "PASS", "DP=42", "5"],
            rows[0],
        )
        self.assertEqual(
            ["chr7", "55019017", "NA", "G", "GA", "1", "EGFR", "frameshift_variant", "HIGH", "50.0", "31", "PASS", "DP=31", "8"],
            rows[3],
        )
        # Multi-ALT record: rows stay in annotation order and share line number.
        self.assertEqual(["chr12", "25245350", "NA", "C", "A", "1", "KRAS", "missense_variant", "MODERATE", "80.0", "60", "PASS", "DP=60", "10"], rows[5])
        self.assertEqual(["chr12", "25245350", "NA", "C", "T", "2", "KRAS", "missense_variant", "MODERATE", "80.0", "60", "PASS", "DP=60", "10"], rows[6])
        # Missing DP on the chr7 stop-gained record.
        self.assertEqual("NA", rows[4][10])
        self.assertEqual([5, 6, 7, 8, 9, 10, 10, 11], [int(row[13]) for row in rows])

    def test_tsv_export_renders_nulls_empty_filter_and_sorted_flag_info(self):
        vcf = single(pos=11856379, ref="A", alt="G", qual=".", filter_value="q10;S50", info="DP=7;ZZ=1;AA=2;FLAG")
        self.sample(vcf=vcf)
        self.start_run()
        _, body = self.service.run_export("run-1", "tsv")
        row = body.decode("utf-8").split("\n")[1].split("\t")
        self.assertEqual(
            ["chr1", "11856379", "NA", "A", "G", "1", "NA", "intergenic_variant", "MODIFIER", "NA", "7", "q10;S50", "AA=2;DP=7;FLAG;ZZ=1", "5"],
            row,
        )

    def test_empty_run_exports_header_or_empty_body(self):
        self.sample()
        self.start_run({"genes": ["NO_SUCH_GENE"]})
        _, jsonl_body = self.service.run_export("run-1", "jsonl")
        _, tsv_body = self.service.run_export("run-1", "tsv")
        self.assertEqual(b"", jsonl_body)
        self.assertEqual((TSV_HEADER + "\n").encode("utf-8"), tsv_body)

    def test_export_validates_format_and_run_existence(self):
        self.sample()
        self.start_run()
        with self.assertRaisesRegex(ValidationError, "format"):
            self.service.run_export("run-1", "csv")
        with self.assertRaises(NotFoundError):
            self.service.run_export("missing", "tsv")
        # An invalid format is a 400 even when the run is also missing.
        with self.assertRaises(ValidationError):
            self.service.run_export("missing", "csv")

    def test_export_does_not_modify_the_run_or_provenance(self):
        self.sample()
        self.start_run()
        before = self.service._document("run-1")
        self.service.run_export("run-1", "tsv")
        self.service.run_export("run-1", "jsonl")
        self.service.run_export("run-1", "vcf")
        after = self.service._document("run-1")
        self.assertEqual(before, after)
        self.assertTrue(self.service.run_provenance("run-1")["verified"])

    def test_vcf_export_header_and_one_record_per_retained_allele(self):
        self.sample()
        self.start_run()
        content_type, body = self.service.run_export("run-1", "vcf")
        self.assertEqual("text/x-variant-call-format; charset=utf-8", content_type)
        text = body.decode("utf-8")
        self.assertNotIn("\r", text)
        lines = text.split("\n")
        self.assertEqual("", lines[-1])  # exactly one trailing LF
        self.assertEqual(
            [
                "##fileformat=VCFv4.2",
                "##reference=GRCh38",
                "##contig=<ID=chr1>",
                "##source=variantrail",
                '##INFO=<ID=ALLELE_INDEX,Number=1,Type=Integer,Description="1-based index of the ALT allele in the source record">',
                '##INFO=<ID=VT_GENE,Number=1,Type=String,Description="Annotated gene symbol, NA when unannotated">',
                '##INFO=<ID=VT_CONSEQUENCE,Number=1,Type=String,Description="Annotated sequence consequence">',
                '##INFO=<ID=VT_IMPACT,Number=1,Type=String,Description="Annotated impact">',
                HEADER,
            ],
            lines[:9],
        )
        rows = lines[9:-1]
        self.assertEqual(8, len(rows))  # seven records, one with two ALTs
        self.assertEqual(
            "chr1\t11856378\trs1\tG\tA\t60.0\tPASS\tDP=42;ALLELE_INDEX=1;VT_GENE=MTHFR;VT_CONSEQUENCE=missense_variant;VT_IMPACT=MODERATE",
            rows[0],
        )
        # Empty original INFO, missing ID.
        self.assertEqual(
            "chr7\t55019017\t.\tG\tT\t50.0\tPASS\tALLELE_INDEX=1;VT_GENE=EGFR;VT_CONSEQUENCE=stop_gained;VT_IMPACT=HIGH",
            rows[4],
        )
        # Multi-ALT record: one row per allele, allele_index order, ALT holds only the current allele.
        self.assertEqual(
            "chr12\t25245350\t.\tC\tA\t80.0\tPASS\tDP=60;ALLELE_INDEX=1;VT_GENE=KRAS;VT_CONSEQUENCE=missense_variant;VT_IMPACT=MODERATE",
            rows[5],
        )
        self.assertEqual(
            "chr12\t25245350\t.\tC\tT\t80.0\tPASS\tDP=60;ALLELE_INDEX=2;VT_GENE=KRAS;VT_CONSEQUENCE=missense_variant;VT_IMPACT=MODERATE",
            rows[6],
        )

    def test_vcf_export_renders_nulls_sorted_info_and_na_gene(self):
        vcf = single(pos=11856379, ref="A", alt="G", qual=".", filter_value="q10;S50", info="DP=7;ZZ=1;AA=2;FLAG")
        self.sample(vcf=vcf)
        self.start_run()
        _, body = self.service.run_export("run-1", "vcf")
        row = body.decode("utf-8").split("\n")[9]
        self.assertEqual(
            "chr1\t11856379\t.\tA\tG\t.\tq10;S50\tAA=2;DP=7;FLAG;ZZ=1;ALLELE_INDEX=1;VT_GENE=NA;VT_CONSEQUENCE=intergenic_variant;VT_IMPACT=MODIFIER",
            row,
        )

    def test_vcf_export_empty_run_returns_header_only(self):
        self.sample()
        self.start_run({"genes": ["NO_SUCH_GENE"]})
        status_body = self.service.run_export("run-1", "vcf")[1]
        text = status_body.decode("utf-8")
        self.assertTrue(text.endswith("\n"))
        self.assertNotIn("\n\n", text)
        lines = text.split("\n")
        self.assertEqual("", lines[-1])
        self.assertEqual(HEADER, lines[-2])
        self.assertEqual(9, len(lines) - 1)  # header lines plus column header only

    def test_vcf_export_is_byte_identical_across_requests(self):
        self.sample()
        self.start_run()
        _, first = self.service.run_export("run-1", "vcf")
        _, second = self.service.run_export("run-1", "vcf")
        self.assertEqual(first, second)

    def test_vcf_export_rejects_reserved_info_keys(self):
        for key in ("ALLELE_INDEX", "VT_GENE", "VT_CONSEQUENCE", "VT_IMPACT"):
            with self.subTest(key=key):
                self.sample(vcf=single(info=f"DP=42;{key}=1"), key=f"s-{key}", sample_id=f"s-{key}")
                self.start_run(key=f"r-{key}", run_id=f"run-{key}", sample_id=f"s-{key}")
                with self.assertRaisesRegex(ValidationError, key):
                    self.service.run_export(f"run-{key}", "vcf")

    def test_vcf_export_validates_format_and_run_existence(self):
        self.sample()
        self.start_run()
        with self.assertRaisesRegex(ValidationError, "jsonl, tsv, vcf"):
            self.service.run_export("run-1", "csv")
        with self.assertRaisesRegex(NotFoundError, "run missing was not found"):
            self.service.run_export("missing", "vcf")


class ExportHttpTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        service = VariantRail(str(Path(self.directory.name) / "test.db"))
        service.create_sample({"id": "s1", "vcf": VCF}, "s1")
        service.create_run("s1", {"id": "run-1"}, "r1")
        Handler.service = service
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

    def test_export_endpoints_serve_rendered_bodies(self):
        status, content_type, body = self.get("/runs/run-1/exports/jsonl")
        self.assertEqual(200, status)
        self.assertEqual("application/x-ndjson; charset=utf-8", content_type)
        self.assertEqual(7, len(body.decode("utf-8").strip().split("\n")))
        status, content_type, body = self.get("/runs/run-1/exports/tsv")
        self.assertEqual(200, status)
        self.assertEqual("text/tab-separated-values; charset=utf-8", content_type)
        self.assertEqual(TSV_HEADER, body.decode("utf-8").split("\n")[0])
        status, content_type, body = self.get("/runs/run-1/exports/vcf")
        self.assertEqual(200, status)
        self.assertEqual("text/x-variant-call-format; charset=utf-8", content_type)
        lines = body.decode("utf-8").split("\n")
        self.assertEqual("##fileformat=VCFv4.2", lines[0])
        self.assertEqual(HEADER, lines[8])
        self.assertEqual(8, len(lines) - 10)  # eight retained alleles, one trailing LF

    def test_export_errors_use_the_standard_error_object(self):
        status, content_type, body = self.get("/runs/missing/exports/tsv")
        self.assertEqual(404, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        self.assertEqual("not_found", json.loads(body)["error"]["code"])
        status, _, body = self.get("/runs/missing/exports/vcf")
        self.assertEqual(404, status)
        self.assertEqual("run missing was not found", json.loads(body)["error"]["message"])
        status, _, body = self.get("/runs/run-1/exports/csv")
        self.assertEqual(400, status)
        error = json.loads(body)["error"]
        self.assertEqual("validation_error", error["code"])
        self.assertEqual("format must be one of jsonl, tsv, vcf", error["message"])
        status, _, body = self.get("/runs/run-1/exports")
        self.assertEqual(404, status)
        self.assertEqual("not_found", json.loads(body)["error"]["code"])


class CompareHttpTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        service = VariantRail(str(Path(self.directory.name) / "test.db"))
        service.create_sample({"id": "s1", "vcf": VCF}, "s1")
        service.create_run("s1", {"id": "run-1"}, "r1")
        service.create_run("s1", {"id": "run-2", "params": {"impacts": ["HIGH"]}}, "r2")
        Handler.service = service
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

    def test_compare_endpoint_returns_the_comparison(self):
        status, content_type, body = self.get("/runs/run-1/compare/run-2")
        self.assertEqual(200, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        payload = json.loads(body)
        self.assertEqual({"left_only": 4, "right_only": 0, "shared": 4, "union": 8}, payload["counts"])
        self.assertEqual("run-1", payload["left"]["run_id"])
        self.assertEqual("run-2", payload["right"]["run_id"])
        self.assertEqual(4, len(payload["shared"]))
        self.assertEqual(4, len(payload["left_only"]))
        self.assertEqual([], payload["right_only"])
        # Repeated requests are byte-for-byte identical.
        _, _, again = self.get("/runs/run-1/compare/run-2")
        self.assertEqual(body, again)

    def test_compare_endpoint_errors_use_the_standard_error_object(self):
        status, content_type, body = self.get("/runs/run-1/compare/missing")
        self.assertEqual(404, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        error = json.loads(body)["error"]
        self.assertEqual("not_found", error["code"])
        self.assertEqual("run missing was not found", error["message"])
        status, _, body = self.get("/runs/missing/compare/run-1")
        self.assertEqual(404, status)
        self.assertEqual("run missing was not found", json.loads(body)["error"]["message"])
        status, _, body = self.get("/runs/run-1/compare")
        self.assertEqual(404, status)
        self.assertEqual("not_found", json.loads(body)["error"]["code"])


class CohortHttpTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        service = VariantRail(str(Path(self.directory.name) / "test.db"))
        service.create_sample({"id": "s1", "vcf": VCF}, "s1")
        service.create_run("s1", {"id": "run-1"}, "r1")
        service.create_run("s1", {"id": "run-2", "params": {"impacts": ["HIGH"]}}, "r2")
        service.create_run("s1", {"id": "run-3", "params": {"genes": ["EGFR"]}}, "r3")
        Handler.service = service
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

    def test_cohort_endpoint_returns_the_cohort(self):
        status, content_type, body = self.get("/runs/run-1/cohort/run-2,run-3")
        self.assertEqual(200, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        payload = json.loads(body)
        self.assertEqual(3, payload["run_count"])
        self.assertEqual(["run-1", "run-2", "run-3"], [run["run_id"] for run in payload["runs"]])
        self.assertEqual([8, 4, 2], [run["allele_count"] for run in payload["runs"]])
        self.assertEqual(8, len(payload["alleles"]))
        self.assertEqual(
            {
                "unique_alleles": 8,
                "core_alleles": 2,
                "variable_alleles": 2,
                "private_alleles": {"run-1": 4, "run-2": 0, "run-3": 0},
            },
            payload["counts"],
        )
        self.assertEqual({"1": 4, "2": 2, "3": 2}, payload["frequency_summary"])
        self.assertIn("gene_summary", payload)
        self.assertIn("impact_summary", payload)
        # Repeated requests are byte-for-byte identical.
        _, _, again = self.get("/runs/run-1/cohort/run-2,run-3")
        self.assertEqual(body, again)
        # The two-run queue uses the cohort shape as well.
        status, _, body = self.get("/runs/run-2/cohort/run-1")
        self.assertEqual(200, status)
        payload = json.loads(body)
        self.assertEqual(2, payload["run_count"])
        self.assertEqual(["run-2", "run-1"], [run["run_id"] for run in payload["runs"]])
        self.assertEqual({"1": 4, "2": 4}, payload["frequency_summary"])

    def test_cohort_endpoint_errors_use_the_standard_error_object(self):
        status, content_type, body = self.get("/runs/run-1/cohort")
        self.assertEqual(400, status)
        self.assertEqual("application/json; charset=utf-8", content_type)
        error = json.loads(body)["error"]
        self.assertEqual("validation_error", error["code"])
        self.assertEqual("cohort requires at least two distinct runs", error["message"])
        status, _, body = self.get("/runs/run-1/cohort/run-2,run-1")
        self.assertEqual(409, status)
        error = json.loads(body)["error"]
        self.assertEqual("conflict", error["code"])
        self.assertEqual("cohort run ids must be distinct", error["message"])
        status, _, body = self.get("/runs/run-1/cohort/run-2,missing")
        self.assertEqual(404, status)
        error = json.loads(body)["error"]
        self.assertEqual("not_found", error["code"])
        self.assertEqual("run missing was not found", error["message"])
        status, _, body = self.get("/runs/missing/cohort/run-1,run-2")
        self.assertEqual(404, status)
        self.assertEqual("run missing was not found", json.loads(body)["error"]["message"])


if __name__ == "__main__":
    unittest.main()

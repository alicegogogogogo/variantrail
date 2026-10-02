import tempfile
import unittest
from pathlib import Path

from variantrail.errors import ConflictError, NotFoundError, ValidationError
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

    # --- exports -----------------------------------------------------------

    def test_jsonl_export_matches_run_variants_line_for_line(self):
        import json

        self.sample()
        run = self.start_run()
        content_type, body = self.service.run_export(run["id"], "jsonl")
        self.assertEqual("application/x-ndjson; charset=utf-8", content_type)
        variants = self.service.run_variants(run["id"])["variants"]
        lines = body.split("\n")
        self.assertEqual(len(variants) + 1, len(lines))
        self.assertEqual("", lines[-1])
        self.assertNotIn("\n\n", body)
        self.assertNotIn("\r", body)
        for line, variant in zip(lines, variants):
            self.assertEqual(json.dumps(variant, ensure_ascii=False, sort_keys=True, separators=(",", ":")), line)
            self.assertEqual(list(sorted(variant)), list(json.loads(line)))

    def test_tsv_export_splits_alleles_and_formats_nulls(self):
        self.sample()
        run = self.start_run()
        content_type, body = self.service.run_export(run["id"], "tsv")
        self.assertEqual("text/tab-separated-values; charset=utf-8", content_type)
        lines = body.split("\n")
        self.assertEqual("", lines[-1])
        self.assertEqual(
            "chrom\tpos\tid\tref\talt\tallele_index\tgene\tconsequence\timpact\tqual\tdp\tfilter\tinfo\tline",
            lines[0],
        )
        self.assertEqual(
            [
                "chr1\t11856378\trs1\tG\tA\t1\tMTHFR\tmissense_variant\tMODERATE\t60.0\t42\tPASS\tDP=42\t5",
                "chr1\t11856378\trs2\tG\tT\t1\tMTHFR\tstop_gained\tHIGH\t20.0\t42\tPASS\tDP=42\t6",
                "chr1\t11856379\tNA\tA\tG\t1\tNA\tintergenic_variant\tMODIFIER\tNA\t5\tq10\tDP=5\t7",
                "chr7\t55019017\tNA\tG\tGA\t1\tEGFR\tframeshift_variant\tHIGH\t50.0\t31\tPASS\tDP=31\t8",
                "chr7\t55019017\tNA\tG\tT\t1\tEGFR\tstop_gained\tHIGH\t50.0\tNA\tPASS\t.\t9",
                "chr12\t25245350\tNA\tC\tA\t1\tKRAS\tmissense_variant\tMODERATE\t80.0\t60\tPASS\tDP=60\t10",
                "chr12\t25245350\tNA\tC\tT\t2\tKRAS\tmissense_variant\tMODERATE\t80.0\t60\tPASS\tDP=60\t10",
                "chr17\t43093445\tNA\tC\tT\t1\tBRCA1\tstop_gained\tHIGH\t99.0\t55\tPASS\tDP=55\t11",
            ],
            lines[1:-1],
        )

    def test_tsv_export_sorts_info_keys_and_keeps_flags_and_filters(self):
        vcf = "\n".join(
            [
                PREAMBLE,
                HEADER,
                "chr1\t11856378\trs1\tG\tA\t60\tq10;lowDP\tDP=42;FLAG;AB=0.5",
            ]
        ) + "\n"
        self.sample(vcf=vcf)
        run = self.start_run()
        _, body = self.service.run_export(run["id"], "tsv")
        row = body.split("\n")[1]
        self.assertEqual(
            "chr1\t11856378\trs1\tG\tA\t1\tMTHFR\tmissense_variant\tMODERATE\t60.0\t42\tq10;lowDP\tAB=0.5;DP=42;FLAG\t5",
            row,
        )

    def test_empty_run_exports_header_or_empty_body(self):
        self.sample()
        run = self.start_run({"genes": ["NOSUCH"]})
        _, jsonl = self.service.run_export(run["id"], "jsonl")
        _, tsv = self.service.run_export(run["id"], "tsv")
        self.assertEqual("", jsonl)
        self.assertEqual(
            "chrom\tpos\tid\tref\talt\tallele_index\tgene\tconsequence\timpact\tqual\tdp\tfilter\tinfo\tline\n",
            tsv,
        )

    def test_export_errors_and_read_only_behaviour(self):
        self.sample()
        run = self.start_run()
        with self.assertRaises(ValidationError):
            self.service.run_export(run["id"], "csv")
        with self.assertRaises(NotFoundError):
            self.service.run_export("missing", "jsonl")
        self.service.run_export(run["id"], "tsv")
        self.assertTrue(self.service.run_provenance(run["id"])["verified"])

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


if __name__ == "__main__":
    unittest.main()

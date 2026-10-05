import copy
import hashlib
import json
import unittest

from variantrail import verify_lineage
from variantrail.lineage import GENESIS, STEP_FIELDS, build_lineage, json_snapshot
from variantrail.pipeline import execute, normalize_params
from variantrail.provenance import canonical, digest_text
from variantrail.vcf import parse_vcf

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
PREAMBLE = "##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1>"

VCF = "\n".join(
    [
        PREAMBLE,
        HEADER,
        "chr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42",
        "chr1\t11856378\trs2\tG\tT\t20\tPASS\tDP=42",
        "chr7\t55019017\t.\tG\tGA\t50\tPASS\tDP=31",
        "chr7\t55019017\t.\tG\tT\t50\tPASS\t.",
        "chr17\t43093445\t.\tC\tT\t99\tPASS\tDP=55",
    ]
) + "\n"


def run(vcf=VCF, raw_params=None, lineage=True):
    parsed = parse_vcf(vcf)
    params = normalize_params(raw_params)
    return execute(parsed, vcf, params, lineage=lineage)


def sha256_body(record):
    body = {field: record[field] for field in STEP_FIELDS}
    return hashlib.sha256(canonical(body).encode("utf-8")).hexdigest()


class LineageManifestTests(unittest.TestCase):
    def test_disabled_run_keeps_the_legacy_result_shape(self):
        legacy = run(lineage=False)
        default = run(lineage=False)
        explicit = run(lineage=False)
        self.assertEqual(legacy, default)
        self.assertEqual(legacy, explicit)
        self.assertEqual({"filter_counts", "provenance", "statistics", "variants"}, set(legacy))
        self.assertNotIn("lineage", legacy)

    def test_enabled_run_adds_the_manifest_and_nothing_else(self):
        result = run()
        self.assertEqual(
            {"filter_counts", "lineage", "provenance", "statistics", "variants"}, set(result)
        )
        lineage = result["lineage"]
        self.assertEqual({"final_sha256", "input", "params", "steps"}, set(lineage))

    def test_input_summary_is_the_content_digest_of_the_vcf(self):
        lineage = run()["lineage"]
        self.assertEqual({"sha256": digest_text(VCF)}, lineage["input"])
        self.assertEqual(lineage["input"]["sha256"], lineage["steps"][0]["input_sha256"])

    def test_params_are_a_complete_call_time_snapshot(self):
        raw = {
            "pass_only": True,
            "min_qual": 30,
            "min_dp": 10,
            "genes": ["MTHFR", "EGFR"],
            "impacts": ["HIGH"],
        }
        lineage = run(raw_params=raw)["lineage"]
        self.assertEqual(
            {
                "genes": ["EGFR", "MTHFR"],
                "impacts": ["HIGH"],
                "min_dp": 10,
                "min_qual": 30.0,
                "pass_only": True,
            },
            lineage["params"],
        )

    def test_snapshot_is_isolated_from_later_mutation(self):
        params = normalize_params({"genes": ["EGFR"], "min_dp": 5})
        result = execute(parse_vcf(VCF), VCF, params, lineage=True)
        params["genes"].append("TP53")
        params["min_dp"] = 999
        params["record_filter"] = {"field": "chrom", "op": "eq", "value": "chr1"}
        params["newly_added"] = True
        manifest = result["lineage"]
        self.assertEqual(["EGFR"], manifest["params"]["genes"])
        self.assertEqual(5, manifest["params"]["min_dp"])
        self.assertNotIn("record_filter", manifest["params"])
        self.assertNotIn("newly_added", manifest["params"])
        self.assertTrue(verify_lineage(copy.deepcopy(manifest)))

    def test_equivalent_values_with_different_dict_key_order_hash_equal(self):
        ordered = build_lineage(
            VCF,
            {"a": 1, "nested": {"x": [1, 2], "y": True}},
            [
                {
                    "name": "ingest",
                    "params": {"b": 2, "a": {"z": 1, "a": 2}},
                    "input_sha256": digest_text(VCF),
                    "output_sha256": "a" * 64,
                }
            ],
        )
        reordered = build_lineage(
            VCF,
            {"nested": {"y": True, "x": [1, 2]}, "a": 1},
            [
                {
                    "name": "ingest",
                    "params": {"a": {"a": 2, "z": 1}, "b": 2},
                    "input_sha256": digest_text(VCF),
                    "output_sha256": "a" * 64,
                }
            ],
        )
        self.assertEqual(ordered, reordered)

    def test_steps_follow_real_execution_order_and_chain_the_genesis(self):
        lineage = run()["lineage"]
        records = lineage["steps"]
        self.assertEqual(["ingest", "filter", "annotate", "summarize"], [r["name"] for r in records])
        self.assertEqual(GENESIS, "0" * 64)
        self.assertEqual(GENESIS, records[0]["previous_sha256"])
        for earlier, later in zip(records, records[1:]):
            self.assertEqual(earlier["sha256"], later["previous_sha256"])
            self.assertEqual(earlier["output_sha256"], later["input_sha256"])
        self.assertEqual(records[-1]["sha256"], lineage["final_sha256"])

    def test_each_record_digest_covers_every_field_but_itself(self):
        records = run()["lineage"]["steps"]
        for record in records:
            self.assertEqual(set(STEP_FIELDS) | {"sha256"}, set(record))
            self.assertEqual(sha256_body(record), record["sha256"])
            self.assertRegex(record["sha256"], r"^[0-9a-f]{64}$")

    def test_record_digest_preserves_json_types_and_array_order(self):
        link = {
            "name": "s",
            "params": {"i": 1, "f": 1.5, "b": True, "n": None, "arr": [3, 1, 2]},
            "input_sha256": "a" * 64,
            "output_sha256": "b" * 64,
        }
        manifest = build_lineage(VCF, {}, [link])
        body = {field: manifest["steps"][0][field] for field in STEP_FIELDS}
        expected_text = json.dumps(body, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
        self.assertEqual(
            hashlib.sha256(expected_text.encode("utf-8")).hexdigest(),
            manifest["steps"][0]["sha256"],
        )
        self.assertEqual([3, 1, 2], manifest["steps"][0]["params"]["arr"])

    def test_filtered_out_everything_still_records_every_executed_step(self):
        result = run(raw_params={"genes": ["NO_SUCH_GENE"]})
        self.assertEqual([], result["variants"])
        records = result["lineage"]["steps"]
        self.assertEqual(["ingest", "filter", "annotate", "summarize"], [r["name"] for r in records])
        self.assertTrue(verify_lineage(copy.deepcopy(result["lineage"])))

    def test_skipped_steps_are_not_reported_as_executed_records(self):
        links = [
            {
                "name": "ingest",
                "params": {"parser_version": "x"},
                "input_sha256": digest_text(VCF),
                "output_sha256": "c" * 64,
            }
        ]
        manifest = build_lineage(VCF, {"min_dp": 0}, links)
        self.assertEqual(["ingest"], [r["name"] for r in manifest["steps"]])
        self.assertEqual(manifest["steps"][-1]["sha256"], manifest["final_sha256"])

    def test_repeated_runs_produce_byte_identical_manifests(self):
        first = run(raw_params={"min_qual": 30, "genes": ["EGFR", "MTHFR"]})["lineage"]
        second = run(raw_params={"min_qual": 30, "genes": ["MTHFR", "EGFR"]})["lineage"]
        self.assertEqual(first, second)

    def test_changing_an_input_changes_the_chain_from_the_first_record(self):
        first = run()["lineage"]
        other_vcf = VCF.replace("DP=42", "DP=43")
        second = run(vcf=other_vcf)["lineage"]
        self.assertNotEqual(first["input"], second["input"])
        self.assertNotEqual(first["steps"][0]["sha256"], second["steps"][0]["sha256"])
        for a, b in zip(first["steps"], second["steps"]):
            self.assertNotEqual(a["sha256"], b["sha256"])
        self.assertNotEqual(first["final_sha256"], second["final_sha256"])

    def test_changing_a_parameter_changes_hashes_from_the_first_affected_record(self):
        first = run(raw_params={"min_dp": 10})["lineage"]
        second = run(raw_params={"min_dp": 11})["lineage"]
        # Ingest runs before filtering and does not see min_dp.
        self.assertEqual(first["steps"][0]["sha256"], second["steps"][0]["sha256"])
        self.assertNotEqual(first["steps"][1]["sha256"], second["steps"][1]["sha256"])
        for a, b in zip(first["steps"][2:], second["steps"][2:]):
            self.assertNotEqual(a["sha256"], b["sha256"])
        self.assertNotEqual(first["final_sha256"], second["final_sha256"])

    def test_changing_gene_allowlist_first_moves_the_annotate_record(self):
        first = run(raw_params={"genes": ["EGFR"]})["lineage"]
        second = run(raw_params={"genes": ["BRCA1"]})["lineage"]
        self.assertEqual(first["steps"][0]["sha256"], second["steps"][0]["sha256"])
        self.assertEqual(first["steps"][1]["sha256"], second["steps"][1]["sha256"])
        self.assertNotEqual(first["steps"][2]["sha256"], second["steps"][2]["sha256"])

    def test_record_filter_is_part_of_the_snapshot_and_chain(self):
        expression = {"all": [{"field": "pos", "op": "gte", "value": 55019017}]}
        lineage = run(raw_params={"record_filter": expression})["lineage"]
        self.assertEqual(expression, lineage["params"]["record_filter"])
        self.assertEqual(
            lineage["params"]["record_filter"],
            lineage["steps"][1]["params"]["record_filter"],
        )
        self.assertTrue(verify_lineage(copy.deepcopy(lineage)))

    def test_non_boolean_lineage_argument_raises_type_error(self):
        parsed = parse_vcf(VCF)
        params = normalize_params({})
        for bad in (1, 0, "true", None, [], object()):
            with self.subTest(bad=bad):
                with self.assertRaises(TypeError):
                    execute(parsed, VCF, params, lineage=bad)

    def test_unjsonable_parameter_fails_before_any_step(self):
        parsed = parse_vcf(VCF)
        params = normalize_params({})
        params["unrepresentable"] = {1, 2, 3}
        with self.assertRaises(TypeError):
            execute(parsed, VCF, params, lineage=True)
        for token in (float("nan"), float("inf")):
            params = normalize_params({})
            params["unrepresentable"] = token
            with self.assertRaises(TypeError):
                execute(parsed, VCF, params, lineage=True)
        # The same value is accepted when lineage stays disabled.
        params = normalize_params({})
        params["unrepresentable"] = {1, 2, 3}
        result = execute(parsed, VCF, params, lineage=False)
        self.assertNotIn("lineage", result)

    def test_json_snapshot_deep_copies_containers(self):
        original = {"genes": ["EGFR"], "nested": {"keep": [True, None]}}
        snapshot = json_snapshot(original)
        snapshot["genes"].append("TP53")
        snapshot["nested"]["keep"].append(False)
        self.assertEqual(["EGFR"], original["genes"])
        self.assertEqual([True, None], original["nested"]["keep"])


class VerifyLineageTests(unittest.TestCase):
    def setUp(self):
        self.manifest = run()["lineage"]

    def verify(self, manifest=None):
        return verify_lineage(copy.deepcopy(manifest if manifest is not None else self.manifest))

    def test_intact_manifest_verifies(self):
        self.assertTrue(verify_lineage(self.manifest))

    def test_verification_does_not_mutate_the_argument(self):
        import json as json_module

        before = json_module.loads(json_module.dumps(self.manifest))
        verify_lineage(self.manifest)
        self.assertEqual(before, self.manifest)

    def test_tampering_with_any_step_field_is_detected(self):
        for index in range(len(self.manifest["steps"])):
            for field, replacement in (
                ("name", "other"),
                ("input_sha256", "f" * 64),
                ("output_sha256", "f" * 64),
                ("previous_sha256", "f" * 64),
            ):
                with self.subTest(index=index, field=field):
                    damaged = copy.deepcopy(self.manifest)
                    damaged["steps"][index][field] = replacement
                    self.assertFalse(verify_lineage(damaged))

    def test_tampering_with_step_params_is_detected(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"][1]["params"]["min_dp"] += 1
        self.assertFalse(verify_lineage(damaged))

    def test_tampering_with_top_level_input_or_params_is_detected(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["input"]["sha256"] = "f" * 64
        self.assertFalse(verify_lineage(damaged))
        damaged = copy.deepcopy(self.manifest)
        damaged["params"]["min_dp"] += 1
        self.assertFalse(verify_lineage(damaged))
        damaged = copy.deepcopy(self.manifest)
        damaged["params"]["extra"] = True
        self.assertFalse(verify_lineage(damaged))

    def test_record_reordering_is_detected(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"] = list(reversed(damaged["steps"]))
        self.assertFalse(verify_lineage(damaged))
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"][0], damaged["steps"][1] = damaged["steps"][1], damaged["steps"][0]
        self.assertFalse(verify_lineage(damaged))

    def test_missing_records_are_detected(self):
        for cut in (1, 2, 3):
            damaged = copy.deepcopy(self.manifest)
            damaged["steps"] = damaged["steps"][:cut]
            self.assertFalse(verify_lineage(damaged))

    def test_wrong_chain_head_value_is_detected(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"][0]["previous_sha256"] = "1" + "0" * 63
        self.assertFalse(verify_lineage(damaged))

    def test_final_digest_mismatch_is_detected(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["final_sha256"] = "f" * 64
        self.assertFalse(verify_lineage(damaged))

    def test_extra_top_level_or_record_fields_fail_instead_of_verifying(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["extra"] = {}
        self.assertFalse(verify_lineage(damaged))
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"][0]["extra"] = 1
        self.assertFalse(verify_lineage(damaged))
        damaged = copy.deepcopy(self.manifest)
        damaged["input"]["length"] = 10
        self.assertFalse(verify_lineage(damaged))

    def test_non_object_raises_value_error(self):
        for bad in (None, 1, 3.2, True, "manifest", [], (), object()):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    verify_lineage(bad)

    def test_missing_required_fields_raise_value_error(self):
        for key in ("input", "params", "steps", "final_sha256"):
            damaged = copy.deepcopy(self.manifest)
            del damaged[key]
            with self.subTest(key=key):
                with self.assertRaises(ValueError):
                    verify_lineage(damaged)
        damaged = copy.deepcopy(self.manifest)
        del damaged["steps"][0]["sha256"]
        with self.assertRaises(ValueError):
            verify_lineage(damaged)
        damaged = copy.deepcopy(self.manifest)
        del damaged["input"]["sha256"]
        with self.assertRaises(ValueError):
            verify_lineage(damaged)

    def test_malformed_digests_raise_value_error(self):
        good = self.manifest["final_sha256"]
        bad_digests = [
            good.upper(),          # uppercase hex
            good[:-1],             # 63 chars
            good + "a",            # 65 chars
            "g" * 64,              # non-hex
            "",                    # empty
            0,                     # wrong type
            None,
        ]
        paths = [
            lambda m, v: m.__setitem__("final_sha256", v),
            lambda m, v: m["input"].__setitem__("sha256", v),
            lambda m, v: m["steps"][0].__setitem__("sha256", v),
            lambda m, v: m["steps"][0].__setitem__("input_sha256", v),
            lambda m, v: m["steps"][0].__setitem__("output_sha256", v),
            lambda m, v: m["steps"][0].__setitem__("previous_sha256", v),
        ]
        for bad in bad_digests:
            for mutate in paths:
                damaged = copy.deepcopy(self.manifest)
                mutate(damaged, bad)
                with self.subTest(bad=bad):
                    with self.assertRaises(ValueError):
                        verify_lineage(damaged)

    def test_step_records_must_be_objects(self):
        damaged = copy.deepcopy(self.manifest)
        damaged["steps"][0] = "ingest"
        with self.assertRaises(ValueError):
            verify_lineage(damaged)


if __name__ == "__main__":
    unittest.main()

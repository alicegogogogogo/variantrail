import copy
import hashlib
import json
import unittest

from variantrail import verify_lineage
from variantrail.lineage import GENESIS_SHA256, canonical
from variantrail.pipeline import execute, normalize_params
from variantrail.vcf import parse_vcf

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
PREAMBLE = "##fileformat=VCFv4.2\n##reference=GRCh38"

VCF = "\n".join(
    [
        PREAMBLE,
        HEADER,
        "chr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42",
        "chr1\t11856379\t.\tA\tG\t20\tPASS\tDP=5",
        "chr7\t55019017\t.\tG\tT\t50\tPASS\tDP=31",
        "chr17\t43093445\t.\tC\tT\t99\tPASS\tDP=55",
    ]
) + "\n"

STEP_NAMES = ["ingest", "filter", "annotate", "summarize"]


def run(vcf=VCF, raw=None, lineage=True, parsed=None):
    params = normalize_params(raw)
    document = parse_vcf(vcf) if parsed is None else parsed
    return execute(document, vcf, params, lineage=lineage), params


def expected_record_digest(entry):
    body = {field: entry[field] for field in ("name", "params", "input_sha256", "output_sha256", "previous_sha256")}
    return hashlib.sha256(canonical(body).encode("utf-8")).hexdigest()


def expected_final_digest(lineage_obj):
    body = {field: lineage_obj[field] for field in ("input_sha256", "params", "steps")}
    return hashlib.sha256(canonical(body).encode("utf-8")).hexdigest()


class LineageManifestTests(unittest.TestCase):
    def test_disabled_run_returns_exactly_the_existing_shape(self):
        result, _ = run(lineage=False)
        self.assertEqual({"filter_counts", "provenance", "statistics", "variants"}, set(result))
        result_default, _ = run(lineage=False)
        # The keyword defaults to disabled.
        parsed = parse_vcf(VCF)
        self.assertEqual(result, execute(parsed, VCF, normalize_params(None)))
        self.assertEqual(result, result_default)

    def test_enabled_run_adds_the_lineage_object(self):
        result, params = run()
        self.assertEqual(
            {"filter_counts", "lineage", "provenance", "statistics", "variants"}, set(result)
        )
        lineage_obj = result["lineage"]
        self.assertEqual({"input_sha256", "params", "steps", "final_sha256"}, set(lineage_obj))

    def test_input_summary_is_the_sha256_of_the_input_text(self):
        result, _ = run()
        self.assertEqual(
            hashlib.sha256(VCF.encode("utf-8")).hexdigest(), result["lineage"]["input_sha256"]
        )

    def test_params_are_the_full_snapshot_of_call_time_values(self):
        result, params = run(raw={"min_qual": 30, "genes": ["EGFR", "MTHFR"]})
        self.assertEqual(params, result["lineage"]["params"])
        self.assertEqual(
            ["EGFR", "MTHFR"], result["lineage"]["params"]["genes"]
        )  # normalized, sorted

    def test_steps_follow_real_execution_order_with_exact_fields(self):
        result, _ = run()
        steps = result["lineage"]["steps"]
        self.assertEqual(STEP_NAMES, [entry["name"] for entry in steps])
        fields = {"name", "params", "input_sha256", "output_sha256", "previous_sha256", "sha256"}
        for entry in steps:
            self.assertEqual(fields, set(entry))
            for key in ("input_sha256", "output_sha256", "previous_sha256", "sha256"):
                self.assertRegex(entry[key], r"[0-9a-f]{64}")

    def test_first_record_uses_the_64_zero_head_and_chain_links_follow(self):
        result, _ = run()
        steps = result["lineage"]["steps"]
        self.assertEqual(GENESIS_SHA256, steps[0]["previous_sha256"])
        for previous, current in zip(steps, steps[1:]):
            self.assertEqual(previous["sha256"], current["previous_sha256"])

    def test_step_inputs_outputs_and_params_are_the_ones_actually_used(self):
        result, _ = run(raw={"record_filter": {"field": "chrom", "op": "eq", "value": "chr1"}})
        provenance = result["provenance"]
        for index, entry in enumerate(result["lineage"]["steps"]):
            self.assertEqual(entry["name"], provenance[index]["name"])
            self.assertEqual(entry["params"], provenance[index]["params"])
            self.assertEqual(entry["input_sha256"], provenance[index]["input_sha256"])
            self.assertEqual(entry["output_sha256"], provenance[index]["output_sha256"])
        self.assertEqual(
            {"min_dp", "min_qual", "pass_only", "record_filter"},
            set(result["lineage"]["steps"][1]["params"]),
        )
        self.assertEqual(
            {"annotation_table_sha256", "genes", "impacts"},
            set(result["lineage"]["steps"][2]["params"]),
        )
        self.assertEqual({}, result["lineage"]["steps"][3]["params"])

    def test_record_digests_cover_every_field_except_themselves(self):
        result, _ = run()
        for entry in result["lineage"]["steps"]:
            self.assertEqual(expected_record_digest(entry), entry["sha256"])

    def test_final_summary_seals_the_whole_manifest(self):
        lineage_obj = run()[0]["lineage"]
        self.assertEqual(expected_final_digest(lineage_obj), lineage_obj["final_sha256"])

    def test_manifest_verifies_and_survives_a_json_round_trip(self):
        lineage_obj = run(raw={"min_dp": 10})[0]["lineage"]
        self.assertTrue(verify_lineage(lineage_obj))
        encoded = json.dumps(lineage_obj)
        self.assertTrue(verify_lineage(json.loads(encoded)))

    def test_repeated_runs_produce_an_identical_manifest(self):
        first = run(raw={"min_qual": 30, "genes": ["EGFR"]})[0]["lineage"]
        second = run(raw={"min_qual": 30, "genes": ["EGFR"]})[0]["lineage"]
        self.assertEqual(first, second)

    def test_dictionary_key_insertion_order_does_not_change_digests(self):
        parsed = parse_vcf(VCF)
        ordered = {
            "genes": ["EGFR"], "impacts": [], "min_dp": 0, "min_qual": 0.0, "pass_only": False
        }
        reversed_order = {
            "pass_only": False, "min_qual": 0.0, "min_dp": 0, "impacts": [], "genes": ["EGFR"]
        }
        first = execute(parsed, VCF, ordered, lineage=True)["lineage"]
        second = execute(parsed, VCF, reversed_order, lineage=True)["lineage"]
        self.assertEqual(first, second)

    def test_snapshot_is_isolated_from_mutation_after_the_call(self):
        params = normalize_params({"genes": ["EGFR"], "impacts": ["HIGH"]})
        result = execute(parse_vcf(VCF), VCF, params, lineage=True)
        params["genes"].append("ZZZ")
        params["impacts"].clear()
        params["min_qual"] = 999.0
        params["new"] = {}
        self.assertEqual(["EGFR"], result["lineage"]["params"]["genes"])
        self.assertEqual(["HIGH"], result["lineage"]["params"]["impacts"])
        self.assertNotIn("new", result["lineage"]["params"])
        self.assertTrue(verify_lineage(result["lineage"]))
        # Step param containers are isolated too.
        result["lineage"]["steps"][1]["params"]["min_dp"] = 1
        self.assertEqual(0, execute(parse_vcf(VCF), VCF, normalize_params(None), lineage=True)
                         ["lineage"]["steps"][1]["params"]["min_dp"])

    def test_empty_input_still_records_every_executed_step(self):
        # execute ingests already-parsed records; an empty run is the empty
        # record list, and every executed step is still recorded.
        result = execute({"records": []}, "", normalize_params(None), lineage=True)
        self.assertEqual(STEP_NAMES, [entry["name"] for entry in result["lineage"]["steps"]])
        self.assertTrue(verify_lineage(result["lineage"]))

    def test_all_records_filtered_still_records_every_executed_step(self):
        result = run(raw={"genes": ["NO_SUCH_GENE"]})[0]
        self.assertEqual([], result["variants"])
        self.assertEqual(STEP_NAMES, [entry["name"] for entry in result["lineage"]["steps"]])
        self.assertTrue(verify_lineage(result["lineage"]))
        result = run(raw={"min_qual": 1000})[0]
        self.assertEqual(STEP_NAMES, [entry["name"] for entry in result["lineage"]["steps"]])
        self.assertTrue(verify_lineage(result["lineage"]))

    def test_changed_input_changes_digests_from_the_first_record(self):
        baseline = run()[0]["lineage"]
        # Same parsed records, different raw text: only the ingest input summary
        # changes, but every chained record must follow via previous_sha256.
        other_text = VCF.replace("##reference=GRCh38", "##reference=GRCh39")
        parsed = parse_vcf(other_text)
        altered = execute(parsed, other_text, normalize_params(None), lineage=True)["lineage"]
        self.assertNotEqual(baseline["input_sha256"], altered["input_sha256"])
        self.assertNotEqual(baseline["steps"][0]["sha256"], altered["steps"][0]["sha256"])
        for left, right in zip(baseline["steps"], altered["steps"]):
            self.assertNotEqual(left["sha256"], right["sha256"])

    def test_changed_parameter_changes_only_from_the_first_affected_record(self):
        baseline = run(raw={"min_dp": 0})[0]["lineage"]
        # min_dp 0 -> 5 changes only the filter step params; no record in the
        # fixture is affected, so outputs stay equal, but the filter record and
        # everything chained after it must differ, while ingest is untouched.
        changed = run(raw={"min_dp": 5})[0]["lineage"]
        self.assertEqual(baseline["steps"][0], changed["steps"][0])
        self.assertNotEqual(baseline["steps"][1]["sha256"], changed["steps"][1]["sha256"])
        self.assertEqual(baseline["steps"][1]["output_sha256"], changed["steps"][1]["output_sha256"])
        for index in (2, 3):
            self.assertNotEqual(baseline["steps"][index]["sha256"], changed["steps"][index]["sha256"])

    def test_annotate_only_parameter_leaves_earlier_records_identical(self):
        baseline = run(raw={"genes": []})[0]["lineage"]
        changed = run(raw={"genes": ["EGFR"]})[0]["lineage"]
        self.assertEqual(baseline["steps"][0], changed["steps"][0])
        self.assertEqual(baseline["steps"][1], changed["steps"][1])
        self.assertNotEqual(baseline["steps"][2]["sha256"], changed["steps"][2]["sha256"])
        self.assertNotEqual(baseline["steps"][3]["sha256"], changed["steps"][3]["sha256"])

    def test_changed_step_output_propagates_through_the_chain(self):
        baseline = run()[0]["lineage"]
        changed = run(raw={"min_qual": 90})[0]["lineage"]
        # Filter keeps fewer records: its output digest changes from step 2.
        self.assertEqual(baseline["steps"][0], changed["steps"][0])
        self.assertNotEqual(
            baseline["steps"][1]["output_sha256"], changed["steps"][1]["output_sha256"]
        )
        for index in (1, 2, 3):
            self.assertNotEqual(baseline["steps"][index]["sha256"], changed["steps"][index]["sha256"])


class VerifyLineageTests(unittest.TestCase):
    def setUp(self):
        self.lineage = run()[0]["lineage"]

    def clone(self, mutate=None):
        clone = copy.deepcopy(self.lineage)
        if mutate is not None:
            mutate(clone)
        return clone

    def test_accepts_an_intact_manifest_without_modifying_it(self):
        snapshot = copy.deepcopy(self.lineage)
        self.assertTrue(verify_lineage(self.lineage))
        self.assertEqual(snapshot, self.lineage)

    def test_tampered_step_params_are_detected(self):
        candidate = self.clone(lambda lin: lin["steps"][1]["params"].__setitem__("min_qual", 1))
        self.assertFalse(verify_lineage(candidate))

    def test_tampered_step_output_is_detected(self):
        candidate = self.clone(lambda lin: lin["steps"][2].__setitem__("output_sha256", GENESIS_SHA256))
        self.assertFalse(verify_lineage(candidate))

    def test_tampered_top_level_params_are_detected_by_the_final_seal(self):
        candidate = self.clone(lambda lin: lin["params"].__setitem__("min_dp", 7))
        # Record digests stay valid (params live only at the top level)...
        self.assertTrue(all(
            expected_record_digest(entry) == entry["sha256"] for entry in candidate["steps"]
        ))
        # ...but the final seal must catch it.
        self.assertFalse(verify_lineage(candidate))

    def test_tampered_input_summary_is_detected(self):
        candidate = self.clone(lambda lin: lin.__setitem__(
            "input_sha256", "a" + lin["input_sha256"][:-1]))
        self.assertFalse(verify_lineage(candidate))

    def test_tampered_final_summary_is_detected(self):
        candidate = self.clone(lambda lin: lin.__setitem__(
            "final_sha256", "f" * 64))
        self.assertFalse(verify_lineage(candidate))

    def test_reordered_records_are_detected(self):
        def reorder(lin):
            lin["steps"][0], lin["steps"][1] = lin["steps"][1], lin["steps"][0]
        self.assertFalse(verify_lineage(self.clone(reorder)))

    def test_missing_record_is_detected(self):
        self.assertFalse(verify_lineage(self.clone(lambda lin: lin["steps"].pop(1))))
        self.assertFalse(verify_lineage(self.clone(lambda lin: lin["steps"].pop())))

    def test_wrong_chain_head_is_detected(self):
        candidate = self.clone(
            lambda lin: lin["steps"][0].__setitem__("previous_sha256", "f" * 64)
        )
        self.assertFalse(verify_lineage(candidate))

    def test_extra_record_field_is_rejected(self):
        candidate = self.clone(lambda lin: lin["steps"][0].__setitem__("extra", 1))
        with self.assertRaises(ValueError):
            verify_lineage(candidate)

    def test_verification_never_mutates_its_argument(self):
        candidate = self.clone(lambda lin: lin["steps"][1]["params"].__setitem__("min_qual", 1))
        before = copy.deepcopy(candidate)
        self.assertFalse(verify_lineage(candidate))
        self.assertEqual(before, candidate)
        malformed = [1, None]
        for value in malformed:
            with self.assertRaises(ValueError):
                verify_lineage(value)

    def test_non_object_raises_value_error(self):
        for value in (None, 1, 1.5, True, "lineage", [], [self.lineage]):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    verify_lineage(value)

    def test_missing_or_unknown_top_level_fields_raise_value_error(self):
        for key in ("input_sha256", "params", "steps", "final_sha256"):
            candidate = self.clone(lambda lin, key=key: lin.pop(key))
            with self.subTest(key=key):
                with self.assertRaises(ValueError):
                    verify_lineage(candidate)
        with self.assertRaises(ValueError):
            verify_lineage(self.clone(lambda lin: lin.__setitem__("extra", 1)))

    def test_malformed_steps_raise_value_error(self):
        cases = [
            lambda lin: lin.__setitem__("steps", []),
            lambda lin: lin.__setitem__("steps", "x"),
            lambda lin: lin.__setitem__("steps", None),
            lambda lin: lin.__setitem__("steps", [1]),
            lambda lin: lin["steps"][0].pop("sha256"),
            lambda lin: lin["steps"][0].pop("name"),
            lambda lin: lin["steps"][0].__setitem__("extra", 1),
        ]
        for mutate in cases:
            candidate = self.clone(mutate)
            with self.subTest(candidate=json.dumps(candidate, default=str)[:80]):
                with self.assertRaises(ValueError):
                    verify_lineage(candidate)

    def test_wrong_but_well_shaped_step_name_breaks_the_record_digest(self):
        candidate = self.clone(lambda lin: lin["steps"][0].__setitem__("name", "extra"))
        self.assertFalse(verify_lineage(candidate))

    def test_non_string_step_name_raises_value_error(self):
        candidate = self.clone(lambda lin: lin["steps"][0].__setitem__("name", 7))
        with self.assertRaises(ValueError):
            verify_lineage(candidate)

    def test_digests_must_be_64_lowercase_hex(self):
        bad_values = [
            "a" * 63,
            "b" * 65,
            "g" * 64,
            "A" * 64,
            "0" * 63 + " ",
            None,
            7,
        ]

        def mutate_for(field, value):
            def mutate(lin):
                if field == "final_sha256":
                    lin["final_sha256"] = value
                elif field == "input_sha256":
                    lin["input_sha256"] = value
                else:
                    lin["steps"][0][field] = value
            return mutate

        for field in ("input_sha256", "final_sha256", "sha256", "previous_sha256", "output_sha256"):
            for value in bad_values:
                candidate = self.clone(mutate_for(field, value))
                with self.subTest(field=field, value=value):
                    with self.assertRaises(ValueError):
                        verify_lineage(candidate)

    def test_non_json_params_shape_raises_value_error(self):
        candidate = self.clone(lambda lin: lin["params"].__setitem__("weird", {1, 2}))
        with self.assertRaises(ValueError):
            verify_lineage(candidate)


class LineageArgumentTests(unittest.TestCase):
    def setUp(self):
        self.parsed = parse_vcf(VCF)

    def test_non_boolean_lineage_argument_raises_type_error(self):
        for value in (0, 1, "true", None, [], {}):
            with self.subTest(value=value):
                with self.assertRaises(TypeError):
                    execute(self.parsed, VCF, normalize_params(None), lineage=value)

    def test_non_json_parameter_value_raises_type_error(self):
        params = normalize_params(None)
        params["weird"] = {1, 2}
        with self.assertRaises(TypeError):
            execute(self.parsed, VCF, params, lineage=True)

        params = normalize_params(None)
        params["weird"] = object()
        with self.assertRaises(TypeError):
            execute(self.parsed, VCF, params, lineage=True)

        params = normalize_params(None)
        params["weird"] = (1, 2)
        with self.assertRaises(TypeError):
            execute(self.parsed, VCF, params, lineage=True)

        params = normalize_params(None)
        params["weird"] = float("nan")
        with self.assertRaises(TypeError):
            execute(self.parsed, VCF, params, lineage=True)

        params = normalize_params(None)
        params["weird"] = {7: 1}  # non-string key
        with self.assertRaises(TypeError):
            execute(self.parsed, VCF, params, lineage=True)

    def test_non_json_parameter_value_is_rejected_without_lineage_disabled_change(self):
        # Without lineage, execute has no JSON-representation requirement.
        params = normalize_params(None)
        params["weird"] = {1, 2}
        result = execute(self.parsed, VCF, params, lineage=False)
        self.assertNotIn("lineage", result)

    def test_non_json_parameter_failure_precedes_step_execution(self):
        params = normalize_params(None)
        params["weird"] = {1, 2}
        with self.assertRaises(TypeError):
            execute({"records": "not even iterated"}, VCF, params, lineage=True)


if __name__ == "__main__":
    unittest.main()

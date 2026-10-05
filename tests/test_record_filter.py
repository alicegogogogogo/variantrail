import unittest

from variantrail.errors import ValidationError
from variantrail.recordfilter import evaluate_record_filter, normalize_record_filter


def record(**overrides):
    base = {
        "alts": ["A", "T"],
        "chrom": "chr1",
        "dp": 42,
        "filter": ["PASS"],
        "id": "rs1",
        "info": {"DP": "42", "AF": "0.5", "FLAG": None, "WORD": "high"},
        "line": 5,
        "pos": 11856378,
        "qual": 60.0,
        "raw": "raw",
        "ref": "G",
    }
    base.update(overrides)
    return base


def matches(expression, rec):
    return evaluate_record_filter(normalize_record_filter(expression), rec)


class RecordFilterConditionTests(unittest.TestCase):
    def test_scalar_eq(self):
        self.assertTrue(matches({"field": "chrom", "op": "eq", "value": "chr1"}, record()))
        self.assertFalse(matches({"field": "chrom", "op": "eq", "value": "chr2"}, record()))
        self.assertTrue(matches({"field": "id", "op": "eq", "value": "rs1"}, record()))
        self.assertTrue(matches({"field": "ref", "op": "eq", "value": "G"}, record()))

    def test_numeric_ordering_on_pos_qual_dp(self):
        for op, value, expected in [
            ("gt", 59, True), ("gt", 60, False), ("gte", 60, True),
            ("lt", 61, True), ("lt", 60, False), ("lte", 60, True),
        ]:
            with self.subTest(op=op, value=value):
                self.assertIs(
                    expected,
                    matches({"field": "qual", "op": op, "value": value}, record()),
                )
        self.assertTrue(matches({"field": "pos", "op": "lte", "value": 11856378}, record()))
        self.assertTrue(matches({"field": "dp", "op": "gte", "value": 42}, record()))

    def test_alt_and_filter_match_any_element_for_eq_and_in(self):
        self.assertTrue(matches({"field": "alt", "op": "eq", "value": "T"}, record()))
        self.assertFalse(matches({"field": "alt", "op": "eq", "value": "C"}, record()))
        self.assertTrue(matches({"field": "alt", "op": "in", "value": ["X", "T"]}, record()))
        self.assertFalse(matches({"field": "alt", "op": "in", "value": ["X", "Y"]}, record()))
        rec = record(filter=["q10", "S50"])
        self.assertTrue(matches({"field": "filter", "op": "eq", "value": "S50"}, rec))
        self.assertTrue(matches({"field": "filter", "op": "in", "value": ["PASS", "q10"]}, rec))
        self.assertFalse(matches({"field": "filter", "op": "eq", "value": "PASS"}, rec))

    def test_in_on_scalar_field_is_membership(self):
        self.assertTrue(matches({"field": "chrom", "op": "in", "value": ["chr2", "chr1"]}, record()))
        self.assertFalse(matches({"field": "chrom", "op": "in", "value": ["chr2"]}, record()))

    def test_exists(self):
        self.assertTrue(matches({"field": "id", "op": "exists"}, record(id="rs1")))
        self.assertFalse(matches({"field": "id", "op": "exists"}, record(id=None)))
        self.assertFalse(matches({"field": "qual", "op": "exists"}, record(qual=None)))
        self.assertFalse(matches({"field": "dp", "op": "exists"}, record(dp=None)))
        self.assertTrue(matches({"field": "info.FLAG", "op": "exists"}, record()))
        self.assertTrue(matches({"field": "info.WORD", "op": "exists"}, record()))
        self.assertFalse(matches({"field": "info.ABSENT", "op": "exists"}, record()))

    def test_info_eq_and_in_compare_strings(self):
        self.assertTrue(matches({"field": "info.WORD", "op": "eq", "value": "high"}, record()))
        self.assertFalse(matches({"field": "info.WORD", "op": "eq", "value": "low"}, record()))
        self.assertTrue(matches({"field": "info.WORD", "op": "in", "value": ["low", "high"]}, record()))

    def test_info_ordering_parses_decimal_and_rejects_garbage(self):
        self.assertTrue(matches({"field": "info.AF", "op": "gt", "value": 0.4}, record()))
        self.assertTrue(matches({"field": "info.AF", "op": "eq", "value": "0.5"}, record()))
        self.assertFalse(matches({"field": "info.AF", "op": "eq", "value": "0.50"}, record()))
        self.assertFalse(matches({"field": "info.AF", "op": "gt", "value": 0.6}, record()))
        self.assertFalse(matches({"field": "info.WORD", "op": "gt", "value": 0}, record()))
        self.assertFalse(matches({"field": "info.FLAG", "op": "gt", "value": 0}, record()))
        self.assertFalse(matches({"field": "info.ABSENT", "op": "gte", "value": 0}, record()))

    def test_missing_or_null_is_false_for_other_ops_and_negatable(self):
        rec = record(id=None, qual=None)
        self.assertFalse(matches({"field": "id", "op": "eq", "value": "rs1"}, rec))
        self.assertFalse(matches({"field": "qual", "op": "lt", "value": 999}, rec))
        self.assertTrue(matches({"not": {"field": "id", "op": "exists"}}, rec))
        self.assertTrue(matches({"not": {"field": "info.ABSENT", "op": "eq", "value": "x"}}, rec))


class RecordFilterLogicalTests(unittest.TestCase):
    def test_all_any_not(self):
        t = {"field": "chrom", "op": "eq", "value": "chr1"}
        f = {"field": "chrom", "op": "eq", "value": "chr2"}
        self.assertTrue(matches({"all": [t, t]}, record()))
        self.assertFalse(matches({"all": [t, f]}, record()))
        self.assertTrue(matches({"any": [f, t]}, record()))
        self.assertFalse(matches({"any": [f, f]}, record()))
        self.assertTrue(matches({"not": f}, record()))
        self.assertFalse(matches({"not": t}, record()))

    def test_nesting(self):
        expression = {
            "all": [
                {"any": [{"field": "pos", "op": "gt", "value": 0}, {"field": "pos", "op": "lt", "value": 0}]},
                {"not": {"field": "id", "op": "exists"}},
            ]
        }
        self.assertTrue(matches(expression, record(id=None)))
        self.assertFalse(matches(expression, record(id="rs1")))


class RecordFilterValidationTests(unittest.TestCase):
    def assert_invalid(self, expression, pattern="record_filter"):
        with self.assertRaisesRegex(ValidationError, pattern):
            normalize_record_filter(expression)

    def test_unknown_field_and_op(self):
        self.assert_invalid({"field": "nope", "op": "eq", "value": "x"}, "unknown field")
        self.assert_invalid({"field": "info.", "op": "eq", "value": "x"}, "unknown field")
        self.assert_invalid({"field": "chrom", "op": "nope", "value": "x"}, "unknown op")

    def test_mixed_and_unknown_node_forms(self):
        self.assert_invalid({"all": [], "op": "eq"})
        self.assert_invalid({"all": [{"field": "chrom", "op": "eq", "value": "x"}],
                             "any": [{"field": "chrom", "op": "eq", "value": "x"}]})
        self.assert_invalid({})
        self.assert_invalid({"field": "chrom"})
        self.assert_invalid({"op": "eq"})
        self.assert_invalid({"value": "x"})
        self.assert_invalid("chrom")
        self.assert_invalid([{"field": "chrom", "op": "eq", "value": "x"}])
        self.assert_invalid({"field": "chrom", "op": "eq", "value": "x", "extra": 1})

    def test_empty_or_wrong_typed_logical_children(self):
        self.assert_invalid({"all": []})
        self.assert_invalid({"any": []})
        self.assert_invalid({"all": "x"})
        self.assert_invalid({"any": {}})
        self.assert_invalid({"not": []})
        self.assert_invalid({"not": {"all": []}})

    def test_value_shape_rules(self):
        self.assert_invalid({"field": "chrom", "op": "exists", "value": "x"}, "exists")
        self.assert_invalid({"field": "chrom", "op": "eq"})
        self.assert_invalid({"field": "chrom", "op": "in", "value": []})
        self.assert_invalid({"field": "chrom", "op": "in", "value": "x"})
        self.assert_invalid({"field": "chrom", "op": "in", "value": [""]})
        self.assert_invalid({"field": "pos", "op": "eq", "value": "x"})
        self.assert_invalid({"field": "pos", "op": "gt", "value": "10"})
        self.assert_invalid({"field": "qual", "op": "lt", "value": True})

    def test_ordering_does_not_apply_to_string_or_list_fields(self):
        for field in ("chrom", "id", "ref", "alt", "filter"):
            with self.subTest(field=field):
                self.assert_invalid({"field": field, "op": "lt", "value": 1})

    def test_depth_limit(self):
        deepest_ok = {"field": "chrom", "op": "eq", "value": "x"}
        deepest_bad = {"field": "chrom", "op": "eq", "value": "x"}
        # 16 nested not-wrappers around a condition is allowed; the 17th is not.
        for _ in range(16):
            deepest_ok = {"not": deepest_ok}
        for _ in range(17):
            deepest_bad = {"not": deepest_bad}
        normalize_record_filter(deepest_ok)
        self.assert_invalid(deepest_bad, "16 levels")


class RecordFilterNormalizationTests(unittest.TestCase):
    def test_in_arrays_are_sorted_but_numbers_stored_as_floats(self):
        normalized = normalize_record_filter({"field": "alt", "op": "in", "value": ["T", "A"]})
        self.assertEqual(["A", "T"], normalized["value"])
        numeric = normalize_record_filter({"field": "pos", "op": "in", "value": [3, 1, 2]})
        self.assertEqual([1.0, 2.0, 3.0], numeric["value"])

    def test_different_in_order_normalizes_to_the_same_expression(self):
        first = normalize_record_filter({"field": "alt", "op": "in", "value": ["A", "T"]})
        second = normalize_record_filter({"field": "alt", "op": "in", "value": ["T", "A"]})
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()

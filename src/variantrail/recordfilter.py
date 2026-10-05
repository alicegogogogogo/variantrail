from __future__ import annotations

import math
import re
from typing import Any

from .errors import ValidationError

MAX_DEPTH = 16
SCALAR_FIELDS = ("chrom", "id", "ref")
NUMERIC_FIELDS = ("pos", "qual", "dp")
LIST_FIELDS = ("alt", "filter")
OPS = ("eq", "in", "lt", "lte", "gt", "gte", "exists")
COMPARISON_OPS = ("lt", "lte", "gt", "gte")
_INFO_KEY = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*$")
_DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$")

# Record attribute holding each field; ``info.<KEY>`` is handled separately.
_FIELD_ATTRIBUTE = {
    "chrom": "chrom",
    "pos": "pos",
    "id": "id",
    "ref": "ref",
    "alt": "alts",
    "qual": "qual",
    "dp": "dp",
    "filter": "filter",
}

_MISSING = object()


def normalize_rule(node: Any, depth: int = 1) -> dict[str, Any]:
    """Validate a record_filter expression and return its normalized form.

    The normalized form keeps exactly the keys of the node's own form, so the
    stored expression is canonical and hashes deterministically.
    """
    if depth > MAX_DEPTH:
        raise ValidationError(f"params.record_filter nesting must not exceed {MAX_DEPTH} levels")
    if not isinstance(node, dict) or not node:
        raise ValidationError("params.record_filter node must be a non-empty object")

    combinators = sorted(set(node) & {"all", "any", "not"})
    if combinators:
        if len(node) != 1:
            raise ValidationError(
                "params.record_filter node must be exactly one of all, any, not or a condition object"
            )
        kind = combinators[0]
        if kind in ("all", "any"):
            children = node[kind]
            if not isinstance(children, list) or not children:
                raise ValidationError(f"params.record_filter.{kind} must be a non-empty array of expressions")
            return {kind: [normalize_rule(child, depth + 1) for child in children]}
        child = node["not"]
        if not isinstance(child, dict):
            raise ValidationError("params.record_filter.not must be a single expression object")
        return {"not": normalize_rule(child, depth + 1)}

    unknown = sorted(set(node) - {"field", "op", "value"})
    if unknown:
        raise ValidationError(f"params.record_filter condition contains unknown keys: {', '.join(unknown)}")
    if "field" not in node or "op" not in node:
        raise ValidationError("params.record_filter condition must contain field and op")

    field = node["field"]
    if not isinstance(field, str) or not _field_known(field):
        raise ValidationError(f"params.record_filter field {field!r} is not supported")
    op = node["op"]
    if not isinstance(op, str) or op not in OPS:
        raise ValidationError(f"params.record_filter op {op!r} is not supported")

    if op == "exists":
        if "value" in node:
            raise ValidationError("params.record_filter.exists does not accept a value")
        return {"field": field, "op": op}

    if "value" not in node:
        raise ValidationError(f"params.record_filter.{op} requires a value")
    value = node["value"]

    if op in COMPARISON_OPS:
        if field not in NUMERIC_FIELDS and not field.startswith("info."):
            raise ValidationError(f"params.record_filter.{op} is not valid for field {field!r}")
        if not _is_number(value):
            raise ValidationError(f"params.record_filter.{op} value must be a finite number")
        return {"field": field, "op": op, "value": value}

    if op == "eq":
        _check_scalar(field, value)
        return {"field": field, "op": op, "value": value}

    # op == "in"
    if not isinstance(value, list) or not value:
        raise ValidationError("params.record_filter.in value must be a non-empty array")
    for item in value:
        _check_scalar(field, item)
    return {"field": field, "op": op, "value": list(value)}


def rule_matches(node: dict[str, Any], record: dict[str, Any]) -> bool:
    """Evaluate a normalized record_filter expression against one record."""
    if "all" in node:
        return all(rule_matches(child, record) for child in node["all"])
    if "any" in node:
        return any(rule_matches(child, record) for child in node["any"])
    if "not" in node:
        return not rule_matches(node["not"], record)

    field = node["field"]
    op = node["op"]
    if field.startswith("info."):
        info = record["info"]
        key = field[len("info."):]
        present = key in info
        value = info[key] if present else _MISSING
        if op == "exists":
            # A bare INFO flag (stored as null) still counts as existing.
            return present
    else:
        value = record[_FIELD_ATTRIBUTE[field]]
        if op == "exists":
            return value is not None

    if value is _MISSING or value is None:
        return False

    expected = node["value"]
    if op == "eq":
        if field in LIST_FIELDS:
            return any(element == expected for element in value)
        return bool(value == expected)
    if op == "in":
        if field in LIST_FIELDS:
            return any(element in expected for element in value)
        return value in expected

    number = _numeric_value(field, value)
    if number is None:
        return False
    if op == "lt":
        return number < expected
    if op == "lte":
        return number <= expected
    if op == "gt":
        return number > expected
    return number >= expected


def _field_known(field: str) -> bool:
    if field in _FIELD_ATTRIBUTE:
        return True
    if field.startswith("info."):
        return bool(_INFO_KEY.match(field[len("info."):]))
    return False


def _is_number(value: Any) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def _check_scalar(field: str, value: Any) -> None:
    """eq/in values must be scalars matching the field's own type."""
    if field in NUMERIC_FIELDS:
        if not _is_number(value):
            raise ValidationError(f"params.record_filter value for field {field!r} must be a finite number")
    elif not isinstance(value, str):
        raise ValidationError(f"params.record_filter value for field {field!r} must be a string")


def _numeric_value(field: str, value: Any) -> float | None:
    """Numeric reading of a present, non-null field value; None when unparseable."""
    if field in NUMERIC_FIELDS:
        return float(value)
    # info.<KEY>: the stored string must fully parse as a finite decimal.
    if not isinstance(value, str) or not _DECIMAL.match(value):
        return None
    number = float(value)
    return number if math.isfinite(number) else None

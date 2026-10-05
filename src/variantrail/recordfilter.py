"""Structured record-level filter expressions.

A node is either a logical expression (``all`` / ``any`` / ``not``) or a
condition object ``{"field", "op", "value"}``. Conditions are evaluated
against parsed VCF records; missing, null or unparseable values make every
non-``exists`` condition false, so ``not`` expresses "missing" tests.
"""

from __future__ import annotations

import math
import re
from typing import Any

from .errors import ValidationError

LOGICAL_OPS = ("all", "any", "not")
CONDITION_OPS = ("eq", "exists", "gt", "gte", "in", "lt", "lte")
ORDERED_OPS = ("gt", "gte", "lt", "lte")
MAX_DEPTH = 16

_SCALAR_FIELDS = {"chrom": "string", "id": "string", "ref": "string"}
_NUMBER_FIELDS = {"dp": "number", "pos": "number", "qual": "number"}
_LIST_FIELDS = {"alt": "list", "filter": "list"}
_DECIMAL = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def normalize_record_filter(raw: Any) -> dict[str, Any]:
    """Validate and canonicalize a record_filter expression.

    Logical arrays keep their order (``all``/``any`` semantics may rely on it
    only for short-circuiting, never for the result); ``in`` arrays are sorted
    because membership is order-independent, and numeric values become floats
    like the other numeric params.
    """
    if not isinstance(raw, dict):
        raise ValidationError("params.record_filter must be an expression object")
    return _normalize_node(raw, 0)


def _normalize_node(node: Any, depth: int) -> dict[str, Any]:
    if not isinstance(node, dict):
        raise ValidationError("record_filter: every node must be an object")
    keys = set(node)
    logical = keys.intersection(LOGICAL_OPS)
    is_condition = "field" in keys or "op" in keys
    if logical and is_condition:
        raise ValidationError("record_filter: a node must not mix all/any/not with a condition")
    if len(logical) > 1:
        raise ValidationError("record_filter: a node must use exactly one of all, any, not")
    if logical:
        (name,) = logical
        if keys != {name}:
            raise ValidationError(f"record_filter: {name} node has unknown keys")
        if depth >= MAX_DEPTH:
            raise ValidationError("record_filter: nesting must not exceed 16 levels")
        child = node[name]
        if name in ("all", "any"):
            if not isinstance(child, list) or not child:
                raise ValidationError(f"record_filter: {name} must be a non-empty array")
            return {name: [_normalize_node(part, depth + 1) for part in child]}
        if not isinstance(child, dict):
            raise ValidationError("record_filter: not must wrap a single expression node")
        return {"not": _normalize_node(child, depth + 1)}
    if not is_condition:
        raise ValidationError(
            "record_filter: every node must be one of all, any, not or a condition object"
        )
    return _normalize_condition(node, keys)


def _normalize_condition(node: dict[str, Any], keys: set[str]) -> dict[str, Any]:
    if "field" not in node or "op" not in node:
        raise ValidationError("record_filter: condition must contain field and op")
    unknown = keys - {"field", "op", "value"}
    if unknown:
        raise ValidationError(f"record_filter: condition has unknown keys: {', '.join(sorted(unknown))}")
    field = node["field"]
    op = node["op"]
    if not isinstance(field, str) or not field:
        raise ValidationError("record_filter: field must be a non-empty string")
    if not isinstance(op, str) or not op:
        raise ValidationError("record_filter: op must be a non-empty string")
    kind = _field_kind(field)
    if kind is None:
        raise ValidationError(f"record_filter: unknown field {field}")
    if op not in CONDITION_OPS:
        raise ValidationError(f"record_filter: unknown op {op}")
    if op == "exists":
        if "value" in node:
            raise ValidationError("record_filter: op exists does not take a value")
        return {"field": field, "op": "exists"}
    if "value" not in node:
        raise ValidationError(f"record_filter: op {op} requires a value")
    value = node["value"]

    if op in ORDERED_OPS:
        if kind not in ("number", "info"):
            raise ValidationError(
                f"record_filter: op {op} only applies to pos, qual, dp or info.<KEY>"
            )
        return {"field": field, "op": op, "value": _number(value, field, op)}

    if op == "eq":
        if kind == "number":
            return {"field": field, "op": "eq", "value": _number(value, field, op)}
        if kind == "list":
            return {"field": field, "op": "eq", "value": _string(value, field, op)}
        # string fields and INFO values are compared verbatim as strings.
        return {"field": field, "op": "eq", "value": _string(value, field, op)}

    # op == "in"
    if not isinstance(value, list) or not value:
        raise ValidationError(f"record_filter: {field} in requires a non-empty array")
    if kind == "number":
        members = [_number(item, field, "in") for item in value]
    else:
        members = [_string(item, field, "in") for item in value]
    return {"field": field, "op": "in", "value": sorted(members)}


def _field_kind(field: str) -> str | None:
    if field in _SCALAR_FIELDS:
        return "string"
    if field in _NUMBER_FIELDS:
        return "number"
    if field in _LIST_FIELDS:
        return "list"
    if field.startswith("info.") and len(field) > len("info."):
        return "info"
    return None


def _number(value: Any, field: str, op: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"record_filter: {field} {op} requires a finite number value")
    return float(value)


def _string(value: Any, field: str, op: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValidationError(f"record_filter: {field} {op} requires a non-empty string value")
    return value


def evaluate_record_filter(expression: dict[str, Any], record: dict[str, Any]) -> bool:
    """Evaluate a normalized expression against one parsed VCF record."""
    if "all" in expression:
        return all(evaluate_record_filter(part, record) for part in expression["all"])
    if "any" in expression:
        return any(evaluate_record_filter(part, record) for part in expression["any"])
    if "not" in expression:
        return not evaluate_record_filter(expression["not"], record)
    return _evaluate_condition(expression, record)


def _evaluate_condition(condition: dict[str, Any], record: dict[str, Any]) -> bool:
    field = condition["field"]
    op = condition["op"]
    present, value = _resolve(field, record)
    if op == "exists":
        # INFO bare flags are stored as null but count as present.
        return present if field.startswith("info.") else present and value is not None
    if not present or value is None:
        return False

    if field.startswith("info."):
        if op == "in":
            return value in condition["value"]
        if op == "eq":
            return value == condition["value"]
        number = _decimal(value)
        if number is None:
            return False
        return _compare(number, op, condition["value"])

    if isinstance(value, list):
        # alt and filter: any element matching satisfies eq and in.
        if op == "eq":
            return condition["value"] in value
        return any(item in condition["value"] for item in value)

    if op == "eq":
        return value == condition["value"]
    if op == "in":
        return value in condition["value"]
    return _compare(value, op, condition["value"])


def _resolve(field: str, record: dict[str, Any]) -> tuple[bool, Any]:
    if field.startswith("info."):
        key = field[len("info."):]
        info = record["info"]
        return key in info, info.get(key)
    if field == "alt":
        value = record["alts"]
    else:
        value = record[field]
    return value is not None, value


def _decimal(text: str | None) -> float | None:
    """Fully parse a finite decimal number; bare flags and junk return None."""
    if text is None or not _DECIMAL.match(text):
        return None
    number = float(text)
    return number if math.isfinite(number) else None


def _compare(actual: float, op: str, expected: float) -> bool:
    if op == "lt":
        return actual < expected
    if op == "lte":
        return actual <= expected
    if op == "gt":
        return actual > expected
    return actual >= expected

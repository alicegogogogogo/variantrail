"""Optional step-level lineage manifest for one pipeline run.

The lineage manifest is a self-contained, hash-chained summary of a run: it
needs no logs or stored files to explain which inputs, parameters and
intermediate steps produced a result. It holds:

- ``input_sha256``  — SHA-256 of the run input text (UTF-8);
- ``params``        — an isolated, JSON-representable snapshot of every
  parameter as it was when the run started;
- ``steps``         — one record per executed step, in real execution order;
- ``final_sha256``  — the final chain summary.

A step record is::

    {"name", "params", "input_sha256", "output_sha256", "previous_sha256",
     "sha256"}

``sha256`` is the SHA-256 of the canonical JSON of all the other fields. The
first record's ``previous_sha256`` is 64 zeroes; every later record covers the
previous one, so a changed input, parameter or step output changes every
digest from the first affected record onward. ``final_sha256`` uses the same
construction one level up: it seals every other manifest field (including the
params snapshot, the input summary and all records), so no top-level field can
be tampered with either.

Canonical JSON is ``json.dumps(value, ensure_ascii=False, separators=(",",
":"), sort_keys=True)``: object keys sorted by Unicode code point, arrays in
order, no insignificant whitespace, numbers/booleans/null in their native JSON
types.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any

GENESIS_SHA256 = "0" * 64

_RECORD_FIELDS = ("input_sha256", "name", "output_sha256", "params", "previous_sha256")
_MANIFEST_FIELDS = ("input_sha256", "params", "steps", "final_sha256")
_HEX64 = re.compile(r"[0-9a-f]{64}\Z")


def canonical(value: Any) -> str:
    """Stable JSON rendering: sorted keys, no insignificant whitespace."""
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def digest(value: Any) -> str:
    """SHA-256 of the canonical JSON rendering of ``value``."""
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def digest_text(text: str) -> str:
    """SHA-256 of raw UTF-8 text (the run input, summarized verbatim)."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def json_snapshot(value: Any) -> Any:
    """Deep-copy ``value`` and prove every leaf is representable as JSON.

    Raises ``TypeError`` for values JSON cannot represent (sets, tuples,
    custom objects, non-string keys, or NaN/Infinity). The returned object is
    fully detached from ``value``: later mutation of the caller's objects can
    never reach the snapshot stored in the manifest.
    """
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return value
    if isinstance(value, float):
        if value != value or value == float("inf") or value == float("-inf"):
            raise TypeError("lineage params must contain only JSON-representable values")
        return value
    if isinstance(value, dict):
        result: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError("lineage params object keys must be strings")
            result[key] = json_snapshot(item)
        return result
    if isinstance(value, list):
        return [json_snapshot(item) for item in value]
    raise TypeError("lineage params must contain only JSON-representable values")


def record(
    name: str,
    params: Any,
    input_sha256: str,
    output_sha256: str,
    previous_sha256: str,
) -> dict[str, Any]:
    """Build one lineage record, hashing every field except its own digest."""
    body: dict[str, Any] = {
        "name": name,
        "params": params,
        "input_sha256": input_sha256,
        "output_sha256": output_sha256,
        "previous_sha256": previous_sha256,
    }
    return {**body, "sha256": digest(body)}


def manifest(vcf_text: str, params_snapshot: Any, steps: list[dict[str, Any]]) -> dict[str, Any]:
    """Assemble a manifest around executed ``steps`` and seal it.

    Only steps that actually executed belong in ``steps``; a step skipped by a
    condition must never be recorded as executed.
    """
    if not steps:
        raise ValueError("lineage requires at least one executed step")
    document: dict[str, Any] = {
        "input_sha256": digest_text(vcf_text),
        "params": params_snapshot,
        "steps": steps,
    }
    return {**document, "final_sha256": digest(document)}


def _is_hex64(value: Any) -> bool:
    return isinstance(value, str) and bool(_HEX64.fullmatch(value))


def _check_json_shape(value: Any, where: str) -> None:
    """Validate a JSON round-trippable shape; malformed input is ValueError."""
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return
    if isinstance(value, float):
        if value != value or value == float("inf") or value == float("-inf"):
            raise ValueError(f"{where} must contain only JSON-representable values")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"{where} must use string object keys")
            _check_json_shape(item, where)
        return
    if isinstance(value, list):
        for item in value:
            _check_json_shape(item, where)
        return
    raise ValueError(f"{where} must contain only JSON-representable values")


def verify_lineage(lineage_obj: Any) -> bool:
    """Recompute a lineage manifest chain and seal; never mutate the input.

    Returns ``True`` only when every record recomputes, links to the previous
    one, starts at the 64-zero head, and the final seal matches. A tampered
    field, reordered or missing record, wrong chain head, or mismatched
    ``final_sha256`` returns ``False``.

    Malformed input — not an object, missing or unknown fields, a non-string
    step name, or a digest that is not 64 lowercase hexadecimal characters —
    raises ``ValueError``.
    """
    if not isinstance(lineage_obj, dict):
        raise ValueError("lineage must be an object")
    if set(lineage_obj) != set(_MANIFEST_FIELDS):
        raise ValueError("lineage object has missing or unknown fields")
    if not _is_hex64(lineage_obj["input_sha256"]):
        raise ValueError("lineage input_sha256 must be 64 lowercase hexadecimal characters")
    if not _is_hex64(lineage_obj["final_sha256"]):
        raise ValueError("lineage final_sha256 must be 64 lowercase hexadecimal characters")

    steps = lineage_obj["steps"]
    if not isinstance(steps, list) or not steps:
        raise ValueError("lineage steps must be a non-empty array")

    previous = GENESIS_SHA256
    for index, entry in enumerate(steps):
        if not isinstance(entry, dict) or set(entry) != set(_RECORD_FIELDS) | {"sha256"}:
            raise ValueError(f"lineage step {index} has missing or unknown fields")
        if not isinstance(entry["name"], str) or not entry["name"]:
            raise ValueError(f"lineage step {index} name must be a non-empty string")
        for field in ("input_sha256", "output_sha256", "previous_sha256", "sha256"):
            if not _is_hex64(entry[field]):
                raise ValueError(f"lineage step {index} {field} must be 64 lowercase hexadecimal characters")
        _check_json_shape(entry["params"], f"lineage step {index} params")
        if entry["previous_sha256"] != previous:
            return False
        body = {field: entry[field] for field in _RECORD_FIELDS}
        if digest(body) != entry["sha256"]:
            return False
        previous = entry["sha256"]

    _check_json_shape(lineage_obj["params"], "lineage params")
    sealed = {
        "input_sha256": lineage_obj["input_sha256"],
        "params": lineage_obj["params"],
        "steps": steps,
    }
    return digest(sealed) == lineage_obj["final_sha256"]

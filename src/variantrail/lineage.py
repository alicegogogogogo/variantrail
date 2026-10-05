"""Optional step-level lineage manifest for a single pipeline run.

When the pipeline run entry point is invoked with ``lineage=True`` it builds a
manifest that lets a caller decide which inputs, parameters and intermediate
steps produced a result without reading logs or stored files:

- ``input`` — the content summary of the input consumed by the run;
- ``params`` — a complete JSON snapshot of every parameter as it was at call
  time, isolated from the caller's parameter object;
- ``steps`` — one record per step that actually executed, in real execution
  order;
- ``final_sha256`` — the digest closing the chain.

Each step record holds its stable ``name``, the parameters the step actually
used, the digests of its input and output, the digest of the previous record
(``"0" * 64`` for the first record) and its own digest. A record digest is the
SHA-256 of the record's canonical JSON over every field except the digest
itself: UTF-8 encoded, object keys sorted by Unicode code point, no
insignificant whitespace, array order preserved, JSON-native number, boolean
and null types.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any

from .provenance import canonical, digest_text

#: Previous-record digest of the first executed step.
GENESIS = "0" * 64

_HEX64 = re.compile(r"[0-9a-f]{64}")

#: Fields of a step record that participate in its digest (every field but the
#: record's own digest), in a stable order.
STEP_FIELDS = ("input_sha256", "name", "output_sha256", "params", "previous_sha256")

_MANIFEST_FIELDS = {"final_sha256", "input", "params", "steps"}
_STEP_RECORD_FIELDS = set(STEP_FIELDS) | {"sha256"}
_INPUT_FIELDS = {"sha256"}


def json_snapshot(value: Any) -> Any:
    """Return a deep copy of ``value`` isolated through a JSON round trip.

    Raises ``TypeError`` for values that cannot be represented as JSON. Taking
    the snapshot at call time means later mutation of the caller's objects can
    never change the returned manifest.
    """
    try:
        return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))
    except (TypeError, ValueError) as error:
        raise TypeError(f"lineage value is not JSON-representable: {error}") from error


def input_summary(vcf_text: str) -> dict[str, str]:
    """Content summary of the VCF text consumed by the run."""
    return {"sha256": digest_text(vcf_text)}


def record_digest(body: dict[str, Any]) -> str:
    """SHA-256 of a step record body over its digest-participating fields."""
    ordered = {field: body[field] for field in STEP_FIELDS}
    return hashlib.sha256(canonical(ordered).encode("utf-8")).hexdigest()


def build_lineage(
    vcf_text: str,
    params_snapshot: Any,
    links: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build a manifest from call-time values and the executed-step links.

    ``links`` is the provenance trail of the steps that really ran, in execution
    order; a step skipped by a conditional simply has no link and therefore
    never appears as an executed record.
    """
    records: list[dict[str, Any]] = []
    previous = GENESIS
    for link in links:
        body = {
            "input_sha256": link["input_sha256"],
            "name": link["name"],
            "output_sha256": link["output_sha256"],
            "params": json_snapshot(link["params"]),
            "previous_sha256": previous,
        }
        record = {**body, "sha256": record_digest(body)}
        records.append(record)
        previous = record["sha256"]
    return {
        "input": input_summary(vcf_text),
        "params": params_snapshot,
        "steps": records,
        "final_sha256": previous,
    }


def _is_hex64(value: Any) -> bool:
    return isinstance(value, str) and bool(_HEX64.fullmatch(value))


def _steps_by_name(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {record["name"]: record for record in manifest["steps"] if isinstance(record, dict)}


def _params_matches_steps(manifest: dict[str, Any]) -> bool:
    """The full parameter snapshot must equal the params the steps used.

    This binds the manifest-level snapshot to the hash-chained step records:
    the filter record carries the record-level parameters and the annotate
    record carries the allele allowlists.
    """
    records = _steps_by_name(manifest)
    filter_record = records.get("filter")
    annotate_record = records.get("annotate")
    if not isinstance(filter_record, dict) or not isinstance(annotate_record, dict):
        return False
    filter_params = filter_record["params"]
    annotate_params = annotate_record["params"]
    if not isinstance(filter_params, dict) or not isinstance(annotate_params, dict):
        return False
    expected = {
        "genes": annotate_params.get("genes"),
        "impacts": annotate_params.get("impacts"),
        "min_dp": filter_params.get("min_dp"),
        "min_qual": filter_params.get("min_qual"),
        "pass_only": filter_params.get("pass_only"),
    }
    if "record_filter" in filter_params:
        expected["record_filter"] = filter_params["record_filter"]
    try:
        return canonical(expected) == canonical(manifest["params"])
    except (TypeError, ValueError):
        return False


def verify_lineage(manifest: Any) -> bool:
    """Recompute a lineage manifest's chain record by record.

    Returns ``True`` only when the whole manifest is intact. A non-object
    argument, a missing required field, or a digest that is not 64 lowercase
    hexadecimal characters raises ``ValueError``; any tampered field, reordered
    or missing record, wrong chain-head value or mismatching ``final_sha256``
    returns ``False``. The argument is never modified.
    """
    if not isinstance(manifest, dict):
        raise ValueError("lineage manifest must be an object")
    if not _MANIFEST_FIELDS <= set(manifest):
        raise ValueError("lineage manifest is missing required fields")
    if set(manifest) != _MANIFEST_FIELDS:
        return False

    summary = manifest["input"]
    if not isinstance(summary, dict):
        raise ValueError("lineage input summary must be an object")
    if "sha256" not in summary:
        raise ValueError("lineage input summary is missing the sha256 field")
    if not _is_hex64(summary["sha256"]):
        raise ValueError("lineage input sha256 must be 64 lowercase hexadecimal characters")
    if set(summary) != _INPUT_FIELDS:
        return False
    if not _is_hex64(manifest["final_sha256"]):
        raise ValueError("lineage final_sha256 must be 64 lowercase hexadecimal characters")

    steps = manifest["steps"]
    if not isinstance(steps, list):
        raise ValueError("lineage steps must be an array")
    if not steps:
        # A real run always executes steps; an empty trail means records were
        # removed, which is chain inconsistency rather than a format error.
        return False

    for index, record in enumerate(steps):
        if not isinstance(record, dict):
            raise ValueError(f"lineage step {index} record must be an object")
        if not _STEP_RECORD_FIELDS <= set(record):
            raise ValueError(f"lineage step {index} record is missing required fields")
        if set(record) != _STEP_RECORD_FIELDS:
            return False
        for field in ("input_sha256", "output_sha256", "previous_sha256", "sha256"):
            if not _is_hex64(record[field]):
                raise ValueError(f"lineage step {index} {field} must be 64 lowercase hexadecimal characters")

    previous = GENESIS
    for record in steps:
        if record["previous_sha256"] != previous:
            return False
        body = {field: record[field] for field in STEP_FIELDS}
        try:
            recomputed = record_digest(body)
        except (TypeError, ValueError):
            return False
        if recomputed != record["sha256"]:
            return False
        previous = record["sha256"]

    if steps[0]["input_sha256"] != summary["sha256"]:
        return False
    if not _params_matches_steps(manifest):
        return False
    return manifest["final_sha256"] == previous

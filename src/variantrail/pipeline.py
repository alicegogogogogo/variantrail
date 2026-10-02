from __future__ import annotations

import math
from typing import Any

from .annotate import annotate, table_sha256
from .errors import ValidationError
from .provenance import digest, digest_text, step

PARSER_VERSION = "variantrail-vcf-1"
ALLOWED_IMPACTS = ("HIGH", "LOW", "MODERATE", "MODIFIER", "UNKNOWN")
RECORD_FILTER_REASONS = ("dp_below_min", "not_pass", "qual_below_min")


def normalize_params(raw: Any) -> dict[str, Any]:
    """Validate run parameters and fill every documented default."""
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValidationError("params must be an object")
    unknown = sorted(set(raw) - {"genes", "impacts", "min_dp", "min_qual", "pass_only"})
    if unknown:
        raise ValidationError(f"params contains unknown fields: {', '.join(unknown)}")

    params: dict[str, Any] = {"genes": [], "impacts": [], "min_dp": 0, "min_qual": 0.0, "pass_only": False}
    if "min_qual" in raw:
        value = raw["min_qual"]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValidationError("params.min_qual must be a non-negative finite number")
        params["min_qual"] = float(value)
    if "min_dp" in raw:
        value = raw["min_dp"]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValidationError("params.min_dp must be a non-negative integer")
        params["min_dp"] = value
    if "pass_only" in raw:
        if not isinstance(raw["pass_only"], bool):
            raise ValidationError("params.pass_only must be a boolean")
        params["pass_only"] = raw["pass_only"]
    for key in ("genes", "impacts"):
        if key not in raw:
            continue
        value = raw[key]
        if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
            raise ValidationError(f"params.{key} must be an array of non-empty strings")
        if len(set(value)) != len(value):
            raise ValidationError(f"params.{key} must not contain duplicate values")
        params[key] = sorted(value)
    for impact in params["impacts"]:
        if impact not in ALLOWED_IMPACTS:
            raise ValidationError(f"params.impacts entry {impact} is not one of {', '.join(ALLOWED_IMPACTS)}")
    return params


def execute(parsed: dict[str, Any], vcf_text: str, params: dict[str, Any]) -> dict[str, Any]:
    """Run ingest -> filter -> annotate -> summarize and record each link."""
    steps: list[dict[str, Any]] = []

    records = parsed["records"]
    steps.append(
        step(
            1,
            "ingest",
            {"parser_version": PARSER_VERSION},
            digest_text(vcf_text),
            digest({"records": records}),
            None,
        )
    )

    counts = dict.fromkeys(RECORD_FILTER_REASONS, 0)
    kept: list[dict[str, Any]] = []
    for record in records:
        reason = record_reason(record, params)
        if reason is None:
            kept.append(record)
        else:
            counts[reason] += 1
    steps.append(
        step(
            2,
            "filter",
            {"min_dp": params["min_dp"], "min_qual": params["min_qual"], "pass_only": params["pass_only"]},
            steps[-1]["output_sha256"],
            digest({"filter_counts": counts, "records": kept}),
            steps[-1]["hash"],
        )
    )

    variants: list[dict[str, Any]] = []
    filtered = 0
    for record in kept:
        annotations = [entry for entry in annotate(record) if allele_selected(entry, params)]
        if not annotations:
            filtered += 1
            continue
        variants.append(variant_document(record, annotations))
    counts = {**counts, "allele_filtered": filtered}
    steps.append(
        step(
            3,
            "annotate",
            {
                "annotation_table_sha256": table_sha256(),
                "genes": params["genes"],
                "impacts": params["impacts"],
            },
            steps[-1]["output_sha256"],
            digest({"filter_counts": counts, "variants": variants}),
            steps[-1]["hash"],
        )
    )

    statistics = summarize(records, variants, counts)
    steps.append(step(4, "summarize", {}, steps[-1]["output_sha256"], digest(statistics), steps[-1]["hash"]))

    return {"filter_counts": counts, "provenance": steps, "statistics": statistics, "variants": variants}


def record_reason(record: dict[str, Any], params: dict[str, Any]) -> str | None:
    """First failing record-level filter, in the documented precedence order."""
    if params["pass_only"] and record["filter"] != ["PASS"]:
        return "not_pass"
    if record["qual"] is None:
        if params["min_qual"] > 0:
            return "qual_below_min"
    elif record["qual"] < params["min_qual"]:
        return "qual_below_min"
    if record["dp"] is None:
        if params["min_dp"] > 0:
            return "dp_below_min"
    elif record["dp"] < params["min_dp"]:
        return "dp_below_min"
    return None


def allele_selected(annotation: dict[str, Any], params: dict[str, Any]) -> bool:
    """Allele-level allowlists, applied after annotation."""
    if params["genes"] and annotation["gene"] not in params["genes"]:
        return False
    if params["impacts"] and annotation["impact"] not in params["impacts"]:
        return False
    return True


def variant_document(record: dict[str, Any], annotations: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "alts": list(record["alts"]),
        "annotations": annotations,
        "chrom": record["chrom"],
        "dp": record["dp"],
        "filter": list(record["filter"]),
        "id": record["id"],
        "info": dict(record["info"]),
        "line": record["line"],
        "pos": record["pos"],
        "qual": record["qual"],
        "raw": record["raw"],
        "ref": record["ref"],
    }


def summarize(
    records: list[dict[str, Any]],
    variants: list[dict[str, Any]],
    counts: dict[str, int],
) -> dict[str, Any]:
    consequence_counts: dict[str, int] = {}
    gene_counts: dict[str, int] = {}
    impact_counts: dict[str, int] = {}
    alleles_kept = 0
    for variant in variants:
        for annotation in variant["annotations"]:
            alleles_kept += 1
            _bump(consequence_counts, annotation["consequence"])
            _bump(gene_counts, annotation["gene"] or "NA")
            _bump(impact_counts, annotation["impact"])
    return {
        "alleles_kept": alleles_kept,
        "alleles_total": sum(len(record["alts"]) for record in records),
        "consequence_counts": _ordered(consequence_counts),
        "filter_counts": _ordered(counts),
        "gene_counts": _ordered(gene_counts),
        "impact_counts": _ordered(impact_counts),
        "records_filtered": len(records) - len(variants),
        "records_kept": len(variants),
        "records_total": len(records),
    }


def _bump(counts: dict[str, int], key: str) -> None:
    counts[key] = counts.get(key, 0) + 1


def _ordered(counts: dict[str, int]) -> dict[str, int]:
    return {key: counts[key] for key in sorted(counts)}

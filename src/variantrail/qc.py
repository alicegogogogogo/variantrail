from __future__ import annotations

from typing import Any

_ACGT = frozenset("ACGT")
_TRANSITIONS = (frozenset("AG"), frozenset("CT"))


def run_qc(document: dict[str, Any]) -> dict[str, Any]:
    """Summarize the variants already stored on a finished run.

    Only stored data is read: no filter step or annotation is re-executed and
    nothing is rewritten. The result is deterministic so repeated requests for
    the same run serialize byte for byte.
    """
    variants = document["variants"]

    total = len(variants)
    pass_count = 0
    missing_qual = 0
    missing_dp = 0
    multiallelic = 0
    qual_values: list[float] = []
    dp_values: list[int] = []
    for variant in variants:
        if variant["filter"] == ["PASS"]:
            pass_count += 1
        if variant["qual"] is None:
            missing_qual += 1
        else:
            qual_values.append(variant["qual"])
        if variant["dp"] is None:
            missing_dp += 1
        else:
            dp_values.append(variant["dp"])
        if len(variant["annotations"]) >= 2:
            multiallelic += 1

    alleles_total = 0
    classes = {"deletion": 0, "insertion": 0, "mnv": 0, "non_acgt_snv": 0, "snv": 0}
    transitions = 0
    transversions = 0
    for variant in variants:
        ref = variant["ref"]
        ref_length = len(ref)
        for annotation in variant["annotations"]:
            alleles_total += 1
            allele = annotation["allele"]
            alt_length = len(allele)
            if ref_length == 1 and alt_length == 1:
                # Each allele lands in exactly one class: an SNV carrying N on
                # either end is non_acgt_snv, never counted as snv.
                if ref in _ACGT and allele in _ACGT:
                    classes["snv"] += 1
                    if frozenset((ref, allele)) in _TRANSITIONS:
                        transitions += 1
                    else:
                        transversions += 1
                else:
                    classes["non_acgt_snv"] += 1
            elif ref_length == alt_length:
                classes["mnv"] += 1
            elif alt_length > ref_length:
                classes["insertion"] += 1
            else:
                classes["deletion"] += 1

    return {
        "run_id": document["id"],
        "sample_id": document["sample_id"],
        "sample_sha256": document["sample_sha256"],
        "provenance_head": document["provenance"][-1]["hash"],
        "records": {
            "total": total,
            "pass": pass_count,
            "missing_qual": missing_qual,
            "missing_dp": missing_dp,
            "multiallelic": multiallelic,
            "qual": _bounds(qual_values),
            "dp": _bounds(dp_values),
        },
        "alleles": {
            "total": alleles_total,
            "snv": classes["snv"],
            "mnv": classes["mnv"],
            "insertion": classes["insertion"],
            "deletion": classes["deletion"],
            "non_acgt_snv": classes["non_acgt_snv"],
            "transitions": transitions,
            "transversions": transversions,
            "ti_tv_ratio": None if transversions == 0 else transitions / transversions,
        },
    }


def _bounds(values: list[Any]) -> dict[str, Any]:
    """Count and min/max the non-null values; nulls were already excluded."""
    if not values:
        return {"count": 0, "min": None, "max": None}
    return {"count": len(values), "min": min(values), "max": max(values)}

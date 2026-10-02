from __future__ import annotations

from typing import Any, Callable

# Buckets an allele of the union can fall into, in response key order.
BUCKETS = ("shared", "left_only", "right_only")
_MISSING_GENE = "NA"


def run_alleles(document: dict[str, Any]) -> dict[tuple[Any, ...], dict[str, Any]]:
    """Index the retained ALT alleles of one run by (chrom, pos, ref, alt).

    An identity that occurs more than once is counted only once.
    """
    alleles: dict[tuple[Any, ...], dict[str, Any]] = {}
    for variant in document["variants"]:
        for annotation in variant["annotations"]:
            identity = (variant["chrom"], variant["pos"], variant["ref"], annotation["allele"])
            if identity in alleles:
                continue
            alleles[identity] = {
                "chrom": variant["chrom"],
                "pos": variant["pos"],
                "ref": variant["ref"],
                "alt": annotation["allele"],
                "gene": annotation["gene"],
                "consequence": annotation["consequence"],
                "impact": annotation["impact"],
            }
    return alleles


def compare_documents(left_document: dict[str, Any], right_document: dict[str, Any]) -> dict[str, Any]:
    """Compare the final retained alleles of two completed, read-only runs."""
    left_alleles = run_alleles(left_document)
    right_alleles = run_alleles(right_document)
    left_ids = set(left_alleles)
    right_ids = set(right_alleles)

    shared = _sorted(left_alleles[i] for i in left_ids & right_ids)
    left_only = _sorted(left_alleles[i] for i in left_ids - right_ids)
    right_only = _sorted(right_alleles[i] for i in right_ids - left_ids)

    return {
        "left": _side(left_document, len(left_alleles)),
        "right": _side(right_document, len(right_alleles)),
        "counts": {
            "shared": len(shared),
            "left_only": len(left_only),
            "right_only": len(right_only),
            "union": len(shared) + len(left_only) + len(right_only),
        },
        "shared": shared,
        "left_only": left_only,
        "right_only": right_only,
        "gene_summary": _group_summary(shared, left_only, right_only, _gene_key),
        "impact_summary": _group_summary(shared, left_only, right_only, lambda allele: allele["impact"]),
    }


def _side(document: dict[str, Any], allele_count: int) -> dict[str, Any]:
    return {
        "run_id": document["id"],
        "sample_id": document["sample_id"],
        "sample_sha256": document["sample_sha256"],
        "allele_count": allele_count,
    }


def _sorted(alleles: Any) -> list[dict[str, Any]]:
    return sorted(alleles, key=lambda allele: (allele["chrom"], allele["pos"], allele["ref"], allele["alt"]))


def _gene_key(allele: dict[str, Any]) -> str:
    return allele["gene"] or _MISSING_GENE


def _group_summary(
    shared: list[dict[str, Any]],
    left_only: list[dict[str, Any]],
    right_only: list[dict[str, Any]],
    key: Callable[[dict[str, Any]], str],
) -> dict[str, dict[str, int]]:
    """Count alleles per group in each of the three buckets, keyed lexicographically."""
    summary: dict[str, dict[str, int]] = {}
    for name, alleles in zip(BUCKETS, (shared, left_only, right_only)):
        for allele in alleles:
            group = key(allele)
            entry = summary.setdefault(group, {bucket: 0 for bucket in BUCKETS})
            entry[name] += 1
    return {group: summary[group] for group in sorted(summary)}

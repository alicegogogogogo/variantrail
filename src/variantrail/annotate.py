from __future__ import annotations

from typing import Any

from .provenance import digest

# Small synthetic demonstration table: (chrom, pos, ref, alt, gene, consequence, impact).
ANNOTATION_TABLE: tuple[tuple[str, int, str, str, str, str, str], ...] = (
    ("chr1", 11856378, "G", "A", "MTHFR", "missense_variant", "MODERATE"),
    ("chr1", 11856378, "G", "T", "MTHFR", "stop_gained", "HIGH"),
    ("chr7", 55019017, "G", "GA", "EGFR", "frameshift_variant", "HIGH"),
    ("chr7", 55019017, "G", "T", "EGFR", "stop_gained", "HIGH"),
    ("chr7", 55273000, "C", "T", "EGFR", "synonymous_variant", "LOW"),
    ("chr11", 5227002, "T", "C", "HBB", "missense_variant", "MODERATE"),
    ("chr12", 25245350, "C", "A", "KRAS", "missense_variant", "MODERATE"),
    ("chr12", 25245350, "C", "T", "KRAS", "missense_variant", "MODERATE"),
    ("chr12", 25245351, "C", "T", "KRAS", "splice_donor_variant", "HIGH"),
    ("chr17", 43093445, "C", "T", "BRCA1", "stop_gained", "HIGH"),
    ("chr17", 43093446, "A", "G", "BRCA1", "missense_variant", "MODERATE"),
    ("chr17", 7676154, "C", "T", "TP53", "missense_variant", "MODERATE"),
    ("chr19", 11200138, "C", "T", "LDLR", "missense_variant", "MODERATE"),
    ("chrX", 31137345, "A", "G", "DMD", "missense_variant", "MODERATE"),
)

UNANNOTATED_CONSEQUENCE = "intergenic_variant"
UNANNOTATED_IMPACT = "MODIFIER"

_LOOKUP = {
    (chrom, pos, ref, alt): (gene, consequence, impact)
    for chrom, pos, ref, alt, gene, consequence, impact in ANNOTATION_TABLE
}


def table_document() -> list[dict[str, Any]]:
    """The annotation table as a deterministic list of records."""
    return [
        {
            "alt": alt,
            "chrom": chrom,
            "consequence": consequence,
            "gene": gene,
            "impact": impact,
            "pos": pos,
            "ref": ref,
        }
        for chrom, pos, ref, alt, gene, consequence, impact in ANNOTATION_TABLE
    ]


def table_sha256() -> str:
    """Content hash of the built-in annotation table."""
    return digest(table_document())


def annotate(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Annotate every ALT allele of one parsed VCF record."""
    annotations: list[dict[str, Any]] = []
    for index, allele in enumerate(record["alts"], start=1):
        entry = _LOOKUP.get((record["chrom"], record["pos"], record["ref"], allele))
        if entry is None:
            gene, consequence, impact = None, UNANNOTATED_CONSEQUENCE, UNANNOTATED_IMPACT
        else:
            gene, consequence, impact = entry
        annotations.append(
            {
                "allele": allele,
                "allele_index": index,
                "consequence": consequence,
                "gene": gene,
                "impact": impact,
            }
        )
    return annotations

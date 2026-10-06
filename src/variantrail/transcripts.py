from __future__ import annotations

from typing import Any

from .annotate import ANNOTATION_TABLE

# Extra transcript candidates beyond the one MANE_SELECT entry per annotation
# table record: (chrom, pos, ref, alt, transcript_id, status, consequence,
# impact). Like the annotation table this is small, synthetic and built in.
EXTRA_TRANSCRIPTS: tuple[tuple[str, int, str, str, str, str, str, str], ...] = (
    ("chr1", 11856378, "G", "A", "MTHFR-002", "CANONICAL", "synonymous_variant", "LOW"),
    ("chr1", 11856378, "G", "A", "MTHFR-003", "OTHER", "intron_variant", "MODIFIER"),
    ("chr7", 55019017, "G", "GA", "EGFR-002", "CANONICAL", "inframe_insertion", "MODERATE"),
)

_STATUS_ORDER = {"MANE_SELECT": 0, "CANONICAL": 1, "OTHER": 2}
_IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "UNKNOWN": 4}


def _candidate(transcript_id: str, status: str, consequence: str, impact: str) -> dict[str, Any]:
    return {
        "consequence": consequence,
        "impact": impact,
        "status": status,
        "transcript_id": transcript_id,
    }


def _build() -> dict[tuple[str, int, str, str], tuple[str, list[dict[str, Any]]]]:
    """The transcript lookup: allele identity -> (gene, candidate list)."""
    table: dict[tuple[str, int, str, str], tuple[str, list[dict[str, Any]]]] = {}
    for chrom, pos, ref, alt, gene, consequence, impact in ANNOTATION_TABLE:
        key = (chrom, pos, ref, alt)
        gene_candidates = table.setdefault(key, (gene, []))
        gene_candidates[1].append(_candidate(f"{gene}-001", "MANE_SELECT", consequence, impact))
    for chrom, pos, ref, alt, transcript_id, status, consequence, impact in EXTRA_TRANSCRIPTS:
        key = (chrom, pos, ref, alt)
        gene_candidates = table.setdefault(key, (transcript_id.rsplit("-", 1)[0], []))
        gene_candidates[1].append(_candidate(transcript_id, status, consequence, impact))
    return table


_LOOKUP = _build()


def transcripts_for(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, list[dict[str, Any]]] | None:
    """``(gene, sorted candidates)`` for one exact allele, or ``None``.

    Candidates sort by status (MANE_SELECT, CANONICAL, OTHER), then impact
    (HIGH, MODERATE, LOW, MODIFIER, UNKNOWN), then transcript id text order,
    so the first entry is always the preferred consequence.
    """
    entry = _LOOKUP.get((chrom, pos, ref, alt))
    if entry is None:
        return None
    gene, candidates = entry
    ordered = sorted(
        candidates,
        key=lambda candidate: (
            _STATUS_ORDER[candidate["status"]],
            _IMPACT_ORDER[candidate["impact"]],
            candidate["transcript_id"],
        ),
    )
    return gene, [dict(candidate) for candidate in ordered]

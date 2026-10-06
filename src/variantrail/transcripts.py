from __future__ import annotations

from typing import Any

from .annotate import ANNOTATION_TABLE

# Extra candidates beyond the per-record MANE_SELECT transcript:
# (chrom, pos, ref, alt, transcript_id, status, consequence, impact).
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


def _build() -> tuple[dict[tuple[str, int, str, str], str], dict[tuple[str, int, str, str], list[dict[str, Any]]]]:
    genes: dict[tuple[str, int, str, str], str] = {}
    table: dict[tuple[str, int, str, str], list[dict[str, Any]]] = {}
    for chrom, pos, ref, alt, gene, consequence, impact in ANNOTATION_TABLE:
        key = (chrom, pos, ref, alt)
        genes[key] = gene
        table.setdefault(key, []).append(_candidate(f"{gene}-001", "MANE_SELECT", consequence, impact))
    for chrom, pos, ref, alt, transcript_id, status, consequence, impact in EXTRA_TRANSCRIPTS:
        table[(chrom, pos, ref, alt)].append(_candidate(transcript_id, status, consequence, impact))
    for candidates in table.values():
        candidates.sort(
            key=lambda candidate: (
                _STATUS_ORDER[candidate["status"]],
                _IMPACT_ORDER[candidate["impact"]],
                candidate["transcript_id"],
            )
        )
    return genes, table


_GENES, _TRANSCRIPTS = _build()


def lookup_transcripts(chrom: str, pos: int, ref: str, alt: str) -> tuple[str | None, list[dict[str, Any]]]:
    """Gene and sorted candidate transcripts for one exact allele.

    Candidates are ordered by status (MANE_SELECT, CANONICAL, OTHER), then by
    impact (HIGH, MODERATE, LOW, MODIFIER, UNKNOWN), then by transcript id.
    With no exact candidate the gene is ``None`` and the list is empty.
    """
    key = (chrom, pos, ref, alt)
    return _GENES.get(key), [dict(candidate) for candidate in _TRANSCRIPTS.get(key, ())]

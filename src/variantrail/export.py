from __future__ import annotations

import json
from typing import Any

from .provenance import canonical

TSV_HEADER = (
    "chrom\tpos\tid\tref\talt\tallele_index\tgene\tconsequence\timpact\tqual\tdp\tfilter\tinfo\tline"
)
_MISSING = "NA"
_EMPTY = "."


def render_jsonl(variants: list[dict[str, Any]]) -> bytes:
    """Render one canonical variant object per line, LF terminated, UTF-8."""
    return "".join(f"{canonical(variant)}\n" for variant in variants).encode("utf-8")


def render_tsv(variants: list[dict[str, Any]]) -> bytes:
    """Render the fixed header plus one row per retained annotation.

    Multi-ALT variants are split into one row per annotation, keeping variant
    file order and annotation order.
    """
    lines = [TSV_HEADER]
    for variant in variants:
        for annotation in variant["annotations"]:
            lines.append(_row(variant, annotation))
    return ("\n".join(lines) + "\n").encode("utf-8")


def _row(variant: dict[str, Any], annotation: dict[str, Any]) -> str:
    return "\t".join(
        [
            variant["chrom"],
            _number(variant["pos"]),
            _na(variant["id"]),
            variant["ref"],
            annotation["allele"],
            _number(annotation["allele_index"]),
            _na(annotation["gene"]),
            annotation["consequence"],
            annotation["impact"],
            _number(variant["qual"]),
            _number(variant["dp"]),
            _filter(variant["filter"]),
            _info(variant["info"]),
            _number(variant["line"]),
        ]
    )


def _na(value: Any) -> str:
    return _MISSING if value is None else value


def _number(value: Any) -> str:
    return _MISSING if value is None else json.dumps(value)


def _filter(values: list[str]) -> str:
    return ";".join(values) if values else _EMPTY


def _info(values: dict[str, str | None]) -> str:
    if not values:
        return _EMPTY
    return ";".join(
        key if value is None else f"{key}={value}"
        for key, value in sorted(values.items())
    )

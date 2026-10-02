from __future__ import annotations

import json
from typing import Any

TSV_HEADER = "chrom\tpos\tid\tref\talt\tallele_index\tgene\tconsequence\timpact\tqual\tdp\tfilter\tinfo\tline"


def to_jsonl(variants: list[dict[str, Any]]) -> str:
    """One compact JSON object per variant, keys sorted, LF-terminated lines."""
    return "".join(
        json.dumps(variant, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
        for variant in variants
    )


def to_tsv(variants: list[dict[str, Any]]) -> str:
    """One row per kept annotation, in variant file order then annotation order."""
    lines = [TSV_HEADER]
    for variant in variants:
        for annotation in variant["annotations"]:
            lines.append(
                "\t".join(
                    [
                        variant["chrom"],
                        _text(variant["pos"]),
                        _text(variant["id"]),
                        variant["ref"],
                        annotation["allele"],
                        _text(annotation["allele_index"]),
                        _text(annotation["gene"]),
                        annotation["consequence"],
                        annotation["impact"],
                        _text(variant["qual"]),
                        _text(variant["dp"]),
                        _filter(variant["filter"]),
                        _info(variant["info"]),
                        _text(variant["line"]),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _text(value: Any) -> str:
    """Null becomes NA; numbers keep their JSON decimal representation."""
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value
    return json.dumps(value)


def _filter(values: list[str]) -> str:
    return ";".join(values) if values else "."


def _info(info: dict[str, str | None]) -> str:
    if not info:
        return "."
    return ";".join(
        key if info[key] is None else f"{key}={info[key]}" for key in sorted(info)
    )

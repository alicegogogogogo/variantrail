from __future__ import annotations

import json
from typing import Any

from .errors import ValidationError
from .provenance import canonical

TSV_HEADER = (
    "chrom\tpos\tid\tref\talt\tallele_index\tgene\tconsequence\timpact\tqual\tdp\tfilter\tinfo\tline"
)
VCF_COLUMN_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
VCF_INFO_DECLARATIONS = (
    '##INFO=<ID=ALLELE_INDEX,Number=1,Type=Integer,Description="Index of the ALT allele in the original record">',
    '##INFO=<ID=VT_GENE,Number=1,Type=String,Description="VariantRail gene annotation, NA when no gene is annotated">',
    '##INFO=<ID=VT_CONSEQUENCE,Number=1,Type=String,Description="VariantRail consequence annotation">',
    '##INFO=<ID=VT_IMPACT,Number=1,Type=String,Description="VariantRail impact annotation">',
)
VCF_RESERVED_INFO_KEYS = ("ALLELE_INDEX", "VT_CONSEQUENCE", "VT_GENE", "VT_IMPACT")
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


def render_vcf(variants: list[dict[str, Any]], meta: list[dict[str, str]]) -> bytes:
    """Render one VCF v4.2 row per retained ALT allele, LF terminated, UTF-8.

    The header opens with a fixed fileformat line, replays the sample's other
    meta lines in upload order, then adds the source and INFO declarations.
    Original INFO keys come first (sorted), followed by the four export keys.
    """
    lines = ["##fileformat=VCFv4.2"]
    lines.extend(entry["raw"] for entry in meta if entry["key"] != "fileformat")
    lines.append("##source=variantrail")
    lines.extend(VCF_INFO_DECLARATIONS)
    lines.append(VCF_COLUMN_HEADER)
    for variant in variants:
        _reject_reserved_info(variant["info"])
        for annotation in variant["annotations"]:
            lines.append(_vcf_row(variant, annotation))
    return ("\n".join(lines) + "\n").encode("utf-8")


def _reject_reserved_info(info: dict[str, str | None]) -> None:
    conflicts = sorted(set(info) & set(VCF_RESERVED_INFO_KEYS))
    if conflicts:
        raise ValidationError(f"variant INFO contains reserved export keys: {', '.join(conflicts)}")


def _vcf_row(variant: dict[str, Any], annotation: dict[str, Any]) -> str:
    info = _info(variant["info"])
    export_info = ";".join(
        [
            f"ALLELE_INDEX={annotation['allele_index']}",
            f"VT_GENE={_na(annotation['gene'])}",
            f"VT_CONSEQUENCE={annotation['consequence']}",
            f"VT_IMPACT={annotation['impact']}",
        ]
    )
    return "\t".join(
        [
            variant["chrom"],
            str(variant["pos"]),
            _dot(variant["id"]),
            variant["ref"],
            annotation["allele"],
            _dot_number(variant["qual"]),
            _filter(variant["filter"]),
            export_info if info == _EMPTY else f"{info};{export_info}",
        ]
    )


def _dot(value: Any) -> str:
    return _EMPTY if value is None else value


def _dot_number(value: Any) -> str:
    return _EMPTY if value is None else json.dumps(value)


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

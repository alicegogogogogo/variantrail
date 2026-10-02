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
    '##INFO=<ID=ALLELE_INDEX,Number=1,Type=Integer,Description="1-based index of the ALT allele in the source record">',
    '##INFO=<ID=VT_GENE,Number=1,Type=String,Description="Annotated gene symbol, NA when unannotated">',
    '##INFO=<ID=VT_CONSEQUENCE,Number=1,Type=String,Description="Annotated sequence consequence">',
    '##INFO=<ID=VT_IMPACT,Number=1,Type=String,Description="Annotated impact">',
)
_RESERVED_INFO_KEYS = ("ALLELE_INDEX", "VT_CONSEQUENCE", "VT_GENE", "VT_IMPACT")
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
    """Render the retained variants as VCF v4.2 text, UTF-8 and LF terminated.

    The header opens with a fresh ##fileformat line, keeps the sample's other
    meta lines in upload order, and declares the four export INFO keys. Each
    retained ALT allele is one record, in variant file order and allele_index
    order. A run with no variants yields just the header.
    """
    lines = ["##fileformat=VCFv4.2"]
    lines.extend(entry["raw"] for entry in meta if entry["key"] != "fileformat")
    lines.append("##source=variantrail")
    lines.extend(VCF_INFO_DECLARATIONS)
    lines.append(VCF_COLUMN_HEADER)
    for variant in variants:
        _check_reserved_info(variant)
        for annotation in sorted(variant["annotations"], key=lambda entry: entry["allele_index"]):
            lines.append(_vcf_row(variant, annotation))
    return ("\n".join(lines) + "\n").encode("utf-8")


def _check_reserved_info(variant: dict[str, Any]) -> None:
    conflicts = sorted(set(variant["info"]) & set(_RESERVED_INFO_KEYS))
    if conflicts:
        raise ValidationError(f"INFO contains reserved export keys: {', '.join(conflicts)}")


def _vcf_row(variant: dict[str, Any], annotation: dict[str, Any]) -> str:
    return "\t".join(
        [
            variant["chrom"],
            _number(variant["pos"]),
            _dot(variant["id"]),
            variant["ref"],
            annotation["allele"],
            _dot_number(variant["qual"]),
            _filter(variant["filter"]),
            _vcf_info(variant, annotation),
        ]
    )


def _vcf_info(variant: dict[str, Any], annotation: dict[str, Any]) -> str:
    entries = [
        key if value is None else f"{key}={value}"
        for key, value in sorted(variant["info"].items())
    ]
    entries.append(f"ALLELE_INDEX={annotation['allele_index']}")
    entries.append(f"VT_GENE={_na(annotation['gene'])}")
    entries.append(f"VT_CONSEQUENCE={annotation['consequence']}")
    entries.append(f"VT_IMPACT={annotation['impact']}")
    return ";".join(entries)


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

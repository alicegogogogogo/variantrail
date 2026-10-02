from __future__ import annotations

import re
from typing import Any

from .errors import ValidationError

_FILEFORMAT = re.compile(r"VCFv4\.[0-9]+$")
_META_KEY = re.compile(r"[A-Za-z][A-Za-z0-9_.-]*$")
_CHROM = re.compile(r"[A-Za-z0-9_.]+$")
_SEQUENCE = re.compile(r"[ACGTN]+$")
_QUAL = re.compile(r"[0-9]+(?:\.[0-9]+)?$")
_FILTER_TOKEN = re.compile(r"[A-Za-z0-9_.]+$")
_INFO_KEY = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*$")
_HEADER = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


def parse_vcf(text: Any) -> dict[str, Any]:
    """Parse the supported VCF subset.

    Every structural or field-level violation raises ``ValidationError`` whose
    message names the offending line. Record order is preserved; ``raw`` keeps
    the original line so anything in the output stays traceable.
    """
    if not isinstance(text, str):
        raise ValidationError("vcf must be a string")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not any(line.strip() for line in lines):
        raise ValidationError("vcf must not be empty")

    meta: list[dict[str, str]] = []
    records: list[dict[str, Any]] = []
    chromosomes: list[str] = []
    seen_chromosomes: set[str] = set()
    seen_variants: dict[tuple[Any, ...], int] = {}
    header_line: int | None = None
    current_chromosome: str | None = None
    last_position = 0

    for number, source in enumerate(lines, start=1):
        line = source[:-1] if source.endswith("\r") else source
        if not line.strip():
            raise ValidationError(f"line {number}: blank lines are not allowed")

        if header_line is None:
            if line.startswith("##"):
                entry = _meta_entry(line, number)
                if not meta and entry["key"] != "fileformat":
                    raise ValidationError(f"line {number}: the first line must be a ##fileformat meta line")
                if entry["key"] == "fileformat" and not _FILEFORMAT.match(entry["value"]):
                    raise ValidationError(f"line {number}: fileformat must be VCFv4.x")
                meta.append(entry)
                continue
            if line.startswith("#"):
                if tuple(line.split("\t")) != _HEADER:
                    raise ValidationError(
                        f"line {number}: the #CHROM header must be exactly "
                        "#CHROM POS ID REF ALT QUAL FILTER INFO separated by tabs"
                    )
                header_line = number
                continue
            raise ValidationError(f"line {number}: data record appears before the #CHROM header line")

        if line.startswith("#"):
            raise ValidationError(f"line {number}: header lines are not allowed after the #CHROM header")

        record = _record(line, number)
        chromosome = record["chrom"]
        if chromosome != current_chromosome:
            if chromosome in seen_chromosomes:
                raise ValidationError(f"line {number}: chromosome {chromosome} appears in more than one block")
            seen_chromosomes.add(chromosome)
            chromosomes.append(chromosome)
            current_chromosome = chromosome
            last_position = 0
        if record["pos"] < last_position:
            raise ValidationError(f"line {number}: records must be sorted by POS within {chromosome}")
        last_position = record["pos"]

        key = (chromosome, record["pos"], record["ref"], tuple(record["alts"]))
        if key in seen_variants:
            raise ValidationError(
                f"line {number}: duplicate variant {chromosome}:{record['pos']} "
                f"{record['ref']}>{','.join(record['alts'])} first appeared at line {seen_variants[key]}"
            )
        seen_variants[key] = number
        records.append(record)

    if header_line is None:
        raise ValidationError("vcf has no #CHROM header line")
    if not records:
        raise ValidationError("vcf has no data records after the #CHROM header")
    return {"chromosomes": chromosomes, "meta": meta, "records": records}


def _meta_entry(line: str, number: int) -> dict[str, str]:
    key, separator, value = line[2:].partition("=")
    if not separator or not _META_KEY.match(key) or not value:
        raise ValidationError(f"line {number}: meta lines must look like ##key=value")
    return {"key": key, "raw": line, "value": value}


def _record(line: str, number: int) -> dict[str, Any]:
    fields = line.split("\t")
    if len(fields) != 8:
        raise ValidationError(f"line {number}: expected 8 tab-separated fields, found {len(fields)}")
    chrom, pos_text, id_text, ref_text, alt_text, qual_text, filter_text, info_text = fields

    if not _CHROM.match(chrom):
        raise ValidationError(f"line {number}: CHROM must be a non-empty chromosome name")
    if not pos_text.isdigit() or int(pos_text) < 1:
        raise ValidationError(f"line {number}: POS must be a positive integer")
    position = int(pos_text)

    identifier = None if id_text == "." else id_text
    if identifier is not None and (not identifier or any(character in identifier for character in " \t;,")):
        raise ValidationError(f"line {number}: ID must be '.' or a value without whitespace, ';' or ','")

    reference = ref_text.upper()
    if not _SEQUENCE.match(reference):
        raise ValidationError(f"line {number}: REF must be a non-empty sequence of A, C, G, T or N")

    alternates = [allele.upper() for allele in alt_text.split(",")]
    if any(not _SEQUENCE.match(allele) for allele in alternates):
        raise ValidationError(f"line {number}: ALT must be a comma-separated list of A, C, G, T or N sequences")
    if len(set(alternates)) != len(alternates):
        raise ValidationError(f"line {number}: ALT must not repeat an allele")

    if qual_text == ".":
        quality: float | None = None
    elif _QUAL.match(qual_text):
        quality = float(qual_text)
    else:
        raise ValidationError(f"line {number}: QUAL must be '.' or a non-negative number")

    if filter_text == ".":
        filters: list[str] = []
    else:
        filters = filter_text.split(";")
        if any(not _FILTER_TOKEN.match(token) for token in filters):
            raise ValidationError(f"line {number}: FILTER must be '.', 'PASS' or ';'-separated tokens")
        if len(set(filters)) != len(filters):
            raise ValidationError(f"line {number}: FILTER must not repeat a token")
        if "PASS" in filters and len(filters) > 1:
            raise ValidationError(f"line {number}: FILTER 'PASS' must not be combined with other tokens")

    info: dict[str, str | None] = {}
    if info_text != ".":
        for entry in info_text.split(";"):
            key, separator, value = entry.partition("=")
            if not _INFO_KEY.match(key):
                raise ValidationError(f"line {number}: INFO entry {entry!r} does not start with a valid key")
            if key in info:
                raise ValidationError(f"line {number}: INFO key {key} is repeated")
            if not separator:
                info[key] = None
            elif not value or any(character.isspace() for character in value):
                raise ValidationError(f"line {number}: INFO {key} must have a non-empty value without whitespace")
            else:
                info[key] = value

    if "DP" not in info:
        depth: int | None = None
    elif info["DP"] is None or not info["DP"].isdigit():
        raise ValidationError(f"line {number}: INFO/DP must be a non-negative integer")
    else:
        depth = int(info["DP"])

    return {
        "alts": alternates,
        "chrom": chrom,
        "dp": depth,
        "filter": filters,
        "id": identifier,
        "info": info,
        "line": number,
        "pos": position,
        "qual": quality,
        "raw": line,
        "ref": reference,
    }

from __future__ import annotations

from typing import Any, Callable

from .errors import ConflictError, NotFoundError, ValidationError
from .export import render_jsonl, render_tsv, render_vcf
from .model import identifier, run_view, sample_view
from .pipeline import execute, normalize_params
from .provenance import digest_text, verify_chain
from .store import Store
from .vcf import parse_vcf


class VariantRail:
    def __init__(self, database: str):
        self.store = Store(database)

    def _idempotent(self, key: str | None, operation: str, action: Callable[[], dict[str, Any]]) -> dict[str, Any]:
        if not key:
            raise ValidationError("Idempotency-Key header is required")
        with self.store.transaction() as connection:
            existing = connection.execute("SELECT operation, response FROM idempotency WHERE key = ?", (key,)).fetchone()
            if existing:
                if existing["operation"] != operation:
                    raise ConflictError("idempotency key was already used for another operation")
                return self.store.decode(existing["response"])
            response = action()
            connection.execute(
                "INSERT INTO idempotency(key, operation, response) VALUES (?, ?, ?)",
                (key, operation, self.store.encode(response)),
            )
            return response

    def create_sample(self, raw: Any, key: str | None) -> dict[str, Any]:
        if not isinstance(raw, dict) or set(raw) != {"id", "vcf"}:
            raise ValidationError("sample must contain exactly id and vcf")
        sample_id = identifier(raw["id"], "sample id")
        vcf_text = raw["vcf"]
        if not isinstance(vcf_text, str):
            raise ValidationError("vcf must be a string")
        parsed = parse_vcf(vcf_text)
        document = {
            "chromosomes": parsed["chromosomes"],
            "id": sample_id,
            "meta": parsed["meta"],
            "records": parsed["records"],
            "sha256": digest_text(vcf_text),
            "vcf": vcf_text,
        }

        def create() -> dict[str, Any]:
            try:
                self.store.connection.execute(
                    "INSERT INTO samples(id, document) VALUES (?, ?)",
                    (sample_id, self.store.encode(document)),
                )
            except Exception as error:
                if "UNIQUE constraint" in str(error):
                    raise ConflictError(f"sample {sample_id} already exists") from error
                raise
            return sample_view(document)

        return self._idempotent(key, f"create-sample:{sample_id}", create)

    def get_sample(self, sample_id: str) -> dict[str, Any]:
        row = self.store.connection.execute("SELECT document FROM samples WHERE id = ?", (sample_id,)).fetchone()
        if not row:
            raise NotFoundError(f"sample {sample_id} was not found")
        return self.store.decode(row["document"])

    def create_run(self, sample_id: str, raw: Any, key: str | None) -> dict[str, Any]:
        if not isinstance(raw, dict) or "id" not in raw or set(raw) - {"id", "params"}:
            raise ValidationError("run must contain exactly id and optional params")
        run_id = identifier(raw["id"], "run id")
        params = normalize_params(raw.get("params"))

        def create() -> dict[str, Any]:
            sample = self.get_sample(sample_id)
            result = execute({"records": sample["records"]}, sample["vcf"], params)
            document = {
                "id": run_id,
                "params": params,
                "provenance": result["provenance"],
                "sample_id": sample_id,
                "sample_sha256": sample["sha256"],
                "statistics": result["statistics"],
                "variants": result["variants"],
            }
            try:
                self.store.connection.execute(
                    "INSERT INTO runs(id, sample_id, document) VALUES (?, ?, ?)",
                    (run_id, sample_id, self.store.encode(document)),
                )
            except Exception as error:
                if "UNIQUE constraint" in str(error):
                    raise ConflictError(f"run {run_id} already exists") from error
                raise
            return run_view(document)

        return self._idempotent(key, f"create-run:{run_id}", create)

    def get_run(self, run_id: str) -> dict[str, Any]:
        return run_view(self._document(run_id))

    def run_variants(self, run_id: str) -> dict[str, Any]:
        variants = self._document(run_id)["variants"]
        return {"count": len(variants), "run_id": run_id, "variants": variants}

    def run_export(self, run_id: str, format: str) -> tuple[str, bytes]:
        """Render the retained variants of a saved run; the run is unchanged."""
        if format == "jsonl":
            content_type = "application/x-ndjson; charset=utf-8"
            variants = self._document(run_id)["variants"]
            body = render_jsonl(variants)
        elif format == "tsv":
            content_type = "text/tab-separated-values; charset=utf-8"
            variants = self._document(run_id)["variants"]
            body = render_tsv(variants)
        elif format == "vcf":
            content_type = "text/x-variant-call-format; charset=utf-8"
            document = self._document(run_id)
            sample = self.get_sample(document["sample_id"])
            body = render_vcf(document["variants"], sample["meta"])
        else:
            raise ValidationError("format must be one of jsonl, tsv, vcf")
        return content_type, body

    def compare_runs(self, run_id: str, other_run_id: str) -> dict[str, Any]:
        """Compare the retained ALT alleles of two finished runs; read-only."""
        left = self._document(run_id)
        right = self._document(other_run_id)
        left_alleles = _alleles(left)
        right_alleles = _alleles(right)
        shared_keys = left_alleles.keys() & right_alleles.keys()
        left_only_keys = left_alleles.keys() - right_alleles.keys()
        right_only_keys = right_alleles.keys() - left_alleles.keys()
        shared = _pick(left_alleles, shared_keys)
        left_only = _pick(left_alleles, left_only_keys)
        right_only = _pick(right_alleles, right_only_keys)
        buckets = (("shared", shared), ("left_only", left_only), ("right_only", right_only))
        return {
            "counts": {
                "left_only": len(left_only),
                "right_only": len(right_only),
                "shared": len(shared),
                "union": len(shared) + len(left_only) + len(right_only),
            },
            "gene_summary": _summary(buckets, lambda allele: allele["gene"] or "NA"),
            "impact_summary": _summary(buckets, lambda allele: allele["impact"]),
            "left": _side(run_id, left, left_alleles),
            "left_only": left_only,
            "right": _side(other_run_id, right, right_alleles),
            "right_only": right_only,
            "shared": shared,
        }

    def cohort_runs(self, run_id: str, others: str | None) -> dict[str, Any]:
        """Compare the retained ALT alleles of two or more runs; read-only."""
        run_ids = [run_id] + (others.split(",") if others else [])
        if len(run_ids) < 2:
            raise ValidationError("cohort requires at least two distinct runs")
        if len(set(run_ids)) != len(run_ids):
            raise ConflictError("cohort run ids must be distinct")
        documents = [self._document(identifier) for identifier in run_ids]
        run_count = len(run_ids)
        per_run = [_alleles(document) for document in documents]

        union: dict[tuple[str, int, str, str], dict[str, Any]] = {}
        for alleles in per_run:
            for key, allele in alleles.items():
                union.setdefault(key, allele)

        rows: list[dict[str, Any]] = []
        frequency: dict[int, int] = {}
        private: dict[str, int] = {identifier: 0 for identifier in run_ids}
        for key in sorted(union):
            present_in = [run_ids[index] for index, alleles in enumerate(per_run) if key in alleles]
            count = len(present_in)
            frequency[count] = frequency.get(count, 0) + 1
            if count == 1:
                private[present_in[0]] += 1
            source = union[key]
            rows.append(
                {
                    "chrom": source["chrom"],
                    "pos": source["pos"],
                    "ref": source["ref"],
                    "alt": source["alt"],
                    "gene": source["gene"],
                    "consequence": source["consequence"],
                    "impact": source["impact"],
                    "run_count": count,
                    "present_in": present_in,
                }
            )

        unique_alleles = len(rows)
        core_alleles = frequency.get(run_count, 0)
        private_alleles = frequency.get(1, 0)
        variable_alleles = unique_alleles - core_alleles - private_alleles

        def cohort_summary(key_of: Callable[[dict[str, Any]], str]) -> dict[str, Any]:
            """Per-key unique/core/variable allele counts, keys sorted."""
            counts: dict[str, dict[str, int]] = {}
            for allele in rows:
                entry = counts.setdefault(
                    key_of(allele),
                    {"unique_alleles": 0, "core_alleles": 0, "variable_alleles": 0},
                )
                entry["unique_alleles"] += 1
                if allele["run_count"] == run_count:
                    entry["core_alleles"] += 1
                elif 1 < allele["run_count"] < run_count:
                    entry["variable_alleles"] += 1
            return {key: counts[key] for key in sorted(counts)}

        return {
            "run_count": run_count,
            "runs": [_side(identifier, document, alleles) for identifier, document, alleles in zip(run_ids, documents, per_run)],
            "alleles": rows,
            "counts": {
                "unique_alleles": unique_alleles,
                "core_alleles": core_alleles,
                "variable_alleles": variable_alleles,
                "private_alleles": {identifier: private[identifier] for identifier in run_ids},
            },
            "frequency_summary": {str(number): frequency.get(number, 0) for number in range(1, run_count + 1)},
            "gene_summary": cohort_summary(lambda allele: allele["gene"] or "NA"),
            "impact_summary": cohort_summary(lambda allele: allele["impact"]),
        }

    def run_provenance(self, run_id: str) -> dict[str, Any]:
        document = self._document(run_id)
        sample = self.get_sample(document["sample_id"])
        reproduced = execute({"records": sample["records"]}, sample["vcf"], document["params"])
        chain_verified = verify_chain(document["provenance"])
        reproduction_verified = (
            [entry["output_sha256"] for entry in reproduced["provenance"]]
            == [entry["output_sha256"] for entry in document["provenance"]]
            and reproduced["statistics"] == document["statistics"]
            and reproduced["variants"] == document["variants"]
        )
        return {
            "chain_verified": chain_verified,
            "head_sha256": document["provenance"][-1]["hash"],
            "reproduction_verified": reproduction_verified,
            "run_id": run_id,
            "steps": document["provenance"],
            "verified": chain_verified and reproduction_verified,
        }

    def _document(self, run_id: str) -> dict[str, Any]:
        row = self.store.connection.execute("SELECT document FROM runs WHERE id = ?", (run_id,)).fetchone()
        if not row:
            raise NotFoundError(f"run {run_id} was not found")
        return self.store.decode(row["document"])


def _alleles(document: dict[str, Any]) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Distinct retained ALT alleles of a run, keyed by (chrom, pos, ref, alt)."""
    alleles: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    for variant in document["variants"]:
        for annotation in variant["annotations"]:
            key = (variant["chrom"], variant["pos"], variant["ref"], annotation["allele"])
            if key not in alleles:
                alleles[key] = {
                    "alt": annotation["allele"],
                    "chrom": variant["chrom"],
                    "consequence": annotation["consequence"],
                    "gene": annotation["gene"],
                    "impact": annotation["impact"],
                    "pos": variant["pos"],
                    "ref": variant["ref"],
                }
    return alleles


def _pick(alleles: dict[tuple[str, int, str, str], dict[str, Any]], keys: Any) -> list[dict[str, Any]]:
    """Alleles for the given identities, sorted by chrom, pos, ref, alt."""
    return [alleles[key] for key in sorted(keys)]


def _side(run_id: str, document: dict[str, Any], alleles: dict[Any, Any]) -> dict[str, Any]:
    return {
        "allele_count": len(alleles),
        "run_id": run_id,
        "sample_id": document["sample_id"],
        "sample_sha256": document["sample_sha256"],
    }


def _summary(buckets: Any, key_of: Callable[[dict[str, Any]], str]) -> dict[str, Any]:
    """Per-key allele counts in each comparison bucket, keys sorted."""
    counts: dict[str, dict[str, int]] = {}
    for bucket, alleles in buckets:
        for allele in alleles:
            entry = counts.setdefault(key_of(allele), {"left_only": 0, "right_only": 0, "shared": 0})
            entry[bucket] += 1
    return {key: counts[key] for key in sorted(counts)}

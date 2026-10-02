from __future__ import annotations

from typing import Any, Callable

from .errors import ConflictError, NotFoundError, ValidationError
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

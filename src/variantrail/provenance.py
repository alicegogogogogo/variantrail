from __future__ import annotations

import hashlib
import json
from typing import Any

STEP_FIELDS = ("step", "name", "params", "input_sha256", "output_sha256", "previous_sha256")


def canonical(value: Any) -> str:
    """Stable JSON rendering: sorted keys, no insignificant whitespace."""
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def digest(value: Any) -> str:
    """SHA-256 of the canonical JSON rendering of ``value``."""
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def digest_text(text: str) -> str:
    """SHA-256 of a raw UTF-8 text payload."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def step(
    number: int,
    name: str,
    params: dict[str, Any],
    input_sha256: str,
    output_sha256: str,
    previous_sha256: str | None,
) -> dict[str, Any]:
    """Build one provenance link; its hash covers the previous link's hash."""
    body = {
        "step": number,
        "name": name,
        "params": params,
        "input_sha256": input_sha256,
        "output_sha256": output_sha256,
        "previous_sha256": previous_sha256,
    }
    return {**body, "hash": digest(body)}


def verify_chain(steps: Any) -> bool:
    """Recompute every link hash and the chain linkage."""
    if not isinstance(steps, list) or not steps:
        return False
    previous: str | None = None
    for number, entry in enumerate(steps, start=1):
        if not isinstance(entry, dict) or set(entry) != set(STEP_FIELDS) | {"hash"}:
            return False
        if entry["step"] != number or entry["previous_sha256"] != previous:
            return False
        body = {field: entry[field] for field in STEP_FIELDS}
        if digest(body) != entry["hash"]:
            return False
        previous = entry["hash"]
    return True

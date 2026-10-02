from __future__ import annotations

import re
from typing import Any

from .errors import ValidationError

_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,99}$")


def identifier(value: Any, field: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.match(value):
        raise ValidationError(
            f"{field} must match [A-Za-z0-9][A-Za-z0-9._:-]{{0,99}}"
        )
    return value


def sample_view(document: dict[str, Any]) -> dict[str, Any]:
    """Public projection of a stored sample; no timestamps, fully deterministic."""
    return {
        "chromosomes": list(document["chromosomes"]),
        "id": document["id"],
        "meta_lines": len(document["meta"]),
        "records": len(document["records"]),
        "sha256": document["sha256"],
    }


def run_view(document: dict[str, Any]) -> dict[str, Any]:
    """Public projection of a stored run."""
    provenance = document["provenance"]
    return {
        "id": document["id"],
        "params": document["params"],
        "provenance_head": provenance[-1]["hash"],
        "provenance_length": len(provenance),
        "sample_id": document["sample_id"],
        "sample_sha256": document["sample_sha256"],
        "statistics": document["statistics"],
        "status": "succeeded",
    }

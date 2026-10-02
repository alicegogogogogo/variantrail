from __future__ import annotations

import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, NamedTuple
from urllib.parse import urlsplit

from .errors import NotFoundError, ValidationError, VariantRailError
from .service import VariantRail


class RawResponse(NamedTuple):
    """A pre-rendered response body with its own Content-Type (e.g. exports)."""

    content_type: str
    body: bytes


class Handler(BaseHTTPRequestHandler):
    service: VariantRail

    def log_message(self, format: str, *args: Any) -> None:
        return

    def _json(self, status: int, value: Any) -> None:
        body = json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode()
        self._respond(status, "application/json; charset=utf-8", body)

    def _respond(self, status: int, content_type: str, body: bytes) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _body(self) -> Any:
        content_type = self.headers.get("Content-Type", "")
        if content_type.split(";", 1)[0].strip().lower() != "application/json":
            raise ValidationError("Content-Type must be application/json")
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length < 0 or length > 5_000_000:
                raise ValueError
            return json.loads(self.rfile.read(length))
        except (ValueError, json.JSONDecodeError) as error:
            raise ValidationError("request body must be valid JSON") from error

    def _dispatch(self) -> tuple[int, Any]:
        parts = [part for part in urlsplit(self.path).path.split("/") if part]
        if self.command == "GET" and parts == ["health"]:
            return 200, {"status": "ok"}
        if self.command == "POST" and parts == ["samples"]:
            return 201, self.service.create_sample(self._body(), self.headers.get("Idempotency-Key"))
        if len(parts) == 3 and parts[0] == "samples" and parts[2] == "runs" and self.command == "POST":
            return 201, self.service.create_run(parts[1], self._body(), self.headers.get("Idempotency-Key"))
        if len(parts) == 2 and parts[0] == "runs" and self.command == "GET":
            return 200, self.service.get_run(parts[1])
        if len(parts) == 3 and parts[0] == "runs" and parts[2] == "variants" and self.command == "GET":
            return 200, self.service.run_variants(parts[1])
        if len(parts) == 4 and parts[0] == "runs" and parts[2] == "exports" and self.command == "GET":
            content_type, body = self.service.run_export(parts[1], parts[3])
            return 200, RawResponse(content_type, body.encode("utf-8"))
        if len(parts) == 3 and parts[0] == "runs" and parts[2] == "provenance" and self.command == "GET":
            return 200, self.service.run_provenance(parts[1])
        raise NotFoundError("route was not found")

    def _handle(self) -> None:
        try:
            status, response = self._dispatch()
            if isinstance(response, RawResponse):
                self._respond(status, response.content_type, response.body)
            else:
                self._json(status, response)
        except VariantRailError as error:
            self._json(error.status, {"error": {"code": error.code, "message": str(error)}})
        except Exception:
            self._json(500, {"error": {"code": "internal_error", "message": "internal server error"}})

    do_GET = _handle
    do_POST = _handle


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the VariantRail HTTP service")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8080, type=int)
    parser.add_argument("--database", default="variantrail.db")
    arguments = parser.parse_args()
    Handler.service = VariantRail(arguments.database)
    server = ThreadingHTTPServer((arguments.host, arguments.port), Handler)
    print(f"VariantRail listening on http://{arguments.host}:{arguments.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()

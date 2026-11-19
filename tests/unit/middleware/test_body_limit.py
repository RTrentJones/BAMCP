"""Request body size limit middleware."""

import pytest
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import PlainTextResponse
from starlette.routing import Route
from starlette.testclient import TestClient

from bamcp import __main__ as entry
from bamcp.config import BAMCPConfig
from bamcp.middleware.bodylimit import AUTH_MAX_REQUEST_BYTES, BodySizeLimitMiddleware


def _app(max_bytes: int = 1000) -> Starlette:
    async def echo(request: Request) -> PlainTextResponse:
        body = await request.body()
        return PlainTextResponse(f"len={len(body)}")

    app = Starlette(
        routes=[
            Route("/mcp", echo, methods=["POST"]),
            Route("/register", echo, methods=["POST"]),
            Route("/token", echo, methods=["POST"]),
        ]
    )
    app.add_middleware(BodySizeLimitMiddleware, max_bytes=max_bytes)
    return app


def _chunks(total: int, size: int = 100):
    sent = 0
    while sent < total:
        n = min(size, total - sent)
        sent += n
        yield b"x" * n


@pytest.mark.unit
class TestBodySizeLimit:
    def test_small_body_passes(self):
        resp = TestClient(_app()).post("/mcp", content=b"x" * 1000)
        assert resp.status_code == 200
        assert resp.text == "len=1000"

    def test_declared_content_length_over_limit_is_413(self):
        resp = TestClient(_app()).post("/mcp", content=b"x" * 1001)
        assert resp.status_code == 413

    def test_streamed_body_without_content_length_is_413(self):
        client = TestClient(_app())
        resp = client.post("/mcp", content=_chunks(5000))
        assert resp.status_code == 413

    def test_streamed_body_under_limit_passes(self):
        resp = TestClient(_app()).post("/mcp", content=_chunks(900))
        assert resp.status_code == 200
        assert resp.text == "len=900"

    def test_lying_content_length_is_cut_off(self):
        """A small declared length with a larger streamed body is still enforced."""
        sent: list[dict] = []
        app = BodySizeLimitMiddleware(_app(max_bytes=10**9), max_bytes=100)

        messages = [
            {"type": "http.request", "body": b"x" * 60, "more_body": True},
            {"type": "http.request", "body": b"x" * 60, "more_body": False},
        ]

        async def receive():
            return messages.pop(0) if messages else {"type": "http.disconnect"}

        async def send(message):
            sent.append(message)

        import anyio

        scope = {
            "type": "http",
            "method": "POST",
            "path": "/mcp",
            "raw_path": b"/mcp",
            "query_string": b"",
            "headers": [(b"content-length", b"10"), (b"host", b"testserver")],
            "server": ("testserver", 80),
            "scheme": "http",
            "root_path": "",
            "http_version": "1.1",
            "client": ("127.0.0.1", 1),
        }
        anyio.run(app, scope, receive, send)
        starts = [m for m in sent if m["type"] == "http.response.start"]
        assert [m["status"] for m in starts] == [413]

    def test_auth_paths_have_smaller_cap(self):
        client = TestClient(_app(max_bytes=10**6))
        body = b"x" * (AUTH_MAX_REQUEST_BYTES + 1)
        assert client.post("/register", content=body).status_code == 413
        assert client.post("/token", content=body).status_code == 413
        assert client.post("/register/", content=body).status_code == 413
        assert client.post("/mcp", content=body).status_code == 200
        assert client.post("/token", content=b"x" * AUTH_MAX_REQUEST_BYTES).status_code == 200

    def test_invalid_content_length_is_400(self):
        resp = TestClient(_app()).post("/mcp", content=b"x", headers={"content-length": "abc"})
        assert resp.status_code == 400

    def test_get_without_body_unaffected(self):
        async def home(request: Request) -> PlainTextResponse:
            return PlainTextResponse("ok")

        app = Starlette(routes=[Route("/", home)])
        app.add_middleware(BodySizeLimitMiddleware, max_bytes=1)
        assert TestClient(app).get("/").status_code == 200


@pytest.mark.unit
def test_wired_from_config_and_413_has_security_headers():
    async def echo(request: Request) -> PlainTextResponse:
        return PlainTextResponse(str(len(await request.body())))

    app = Starlette(routes=[Route("/mcp", echo, methods=["POST"])])
    entry._add_http_middleware(app, BAMCPConfig(max_request_bytes=50, rate_limit=1))
    client = TestClient(app)
    resp = client.post("/mcp", content=b"x" * 51)
    assert resp.status_code == 413
    assert resp.headers.get("x-content-type-options") == "nosniff"
    # The rejected request did not consume the single rate-limit slot.
    assert client.post("/mcp", content=b"x" * 50).status_code == 200


@pytest.mark.unit
def test_max_request_bytes_env(monkeypatch):
    monkeypatch.setenv("BAMCP_MAX_REQUEST_BYTES", "2048")
    assert BAMCPConfig.from_env().max_request_bytes == 2048
    monkeypatch.delenv("BAMCP_MAX_REQUEST_BYTES")
    assert BAMCPConfig.from_env().max_request_bytes == 1024 * 1024
    with pytest.raises(ValueError, match="max_request_bytes"):
        BAMCPConfig(max_request_bytes=0)

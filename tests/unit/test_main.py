"""Unit tests for the bamcp.__main__ entry point."""

import os

import pytest
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import PlainTextResponse
from starlette.routing import Route
from starlette.testclient import TestClient

from bamcp import __main__ as entry
from bamcp.config import BAMCPConfig


def _app() -> Starlette:
    async def home(request: Request) -> PlainTextResponse:
        return PlainTextResponse("ok")

    return Starlette(routes=[Route("/", home)])


class TestMiddlewareOrder:
    @pytest.mark.unit
    def test_trusted_host_is_outermost_ratelimit_innermost(self):
        from starlette.middleware.trustedhost import TrustedHostMiddleware

        from bamcp.middleware.ratelimit import RateLimitMiddleware
        from bamcp.middleware.security import SecurityHeadersMiddleware

        app = _app()
        entry._add_http_middleware(app, BAMCPConfig(trusted_hosts=["testserver"]))
        # user_middleware[0] is the outermost layer.
        assert [m.cls for m in app.user_middleware] == [
            TrustedHostMiddleware,
            SecurityHeadersMiddleware,
            RateLimitMiddleware,
        ]

    @pytest.mark.unit
    def test_rate_limited_response_carries_security_headers(self):
        app = _app()
        entry._add_http_middleware(app, BAMCPConfig(rate_limit=1))
        client = TestClient(app)
        assert client.get("/").status_code == 200
        resp = client.get("/")
        assert resp.status_code == 429
        assert resp.headers.get("x-content-type-options") == "nosniff"

    @pytest.mark.unit
    def test_bad_host_rejected_before_rate_limit(self):
        app = _app()
        entry._add_http_middleware(app, BAMCPConfig(rate_limit=1, trusted_hosts=["good.example"]))
        client = TestClient(app)
        for _ in range(3):
            assert client.get("/").status_code == 400  # Host: testserver is not allowed
        # The rejected requests never consumed the (single) rate-limit slot.
        assert client.get("/", headers={"host": "good.example"}).status_code == 200


class TestMainConfigErrors:
    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("var", "value"),
        [("BAMCP_TRANSPORT", "carrier-pigeon"), ("BAMCP_PORT", "abc"), ("BAMCP_AUTH_ENABLED", "?")],
    )
    def test_invalid_env_exits_with_named_variable(self, monkeypatch, capsys, var, value):
        for key in list(os.environ):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)
        monkeypatch.setenv(var, value)
        with pytest.raises(SystemExit) as exc:
            entry.main()
        assert exc.value.code == 1
        err = capsys.readouterr().err
        assert err.startswith("Invalid BAMCP configuration:")
        assert value in err

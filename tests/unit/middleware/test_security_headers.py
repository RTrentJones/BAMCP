"""Unit tests for security headers middleware."""

import pytest
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import PlainTextResponse
from starlette.routing import Route
from starlette.testclient import TestClient

from bamcp.middleware.security import SecurityHeadersMiddleware


def _make_app() -> Starlette:
    """Create a test Starlette app with security headers."""

    async def homepage(request: Request) -> PlainTextResponse:
        return PlainTextResponse("ok")

    app = Starlette(routes=[Route("/", homepage)])
    app.add_middleware(SecurityHeadersMiddleware)
    return app


class TestSecurityHeadersMiddleware:
    """Tests for SecurityHeadersMiddleware."""

    @pytest.mark.unit
    def test_x_content_type_options(self):
        client = TestClient(_make_app())
        resp = client.get("/")
        assert resp.headers["X-Content-Type-Options"] == "nosniff"

    @pytest.mark.unit
    def test_x_frame_options(self):
        client = TestClient(_make_app())
        resp = client.get("/")
        assert resp.headers["X-Frame-Options"] == "DENY"

    @pytest.mark.unit
    def test_referrer_policy(self):
        client = TestClient(_make_app())
        resp = client.get("/")
        assert resp.headers["Referrer-Policy"] == "strict-origin-when-cross-origin"

    @pytest.mark.unit
    def test_permissions_policy(self):
        client = TestClient(_make_app())
        resp = client.get("/")
        assert "camera=()" in resp.headers["Permissions-Policy"]
        assert "microphone=()" in resp.headers["Permissions-Policy"]

    @pytest.mark.unit
    def test_content_security_policy(self):
        client = TestClient(_make_app())
        resp = client.get("/")
        csp = resp.headers["Content-Security-Policy"]
        assert "default-src 'self'" in csp
        assert "unsafe-inline" not in csp
        assert "frame-ancestors 'none'" in csp

    @pytest.mark.unit
    def test_does_not_override_existing_headers(self):
        """Middleware should use setdefault, not overwrite."""

        async def homepage(request: Request) -> PlainTextResponse:
            return PlainTextResponse("ok", headers={"X-Frame-Options": "SAMEORIGIN"})

        app = Starlette(routes=[Route("/", homepage)])
        app.add_middleware(SecurityHeadersMiddleware)
        client = TestClient(app)
        resp = client.get("/")
        # Should keep the handler's value, not overwrite with DENY
        assert resp.headers["X-Frame-Options"] == "SAMEORIGIN"

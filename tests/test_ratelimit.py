"""Unit tests for rate limiting middleware."""

import pytest
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import PlainTextResponse
from starlette.routing import Route
from starlette.testclient import TestClient

from bamcp.ratelimit import RateLimitMiddleware


def _make_app(requests_per_minute: int = 5) -> Starlette:
    """Create a test Starlette app with rate limiting."""

    async def homepage(request: Request) -> PlainTextResponse:
        return PlainTextResponse("ok")

    app = Starlette(routes=[Route("/", homepage)])
    app.add_middleware(RateLimitMiddleware, requests_per_minute=requests_per_minute)
    return app


class TestRateLimitMiddleware:
    """Tests for RateLimitMiddleware."""

    @pytest.mark.unit
    def test_allows_requests_under_limit(self):
        client = TestClient(_make_app(requests_per_minute=5))
        for _ in range(5):
            resp = client.get("/")
            assert resp.status_code == 200

    @pytest.mark.unit
    def test_blocks_requests_over_limit(self):
        client = TestClient(_make_app(requests_per_minute=3))
        # Use up the limit
        for _ in range(3):
            resp = client.get("/")
            assert resp.status_code == 200

        # Next request should be rate limited
        resp = client.get("/")
        assert resp.status_code == 429
        assert "Retry-After" in resp.headers

    @pytest.mark.unit
    def test_retry_after_header(self):
        client = TestClient(_make_app(requests_per_minute=1))
        client.get("/")  # Use the one allowed request

        resp = client.get("/")
        assert resp.status_code == 429
        retry_after = int(resp.headers["Retry-After"])
        assert retry_after > 0
        assert retry_after <= 61

    @pytest.mark.unit
    def test_different_ips_tracked_separately(self):
        client = TestClient(_make_app(requests_per_minute=2))

        # Requests from different IPs (via X-Forwarded-For)
        for _ in range(2):
            resp = client.get("/", headers={"X-Forwarded-For": "1.2.3.4"})
            assert resp.status_code == 200

        # IP 1.2.3.4 is now rate limited
        resp = client.get("/", headers={"X-Forwarded-For": "1.2.3.4"})
        assert resp.status_code == 429

        # Different IP should still work
        resp = client.get("/", headers={"X-Forwarded-For": "5.6.7.8"})
        assert resp.status_code == 200

    @pytest.mark.unit
    def test_x_forwarded_for_first_ip(self):
        """Should use the first IP in X-Forwarded-For chain."""
        client = TestClient(_make_app(requests_per_minute=1))

        resp = client.get("/", headers={"X-Forwarded-For": "1.2.3.4, 10.0.0.1"})
        assert resp.status_code == 200

        # Same first IP should be rate limited
        resp = client.get("/", headers={"X-Forwarded-For": "1.2.3.4, 10.0.0.2"})
        assert resp.status_code == 429

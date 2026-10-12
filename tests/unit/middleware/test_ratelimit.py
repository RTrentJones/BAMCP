"""Unit tests for rate limiting middleware."""

from types import SimpleNamespace

import pytest
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import PlainTextResponse
from starlette.routing import Route
from starlette.testclient import TestClient

from bamcp.middleware.ratelimit import RateLimitMiddleware


def _make_app(requests_per_minute: int = 5) -> Starlette:
    """Create a test Starlette app with rate limiting."""

    async def homepage(request: Request) -> PlainTextResponse:
        return PlainTextResponse("ok")

    app = Starlette(routes=[Route("/", homepage)])
    app.add_middleware(RateLimitMiddleware, requests_per_minute=requests_per_minute)
    return app


def _mw(**kwargs) -> RateLimitMiddleware:
    """A middleware instance not attached to an app, for unit-testing helpers."""
    return RateLimitMiddleware(app=lambda *a, **k: None, **kwargs)


def _req(peer: str, xff: str | None = None, cf: str | None = None) -> SimpleNamespace:
    """A minimal stand-in exposing only what _get_client_ip reads."""
    headers = {"x-forwarded-for": xff} if xff is not None else {}
    if cf is not None:
        headers["cf-connecting-ip"] = cf
    return SimpleNamespace(client=SimpleNamespace(host=peer), headers=headers)


class TestRateLimitMiddleware:
    """End-to-end limit behavior (single peer bucket; TestClient peer is untrusted)."""

    @pytest.mark.unit
    def test_allows_requests_under_limit(self):
        client = TestClient(_make_app(requests_per_minute=5))
        for _ in range(5):
            assert client.get("/").status_code == 200

    @pytest.mark.unit
    def test_blocks_requests_over_limit(self):
        client = TestClient(_make_app(requests_per_minute=3))
        for _ in range(3):
            assert client.get("/").status_code == 200
        resp = client.get("/")
        assert resp.status_code == 429
        assert "Retry-After" in resp.headers

    @pytest.mark.unit
    def test_retry_after_header(self):
        client = TestClient(_make_app(requests_per_minute=1))
        client.get("/")
        resp = client.get("/")
        assert resp.status_code == 429
        retry_after = int(resp.headers["Retry-After"])
        assert 0 < retry_after <= 61


class TestProxyAwareClientIP:
    """X-Forwarded-For is honored only behind a trusted proxy — the anti-spoofing fix."""

    @pytest.mark.unit
    def test_untrusted_peer_ignores_xff(self):
        """An untrusted client cannot spoof its identity via XFF."""
        mw = _mw()  # default trusted proxies = loopback
        # Peer is a public IP (not a trusted proxy) → XFF ignored, real socket IP used.
        assert mw._get_client_ip(_req("203.0.113.9", xff="1.2.3.4")) == "203.0.113.9"

    @pytest.mark.unit
    def test_trusted_proxy_honors_xff(self):
        mw = _mw(trusted_proxies=["10.0.0.1"])
        assert mw._get_client_ip(_req("10.0.0.1", xff="1.2.3.4")) == "1.2.3.4"

    @pytest.mark.unit
    def test_loopback_proxy_trusted_by_default(self):
        mw = _mw()  # cloudflared sidecar connects from loopback
        assert mw._get_client_ip(_req("127.0.0.1", xff="1.2.3.4")) == "1.2.3.4"

    @pytest.mark.unit
    def test_returns_rightmost_untrusted_hop(self):
        """With a proxy chain, the client is the rightmost hop that is not a trusted proxy."""
        mw = _mw(trusted_proxies=["127.0.0.1", "10.0.0.0/8"])
        # client -> edge(1.2.3.4) -> proxy(10.0.0.5) -> local proxy(127.0.0.1)
        ip = mw._get_client_ip(_req("127.0.0.1", xff="1.2.3.4, 10.0.0.5"))
        assert ip == "1.2.3.4"

    @pytest.mark.unit
    def test_no_xff_uses_peer(self):
        mw = _mw()
        assert mw._get_client_ip(_req("127.0.0.1")) == "127.0.0.1"


class TestSpoofingIsSeparatedByRealIP:
    """Distinct real clients get distinct buckets; a spoofer behind no proxy shares one."""

    @pytest.mark.unit
    def test_distinct_clients_behind_trusted_proxy(self):
        mw = _mw(requests_per_minute=2, trusted_proxies=["10.0.0.1"])
        import time

        now = time.monotonic()
        mw._requests["1.2.3.4"] = [now, now]  # this client is at its limit
        # A different forwarded client is unaffected.
        assert mw._get_client_ip(_req("10.0.0.1", xff="5.6.7.8")) == "5.6.7.8"


class TestBoundedTable:
    """The tracking table cannot grow without bound (memory-exhaustion defense)."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_table_capped_by_max_tracked_ips(self):
        mw = _mw(requests_per_minute=1000, trusted_proxies=["10.0.0.1"], max_tracked_ips=50)

        async def _next(_request):
            return PlainTextResponse("ok")

        # 500 distinct forwarded clients through the trusted proxy.
        for i in range(500):
            await mw.dispatch(_req("10.0.0.1", xff=f"198.51.{i // 256}.{i % 256}"), _next)
        assert len(mw._requests) <= 50

    @pytest.mark.unit
    def test_expired_entries_swept(self):
        import time

        mw = _mw()
        stale = time.monotonic() - 120  # older than the 60s window
        mw._requests["9.9.9.9"] = [stale]
        mw._sweep_expired(time.monotonic())
        assert "9.9.9.9" not in mw._requests


class TestCloudflareAndValidation:
    """CF-Connecting-IP preference and IP validation of the chosen hop."""

    @pytest.mark.unit
    def test_cf_connecting_ip_preferred_behind_loopback(self):
        mw = _mw()
        req = _req("127.0.0.1", xff="9.9.9.9", cf="203.0.113.7")
        assert mw._get_client_ip(req) == "203.0.113.7"

    @pytest.mark.unit
    def test_cf_connecting_ip_ignored_from_untrusted_peer(self):
        mw = _mw()
        assert mw._get_client_ip(_req("198.51.100.1", cf="203.0.113.7")) == "198.51.100.1"

    @pytest.mark.unit
    def test_invalid_cf_header_falls_back_to_xff(self):
        mw = _mw()
        assert mw._get_client_ip(_req("127.0.0.1", xff="5.6.7.8", cf="not-an-ip")) == "5.6.7.8"

    @pytest.mark.unit
    def test_invalid_xff_hop_falls_back_to_peer(self):
        mw = _mw()
        assert mw._get_client_ip(_req("127.0.0.1", xff="garbage-token")) == "127.0.0.1"

    @pytest.mark.unit
    def test_ipv6_cf_header_is_normalized(self):
        mw = _mw()
        assert mw._get_client_ip(_req("::1", cf="2001:DB8:0:0::1")) == "2001:db8::1"

    @pytest.mark.unit
    def test_default_limit_is_120(self):
        assert _mw().requests_per_minute == 120


class TestAmortizedSweep:
    """Sweeping pops stale entries off the LRU front and stops at the first live one."""

    @pytest.mark.unit
    def test_sweep_stops_at_first_live_entry(self):
        import time

        mw = _mw()
        now = time.monotonic()
        mw._requests["1.1.1.1"] = [now - 200]
        mw._requests["2.2.2.2"] = [now - 120]
        mw._requests["3.3.3.3"] = [now - 1]
        mw._requests["4.4.4.4"] = [now - 300]  # behind a live entry: left for later / the cap
        mw._sweep_expired(now)
        assert list(mw._requests) == ["3.3.3.3", "4.4.4.4"]

    @pytest.mark.unit
    def test_sweep_handles_empty_table(self):
        mw = _mw()
        mw._sweep_expired(0.0)
        assert len(mw._requests) == 0

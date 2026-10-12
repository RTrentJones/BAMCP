"""Tests for the /__version deploy-identity route (Greenlight's artifact-identity gate)."""

import pytest
from starlette.testclient import TestClient

from bamcp.config import BAMCPConfig
from bamcp.server import create_server

SHA = "0123456789abcdef0123456789abcdef01234567"


def _client(config: BAMCPConfig) -> TestClient:
    return TestClient(create_server(config).streamable_http_app())


class TestVersionRoute:
    @pytest.mark.unit
    @pytest.mark.parametrize("path", ["/__version", "/mcp/__version"])
    def test_serves_baked_sha(self, monkeypatch, path):
        monkeypatch.setenv("GREENLIGHT_SHA", SHA)
        resp = _client(BAMCPConfig()).get(path)
        assert resp.status_code == 200
        assert resp.json() == {"sha": SHA}
        assert resp.headers["Cache-Control"] == "no-store"

    @pytest.mark.unit
    @pytest.mark.parametrize("value", [None, ""])
    def test_unset_sha_reports_null(self, monkeypatch, value):
        """A build without the arg must say so (null), never invent an identity."""
        if value is None:
            monkeypatch.delenv("GREENLIGHT_SHA", raising=False)
        else:
            monkeypatch.setenv("GREENLIGHT_SHA", value)
        resp = _client(BAMCPConfig()).get("/mcp/__version")
        assert resp.status_code == 200
        assert resp.json() == {"sha": None}

    @pytest.mark.unit
    def test_public_when_auth_enabled(self, monkeypatch):
        """Identity is served without a token, while the MCP endpoint itself stays protected."""
        monkeypatch.setenv("GREENLIGHT_SHA", SHA)
        client = _client(BAMCPConfig(auth_enabled=True, verify_token="secret"))
        assert client.get("/mcp/__version").json() == {"sha": SHA}
        assert client.post("/mcp", json={}).status_code == 401

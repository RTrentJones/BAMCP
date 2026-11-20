"""Size caps on dynamically registered OAuth client metadata."""

import pytest
from mcp.server.auth.provider import RegistrationError
from mcp.server.auth.routes import create_auth_routes
from mcp.server.auth.settings import ClientRegistrationOptions
from mcp.shared.auth import OAuthClientInformationFull
from pydantic import AnyHttpUrl, AnyUrl
from starlette.applications import Starlette
from starlette.testclient import TestClient

from bamcp.middleware.auth import (
    MAX_CLIENT_METADATA_BYTES,
    MAX_CLIENT_NAME_LENGTH,
    MAX_REDIRECT_URI_LENGTH,
    MAX_REDIRECT_URIS,
    BAMCPAuthProvider,
)


def _client(**overrides) -> OAuthClientInformationFull:
    fields = {
        "client_id": "c1",
        "client_name": "ok",
        "redirect_uris": [AnyUrl("http://localhost:3000/cb")],
        "grant_types": ["authorization_code", "refresh_token"],
        "response_types": ["code"],
        "token_endpoint_auth_method": "none",
    }
    fields.update(overrides)
    return OAuthClientInformationFull(**fields)


@pytest.fixture
def provider():
    return BAMCPAuthProvider(token_expiry=3600, allow_dynamic_registration=True)


@pytest.mark.unit
class TestRegistrationCaps:
    async def test_normal_client_registers(self, provider):
        await provider.register_client(_client())
        assert await provider.get_client("c1") is not None

    async def test_too_many_redirect_uris(self, provider):
        uris = [AnyUrl(f"http://localhost/cb{i}") for i in range(MAX_REDIRECT_URIS + 1)]
        with pytest.raises(RegistrationError) as exc:
            await provider.register_client(_client(redirect_uris=uris))
        assert exc.value.error == "invalid_redirect_uri"
        assert await provider.get_client("c1") is None

    async def test_redirect_uri_too_long(self, provider):
        uri = AnyUrl("http://localhost/" + "a" * MAX_REDIRECT_URI_LENGTH)
        with pytest.raises(RegistrationError):
            await provider.register_client(_client(redirect_uris=[uri]))

    async def test_client_name_too_long(self, provider):
        with pytest.raises(RegistrationError) as exc:
            await provider.register_client(_client(client_name="n" * (MAX_CLIENT_NAME_LENGTH + 1)))
        assert exc.value.error == "invalid_client_metadata"

    async def test_total_size_cap(self, provider):
        # Each field is individually within bounds, but the record is too large overall.
        with pytest.raises(RegistrationError) as exc:
            await provider.register_client(
                _client(
                    client_uri=AnyHttpUrl("http://localhost/" + "u" * 8000),
                    policy_uri=AnyHttpUrl("http://localhost/" + "p" * 8000),
                    tos_uri=AnyHttpUrl("http://localhost/" + "t" * 8000),
                )
            )
        assert exc.value.error == "invalid_client_metadata"

    async def test_total_cap_is_4_kib(self, provider):
        assert MAX_CLIENT_METADATA_BYTES == 4 * 1024
        # ~5 KiB of individually in-bounds fields: accepted under the old 16 KiB cap.
        with pytest.raises(RegistrationError) as exc:
            await provider.register_client(
                _client(
                    redirect_uris=[
                        AnyUrl(f"https://app.example/cb{i}/" + "r" * 1000) for i in range(5)
                    ]
                )
            )
        assert exc.value.error == "invalid_client_metadata"

    async def test_realistic_hosted_client_fits(self, provider):
        # A claude.ai / desktop-style registration with several callbacks and branding URIs.
        client = _client(
            client_id="x" * 43,
            client_secret="s" * 43,
            client_id_issued_at=1_760_000_000,
            client_secret_expires_at=0,
            client_name="Claude Desktop (bamcp)",
            redirect_uris=[
                AnyUrl("https://claude.ai/api/mcp/auth_callback"),
                AnyUrl("https://claude.com/api/mcp/auth_callback"),
                AnyUrl("http://localhost:6274/oauth/callback"),
                AnyUrl("http://localhost:6274/oauth/callback/debug"),
                AnyUrl("http://127.0.0.1:33418/callback"),
            ],
            token_endpoint_auth_method="client_secret_post",
            scope="bamcp:read",
            client_uri=AnyHttpUrl("https://claude.ai"),
            logo_uri=AnyHttpUrl("https://claude.ai/images/claude_app_icon.png"),
            software_id="claude-desktop",
            software_version="1.0.0",
        )
        assert len(client.model_dump_json().encode()) < MAX_CLIENT_METADATA_BYTES // 2
        await provider.register_client(client)


@pytest.mark.unit
def test_register_endpoint_returns_400(provider):
    routes = create_auth_routes(
        provider,
        AnyHttpUrl("http://localhost:8000"),
        client_registration_options=ClientRegistrationOptions(enabled=True),
    )
    client = TestClient(Starlette(routes=routes))
    body = {
        "redirect_uris": [f"http://localhost/cb{i}" for i in range(MAX_REDIRECT_URIS + 1)],
        "token_endpoint_auth_method": "none",
    }
    resp = client.post("/register", json=body)
    assert resp.status_code == 400
    assert resp.json()["error"] == "invalid_redirect_uri"
    ok = client.post("/register", json={**body, "redirect_uris": ["http://localhost/cb"]})
    assert ok.status_code == 201

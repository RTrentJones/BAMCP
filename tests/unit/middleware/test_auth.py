"""Unit tests for bamcp.auth module."""

import time

import pytest
from mcp.server.auth.provider import AuthorizationParams
from mcp.shared.auth import OAuthClientInformationFull
from pydantic import AnyUrl

from bamcp.config import BAMCPConfig
from bamcp.middleware.auth import BAMCPAuthProvider, build_auth_settings


@pytest.fixture
def provider():
    """Auth provider with 1-hour token expiry.

    Dynamic registration is enabled here because these tests exercise the
    authorization-code flow mechanics, which begin with register_client. Prod
    defaults to registration OFF — see TestDynamicRegistrationDisabled.
    """
    return BAMCPAuthProvider(token_expiry=3600, allow_dynamic_registration=True)


@pytest.fixture
def client_info():
    """A registered OAuth client."""
    return OAuthClientInformationFull(
        client_id="test-client",
        client_name="Test Client",
        redirect_uris=[AnyUrl("http://localhost:3000/callback")],
        grant_types=["authorization_code", "refresh_token"],
        response_types=["code"],
        token_endpoint_auth_method="none",
    )


@pytest.fixture
def auth_params():
    """Standard authorization parameters."""
    return AuthorizationParams(
        state="test-state",
        scopes=["read"],
        code_challenge="test-challenge",
        redirect_uri=AnyUrl("http://localhost:3000/callback"),
        redirect_uri_provided_explicitly=True,
    )


class TestBuildAuthSettings:
    """Tests for build_auth_settings helper."""

    @pytest.mark.unit
    def test_builds_settings(self):
        config = BAMCPConfig(
            issuer_url="http://localhost:8000",
            resource_server_url="http://localhost:8000",
            required_scopes=["read", "write"],
        )
        settings = build_auth_settings(config)
        assert str(settings.issuer_url) == "http://localhost:8000/"
        assert settings.required_scopes == ["read", "write"]
        assert settings.client_registration_options is not None
        # Dynamic registration follows the config flag, which defaults OFF.
        assert settings.client_registration_options.enabled is False
        assert settings.revocation_options is not None
        assert settings.revocation_options.enabled is True

    @pytest.mark.unit
    def test_registration_enabled_when_configured(self):
        config = BAMCPConfig(
            issuer_url="http://localhost:8000",
            resource_server_url="http://localhost:8000",
            allow_dynamic_registration=True,
        )
        settings = build_auth_settings(config)
        assert settings.client_registration_options.enabled is True

    @pytest.mark.unit
    def test_builds_settings_no_scopes(self):
        config = BAMCPConfig(
            issuer_url="http://localhost:8000",
            resource_server_url="http://localhost:8000",
        )
        settings = build_auth_settings(config)
        assert settings.required_scopes == []


class TestBAMCPAuthProvider:
    """Tests for the in-memory OAuth provider."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_register_and_get_client(self, provider, client_info):
        await provider.register_client(client_info)
        result = await provider.get_client("test-client")
        assert result is not None
        assert result.client_id == "test-client"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_get_unknown_client(self, provider):
        result = await provider.get_client("nonexistent")
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_authorize_returns_redirect(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        assert "http://localhost:3000/callback" in redirect
        assert "code=" in redirect
        assert "state=test-state" in redirect

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_full_auth_code_flow(self, provider, client_info, auth_params):
        """Test the complete authorization code -> token exchange flow."""
        await provider.register_client(client_info)

        # Authorize and extract code
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]

        # Load and exchange code
        auth_code = await provider.load_authorization_code(client_info, code)
        assert auth_code is not None
        assert auth_code.client_id == "test-client"

        token = await provider.exchange_authorization_code(client_info, auth_code)
        assert token.access_token
        assert token.refresh_token
        assert token.token_type == "Bearer"
        assert token.expires_in == 3600

        # Code should be consumed (single-use)
        reused = await provider.load_authorization_code(client_info, code)
        assert reused is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_access_token_validation(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]
        auth_code = await provider.load_authorization_code(client_info, code)
        token = await provider.exchange_authorization_code(client_info, auth_code)

        # Validate the access token
        access = await provider.load_access_token(token.access_token)
        assert access is not None
        assert access.client_id == "test-client"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_expired_token_rejected(self, provider, client_info, auth_params):
        """Expired tokens should not validate."""
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]
        auth_code = await provider.load_authorization_code(client_info, code)
        token = await provider.exchange_authorization_code(client_info, auth_code)

        # Manually expire the token
        provider._access_tokens[token.access_token].expires_at = int(time.time()) - 1

        access = await provider.load_access_token(token.access_token)
        assert access is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_refresh_token_flow(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]
        auth_code = await provider.load_authorization_code(client_info, code)
        token = await provider.exchange_authorization_code(client_info, auth_code)

        # Use refresh token
        rt = await provider.load_refresh_token(client_info, token.refresh_token)
        assert rt is not None
        new_token = await provider.exchange_refresh_token(client_info, rt, ["read"])
        assert new_token.access_token
        assert new_token.access_token != token.access_token

        # Old refresh token should be consumed
        reused = await provider.load_refresh_token(client_info, token.refresh_token)
        assert reused is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_revoke_access_token(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]
        auth_code = await provider.load_authorization_code(client_info, code)
        token = await provider.exchange_authorization_code(client_info, auth_code)

        access = await provider.load_access_token(token.access_token)
        assert access is not None
        await provider.revoke_token(access)

        revoked = await provider.load_access_token(token.access_token)
        assert revoked is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_revoke_refresh_token(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]
        auth_code = await provider.load_authorization_code(client_info, code)
        token = await provider.exchange_authorization_code(client_info, auth_code)

        rt = await provider.load_refresh_token(client_info, token.refresh_token)
        assert rt is not None
        await provider.revoke_token(rt)

        revoked = await provider.load_refresh_token(client_info, token.refresh_token)
        assert revoked is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_expired_auth_code_rejected(self, provider, client_info):
        """Auth codes with expired timestamps should be rejected."""
        await provider.register_client(client_info)
        params = AuthorizationParams(
            state=None,
            scopes=[],
            code_challenge="c",
            redirect_uri=AnyUrl("http://localhost:3000/callback"),
            redirect_uri_provided_explicitly=True,
        )
        redirect = await provider.authorize(client_info, params)
        code = redirect.split("code=")[1].split("&")[0]

        # Manually expire the code
        provider._auth_codes[code].expires_at = time.time() - 1

        result = await provider.load_authorization_code(client_info, code)
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_wrong_client_cannot_load_code(self, provider, client_info, auth_params):
        await provider.register_client(client_info)
        redirect = await provider.authorize(client_info, auth_params)
        code = redirect.split("code=")[1].split("&")[0]

        other = OAuthClientInformationFull(
            client_id="other-client",
            redirect_uris=[AnyUrl("http://localhost:3000/callback")],
            grant_types=["authorization_code"],
            response_types=["code"],
            token_endpoint_auth_method="none",
        )
        result = await provider.load_authorization_code(other, code)
        assert result is None


class TestDynamicRegistrationDisabled:
    """Adversarial: with the locked-down default, an anonymous client cannot get access.

    The reviewed 'OAuth' was a token-minting service — open dynamic registration plus a
    no-consent authorize let anyone register and mint a token. With registration OFF
    (the default), self-registration is refused, so the authorization-code flow is
    unreachable and the only accepted credential is the configured service token.
    """

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_anonymous_client_cannot_self_register(self, client_info):
        provider = BAMCPAuthProvider()  # default: registration disabled
        with pytest.raises(ValueError, match="registration is disabled"):
            await provider.register_client(client_info)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_unregistered_client_is_unknown(self, client_info):
        """Without registration, the client is never known → SDK rejects before authorize."""
        provider = BAMCPAuthProvider()
        assert await provider.get_client(client_info.client_id) is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_forged_bearer_is_rejected_without_service_token(self):
        """No service token configured → no bearer is ever valid."""
        provider = BAMCPAuthProvider()
        assert await provider.load_access_token("anything-a-client-made-up") is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_only_the_service_token_is_accepted(self):
        provider = BAMCPAuthProvider(verify_token="svc-secret", verify_scopes=["bamcp:read"])
        assert await provider.load_access_token("svc-secret") is not None
        assert await provider.load_access_token("svc-secre") is None  # near-miss rejected

    @pytest.mark.unit
    def test_config_default_locks_down_registration(self, monkeypatch):
        monkeypatch.setenv("BAMCP_AUTH_ENABLED", "true")
        monkeypatch.delenv("BAMCP_ALLOW_DYNAMIC_REGISTRATION", raising=False)
        monkeypatch.delenv("BAMCP_REQUIRED_SCOPES", raising=False)
        cfg = BAMCPConfig.from_env()
        assert cfg.allow_dynamic_registration is False
        # Auth-on with no explicit scopes → a scope is required, never scope-less.
        assert cfg.required_scopes == ["bamcp:read"]


class TestServiceToken:
    """M2M service token (CI verify) — stateless, restart-proof, read scope."""

    @pytest.mark.asyncio
    async def test_service_token_is_accepted_with_scopes_and_no_expiry(self):
        provider = BAMCPAuthProvider(verify_token="svc-secret", verify_scopes=["read"])
        at = await provider.load_access_token("svc-secret")
        assert at is not None
        assert at.scopes == ["read"]
        assert at.expires_at is None  # never expires → survives restarts

    @pytest.mark.asyncio
    async def test_wrong_token_is_rejected(self):
        provider = BAMCPAuthProvider(verify_token="svc-secret", verify_scopes=["read"])
        assert await provider.load_access_token("not-the-token") is None

    @pytest.mark.asyncio
    async def test_disabled_when_unset(self):
        provider = BAMCPAuthProvider()  # no verify_token
        assert await provider.load_access_token("svc-secret") is None

    @pytest.mark.asyncio
    async def test_config_reads_bamcp_verify_token(self, monkeypatch):
        monkeypatch.setenv("BAMCP_VERIFY_TOKEN", "from-env")
        cfg = BAMCPConfig.from_env()
        assert cfg.verify_token == "from-env"


class TestServiceTokenNonAscii:
    """A latin-1-decoded non-ASCII bearer must be rejected, never raise (was a 500)."""

    @pytest.mark.asyncio
    @pytest.mark.parametrize("bearer", ["\xe9", "svc-secr\xe9t", "☃", "\ud800"])
    async def test_non_ascii_bearer_is_rejected(self, bearer):
        provider = BAMCPAuthProvider(verify_token="svc-secret", verify_scopes=["read"])
        assert await provider.load_access_token(bearer) is None

    @pytest.mark.asyncio
    async def test_non_ascii_service_token_still_matches(self):
        provider = BAMCPAuthProvider(verify_token="s\xe9cret", verify_scopes=["read"])
        assert await provider.load_access_token("s\xe9cret") is not None
        assert await provider.load_access_token("secret") is None


def _client(cid: str) -> OAuthClientInformationFull:
    return OAuthClientInformationFull(
        client_id=cid,
        redirect_uris=[AnyUrl("http://localhost:3000/callback")],
        grant_types=["authorization_code", "refresh_token"],
        response_types=["code"],
        token_endpoint_auth_method="none",
    )


async def _mint(provider, client, params):
    redirect = await provider.authorize(client, params)
    code = redirect.split("code=")[1].split("&")[0]
    auth_code = await provider.load_authorization_code(client, code)
    return await provider.exchange_authorization_code(client, auth_code)


class TestBoundedStores:
    """The open demo lets anyone register and mint tokens — memory must stay bounded."""

    @pytest.mark.asyncio
    async def test_client_cap_evicts_oldest(self):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True, max_clients=3)
        for i in range(5):
            await provider.register_client(_client(f"c{i}"))
        assert list(provider._clients) == ["c2", "c3", "c4"]
        assert await provider.get_client("c0") is None
        assert await provider.get_client("c4") is not None

    @pytest.mark.asyncio
    async def test_reregistration_refreshes_client_slot(self):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True, max_clients=2)
        await provider.register_client(_client("a"))
        await provider.register_client(_client("b"))
        await provider.register_client(_client("a"))
        await provider.register_client(_client("c"))
        assert list(provider._clients) == ["a", "c"]

    @pytest.mark.asyncio
    async def test_token_caps_evict_oldest(self, client_info, auth_params):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True, max_tokens=2)
        await provider.register_client(client_info)
        first = await _mint(provider, client_info, auth_params)
        for _ in range(3):
            await _mint(provider, client_info, auth_params)
        assert len(provider._access_tokens) == 2
        assert len(provider._refresh_tokens) == 2
        assert len(provider._refresh_to_access) <= 2
        assert await provider.load_access_token(first.access_token) is None
        for _ in range(3):
            await provider.authorize(client_info, auth_params)
        assert len(provider._auth_codes) == 2

    @pytest.mark.asyncio
    async def test_refresh_tokens_have_a_ttl(self, client_info, auth_params, monkeypatch):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True, refresh_token_expiry=100)
        await provider.register_client(client_info)
        token = await _mint(provider, client_info, auth_params)
        rt = await provider.load_refresh_token(client_info, token.refresh_token)
        assert rt is not None and rt.expires_at is not None
        real = time.time()
        monkeypatch.setattr(time, "time", lambda: real + 101)
        assert await provider.load_refresh_token(client_info, token.refresh_token) is None
        assert token.refresh_token not in provider._refresh_tokens

    @pytest.mark.asyncio
    async def test_rotation_revokes_old_access_token(self, client_info, auth_params):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True)
        await provider.register_client(client_info)
        token = await _mint(provider, client_info, auth_params)
        rt = await provider.load_refresh_token(client_info, token.refresh_token)
        new = await provider.exchange_refresh_token(client_info, rt, [])
        assert await provider.load_access_token(token.access_token) is None
        assert await provider.load_access_token(new.access_token) is not None
        assert token.refresh_token not in provider._refresh_to_access

    @pytest.mark.asyncio
    async def test_revoking_refresh_revokes_paired_access(self, client_info, auth_params):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True)
        await provider.register_client(client_info)
        token = await _mint(provider, client_info, auth_params)
        rt = await provider.load_refresh_token(client_info, token.refresh_token)
        await provider.revoke_token(rt)
        assert await provider.load_access_token(token.access_token) is None

    @pytest.mark.asyncio
    async def test_periodic_sweep_drops_expired_entries(
        self, client_info, auth_params, monkeypatch
    ):
        provider = BAMCPAuthProvider(
            allow_dynamic_registration=True,
            token_expiry=10,
            refresh_token_expiry=20,
            sweep_interval=60,
        )
        await provider.register_client(client_info)
        await _mint(provider, client_info, auth_params)
        await provider.authorize(client_info, auth_params)  # an unused code
        assert provider._access_tokens and provider._auth_codes and provider._refresh_tokens

        real_wall, real_mono = time.time(), time.monotonic()
        monkeypatch.setattr(time, "time", lambda: real_wall + 1000)
        monkeypatch.setattr(time, "monotonic", lambda: real_mono + 61)
        await provider.register_client(_client("trigger"))
        assert provider._auth_codes == {}
        assert provider._access_tokens == {}
        assert provider._refresh_tokens == {}
        assert provider._refresh_to_access == {}

    @pytest.mark.asyncio
    async def test_sweep_is_rate_limited(self, client_info, auth_params, monkeypatch):
        provider = BAMCPAuthProvider(allow_dynamic_registration=True, token_expiry=10)
        await provider.register_client(client_info)
        await _mint(provider, client_info, auth_params)
        real_wall = time.time()
        monkeypatch.setattr(time, "time", lambda: real_wall + 1000)
        # Within sweep_interval of construction: the expired token stays until lookup.
        await provider.register_client(_client("other"))
        assert len(provider._access_tokens) == 1

    @pytest.mark.asyncio
    async def test_sweep_drops_orphaned_pairings(self):
        provider = BAMCPAuthProvider()
        provider._refresh_to_access["gone"] = "x"
        provider._sweep_expired(time.time())
        assert provider._refresh_to_access == {}

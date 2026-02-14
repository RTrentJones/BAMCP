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
    """Auth provider with 1-hour token expiry."""
    return BAMCPAuthProvider(token_expiry=3600)


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
        assert settings.client_registration_options.enabled is True
        assert settings.revocation_options is not None
        assert settings.revocation_options.enabled is True

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

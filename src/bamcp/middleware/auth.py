"""OAuth 2.0 Authorization Server provider for BAMCP.

Implements the MCP SDK's OAuthAuthorizationServerProvider protocol with
in-memory storage.  Suitable for single-instance deployments.  For
multi-instance / HA setups, swap the dicts for a shared store (Redis, DB).
"""

from __future__ import annotations

import secrets
import time
from typing import cast
from urllib.parse import urlencode

from mcp.server.auth.provider import (
    AccessToken,
    AuthorizationCode,
    AuthorizationParams,
    OAuthAuthorizationServerProvider,
    RefreshToken,
)
from mcp.server.auth.settings import (
    AuthSettings,
    ClientRegistrationOptions,
    RevocationOptions,
)
from mcp.shared.auth import OAuthClientInformationFull, OAuthToken
from pydantic import AnyHttpUrl

from ..config import BAMCPConfig


def build_auth_settings(config: BAMCPConfig) -> AuthSettings:
    """Build AuthSettings from BAMCPConfig."""
    return AuthSettings(
        issuer_url=cast(AnyHttpUrl, config.issuer_url),
        resource_server_url=cast(AnyHttpUrl, config.resource_server_url),
        required_scopes=config.required_scopes or [],
        client_registration_options=ClientRegistrationOptions(
            enabled=True,
            valid_scopes=config.required_scopes or [],
            default_scopes=config.required_scopes or [],
        ),
        revocation_options=RevocationOptions(enabled=True),
    )


class BAMCPAuthProvider(
    OAuthAuthorizationServerProvider[AuthorizationCode, RefreshToken, AccessToken],
):
    """In-memory OAuth 2.0 authorization server for BAMCP."""

    def __init__(self, token_expiry: int = 3600) -> None:
        self.token_expiry = token_expiry
        self._clients: dict[str, OAuthClientInformationFull] = {}
        self._auth_codes: dict[str, AuthorizationCode] = {}
        self._access_tokens: dict[str, AccessToken] = {}
        self._refresh_tokens: dict[str, RefreshToken] = {}

    # -- Client management ---------------------------------------------------

    async def get_client(self, client_id: str) -> OAuthClientInformationFull | None:
        return self._clients.get(client_id)

    async def register_client(self, client_info: OAuthClientInformationFull) -> None:
        if client_info.client_id is not None:
            self._clients[client_info.client_id] = client_info

    # -- Authorization -------------------------------------------------------

    async def authorize(
        self,
        client: OAuthClientInformationFull,
        params: AuthorizationParams,
    ) -> str:
        """Generate an authorization code and return the redirect URI."""
        code = secrets.token_urlsafe(32)
        client_id = client.client_id or ""
        self._auth_codes[code] = AuthorizationCode(
            code=code,
            scopes=params.scopes or [],
            expires_at=time.time() + 300,  # 5-minute code lifetime
            client_id=client_id,
            code_challenge=params.code_challenge,
            redirect_uri=params.redirect_uri,
            redirect_uri_provided_explicitly=params.redirect_uri_provided_explicitly,
        )

        query = urlencode({"code": code, **({"state": params.state} if params.state else {})})
        return f"{params.redirect_uri}?{query}"

    # -- Code exchange -------------------------------------------------------

    async def load_authorization_code(
        self,
        client: OAuthClientInformationFull,
        authorization_code: str,
    ) -> AuthorizationCode | None:
        code_obj = self._auth_codes.get(authorization_code)
        if code_obj is None:
            return None
        if code_obj.client_id != (client.client_id or ""):
            return None
        if code_obj.expires_at < time.time():
            self._auth_codes.pop(authorization_code, None)
            return None
        return code_obj

    async def exchange_authorization_code(
        self,
        client: OAuthClientInformationFull,
        authorization_code: AuthorizationCode,
    ) -> OAuthToken:
        self._auth_codes.pop(authorization_code.code, None)
        return self._issue_token(client.client_id or "", authorization_code.scopes)

    # -- Refresh token -------------------------------------------------------

    async def load_refresh_token(
        self,
        client: OAuthClientInformationFull,
        refresh_token: str,
    ) -> RefreshToken | None:
        rt = self._refresh_tokens.get(refresh_token)
        if rt is None:
            return None
        if rt.client_id != (client.client_id or ""):
            return None
        return rt

    async def exchange_refresh_token(
        self,
        client: OAuthClientInformationFull,
        refresh_token: RefreshToken,
        scopes: list[str],
    ) -> OAuthToken:
        self._refresh_tokens.pop(refresh_token.token, None)
        effective_scopes = scopes if scopes else refresh_token.scopes
        return self._issue_token(client.client_id or "", effective_scopes)

    # -- Access token verification -------------------------------------------

    async def load_access_token(self, token: str) -> AccessToken | None:
        at = self._access_tokens.get(token)
        if at is None:
            return None
        if at.expires_at is not None and at.expires_at < int(time.time()):
            self._access_tokens.pop(token, None)
            return None
        return at

    # -- Revocation ----------------------------------------------------------

    async def revoke_token(
        self,
        token: AccessToken | RefreshToken,
    ) -> None:
        if isinstance(token, AccessToken):
            self._access_tokens.pop(token.token, None)
        else:
            self._refresh_tokens.pop(token.token, None)

    # -- Helpers -------------------------------------------------------------

    def _issue_token(self, client_id: str, scopes: list[str]) -> OAuthToken:
        access = secrets.token_urlsafe(32)
        refresh = secrets.token_urlsafe(32)
        now = int(time.time())

        self._access_tokens[access] = AccessToken(
            token=access,
            client_id=client_id,
            scopes=scopes,
            expires_at=now + self.token_expiry,
        )
        self._refresh_tokens[refresh] = RefreshToken(
            token=refresh,
            client_id=client_id,
            scopes=scopes,
        )

        return OAuthToken(
            access_token=access,
            token_type="Bearer",  # noqa: S106
            expires_in=self.token_expiry,
            refresh_token=refresh,
            scope=" ".join(scopes) if scopes else None,
        )

"""OAuth 2.0 Authorization Server provider for BAMCP.

Implements the MCP SDK's OAuthAuthorizationServerProvider protocol with
in-memory storage.  Suitable for single-instance deployments.  For
multi-instance / HA setups, swap the dicts for a shared store (Redis, DB).
"""

from __future__ import annotations

import secrets
import time
from typing import TypeVar, cast
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

_K = TypeVar("_K")
_V = TypeVar("_V")


def build_auth_settings(config: BAMCPConfig) -> AuthSettings:
    """Build AuthSettings from BAMCPConfig.

    Dynamic client registration is gated on ``config.allow_dynamic_registration``
    (default OFF). With it off, the ``/register`` endpoint is not exposed, so an
    anonymous client cannot self-register and mint a token — the only accepted
    credential is the configured service token (``verify_token``) or a client
    registered out-of-band. This is the honest access-control posture; an open
    registration endpoint with a no-consent authorize is not access control.
    """
    return AuthSettings(
        issuer_url=cast(AnyHttpUrl, config.issuer_url),
        resource_server_url=cast(AnyHttpUrl, config.resource_server_url),
        required_scopes=config.required_scopes or [],
        client_registration_options=ClientRegistrationOptions(
            enabled=config.allow_dynamic_registration,
            valid_scopes=config.required_scopes or [],
            default_scopes=config.required_scopes or [],
        ),
        revocation_options=RevocationOptions(enabled=True),
    )


# Resource bounds for the in-memory stores. The demo is open (dynamic registration +
# auto-approve), so anyone can register clients and mint tokens; without these caps and the
# periodic sweep, expired entries would only ever be dropped on lookup and the dicts would grow
# for the life of the process.
DEFAULT_REFRESH_TOKEN_EXPIRY = 30 * 24 * 3600  # 30 days
DEFAULT_MAX_CLIENTS = 10_000
DEFAULT_MAX_TOKENS = 50_000  # per store: auth codes, access tokens, refresh tokens
DEFAULT_SWEEP_INTERVAL = 60.0  # seconds between full expiry sweeps
AUTH_CODE_LIFETIME = 300  # 5-minute code lifetime


def _to_bytes(value: str) -> bytes:
    """Encode a bearer for constant-time comparison.

    ``secrets.compare_digest`` raises TypeError on non-ASCII ``str`` — and Starlette decodes
    header bytes as latin-1, so ``Authorization: Bearer \xe9`` would surface as a 500. Comparing
    bytes accepts any input; ``surrogatepass`` keeps even a lone surrogate from raising.
    """
    return value.encode("utf-8", "surrogatepass")


def _evict_oldest(store: dict[_K, _V], cap: int) -> list[_V]:
    """Drop insertion-oldest entries until ``len(store) <= cap``; return the evicted values."""
    evicted: list[_V] = []
    while len(store) > cap:
        key = next(iter(store))
        evicted.append(store.pop(key))
    return evicted


class BAMCPAuthProvider(
    OAuthAuthorizationServerProvider[AuthorizationCode, RefreshToken, AccessToken],
):
    """In-memory OAuth 2.0 authorization server for BAMCP.

    Memory is bounded: expired codes/tokens are swept at most every ``sweep_interval`` seconds
    (amortized over mutating calls, so no request pays a full scan more than once a minute),
    refresh tokens carry a TTL, each store is capped (oldest evicted first), and a refresh-token
    rotation revokes the access token issued alongside the old refresh token.
    """

    def __init__(
        self,
        token_expiry: int = 3600,
        verify_token: str | None = None,
        verify_scopes: list[str] | None = None,
        allow_dynamic_registration: bool = False,
        *,
        refresh_token_expiry: int = DEFAULT_REFRESH_TOKEN_EXPIRY,
        max_clients: int = DEFAULT_MAX_CLIENTS,
        max_tokens: int = DEFAULT_MAX_TOKENS,
        sweep_interval: float = DEFAULT_SWEEP_INTERVAL,
    ) -> None:
        self.token_expiry = token_expiry
        self.refresh_token_expiry = refresh_token_expiry
        self.max_clients = max_clients
        self.max_tokens = max_tokens
        self.sweep_interval = sweep_interval
        # Optional M2M service token (CI verify): a single configured bearer accepted statelessly,
        # so it survives restarts (the OAuth stores below are in-memory and reset on every restart).
        self._verify_token = verify_token or None
        self._verify_token_bytes = _to_bytes(verify_token) if verify_token else None
        self._verify_scopes = verify_scopes or []
        # Dynamic client registration (default OFF). When off, register_client refuses, so no
        # anonymous client can self-register and mint a token — access is the service token only.
        self._allow_dynamic_registration = allow_dynamic_registration
        # Plain dicts preserve insertion order, so the first key is always the oldest entry.
        self._clients: dict[str, OAuthClientInformationFull] = {}
        self._auth_codes: dict[str, AuthorizationCode] = {}
        self._access_tokens: dict[str, AccessToken] = {}
        self._refresh_tokens: dict[str, RefreshToken] = {}
        # refresh token -> the access token issued with it (revoked together on rotation).
        self._refresh_to_access: dict[str, str] = {}
        self._last_sweep = time.monotonic()

    # -- Housekeeping --------------------------------------------------------

    def _maybe_sweep(self) -> None:
        """Run a full expiry sweep if ``sweep_interval`` has elapsed since the last one."""
        mono = time.monotonic()
        if mono - self._last_sweep >= self.sweep_interval:
            self._last_sweep = mono
            self._sweep_expired(time.time())

    def _sweep_expired(self, now: float) -> None:
        """Drop every expired auth code, access token and refresh token."""
        for code, ac in list(self._auth_codes.items()):
            if ac.expires_at < now:
                del self._auth_codes[code]
        for tok, at in list(self._access_tokens.items()):
            if at.expires_at is not None and at.expires_at < now:
                del self._access_tokens[tok]
        for tok, rt in list(self._refresh_tokens.items()):
            if rt.expires_at is not None and rt.expires_at < now:
                self._drop_refresh(tok)
        # Pairings whose refresh token is gone (e.g. evicted) are dead weight.
        for tok in [t for t in self._refresh_to_access if t not in self._refresh_tokens]:
            del self._refresh_to_access[tok]

    def _drop_refresh(self, token: str, *, revoke_access: bool = False) -> None:
        self._refresh_tokens.pop(token, None)
        access = self._refresh_to_access.pop(token, None)
        if revoke_access and access is not None:
            self._access_tokens.pop(access, None)

    # -- Client management ---------------------------------------------------

    async def get_client(self, client_id: str) -> OAuthClientInformationFull | None:
        return self._clients.get(client_id)

    async def register_client(self, client_info: OAuthClientInformationFull) -> None:
        # Refuse self-registration when dynamic registration is disabled. Defense in depth:
        # the SDK also withholds the /register endpoint, but this guarantees no client is
        # admitted through this provider regardless of how it is reached.
        if not self._allow_dynamic_registration:
            raise ValueError("Dynamic client registration is disabled")
        self._maybe_sweep()
        if client_info.client_id is not None:
            # Re-registering moves the client to the newest slot.
            self._clients.pop(client_info.client_id, None)
            self._clients[client_info.client_id] = client_info
            _evict_oldest(self._clients, self.max_clients)

    # -- Authorization -------------------------------------------------------

    async def authorize(
        self,
        client: OAuthClientInformationFull,
        params: AuthorizationParams,
    ) -> str:
        """Generate an authorization code and return the redirect URI."""
        self._maybe_sweep()
        code = secrets.token_urlsafe(32)
        client_id = client.client_id or ""
        self._auth_codes[code] = AuthorizationCode(
            code=code,
            scopes=params.scopes or [],
            expires_at=time.time() + AUTH_CODE_LIFETIME,
            client_id=client_id,
            code_challenge=params.code_challenge,
            redirect_uri=params.redirect_uri,
            redirect_uri_provided_explicitly=params.redirect_uri_provided_explicitly,
        )
        _evict_oldest(self._auth_codes, self.max_tokens)

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
        if rt.expires_at is not None and rt.expires_at < time.time():
            self._drop_refresh(refresh_token)
            return None
        return rt

    async def exchange_refresh_token(
        self,
        client: OAuthClientInformationFull,
        refresh_token: RefreshToken,
        scopes: list[str],
    ) -> OAuthToken:
        # Rotation: the old refresh token and the access token minted with it both die.
        self._drop_refresh(refresh_token.token, revoke_access=True)
        effective_scopes = scopes if scopes else refresh_token.scopes
        return self._issue_token(client.client_id or "", effective_scopes)

    # -- Access token verification -------------------------------------------

    async def load_access_token(self, token: str) -> AccessToken | None:
        # M2M service token (CI verify): accepted statelessly with the required scopes and no
        # expiry, so it works right after a restart. Constant-time compare avoids a timing oracle;
        # bytes (not str) so a non-ASCII bearer is simply rejected instead of raising.
        if self._verify_token_bytes is not None and secrets.compare_digest(
            _to_bytes(token), self._verify_token_bytes
        ):
            return AccessToken(
                token=token,
                client_id="bamcp-verify",
                scopes=list(self._verify_scopes),
                expires_at=None,
            )
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
            # RFC 7009 §2.1: revoking a refresh token may also invalidate its access tokens.
            self._drop_refresh(token.token, revoke_access=True)

    # -- Helpers -------------------------------------------------------------

    def _issue_token(self, client_id: str, scopes: list[str]) -> OAuthToken:
        self._maybe_sweep()
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
            expires_at=now + self.refresh_token_expiry,
        )
        self._refresh_to_access[refresh] = access
        _evict_oldest(self._access_tokens, self.max_tokens)
        for rt in _evict_oldest(self._refresh_tokens, self.max_tokens):
            self._refresh_to_access.pop(rt.token, None)

        return OAuthToken(
            access_token=access,
            token_type="Bearer",  # noqa: S106
            expires_in=self.token_expiry,
            refresh_token=refresh,
            scope=" ".join(scopes) if scopes else None,
        )

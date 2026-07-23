"""IP-based rate limiting middleware for Starlette/ASGI."""

from __future__ import annotations

import ipaddress
import time
from collections import OrderedDict

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response

# Loopback is trusted by default: the cloudflared tunnel sidecar is co-located and connects
# from localhost, so it is the immediate proxy in the deployed topology.
_DEFAULT_TRUSTED_PROXIES = ("127.0.0.1", "::1")


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Sliding-window rate limiter keyed by client IP.

    Tracks request timestamps per IP in a sliding window and returns HTTP 429 with a
    Retry-After header when the limit is exceeded.

    Two hardening properties over a naive limiter:

    * **Proxy-aware client IP.** ``X-Forwarded-For`` is honored only when the direct peer is a
      configured trusted proxy; otherwise the header is ignored and the real socket IP is used.
      An untrusted client therefore cannot spoof rotating identities to bypass per-IP limits or
      balloon the tracking table.
    * **Bounded memory.** The per-IP table is an LRU capped at ``max_tracked_ips``; a global
      sweep drops fully-expired entries each request, and the least-recently-seen IP is evicted
      when the cap is reached — so churn/spoofing cannot grow it without bound.

    Args:
        app: The ASGI application.
        requests_per_minute: Maximum requests allowed per minute per IP.
        trusted_proxies: IPs/CIDRs of reverse proxies whose XFF header is trusted. Defaults to
            loopback.
        max_tracked_ips: Hard ceiling on distinct IPs tracked at once.
    """

    def __init__(
        self,
        app,  # noqa: ANN001
        requests_per_minute: int = 60,
        trusted_proxies: list[str] | None = None,
        max_tracked_ips: int = 10_000,
    ):
        super().__init__(app)
        self.requests_per_minute = requests_per_minute
        self.window_seconds = 60
        self.max_tracked_ips = max_tracked_ips
        self._trusted_networks = self._parse_networks(
            trusted_proxies if trusted_proxies is not None else list(_DEFAULT_TRUSTED_PROXIES)
        )
        # IP -> list of request timestamps, ordered by last activity (LRU: newest at the end).
        self._requests: OrderedDict[str, list[float]] = OrderedDict()

    @staticmethod
    def _parse_networks(entries: list[str]) -> list[ipaddress._BaseNetwork]:
        nets: list[ipaddress._BaseNetwork] = []
        for entry in entries:
            try:
                nets.append(ipaddress.ip_network(entry, strict=False))
            except ValueError:
                continue
        return nets

    def _is_trusted_proxy(self, ip: str) -> bool:
        try:
            addr = ipaddress.ip_address(ip)
        except ValueError:
            return False
        return any(addr in net for net in self._trusted_networks)

    def _get_client_ip(self, request: Request) -> str:
        """Resolve the client IP, honoring X-Forwarded-For only behind a trusted proxy."""
        peer = request.client.host if request.client else "unknown"

        # Only trust XFF if the direct peer is a trusted proxy. Walk the header right-to-left and
        # return the first hop that is not itself a trusted proxy — the real client.
        if self._is_trusted_proxy(peer):
            forwarded = request.headers.get("x-forwarded-for")
            if forwarded:
                for hop in reversed([h.strip() for h in forwarded.split(",") if h.strip()]):
                    if not self._is_trusted_proxy(hop):
                        return hop
        return peer

    def _sweep_expired(self, now: float) -> None:
        """Drop entries whose entire window has expired (bounds memory under churn)."""
        cutoff = now - self.window_seconds
        for ip in list(self._requests.keys()):
            self._requests[ip] = [t for t in self._requests[ip] if t > cutoff]
            if not self._requests[ip]:
                del self._requests[ip]

    def _prune_within_window(self, ip: str, now: float) -> list[float]:
        cutoff = now - self.window_seconds
        timestamps = [t for t in self._requests.get(ip, []) if t > cutoff]
        if timestamps:
            self._requests[ip] = timestamps
            self._requests.move_to_end(ip)
        else:
            self._requests.pop(ip, None)
        return timestamps

    def _enforce_cap(self) -> None:
        while len(self._requests) > self.max_tracked_ips:
            # Evict the least-recently-seen IP (front of the OrderedDict).
            self._requests.popitem(last=False)

    async def dispatch(self, request: Request, call_next):  # noqa: ANN001, ANN201
        ip = self._get_client_ip(request)
        now = time.monotonic()

        self._sweep_expired(now)
        timestamps = self._prune_within_window(ip, now)

        if len(timestamps) >= self.requests_per_minute:
            oldest = timestamps[0]
            retry_after = int(self.window_seconds - (now - oldest)) + 1
            return Response(
                content="Rate limit exceeded",
                status_code=429,
                headers={"Retry-After": str(max(retry_after, 1))},
            )

        timestamps.append(now)
        self._requests[ip] = timestamps
        self._requests.move_to_end(ip)
        self._enforce_cap()
        return await call_next(request)

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
      balloon the tracking table. Behind the proxy, ``CF-Connecting-IP`` is preferred and any
      chosen header value must be a valid IP.
    * **Bounded memory.** The per-IP table is an LRU capped at ``max_tracked_ips``; fully-expired
      entries are popped off the LRU front (amortized), and the least-recently-seen IP is
      evicted when the cap is reached — so churn/spoofing cannot grow it without bound.

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
        requests_per_minute: int = 120,
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

    @staticmethod
    def _valid_ip(value: str | None) -> str | None:
        """Return ``value`` normalized if it parses as an IP address, else None."""
        if not value:
            return None
        try:
            return str(ipaddress.ip_address(value.strip()))
        except ValueError:
            return None

    def _get_client_ip(self, request: Request) -> str:
        """Resolve the client IP, honoring proxy headers only behind a trusted proxy.

        Behind a trusted peer (default: loopback, i.e. the co-located cloudflared sidecar):

        1. ``CF-Connecting-IP`` — Cloudflare sets it to the single real client IP and overwrites
           any client-supplied value at the edge, so it is the most reliable source.
        2. ``X-Forwarded-For`` walked right-to-left; the first hop that is not itself a trusted
           proxy is the client.

        Whatever header value is chosen must parse as a real IP address; garbage (or a spoofed
        non-IP token that would mint a fresh bucket per request) falls back to the socket peer.
        """
        peer = request.client.host if request.client else "unknown"
        if not self._is_trusted_proxy(peer):
            return peer

        cf_ip = self._valid_ip(request.headers.get("cf-connecting-ip"))
        if cf_ip is not None:
            return cf_ip

        forwarded = request.headers.get("x-forwarded-for")
        if forwarded:
            for hop in reversed([h.strip() for h in forwarded.split(",") if h.strip()]):
                if not self._is_trusted_proxy(hop):
                    return self._valid_ip(hop) or peer
        return peer

    def _sweep_expired(self, now: float) -> None:
        """Drop fully-expired entries from the LRU front (amortized O(1) per request).

        The table is ordered by last activity, so stale IPs collect at the front: pop while the
        front entry's newest timestamp is outside the window and stop at the first live one.
        Each entry is popped at most once, so the cost is amortized over the requests that
        created it — no full 10k-IP scan per request. (Ordering is approximate — a throttled
        request refreshes recency without a new timestamp — so a few stale entries may linger
        behind a live one; the ``max_tracked_ips`` cap still bounds the table.)
        """
        cutoff = now - self.window_seconds
        while self._requests:
            ip, timestamps = next(iter(self._requests.items()))
            if timestamps and timestamps[-1] > cutoff:
                break
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

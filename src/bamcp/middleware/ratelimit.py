"""IP-based rate limiting middleware for Starlette/ASGI."""

from __future__ import annotations

import time
from collections import defaultdict

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Sliding-window rate limiter keyed by client IP.

    Tracks request timestamps per IP in a sliding window and returns
    HTTP 429 with Retry-After header when the limit is exceeded.

    Args:
        app: The ASGI application.
        requests_per_minute: Maximum requests allowed per minute per IP.
    """

    def __init__(self, app, requests_per_minute: int = 60):  # noqa: ANN001
        super().__init__(app)
        self.requests_per_minute = requests_per_minute
        self.window_seconds = 60
        # IP -> list of request timestamps
        self._requests: dict[str, list[float]] = defaultdict(list)

    def _get_client_ip(self, request: Request) -> str:
        """Extract client IP, respecting X-Forwarded-For behind a proxy."""
        forwarded = request.headers.get("x-forwarded-for")
        if forwarded:
            # Take the first IP (original client)
            return forwarded.split(",")[0].strip()
        return request.client.host if request.client else "unknown"

    def _cleanup_old_entries(self, ip: str, now: float) -> None:
        """Remove timestamps outside the current window."""
        cutoff = now - self.window_seconds
        self._requests[ip] = [t for t in self._requests[ip] if t > cutoff]
        if not self._requests[ip]:
            del self._requests[ip]

    async def dispatch(self, request: Request, call_next):  # noqa: ANN001, ANN201
        ip = self._get_client_ip(request)
        now = time.monotonic()

        self._cleanup_old_entries(ip, now)

        timestamps = self._requests.get(ip, [])
        if len(timestamps) >= self.requests_per_minute:
            # Calculate retry-after from oldest request in window
            oldest = timestamps[0]
            retry_after = int(self.window_seconds - (now - oldest)) + 1
            return Response(
                content="Rate limit exceeded",
                status_code=429,
                headers={"Retry-After": str(max(retry_after, 1))},
            )

        self._requests[ip].append(now)
        return await call_next(request)

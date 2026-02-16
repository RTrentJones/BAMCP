"""Security headers middleware for Starlette/ASGI."""

from __future__ import annotations

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Add standard security headers to all HTTP responses.

    Headers set:
    - X-Content-Type-Options: nosniff
    - X-Frame-Options: DENY
    - Referrer-Policy: strict-origin-when-cross-origin
    - Permissions-Policy: restrictive defaults
    - Content-Security-Policy: basic policy
    """

    async def dispatch(self, request: Request, call_next):  # noqa: ANN001, ANN201
        response = await call_next(request)

        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("Referrer-Policy", "strict-origin-when-cross-origin")
        response.headers.setdefault(
            "Permissions-Policy", "camera=(), microphone=(), geolocation=()"
        )
        csp = (
            "default-src 'self';"
            " script-src 'self';"
            " style-src 'self';"
            " connect-src 'self';"
            " img-src 'self' data:;"
            " frame-ancestors 'none'"
        )
        response.headers.setdefault("Content-Security-Policy", csp)

        return response

"""Entry point for running BAMCP as a module: python -m bamcp."""

import atexit
import sys
from typing import Literal

from .config import BAMCPConfig
from .core.tools import get_cache
from .server import create_server

Transport = Literal["stdio", "sse", "streamable-http"]


def _add_http_middleware(app, config: BAMCPConfig):  # noqa: ANN001, ANN201
    """Add security middleware to a Starlette app for HTTP transports."""
    from starlette.middleware.trustedhost import TrustedHostMiddleware

    from .middleware.ratelimit import RateLimitMiddleware
    from .middleware.security import SecurityHeadersMiddleware

    # Order matters. Starlette's add_middleware() PREPENDS, so the LAST one added is the
    # OUTERMOST layer and sees the request first. The resulting request path is:
    #
    #   TrustedHost (if configured) -> SecurityHeaders -> RateLimit -> app
    #
    # * TrustedHost is outermost: a bad Host header (DNS rebinding) is rejected before any
    #   other work, and never consumes a rate-limit slot.
    # * SecurityHeaders wraps RateLimit, so 429 responses carry the security headers too.
    # * RateLimit is innermost, directly in front of the MCP app.
    #
    # Added innermost-first below.
    app.add_middleware(
        RateLimitMiddleware,
        requests_per_minute=config.rate_limit,
        trusted_proxies=config.rate_limit_trusted_proxies,
        max_tracked_ips=config.rate_limit_max_tracked_ips,
    )

    app.add_middleware(SecurityHeadersMiddleware)

    if config.trusted_hosts:
        app.add_middleware(TrustedHostMiddleware, allowed_hosts=config.trusted_hosts)


def main() -> None:
    """Run the BAMCP MCP server."""
    # BAMCPConfig validates everything (incl. the transport, in __post_init__) and names the
    # offending variable; surface that as a one-line error, not a traceback.
    try:
        config = BAMCPConfig.from_env()
    except ValueError as exc:
        print(f"Invalid BAMCP configuration: {exc}", file=sys.stderr)
        sys.exit(1)

    from .middleware.telemetry import configure_telemetry

    configure_telemetry(
        enabled=config.telemetry_enabled,
        jsonl_path=config.telemetry_path,
        otel_enabled=config.telemetry_otel_enabled,
    )

    # Register cleanup on exit to remove this session's cache files
    def cleanup_on_exit() -> None:
        cache = get_cache(config)
        removed = cache.cleanup_session()
        if removed > 0:
            print(f"Cleaned up {removed} cached index files", file=sys.stderr)

    atexit.register(cleanup_on_exit)

    transport: Transport = config.transport  # type: ignore[assignment]
    server = create_server(config)

    if transport == "stdio":
        server.run(transport="stdio")
    else:
        # For HTTP transports, get the Starlette app and add middleware
        import anyio
        import uvicorn

        app = server.sse_app() if transport == "sse" else server.streamable_http_app()

        _add_http_middleware(app, config)

        async def _serve() -> None:
            uvi_config = uvicorn.Config(
                app,
                host=config.host,
                port=config.port,
                log_level="info",
            )
            uvi_server = uvicorn.Server(uvi_config)
            await uvi_server.serve()

        anyio.run(_serve)


if __name__ == "__main__":
    main()

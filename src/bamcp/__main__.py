"""Entry point for running BAMCP as a module: python -m bamcp."""

import atexit
import sys
from typing import Literal, get_args

from .config import BAMCPConfig
from .core.tools import get_cache
from .server import create_server

Transport = Literal["stdio", "sse", "streamable-http"]
VALID_TRANSPORTS: tuple[str, ...] = get_args(Transport)


def _add_http_middleware(app, config: BAMCPConfig):  # noqa: ANN001, ANN201
    """Add security middleware to a Starlette app for HTTP transports."""
    from starlette.middleware.trustedhost import TrustedHostMiddleware

    from .middleware.ratelimit import RateLimitMiddleware
    from .middleware.security import SecurityHeadersMiddleware

    # Order matters: outermost middleware runs first
    # 1. Rate limiting (reject before doing work)
    app.add_middleware(RateLimitMiddleware, requests_per_minute=config.rate_limit)

    # 2. Security headers (add to all responses)
    app.add_middleware(SecurityHeadersMiddleware)

    # 3. Trusted host validation (DNS rebinding protection)
    if config.trusted_hosts:
        app.add_middleware(TrustedHostMiddleware, allowed_hosts=config.trusted_hosts)


def main() -> None:
    """Run the BAMCP MCP server."""
    config = BAMCPConfig.from_env()

    if config.transport not in VALID_TRANSPORTS:
        print(
            f"Invalid BAMCP_TRANSPORT={config.transport!r}. "
            f"Must be one of: {', '.join(VALID_TRANSPORTS)}",
            file=sys.stderr,
        )
        sys.exit(1)

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

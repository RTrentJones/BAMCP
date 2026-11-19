"""Entry point for running BAMCP as a module: python -m bamcp."""

import atexit
import logging
import os
import sys
from pathlib import Path
from typing import Literal

from .config import BAMCPConfig
from .core.tools import get_cache
from .core.validation import warn_if_remote_unrestricted
from .server import create_server

Transport = Literal["stdio", "sse", "streamable-http"]

logger = logging.getLogger(__name__)


def _configure_htslib_reference_lookup(config: BAMCPConfig) -> dict[str, str]:
    """Keep htslib's CRAM reference lookup (REF_PATH / REF_CACHE) local to this server.

    A CRAM decoded without an explicit reference makes htslib look the reference up by the
    header's M5 checksum along ``REF_PATH`` — older htslib builds default that to the EBI MD5
    service, and any ``http(s)://`` entry is fetched with libcurl, bypassing BAMCP's SSRF
    checks and host allow-list. When the operator has not set them, point both variables at a
    local directory under the BAMCP cache dir: a search path with no URL entry never goes to
    the network, so a reference-less CRAM fails fast instead. An operator-set value is left
    alone (logged with a warning when it names a URL). Returns the values in effect.
    """
    local = str(Path(config.cache_dir) / "hts-ref" / "%2s" / "%2s" / "%s")
    applied: dict[str, str] = {}
    for name in ("REF_PATH", "REF_CACHE"):
        current = os.environ.get(name)
        if current:
            if "://" in current or "URL=" in current:
                logger.warning(
                    "%s=%s lets htslib fetch CRAM references from the network outside "
                    "BAMCP's remote-host allow-list",
                    name,
                    current,
                )
            applied[name] = current
            continue
        os.environ[name] = local
        applied[name] = local
    logger.info(
        "htslib CRAM reference lookup: REF_PATH=%s REF_CACHE=%s",
        applied["REF_PATH"],
        applied["REF_CACHE"],
    )
    return applied


def _add_http_middleware(app, config: BAMCPConfig):  # noqa: ANN001, ANN201
    """Add security middleware to a Starlette app for HTTP transports."""
    from starlette.middleware.trustedhost import TrustedHostMiddleware

    from .middleware.bodylimit import BodySizeLimitMiddleware
    from .middleware.ratelimit import RateLimitMiddleware
    from .middleware.security import SecurityHeadersMiddleware

    # Order matters. Starlette's add_middleware() PREPENDS, so the LAST one added is the
    # OUTERMOST layer and sees the request first. The resulting request path is:
    #
    #   TrustedHost (if configured) -> SecurityHeaders -> BodySizeLimit -> RateLimit -> app
    #
    # * TrustedHost is outermost: a bad Host header (DNS rebinding) is rejected before any
    #   other work, and never consumes a rate-limit slot.
    # * SecurityHeaders only decorates responses (it reads no body), so BodySizeLimit is
    #   effectively the first gate after the Host check, and its 413s carry the security
    #   headers too. An oversized Content-Length is refused before any body is read and
    #   without consuming a rate-limit slot; a streamed body is cut off at the limit.
    # * SecurityHeaders also wraps RateLimit, so 429 responses carry the headers.
    # * RateLimit is innermost, directly in front of the MCP app.
    #
    # Added innermost-first below.
    app.add_middleware(
        RateLimitMiddleware,
        requests_per_minute=config.rate_limit,
        trusted_proxies=config.rate_limit_trusted_proxies,
        max_tracked_ips=config.rate_limit_max_tracked_ips,
    )

    app.add_middleware(BodySizeLimitMiddleware, max_bytes=config.max_request_bytes)

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

    warn_if_remote_unrestricted(config)
    _configure_htslib_reference_lookup(config)

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
            from .core.prewarm import start_prewarm

            # Best-effort background warm-up of the demo datasets; never blocks readiness and
            # never raises (HTTP only — stdio never reaches this branch). Keep a reference so
            # the task isn't garbage-collected mid-flight.
            prewarm_task = start_prewarm(config)  # noqa: F841
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

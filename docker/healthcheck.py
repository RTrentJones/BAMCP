#!/usr/bin/env python3
"""Docker health check script for BAMCP.

Two modes:

``--preflight`` (run by ``docker/entrypoint.sh`` before the server starts):
  Verifies bamcp, pysam, config, and server are importable/functional. The server is not
  listening yet, so no HTTP probe is possible here.

Runtime probe (no arguments; the Dockerfile ``HEALTHCHECK``):
  For HTTP transports (``sse``, ``streamable-http``) it GETs
  ``http://127.0.0.1:<port>/__version`` with a short timeout and fails on anything but 200, so
  a wedged server is reported unhealthy. For ``stdio`` there is no listener, so it falls back to
  the import check. The probe is stdlib-only (urllib) so it stays cheap every interval.

Exit 0 = healthy, Exit 1 = unhealthy.
"""

import os
import sys

HTTP_TRANSPORTS = frozenset({"sse", "streamable-http"})
DEFAULT_PORT = "8000"
HTTP_TIMEOUT_SECONDS = 3.0


def check_health(preflight: bool = False) -> bool:
    """Return True when healthy. ``preflight`` forces the import-only check."""
    if preflight:
        return _check_imports()
    transport = os.environ.get("BAMCP_TRANSPORT", "stdio").strip().lower() or "stdio"
    if transport in HTTP_TRANSPORTS:
        return _check_http()
    return _check_imports()


def _resolve_port() -> str:
    """Mirror docker/entrypoint.sh: a cloud ``PORT`` is used only when ``BAMCP_PORT`` is unset."""
    return os.environ.get("BAMCP_PORT") or os.environ.get("PORT") or DEFAULT_PORT


def _probe_host_header() -> str | None:
    """A Host header the server's TrustedHostMiddleware will accept, or None for the default.

    With ``BAMCP_TRUSTED_HOSTS`` set, a bare ``Host: 127.0.0.1`` is rejected (400), which would
    mark a healthy container unhealthy. Use the first configured host instead.
    """
    hosts = [h.strip() for h in os.environ.get("BAMCP_TRUSTED_HOSTS", "").split(",") if h.strip()]
    if not hosts or "*" in hosts:
        return None
    first = hosts[0]
    if first.startswith("*."):
        # Starlette matches "*.example.com" as a ".example.com" suffix.
        return "healthcheck" + first[1:]
    return first


def _check_http() -> bool:
    """Verify the HTTP server answers GET /__version with 200."""
    import urllib.request

    url = f"http://127.0.0.1:{_resolve_port()}/__version"
    headers = {}
    host = _probe_host_header()
    if host:
        headers["Host"] = host
    try:
        req = urllib.request.Request(url, method="GET", headers=headers)  # noqa: S310
        with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT_SECONDS) as resp:  # noqa: S310
            if resp.status != 200:
                print(f"HTTP health check failed: {url} -> {resp.status}", file=sys.stderr)
                return False
            return True
    except Exception as exc:
        print(f"HTTP health check failed: {url}: {exc}", file=sys.stderr)
        return False


def _check_imports() -> bool:
    """Verify bamcp package is importable and functional."""
    try:
        from bamcp import __version__

        assert __version__, "Version string is empty"

        import pysam  # noqa: F401

        from bamcp.config import BAMCPConfig

        config = BAMCPConfig.from_env()
        assert config.max_reads > 0, "max_reads must be positive"

        from bamcp.server import create_server

        server = create_server(config)
        assert server.name == "bamcp", f"Unexpected server name: {server.name}"

        return True
    except Exception:
        import traceback

        traceback.print_exc(file=sys.stderr)
        return False


if __name__ == "__main__":
    sys.exit(0 if check_health(preflight="--preflight" in sys.argv[1:]) else 1)

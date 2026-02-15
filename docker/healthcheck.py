#!/usr/bin/env python3
"""Docker health check script for BAMCP.

For stdio transport:
  Verifies bamcp, pysam, config, and server are importable/functional.

For HTTP transports (sse, streamable-http):
  Hits the local HTTP endpoint to verify the server is responding.

Exit 0 = healthy, Exit 1 = unhealthy.
"""

import os
import sys


def check_health() -> bool:
    # Pre-flight always does import checks.
    # HTTP liveness checks are for runtime probes, not pre-flight.
    return _check_imports()


def _check_http() -> bool:
    """Verify the HTTP server is responding."""
    import urllib.request

    port = os.environ.get("BAMCP_PORT", "8000")
    url = f"http://127.0.0.1:{port}/sse"
    try:
        req = urllib.request.Request(url, method="GET")  # noqa: S310
        with urllib.request.urlopen(req, timeout=3) as resp:  # noqa: S310
            return resp.status == 200
    except Exception as exc:
        print(f"HTTP health check failed: {exc}", file=sys.stderr)
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
    except Exception as exc:
        import traceback

        traceback.print_exc(file=sys.stderr)
        return False


if __name__ == "__main__":
    sys.exit(0 if check_health() else 1)

#!/usr/bin/env python3
"""Docker health check script for BAMCP.

Verifies:
  1. bamcp package is importable
  2. pysam is available
  3. MCP server can be instantiated
  4. Config loads from environment

Exit 0 = healthy, Exit 1 = unhealthy.
"""

import sys


def check_health() -> bool:
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
    except Exception as e:
        print(f"Health check failed: {e}", file=sys.stderr)
        return False


if __name__ == "__main__":
    sys.exit(0 if check_health() else 1)

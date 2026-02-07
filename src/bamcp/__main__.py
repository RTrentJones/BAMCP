"""Entry point for running BAMCP as a module: python -m bamcp."""

import atexit
import sys
from typing import Literal, get_args

from .config import BAMCPConfig
from .server import create_server
from .tools import get_cache

Transport = Literal["stdio", "sse", "streamable-http"]
VALID_TRANSPORTS: tuple[str, ...] = get_args(Transport)


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
    server.run(transport=transport)


if __name__ == "__main__":
    main()

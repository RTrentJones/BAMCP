"""Entry point for running BAMCP as a module: python -m bamcp."""

import sys
from typing import Literal, get_args

from .config import BAMCPConfig
from .server import create_server

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

    transport: Transport = config.transport  # type: ignore[assignment]
    server = create_server(config)
    server.run(transport=transport)


if __name__ == "__main__":
    main()

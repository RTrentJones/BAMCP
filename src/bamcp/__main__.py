"""Entry point for running BAMCP as a module: python -m bamcp."""

import asyncio
import sys

from mcp.server.stdio import stdio_server

from .config import BAMCPConfig
from .server import create_server


async def main() -> None:
    """Run the BAMCP MCP server over stdio."""
    config = BAMCPConfig.from_env()
    app = create_server(config)

    async with stdio_server() as (read_stream, write_stream):
        await app.run(read_stream, write_stream, app.create_initialization_options())


if __name__ == "__main__":
    asyncio.run(main())

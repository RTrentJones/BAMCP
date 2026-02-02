"""MCP server setup for BAMCP using FastMCP."""

from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from .config import BAMCPConfig
from .resources import get_viewer_html
from .tools import (
    handle_browse_region,
    handle_get_coverage,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
)


def create_server(config: BAMCPConfig | None = None) -> FastMCP:
    """Create and configure the BAMCP MCP server."""
    if config is None:
        config = BAMCPConfig.from_env()

    kwargs: dict = {
        "name": "bamcp",
        "host": config.host,
        "port": config.port,
    }

    if config.auth_enabled:
        from .auth import BAMCPAuthProvider, build_auth_settings

        kwargs["auth_server_provider"] = BAMCPAuthProvider(
            token_expiry=config.token_expiry,
        )
        kwargs["auth"] = build_auth_settings(config)

    mcp = FastMCP(**kwargs)

    # -- Tools ---------------------------------------------------------------
    # Thin wrappers delegate to the existing handlers in tools.py.
    # FastMCP derives the JSON-Schema from the function signature.

    @mcp.tool(
        description="View aligned reads in a genomic region with interactive visualization",
    )
    async def browse_region(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> str:
        result = await handle_browse_region(
            {"file_path": file_path, "region": region, "reference": reference},
            config,
        )
        return str(result["content"][0]["text"])

    @mcp.tool(description="Detect and return variants in a genomic region")
    async def get_variants(
        file_path: str,
        region: str,
        reference: str | None = None,
        min_vaf: float | None = None,
        min_depth: int | None = None,
    ) -> str:
        args: dict = {"file_path": file_path, "region": region}
        if reference is not None:
            args["reference"] = reference
        if min_vaf is not None:
            args["min_vaf"] = min_vaf
        if min_depth is not None:
            args["min_depth"] = min_depth
        result = await handle_get_variants(args, config)
        return str(result["content"][0]["text"])

    @mcp.tool(description="Calculate depth of coverage statistics for a region")
    async def get_coverage(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> str:
        args: dict = {"file_path": file_path, "region": region}
        if reference is not None:
            args["reference"] = reference
        result = await handle_get_coverage(args, config)
        return str(result["content"][0]["text"])

    @mcp.tool(description="List chromosomes/contigs in a BAM/CRAM file header")
    async def list_contigs(
        file_path: str,
        reference: str | None = None,
    ) -> str:
        args: dict = {"file_path": file_path}
        if reference is not None:
            args["reference"] = reference
        result = await handle_list_contigs(args, config)
        return str(result["content"][0]["text"])

    @mcp.tool(
        description="Jump to a specific genomic position and view surrounding reads",
    )
    async def jump_to(
        file_path: str,
        position: int,
        contig: str | None = None,
        window: int | None = None,
        reference: str | None = None,
    ) -> str:
        args: dict = {"file_path": file_path, "position": position}
        if contig is not None:
            args["contig"] = contig
        if window is not None:
            args["window"] = window
        if reference is not None:
            args["reference"] = reference
        result = await handle_jump_to(args, config)
        return str(result["content"][0]["text"])

    # -- Resources -----------------------------------------------------------

    @mcp.resource(
        "ui://bamcp/viewer",
        name="BAMCP Alignment Viewer",
        description="Interactive BAM/CRAM alignment visualization",
        mime_type="text/html+mcp",
    )
    def viewer() -> str:
        return get_viewer_html()

    return mcp

"""MCP server setup for BAMCP using FastMCP."""

from __future__ import annotations

from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.types import CallToolResult, TextContent

from .config import BAMCPConfig
from .resources import get_viewer_html
from .tools import (
    handle_browse_region,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_visualize_region,
)

_VIEWER_URI = "ui://bamcp/viewer"
_VIEWER_META: dict = {"ui": {"resourceUri": _VIEWER_URI}}

# MCP Apps extension capability (SEP-1865)
_MCP_APPS_EXPERIMENTAL: dict[str, dict[str, Any]] = {"io.modelcontextprotocol/ui": {}}


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

    # Advertise MCP Apps extension capability (SEP-1865)
    # Patch create_initialization_options to include experimental capabilities
    _original_create_init_opts = mcp._mcp_server.create_initialization_options

    def _patched_create_init_opts(
        notification_options: Any = None,
        experimental_capabilities: dict[str, dict[str, Any]] | None = None,
    ) -> Any:
        merged = {**_MCP_APPS_EXPERIMENTAL}
        if experimental_capabilities:
            merged.update(experimental_capabilities)
        return _original_create_init_opts(notification_options, merged)

    mcp._mcp_server.create_initialization_options = _patched_create_init_opts

    # -- Tools ---------------------------------------------------------------
    # Thin wrappers delegate to the existing handlers in tools.py.
    # FastMCP derives the JSON-Schema from the function signature.

    @mcp.tool(
        description="View aligned reads in a genomic region with interactive visualization",
        meta=_VIEWER_META,
    )
    async def browse_region(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> CallToolResult:
        result = await handle_browse_region(
            {"file_path": file_path, "region": region, "reference": reference},
            config,
        )
        payload = result.get("_meta", {}).get("ui/init", {})
        return CallToolResult(
            content=[TextContent(type="text", text=result["content"][0]["text"])],
            structuredContent=payload or None,
        )

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
        meta=_VIEWER_META,
    )
    async def jump_to(
        file_path: str,
        position: int,
        contig: str | None = None,
        window: int | None = None,
        reference: str | None = None,
    ) -> CallToolResult:
        args: dict = {"file_path": file_path, "position": position}
        if contig is not None:
            args["contig"] = contig
        if window is not None:
            args["window"] = window
        if reference is not None:
            args["reference"] = reference
        result = await handle_jump_to(args, config)
        payload = result.get("_meta", {}).get("ui/init", {})
        return CallToolResult(
            content=[TextContent(type="text", text=result["content"][0]["text"])],
            structuredContent=payload or None,
        )

    # -- Phase 1: MCP Apps Tools --------------------------------------------

    @mcp.tool(
        description="Visualize aligned reads in a genomic region with interactive MCP Apps viewer",
        meta=_VIEWER_META,
    )
    async def visualize_region(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> CallToolResult:
        result = await handle_visualize_region(
            {"file_path": file_path, "region": region, "reference": reference},
            config,
        )
        payload = result.get("_meta", {}).get("ui/init", {})
        return CallToolResult(
            content=[TextContent(type="text", text=result["content"][0]["text"])],
            structuredContent=payload or None,
        )

    @mcp.tool(
        description=(
            "Get a text summary of a genomic region "
            "(read count, coverage, variants) for LLM reasoning"
        ),
    )
    async def get_region_summary(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> str:
        args: dict = {"file_path": file_path, "region": region}
        if reference is not None:
            args["reference"] = reference
        result = await handle_get_region_summary(args, config)
        return str(result["content"][0]["text"])

    # -- Phase 2: External Database Tools ------------------------------------

    @mcp.tool(
        description=(
            "Look up a variant in ClinVar for clinical significance, review status, "
            "and associated conditions. Research-grade information, not for clinical use."
        ),
    )
    async def lookup_clinvar(
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> str:
        result = await handle_lookup_clinvar(
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
            config,
        )
        return str(result["content"][0]["text"])

    @mcp.tool(
        description=(
            "Look up a variant in gnomAD for population allele frequency data. "
            "Research-grade information, not for clinical use."
        ),
    )
    async def lookup_gnomad(
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> str:
        result = await handle_lookup_gnomad(
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
            config,
        )
        return str(result["content"][0]["text"])

    # -- Resources -----------------------------------------------------------

    @mcp.resource(
        "ui://bamcp/viewer",
        name="BAMCP Alignment Viewer",
        description="Interactive BAM/CRAM alignment visualization",
        mime_type="text/html;profile=mcp-app",
    )
    def viewer() -> str:
        return get_viewer_html()

    return mcp

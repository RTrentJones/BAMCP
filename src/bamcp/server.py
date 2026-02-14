"""MCP server setup for BAMCP using FastMCP."""

from __future__ import annotations

from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.types import CallToolResult, TextContent

from .analysis.curation import handle_get_variant_curation_summary
from .config import BAMCPConfig
from .core.tools import (
    get_cache,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_visualize_region,
)
from .resources import get_viewer_html

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
        from .middleware.auth import BAMCPAuthProvider, build_auth_settings

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

    mcp._mcp_server.create_initialization_options = _patched_create_init_opts  # type: ignore[method-assign]

    # -- Tools ---------------------------------------------------------------
    # Thin wrappers delegate to the existing handlers in tools.py.
    # FastMCP derives the JSON-Schema from the function signature.

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

    @mcp.tool(
        description=(
            "List chromosomes/contigs in a BAM/CRAM file and detect genome build. "
            "Use this first on new BAM files to determine whether reference is GRCh37 or GRCh38. "
            "Returns detected build, confidence, and suggested public reference URL if available."
        ),
    )
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

    # -- Phase 3: Curation Tools ---------------------------------------------

    @mcp.tool(
        description=(
            "Get a detailed curation summary for a specific variant with artifact risk "
            "assessment, quality metrics, and recommendations for clinical interpretation. "
            "Use this for in-depth variant review by genetic curators."
        ),
    )
    async def get_variant_curation_summary(
        file_path: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        window: int = 50,
        reference: str | None = None,
    ) -> str:
        args: dict = {
            "file_path": file_path,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "window": window,
        }
        if reference is not None:
            args["reference"] = reference
        result = await handle_get_variant_curation_summary(args, config)
        return str(result["content"][0]["text"])

    @mcp.tool(
        description=(
            "Search for a gene by symbol and return its genomic coordinates. "
            "Use this to navigate to a gene's location by name (e.g., BRCA1, TP53)."
        ),
    )
    async def search_gene(symbol: str) -> str:
        import json

        from .clients.genes import GeneClient

        client = GeneClient(
            api_key=config.ncbi_api_key,
            genome_build=config.genome_build,
        )

        try:
            result = await client.search(symbol)
        finally:
            await client.close()

        if result is None:
            return json.dumps({"error": f"Gene '{symbol}' not found"})

        return json.dumps(
            {
                "symbol": result.symbol,
                "name": result.name,
                "region": f"{result.chrom}:{result.start}-{result.end}",
                "strand": result.strand,
            }
        )

    # -- Cache Management Tool -----------------------------------------------

    @mcp.tool(description="Clean up this session's BAM index cache files")
    async def cleanup_cache() -> str:
        """Remove cached index files downloaded during this session.

        Each session stores its cache files in an isolated subdirectory,
        so this only affects the current session's files.

        Returns:
            Summary of cleanup action.
        """
        cache = get_cache(config)
        removed = cache.cleanup_session()
        return f"Removed {removed} cached index files from session {cache.session_id}"

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

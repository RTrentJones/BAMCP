"""MCP server setup for BAMCP using FastMCP."""

from __future__ import annotations

import os
from collections.abc import Awaitable, Callable
from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.types import CallToolResult, TextContent
from starlette.requests import Request
from starlette.responses import JSONResponse

from .analysis.acmg import handle_classify_variant
from .analysis.curation import handle_get_variant_curation_summary
from .config import BAMCPConfig
from .core.errors import sanitize_tool_exception
from .core.tools import (
    handle_cleanup_cache,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_scan_variants,
    handle_search_gene,
    handle_visualize_region,
)
from .demos import resolve_demo_alias
from .prompts import SERVER_INSTRUCTIONS, register_demo_surface
from .resources import get_viewer_html
from .tool_specs import AUTO_REFERENCE_NOTE, DETECTION_DEFAULTS_NOTE, REGION_FORMAT_NOTE

_VIEWER_URI = "ui://bamcp/viewer"
_VIEWER_META: dict = {"ui": {"resourceUri": _VIEWER_URI}}

_Handler = Callable[[dict[str, Any], BAMCPConfig], Awaitable[dict]]


def _args(**kwargs: Any) -> dict[str, Any]:
    """Tool arguments with unset optionals dropped, so handlers apply their own defaults."""
    return {k: v for k, v in kwargs.items() if v is not None}


def _text(result: dict) -> str:
    return str(result["content"][0]["text"])


def _viewer_result(result: dict) -> CallToolResult:
    """Text summary for the LLM + the viewer payload as structuredContent."""
    payload = result.get("_meta", {}).get("ui/init", {})
    return CallToolResult(
        content=[TextContent(type="text", text=_text(result))],
        structuredContent=payload or None,
    )


def create_server(config: BAMCPConfig | None = None) -> FastMCP:
    """Create and configure the BAMCP MCP server."""
    if config is None:
        config = BAMCPConfig.from_env()

    kwargs: dict = {
        "name": "bamcp",
        "instructions": SERVER_INSTRUCTIONS,
        "host": config.host,
        "port": config.port,
    }

    if config.auth_enabled:
        from .middleware.auth import BAMCPAuthProvider, build_auth_settings

        # Fail loud on an unusable auth config: with dynamic registration disabled (the default)
        # and no service token, no client could ever obtain a credential — the server would be
        # reachable but reject every request. Refuse to start instead of running a dead server.
        if not config.verify_token and not config.allow_dynamic_registration:
            raise ValueError(
                "Auth is enabled but there is no way to obtain a credential: set "
                "BAMCP_VERIFY_TOKEN (service token) or BAMCP_ALLOW_DYNAMIC_REGISTRATION=true."
            )

        kwargs["auth_server_provider"] = BAMCPAuthProvider(
            token_expiry=config.token_expiry,
            verify_token=config.verify_token,
            verify_scopes=config.required_scopes or [],
            allow_dynamic_registration=config.allow_dynamic_registration,
        )
        kwargs["auth"] = build_auth_settings(config)

    mcp = FastMCP(**kwargs)

    async def _run(tool: str, handler: _Handler, args: dict[str, Any]) -> dict:
        """Run a handler under the one error envelope: any failure surfaces as an MCP error
        result (isError=true) with a sanitized message — never raw pysam/OS text that could
        carry operator paths (see core/errors.py)."""
        try:
            return await handler(resolve_demo_alias(args), config)
        except Exception as e:  # noqa: BLE001 — every failure goes through the sanitizer
            raise sanitize_tool_exception(e, tool, config) from e

    # -- Tools ---------------------------------------------------------------
    # Thin wrappers delegate to the handlers in core/tools.py (and analysis/*) through _run.
    # FastMCP derives each JSON-Schema from the function signature; tool_specs.py mirrors those
    # for the eval harness, and a parity test keeps the two in lockstep.

    @mcp.tool(
        description=(
            "Detect and return candidate variants in a genomic region. "
            + REGION_FORMAT_NOTE
            + " Each variant carries artifact_risk, confidence, concerns and a plain-language "
            "interpretation (with a zygosity hint). "
            + DETECTION_DEFAULTS_NOTE
            + " "
            + AUTO_REFERENCE_NOTE
            + " variant_source: 'auto' (BAMCP candidates, VCF overlaid if given), "
            "'vcf' (the VCF is authoritative — BAMCP attaches read-level evidence at "
            "each site; requires vcf_path), or 'bamcp' (BAMCP candidates only)."
        )
    )
    async def get_variants(
        file_path: str,
        region: str,
        reference: str | None = None,
        vcf_path: str | None = None,
        variant_source: str = "auto",
        min_vaf: float | None = None,
        min_depth: int | None = None,
        min_alt_count: int | None = None,
    ) -> str:
        args = _args(
            file_path=file_path,
            region=region,
            reference=reference,
            vcf_path=vcf_path,
            variant_source=variant_source,
            min_vaf=min_vaf,
            min_depth=min_depth,
            min_alt_count=min_alt_count,
        )
        return _text(await _run("get_variants", handle_get_variants, args))

    @mcp.tool(
        description="Calculate depth of coverage statistics for a region. " + REGION_FORMAT_NOTE
    )
    async def get_coverage(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> str:
        args = _args(file_path=file_path, region=region, reference=reference)
        return _text(await _run("get_coverage", handle_get_coverage, args))

    @mcp.tool(
        description=(
            "Detect a BAM/CRAM's genome build (GRCh37 vs GRCh38) and contig naming style, and "
            "list its primary chromosomes. Use this first on new BAM files. The summary (build, "
            "style, suggested public reference URL, contig counts) comes first; alt/decoy/HLA/"
            "unplaced contigs are counted, not listed, unless include_all=true."
        ),
    )
    async def list_contigs(
        file_path: str,
        reference: str | None = None,
        include_all: bool = False,
    ) -> str:
        args = _args(file_path=file_path, reference=reference, include_all=include_all)
        return _text(await _run("list_contigs", handle_list_contigs, args))

    @mcp.tool(
        description=(
            "Jump to a specific genomic position and view surrounding reads. "
            "position is 1-based (VCF/samtools convention). If contig is omitted, the "
            "BAM's first contig is used."
        ),
        meta=_VIEWER_META,
    )
    async def jump_to(
        file_path: str,
        position: int,
        contig: str | None = None,
        window: int | None = None,
        reference: str | None = None,
        vcf_path: str | None = None,
        variant_source: str = "auto",
    ) -> CallToolResult:
        args = _args(
            file_path=file_path,
            position=position,
            contig=contig,
            window=window,
            reference=reference,
            vcf_path=vcf_path,
            variant_source=variant_source,
        )
        return _viewer_result(await _run("jump_to", handle_jump_to, args))

    # -- Phase 1: MCP Apps Tools --------------------------------------------

    @mcp.tool(
        description=(
            "Visualize aligned reads in a genomic region with interactive MCP Apps viewer. "
            + REGION_FORMAT_NOTE
        ),
        meta=_VIEWER_META,
    )
    async def visualize_region(
        file_path: str,
        region: str,
        reference: str | None = None,
        vcf_path: str | None = None,
        variant_source: str = "auto",
    ) -> CallToolResult:
        args = _args(
            file_path=file_path,
            region=region,
            reference=reference,
            vcf_path=vcf_path,
            variant_source=variant_source,
        )
        return _viewer_result(await _run("visualize_region", handle_visualize_region, args))

    @mcp.tool(
        description=(
            "Get a text summary of a genomic region "
            "(read count, coverage, candidate variants each with a plain-language "
            "interpretation) for LLM reasoning. " + REGION_FORMAT_NOTE + " " + AUTO_REFERENCE_NOTE
        ),
    )
    async def get_region_summary(
        file_path: str,
        region: str,
        reference: str | None = None,
        vcf_path: str | None = None,
        variant_source: str = "auto",
        min_alt_count: int | None = None,
    ) -> str:
        args = _args(
            file_path=file_path,
            region=region,
            reference=reference,
            vcf_path=vcf_path,
            variant_source=variant_source,
            min_alt_count=min_alt_count,
        )
        return _text(await _run("get_region_summary", handle_get_region_summary, args))

    # -- Phase 2: External Database Tools ------------------------------------

    @mcp.tool(
        description=(
            "Look up a variant in ClinVar for clinical significance, review status, "
            "and associated conditions. Positions are 1-based in `build` (GRCh37/hg19 or "
            "GRCh38/hg38; default: server build) — pass the BAM's build from list_contigs. "
            "Research-grade information, not for clinical use."
        ),
    )
    async def lookup_clinvar(
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        build: str | None = None,
    ) -> str:
        args = _args(chrom=chrom, pos=pos, ref=ref, alt=alt, build=build)
        return _text(await _run("lookup_clinvar", handle_lookup_clinvar, args))

    @mcp.tool(
        description=(
            "Look up a variant in gnomAD for population allele frequency data. Positions are "
            "1-based in `build` (GRCh37 queries gnomAD v2.1, GRCh38 queries v4; default: server "
            "build) — pass the BAM's build from list_contigs. "
            "Research-grade information, not for clinical use."
        ),
    )
    async def lookup_gnomad(
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        build: str | None = None,
    ) -> str:
        args = _args(chrom=chrom, pos=pos, ref=ref, alt=alt, build=build)
        return _text(await _run("lookup_gnomad", handle_lookup_gnomad, args))

    # -- Phase 3: Curation Tools ---------------------------------------------

    @mcp.tool(
        description=(
            "Get a detailed curation summary for a specific variant with artifact risk "
            "assessment, quality metrics, and recommendations for clinical interpretation. "
            "Use this for in-depth variant review by genetic curators. Set format='rubric' "
            "for machine-scorable JSON; include_clinical=true adds ClinVar/gnomAD context."
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
        format: str = "text",
        include_clinical: bool = False,
    ) -> str:
        args = _args(
            file_path=file_path,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            window=window,
            reference=reference,
            format=format,
            include_clinical=include_clinical,
        )
        return _text(
            await _run("get_variant_curation_summary", handle_get_variant_curation_summary, args)
        )

    @mcp.tool(
        description=(
            "Assemble a multi-source evidence package (BAM observation + ClinVar + gnomAD) for "
            "a variant and return an ACMG reasoning scaffold. Apply each applicable criterion, "
            "then give a classification (Pathogenic / Likely Pathogenic / VUS / Likely Benign / "
            "Benign) with confidence. Research-grade, not a clinical determination. Positions "
            "are 1-based."
        ),
    )
    async def classify_variant(
        file_path: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        gene: str | None = None,
        window: int = 50,
        reference: str | None = None,
    ) -> CallToolResult:
        args = _args(
            file_path=file_path,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gene=gene,
            window=window,
            reference=reference,
        )
        result = await _run("classify_variant", handle_classify_variant, args)
        return CallToolResult(
            content=[TextContent(type="text", text=_text(result))],
            structuredContent=result.get("structuredContent"),
        )

    @mcp.tool(
        description=(
            "Search for a gene by symbol and return its genomic coordinates (e.g. BRCA1, TP53). "
            "Coordinates are build-specific: pass the BAM's file_path so the build is "
            "auto-detected and the coordinates line up in that BAM (a GRCh38 locus in a GRCh37 "
            "BAM is ~0.5-2 Mb off, silently). Optionally set build (GRCh37/GRCh38) to override. "
            "The response states the coordinates' genome_build and warns on any mismatch. "
            "The returned region is 1-based inclusive, ready to pass to visualize_region."
        ),
    )
    async def search_gene(
        symbol: str,
        file_path: str | None = None,
        build: str | None = None,
    ) -> str:
        args = _args(symbol=symbol, file_path=file_path, build=build)
        return _text(await _run("search_gene", handle_search_gene, args))

    # -- Scanning Tools -------------------------------------------------------

    @mcp.tool(
        description=(
            "Scan an entire contig for candidate variants using fast coverage-based detection. "
            "Requires a reference genome (" + AUTO_REFERENCE_NOTE + ") Returns up to 500 "
            "candidate variants ranked by VAF. If contig is omitted, the BAM's first contig is "
            "scanned. " + DETECTION_DEFAULTS_NOTE
        ),
    )
    async def scan_variants(
        file_path: str,
        contig: str | None = None,
        reference: str | None = None,
        min_vaf: float | None = None,
        min_depth: int | None = None,
        min_alt_count: int | None = None,
    ) -> str:
        args = _args(
            file_path=file_path,
            contig=contig,
            reference=reference,
            min_vaf=min_vaf,
            min_depth=min_depth,
            min_alt_count=min_alt_count,
        )
        return _text(await _run("scan_variants", handle_scan_variants, args))

    # -- Cache Management Tool -----------------------------------------------

    @mcp.tool(
        description=(
            "Clear the server's downloaded BAM index files and region cache (single-client "
            "stdio servers only; shared HTTP servers expire their cache by TTL)"
        )
    )
    async def cleanup_cache() -> str:
        return _text(await _run("cleanup_cache", handle_cleanup_cache, {}))

    # -- Resources -----------------------------------------------------------

    @mcp.resource(
        "ui://bamcp/viewer",
        name="BAMCP Alignment Viewer",
        description="Interactive BAM/CRAM alignment visualization",
        mime_type="text/html;profile=mcp-app",
    )
    def viewer() -> str:
        return get_viewer_html()

    register_demo_surface(mcp)

    # -- Deploy identity -----------------------------------------------------
    # Greenlight's verify gate compares the commit the image was built from (GREENLIGHT_SHA, baked
    # in at docker build) against the sha CI just shipped, so a stale container can't pass as the
    # new build. Served unauthenticated (a commit sha is public in this repo) at the root and under
    # the MCP path, because the verify base URL is <host>/mcp. A build without the arg reports null.

    async def deploy_version(request: Request) -> JSONResponse:
        return JSONResponse(
            {"sha": os.environ.get("GREENLIGHT_SHA") or None},
            headers={"Cache-Control": "no-store"},
        )

    mcp.custom_route("/__version", methods=["GET"])(deploy_version)
    mcp.custom_route("/mcp/__version", methods=["GET"])(deploy_version)

    return mcp

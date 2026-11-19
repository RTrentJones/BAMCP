"""MCP server setup for BAMCP using FastMCP."""

from __future__ import annotations

import os
from collections.abc import Awaitable, Callable
from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.server.fastmcp.exceptions import ToolError
from mcp.types import CallToolResult, TextContent
from pydantic import ValidationError
from starlette.requests import Request
from starlette.responses import JSONResponse

from . import __version__
from .analysis.acmg import handle_classify_variant
from .analysis.curation import handle_get_variant_curation_summary
from .config import BAMCPConfig
from .core.errors import format_argument_error, sanitize_tool_exception
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
from .tool_specs import TOOL_SPECS

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


def _install_argument_error_mapping(mcp: FastMCP) -> None:
    """Map schema-validation failures to a short message naming the expected parameters.

    Argument validation runs inside FastMCP (before ``_run``), which renders the raw pydantic
    dump (with an errors.pydantic.dev URL) as the tool error. Wrap the tool manager's
    ``call_tool`` so that one failure mode gets the same client-friendly envelope as the rest.
    """
    manager = mcp._tool_manager
    original = manager.call_tool

    async def call_tool(
        name: str, arguments: dict[str, Any], context: Any = None, convert_result: bool = False
    ) -> Any:
        try:
            return await original(name, arguments, context=context, convert_result=convert_result)
        except ToolError as e:
            cause = e.__cause__
            tool = manager.get_tool(name)
            if isinstance(cause, ValidationError) and tool is not None:
                raise ToolError(
                    format_argument_error(name, cause, tool.parameters, arguments)
                ) from None
            raise

    manager.call_tool = call_tool  # type: ignore[method-assign]


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
    # serverInfo.version: BAMCP's own version, not the mcp SDK's (FastMCP's default fallback).
    mcp._mcp_server.version = __version__

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
    # for the eval harness, and a parity test keeps the two in lockstep. Descriptions have one
    # source: TOOL_SPECS (the parity test asserts the registered text equals it).

    @mcp.tool(description=TOOL_SPECS["get_variants"]["description"])
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

    @mcp.tool(description=TOOL_SPECS["get_coverage"]["description"])
    async def get_coverage(
        file_path: str,
        region: str,
        reference: str | None = None,
    ) -> str:
        args = _args(file_path=file_path, region=region, reference=reference)
        return _text(await _run("get_coverage", handle_get_coverage, args))

    @mcp.tool(description=TOOL_SPECS["list_contigs"]["description"])
    async def list_contigs(
        file_path: str,
        reference: str | None = None,
        include_all: bool = False,
    ) -> str:
        args = _args(file_path=file_path, reference=reference, include_all=include_all)
        return _text(await _run("list_contigs", handle_list_contigs, args))

    @mcp.tool(description=TOOL_SPECS["jump_to"]["description"], meta=_VIEWER_META)
    async def jump_to(
        file_path: str,
        contig: str,
        position: int,
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

    @mcp.tool(description=TOOL_SPECS["visualize_region"]["description"], meta=_VIEWER_META)
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

    @mcp.tool(description=TOOL_SPECS["get_region_summary"]["description"])
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

    @mcp.tool(description=TOOL_SPECS["lookup_clinvar"]["description"])
    async def lookup_clinvar(
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        build: str | None = None,
    ) -> str:
        args = _args(chrom=chrom, pos=pos, ref=ref, alt=alt, build=build)
        return _text(await _run("lookup_clinvar", handle_lookup_clinvar, args))

    @mcp.tool(description=TOOL_SPECS["lookup_gnomad"]["description"])
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

    @mcp.tool(description=TOOL_SPECS["get_variant_curation_summary"]["description"])
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

    @mcp.tool(description=TOOL_SPECS["classify_variant"]["description"])
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

    @mcp.tool(description=TOOL_SPECS["search_gene"]["description"])
    async def search_gene(
        symbol: str,
        file_path: str | None = None,
        build: str | None = None,
    ) -> str:
        args = _args(symbol=symbol, file_path=file_path, build=build)
        return _text(await _run("search_gene", handle_search_gene, args))

    # -- Scanning Tools -------------------------------------------------------

    @mcp.tool(description=TOOL_SPECS["scan_variants"]["description"])
    async def scan_variants(
        file_path: str,
        contig: str | None = None,
        region: str | None = None,
        reference: str | None = None,
        min_vaf: float | None = None,
        min_depth: int | None = None,
        min_alt_count: int | None = None,
    ) -> str:
        args = _args(
            file_path=file_path,
            contig=contig,
            region=region,
            reference=reference,
            min_vaf=min_vaf,
            min_depth=min_depth,
            min_alt_count=min_alt_count,
        )
        return _text(await _run("scan_variants", handle_scan_variants, args))

    # -- Cache Management Tool -----------------------------------------------

    @mcp.tool(description=TOOL_SPECS["cleanup_cache"]["description"])
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
    _install_argument_error_mapping(mcp)

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

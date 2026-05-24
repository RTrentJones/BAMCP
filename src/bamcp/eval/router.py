"""Tool routing for the eval harness.

A router takes a tool name + arguments and dispatches to the correct BAMCP
handler. Two implementations:

- ``InProcessRouter``: imports the handlers and calls them directly. Used by
  the default ``--router in-process`` mode in Phase 1.
- ``MCPStdioRouter``: launches ``python -m bamcp`` as a stdio subprocess and
  routes calls through the MCP protocol. Stub for Phase 2.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Protocol

from ..analysis.curation import handle_get_variant_curation_summary
from ..config import BAMCPConfig
from ..core.tools import (
    get_gene_client,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_scan_variants,
    handle_visualize_region,
)


@dataclass
class RouterResult:
    """Outcome of a tool call routed through the harness."""

    ok: bool
    text: str  # Tool's content[0].text — what the LLM sees
    error: str | None = None
    # Surfaced when a tool's response carried ``_meta.ui/init`` (the MCP Apps
    # viewer payload). The visual-eval runner uses this to drive screenshot
    # capture; text-only runs ignore it.
    ui_payload: dict | None = None


class ToolRouter(Protocol):
    """Dispatches a tool name + arguments to the right BAMCP handler."""

    def list_tools(self) -> list[str]: ...

    async def call(self, name: str, arguments: dict[str, Any]) -> RouterResult: ...


# Map MCP-exposed tool names to their async handler functions.
_HANDLERS: dict[str, Callable] = {
    "get_variants": handle_get_variants,
    "get_coverage": handle_get_coverage,
    "list_contigs": handle_list_contigs,
    "jump_to": handle_jump_to,
    "visualize_region": handle_visualize_region,
    "get_region_summary": handle_get_region_summary,
    "lookup_clinvar": handle_lookup_clinvar,
    "lookup_gnomad": handle_lookup_gnomad,
    "scan_variants": handle_scan_variants,
    "get_variant_curation_summary": handle_get_variant_curation_summary,
}


class InProcessRouter:
    """Routes tool calls by importing and invoking handlers directly.

    Useful for fast deterministic evals and for tests. Skips the MCP
    serialization layer — exercise that path in Phase 2 with MCPStdioRouter.
    """

    def __init__(self, config: BAMCPConfig | None = None) -> None:
        self.config = config or BAMCPConfig()

    def list_tools(self) -> list[str]:
        names = list(_HANDLERS.keys())
        names.append("search_gene")
        names.append("cleanup_cache")
        return names

    async def call(self, name: str, arguments: dict[str, Any]) -> RouterResult:
        if name == "search_gene":
            return await self._call_search_gene(arguments)
        if name == "cleanup_cache":
            from ..core.tools import get_cache

            removed = get_cache(self.config).cleanup_session()
            return RouterResult(ok=True, text=f"Removed {removed} cached index files")
        handler = _HANDLERS.get(name)
        if handler is None:
            return RouterResult(ok=False, text="", error=f"Unknown tool: {name}")
        try:
            result = await handler(arguments, self.config)
        except Exception as e:  # noqa: BLE001 — surface failure to LLM
            return RouterResult(ok=False, text="", error=f"{type(e).__name__}: {e}")
        content = result.get("content") if isinstance(result, dict) else None
        if isinstance(content, list) and content:
            text = content[0].get("text", "") if isinstance(content[0], dict) else ""
        else:
            text = json.dumps(result, default=str)
        ui_payload = None
        if isinstance(result, dict):
            meta = result.get("_meta")
            if isinstance(meta, dict):
                candidate = meta.get("ui/init")
                if isinstance(candidate, dict):
                    ui_payload = candidate
        return RouterResult(ok=True, text=str(text), ui_payload=ui_payload)

    async def _call_search_gene(self, arguments: dict[str, Any]) -> RouterResult:
        symbol = arguments.get("symbol", "")
        if not symbol:
            return RouterResult(ok=False, text="", error="search_gene requires 'symbol'")
        try:
            client = get_gene_client(self.config)
            result = await client.search(symbol)
        except Exception as e:  # noqa: BLE001
            return RouterResult(ok=False, text="", error=f"{type(e).__name__}: {e}")
        if result is None:
            return RouterResult(ok=True, text=json.dumps({"error": f"Gene {symbol!r} not found"}))
        return RouterResult(
            ok=True,
            text=json.dumps(
                {
                    "symbol": result.symbol,
                    "name": result.name,
                    "region": f"{result.chrom}:{result.start}-{result.end}",
                    "strand": result.strand,
                }
            ),
        )


def tool_descriptors(router: ToolRouter) -> list:
    """Build provider-agnostic tool descriptors from the router's tool list.

    Each tool exposes a permissive JSON schema (object with arbitrary string
    properties) since the eval harness doesn't have access to FastMCP's
    inferred schemas. Real providers handle this gracefully.
    """
    from .providers import ToolDescriptor

    out: list[ToolDescriptor] = []
    for name in router.list_tools():
        out.append(
            ToolDescriptor(
                name=name,
                description=f"BAMCP tool: {name}",
                input_schema={
                    "type": "object",
                    "additionalProperties": True,
                },
            )
        )
    return out

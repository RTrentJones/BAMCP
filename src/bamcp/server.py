"""MCP server setup for BAMCP."""

from mcp.server import Server
from mcp.types import Resource, TextContent, Tool

from .config import BAMCPConfig
from .resources import get_viewer_html
from .tools import (
    handle_browse_region,
    handle_get_coverage,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
)

TOOLS = [
    Tool(
        name="browse_region",
        description="View aligned reads in a genomic region with interactive visualization",
        inputSchema={
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                "region": {
                    "type": "string",
                    "description": "Genomic region (e.g., chr1:1000-2000)",
                },
                "reference": {
                    "type": "string",
                    "description": "Reference FASTA path (required for CRAM)",
                },
            },
            "required": ["file_path", "region"],
        },
    ),
    Tool(
        name="get_variants",
        description="Detect and return variants in a genomic region",
        inputSchema={
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                "region": {
                    "type": "string",
                    "description": "Genomic region (e.g., chr1:1000-2000)",
                },
                "reference": {"type": "string", "description": "Reference FASTA path"},
                "min_vaf": {
                    "type": "number",
                    "description": "Minimum variant allele frequency",
                },
                "min_depth": {"type": "integer", "description": "Minimum read depth"},
            },
            "required": ["file_path", "region"],
        },
    ),
    Tool(
        name="get_coverage",
        description="Calculate depth of coverage statistics for a region",
        inputSchema={
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                "region": {
                    "type": "string",
                    "description": "Genomic region (e.g., chr1:1000-2000)",
                },
                "reference": {"type": "string", "description": "Reference FASTA path"},
            },
            "required": ["file_path", "region"],
        },
    ),
    Tool(
        name="list_contigs",
        description="List chromosomes/contigs in a BAM/CRAM file header",
        inputSchema={
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                "reference": {"type": "string", "description": "Reference FASTA path"},
            },
            "required": ["file_path"],
        },
    ),
    Tool(
        name="jump_to",
        description="Jump to a specific genomic position and view surrounding reads",
        inputSchema={
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                "position": {"type": "integer", "description": "Genomic position to center on"},
                "contig": {
                    "type": "string",
                    "description": "Chromosome/contig name (default: chr1)",
                },
                "window": {
                    "type": "integer",
                    "description": "Window size in bp around position",
                },
                "reference": {"type": "string", "description": "Reference FASTA path"},
            },
            "required": ["file_path", "position"],
        },
    ),
]

TOOL_HANDLERS = {
    "browse_region": handle_browse_region,
    "get_variants": handle_get_variants,
    "get_coverage": handle_get_coverage,
    "list_contigs": handle_list_contigs,
    "jump_to": handle_jump_to,
}


def create_server(config: BAMCPConfig | None = None) -> Server:
    """Create and configure the BAMCP MCP server."""
    if config is None:
        config = BAMCPConfig.from_env()

    app = Server("bamcp")

    @app.list_tools()
    async def list_tools() -> list[Tool]:
        return TOOLS

    @app.call_tool()
    async def call_tool(name: str, arguments: dict) -> list[TextContent]:
        handler = TOOL_HANDLERS.get(name)
        if handler is None:
            raise ValueError(f"Unknown tool: {name}")

        result = await handler(arguments, config)
        return [TextContent(type="text", text=result["content"][0]["text"])]

    @app.list_resources()
    async def list_resources() -> list[Resource]:
        return [
            Resource(
                uri="ui://bamcp/viewer",
                name="BAMCP Alignment Viewer",
                mimeType="text/html+mcp",
                description="Interactive BAM/CRAM alignment visualization",
            )
        ]

    @app.read_resource()
    async def read_resource(uri: str) -> str:
        if str(uri) == "ui://bamcp/viewer":
            return get_viewer_html()
        raise ValueError(f"Unknown resource: {uri}")

    return app

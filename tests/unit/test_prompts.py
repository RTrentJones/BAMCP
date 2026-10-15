"""Newcomer MCP surface: instructions, resources, prompts, and the demo alias through _run."""

from __future__ import annotations

import json

import pytest
from mcp.shared.memory import create_connected_server_and_client_session

from bamcp.config import BAMCPConfig
from bamcp.demos import DEMO_DATASETS
from bamcp.prompts import SERVER_INSTRUCTIONS
from bamcp.server import create_server

EXPECTED_TOOLS = {
    "get_variants",
    "get_coverage",
    "list_contigs",
    "jump_to",
    "visualize_region",
    "get_region_summary",
    "lookup_clinvar",
    "lookup_gnomad",
    "get_variant_curation_summary",
    "classify_variant",
    "search_gene",
    "scan_variants",
    "cleanup_cache",
}


def _client(config: BAMCPConfig | None = None):
    server = create_server(config or BAMCPConfig())
    return create_connected_server_and_client_session(server._mcp_server)


@pytest.mark.unit
def test_instructions_cover_newcomer_essentials():
    text = SERVER_INSTRUCTIONS
    assert "demo:HG002" in text and "demo:HG00096" in text
    assert "list_contigs first" in text
    assert "1-based" in text
    assert "not for clinical use" in text
    assert create_server(BAMCPConfig()).instructions == SERVER_INSTRUCTIONS


@pytest.mark.integration
async def test_tool_list_is_unchanged():
    async with _client() as client:
        tools = await client.list_tools()
    assert {t.name for t in tools.tools} == EXPECTED_TOOLS
    assert len(tools.tools) == 13


@pytest.mark.integration
async def test_resources_listed_and_readable():
    async with _client() as client:
        listed = await client.list_resources()
        uris = {str(r.uri) for r in listed.resources}
        assert {"bamcp://demos", "bamcp://glossary", "bamcp://guide"} <= uris
        assert "ui://bamcp/viewer" in uris

        demos = await client.read_resource("bamcp://demos")
        data = json.loads(demos.contents[0].text)
        assert [d["alias"] for d in data["datasets"]] == ["demo:HG002", "demo:HG00096"]

        glossary = await client.read_resource("bamcp://glossary")
        assert "VAF" in glossary.contents[0].text
        assert "Do not conclude" in glossary.contents[0].text

        guide = await client.read_resource("bamcp://guide")
        assert "list_contigs" in guide.contents[0].text
        assert "demo:HG002" in guide.contents[0].text


@pytest.mark.integration
async def test_prompts_listed_and_rendered():
    async with _client() as client:
        listed = await client.list_prompts()
        names = {p.name for p in listed.prompts}
        assert {"explore_demo", "check_variant", "explain_terms"} <= names

        tour = await client.get_prompt("explore_demo", {"dataset": "demo:HG00096"})
        text = tour.messages[0].content.text
        assert "demo:HG00096" in text and "GRCh37" in text
        assert DEMO_DATASETS["HG00096"].reference_url in text

        default_tour = await client.get_prompt("explore_demo", {})
        text = default_tour.messages[0].content.text
        assert "demo:HG002" in text
        # Gene-scale loci are surveyed with get_coverage, never read whole by get_variants /
        # visualize_region (81 kb: ~1 min and a timeout); small loci go to the read tools.
        small, _, genes = text.partition("Gene-scale loci")
        assert "chr17:43044295-43125483" in genes and "chr17:43044295-43125483" not in small
        assert "get_coverage" in genes and "do NOT call get_variants" in genes
        assert "chr17:43099864-43099964" in small
        assert "even if get_variants did not report it" in text

        check = await client.get_prompt(
            "check_variant", {"file_path": "demo:HG002", "region": "chr17:43099914"}
        )
        text = check.messages[0].content.text
        assert "chr17:43099914-43099914" in text and "list_contigs" in text

        terms = await client.get_prompt("explain_terms", {"terms": "VAF, MAPQ"})
        text = terms.messages[0].content.text
        assert "VAF" in text and "MAPQ" in text and "Homopolymer" not in text

        everything = await client.get_prompt("explain_terms", {"terms": "zzz-unknown"})
        assert "Homopolymer" in everything.messages[0].content.text


@pytest.mark.integration
async def test_unknown_demo_alias_is_clear_tool_error():
    async with _client() as client:
        result = await client.call_tool("list_contigs", {"file_path": "demo:NOPE"})
    assert result.isError
    text = result.content[0].text
    assert "Unknown demo dataset" in text and "demo:HG002" in text


@pytest.mark.integration
async def test_demo_alias_resolves_before_validation():
    """The resolved URL goes through normal remote-URL validation (hermetic: allow-list miss)."""
    config = BAMCPConfig(allow_remote_files=True, allowed_remote_hosts=["example.org"])
    async with _client(config) as client:
        result = await client.call_tool("list_contigs", {"file_path": "demo:HG002"})
    assert result.isError
    assert "ftp-trace.ncbi.nlm.nih.gov" in result.content[0].text


@pytest.mark.network
async def test_demo_hg002_real_het_snv_live():
    demo = DEMO_DATASETS["HG002"]
    config = BAMCPConfig(allow_remote_files=True)
    async with _client(config) as client:
        result = await client.call_tool(
            "get_variants",
            {
                "file_path": demo.alias,
                "region": "chr17:43099900-43099930",
                "reference": demo.reference_url,
            },
        )
    assert not result.isError, result.content[0].text
    variants = json.loads(result.content[0].text)["variants"]
    assert any(v["position"] == 43099914 and v["alt"] == "A" for v in variants)

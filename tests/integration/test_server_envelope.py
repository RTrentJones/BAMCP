"""End-to-end (in-memory MCP client) checks of the tool error envelope and server wiring."""

from __future__ import annotations

import json
import os

import pytest
from mcp.shared.memory import create_connected_server_and_client_session

from bamcp.config import BAMCPConfig
from bamcp.core.tools import get_services
from bamcp.server import create_server

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "fixtures")
SMALL_BAM = os.path.join(FIXTURES, "small.bam")
REF = os.path.join(FIXTURES, "ref.fa")


async def _call(config: BAMCPConfig, name: str, args: dict):
    server = create_server(config)
    async with create_connected_server_and_client_session(server._mcp_server) as client:
        return await client.call_tool(name, args)


def _text(result) -> str:
    return result.content[0].text


@pytest.mark.integration
async def test_missing_operator_reference_is_error_without_path():
    result = await _call(
        BAMCPConfig(reference="/opt/secret/x.fa"),
        "get_variants",
        {"file_path": SMALL_BAM, "region": "chr1:1-200"},
    )
    assert result.isError
    assert "/opt/secret" not in _text(result)


@pytest.mark.integration
async def test_jump_to_unknown_contig_is_clean_error():
    result = await _call(
        BAMCPConfig(), "jump_to", {"file_path": SMALL_BAM, "position": 100, "contig": "chrZZ"}
    )
    assert result.isError
    assert "KeyError" not in _text(result)
    assert "chrZZ" in _text(result)


@pytest.mark.integration
async def test_jump_to_validates_region_like_visualize_region():
    result = await _call(
        BAMCPConfig(), "jump_to", {"file_path": SMALL_BAM, "position": 100, "contig": "chr 1"}
    )
    assert result.isError
    assert "Invalid region format" in _text(result)


@pytest.mark.integration
async def test_jump_to_defaults_to_first_bam_contig():
    result = await _call(BAMCPConfig(), "jump_to", {"file_path": SMALL_BAM, "position": 100})
    assert not result.isError
    assert result.structuredContent["contig"] == "chr1"


@pytest.mark.integration
async def test_scan_variants_without_reference_is_error():
    result = await _call(BAMCPConfig(), "scan_variants", {"file_path": SMALL_BAM})
    assert result.isError
    assert "Reference genome required" in _text(result)


@pytest.mark.integration
async def test_scan_variants_without_contig_scans_first_contig():
    result = await _call(BAMCPConfig(reference=REF), "scan_variants", {"file_path": SMALL_BAM})
    assert not result.isError
    assert json.loads(_text(result))["contig"] == "chr1"


@pytest.mark.integration
async def test_curation_summary_exposes_include_clinical(monkeypatch):
    seen: dict = {}

    async def fake(args, config):
        seen.update(args)
        return {"content": [{"type": "text", "text": "ok"}]}

    monkeypatch.setattr("bamcp.server.handle_get_variant_curation_summary", fake)
    result = await _call(
        BAMCPConfig(),
        "get_variant_curation_summary",
        {
            "file_path": SMALL_BAM,
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "C",
            "include_clinical": True,
        },
    )
    assert not result.isError
    assert seen["include_clinical"] is True


@pytest.mark.integration
async def test_unset_optionals_are_not_forwarded(monkeypatch):
    seen: dict = {}

    async def fake(args, config):
        seen.update(args)
        return {"content": [{"type": "text", "text": "ok"}]}

    monkeypatch.setattr("bamcp.server.handle_get_variants", fake)
    await _call(BAMCPConfig(), "get_variants", {"file_path": SMALL_BAM, "region": "chr1:1-2"})
    assert seen == {"file_path": SMALL_BAM, "region": "chr1:1-2", "variant_source": "auto"}


class _Closable:
    def __init__(self):
        self.closed = False

    async def close(self):
        self.closed = True


@pytest.mark.integration
async def test_cleanup_cache_never_closes_shared_clients(tmp_path):
    config = BAMCPConfig(cache_dir=str(tmp_path))
    services = get_services(config)
    client = _Closable()
    services._clinvar = client  # type: ignore[assignment]
    services.region_cache[("k",)] = (0.0, object(), 1)  # type: ignore[assignment]
    services.cache.cache_dir.mkdir(parents=True, exist_ok=True)
    (services.cache.cache_dir / "x.bai").write_bytes(b"x")

    result = await _call(config, "cleanup_cache", {})

    assert not result.isError
    assert "Removed 1 cached index files" in _text(result)
    assert services._clinvar is client and not client.closed
    assert not services.region_cache


@pytest.mark.integration
async def test_cleanup_cache_is_noop_on_shared_http_server(tmp_path):
    config = BAMCPConfig(cache_dir=str(tmp_path), transport="streamable-http")
    services = get_services(config)
    services.cache.cache_dir.mkdir(parents=True, exist_ok=True)
    index = services.cache.cache_dir / "x.bai"
    index.write_bytes(b"x")

    result = await _call(config, "cleanup_cache", {})

    assert not result.isError
    assert "disabled on shared" in _text(result)
    assert index.exists()
    assert "this session" not in _text(result)

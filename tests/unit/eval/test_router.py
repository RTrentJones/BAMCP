"""Unit tests for the in-process tool router."""

from __future__ import annotations

import json

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval.router import InProcessRouter, tool_descriptors


@pytest.fixture
def cfg(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path)


@pytest.mark.unit
def test_list_tools_includes_all_handlers():
    router = InProcessRouter()
    tools = router.list_tools()
    assert "get_variants" in tools
    assert "get_variant_curation_summary" in tools
    assert "search_gene" in tools
    assert "cleanup_cache" in tools


@pytest.mark.unit
async def test_unknown_tool_returns_error_result():
    router = InProcessRouter()
    result = await router.call("does_not_exist", {})
    assert result.ok is False
    assert "Unknown tool" in (result.error or "")


@pytest.mark.unit
async def test_known_tool_returns_text(small_bam_path, cfg):
    router = InProcessRouter(config=cfg)
    result = await router.call(
        "get_coverage", {"file_path": small_bam_path, "region": "chr1:90-200"}
    )
    assert result.ok is True
    # Should be JSON-parseable
    parsed = json.loads(result.text)
    assert "mean" in parsed
    assert "region" in parsed


@pytest.mark.unit
async def test_search_gene_missing_symbol():
    router = InProcessRouter()
    result = await router.call("search_gene", {})
    assert result.ok is False
    assert "symbol" in (result.error or "").lower()


@pytest.mark.unit
def test_tool_descriptors_built_from_router_list():
    router = InProcessRouter()
    descs = tool_descriptors(router)
    names = [d.name for d in descs]
    assert "get_variants" in names
    # Every descriptor exposes a schema dict.
    for d in descs:
        assert isinstance(d.input_schema, dict)
        assert d.input_schema["type"] == "object"


@pytest.mark.unit
async def test_cleanup_cache_path():
    """The router's cleanup_cache shortcut returns a count of removed files."""
    router = InProcessRouter()
    result = await router.call("cleanup_cache", {})
    assert result.ok is True
    assert "cached index files" in result.text


@pytest.mark.unit
async def test_search_gene_caught_exception(monkeypatch):
    """search_gene wraps client failures as RouterResult errors."""
    from bamcp.eval import router as router_mod

    class BoomClient:
        async def search(self, symbol):
            raise RuntimeError("network down")

    monkeypatch.setattr(router_mod, "get_gene_client", lambda config: BoomClient())
    result = await InProcessRouter().call("search_gene", {"symbol": "TP53"})
    assert result.ok is False
    assert "network down" in (result.error or "")


@pytest.mark.unit
async def test_search_gene_not_found(monkeypatch):
    from bamcp.eval import router as router_mod

    class EmptyClient:
        async def search(self, symbol):
            return None

    monkeypatch.setattr(router_mod, "get_gene_client", lambda config: EmptyClient())
    result = await InProcessRouter().call("search_gene", {"symbol": "BOGUS"})
    assert result.ok is True
    assert "not found" in result.text


@pytest.mark.unit
async def test_search_gene_returns_region(monkeypatch):
    from bamcp.eval import router as router_mod

    class Hit:
        symbol = "TP53"
        name = "tumor protein p53"
        chrom = "chr17"
        start = 7670000
        end = 7680000
        strand = "-"

    class HitClient:
        async def search(self, symbol):
            return Hit()

    monkeypatch.setattr(router_mod, "get_gene_client", lambda config: HitClient())
    result = await InProcessRouter().call("search_gene", {"symbol": "TP53"})
    assert result.ok is True
    parsed = json.loads(result.text)
    assert parsed["symbol"] == "TP53"
    assert parsed["region"] == "chr17:7670000-7680000"


@pytest.mark.unit
async def test_handler_exception_surfaces_as_error_result():
    """Internal handler crashes become RouterResult errors instead of raising."""
    from bamcp.eval import router as router_mod

    async def boom(args, config):  # noqa: ANN001
        raise RuntimeError("handler bug")

    # Patch one entry in the handler table.
    monkeypatch_handlers = dict(router_mod._HANDLERS)
    monkeypatch_handlers["get_variants"] = boom
    router_mod._HANDLERS, original = monkeypatch_handlers, router_mod._HANDLERS
    try:
        router = InProcessRouter()
        result = await router.call("get_variants", {"file_path": "x", "region": "chr1:1-2"})
        assert result.ok is False
        assert "handler bug" in (result.error or "")
    finally:
        router_mod._HANDLERS = original

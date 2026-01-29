"""Unit tests for bamcp.server module."""

import pytest

from bamcp.config import BAMCPConfig
from bamcp.server import TOOLS, TOOL_HANDLERS, create_server


class TestToolDefinitions:
    """Tests for the tool definitions."""

    @pytest.mark.unit
    def test_tool_count(self):
        assert len(TOOLS) == 4

    @pytest.mark.unit
    def test_tool_names(self):
        names = {t.name for t in TOOLS}
        assert names == {"browse_region", "get_variants", "get_coverage", "list_contigs"}

    @pytest.mark.unit
    def test_browse_region_schema(self):
        tool = next(t for t in TOOLS if t.name == "browse_region")
        schema = tool.inputSchema
        assert "file_path" in schema["properties"]
        assert "region" in schema["properties"]
        assert "reference" in schema["properties"]
        assert set(schema["required"]) == {"file_path", "region"}

    @pytest.mark.unit
    def test_get_variants_schema(self):
        tool = next(t for t in TOOLS if t.name == "get_variants")
        schema = tool.inputSchema
        assert "file_path" in schema["properties"]
        assert "region" in schema["properties"]
        assert "min_vaf" in schema["properties"]
        assert "min_depth" in schema["properties"]

    @pytest.mark.unit
    def test_get_coverage_schema(self):
        tool = next(t for t in TOOLS if t.name == "get_coverage")
        schema = tool.inputSchema
        assert set(schema["required"]) == {"file_path", "region"}

    @pytest.mark.unit
    def test_list_contigs_schema(self):
        tool = next(t for t in TOOLS if t.name == "list_contigs")
        schema = tool.inputSchema
        assert schema["required"] == ["file_path"]

    @pytest.mark.unit
    def test_all_tools_have_descriptions(self):
        for tool in TOOLS:
            assert tool.description, f"Tool {tool.name} missing description"
            assert len(tool.description) > 10

    @pytest.mark.unit
    def test_handler_map_matches_tools(self):
        tool_names = {t.name for t in TOOLS}
        handler_names = set(TOOL_HANDLERS.keys())
        assert tool_names == handler_names


class TestCreateServer:
    """Tests for the create_server function."""

    @pytest.mark.unit
    def test_creates_server(self):
        config = BAMCPConfig()
        server = create_server(config)
        assert server is not None
        assert server.name == "bamcp"

    @pytest.mark.unit
    def test_creates_server_default_config(self, monkeypatch):
        """Should create server with default config from env."""
        import os
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        server = create_server()
        assert server is not None

    @pytest.mark.unit
    def test_server_has_name(self):
        server = create_server(BAMCPConfig())
        assert server.name == "bamcp"

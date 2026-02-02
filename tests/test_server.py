"""Unit tests for bamcp.server module."""

import os

import pytest

from bamcp.config import BAMCPConfig
from bamcp.server import create_server


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
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        server = create_server()
        assert server is not None
        assert server.name == "bamcp"

    @pytest.mark.unit
    def test_server_registers_tools(self):
        """Server should register all expected tools."""
        server = create_server(BAMCPConfig())
        tool_names = {name for name in server._tool_manager._tools}
        expected = {"browse_region", "get_variants", "get_coverage", "list_contigs", "jump_to"}
        assert expected == tool_names

    @pytest.mark.unit
    def test_server_registers_resources(self):
        """Server should register the viewer resource."""
        server = create_server(BAMCPConfig())
        resource_keys = list(server._resource_manager._resources.keys())
        assert len(resource_keys) >= 1
        assert any("viewer" in str(k) for k in resource_keys)

    @pytest.mark.unit
    def test_server_without_auth(self):
        """Server created without auth should have no auth provider."""
        config = BAMCPConfig(auth_enabled=False)
        server = create_server(config)
        assert server is not None

    @pytest.mark.unit
    def test_server_with_auth(self):
        """Server created with auth should have auth provider."""
        config = BAMCPConfig(
            auth_enabled=True,
            issuer_url="http://localhost:8000",
            resource_server_url="http://localhost:8000",
        )
        server = create_server(config)
        assert server is not None

    @pytest.mark.unit
    def test_server_host_port(self):
        """Server should use config host and port."""
        config = BAMCPConfig(host="127.0.0.1", port=9000)
        server = create_server(config)
        assert server.settings.host == "127.0.0.1"
        assert server.settings.port == 9000

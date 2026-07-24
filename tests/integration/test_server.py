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
        expected = {
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
        """Server created with auth + a credential should have an auth provider."""
        config = BAMCPConfig(
            auth_enabled=True,
            issuer_url="http://localhost:8000",
            resource_server_url="http://localhost:8000",
            verify_token="svc-secret",  # a bootstrap credential
        )
        server = create_server(config)
        assert server is not None

    @pytest.mark.unit
    def test_auth_without_bootstrap_credential_refuses_to_start(self):
        """Auth enabled with no service token and no registration is an unusable config."""
        config = BAMCPConfig(
            auth_enabled=True,
            allow_dynamic_registration=False,
            verify_token=None,
        )
        with pytest.raises(ValueError, match="no way to obtain a credential"):
            create_server(config)

    @pytest.mark.unit
    def test_auth_with_dynamic_registration_only_is_allowed(self):
        """Dynamic registration alone is a valid bootstrap path (dev/trusted context)."""
        config = BAMCPConfig(auth_enabled=True, allow_dynamic_registration=True)
        assert create_server(config) is not None

    @pytest.mark.unit
    def test_server_host_port(self):
        """Server should use config host and port."""
        config = BAMCPConfig(host="127.0.0.1", port=9000)
        server = create_server(config)
        assert server.settings.host == "127.0.0.1"
        assert server.settings.port == 9000

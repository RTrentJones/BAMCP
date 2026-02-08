"""Unit tests for bamcp.config module."""

import os

import pytest

from bamcp.config import BAMCPConfig
from bamcp.constants import (
    DEFAULT_CACHE_TTL_SECONDS,
    DEFAULT_HOST,
    DEFAULT_MAX_READS,
    DEFAULT_MIN_DEPTH,
    DEFAULT_MIN_MAPQ,
    DEFAULT_MIN_VAF,
    DEFAULT_PORT,
    DEFAULT_TOKEN_EXPIRY_SECONDS,
    DEFAULT_WINDOW_SIZE,
)


class TestBAMCPConfig:
    """Tests for BAMCPConfig dataclass."""

    @pytest.mark.unit
    def test_default_values(self):
        """Config defaults should match documented values."""
        config = BAMCPConfig()
        assert config.reference is None
        assert config.max_reads == DEFAULT_MAX_READS
        assert config.default_window == DEFAULT_WINDOW_SIZE
        assert config.min_vaf == DEFAULT_MIN_VAF
        assert config.min_depth == DEFAULT_MIN_DEPTH
        assert config.min_mapq == DEFAULT_MIN_MAPQ
        assert config.transport == "stdio"
        assert config.host == DEFAULT_HOST
        assert config.port == DEFAULT_PORT
        assert config.auth_enabled is False
        assert config.token_expiry == DEFAULT_TOKEN_EXPIRY_SECONDS

    @pytest.mark.unit
    def test_custom_values(self):
        """Config should accept custom values."""
        config = BAMCPConfig(
            reference="/path/to/ref.fa",
            max_reads=5000,
            default_window=250,
            min_vaf=0.05,
            min_depth=20,
            min_mapq=10,
        )
        assert config.reference == "/path/to/ref.fa"
        assert config.max_reads == 5000
        assert config.default_window == 250
        assert config.min_vaf == 0.05
        assert config.min_depth == 20
        assert config.min_mapq == 10

    @pytest.mark.unit
    def test_from_env_defaults(self, monkeypatch):
        """from_env with no env vars should return defaults."""
        # Clear any existing BAMCP_ env vars
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        config = BAMCPConfig.from_env()
        assert config.reference is None
        assert config.max_reads == DEFAULT_MAX_READS
        assert config.default_window == DEFAULT_WINDOW_SIZE
        assert config.min_vaf == DEFAULT_MIN_VAF
        assert config.min_depth == DEFAULT_MIN_DEPTH
        assert config.min_mapq == DEFAULT_MIN_MAPQ

    @pytest.mark.unit
    def test_from_env_custom(self, monkeypatch):
        """from_env should read environment variables."""
        monkeypatch.setenv("BAMCP_REFERENCE", "/data/hg38.fa")
        monkeypatch.setenv("BAMCP_MAX_READS", "5000")
        monkeypatch.setenv("BAMCP_DEFAULT_WINDOW", "1000")
        monkeypatch.setenv("BAMCP_MIN_VAF", "0.05")
        monkeypatch.setenv("BAMCP_MIN_DEPTH", "20")
        monkeypatch.setenv("BAMCP_MIN_MAPQ", "30")

        config = BAMCPConfig.from_env()
        assert config.reference == "/data/hg38.fa"
        assert config.max_reads == 5000
        assert config.default_window == 1000
        assert config.min_vaf == 0.05
        assert config.min_depth == 20
        assert config.min_mapq == 30

    @pytest.mark.unit
    def test_from_env_partial(self, monkeypatch):
        """from_env should handle partial env var overrides."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        monkeypatch.setenv("BAMCP_MAX_READS", "500")

        config = BAMCPConfig.from_env()
        assert config.max_reads == 500
        assert config.reference is None  # default
        assert config.min_vaf == DEFAULT_MIN_VAF  # default

    @pytest.mark.unit
    def test_from_env_transport_fields(self, monkeypatch):
        """from_env should read transport-related env vars."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        monkeypatch.setenv("BAMCP_TRANSPORT", "sse")
        monkeypatch.setenv("BAMCP_HOST", "127.0.0.1")
        monkeypatch.setenv("BAMCP_PORT", "9000")

        config = BAMCPConfig.from_env()
        assert config.transport == "sse"
        assert config.host == "127.0.0.1"
        assert config.port == 9000

    @pytest.mark.unit
    def test_from_env_auth_fields(self, monkeypatch):
        """from_env should read auth-related env vars."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        monkeypatch.setenv("BAMCP_AUTH_ENABLED", "true")
        monkeypatch.setenv("BAMCP_ISSUER_URL", "https://auth.example.com")
        monkeypatch.setenv("BAMCP_RESOURCE_SERVER_URL", "https://mcp.example.com")
        monkeypatch.setenv("BAMCP_REQUIRED_SCOPES", "read,write")
        monkeypatch.setenv("BAMCP_TOKEN_EXPIRY", "7200")

        config = BAMCPConfig.from_env()
        assert config.auth_enabled is True
        assert config.issuer_url == "https://auth.example.com"
        assert config.resource_server_url == "https://mcp.example.com"
        assert config.required_scopes == ["read", "write"]
        assert config.token_expiry == 7200

    @pytest.mark.unit
    def test_from_env_auth_disabled_by_default(self, monkeypatch):
        """Auth should be disabled when env var is empty."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        config = BAMCPConfig.from_env()
        assert config.auth_enabled is False
        assert config.required_scopes is None

    @pytest.mark.unit
    def test_from_env_cache_defaults(self, monkeypatch, tmp_path):
        """Cache should use default directory and TTL."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        config = BAMCPConfig.from_env()
        assert config.cache_dir.endswith(".cache/bamcp")
        assert config.cache_ttl == DEFAULT_CACHE_TTL_SECONDS

    @pytest.mark.unit
    def test_from_env_cache_custom(self, monkeypatch, tmp_path):
        """Cache should respect custom env vars."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        cache_dir = str(tmp_path / "custom_cache")
        monkeypatch.setenv("BAMCP_CACHE_DIR", cache_dir)
        monkeypatch.setenv("BAMCP_CACHE_TTL", "3600")

        config = BAMCPConfig.from_env()
        assert config.cache_dir == cache_dir
        assert config.cache_ttl == 3600

    @pytest.mark.unit
    def test_from_env_creates_cache_directory(self, monkeypatch, tmp_path):
        """from_env should create the cache directory if it doesn't exist."""
        for key in list(os.environ.keys()):
            if key.startswith("BAMCP_"):
                monkeypatch.delenv(key, raising=False)

        cache_dir = tmp_path / "new_cache"
        assert not cache_dir.exists()

        monkeypatch.setenv("BAMCP_CACHE_DIR", str(cache_dir))
        BAMCPConfig.from_env()

        assert cache_dir.exists()

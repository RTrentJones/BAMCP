"""Unit tests for bamcp.config module."""

import os
import pytest

from bamcp.config import BAMCPConfig


class TestBAMCPConfig:
    """Tests for BAMCPConfig dataclass."""

    @pytest.mark.unit
    def test_default_values(self):
        """Config defaults should match documented values."""
        config = BAMCPConfig()
        assert config.reference is None
        assert config.max_reads == 10000
        assert config.default_window == 500
        assert config.min_vaf == 0.1
        assert config.min_depth == 10
        assert config.min_mapq == 0

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
        assert config.max_reads == 10000
        assert config.default_window == 500
        assert config.min_vaf == 0.1
        assert config.min_depth == 10
        assert config.min_mapq == 0

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
        assert config.min_vaf == 0.1  # default

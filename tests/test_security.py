"""Unit tests for security validation."""

import pytest
from bamcp.config import BAMCPConfig
from bamcp.tools import validate_path


class TestValidatePath:
    """Tests for validate_path function."""

    @pytest.mark.unit
    def test_remote_files_disabled_by_default(self):
        config = BAMCPConfig(allow_remote_files=False)
        with pytest.raises(ValueError, match="Remote files are disabled"):
            validate_path("https://example.com/file.bam", config)

    @pytest.mark.unit
    def test_remote_files_allowed(self):
        config = BAMCPConfig(allow_remote_files=True)
        # Should not raise
        validate_path("https://example.com/file.bam", config)

    @pytest.mark.unit
    def test_remote_files_invalid_scheme(self):
        config = BAMCPConfig(allow_remote_files=True)
        with pytest.raises(ValueError, match="Scheme not supported for remote file"):
            validate_path("ftp://example.com/file.bam", config)

    @pytest.mark.unit
    def test_local_files_no_restriction(self, tmp_path):
        config = BAMCPConfig(allowed_directories=None)
        f = tmp_path / "test.bam"
        f.touch()
        # Should not raise
        validate_path(str(f), config)

    @pytest.mark.unit
    def test_local_files_allowed_directory(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        config = BAMCPConfig(allowed_directories=[str(allowed)])

        f = allowed / "test.bam"
        f.touch()

        # Should not raise
        validate_path(str(f), config)

    @pytest.mark.unit
    def test_local_files_denied_directory(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        forbidden = tmp_path / "forbidden"
        forbidden.mkdir()

        config = BAMCPConfig(allowed_directories=[str(allowed)])

        f = forbidden / "test.bam"
        f.touch()

        with pytest.raises(ValueError, match="Path .* is not in allowed directories"):
            validate_path(str(f), config)

    @pytest.mark.unit
    def test_path_traversal_attempt(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        secret = tmp_path / "secret.txt"
        secret.touch()

        config = BAMCPConfig(allowed_directories=[str(allowed)])

        # Try to break out of allowed directory using ..
        traversal = allowed / "../secret.txt"

        with pytest.raises(ValueError, match="Path .* is not in allowed directories"):
            validate_path(str(traversal), config)

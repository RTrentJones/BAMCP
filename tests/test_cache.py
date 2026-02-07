"""Unit tests for bamcp.cache module."""

import time
from pathlib import Path

import pytest

from bamcp.cache import BAMIndexCache


class TestBAMIndexCache:
    """Tests for BAMIndexCache."""

    @pytest.mark.unit
    def test_init_creates_session_directory(self, tmp_path: Path):
        """Cache should create a session-specific subdirectory."""
        cache = BAMIndexCache(str(tmp_path), session_id="test123")

        assert cache.cache_dir.exists()
        assert cache.cache_dir.name == "test123"
        assert cache.cache_dir.parent == tmp_path

    @pytest.mark.unit
    def test_init_generates_session_id(self, tmp_path: Path):
        """Cache should generate a session ID if not provided."""
        cache = BAMIndexCache(str(tmp_path))

        assert cache.session_id is not None
        assert len(cache.session_id) == 8  # First 8 chars of UUID

    @pytest.mark.unit
    def test_get_index_path_local_file_returns_none(self, tmp_path: Path):
        """Local files should return None (let pysam find the index)."""
        cache = BAMIndexCache(str(tmp_path))

        assert cache.get_index_path("/path/to/local.bam") is None
        assert cache.get_index_path("relative/path.bam") is None
        assert cache.get_index_path("./local.bam") is None

    @pytest.mark.unit
    def test_get_index_path_remote_bam(self, tmp_path: Path):
        """Remote BAM URLs should return a cache path in session dir."""
        cache = BAMIndexCache(str(tmp_path), session_id="sess1")

        url = "https://example.com/data/sample.bam"
        index_path = cache.get_index_path(url)

        assert index_path is not None
        assert "sess1" in index_path
        assert index_path.endswith(".bai")
        assert "sample.bam" in index_path

    @pytest.mark.unit
    def test_get_index_path_remote_cram(self, tmp_path: Path):
        """Remote CRAM URLs should return a .crai cache path."""
        cache = BAMIndexCache(str(tmp_path))

        url = "https://example.com/data/sample.cram"
        index_path = cache.get_index_path(url)

        assert index_path is not None
        assert index_path.endswith(".crai")
        assert "sample.cram" in index_path

    @pytest.mark.unit
    def test_get_index_path_http_url(self, tmp_path: Path):
        """HTTP URLs (not just HTTPS) should work."""
        cache = BAMIndexCache(str(tmp_path))

        url = "http://example.com/data/sample.bam"
        index_path = cache.get_index_path(url)

        assert index_path is not None
        assert index_path.endswith(".bai")

    @pytest.mark.unit
    def test_get_index_path_unique_per_url(self, tmp_path: Path):
        """Different URLs should produce different cache paths."""
        cache = BAMIndexCache(str(tmp_path))

        url1 = "https://example.com/data/sample.bam"
        url2 = "https://other.com/data/sample.bam"

        path1 = cache.get_index_path(url1)
        path2 = cache.get_index_path(url2)

        assert path1 != path2

    @pytest.mark.unit
    def test_get_index_path_deterministic(self, tmp_path: Path):
        """Same URL should always produce the same cache path."""
        cache = BAMIndexCache(str(tmp_path))

        url = "https://example.com/data/sample.bam"

        path1 = cache.get_index_path(url)
        path2 = cache.get_index_path(url)

        assert path1 == path2

    @pytest.mark.unit
    def test_sessions_isolated(self, tmp_path: Path):
        """Different sessions should have different cache directories."""
        cache1 = BAMIndexCache(str(tmp_path), session_id="session1")
        cache2 = BAMIndexCache(str(tmp_path), session_id="session2")

        url = "https://example.com/data/sample.bam"
        path1 = cache1.get_index_path(url)
        path2 = cache2.get_index_path(url)

        assert path1 != path2
        assert "session1" in path1
        assert "session2" in path2

    @pytest.mark.unit
    def test_is_valid_nonexistent_file(self, tmp_path: Path):
        """Nonexistent files should be invalid."""
        cache = BAMIndexCache(str(tmp_path))

        assert not cache.is_valid(str(tmp_path / "nonexistent.bai"))

    @pytest.mark.unit
    def test_is_valid_fresh_file(self, tmp_path: Path):
        """Recently created files should be valid."""
        cache = BAMIndexCache(str(tmp_path), ttl_seconds=3600)

        cache_file = cache.cache_dir / "test.bai"
        cache_file.write_bytes(b"test content")

        assert cache.is_valid(str(cache_file))

    @pytest.mark.unit
    def test_is_valid_expired_file(self, tmp_path: Path):
        """Files older than TTL should be invalid."""
        cache = BAMIndexCache(str(tmp_path), ttl_seconds=1)

        cache_file = cache.cache_dir / "test.bai"
        cache_file.write_bytes(b"test content")

        # Wait for expiration
        time.sleep(1.1)

        assert not cache.is_valid(str(cache_file))

    @pytest.mark.unit
    def test_cleanup_session_removes_files(self, tmp_path: Path):
        """cleanup_session should remove all files in session directory."""
        cache = BAMIndexCache(str(tmp_path), session_id="mysession")

        # Create files in session directory
        (cache.cache_dir / "file1.bai").write_bytes(b"content1")
        (cache.cache_dir / "file2.bai").write_bytes(b"content2")

        removed = cache.cleanup_session()

        assert removed == 2
        assert not cache.cache_dir.exists()  # Directory should be removed too

    @pytest.mark.unit
    def test_cleanup_session_empty(self, tmp_path: Path):
        """cleanup_session on empty session should return zero."""
        cache = BAMIndexCache(str(tmp_path))

        removed = cache.cleanup_session()

        assert removed == 0

    @pytest.mark.unit
    def test_cleanup_session_only_affects_own_session(self, tmp_path: Path):
        """cleanup_session should not affect other sessions."""
        cache1 = BAMIndexCache(str(tmp_path), session_id="session1")
        cache2 = BAMIndexCache(str(tmp_path), session_id="session2")

        # Create files in both sessions
        (cache1.cache_dir / "file1.bai").write_bytes(b"content1")
        (cache2.cache_dir / "file2.bai").write_bytes(b"content2")

        # Cleanup only session1
        removed = cache1.cleanup_session()

        assert removed == 1
        assert not cache1.cache_dir.exists()
        assert cache2.cache_dir.exists()
        assert (cache2.cache_dir / "file2.bai").exists()

    @pytest.mark.unit
    def test_cleanup_expired_sessions_removes_old(self, tmp_path: Path):
        """cleanup_expired_sessions should remove expired session directories."""
        # Create an old session with expired files
        old_session = tmp_path / "old_session"
        old_session.mkdir()
        old_file = old_session / "old.bai"
        old_file.write_bytes(b"old content")
        time.sleep(1.1)

        # Create current session
        cache = BAMIndexCache(str(tmp_path), ttl_seconds=1, session_id="current")

        removed, remaining = cache.cleanup_expired_sessions()

        assert removed == 1
        assert remaining == 0
        assert not old_session.exists()
        assert cache.cache_dir.exists()  # Current session untouched

    @pytest.mark.unit
    def test_cleanup_expired_sessions_keeps_active(self, tmp_path: Path):
        """cleanup_expired_sessions should keep sessions with valid files."""
        # Create another session with fresh files
        other_session = tmp_path / "other"
        other_session.mkdir()
        (other_session / "fresh.bai").write_bytes(b"fresh")

        # Create current session
        cache = BAMIndexCache(str(tmp_path), ttl_seconds=3600, session_id="current")

        removed, remaining = cache.cleanup_expired_sessions()

        assert removed == 0
        assert remaining == 1
        assert other_session.exists()

    @pytest.mark.unit
    def test_cleanup_expired_sessions_ignores_current(self, tmp_path: Path):
        """cleanup_expired_sessions should never remove current session."""
        cache = BAMIndexCache(str(tmp_path), ttl_seconds=0, session_id="current")

        # Create a file that would be "expired"
        (cache.cache_dir / "test.bai").write_bytes(b"content")
        time.sleep(0.1)

        removed, remaining = cache.cleanup_expired_sessions()

        # Current session should not be removed even with TTL=0
        assert cache.cache_dir.exists()

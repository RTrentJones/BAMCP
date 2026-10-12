"""Unit tests for bamcp.cache module."""

import os
import time
from pathlib import Path

import pytest

from bamcp.constants import CACHE_SESSION_ID_LENGTH, DEFAULT_CACHE_TTL_SECONDS
from bamcp.core.cache import BAMIndexCache


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
        assert len(cache.session_id) == CACHE_SESSION_ID_LENGTH

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
    def test_get_index_path_signed_uppercase_cram(self, tmp_path: Path):
        """A signed (query-string) or upper-case CRAM URL still maps to a .crai."""
        cache = BAMIndexCache(str(tmp_path))

        index_path = cache.get_index_path("https://example.com/d/S1.CRAM?X-Amz-Signature=a/b.bam")

        assert index_path is not None
        assert index_path.endswith("_S1.CRAM.crai")
        # The query string never leaks into the filename.
        assert "Signature" not in index_path
        assert Path(index_path).parent == cache.cache_dir

    @pytest.mark.unit
    def test_get_index_path_long_name_truncated(self, tmp_path: Path):
        """An over-long URL basename is truncated so the filename stays valid."""
        cache = BAMIndexCache(str(tmp_path))

        index_path = cache.get_index_path("https://example.com/" + "a" * 400 + ".bam")

        assert index_path is not None
        assert len(Path(index_path).name) < 255

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

    @pytest.mark.unit
    def test_default_ttl_uses_shared_constant(self, tmp_path: Path):
        """Default cache TTL should come from shared constants."""
        cache = BAMIndexCache(str(tmp_path))
        assert cache.ttl_seconds == DEFAULT_CACHE_TTL_SECONDS


def _write(path: Path, size: int, age_seconds: float) -> Path:
    path.write_bytes(b"x" * size)
    ts = time.time() - age_seconds
    os.utime(path, (ts, ts))
    return path


class TestIndexCacheBudget:
    """The session's index cache is bounded by a total byte budget (LRU by mtime)."""

    @pytest.mark.unit
    def test_evicts_oldest_until_within_budget(self, tmp_path: Path):
        cache = BAMIndexCache(str(tmp_path), session_id="s", max_total_bytes=250)
        oldest = _write(cache.cache_dir / "a.bai", 100, 30)
        middle = _write(cache.cache_dir / "b.bai", 100, 20)
        newest = _write(cache.cache_dir / "c.bai", 100, 10)

        assert cache.enforce_budget() == 1
        assert not oldest.exists()
        assert middle.exists() and newest.exists()

    @pytest.mark.unit
    def test_keep_file_never_evicted(self, tmp_path: Path):
        cache = BAMIndexCache(str(tmp_path), session_id="s", max_total_bytes=50)
        old = _write(cache.cache_dir / "a.bai", 40, 30)
        kept = _write(cache.cache_dir / "big.bai", 100, 60)  # oldest AND over budget alone

        cache.enforce_budget(keep=str(kept))

        assert kept.exists()
        assert not old.exists()

    @pytest.mark.unit
    def test_expired_files_dropped_and_tmp_ignored(self, tmp_path: Path):
        cache = BAMIndexCache(str(tmp_path), session_id="s", ttl_seconds=100)
        expired = _write(cache.cache_dir / "old.bai", 10, 200)
        fresh = _write(cache.cache_dir / "new.bai", 10, 1)
        in_flight = _write(cache.cache_dir / ".index-abc.tmp", 10, 200)

        assert cache.enforce_budget() == 1
        assert not expired.exists()
        assert fresh.exists() and in_flight.exists()

    @pytest.mark.unit
    def test_missing_cache_dir_is_noop(self, tmp_path: Path):
        cache = BAMIndexCache(str(tmp_path), session_id="s")
        cache.cache_dir.rmdir()
        assert cache.enforce_budget() == 0

    @pytest.mark.unit
    def test_unlink_failure_tolerated(self, tmp_path: Path):
        assert BAMIndexCache._unlink(tmp_path / "does-not-exist") is False

    @pytest.mark.unit
    def test_cleanup_expired_sessions_missing_base_dir(self, tmp_path: Path):
        cache = BAMIndexCache(str(tmp_path / "base"), session_id="s")
        cache.cache_dir.rmdir()
        cache.base_cache_dir.rmdir()
        assert cache.cleanup_expired_sessions() == (0, 0)


class TestServicesCacheWiring:
    """BAMCPServices wires the budget from config and sweeps expired sessions at init."""

    @pytest.mark.unit
    def test_budget_from_config_and_expired_sweep(self, tmp_path: Path):
        from bamcp.config import BAMCPConfig
        from bamcp.core.tools import BAMCPServices

        stale = tmp_path / "stale_session"
        stale.mkdir()
        _write(stale / "old.bai", 10, 10_000)

        config = BAMCPConfig(cache_dir=str(tmp_path), cache_ttl=60, max_index_cache_bytes=1234)
        services = BAMCPServices.from_config(config)

        assert services.cache.max_total_bytes == 1234
        assert not stale.exists()

    @pytest.mark.unit
    def test_config_budget_from_env_and_validation(self, monkeypatch):
        from bamcp.config import BAMCPConfig
        from bamcp.constants import DEFAULT_MAX_INDEX_CACHE_BYTES

        assert BAMCPConfig().max_index_cache_bytes == DEFAULT_MAX_INDEX_CACHE_BYTES
        monkeypatch.setenv("BAMCP_MAX_INDEX_CACHE_BYTES", "4096")
        assert BAMCPConfig.from_env().max_index_cache_bytes == 4096
        with pytest.raises(ValueError, match="max_index_cache_bytes must be at least 1"):
            BAMCPConfig(max_index_cache_bytes=0)

"""Cache management for remote BAM index files."""

from __future__ import annotations

import contextlib
import hashlib
import uuid
from datetime import datetime
from pathlib import Path
from urllib.parse import urlparse

from ..constants import (
    CACHE_SESSION_ID_LENGTH,
    DEFAULT_CACHE_TTL_SECONDS,
    DEFAULT_MAX_INDEX_CACHE_BYTES,
    REMOTE_FILE_SCHEMES,
)
from .validation import is_cram

# Cap on the URL-derived part of a cache filename (filesystems limit names to 255 bytes).
_MAX_BASENAME_CHARS = 100


class BAMIndexCache:
    """File-based cache for remote BAM index files.

    When pysam opens a remote BAM via URL, it downloads the index file.
    This cache redirects those downloads to a session-specific subdirectory
    to avoid conflicts between concurrent sessions.

    Each cache instance gets its own session ID, isolating its files from
    other sessions.
    """

    def __init__(
        self,
        cache_dir: str,
        ttl_seconds: int = DEFAULT_CACHE_TTL_SECONDS,
        session_id: str | None = None,
        max_total_bytes: int = DEFAULT_MAX_INDEX_CACHE_BYTES,
    ):
        """Initialize the cache.

        Args:
            cache_dir: Base directory for cached index files.
            ttl_seconds: Time-to-live in seconds (default: 24 hours).
            session_id: Unique session identifier. If None, generates a UUID.
            max_total_bytes: On-disk budget for this session's cached indexes;
                :meth:`enforce_budget` evicts the oldest files beyond it.
        """
        self.base_cache_dir = Path(cache_dir)
        self.ttl_seconds = ttl_seconds
        self.max_total_bytes = max_total_bytes
        self.session_id = session_id or str(uuid.uuid4())[:CACHE_SESSION_ID_LENGTH]
        self.cache_dir = self.base_cache_dir / self.session_id
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_index_path(self, bam_path: str) -> str | None:
        """Get cache path for a remote BAM's index file.

        Returns None for local files (let pysam find the index normally).

        Args:
            bam_path: Path or URL to the BAM/CRAM file.

        Returns:
            Cache path for the index file, or None for local files.
        """
        if not bam_path.startswith(REMOTE_FILE_SCHEMES):
            return None

        # Hash URL + basename for unique filename (not used for security). The
        # basename comes from the URL *path* so a signed URL's query string can't
        # leak '/'-split fragments or an over-long name into the cache filename.
        hash_prefix = hashlib.md5(bam_path.encode()).hexdigest()[:8]  # noqa: S324
        basename = (Path(urlparse(bam_path).path).name or "remote")[:_MAX_BASENAME_CHARS]
        suffix = ".crai" if is_cram(bam_path) else ".bai"
        return str(self.cache_dir / f"{hash_prefix}_{basename}{suffix}")

    def get_reference_index_paths(self, reference_url: str) -> tuple[str, str | None] | None:
        """Cache paths for a remote FASTA's ``.fai`` (and ``.gzi`` when bgzip-compressed).

        Returns None for local references (htslib reads their index in place).
        """
        if not reference_url.startswith(REMOTE_FILE_SCHEMES):
            return None
        hash_prefix = hashlib.md5(reference_url.encode()).hexdigest()[:8]  # noqa: S324
        url_path = urlparse(reference_url).path
        basename = (Path(url_path).name or "reference")[:_MAX_BASENAME_CHARS]
        stem = str(self.cache_dir / f"{hash_prefix}_{basename}")
        gzi = f"{stem}.gzi" if url_path.lower().endswith((".gz", ".bgz")) else None
        return f"{stem}.fai", gzi

    def is_valid(self, cache_path: str) -> bool:
        """Check if cached index exists and hasn't expired.

        Args:
            cache_path: Path to the cached index file.

        Returns:
            True if cache exists and is within TTL.
        """
        path = Path(cache_path)
        if not path.exists():
            return False
        age = datetime.now().timestamp() - path.stat().st_mtime
        return age < self.ttl_seconds

    def enforce_budget(self, keep: str | None = None) -> int:
        """Bound this session's cache: drop expired files, then evict oldest-first.

        Removes expired index files, then deletes the oldest (by mtime) remaining
        files until the session's total size is within ``max_total_bytes``. The
        file named by ``keep`` (the one just downloaded) is never evicted, even if
        it alone exceeds the budget. In-flight ``.index-*.tmp`` downloads are left
        alone. Filesystem races (a file vanishing concurrently) are tolerated.

        Returns:
            Number of files removed.
        """
        keep_path = Path(keep).resolve() if keep else None
        now = datetime.now().timestamp()
        entries: list[tuple[float, int, Path]] = []
        removed = 0
        try:
            candidates = [p for p in self.cache_dir.iterdir() if not p.name.endswith(".tmp")]
        except OSError:
            return 0
        for path in candidates:
            try:
                st = path.stat()
            except OSError:
                continue
            if not path.is_file():
                continue
            is_keep = keep_path is not None and path.resolve() == keep_path
            if not is_keep and now - st.st_mtime >= self.ttl_seconds:
                if self._unlink(path):
                    removed += 1
                continue
            entries.append((st.st_mtime, st.st_size, path))

        total = sum(size for _, size, _ in entries)
        for _mtime, size, path in sorted(entries, key=lambda e: e[0]):
            if total <= self.max_total_bytes:
                break
            if keep_path is not None and path.resolve() == keep_path:
                continue
            if self._unlink(path):
                removed += 1
            total -= size
        return removed

    @staticmethod
    def _unlink(path: Path) -> bool:
        try:
            path.unlink()
        except OSError:  # incl. FileNotFoundError: already gone is fine
            return False
        return True

    def cleanup_session(self) -> int:
        """Remove the entire session cache directory.

        Removes all files downloaded during this session. Safe to call
        without affecting other sessions since each session has its own
        subdirectory.

        Returns:
            Number of files removed.
        """
        removed = 0
        if self.cache_dir.exists():
            for cache_file in self.cache_dir.glob("*"):
                if cache_file.is_file():
                    cache_file.unlink()
                    removed += 1
            # Remove the empty session directory
            with contextlib.suppress(OSError):
                self.cache_dir.rmdir()
        return removed

    def cleanup_expired_sessions(self) -> tuple[int, int]:
        """Remove expired session directories (older than TTL).

        Cleans up session directories from other sessions that have expired.
        Does not remove the current session's directory.

        Returns:
            Tuple of (sessions_removed, sessions_remaining).
        """
        removed = 0
        remaining = 0
        try:
            session_dirs = list(self.base_cache_dir.iterdir())
        except OSError:
            return 0, 0
        for session_dir in session_dirs:
            if session_dir.is_dir() and session_dir.name != self.session_id:
                # Check if any file in the session is still valid
                files = list(session_dir.glob("*"))
                if not files:
                    # Empty directory, remove it
                    try:
                        session_dir.rmdir()
                        removed += 1
                    except OSError:
                        remaining += 1
                elif all(not self.is_valid(str(f)) for f in files if f.is_file()):
                    # All files expired, remove them and the directory
                    for f in files:
                        if f.is_file():
                            f.unlink(missing_ok=True)
                    try:
                        session_dir.rmdir()
                        removed += 1
                    except OSError:
                        remaining += 1
                else:
                    remaining += 1
        return removed, remaining

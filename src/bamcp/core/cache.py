"""Cache management for remote BAM index files."""

from __future__ import annotations

import contextlib
import hashlib
import uuid
from datetime import datetime
from pathlib import Path

from ..constants import (
    CACHE_SESSION_ID_LENGTH,
    DEFAULT_CACHE_TTL_SECONDS,
    REMOTE_FILE_SCHEMES,
)


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
    ):
        """Initialize the cache.

        Args:
            cache_dir: Base directory for cached index files.
            ttl_seconds: Time-to-live in seconds (default: 24 hours).
            session_id: Unique session identifier. If None, generates a UUID.
        """
        self.base_cache_dir = Path(cache_dir)
        self.ttl_seconds = ttl_seconds
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

        # Hash URL + basename for unique filename (not used for security)
        hash_prefix = hashlib.md5(bam_path.encode()).hexdigest()[:8]  # noqa: S324
        basename = Path(bam_path).name
        suffix = ".crai" if bam_path.endswith(".cram") else ".bai"
        return str(self.cache_dir / f"{hash_prefix}_{basename}{suffix}")

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
        for session_dir in self.base_cache_dir.iterdir():
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
                            f.unlink()
                    try:
                        session_dir.rmdir()
                        removed += 1
                    except OSError:
                        remaining += 1
                else:
                    remaining += 1
        return removed, remaining

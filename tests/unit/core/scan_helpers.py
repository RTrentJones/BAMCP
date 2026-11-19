"""Run the production contig scanner (``tools._scan_contig_progressive``) synchronously.

The scan used to have a second, test-only implementation in ``parsers``; tests now drive the
code path ``scan_variants`` actually runs, with no deadline.
"""

from __future__ import annotations

import time

from bamcp.constants import SCAN_VARIANTS_CHUNK_SIZE, SCAN_VARIANTS_MAX_REGION
from bamcp.core import tools


def scan(
    bam_path: str,
    contig: str,
    reference_path: str,
    *,
    chunk_size: int = SCAN_VARIANTS_CHUNK_SIZE,
    min_vaf: float = 0.1,
    min_depth: int = 3,
    min_mapq: int = 0,
    min_baseq: int = 0,
    max_region: int = SCAN_VARIANTS_MAX_REGION,
    min_alt_count: int | None = None,
    index_filename: str | None = None,
) -> tools._ScanProgress:
    """Scan ``contig`` to completion; returns the filled-in progress record."""
    progress = tools._ScanProgress()
    tools._scan_contig_progressive(
        bam_path,
        contig,
        reference_path,
        progress,
        time.monotonic() + 3600,
        chunk_size=chunk_size,
        min_vaf=min_vaf,
        min_depth=min_depth,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        max_region=max_region,
        index_filename=index_filename,
        min_alt_count=min_alt_count,
    )
    assert not progress.deadline_hit
    return progress

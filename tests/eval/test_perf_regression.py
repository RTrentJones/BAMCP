"""Timing regression gate on a fixed fixture.

Guards against gross performance regressions in the core parse+detect path — an
accidental O(n^2) loop, a full-file rescan, a dropped short-circuit — rather than
micro-optimizations. The bound is deliberately generous (roughly 25-40x the
observed local runtime) so it stays robust on slow, shared CI runners while still
tripping on a catastrophic slowdown. We time the *uncached* ``fetch_region`` (the
in-process region cache would otherwise make repeated calls trivial) and take the
minimum of several runs, the statistic least perturbed by runner noise.
"""

from __future__ import annotations

import time

import pytest

from bamcp.core.parsers import fetch_region

# Full chr1 of the comprehensive fixture: ~230 reads, SNV + indel detection.
_REGION = "chr1:1-10000"
_RUNS = 5
# Generous ceiling: the operation runs in ~20ms locally; a real regression makes
# it many times slower. If the *fastest* of several runs exceeds this, something
# is badly wrong (not just a busy runner).
_MAX_SECONDS = 0.75


@pytest.mark.integration
def test_fetch_region_within_time_budget(comprehensive_bam_path, comprehensive_ref_fasta_path):
    bam = str(comprehensive_bam_path)
    ref = str(comprehensive_ref_fasta_path)

    fetch_region(bam, _REGION, ref)  # warm the OS/file caches

    best = min(_timed(lambda: fetch_region(bam, _REGION, ref)) for _ in range(_RUNS))
    assert best < _MAX_SECONDS, (
        f"fetch_region({_REGION}) best-of-{_RUNS} took {best * 1000:.0f}ms "
        f"(> {_MAX_SECONDS * 1000:.0f}ms) — likely a performance regression"
    )


@pytest.mark.integration
def test_fetch_region_scales_sublinearly_with_empty_tail(
    comprehensive_bam_path, comprehensive_ref_fasta_path
):
    """A region with no reads must not cost more than the read-dense one.

    Catches a regression that makes cost scale with region *width* rather than
    with the reads/coverage actually present (e.g. a per-base Python loop over
    the whole span).
    """
    bam = str(comprehensive_bam_path)
    ref = str(comprehensive_ref_fasta_path)

    # chr1:5000-9000 is reference-only (no planted reads) but wide.
    empty = min(_timed(lambda: fetch_region(bam, "chr1:5000-9000", ref)) for _ in range(_RUNS))
    dense = min(_timed(lambda: fetch_region(bam, "chr1:1000-3100", ref)) for _ in range(_RUNS))

    # The empty-but-wide region should be no more than a small multiple of the
    # read-dense one; a large blow-up means cost tracks width, not content.
    assert empty < max(dense * 5, 0.5), (
        f"empty wide region ({empty * 1000:.0f}ms) far exceeds dense region "
        f"({dense * 1000:.0f}ms) — cost may scale with region width"
    )


def _timed(fn) -> float:
    start = time.perf_counter()
    fn()
    return time.perf_counter() - start

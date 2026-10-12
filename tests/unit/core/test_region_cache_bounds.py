"""Region-cache memory bound, local-file staleness, and parse-gate loop binding."""

from __future__ import annotations

import asyncio
import os
import shutil
import threading
from collections import OrderedDict

import pytest

import bamcp.core.tools as tm
from bamcp.config import BAMCPConfig
from bamcp.core.parsers import AlignedRead, RegionData


def _region(n_reads: int, read_len: int = 100, cov: int = 0) -> RegionData:
    reads = [
        AlignedRead(
            name=f"r{i}",
            sequence="A" * read_len,
            qualities=[30] * read_len,
            cigar=f"{read_len}M",
            position=i,
            end_position=i + read_len,
            mapping_quality=60,
            is_reverse=False,
        )
        for i in range(n_reads)
    ]
    return RegionData("chr1", 0, cov, reads, [0] * cov, [])


@pytest.mark.unit
def test_approx_bytes_scales_with_payload():
    small = tm._approx_region_bytes(_region(10))
    big = tm._approx_region_bytes(_region(10_000, cov=1_000_000))
    assert big > 1000 * small
    # A full-mode 1 Mb / 10k-read entry is tens of MB — the reason a count cap is not enough.
    assert big > 20 * 1024 * 1024


@pytest.mark.unit
def test_oversized_entry_is_served_uncached(monkeypatch):
    monkeypatch.setattr(tm, "REGION_CACHE_MAX_BYTES", 100_000)
    cache: OrderedDict = OrderedDict()
    data = _region(1000)
    assert tm._set_cached_region(cache, ("big",), data, ttl=3600) is data
    assert ("big",) not in cache


@pytest.mark.unit
def test_byte_budget_evicts_lru(monkeypatch):
    one = tm._approx_region_bytes(_region(20))
    monkeypatch.setattr(tm, "REGION_CACHE_MAX_BYTES", one * 4 + 10)
    cache: OrderedDict = OrderedDict()
    for i in range(6):
        tm._set_cached_region(cache, (f"k{i}",), _region(20), ttl=3600)
    assert len(cache) == 4  # count cap (128) not hit; the byte budget evicted the oldest
    assert ("k0",) not in cache and ("k5",) in cache
    assert sum(entry[2] for entry in cache.values()) <= one * 4 + 10


@pytest.mark.unit
def test_cache_key_changes_when_local_file_is_replaced(tmp_path, small_bam_path):
    bam = tmp_path / "x.bam"
    shutil.copy(small_bam_path, bam)
    cfg = BAMCPConfig()
    key = tm._cache_key(str(bam), "chr1:1-2", None, cfg, "full", None, None, None)
    st = os.stat(bam)
    os.utime(bam, ns=(st.st_atime_ns, st.st_mtime_ns + 5_000_000_000))
    assert tm._cache_key(str(bam), "chr1:1-2", None, cfg, "full", None, None, None) != key


@pytest.mark.unit
def test_fingerprint_ignores_remote_and_missing():
    assert tm._local_fingerprint("https://example.org/x.bam") is None
    assert tm._local_fingerprint("/does/not/exist.bam") is None
    assert tm._local_fingerprint(None) is None


@pytest.mark.unit
def test_parse_gate_works_across_event_loops():
    """A module-level Semaphore binds to the first contending loop; the per-loop gate must not."""
    release = threading.Event()

    def _blocking(i):
        release.wait(timeout=5)
        return i

    async def contend():
        gate = tm._get_parse_gate()
        # Fill every permit so the next acquire must wait — which is what binds a loop.
        tasks = [
            asyncio.create_task(tm._run_in_parse_pool(_blocking, i))
            for i in range(tm.PARSE_MAX_CONCURRENCY + 1)
        ]
        await asyncio.sleep(0.05)
        assert gate.locked()
        release.set()
        out = await asyncio.gather(*tasks)
        release.clear()
        return out

    assert len(asyncio.run(contend())) == tm.PARSE_MAX_CONCURRENCY + 1
    assert len(asyncio.run(contend())) == tm.PARSE_MAX_CONCURRENCY + 1  # second loop: no error

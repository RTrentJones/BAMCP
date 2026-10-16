"""scan_variants runs against a deadline and returns partial results instead of timing out."""

import asyncio
import json
import time

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core import tools
from bamcp.core.errors import ToolInputError
from bamcp.core.parsers import scan_contig

pytestmark = pytest.mark.unit


@pytest.fixture
def ref_config(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path)


def _payload(result):
    return json.loads(result["content"][0]["text"])


async def test_full_scan_matches_scan_contig(small_bam_path, ref_fasta_path, ref_config):
    """The progressive scanner finds exactly what the whole-contig scanner finds."""
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    expected = scan_contig(
        small_bam_path,
        "chr1",
        ref_fasta_path,
        min_vaf=ref_config.min_vaf,
        min_depth=ref_config.min_depth,
        min_mapq=ref_config.min_mapq,
        min_baseq=ref_config.min_baseq,
        min_alt_count=ref_config.min_alt_count,
    )
    got = {(v["position"], v["alt"], v["vaf"]) for v in payload["variants"]}
    want = {(v["position"] + 1, v["alt"], v["vaf"]) for v in expected.variants}
    assert got == want
    assert payload["truncated"] is False
    assert "truncation_reason" not in payload
    assert payload["scanned_region"] == f"chr1:1-{expected.scanned_end}"
    assert payload["applied_filters"]["truncated_by_deadline"] is False


async def test_deadline_returns_partial_results(small_bam_path, ref_config, monkeypatch):
    """A slow chunk past the soft deadline yields truncated:true + what was scanned."""
    monkeypatch.setattr(tools, "SCAN_VARIANTS_CHUNK_SIZE", 100)
    monkeypatch.setattr(tools, "SCAN_VARIANTS_DEADLINE_SECONDS", 0.15)
    real_chunk = tools._scan_chunk
    calls = []

    def slow_chunk(*args, **kwargs):
        calls.append(args[3])
        time.sleep(0.06)
        return real_chunk(*args, **kwargs)

    monkeypatch.setattr(tools, "_scan_chunk", slow_chunk)
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    assert payload["truncated"] is True
    assert "time budget" in payload["truncation_reason"]
    assert payload["applied_filters"]["truncated_by_deadline"] is True
    scanned_end = int(payload["scanned_region"].rsplit("-", 1)[1])
    assert payload["scanned_region"].startswith("chr1:1-")
    assert scanned_end == len(calls) * 100
    assert scanned_end < payload["applied_filters"]["contig_length"]
    # Only candidates inside the scanned region are reported (1-based positions).
    assert all(v["position"] <= scanned_end for v in payload["variants"])
    # The remainder is named exactly, as a region scan_variants accepts.
    contig_length = payload["applied_filters"]["contig_length"]
    assert payload["next_region"] == f"chr1:{scanned_end + 1}-{contig_length}"
    assert f"region='{payload['next_region']}'" in payload["truncation_reason"]


async def test_region_resumes_a_truncated_scan(small_bam_path, ref_config):
    """Scanning the remainder region covers exactly the bases a full scan would."""
    full = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    first = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "region": "chr1:1-300"}, ref_config
        )
    )
    rest_end = full["applied_filters"]["contig_length"]
    rest = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "region": f"chr1:301-{rest_end}"}, ref_config
        )
    )
    assert first["scanned_region"] == "chr1:1-300"
    assert rest["scanned_region"] == f"chr1:301-{rest_end}"
    assert all(v["position"] <= 300 for v in first["variants"])
    assert all(v["position"] > 300 for v in rest["variants"])
    key = lambda v: (v["position"], v["alt"])  # noqa: E731
    assert sorted(map(key, first["variants"] + rest["variants"])) == sorted(
        map(key, full["variants"])
    )


async def test_scan_output_sorted_by_position_with_interpretation(small_bam_path, ref_config):
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1", "min_alt_count": 1, "min_vaf": 0.05},
            ref_config,
        )
    )
    positions = [v["position"] for v in payload["variants"]]
    assert positions == sorted(positions)
    assert all("interpretation" in v for v in payload["variants"])


async def test_region_and_contig_must_agree(small_bam_path, ref_config):
    with pytest.raises(ToolInputError, match="does not match"):
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr2", "region": "chr1:1-100"}, ref_config
        )


async def test_hard_timeout_returns_partial_results(small_bam_path, ref_config, monkeypatch):
    """A single chunk hung past the hard timeout still answers with the completed chunks."""
    monkeypatch.setattr(tools, "SCAN_VARIANTS_CHUNK_SIZE", 100)
    monkeypatch.setattr(tools, "SCAN_VARIANTS_DEADLINE_SECONDS", 30.0)
    monkeypatch.setattr(tools, "SCAN_VARIANTS_TIMEOUT_SECONDS", 0.2)
    real_chunk = tools._scan_chunk

    def hang_second(*args, **kwargs):
        if args[3] >= 100:
            time.sleep(0.4)
        return real_chunk(*args, **kwargs)

    monkeypatch.setattr(tools, "_scan_chunk", hang_second)
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    assert payload["truncated"] is True
    assert payload["scanned_region"] == "chr1:1-100"
    assert "0.2s" in payload["truncation_reason"]
    await asyncio.sleep(0.5)  # let the abandoned worker see the stop flag and exit


async def test_nothing_scanned_is_an_error(small_bam_path, ref_config, monkeypatch):
    """No header read at all → no partial answer to give → an error, not an empty result."""
    monkeypatch.setattr(tools, "SCAN_VARIANTS_TIMEOUT_SECONDS", 0.05)

    def stuck(*args, **kwargs):
        time.sleep(0.2)

    monkeypatch.setattr(tools, "_scan_contig_progressive", stuck)
    with pytest.raises(ToolInputError, match="before any data"):
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    await asyncio.sleep(0.25)


async def test_min_alt_count_is_honored_and_echoed(small_bam_path, ref_config):
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1", "min_alt_count": 10_000}, ref_config
        )
    )
    assert payload["variants"] == []
    assert payload["applied_filters"]["min_alt_count"] == 10_000
    assert payload["truncated"] is False


async def test_region_start_beyond_contig_is_an_error(small_bam_path, ref_config):
    """A start past the contig end must not read as a successful, variant-free scan."""
    full = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    past = full["applied_filters"]["contig_length"] + 100
    with pytest.raises(ToolInputError, match="beyond the end of chr1"):
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "region": f"chr1:{past}-{past + 50}"}, ref_config
        )


async def test_max_region_cap_is_reported_with_continuation(
    small_bam_path, ref_config, monkeypatch
):
    """A scan capped by SCAN_VARIANTS_MAX_REGION is truncated, with a next_region through the
    caller's requested stop — never truncated:false over an unscanned remainder."""
    monkeypatch.setattr(tools, "SCAN_VARIANTS_MAX_REGION", 300)
    payload = _payload(
        await tools.handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, ref_config
        )
    )
    contig_length = payload["applied_filters"]["contig_length"]
    assert payload["truncated"] is True
    assert payload["applied_filters"]["truncated_by_max_region"] is True
    assert payload["applied_filters"]["truncated_by_deadline"] is False
    assert payload["scanned_region"] == "chr1:1-300"
    assert payload["next_region"] == f"chr1:301-{contig_length}"
    assert "at most 300 bp" in payload["truncation_reason"]


async def test_alt_count_floor_applied_before_deletion_pileups(
    small_bam_path, ref_config, monkeypatch
):
    """Candidates under min_alt_count must never reach the (slow) deletion-depth pileup."""
    from bamcp.core import parsers

    seen: list[int] = []
    real = parsers._deletion_depth_at

    def spy(samfile, contig, positions, min_mapq):
        seen.extend(positions)
        return real(samfile, contig, positions, min_mapq)

    monkeypatch.setattr(parsers, "_deletion_depth_at", spy)
    await tools.handle_scan_variants(
        {"file_path": small_bam_path, "contig": "chr1", "min_alt_count": 10_000}, ref_config
    )
    assert seen == []

"""Unit tests for the GIAB generator's pure helpers (no network, no pysam).

The download / VCF-stream / read-sim paths are network-gated and covered by
running fetch_giab.py manually (see GIAB_RESULTS.md). Here we lock the logic
that is easy to get subtly wrong: BED interval lookup, negative-region
selection, coordinate parsing, and the error-rate → Phred mapping.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

# Load the module by path — datasets/giab is not an importable package.
# Resolve relative to this file so collection works from any CWD, and register
# in sys.modules before exec so @dataclass can resolve the module.
_PATH = Path(__file__).resolve().parents[2] / "eval/datasets/giab/fetch_giab.py"
_spec = importlib.util.spec_from_file_location("giab_fetch", _PATH)
giab = importlib.util.module_from_spec(_spec)
sys.modules["giab_fetch"] = giab
_spec.loader.exec_module(giab)


def test_load_bed_intervals_clips_to_region(tmp_path):
    bed = tmp_path / "hc.bed"
    bed.write_text(
        "chr20\t900\t1100\tx\n"  # overlaps left edge
        "chr20\t2000\t2500\tx\n"  # fully inside
        "chr20\t9000\t9999\tx\n"  # outside region
        "chr21\t1000\t2000\tx\n",  # wrong chrom
        encoding="utf-8",
    )
    intervals = giab.load_bed_intervals(bed, "chr20", 1000, 3000)
    assert intervals == [(1000, 1100), (2000, 2500)]


def test_in_intervals_boundaries():
    intervals = [(1000, 1100), (2000, 2500)]
    starts = [iv[0] for iv in intervals]
    assert giab._in_intervals(1000, starts, intervals)  # inclusive start
    assert not giab._in_intervals(1100, starts, intervals)  # exclusive end
    assert giab._in_intervals(2499, starts, intervals)
    assert not giab._in_intervals(1500, starts, intervals)  # in the gap
    assert not giab._in_intervals(500, starts, intervals)  # before all


def test_pick_negative_regions_avoids_truth():
    intervals = [(1000, 10000)]
    truth = [giab.TruthSNV(pos0=5000, ref="A", alt="T", zygosity="het")]
    negs = giab._pick_negative_regions(intervals, truth, want=2, span=1000)
    assert len(negs) == 2
    # None of the chosen windows may contain the truth SNV at 5000.
    for start, end in negs:
        assert not (start <= 5000 < end)


def test_phred_for_error_rate():
    assert giab._phred_for_error_rate(0.0) == 40
    assert giab._phred_for_error_rate(0.001) == 30  # Q30
    assert giab._phred_for_error_rate(0.01) == 20  # Q20
    # Clamped into [2, 40].
    assert giab._phred_for_error_rate(0.9) == 2
    assert giab._phred_for_error_rate(1e-9) == 40


def test_parse_region():
    assert giab._parse_region("chr20:1000000-1060000") == ("chr20", 1000000, 1060000)


def test_truthsnv_is_zero_based_hashable():
    t = giab.TruthSNV(pos0=1002294, ref="G", alt="A", zygosity="het")
    assert t.pos0 == 1002294
    assert {t}  # frozen/hashable

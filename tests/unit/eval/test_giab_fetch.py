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


def test_truthsnv_is_hashable():
    t = giab.TruthSNV(pos0=1002294, ref="G", alt="A", zygosity="het")
    assert {t}  # frozen/hashable


def test_truthsnv_exposes_both_coordinate_bases():
    """The two bases are distinct and must not be conflated.

    ``pos0`` is simulation space (pysam, 0-based). ``pos1`` is what goes in the
    manifest, because the scorer matches it against ``get_variants`` output,
    which ``core/tools.py::_one_based`` converts to 1-based.
    """
    t = giab.TruthSNV(pos0=1002294, ref="G", alt="A", zygosity="het")
    assert t.pos0 == 1002294
    assert t.pos1 == 1002295
    assert t.pos1 == t.pos0 + 1


def _write_manifest(tmp_path, truth, negatives=()):
    import yaml

    out = tmp_path / "manifest.yaml"
    giab.write_manifest(
        out,
        "chr20",
        1000000,
        1060000,
        list(truth),
        list(negatives),
        tmp_path / "chr20.fa",
        tmp_path / "reads.bam",
        depth=30,
        error_rate=0.001,
        seed=1,
    )
    return yaml.safe_load(out.read_text(encoding="utf-8"))


def test_manifest_positions_are_one_based(tmp_path):
    """The regression guard for the bug this test file previously locked in.

    A 0-based manifest is not a loud failure — the scorer's exact
    ``(chrom, pos, ref, alt)`` match simply finds nothing and reports recall
    0.0, which reads as "the caller regressed" rather than "the manifest is in
    the wrong coordinate space".
    """
    truth = [
        giab.TruthSNV(pos0=1002294, ref="G", alt="A", zygosity="het"),
        giab.TruthSNV(pos0=1002466, ref="G", alt="C", zygosity="hom"),
    ]
    doc = _write_manifest(tmp_path, truth)

    positions = [s["pos"] for s in doc["variant_sites"]]
    assert positions == [1002295, 1002467]
    assert all(s["pos"] == t.pos1 for s, t in zip(doc["variant_sites"], truth, strict=True))


def test_manifest_negative_region_note_is_quoted(tmp_path):
    """An unquoted comma would split the note into a second, null-valued key."""
    doc = _write_manifest(
        tmp_path,
        [giab.TruthSNV(pos0=1002294, ref="G", alt="A", zygosity="het")],
        negatives=[(1000000, 1002000)],
    )

    (negative,) = doc["negative_regions"]
    assert set(negative) == {"region", "note"}
    assert negative["note"] == "high-confidence, no truth SNV"

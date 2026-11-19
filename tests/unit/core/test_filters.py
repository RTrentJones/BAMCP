"""The one candidate-filter predicate (core.filters) and every tool-layer filter built on it."""

from __future__ import annotations

from bamcp.config import BAMCPConfig
from bamcp.core.filters import failed_candidate_filters, passes_candidate_filters
from bamcp.core.parsers import RegionData
from bamcp.core.tools import _apply_alt_count, _select_get_variants, site_filter_check


def _cand(alt_count: int | None, *, vaf: float = 0.5, depth: int = 20, source: str = "bamcp"):
    v = {"position": 10, "ref": "A", "alt": "G", "vaf": vaf, "depth": depth, "source": source}
    if alt_count is not None:
        v["alt_count"] = alt_count
    return v


def _data(variants: list[dict]) -> RegionData:
    return RegionData(
        contig="chr1", start=0, end=100, reads=[], coverage=[], variants=variants,
        reference_sequence=None,
    )  # fmt: skip


class TestPredicate:
    def test_unknown_alt_count_is_not_dropped_by_min_alt_count(self):
        assert passes_candidate_filters(_cand(None), min_alt_count=5)
        v = {**_cand(0), "alt_count": None}  # explicit None, as well as a missing key
        assert passes_candidate_filters(v, min_alt_count=5)

    def test_known_alt_count_below_floor_fails(self):
        assert not passes_candidate_filters(_cand(2), min_alt_count=3)
        assert passes_candidate_filters(_cand(3), min_alt_count=3)

    def test_failed_lists_every_floor_in_order(self):
        v = _cand(1, vaf=0.01, depth=2)
        assert failed_candidate_filters(v, min_vaf=0.1, min_depth=5, min_alt_count=2) == [
            "min_vaf",
            "min_depth",
            "min_alt_count",
        ]

    def test_none_threshold_is_not_applied(self):
        assert failed_candidate_filters(_cand(1, vaf=0.0, depth=0)) == []

    def test_vcf_records_always_pass(self):
        assert passes_candidate_filters(
            _cand(0, vaf=0.0, depth=0, source="vcf"), min_vaf=0.5, min_depth=10, min_alt_count=3
        )


class TestCallersAgree:
    """get_variants' selection, the viewer/summary alt-read floor and filter_check agree."""

    def test_alt_count_floor_keeps_unknown_drops_low(self):
        data = _data([_cand(None), _cand(1), _cand(5)])
        kept = _apply_alt_count(data, 3).variants
        assert [v.get("alt_count") for v in kept] == [None, 5]

    def test_get_variants_selection_matches(self):
        data = _data([_cand(None), _cand(1), _cand(5)])
        kept, total = _select_get_variants(data, 0.0, 1, 3)
        assert total == 2
        assert sorted(v.get("alt_count") or -1 for v in kept) == [-1, 5]

    def test_site_filter_check_unknown_alt_count(self):
        config = BAMCPConfig(min_alt_count=3, min_vaf=0.1, min_depth=5)
        check = site_filter_check(_cand(None), config)
        assert check["passes_filters"]
        assert check["at_filters"]["alt_count"] is None

    def test_site_filter_check_low_alt_count_and_baseq(self):
        config = BAMCPConfig(min_alt_count=3, min_vaf=0.1, min_depth=5)
        check = site_filter_check(
            _cand(1), config, unfiltered={"alt_count": 4, "depth": 20, "vaf": 0.2}
        )
        assert check["failed_filters"] == ["min_baseq", "min_alt_count"]
        assert check["low_baseq_alt_reads_removed"] == 3

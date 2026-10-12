"""applied_filters must make every empty/partial result self-describing.

Detection disabled (no reference / VCF-authoritative), indel detection skipped on a
large region, max_reads truncation, the fixed read exclusions, and scan_variants'
scope (scanned range, SNV-only, result cap) are all echoed — so a variant's absence
is never indistinguishable from "filtered out" or "not looked for".
"""

from __future__ import annotations

import json

import pysam
import pytest

from bamcp.config import BAMCPConfig
from bamcp.constants import INDEL_DETECTION_MAX_REGION
from bamcp.core.parsers import fetch_region, scan_contig
from bamcp.core.serialization import build_enhanced_variants, serialize_region_data
from bamcp.core.tools import (
    _format_filters_line,
    handle_get_region_summary,
    handle_get_variants,
    handle_scan_variants,
    handle_visualize_region,
)

SEQ = "ACGT" * 250  # 1000 bp


def _ref(tmp_path) -> str:
    path = str(tmp_path / "ref.fa")
    with open(path, "w") as fh:
        fh.write(">chr1\n" + SEQ + "\n")
    pysam.faidx(path)
    return path


def _bam(tmp_path, reads: list[dict]) -> str:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(SEQ)}]}
    path = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, r in enumerate(sorted(reads, key=lambda r: r["start"])):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = r["seq"]
            a.flag = r.get("flag", 0)
            a.reference_id = 0
            a.reference_start = r["start"]
            a.mapping_quality = 60
            a.cigartuples = r["cigar"]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["seq"]))
            out.write(a)
    pysam.index(path)
    return path


def _alt(pos: int) -> str:
    return "C" if SEQ[pos] != "C" else "G"


def _snv_read(start: int, length: int, pos: int | None = None, flag: int = 0) -> dict:
    seq = list(SEQ[start : start + length])
    if pos is not None:
        seq[pos - start] = _alt(pos)
    return {"start": start, "seq": "".join(seq), "cigar": [(0, length)], "flag": flag}


class TestDetectionStatus:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_get_variants_without_reference_says_detection_disabled(self, small_bam_path):
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, BAMCPConfig()
        )
        f = json.loads(result["content"][0]["text"])["applied_filters"]
        assert f["snv_detection"] == "disabled: no reference"
        assert f["indel_detection"] is False

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_get_variants_with_reference_says_enabled(self, small_bam_path, ref_fasta_path):
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:101-200"},
            BAMCPConfig(reference=ref_fasta_path),
        )
        f = json.loads(result["content"][0]["text"])["applied_filters"]
        assert f["snv_detection"] == "enabled"
        assert f["indel_detection"] is True
        assert f["excluded_reads"] == [
            "unmapped",
            "secondary",
            "supplementary",
            "duplicate",
            "qcfail",
        ]
        assert f["overlapping_mates"] == "both_counted"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_large_region_echoes_indel_detection_off(self, tmp_path):
        ref = str(tmp_path / "big.fa")
        length = INDEL_DETECTION_MAX_REGION + 5000
        with open(ref, "w") as fh:
            fh.write(">chr1\n" + "A" * length + "\n")
        pysam.faidx(ref)
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}]}
        bam = str(tmp_path / "big.bam")
        with pysam.AlignmentFile(bam, "wb", header=header):
            pass
        pysam.index(bam)
        region = f"chr1:1-{INDEL_DETECTION_MAX_REGION + 1}"
        result = await handle_get_variants(
            {"file_path": bam, "region": region}, BAMCPConfig(reference=ref)
        )
        f = json.loads(result["content"][0]["text"])["applied_filters"]
        assert f["snv_detection"] == "enabled"
        assert f["indel_detection"] is False
        assert "exceeds" in f["indel_detection_note"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_viewer_and_summary_text_say_detection_disabled(self, small_bam_path):
        cfg = BAMCPConfig()
        viz = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, cfg
        )
        assert viz["_meta"]["ui/init"]["applied_filters"]["snv_detection"] == (
            "disabled: no reference"
        )
        assert "SNV detection disabled: no reference" in viz["content"][0]["text"]
        summary = await handle_get_region_summary(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, cfg
        )
        assert "SNV detection disabled: no reference" in summary["content"][0]["text"]

    @pytest.mark.unit
    def test_filters_line_mentions_exclusions_and_truncation(self):
        line = _format_filters_line(
            {
                "min_vaf": 0.1,
                "min_depth": 5,
                "min_mapq": 0,
                "min_baseq": 0,
                "max_reads": 10,
                "reads_truncated": True,
                "excluded_reads": ["duplicate", "qcfail"],
                "snv_detection": "enabled",
                "indel_detection": False,
                "indel_detection_note": "skipped: too large",
            }
        )
        assert "TRUNCATED" in line
        assert "excluded reads: duplicate/qcfail" in line
        assert "indel detection skipped: too large" in line
        assert "SNV detection" not in line  # enabled is the unremarkable default


class TestReadsTruncated:
    def _setup(self, tmp_path):
        # Five reads at 100 fill the cap; the alt site at 320 is only seen past it.
        reads = [_snv_read(100, 50) for _ in range(5)]
        reads += [_snv_read(300, 50, pos=320, flag=0) for _ in range(3)]
        reads += [_snv_read(300, 50, pos=320, flag=16) for _ in range(2)]
        reads += [_snv_read(300, 50) for _ in range(5)]
        # A deletion past the cap too (0-based 330..332 deleted).
        del_seq = SEQ[300:330] + SEQ[333:353]
        reads += [{"start": 300, "seq": del_seq, "cigar": [(0, 30), (2, 3), (0, 20)]}] * 5
        return _bam(tmp_path, reads), _ref(tmp_path)

    @pytest.mark.unit
    def test_truncated_candidates_scored_from_pileup(self, tmp_path):
        bam, ref = self._setup(tmp_path)
        data = fetch_region(bam, "chr1:1-1000", ref, max_reads=5, min_vaf=0.1, min_depth=1)
        assert data.reads_truncated is True
        assert len(data.reads) == 5
        enhanced, evidence = build_enhanced_variants(data)
        [snv] = [v for v in enhanced if v["position"] == 320]
        assert snv["evidence_source"] == "pileup"
        # Real strand evidence, not the zeros the capped read list would give.
        assert (snv["strand_forward"], snv["strand_reverse"]) == (3, 2)
        assert evidence[f"320:{SEQ[320]}>{_alt(320)}"]["forward_count"] == 3
        assert not any(k.startswith("_alt_") for k in snv)
        [dele] = [v for v in enhanced if v.get("indel_type") == "del"]
        assert dele["evidence_incomplete"] is True

    @pytest.mark.unit
    def test_untruncated_region_is_unchanged(self, tmp_path):
        bam, ref = self._setup(tmp_path)
        data = fetch_region(bam, "chr1:1-1000", ref, max_reads=1000, min_vaf=0.1, min_depth=1)
        assert data.reads_truncated is False
        payload = serialize_region_data(data)
        [snv] = [v for v in payload["variants"] if v["position"] == 320]
        assert "evidence_source" not in snv
        assert (snv["strand_forward"], snv["strand_reverse"]) == (3, 2)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_visualize_region_echoes_reads_truncated(self, tmp_path):
        bam, ref = self._setup(tmp_path)
        cfg = BAMCPConfig(reference=ref, max_reads=5, min_vaf=0.1, min_depth=1)
        viz = await handle_visualize_region({"file_path": bam, "region": "chr1:1-1000"}, cfg)
        assert viz["_meta"]["ui/init"]["applied_filters"]["reads_truncated"] is True
        assert "TRUNCATED" in viz["content"][0]["text"]
        cfg_big = BAMCPConfig(reference=ref, max_reads=1000, min_vaf=0.1, min_depth=1)
        viz2 = await handle_visualize_region({"file_path": bam, "region": "chr1:1-1000"}, cfg_big)
        assert viz2["_meta"]["ui/init"]["applied_filters"]["reads_truncated"] is False


class TestScanScope:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_scan_echoes_scope(self, small_bam_path, ref_fasta_path):
        result = await handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"},
            BAMCPConfig(reference=ref_fasta_path),
        )
        f = json.loads(result["content"][0]["text"])["applied_filters"]
        assert f["scanned_region"] == "chr1:1-1000"
        assert f["contig_length"] == 1000
        assert f["truncated_by_max_region"] is False
        assert f["indel_detection"] is False
        assert f["results_capped"] is False
        assert f["max_results"] == 500
        assert f["total_candidates"] >= 0

    @pytest.mark.unit
    def test_scan_reports_max_region_truncation(self, small_bam_path, ref_fasta_path):
        scan = scan_contig(small_bam_path, "chr1", ref_fasta_path, max_region=500)
        assert (scan.scanned_end, scan.truncated_by_max_region) == (500, True)

    @pytest.mark.unit
    def test_scan_vaf_counts_deletion_depth_like_get_variants(self, tmp_path):
        """2 alt + 2 ref + 4 reads deleting the site: VAF is 2/8, not 2/4."""
        site = 150
        reads = [_snv_read(100, 100, pos=site) for _ in range(2)]
        reads += [_snv_read(100, 100) for _ in range(2)]
        del_seq = SEQ[100:148] + SEQ[153:200]
        reads += [{"start": 100, "seq": del_seq, "cigar": [(0, 48), (2, 5), (0, 47)]}] * 4
        bam, ref = _bam(tmp_path, reads), _ref(tmp_path)

        scan = scan_contig(bam, "chr1", ref, min_vaf=0.1, min_depth=1)
        [v] = [v for v in scan.variants if v["position"] == site]
        assert (v["depth"], v["alt_count"], v["vaf"]) == (8, 2, 0.25)

        region_data = fetch_region(bam, "chr1:101-200", ref, min_vaf=0.1, min_depth=1)
        [g] = [x for x in region_data.variants if x["position"] == site]
        assert (g["depth"], g["vaf"]) == (v["depth"], v["vaf"])

        # The deletion depth can push a candidate below min_vaf: re-filtered.
        strict = scan_contig(bam, "chr1", ref, min_vaf=0.3, min_depth=1)
        assert [x for x in strict.variants if x["position"] == site] == []

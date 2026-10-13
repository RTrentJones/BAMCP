"""Candidate-detection correctness on purpose-built BAMs.

Covers: N/IUPAC reference bases, spliced (N-CIGAR) reads, the overlapping-mates
policy, and deletions under a positive min_baseq.
"""

from __future__ import annotations

import pysam
import pytest

from bamcp.core.parsers import (
    OVERLAPPING_MATES_POLICY,
    annotate_vcf_snv_support,
    count_coverage_arrays,
    detect_indels,
    fetch_region,
)

CONTIG_LEN = 1000


def _ref(tmp_path, seq: str) -> str:
    path = str(tmp_path / "ref.fa")
    with open(path, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(path)
    return path


def _bam(tmp_path, reads: list[dict], name: str = "t.bam") -> str:
    """Write reads (dicts: start, seq, cigar, flag?, qual?, mate?) as a sorted, indexed BAM."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": CONTIG_LEN}]}
    path = str(tmp_path / name)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, r in enumerate(sorted(reads, key=lambda r: r["start"])):
            a = pysam.AlignedSegment()
            a.query_name = r.get("name", f"r{i}")
            a.query_sequence = r["seq"]
            a.flag = r.get("flag", 0)
            a.reference_id = 0
            a.reference_start = r["start"]
            a.mapping_quality = 60
            a.cigartuples = r["cigar"]
            a.query_qualities = pysam.qualitystring_to_array(r.get("qual", "I" * len(r["seq"])))
            if "mate" in r:
                a.next_reference_id = 0
                a.next_reference_start = r["mate"]
                a.template_length = r.get("tlen", 0)
            out.write(a)
    pysam.index(path)
    return path


SEQ = "ACGT" * (CONTIG_LEN // 4)


class TestAmbiguousReference:
    @pytest.mark.unit
    def test_n_reference_base_yields_no_snv(self, tmp_path):
        ref = list(SEQ)
        ref[150] = "N"
        ref_path = _ref(tmp_path, "".join(ref))
        reads = [{"start": 100, "seq": SEQ[100:200], "cigar": [(0, 100)]} for _ in range(5)]
        bam = _bam(tmp_path, reads)
        data = fetch_region(bam, "chr1:101-200", ref_path, min_vaf=0.1, min_depth=1)
        # Every read carries SEQ[150] where the reference says N: not a 100%-VAF variant.
        assert [v for v in data.variants if v["position"] == 150] == []

    @pytest.mark.unit
    def test_iupac_anchor_yields_no_indel(self, tmp_path):
        ref = list(SEQ)
        ref[149] = "R"
        ref_path = _ref(tmp_path, "".join(ref))
        ins_seq = SEQ[100:150] + "TT" + SEQ[150:198]
        reads = [{"start": 100, "seq": ins_seq, "cigar": [(0, 50), (1, 2), (0, 48)]}] * 4
        bam = _bam(tmp_path, reads)
        data = fetch_region(bam, "chr1:101-200", ref_path, min_vaf=0.1, min_depth=1)
        assert [v for v in data.variants if v["variant_kind"] == "indel"] == []

    @pytest.mark.unit
    def test_deletion_spanning_n_is_skipped(self, tmp_path):
        ref = list(SEQ)
        ref[151] = "N"
        ref_path = _ref(tmp_path, "".join(ref))
        del_seq = SEQ[100:150] + SEQ[153:200]
        reads = [{"start": 100, "seq": del_seq, "cigar": [(0, 50), (2, 3), (0, 47)]}] * 4
        bam = _bam(tmp_path, reads)
        data = fetch_region(bam, "chr1:101-200", ref_path, min_vaf=0.1, min_depth=1)
        assert [v for v in data.variants if v["variant_kind"] == "indel"] == []


class TestSplicedReads:
    @pytest.mark.unit
    def test_refskip_reads_do_not_inflate_the_snv_denominator(self, tmp_path):
        ref_path = _ref(tmp_path, SEQ)
        alt = "C" if SEQ[150] != "C" else "G"
        alt_seq = SEQ[100:150] + alt + SEQ[151:200]
        reads = [{"start": 100, "seq": SEQ[100:200], "cigar": [(0, 100)]} for _ in range(4)]
        reads.append({"start": 100, "seq": alt_seq, "cigar": [(0, 100)]})
        # Five spliced reads skip clean over position 150 (50M200N50M).
        spliced = SEQ[100:140] + SEQ[340:400]
        reads += [{"start": 100, "seq": spliced, "cigar": [(0, 40), (3, 200), (0, 60)]}] * 5
        bam = _bam(tmp_path, reads)

        _, gap_depth = detect_indels(
            bam, "chr1:101-200", ref_path, SEQ[100:200], min_vaf=0.1, min_depth=1
        )
        assert gap_depth[150 - 100] == 0  # a skip is not a deletion

        data = fetch_region(bam, "chr1:101-200", ref_path, min_vaf=0.1, min_depth=1)
        [snv] = [v for v in data.variants if v["position"] == 150]
        assert (snv["depth"], snv["alt_count"], snv["vaf"]) == (5, 1, 0.2)


class TestOverlappingMates:
    def _pair(self, tmp_path):
        """One fragment whose mates both cover 0-based 130 and both carry the alt."""
        alt = "C" if SEQ[130] != "C" else "G"
        r1 = SEQ[100:130] + alt + SEQ[131:150]
        r2 = SEQ[120:130] + alt + SEQ[131:170]
        reads = [
            {
                "name": "frag",
                "start": 100,
                "seq": r1,
                "cigar": [(0, 50)],
                "flag": 99,
                "mate": 120,
                "tlen": 70,
            },
            {
                "name": "frag",
                "start": 120,
                "seq": r2,
                "cigar": [(0, 50)],
                "flag": 147,
                "mate": 100,
                "tlen": -70,
            },
        ]
        return _bam(tmp_path, reads), alt

    @pytest.mark.unit
    def test_policy_is_both_mates_counted(self):
        assert OVERLAPPING_MATES_POLICY == "both_counted"

    @pytest.mark.unit
    @pytest.mark.parametrize("min_baseq", [0, 20])
    def test_vcf_support_depth_matches_count_coverage(self, tmp_path, min_baseq):
        bam, alt = self._pair(tmp_path)
        *_, counts = count_coverage_arrays(bam, "chr1:131-131", min_baseq=min_baseq)
        cc_depth = sum(int(c[0]) for c in counts)
        v = {"source": "vcf", "position": 130, "ref": SEQ[130], "alt": alt}
        annotate_vcf_snv_support(bam, "chr1:101-200", None, [v], min_baseq=min_baseq)
        assert cc_depth == 2
        assert v["read_depth"] == cc_depth
        assert v["strand_forward"] + v["strand_reverse"] == 2
        # True base qualities — pysam's overlap handling used to rewrite them to 80/0.
        assert v["_alt_qualities"] == [40, 40]

    @pytest.mark.unit
    def test_indel_depth_counts_both_mates(self, tmp_path):
        ref_path = _ref(tmp_path, SEQ)
        bam, _ = self._pair(tmp_path)
        data = fetch_region(bam, "chr1:101-200", ref_path, min_vaf=0.1, min_depth=1)
        [snv] = [v for v in data.variants if v["position"] == 130]
        assert snv["depth"] == 2


class TestDeletionsUnderMinBaseq:
    @pytest.mark.unit
    def test_deletion_reads_count_toward_depth_with_min_baseq(self, tmp_path):
        ref_path = _ref(tmp_path, SEQ)
        del_seq = SEQ[100:120] + SEQ[123:180]
        reads = [{"start": 100, "seq": del_seq, "cigar": [(0, 20), (2, 3), (0, 57)]}] * 4
        reads += [{"start": 100, "seq": SEQ[100:180], "cigar": [(0, 80)]}] * 4
        bam = _bam(tmp_path, reads)
        indels, gap_depth = detect_indels(
            bam, "chr1:101-200", ref_path, SEQ[100:200], min_baseq=20, min_vaf=0.1, min_depth=1
        )
        assert gap_depth[120 - 100 : 123 - 100] == [4, 4, 4]
        [d] = [v for v in indels if v["indel_type"] == "del"]
        assert (d["position"], d["depth"], d["alt_count"], d["vaf"]) == (119, 8, 4, 0.5)

    @pytest.mark.unit
    def test_low_quality_matched_bases_excluded_like_count_coverage(self, tmp_path):
        ref_path = _ref(tmp_path, SEQ)
        del_seq = SEQ[100:120] + SEQ[123:180]
        reads = [{"start": 100, "seq": del_seq, "cigar": [(0, 20), (2, 3), (0, 57)]}] * 2
        # Two reference reads with a low-quality base at the anchor column (0-based 119).
        lowq = "I" * 19 + "#" + "I" * 60
        reads += [{"start": 100, "seq": SEQ[100:180], "cigar": [(0, 80)], "qual": lowq}] * 2
        bam = _bam(tmp_path, reads)
        indels, _ = detect_indels(
            bam, "chr1:101-200", ref_path, SEQ[100:200], min_baseq=20, min_vaf=0.1, min_depth=1
        )
        [d] = [v for v in indels if v["indel_type"] == "del"]
        # Depth at the anchor = the 2 deletion reads' anchor bases (the low-quality ref
        # bases are dropped, as count_coverage drops them).
        assert d["depth"] == 2

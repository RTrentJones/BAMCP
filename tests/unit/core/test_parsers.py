"""Unit tests for bamcp.parsers module."""

import pytest

from bamcp.core.parsers import (
    AlignedRead,
    RegionData,
    SoftClip,
    detect_variants,
    extract_soft_clips,
    fetch_region,
    parse_region,
    scan_variants_chunked,
)


class TestParseRegion:
    """Tests for the parse_region function."""

    @pytest.mark.unit
    def test_standard_format(self):
        contig, start, end = parse_region("chr1:1000-2000")
        assert contig == "chr1"
        assert start == 1000
        assert end == 2000

    @pytest.mark.unit
    def test_no_chr_prefix(self):
        contig, start, end = parse_region("1:1000-2000")
        assert contig == "1"
        assert start == 1000
        assert end == 2000

    @pytest.mark.unit
    def test_commas_in_coordinates(self):
        contig, start, end = parse_region("chr1:1,000-2,000")
        assert contig == "chr1"
        assert start == 1000
        assert end == 2000

    @pytest.mark.unit
    def test_large_coordinates(self):
        contig, start, end = parse_region("chr17:7,577,000-7,577,500")
        assert contig == "chr17"
        assert start == 7577000
        assert end == 7577500

    @pytest.mark.unit
    def test_invalid_format_no_colon(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1_1000_2000")

    @pytest.mark.unit
    def test_invalid_format_no_dash(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1:1000")

    @pytest.mark.unit
    def test_invalid_format_non_numeric(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1:abc-def")

    @pytest.mark.unit
    def test_negative_start(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1:-10-100")

    @pytest.mark.unit
    def test_end_less_than_start(self):
        with pytest.raises(ValueError, match="must be greater than start"):
            parse_region("chr1:2000-1000")

    @pytest.mark.unit
    def test_end_equals_start(self):
        with pytest.raises(ValueError, match="must be greater than start"):
            parse_region("chr1:1000-1000")


class TestDetectVariants:
    """Tests for the detect_variants function."""

    @pytest.mark.unit
    def test_no_ref_seq(self):
        """No variants without reference sequence."""
        # A=10
        coverage_counts = ([10], [0], [0], [0])
        variants = detect_variants(coverage_counts, None, "chr1", 100)
        assert variants == []

    @pytest.mark.unit
    def test_no_variants(self):
        """No variants when all bases match reference."""
        # 4 positions: A, C, G, T
        cov_A = [20, 0, 0, 0]
        cov_C = [0, 20, 0, 0]
        cov_G = [0, 0, 20, 0]
        cov_T = [0, 0, 0, 20]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "ACGT"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100)
        assert variants == []

    @pytest.mark.unit
    def test_snv_detected(self):
        """SNV above VAF threshold should be detected."""
        # Position 0: ref=A, 15A + 5T = VAF 0.25
        cov_A = [15]
        cov_C = [0]
        cov_G = [0]
        cov_T = [5]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "A"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100, min_vaf=0.1, min_depth=10)
        assert len(variants) == 1
        assert variants[0]["position"] == 100
        assert variants[0]["ref"] == "A"
        assert variants[0]["alt"] == "T"
        assert variants[0]["vaf"] == 0.25
        assert variants[0]["depth"] == 20
        assert variants[0]["alt_count"] == 5

    @pytest.mark.unit
    def test_snv_below_vaf_threshold(self):
        """SNV below VAF threshold should not be detected."""
        # VAF = 0.05 (1/20)
        cov_A = [19]
        cov_C = [0]
        cov_G = [0]
        cov_T = [1]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "A"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100, min_vaf=0.1, min_depth=10)
        assert len(variants) == 0

    @pytest.mark.unit
    def test_low_depth_filtered(self):
        """Positions with low depth should be filtered."""
        # depth = 6
        cov_A = [3]
        cov_C = [0]
        cov_G = [0]
        cov_T = [3]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "A"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100, min_vaf=0.1, min_depth=10)
        assert len(variants) == 0

    @pytest.mark.unit
    def test_empty_counts(self):
        """Empty counts should be handled."""
        coverage_counts = ([], [], [], [])
        ref_seq = "A"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100)
        assert variants == []

    @pytest.mark.unit
    def test_multiple_alts(self):
        """Multiple alternate alleles at one position."""
        # A=10, T=5, G=5
        cov_A = [10]
        cov_C = [0]
        cov_G = [5]
        cov_T = [5]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "A"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 100, min_vaf=0.1, min_depth=10)
        assert len(variants) == 2
        alts = {v["alt"] for v in variants}
        assert alts == {"T", "G"}

    @pytest.mark.unit
    def test_correct_positions(self):
        """Variants should have correct genomic positions."""
        # Pos 0: A (ref=A)
        # Pos 1: C=10, G=10 (ref=C) -> Alt G
        # Pos 2: G (ref=G)
        cov_A = [20, 0, 0]
        cov_C = [0, 10, 0]
        cov_G = [0, 10, 20]
        cov_T = [0, 0, 0]
        coverage_counts = (cov_A, cov_C, cov_G, cov_T)

        ref_seq = "ACG"
        variants = detect_variants(coverage_counts, ref_seq, "chr1", 500, min_vaf=0.1, min_depth=10)
        assert len(variants) == 1
        assert variants[0]["position"] == 501
        assert variants[0]["ref"] == "C"
        assert variants[0]["alt"] == "G"


class TestFetchRegion:
    """Tests for the fetch_region function using real BAM fixtures."""

    @pytest.mark.unit
    def test_basic_fetch(self, small_bam_path):
        """Should fetch reads from a BAM file."""
        data = fetch_region(small_bam_path, "chr1:90-200")
        assert isinstance(data, RegionData)
        assert data.contig == "chr1"
        assert data.start == 90
        assert data.end == 200
        assert len(data.reads) > 0

    @pytest.mark.unit
    def test_reads_have_correct_fields(self, small_bam_path):
        """Reads should have all expected fields."""
        data = fetch_region(small_bam_path, "chr1:90-200")
        read = data.reads[0]
        assert isinstance(read, AlignedRead)
        assert isinstance(read.name, str)
        assert isinstance(read.sequence, str)
        assert isinstance(read.qualities, list)
        assert isinstance(read.cigar, str)
        assert isinstance(read.position, int)
        assert isinstance(read.end_position, int)
        assert isinstance(read.mapping_quality, int)
        assert isinstance(read.is_reverse, bool)
        assert isinstance(read.mismatches, list)

    @pytest.mark.unit
    def test_coverage_array(self, small_bam_path):
        """Coverage array should have correct length."""
        data = fetch_region(small_bam_path, "chr1:100-200")
        assert len(data.coverage) == 100  # end - start
        # Should have non-zero coverage where reads are
        assert sum(data.coverage) > 0

    @pytest.mark.unit
    def test_coverage_values(self, small_bam_path):
        """Coverage values should be non-negative integers."""
        data = fetch_region(small_bam_path, "chr1:100-200")
        for cov in data.coverage:
            assert isinstance(cov, int)
            assert cov >= 0

    @pytest.mark.unit
    def test_min_mapq_filter(self, small_bam_path):
        """Reads below min_mapq should be filtered."""
        data_no_filter = fetch_region(small_bam_path, "chr1:90-200", min_mapq=0)
        data_filtered = fetch_region(small_bam_path, "chr1:90-200", min_mapq=30)

        # read6_lowq has MAPQ 5, should be excluded with min_mapq=30
        names_no_filter = {r.name for r in data_no_filter.reads}
        names_filtered = {r.name for r in data_filtered.reads}

        assert "read6_lowq" in names_no_filter
        assert "read6_lowq" not in names_filtered

    @pytest.mark.unit
    def test_secondary_reads_filtered(self, small_bam_path):
        """Secondary alignments should be filtered."""
        data = fetch_region(small_bam_path, "chr1:390-450")
        names = {r.name for r in data.reads}
        assert "read7_secondary" not in names

    @pytest.mark.unit
    def test_secondary_reads_excluded_from_coverage(self, small_bam_path):
        """Secondary/supplementary reads should not inflate coverage counts."""
        data = fetch_region(small_bam_path, "chr1:390-450")
        # Coverage at any position should not exceed the number of primary reads
        # covering that position
        for i, cov in enumerate(data.coverage):
            pos = data.start + i
            primary_count = sum(1 for r in data.reads if r.position <= pos < r.end_position)
            assert cov <= primary_count + 1, (
                f"Coverage {cov} at pos {pos} exceeds primary read count {primary_count}"
            )

    @pytest.mark.unit
    def test_max_reads_limit(self, small_bam_path):
        """Should respect max_reads limit."""
        data = fetch_region(small_bam_path, "chr1:0-1000", max_reads=3)
        assert len(data.reads) <= 3

    @pytest.mark.unit
    def test_empty_region(self, small_bam_path):
        """Region with no reads should return empty lists."""
        data = fetch_region(small_bam_path, "chr1:900-950")
        assert len(data.reads) == 0
        assert all(c == 0 for c in data.coverage)

    @pytest.mark.unit
    def test_with_reference(self, small_bam_path, ref_fasta_path):
        """Should compute mismatches when reference is provided."""
        data = fetch_region(small_bam_path, "chr1:90-200", reference_path=ref_fasta_path)
        assert data.reference_sequence is not None
        assert len(data.reference_sequence) == 110  # 200 - 90

    @pytest.mark.unit
    def test_without_reference(self, small_bam_path):
        """Without reference, mismatches and variants should be empty."""
        data = fetch_region(small_bam_path, "chr1:90-200")
        assert data.reference_sequence is None
        assert data.variants == []

    @pytest.mark.unit
    def test_different_contigs(self, small_bam_path):
        """Should fetch reads from different contigs."""
        data = fetch_region(small_bam_path, "chr2:40-100")
        names = {r.name for r in data.reads}
        assert "read8_chr2" in names

    @pytest.mark.unit
    def test_forward_reverse_strand(self, small_bam_path):
        """Should correctly identify forward and reverse reads."""
        data = fetch_region(small_bam_path, "chr1:90-200")
        forward = [r for r in data.reads if not r.is_reverse]
        reverse = [r for r in data.reads if r.is_reverse]
        assert len(forward) > 0
        assert len(reverse) > 0

    @pytest.mark.unit
    def test_read_cigar(self, small_bam_path):
        """Reads should have CIGAR strings."""
        data = fetch_region(small_bam_path, "chr1:190-260")
        for read in data.reads:
            assert read.cigar != ""

    @pytest.mark.unit
    def test_insertion_cigar(self, small_bam_path):
        """Read with insertion should have correct CIGAR."""
        data = fetch_region(small_bam_path, "chr1:190-260")
        read4 = next((r for r in data.reads if r.name == "read4"), None)
        if read4:
            assert "I" in read4.cigar

    @pytest.mark.unit
    def test_deletion_cigar(self, small_bam_path):
        """Read with deletion should have correct CIGAR."""
        data = fetch_region(small_bam_path, "chr1:290-380")
        read5 = next((r for r in data.reads if r.name == "read5"), None)
        if read5:
            assert "D" in read5.cigar

    @pytest.mark.unit
    def test_invalid_region_format(self, small_bam_path):
        """Should raise ValueError for invalid region."""
        with pytest.raises(ValueError):
            fetch_region(small_bam_path, "invalid_region")

    @pytest.mark.unit
    def test_nonexistent_file(self):
        """Should raise error for nonexistent BAM file."""
        with pytest.raises((FileNotFoundError, OSError)):
            fetch_region("/nonexistent/file.bam", "chr1:100-200")

    @pytest.mark.unit
    def test_empty_bam(self, empty_bam_path):
        """Should return empty data for BAM with no reads."""
        data = fetch_region(empty_bam_path, "chr1:100-200")
        assert len(data.reads) == 0
        assert all(c == 0 for c in data.coverage)


class TestAlignedRead:
    """Tests for the AlignedRead dataclass."""

    @pytest.mark.unit
    def test_default_mismatches(self):
        """Mismatches should default to empty list."""
        read = AlignedRead(
            name="test",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
        )
        assert read.mismatches == []

    @pytest.mark.unit
    def test_all_fields(self):
        """All fields should be accessible."""
        read = AlignedRead(
            name="read1",
            sequence="ACGT",
            qualities=[30, 31, 32, 33],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=True,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        assert read.name == "read1"
        assert read.sequence == "ACGT"
        assert read.qualities == [30, 31, 32, 33]
        assert read.cigar == "4M"
        assert read.position == 100
        assert read.end_position == 104
        assert read.mapping_quality == 60
        assert read.is_reverse is True
        assert len(read.mismatches) == 1


class TestRegionData:
    """Tests for the RegionData dataclass."""

    @pytest.mark.unit
    def test_region_data_creation(self):
        """RegionData should store all fields."""
        data = RegionData(
            contig="chr1",
            start=100,
            end=200,
            reads=[],
            coverage=[0] * 100,
            variants=[],
            reference_sequence="ACGT" * 25,
        )
        assert data.contig == "chr1"
        assert data.start == 100
        assert data.end == 200
        assert len(data.coverage) == 100
        assert data.reference_sequence is not None

    @pytest.mark.unit
    def test_region_data_optional_ref(self):
        """reference_sequence should default to None."""
        data = RegionData(
            contig="chr1",
            start=100,
            end=200,
            reads=[],
            coverage=[],
            variants=[],
        )
        assert data.reference_sequence is None


class TestSoftClip:
    """Tests for the SoftClip dataclass."""

    @pytest.mark.unit
    def test_soft_clip_creation(self):
        """SoftClip should store all fields."""
        clip = SoftClip(
            position=100,
            length=10,
            sequence="ACGTACGTAC",
            side="left",
        )
        assert clip.position == 100
        assert clip.length == 10
        assert clip.sequence == "ACGTACGTAC"
        assert clip.side == "left"

    @pytest.mark.unit
    def test_soft_clip_optional_sequence(self):
        """sequence can be None."""
        clip = SoftClip(
            position=100,
            length=10,
            sequence=None,
            side="right",
        )
        assert clip.sequence is None


class TestExtractSoftClips:
    """Tests for the extract_soft_clips function."""

    @pytest.mark.unit
    def test_no_soft_clips(self, small_bam_path):
        """Reads without soft clips should return empty list."""
        import pysam

        with pysam.AlignmentFile(small_bam_path, "rb") as samfile:
            for read in samfile.fetch("chr1", 100, 200):
                if read.cigarstring and "S" not in read.cigarstring:
                    clips = extract_soft_clips(read)
                    assert clips == []
                    return
        # If no reads without S, skip test
        pytest.skip("No reads without soft clips found")

    @pytest.mark.unit
    def test_aligned_read_has_soft_clips_field(self):
        """AlignedRead should have soft_clips field."""
        read = AlignedRead(
            name="test",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
        )
        assert hasattr(read, "soft_clips")
        assert read.soft_clips == []

    @pytest.mark.unit
    def test_aligned_read_with_soft_clips(self):
        """AlignedRead should store soft_clips."""
        clips = [
            SoftClip(position=95, length=5, sequence="NNNNN", side="left"),
        ]
        read = AlignedRead(
            name="test",
            sequence="NNNNNACGT",
            qualities=[30] * 9,
            cigar="5S4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            soft_clips=clips,
        )
        assert len(read.soft_clips) == 1
        assert read.soft_clips[0].side == "left"
        assert read.soft_clips[0].length == 5


class TestScanVariantsChunked:
    """Tests for the scan_variants_chunked function."""

    @pytest.mark.unit
    def test_scans_contig(self, small_bam_path, ref_fasta_path):
        """Should find variants across chr1."""
        variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, chunk_size=500
        )
        assert isinstance(variants, list)
        # All variants should be on chr1
        for v in variants:
            assert v["contig"] == "chr1"
            assert "position" in v
            assert "ref" in v
            assert "alt" in v
            assert "vaf" in v

    @pytest.mark.unit
    def test_sorted_by_vaf(self, small_bam_path, ref_fasta_path):
        """Results should be sorted by VAF descending."""
        variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, chunk_size=500
        )
        if len(variants) >= 2:
            for i in range(len(variants) - 1):
                assert variants[i]["vaf"] >= variants[i + 1]["vaf"]

    @pytest.mark.unit
    def test_nonexistent_contig(self, small_bam_path, ref_fasta_path):
        """Should return empty list for nonexistent contig."""
        variants = scan_variants_chunked(
            small_bam_path, "chrZ", ref_fasta_path
        )
        assert variants == []

    @pytest.mark.unit
    def test_respects_min_vaf(self, small_bam_path, ref_fasta_path):
        """Higher min_vaf should return fewer or equal variants."""
        all_variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, min_vaf=0.05
        )
        strict_variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, min_vaf=0.5
        )
        assert len(strict_variants) <= len(all_variants)

    @pytest.mark.unit
    def test_respects_min_depth(self, small_bam_path, ref_fasta_path):
        """Higher min_depth should return fewer or equal variants."""
        all_variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, min_depth=1
        )
        strict_variants = scan_variants_chunked(
            small_bam_path, "chr1", ref_fasta_path, min_depth=10
        )
        assert len(strict_variants) <= len(all_variants)

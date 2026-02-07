"""Unit tests for bamcp.reference module."""

import pytest

from bamcp.reference import (
    GENOME_BUILDS,
    detect_genome_build,
    get_public_reference_url,
    normalize_build_name,
)


class TestNormalizeBuildName:
    """Tests for normalize_build_name function."""

    @pytest.mark.unit
    def test_canonical_grch38(self):
        """GRCh38 should normalize to itself."""
        assert normalize_build_name("GRCh38") == "GRCh38"

    @pytest.mark.unit
    def test_canonical_grch37(self):
        """GRCh37 should normalize to itself."""
        assert normalize_build_name("GRCh37") == "GRCh37"

    @pytest.mark.unit
    def test_hg38_alias(self):
        """hg38 should normalize to GRCh38."""
        assert normalize_build_name("hg38") == "GRCh38"

    @pytest.mark.unit
    def test_hg19_alias(self):
        """hg19 should normalize to GRCh37."""
        assert normalize_build_name("hg19") == "GRCh37"

    @pytest.mark.unit
    def test_b37_alias(self):
        """b37 should normalize to GRCh37."""
        assert normalize_build_name("b37") == "GRCh37"

    @pytest.mark.unit
    def test_case_insensitive(self):
        """Normalization should be case insensitive."""
        assert normalize_build_name("GRCH38") == "GRCh38"
        assert normalize_build_name("HG38") == "GRCh38"
        assert normalize_build_name("grch37") == "GRCh37"

    @pytest.mark.unit
    def test_unknown_build(self):
        """Unknown builds should return None."""
        assert normalize_build_name("unknown") is None
        assert normalize_build_name("mm10") is None
        assert normalize_build_name("") is None


class TestGetPublicReferenceUrl:
    """Tests for get_public_reference_url function."""

    @pytest.mark.unit
    def test_grch38_url(self):
        """GRCh38 should return UCSC hg38 URL."""
        url = get_public_reference_url("GRCh38")
        assert url is not None
        assert "hg38" in url
        assert url.startswith("https://")

    @pytest.mark.unit
    def test_grch37_url(self):
        """GRCh37 should return UCSC hg19 URL."""
        url = get_public_reference_url("GRCh37")
        assert url is not None
        assert "hg19" in url
        assert url.startswith("https://")

    @pytest.mark.unit
    def test_alias_resolution(self):
        """Aliases should work for URL lookup."""
        assert get_public_reference_url("hg38") == get_public_reference_url("GRCh38")
        assert get_public_reference_url("hg19") == get_public_reference_url("GRCh37")

    @pytest.mark.unit
    def test_unknown_build_returns_none(self):
        """Unknown builds should return None."""
        assert get_public_reference_url("unknown") is None
        assert get_public_reference_url("mm10") is None


class TestDetectGenomeBuild:
    """Tests for detect_genome_build function."""

    @pytest.mark.unit
    def test_grch38_chr1_length(self):
        """GRCh38 should be detected from chr1 length."""
        contigs = [
            {"name": "chr1", "length": 248956422},
            {"name": "chr2", "length": 242193529},
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "GRCh38"
        assert result["confidence"] == "high"
        assert len(result["evidence"]) > 0
        assert "248,956,422" in result["evidence"][0]

    @pytest.mark.unit
    def test_grch37_chr1_length(self):
        """GRCh37 should be detected from chr1 length."""
        contigs = [
            {"name": "chr1", "length": 249250621},
            {"name": "chr2", "length": 243199373},
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "GRCh37"
        assert result["confidence"] == "high"

    @pytest.mark.unit
    def test_ncbi_style_naming(self):
        """Detection should work with NCBI-style naming (no chr prefix)."""
        contigs = [
            {"name": "1", "length": 248956422},
            {"name": "2", "length": 242193529},
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "GRCh38"
        assert result["confidence"] == "high"

    @pytest.mark.unit
    def test_empty_contigs(self):
        """Empty contigs should return unknown with low confidence."""
        result = detect_genome_build([])

        assert result["build"] == "unknown"
        assert result["confidence"] == "low"

    @pytest.mark.unit
    def test_no_chr1(self):
        """Missing chr1 should return unknown."""
        contigs = [
            {"name": "chr2", "length": 242193529},
            {"name": "chr3", "length": 198295559},
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "unknown"
        assert result["confidence"] == "low"
        assert any("chr1" in e for e in result["evidence"])

    @pytest.mark.unit
    def test_unknown_chr1_length(self):
        """Unknown chr1 length should return unknown with medium confidence."""
        contigs = [
            {"name": "chr1", "length": 100000000},  # Not a known human build
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "unknown"
        assert result["confidence"] == "medium"

    @pytest.mark.unit
    def test_small_test_bam(self):
        """Small test BAM with arbitrary lengths should return unknown."""
        # This matches the test fixtures in tests/fixtures/
        contigs = [
            {"name": "chr1", "length": 1000},
            {"name": "chr2", "length": 500},
        ]
        result = detect_genome_build(contigs)

        assert result["build"] == "unknown"
        assert result["confidence"] == "medium"


class TestGenomeBuildsRegistry:
    """Tests for the GENOME_BUILDS registry."""

    @pytest.mark.unit
    def test_grch38_has_required_fields(self):
        """GRCh38 entry should have all required fields."""
        assert "GRCh38" in GENOME_BUILDS
        grch38 = GENOME_BUILDS["GRCh38"]

        assert "aliases" in grch38
        assert "chr1_length" in grch38
        assert "fasta_url" in grch38
        assert grch38["chr1_length"] == 248956422

    @pytest.mark.unit
    def test_grch37_has_required_fields(self):
        """GRCh37 entry should have all required fields."""
        assert "GRCh37" in GENOME_BUILDS
        grch37 = GENOME_BUILDS["GRCh37"]

        assert "aliases" in grch37
        assert "chr1_length" in grch37
        assert "fasta_url" in grch37
        assert grch37["chr1_length"] == 249250621

    @pytest.mark.unit
    def test_chr1_lengths_differ(self):
        """GRCh37 and GRCh38 chr1 lengths should be different."""
        grch37_len = GENOME_BUILDS["GRCh37"]["chr1_length"]
        grch38_len = GENOME_BUILDS["GRCh38"]["chr1_length"]

        assert grch37_len != grch38_len
        assert grch37_len > grch38_len  # GRCh37 is longer

    @pytest.mark.unit
    def test_fasta_urls_are_https(self):
        """All FASTA URLs should use HTTPS."""
        for build_name, info in GENOME_BUILDS.items():
            assert info["fasta_url"].startswith("https://"), f"{build_name} URL not HTTPS"

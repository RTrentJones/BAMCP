"""Unit tests for security validation."""

import pytest

from bamcp.config import BAMCPConfig
from bamcp.tools import ALLELE_PATTERN, CHROM_PATTERN, validate_lookup_inputs, validate_path


class TestValidatePath:
    """Tests for validate_path function."""

    @pytest.mark.unit
    def test_remote_files_disabled_by_default(self):
        config = BAMCPConfig(allow_remote_files=False)
        with pytest.raises(ValueError, match="Remote files are disabled"):
            validate_path("https://example.com/file.bam", config)

    @pytest.mark.unit
    def test_remote_files_allowed(self):
        config = BAMCPConfig(allow_remote_files=True)
        # Should not raise
        validate_path("https://example.com/file.bam", config)

    @pytest.mark.unit
    def test_remote_files_invalid_scheme(self):
        config = BAMCPConfig(allow_remote_files=True)
        with pytest.raises(ValueError, match="Scheme not supported for remote file"):
            validate_path("ftp://example.com/file.bam", config)

    @pytest.mark.unit
    def test_local_files_no_restriction(self, tmp_path):
        config = BAMCPConfig(allowed_directories=None)
        f = tmp_path / "test.bam"
        f.touch()
        # Should not raise
        validate_path(str(f), config)

    @pytest.mark.unit
    def test_local_files_allowed_directory(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        config = BAMCPConfig(allowed_directories=[str(allowed)])

        f = allowed / "test.bam"
        f.touch()

        # Should not raise
        validate_path(str(f), config)

    @pytest.mark.unit
    def test_local_files_denied_directory(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        forbidden = tmp_path / "forbidden"
        forbidden.mkdir()

        config = BAMCPConfig(allowed_directories=[str(allowed)])

        f = forbidden / "test.bam"
        f.touch()

        with pytest.raises(ValueError, match="Path .* is not in allowed directories"):
            validate_path(str(f), config)

    @pytest.mark.unit
    def test_path_traversal_attempt(self, tmp_path):
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        secret = tmp_path / "secret.txt"
        secret.touch()

        config = BAMCPConfig(allowed_directories=[str(allowed)])

        # Try to break out of allowed directory using ..
        traversal = allowed / "../secret.txt"

        with pytest.raises(ValueError, match="Path .* is not in allowed directories"):
            validate_path(str(traversal), config)


class TestChromosomePattern:
    """Tests for chromosome validation pattern."""

    @pytest.mark.unit
    def test_valid_autosomal_chromosomes(self):
        """Standard numbered chromosomes should be valid."""
        for num in range(1, 23):
            assert CHROM_PATTERN.match(str(num))
            assert CHROM_PATTERN.match(f"chr{num}")

    @pytest.mark.unit
    def test_valid_sex_chromosomes(self):
        """X, Y chromosomes should be valid."""
        for chrom in ["X", "Y", "chrX", "chrY"]:
            assert CHROM_PATTERN.match(chrom)

    @pytest.mark.unit
    def test_valid_mitochondrial(self):
        """MT and M should be valid."""
        for chrom in ["M", "MT", "chrM", "chrMT"]:
            assert CHROM_PATTERN.match(chrom)

    @pytest.mark.unit
    def test_case_insensitive(self):
        """Chromosome pattern should be case insensitive."""
        assert CHROM_PATTERN.match("chr1")
        assert CHROM_PATTERN.match("CHR1")
        assert CHROM_PATTERN.match("chrX")
        assert CHROM_PATTERN.match("CHRX")

    @pytest.mark.unit
    def test_invalid_chromosome_names(self):
        """Invalid chromosome names should not match."""
        invalid = ["chrABC", "chr1A", "chromosome1", "chr_1", ""]
        for chrom in invalid:
            assert not CHROM_PATTERN.match(chrom), f"{chrom} should be invalid"


class TestAllelePattern:
    """Tests for allele validation pattern."""

    @pytest.mark.unit
    def test_valid_single_nucleotide(self):
        """Single nucleotides should be valid."""
        for base in "ACGTN":
            assert ALLELE_PATTERN.match(base)

    @pytest.mark.unit
    def test_valid_multiple_nucleotides(self):
        """Multiple nucleotides (insertions, deletions) should be valid."""
        assert ALLELE_PATTERN.match("ACGT")
        assert ALLELE_PATTERN.match("AAAAAAA")
        assert ALLELE_PATTERN.match("TGCATGCA")

    @pytest.mark.unit
    def test_case_insensitive(self):
        """Allele pattern should be case insensitive."""
        assert ALLELE_PATTERN.match("acgt")
        assert ALLELE_PATTERN.match("ACGT")
        assert ALLELE_PATTERN.match("AcGt")

    @pytest.mark.unit
    def test_invalid_alleles(self):
        """Invalid allele sequences should not match."""
        invalid = ["", "X", "U", "ACGU", "A-T", "A T", "123"]
        for allele in invalid:
            assert not ALLELE_PATTERN.match(allele), f"{allele} should be invalid"


class TestValidateLookupInputs:
    """Tests for validate_lookup_inputs function."""

    @pytest.mark.unit
    def test_valid_input(self):
        """Valid chromosome, position, ref, alt should pass validation."""
        # Should not raise
        validate_lookup_inputs("chr17", 7674220, "G", "A")

    @pytest.mark.unit
    def test_valid_input_without_chr_prefix(self):
        """Chromosome without chr prefix should be valid."""
        validate_lookup_inputs("17", 7674220, "G", "A")

    @pytest.mark.unit
    def test_invalid_chromosome(self):
        """Invalid chromosome should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid chromosome"):
            validate_lookup_inputs("chrABC", 100, "A", "T")

    @pytest.mark.unit
    def test_negative_position(self):
        """Negative position should raise ValueError."""
        with pytest.raises(ValueError, match="Position must be positive"):
            validate_lookup_inputs("chr1", -1, "A", "T")

    @pytest.mark.unit
    def test_zero_position(self):
        """Zero position should raise ValueError."""
        with pytest.raises(ValueError, match="Position must be positive"):
            validate_lookup_inputs("chr1", 0, "A", "T")

    @pytest.mark.unit
    def test_invalid_ref_allele(self):
        """Invalid reference allele should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid reference allele"):
            validate_lookup_inputs("chr1", 100, "X", "A")

    @pytest.mark.unit
    def test_invalid_alt_allele(self):
        """Invalid alternate allele should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid alternate allele"):
            validate_lookup_inputs("chr1", 100, "A", "X")

    @pytest.mark.unit
    def test_empty_ref_allele(self):
        """Empty reference allele should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid reference allele"):
            validate_lookup_inputs("chr1", 100, "", "A")

    @pytest.mark.unit
    def test_empty_alt_allele(self):
        """Empty alternate allele should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid alternate allele"):
            validate_lookup_inputs("chr1", 100, "A", "")

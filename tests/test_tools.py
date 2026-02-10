"""Unit tests for bamcp.tools module."""

import json

import pytest

from bamcp.config import BAMCPConfig
from bamcp.constants import (
    DEFAULT_CONTIG,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
    LOW_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MIN_MEAN_QUALITY,
    LOW_CONFIDENCE_MIN_VAF,
    VIEWER_RESOURCE_URI,
)
from bamcp.parsers import AlignedRead, RegionData
from bamcp.tools import (
    _serialize_region_data,
    handle_browse_region,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_visualize_region,
)


@pytest.fixture
def config():
    """Default test config."""
    return BAMCPConfig()


@pytest.fixture
def config_with_ref(ref_fasta_path):
    """Config with reference path."""
    return BAMCPConfig(reference=ref_fasta_path)


class TestSerializeRegionData:
    """Tests for _serialize_region_data helper."""

    @pytest.mark.unit
    def test_empty_region(self):
        data = RegionData(
            contig="chr1", start=100, end=200, reads=[], coverage=[0] * 100, variants=[]
        )
        result = _serialize_region_data(data)
        assert result["contig"] == "chr1"
        assert result["start"] == 100
        assert result["end"] == 200
        assert result["reads"] == []
        assert len(result["coverage"]) == 100
        assert result["variants"] == []
        assert result["reference_sequence"] is None

    @pytest.mark.unit
    def test_with_reads(self):
        read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=200,
            reads=[read],
            coverage=[1] * 100,
            variants=[],
        )
        # Default is compact=True, which omits sequence
        result = _serialize_region_data(data)
        assert len(result["reads"]) == 1
        r = result["reads"][0]
        assert r["name"] == "r1"
        assert "sequence" not in r  # Compact mode omits sequence
        assert r["cigar"] == "4M"
        assert r["position"] == 100
        assert r["end_position"] == 104
        assert r["mapping_quality"] == 60
        assert "qualities" not in r  # Qualities are never serialized
        assert r["is_reverse"] is False
        assert len(r["mismatches"]) == 1

    @pytest.mark.unit
    def test_non_compact_includes_sequence(self):
        """Non-compact mode includes sequence but never qualities."""
        read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=200,
            reads=[read],
            coverage=[1] * 100,
            variants=[],
        )
        result = _serialize_region_data(data, compact=False)
        r = result["reads"][0]
        assert r["sequence"] == "ACGT"
        assert "qualities" not in r  # Qualities never serialized

    @pytest.mark.unit
    def test_serialization_is_json_compatible(self):
        data = RegionData(
            contig="chr1",
            start=0,
            end=10,
            reads=[],
            coverage=[0] * 10,
            variants=[
                {
                    "contig": "chr1",
                    "position": 5,
                    "ref": "A",
                    "alt": "T",
                    "vaf": 0.5,
                    "depth": 20,
                    "alt_count": 10,
                }
            ],
            reference_sequence="ACGTACGTAC",
        )
        result = _serialize_region_data(data)
        # Should be JSON serializable
        json_str = json.dumps(result)
        parsed = json.loads(json_str)
        assert parsed == result


class TestHandleBrowseRegion:
    """Tests for handle_browse_region tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_ui_metadata(self, small_bam_path, config):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        assert "content" in result
        assert "_meta" in result
        assert result["_meta"]["ui/resourceUri"] == VIEWER_RESOURCE_URI
        assert "ui/init" in result["_meta"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_content_is_summary(self, small_bam_path, config):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        # Content is now a summary string, not JSON payload
        text = result["content"][0]["text"]
        assert "Region chr1:90-200" in text
        assert "reads" in text
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert "contig" in payload
        assert "reads" in payload
        assert "coverage" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_with_reference(self, small_bam_path, ref_fasta_path, config):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200", "reference": ref_fasta_path},
            config,
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["reference_sequence"] is not None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_config_reference_fallback(self, small_bam_path, config_with_ref):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["reference_sequence"] is not None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_invalid_region(self, small_bam_path, config):
        with pytest.raises(ValueError):
            await handle_browse_region({"file_path": small_bam_path, "region": "invalid"}, config)


class TestHandleGetVariants:
    """Tests for handle_get_variants tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_variants_structure(self, small_bam_path, config):
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        text = result["content"][0]["text"]
        payload = json.loads(text)
        assert "variants" in payload
        assert "count" in payload
        assert isinstance(payload["variants"], list)
        assert payload["count"] == len(payload["variants"])

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_custom_vaf_filter(self, small_bam_path, config_with_ref):
        result = await handle_get_variants(
            {
                "file_path": small_bam_path,
                "region": "chr1:90-200",
                "min_vaf": 0.5,
                "min_depth": 1,
            },
            config_with_ref,
        )
        payload = json.loads(result["content"][0]["text"])
        for v in payload["variants"]:
            assert v["vaf"] >= 0.5

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_no_meta_key(self, small_bam_path, config):
        """get_variants should not include UI metadata."""
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        assert "_meta" not in result


class TestHandleGetCoverage:
    """Tests for handle_get_coverage tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_coverage_stats(self, small_bam_path, config):
        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert "region" in payload
        assert "mean" in payload
        assert "min" in payload
        assert "max" in payload
        assert "median" in payload
        assert "bases_covered" in payload
        assert "total_bases" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_coverage_values(self, small_bam_path, config):
        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["total_bases"] == 110
        assert payload["min"] >= 0
        assert payload["max"] >= payload["min"]
        assert payload["mean"] >= 0

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_empty_region_coverage(self, small_bam_path, config):
        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:900-950"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["mean"] == 0
        assert payload["max"] == 0
        assert payload["bases_covered"] == 0

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_region_format_in_output(self, small_bam_path, config):
        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["region"] == "chr1:100-200"


class TestHandleListContigs:
    """Tests for handle_list_contigs tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_contigs(self, small_bam_path, config):
        result = await handle_list_contigs({"file_path": small_bam_path}, config)
        payload = json.loads(result["content"][0]["text"])
        assert "contigs" in payload
        contigs = payload["contigs"]
        assert len(contigs) == 2

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_contig_format(self, small_bam_path, config):
        result = await handle_list_contigs({"file_path": small_bam_path}, config)
        payload = json.loads(result["content"][0]["text"])
        for contig in payload["contigs"]:
            assert "name" in contig
            assert "length" in contig
            assert isinstance(contig["name"], str)
            assert isinstance(contig["length"], int)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_expected_contigs(self, small_bam_path, config):
        result = await handle_list_contigs({"file_path": small_bam_path}, config)
        payload = json.loads(result["content"][0]["text"])
        names = {c["name"] for c in payload["contigs"]}
        assert "chr1" in names
        assert "chr2" in names

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_contig_lengths(self, small_bam_path, config):
        result = await handle_list_contigs({"file_path": small_bam_path}, config)
        payload = json.loads(result["content"][0]["text"])
        contig_map = {c["name"]: c["length"] for c in payload["contigs"]}
        assert contig_map["chr1"] == 1000
        assert contig_map["chr2"] == 500

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_nonexistent_file(self, config):
        with pytest.raises((FileNotFoundError, OSError)):
            await handle_list_contigs({"file_path": "/nonexistent.bam"}, config)


class TestHandleJumpTo:
    """Tests for handle_jump_to tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_ui_metadata(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1"}, config
        )
        assert "_meta" in result
        assert result["_meta"]["ui/resourceUri"] == VIEWER_RESOURCE_URI

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_centers_on_position(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1"}, config
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["contig"] == DEFAULT_CONTIG
        # Position 150 should be within the returned start-end range
        assert payload["start"] <= 150 <= payload["end"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_custom_window(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 100},
            config,
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        span = payload["end"] - payload["start"]
        assert span == 100

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_window_from_config(self, small_bam_path):
        config = BAMCPConfig(default_window=200)
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 300, "contig": "chr1"}, config
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        span = payload["end"] - payload["start"]
        assert span == 200

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_contig(self, small_bam_path, config):
        """Without contig arg, defaults to configured constant."""
        result = await handle_jump_to({"file_path": small_bam_path, "position": 150}, config)
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["contig"] == DEFAULT_CONTIG

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_reads_returned(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 200},
            config,
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert "reads" in payload
        assert len(payload["reads"]) > 0


class TestHandleVisualizeRegion:
    """Tests for handle_visualize_region tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_ui_metadata(self, small_bam_path, config):
        result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        assert "content" in result
        assert "_meta" in result
        assert result["_meta"]["ui/resourceUri"] == VIEWER_RESOURCE_URI
        assert "ui/init" in result["_meta"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_content_is_summary(self, small_bam_path, config):
        result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        # Content is now a summary string, not JSON payload
        text = result["content"][0]["text"]
        assert "Region chr1:90-200" in text
        assert "reads" in text
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert "contig" in payload
        assert "reads" in payload
        assert "coverage" in payload


class TestHandleGetRegionSummary:
    """Tests for handle_get_region_summary tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_text_summary(self, small_bam_path, config):
        result = await handle_get_region_summary(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        assert "content" in result
        assert "_meta" not in result
        text = result["content"][0]["text"]
        assert "Region:" in text
        assert "Reads:" in text
        assert "Coverage:" in text
        assert "Variants detected:" in text

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_summary_includes_stats(self, small_bam_path, config):
        result = await handle_get_region_summary(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        text = result["content"][0]["text"]
        assert "chr1:90-200" in text
        assert "mean=" in text
        assert "max=" in text


class TestHandleLookupClinvar:
    """Tests for handle_lookup_clinvar tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_annotation(self, config, monkeypatch):
        """Mock ClinVarClient to return a known result."""
        from bamcp.clinvar import ClinVarResult

        async def mock_lookup(self, chrom, pos, ref, alt):
            return ClinVarResult(
                variation_id=12345,
                clinical_significance="Pathogenic",
                review_status="criteria provided, multiple submitters, no conflicts",
                stars=2,
                conditions=["Li-Fraumeni syndrome"],
                last_evaluated="2023-01-15",
                gene="TP53",
                variant_name="NM_000546.6(TP53):c.743G>A",
            )

        from bamcp import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", mock_lookup)

        result = await handle_lookup_clinvar(
            {"chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["clinical_significance"] == "Pathogenic"
        assert payload["gene"] == "TP53"
        assert payload["stars"] == 2
        assert "disclaimer" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_includes_disclaimer(self, config, monkeypatch):
        from bamcp.clinvar import ClinVarResult

        async def mock_lookup(self, chrom, pos, ref, alt):
            return ClinVarResult(
                variation_id=1,
                clinical_significance="Benign",
                review_status="no assertion criteria provided",
                stars=0,
                conditions=[],
                last_evaluated=None,
                gene=None,
                variant_name=None,
            )

        from bamcp import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", mock_lookup)

        result = await handle_lookup_clinvar(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert (
            "not for clinical" in payload["disclaimer"].lower()
            or "not intended" in payload["disclaimer"].lower()
        )

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_not_found(self, config, monkeypatch):
        async def mock_lookup(self, chrom, pos, ref, alt):
            return None

        from bamcp import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", mock_lookup)

        result = await handle_lookup_clinvar(
            {"chrom": "chr99", "pos": 1, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["found"] is False
        assert "disclaimer" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_error_handling(self, config, monkeypatch):
        async def mock_lookup(self, chrom, pos, ref, alt):
            raise ConnectionError("API unavailable")

        from bamcp import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", mock_lookup)

        result = await handle_lookup_clinvar(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert "error" in payload


class TestHandleLookupGnomad:
    """Tests for handle_lookup_gnomad tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_frequency(self, config, monkeypatch):
        from bamcp.gnomad import GnomadResult, PopulationFrequency

        async def mock_lookup(self, chrom, pos, ref, alt):
            return GnomadResult(
                variant_id="17-7674220-G-A",
                global_af=3.28e-05,
                ac=5,
                an=152312,
                homozygote_count=0,
                populations=[
                    PopulationFrequency(id="afr", ac=2, an=41406, homozygote_count=0, af=4.83e-05),
                ],
                filters=["PASS"],
                source="genome",
            )

        from bamcp import gnomad

        monkeypatch.setattr(gnomad.GnomadClient, "lookup", mock_lookup)

        result = await handle_lookup_gnomad(
            {"chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["global_af"] == pytest.approx(3.28e-05)
        assert payload["source"] == "genome"
        assert "disclaimer" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_not_found(self, config, monkeypatch):
        async def mock_lookup(self, chrom, pos, ref, alt):
            return None

        from bamcp import gnomad

        monkeypatch.setattr(gnomad.GnomadClient, "lookup", mock_lookup)

        result = await handle_lookup_gnomad(
            {"chrom": "chr99", "pos": 1, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["found"] is False
        assert "disclaimer" in payload


class TestLowConfidenceThresholds:
    """Tests for low-confidence variant labeling thresholds."""

    @pytest.mark.unit
    def test_marks_low_confidence_using_shared_thresholds(self):
        read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[LOW_CONFIDENCE_MIN_MEAN_QUALITY - 1] * 4,
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=120,
            reads=[read],
            coverage=[1] * 20,
            variants=[
                {
                    "contig": "chr1",
                    "position": 101,
                    "ref": "C",
                    "alt": "T",
                    "vaf": LOW_CONFIDENCE_MIN_VAF - 0.01,
                    "depth": LOW_CONFIDENCE_MIN_DEPTH - 1,
                    "alt_count": 1,
                }
            ],
        )

        result = _serialize_region_data(data)

        assert result["variants"][0]["is_low_confidence"] is True

    @pytest.mark.unit
    def test_keeps_high_confidence_when_all_thresholds_pass(self):
        forward_read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[LOW_CONFIDENCE_MIN_MEAN_QUALITY + 10] * 4,
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        reverse_read = AlignedRead(
            name="r2",
            sequence="ACGT",
            qualities=[LOW_CONFIDENCE_MIN_MEAN_QUALITY + 10] * 4,
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=True,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=120,
            reads=[forward_read, reverse_read],
            coverage=[2] * 20,
            variants=[
                {
                    "contig": "chr1",
                    "position": 101,
                    "ref": "C",
                    "alt": "T",
                    "vaf": LOW_CONFIDENCE_MIN_VAF + 0.1,
                    "depth": LOW_CONFIDENCE_MIN_DEPTH + 10,
                    "alt_count": 2,
                }
            ],
        )

        result = _serialize_region_data(data)

        assert (
            result["variant_evidence"]["101:C>T"]["strand_bias"] <= LOW_CONFIDENCE_MAX_STRAND_BIAS
        )
        assert result["variants"][0]["is_low_confidence"] is False


class TestBinValues:
    """Tests for _bin_values helper function."""

    @pytest.mark.unit
    def test_empty_values(self):
        from bamcp.tools import _bin_values

        result = _bin_values([], [0, 10, 20, 30])
        assert result == [0, 0, 0, 0]

    @pytest.mark.unit
    def test_single_bin(self):
        from bamcp.tools import _bin_values

        result = _bin_values([5, 15, 25], [0, 10, 20, 30])
        assert result == [1, 1, 1, 0]

    @pytest.mark.unit
    def test_all_in_last_bin(self):
        from bamcp.tools import _bin_values

        result = _bin_values([35, 40, 100], [0, 10, 20, 30])
        assert result == [0, 0, 0, 3]

    @pytest.mark.unit
    def test_quality_histogram_bins(self):
        from bamcp.constants import QUALITY_HISTOGRAM_BINS
        from bamcp.tools import _bin_values

        values = [5, 15, 25, 35, 45, 50]
        result = _bin_values(values, QUALITY_HISTOGRAM_BINS)
        assert len(result) == 5
        assert result == [1, 1, 1, 1, 2]  # 45 and 50 go in last bin


class TestDetectHomopolymer:
    """Tests for _detect_homopolymer helper function."""

    @pytest.mark.unit
    def test_single_base(self):
        from bamcp.tools import _detect_homopolymer

        result = _detect_homopolymer("ACGT", 0)
        assert result == 1

    @pytest.mark.unit
    def test_homopolymer_run(self):
        from bamcp.tools import _detect_homopolymer

        result = _detect_homopolymer("AAAAAACGT", 3)
        assert result == 6

    @pytest.mark.unit
    def test_position_at_end(self):
        from bamcp.tools import _detect_homopolymer

        result = _detect_homopolymer("CGTTTTT", 5)
        assert result == 5

    @pytest.mark.unit
    def test_empty_sequence(self):
        from bamcp.tools import _detect_homopolymer

        result = _detect_homopolymer("", 0)
        assert result == 0

    @pytest.mark.unit
    def test_position_out_of_range(self):
        from bamcp.tools import _detect_homopolymer

        result = _detect_homopolymer("ACGT", 10)
        assert result == 0


class TestComputeArtifactRisk:
    """Tests for _compute_artifact_risk function."""

    @pytest.mark.unit
    def test_low_risk(self):
        from bamcp.tools import _compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],  # All in middle
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],  # High MAPQ
        }
        result = _compute_artifact_risk(variant, evidence, None, 0)
        assert result["artifact_likelihood"] == "low"
        assert len(result["risks"]) == 0

    @pytest.mark.unit
    def test_high_strand_bias(self):
        from bamcp.tools import _compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.98,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = _compute_artifact_risk(variant, evidence, None, 0)
        assert result["artifact_likelihood"] in ["medium", "high"]
        assert any(r["type"] == "strand_bias" for r in result["risks"])

    @pytest.mark.unit
    def test_near_end_bias(self):
        from bamcp.tools import _compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [10, 10, 1, 0, 0, 0],  # Mostly near ends
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = _compute_artifact_risk(variant, evidence, None, 0)
        assert any(r["type"] == "read_position_bias" for r in result["risks"])

    @pytest.mark.unit
    def test_homopolymer_detection(self):
        from bamcp.tools import _compute_artifact_risk

        variant = {"position": 105, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        # Position 5 in reference is in an AAAAA run
        ref_seq = "ACGTAAAAACGT"
        result = _compute_artifact_risk(variant, evidence, ref_seq, 100)
        assert any(r["type"] == "homopolymer" for r in result["risks"])

    @pytest.mark.unit
    def test_low_depth_risk(self):
        from bamcp.tools import _compute_artifact_risk

        variant = {"position": 100, "depth": 3, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = _compute_artifact_risk(variant, evidence, None, 0)
        assert any(r["type"] == "low_depth" for r in result["risks"])


class TestVariantEvidenceHistograms:
    """Tests for histogram data in variant evidence."""

    @pytest.mark.unit
    def test_evidence_includes_histograms(self):
        read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=120,
            reads=[read],
            coverage=[1] * 20,
            variants=[
                {
                    "contig": "chr1",
                    "position": 101,
                    "ref": "C",
                    "alt": "T",
                    "vaf": 0.5,
                    "depth": 20,
                    "alt_count": 1,
                }
            ],
        )

        result = _serialize_region_data(data)
        evidence = result["variant_evidence"]["101:C>T"]

        assert "quality_histogram" in evidence
        assert "position_histogram" in evidence
        assert "mapq_histogram" in evidence
        assert "artifact_risk" in evidence

    @pytest.mark.unit
    def test_evidence_artifact_risk_structure(self):
        read = AlignedRead(
            name="r1",
            sequence="ACGT",
            qualities=[30, 30, 30, 30],
            cigar="4M",
            position=100,
            end_position=104,
            mapping_quality=60,
            is_reverse=False,
            mismatches=[{"pos": 101, "ref": "C", "alt": "T"}],
        )
        data = RegionData(
            contig="chr1",
            start=100,
            end=120,
            reads=[read],
            coverage=[1] * 20,
            variants=[
                {
                    "contig": "chr1",
                    "position": 101,
                    "ref": "C",
                    "alt": "T",
                    "vaf": 0.5,
                    "depth": 20,
                    "alt_count": 1,
                }
            ],
        )

        result = _serialize_region_data(data)
        evidence = result["variant_evidence"]["101:C>T"]
        artifact_risk = evidence["artifact_risk"]

        assert "risks" in artifact_risk
        assert "risk_score" in artifact_risk
        assert "artifact_likelihood" in artifact_risk
        assert artifact_risk["artifact_likelihood"] in ["low", "medium", "high"]


class TestHandleGetVariantCurationSummary:
    """Tests for handle_get_variant_curation_summary tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_formatted_summary(self, small_bam_path, config_with_ref):
        from bamcp.tools import handle_get_variant_curation_summary

        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "chr1",
                "pos": 105,
                "ref": "A",
                "alt": "T",
                "window": 50,
            },
            config_with_ref,
        )
        text = result["content"][0]["text"]
        # Should either find variant or report not found
        assert "Variant" in text or "not found" in text.lower()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_invalid_chromosome_rejected(self, small_bam_path, config):
        from bamcp.tools import handle_get_variant_curation_summary

        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "invalid_chr",
                "pos": 100,
                "ref": "A",
                "alt": "T",
            },
            config,
        )
        text = result["content"][0]["text"]
        assert "Invalid" in text

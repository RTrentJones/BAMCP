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
from bamcp.core.parsers import AlignedRead, RegionData
from bamcp.core.serialization import serialize_region_data
from bamcp.core.tools import (
    close_external_clients,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_scan_variants,
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
    """Tests for serialize_region_data helper."""

    @pytest.mark.unit
    def test_empty_region(self):
        data = RegionData(
            contig="chr1", start=100, end=200, reads=[], coverage=[0] * 100, variants=[]
        )
        result = serialize_region_data(data)
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
        # Default is compact=False, which includes sequence for zoomed views
        result = serialize_region_data(data)
        assert len(result["reads"]) == 1
        r = result["reads"][0]
        assert r["name"] == "r1"
        assert r["sequence"] == "ACGT"  # Default includes sequence
        assert r["cigar"] == "4M"
        assert r["position"] == 100
        assert r["end_position"] == 104
        assert r["mapping_quality"] == 60
        assert r["qualities"] == [30, 30, 30, 30]  # Qualities ride with sequence
        assert r["is_reverse"] is False
        assert len(r["mismatches"]) == 1

    @pytest.mark.unit
    def test_compact_mode_omits_sequence_and_qualities(self):
        """Compact mode omits both sequence and qualities to keep payloads small."""
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
        result = serialize_region_data(data, compact=True)
        r = result["reads"][0]
        assert "sequence" not in r  # Compact omits sequence
        assert "qualities" not in r  # And the matching qualities

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
        result = serialize_region_data(data)
        # Should be JSON serializable
        json_str = json.dumps(result)
        parsed = json.loads(json_str)
        assert parsed == result


class TestHandleBrowseRegion:
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

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_with_reference(self, small_bam_path, ref_fasta_path, config):
        result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200", "reference": ref_fasta_path},
            config,
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["reference_sequence"] is not None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_config_reference_fallback(self, small_bam_path, config_with_ref):
        result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        # Payload is in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["reference_sequence"] is not None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_invalid_region(self, small_bam_path, config):
        with pytest.raises(ValueError):
            await handle_visualize_region(
                {"file_path": small_bam_path, "region": "invalid"}, config
            )


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
        assert payload["coordinate_system"] == "1-based-inclusive"
        assert payload["variant_type"] == "candidate"
        assert "disclaimer" in payload
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

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_positions_are_one_based(self, small_bam_path, config_with_ref):
        """Regression: get_variants must report 1-based positions.

        Internally variants are detected in pysam's 0-based coordinates; every
        genomics consumer a caller reaches next (VCF, dbSNP, ClinVar, gnomAD, IGV)
        is 1-based. Leaking the 0-based value made lookup_clinvar / lookup_gnomad
        miss by one. The LLM-facing position must be the internal one + 1.
        """
        from bamcp.core.tools import _fetch_region_with_timeout

        region = "chr1:490-600"
        data = await _fetch_region_with_timeout(
            small_bam_path,
            region,
            config_with_ref.reference,
            config_with_ref,
            min_vaf=0.05,
            min_depth=1,
        )
        raw_positions = sorted(v["position"] for v in data.variants)  # 0-based
        assert raw_positions, "fixture region should contain at least one variant"

        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": region, "min_vaf": 0.05, "min_depth": 1},
            config_with_ref,
        )
        tool_positions = sorted(
            v["position"] for v in json.loads(result["content"][0]["text"])["variants"]
        )
        assert tool_positions == [p + 1 for p in raw_positions]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_position_round_trips_into_curation(self, small_bam_path, config_with_ref):
        """A position from get_variants must be found by get_variant_curation_summary.

        Guards the convention pairing: both speak 1-based, so the exact
        detect -> curate/lookup chain a caller follows resolves the variant.
        """
        from bamcp.analysis.curation import handle_get_variant_curation_summary

        result = await handle_get_variants(
            {
                "file_path": small_bam_path,
                "region": "chr1:490-600",
                "min_vaf": 0.05,
                "min_depth": 1,
            },
            config_with_ref,
        )
        variants = json.loads(result["content"][0]["text"])["variants"]
        assert variants, "need a variant to curate"
        v = variants[0]

        summary = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": v["contig"],
                "pos": v["position"],
                "ref": v["ref"],
                "alt": v["alt"],
                "window": 50,
            },
            config_with_ref,
        )
        text = summary["content"][0]["text"]
        assert "not found" not in text.lower(), text[:200]


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
        assert "Candidate variants detected:" in text

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
        from bamcp.clients.clinvar import ClinVarResult

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

        from bamcp.clients import clinvar

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
        from bamcp.clients.clinvar import ClinVarResult

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

        from bamcp.clients import clinvar

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

        from bamcp.clients import clinvar

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

        from bamcp.clients import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", mock_lookup)

        result = await handle_lookup_clinvar(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert "error" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_disabled_by_config(self, monkeypatch):
        async def fail_lookup(self, chrom, pos, ref, alt):
            raise AssertionError("disabled ClinVar should not call API client")

        from bamcp.clients import clinvar

        monkeypatch.setattr(clinvar.ClinVarClient, "lookup", fail_lookup)
        config = BAMCPConfig(clinvar_enabled=False)

        result = await handle_lookup_clinvar(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["error"] == "ClinVar lookup is disabled"
        assert "disclaimer" in payload


class TestHandleLookupGnomad:
    """Tests for handle_lookup_gnomad tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_frequency(self, config, monkeypatch):
        from bamcp.clients.gnomad import GnomadResult, PopulationFrequency

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

        from bamcp.clients import gnomad

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

        from bamcp.clients import gnomad

        monkeypatch.setattr(gnomad.GnomadClient, "lookup", mock_lookup)

        result = await handle_lookup_gnomad(
            {"chrom": "chr99", "pos": 1, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["found"] is False
        assert "disclaimer" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_disabled_by_config(self, monkeypatch):
        async def fail_lookup(self, chrom, pos, ref, alt):
            raise AssertionError("disabled gnomAD should not call API client")

        from bamcp.clients import gnomad

        monkeypatch.setattr(gnomad.GnomadClient, "lookup", fail_lookup)
        config = BAMCPConfig(gnomad_enabled=False)

        result = await handle_lookup_gnomad(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["error"] == "gnomAD lookup is disabled"
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

        result = serialize_region_data(data)

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

        result = serialize_region_data(data)

        assert (
            result["variant_evidence"]["101:C>T"]["strand_bias"] <= LOW_CONFIDENCE_MAX_STRAND_BIAS
        )
        assert result["variants"][0]["is_low_confidence"] is False


class TestBinValues:
    """Tests for bin_values helper function."""

    @pytest.mark.unit
    def test_empty_values(self):
        from bamcp.analysis.evidence import bin_values

        result = bin_values([], [0, 10, 20, 30])
        assert result == [0, 0, 0, 0]

    @pytest.mark.unit
    def test_single_bin(self):
        from bamcp.analysis.evidence import bin_values

        result = bin_values([5, 15, 25], [0, 10, 20, 30])
        assert result == [1, 1, 1, 0]

    @pytest.mark.unit
    def test_all_in_last_bin(self):
        from bamcp.analysis.evidence import bin_values

        result = bin_values([35, 40, 100], [0, 10, 20, 30])
        assert result == [0, 0, 0, 3]

    @pytest.mark.unit
    def test_quality_histogram_bins(self):
        from bamcp.analysis.evidence import bin_values
        from bamcp.constants import QUALITY_HISTOGRAM_BINS

        values = [5, 15, 25, 35, 45, 50]
        result = bin_values(values, QUALITY_HISTOGRAM_BINS)
        assert len(result) == 5
        assert result == [1, 1, 1, 1, 2]  # 45 and 50 go in last bin


class TestDetectHomopolymer:
    """Tests for detect_homopolymer helper function."""

    @pytest.mark.unit
    def test_single_base(self):
        from bamcp.analysis.evidence import detect_homopolymer

        result = detect_homopolymer("ACGT", 0)
        assert result == 1

    @pytest.mark.unit
    def test_homopolymer_run(self):
        from bamcp.analysis.evidence import detect_homopolymer

        result = detect_homopolymer("AAAAAACGT", 3)
        assert result == 6

    @pytest.mark.unit
    def test_position_at_end(self):
        from bamcp.analysis.evidence import detect_homopolymer

        result = detect_homopolymer("CGTTTTT", 5)
        assert result == 5

    @pytest.mark.unit
    def test_empty_sequence(self):
        from bamcp.analysis.evidence import detect_homopolymer

        result = detect_homopolymer("", 0)
        assert result == 0

    @pytest.mark.unit
    def test_position_out_of_range(self):
        from bamcp.analysis.evidence import detect_homopolymer

        result = detect_homopolymer("ACGT", 10)
        assert result == 0


class TestComputeArtifactRisk:
    """Tests for compute_artifact_risk function."""

    @pytest.mark.unit
    def test_low_risk(self):
        from bamcp.analysis.evidence import compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],  # All in middle
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],  # High MAPQ
        }
        result = compute_artifact_risk(variant, evidence, None, 0)
        assert result["artifact_likelihood"] == "low"
        assert len(result["risks"]) == 0

    @pytest.mark.unit
    def test_high_strand_bias(self):
        from bamcp.analysis.evidence import compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.98,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = compute_artifact_risk(variant, evidence, None, 0)
        assert result["artifact_likelihood"] in ["medium", "high"]
        assert any(r["type"] == "strand_bias" for r in result["risks"])

    @pytest.mark.unit
    def test_near_end_bias(self):
        from bamcp.analysis.evidence import compute_artifact_risk

        variant = {"position": 100, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [10, 10, 1, 0, 0, 0],  # Mostly near ends
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = compute_artifact_risk(variant, evidence, None, 0)
        assert any(r["type"] == "read_position_bias" for r in result["risks"])

    @pytest.mark.unit
    def test_homopolymer_detection(self):
        from bamcp.analysis.evidence import compute_artifact_risk

        variant = {"position": 105, "depth": 50, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        # Position 5 in reference is in an AAAAA run
        ref_seq = "ACGTAAAAACGT"
        result = compute_artifact_risk(variant, evidence, ref_seq, 100)
        assert any(r["type"] == "homopolymer" for r in result["risks"])

    @pytest.mark.unit
    def test_low_depth_risk(self):
        from bamcp.analysis.evidence import compute_artifact_risk

        variant = {"position": 100, "depth": 3, "vaf": 0.5}
        evidence = {
            "strand_bias": 0.1,
            "position_histogram": [0, 0, 5, 5, 5, 5],
            "mapq_histogram": [0, 0, 0, 0, 5, 10, 20],
        }
        result = compute_artifact_risk(variant, evidence, None, 0)
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

        result = serialize_region_data(data)
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

        result = serialize_region_data(data)
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
        from bamcp.analysis.curation import handle_get_variant_curation_summary

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
        from bamcp.analysis.curation import handle_get_variant_curation_summary

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


class TestHandleScanVariants:
    """Tests for handle_scan_variants tool handler."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_requires_reference(self, small_bam_path, config):
        """Should return error when no reference is configured."""
        result = await handle_scan_variants({"file_path": small_bam_path, "contig": "chr1"}, config)
        payload = json.loads(result["content"][0]["text"])
        assert "error" in payload
        assert "reference" in payload["error"].lower()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_returns_variants(self, small_bam_path, config_with_ref):
        """Should return ranked variants for chr1."""
        result = await handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, config_with_ref
        )
        payload = json.loads(result["content"][0]["text"])
        assert "variants" in payload
        assert "contig" in payload
        assert payload["coordinate_system"] == "1-based-inclusive"
        assert payload["variant_type"] == "candidate"
        assert "disclaimer" in payload
        assert payload["contig"] == "chr1"
        assert isinstance(payload["variants"], list)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_variants_sorted_by_vaf(self, small_bam_path, config_with_ref):
        """Returned variants should be sorted by VAF descending."""
        result = await handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chr1"}, config_with_ref
        )
        payload = json.loads(result["content"][0]["text"])
        variants = payload["variants"]
        if len(variants) >= 2:
            for i in range(len(variants) - 1):
                assert variants[i]["vaf"] >= variants[i + 1]["vaf"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_nonexistent_contig_empty(self, small_bam_path, config_with_ref):
        """Should return empty list for nonexistent contig."""
        result = await handle_scan_variants(
            {"file_path": small_bam_path, "contig": "chrZ"}, config_with_ref
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["variants"] == []
        assert payload["count"] == 0


class TestMediumRoadmapCoverage:
    """Targeted coverage for medium-roadmap fast paths, caching, and cleanup."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_region_coverage_mode_uses_coverage_fast_path(self, config, monkeypatch):
        """Coverage mode should avoid the full read-serialization parser."""
        from bamcp.core import tools as tools_module

        calls = []

        async def mock_ensure_cached_index(file_path, cfg):
            return None

        def mock_fetch_coverage_only(*args):
            calls.append(args)
            return RegionData("chr1", 10, 20, reads=[], coverage=[3] * 10, variants=[])

        def fail_fetch_region(*args):
            raise AssertionError("coverage mode should not call fetch_region")

        monkeypatch.setattr(tools_module, "_ensure_cached_index", mock_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_coverage_only", mock_fetch_coverage_only)
        monkeypatch.setattr(tools_module, "fetch_region", fail_fetch_region)

        data = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, config, mode="coverage"
        )

        assert data.coverage == [3] * 10
        assert len(calls) == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_region_cache_reuses_identical_query(self, config, monkeypatch):
        """Identical region queries should hit the in-process region cache."""
        from bamcp.core import tools as tools_module

        call_count = 0

        async def mock_ensure_cached_index(file_path, cfg):
            return None

        def mock_fetch_coverage_only(*args):
            nonlocal call_count
            call_count += 1
            return RegionData("chr1", 10, 20, reads=[], coverage=[call_count] * 10, variants=[])

        monkeypatch.setattr(tools_module, "_ensure_cached_index", mock_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_coverage_only", mock_fetch_coverage_only)

        cfg = BAMCPConfig(cache_ttl=60)
        first = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, cfg, mode="coverage"
        )
        second = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, cfg, mode="coverage"
        )

        # Readless queries are tile-cached and sliced per call, so the objects
        # differ but the underlying tile is computed only once.
        assert first.coverage == second.coverage
        assert call_count == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_region_cache_separates_query_mode(self, config, monkeypatch):
        """Coverage and full modes must not share cached RegionData."""
        from bamcp.core import tools as tools_module

        async def mock_ensure_cached_index(file_path, cfg):
            return None

        def mock_fetch_coverage_only(*args):
            return RegionData("chr1", 10, 20, reads=[], coverage=[1] * 10, variants=[])

        def mock_fetch_region(*args):
            return RegionData(
                "chr1",
                10,
                20,
                reads=[
                    AlignedRead(
                        name="r1",
                        sequence="A",
                        qualities=[30],
                        cigar="1M",
                        position=10,
                        end_position=11,
                        mapping_quality=60,
                        is_reverse=False,
                    )
                ],
                coverage=[2] * 10,
                variants=[],
            )

        monkeypatch.setattr(tools_module, "_ensure_cached_index", mock_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_coverage_only", mock_fetch_coverage_only)
        monkeypatch.setattr(tools_module, "fetch_region", mock_fetch_region)

        cfg = BAMCPConfig(cache_ttl=60)
        coverage = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, cfg, mode="coverage"
        )
        full = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, cfg, mode="full"
        )

        assert coverage.reads == []
        assert len(full.reads) == 1

    @pytest.mark.unit
    def test_tile_bounds_snap_outward(self):
        from bamcp.core.tools import _tile_bounds

        assert _tile_bounds(100, 200, 4096) == (0, 4096)
        assert _tile_bounds(0, 4096, 4096) == (0, 4096)
        assert _tile_bounds(4096, 4100, 4096) == (4096, 8192)
        assert _tile_bounds(4095, 4097, 4096) == (0, 8192)

    @pytest.mark.unit
    def test_slice_readless_region_data(self):
        from bamcp.core.tools import _slice_readless_region_data

        data = RegionData(
            "chr1",
            0,
            10,
            reads=[],
            coverage=list(range(10)),
            variants=[
                {"position": 2, "ref": "A", "alt": "T"},
                {"position": 7, "ref": "C", "alt": "G"},
                {"position": 9, "ref": "A", "alt": "C"},
            ],
            reference_sequence="ACGTACGTAC",
        )
        sliced = _slice_readless_region_data(data, 2, 8)
        assert (sliced.start, sliced.end) == (2, 8)
        assert sliced.coverage == [2, 3, 4, 5, 6, 7]
        assert sliced.reference_sequence == "GTACGT"  # ref[2:8]
        assert [v["position"] for v in sliced.variants] == [2, 7]  # 2 <= pos < 8
        assert sliced.reads == []

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_readless_tile_reused_across_overlapping_regions(self, monkeypatch):
        """Two sub-regions in the same tile compute once and slice correctly."""
        from bamcp.core import tools as tools_module

        call_count = 0

        async def mock_ensure_cached_index(file_path, cfg):
            return None

        def mock_fetch_coverage_only(*args):
            nonlocal call_count
            call_count += 1
            # Tile region is chr1:0-4096 -> coverage indexed by absolute position.
            return RegionData("chr1", 0, 4096, reads=[], coverage=list(range(4096)), variants=[])

        monkeypatch.setattr(tools_module, "_ensure_cached_index", mock_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_coverage_only", mock_fetch_coverage_only)

        cfg = BAMCPConfig(cache_ttl=60)
        a = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:100-200", None, cfg, mode="coverage"
        )
        b = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:300-400", None, cfg, mode="coverage"
        )

        assert call_count == 1  # one tile served both sub-regions
        assert (a.start, a.end) == (100, 200)
        assert a.coverage == list(range(100, 200))
        assert b.coverage == list(range(300, 400))

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_tiled_coverage_matches_direct_fetch(self, small_bam_path, config):
        """Tiled+sliced coverage must equal a direct exact fetch of the same region."""
        from bamcp.core.parsers import fetch_coverage_only
        from bamcp.core.tools import _fetch_region_with_timeout

        region = "chr1:100-160"
        tiled = await _fetch_region_with_timeout(
            small_bam_path, region, None, config, mode="coverage"
        )
        direct = fetch_coverage_only(
            small_bam_path, region, None, config.min_mapq, config.min_baseq, None
        )
        assert (tiled.start, tiled.end) == (100, 160)
        assert tiled.coverage == direct.coverage

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_region_overlays_vcf_and_deduplicates(self, config, monkeypatch, tmp_path):
        """VCF variants should be overlaid without duplicating matching candidate variants."""
        from bamcp.core import tools as tools_module

        vcf_path = tmp_path / "calls.vcf.gz"
        vcf_path.touch()

        async def mock_ensure_cached_index(file_path, cfg):
            return None

        def mock_fetch_region(*args):
            return RegionData(
                "chr1",
                10,
                20,
                reads=[],
                coverage=[1] * 10,
                variants=[
                    {
                        "contig": "chr1",
                        "position": 12,
                        "ref": "A",
                        "alt": "T",
                        "vaf": 0.5,
                        "depth": 10,
                        "alt_count": 5,
                    }
                ],
            )

        def mock_load_vcf_variants(path, region):
            assert path == str(vcf_path)
            assert region == "chr1:10-20"
            return [
                {
                    "contig": "chr1",
                    "position": 12,
                    "ref": "A",
                    "alt": "T",
                    "vaf": 0.5,
                    "depth": 10,
                    "alt_count": 5,
                    "source": "vcf",
                },
                {
                    "contig": "chr1",
                    "position": 15,
                    "ref": "G",
                    "alt": "C",
                    "vaf": 0.0,
                    "depth": 0,
                    "alt_count": 0,
                    "source": "vcf",
                },
            ]

        monkeypatch.setattr(tools_module, "_ensure_cached_index", mock_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_region", mock_fetch_region)
        monkeypatch.setattr(tools_module, "load_vcf_variants", mock_load_vcf_variants)

        data = await tools_module._fetch_region_with_timeout(
            "sample.bam", "chr1:10-20", None, config, vcf_path=str(vcf_path)
        )

        assert [(v["position"], v["ref"], v["alt"]) for v in data.variants] == [
            (12, "A", "T"),
            (15, "G", "C"),
        ]
        assert data.variants[1]["source"] == "vcf"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_region_uses_vcf_as_primary_variant_source_without_reference(
        self, config, monkeypatch, tmp_path
    ):
        """Variant mode with VCF and no reference should avoid BAM parsing entirely."""
        from bamcp.core import tools as tools_module

        vcf_path = tmp_path / "primary.vcf"
        vcf_path.touch()

        def fail_ensure_cached_index(file_path, cfg):
            raise AssertionError("VCF-primary variant mode should not fetch BAM indexes")

        def fail_fetch_region(*args):
            raise AssertionError("VCF-primary variant mode should not parse BAM reads")

        def mock_load_vcf_variants(path, region):
            assert path == str(vcf_path)
            assert region == "chr3:100-200"
            return [
                {
                    "contig": "chr3",
                    "position": 150,
                    "ref": "N",
                    "alt": "<DUP>",
                    "variant_kind": "sv",
                    "vaf": 0.0,
                    "depth": 0,
                    "alt_count": 0,
                    "source": "vcf",
                    "sample_names": ["sample1"],
                    "samples": {"sample1": {"GT": [0, 1]}},
                    "sv_type": "DUP",
                    "sv_end": 300,
                    "sv_len": 150,
                }
            ]

        monkeypatch.setattr(tools_module, "_ensure_cached_index", fail_ensure_cached_index)
        monkeypatch.setattr(tools_module, "fetch_region", fail_fetch_region)
        monkeypatch.setattr(tools_module, "load_vcf_variants", mock_load_vcf_variants)

        data = await tools_module._fetch_region_with_timeout(
            "sample.bam",
            "chr3:100-200",
            None,
            config,
            mode="variants",
            vcf_path=str(vcf_path),
        )

        assert data.reads == []
        assert data.coverage == []
        assert data.variants[0]["variant_kind"] == "sv"
        assert data.variants[0]["sample_names"] == ["sample1"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_download_index_streaming_revalidates_and_size_caps(
        self, config, monkeypatch, tmp_path
    ):
        """Streaming index download should validate derived URLs and enforce size caps."""
        from bamcp.core import tools as tools_module

        validated = []

        def mock_validate_remote_url(url, cfg):
            validated.append(url)

        class ResponseStub:
            status_code = 200
            headers = {"content-length": "5"}

            async def aiter_bytes(self):
                yield b"abc"
                yield b"de"

        class StreamContext:
            async def __aenter__(self):
                return ResponseStub()

            async def __aexit__(self, exc_type, exc, tb):
                return False

        class ClientStub:
            def stream(self, method, url):
                assert method == "GET"
                assert url == "https://example.com/sample.bam.bai"
                return StreamContext()

        monkeypatch.setattr(tools_module, "validate_remote_url", mock_validate_remote_url)
        index_path = tmp_path / "sample.bai"

        result = await tools_module._download_index_streaming(
            ClientStub(),
            "https://example.com/sample.bam.bai",
            str(index_path),
            config,
        )

        assert result == str(index_path)
        assert index_path.read_bytes() == b"abcde"
        assert validated == ["https://example.com/sample.bam.bai"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_download_index_streaming_rejects_oversized_content_length(
        self, monkeypatch, tmp_path
    ):
        """Content-Length above the configured cap should abort before writing."""
        from bamcp.core import tools as tools_module

        class ResponseStub:
            status_code = 200
            headers = {"content-length": "6"}

            async def aiter_bytes(self):
                yield b"abcdef"

        class StreamContext:
            async def __aenter__(self):
                return ResponseStub()

            async def __aexit__(self, exc_type, exc, tb):
                return False

        class ClientStub:
            def stream(self, method, url):
                return StreamContext()

        monkeypatch.setattr(tools_module, "validate_remote_url", lambda url, cfg: None)

        with pytest.raises(ValueError, match="exceeds configured maximum"):
            await tools_module._download_index_streaming(
                ClientStub(),
                "https://example.com/sample.bam.bai",
                str(tmp_path / "sample.bai"),
                BAMCPConfig(max_remote_index_bytes=5),
            )

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_close_external_clients_closes_clients_and_clears_cache(self, config):
        """External client cleanup should close clients and clear the region cache."""
        from bamcp.core.tools import get_services

        closed = []

        class ClientStub:
            def __init__(self, name):
                self.name = name

            async def close(self):
                closed.append(self.name)

        services = get_services(config)
        services._clinvar = ClientStub("clinvar")
        services._gnomad = ClientStub("gnomad")
        services._genes = ClientStub("gene")
        services.region_cache[("coverage", "sample.bam")] = (
            0,
            RegionData("chr1", 0, 1, [], [0], []),
        )

        await close_external_clients(config)

        assert sorted(closed) == ["clinvar", "gene", "gnomad"]
        assert services._clinvar is None
        assert services._gnomad is None
        assert services._genes is None
        assert services.region_cache == {}

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_services_isolated_per_config(self):
        """Two servers built from different configs must not share cache or clients."""
        from bamcp.core.tools import get_services

        cfg_a = BAMCPConfig(gnomad_dataset="gnomad_r4")
        cfg_b = BAMCPConfig(gnomad_dataset="gnomad_r2_1")

        services_a = get_services(cfg_a)
        services_b = get_services(cfg_b)

        assert services_a is not services_b
        assert services_a.cache is not services_b.cache
        # gnomAD clients are built lazily from each config's dataset.
        assert services_a.gnomad() is not services_b.gnomad()
        assert services_a.gnomad().dataset == "gnomad_r4"
        assert services_b.gnomad().dataset == "gnomad_r2_1"

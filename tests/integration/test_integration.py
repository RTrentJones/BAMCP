"""Integration tests for BAMCP - testing full tool call workflows."""

import json

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core.tools import (
    handle_get_coverage,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_visualize_region,
)


@pytest.fixture
def config_with_ref(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path, min_depth=1, min_vaf=0.05)


class TestBrowseRegionIntegration:
    """Integration tests for browse_region end-to-end flow."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_full_browse_region_flow(self, small_bam_path, config_with_ref):
        """Test complete browse_region: fetch -> serialize -> verify all data."""
        result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )

        # Verify structure
        assert result["_meta"]["ui/resourceUri"] == "ui://bamcp/viewer"

        # Content should be a summary string (not JSON payload)
        content_text = result["content"][0]["text"]
        assert "Region chr1:90-200" in content_text
        assert "reads" in content_text
        assert "variants" in content_text

        # Payload is now only in _meta.ui/init (not duplicated in content)
        payload = result["_meta"]["ui/init"]

        # Verify region
        assert payload["contig"] == "chr1"
        assert payload["start"] == 90
        assert payload["end"] == 200

        # Verify reads are present
        assert len(payload["reads"]) > 0

        # Verify read structure (compact mode is now default)
        read = payload["reads"][0]
        required_fields = {
            "name",
            "cigar",
            "position",
            "end_position",
            "mapping_quality",
            "is_reverse",
            "mismatches",
        }
        # In compact mode, sequence and qualities are omitted
        assert required_fields.issubset(set(read.keys()))
        assert "qualities" not in read  # Qualities are never serialized (unused by frontend)

        # Verify coverage length
        assert len(payload["coverage"]) == 110

        # Verify reference
        assert payload["reference_sequence"] is not None
        assert len(payload["reference_sequence"]) == 110

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_browse_then_variant_consistency(self, small_bam_path, config_with_ref):
        """browse_region and get_variants should return consistent variant data."""
        browse_result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:490-600"}, config_with_ref
        )
        variant_result = await handle_get_variants(
            {
                "file_path": small_bam_path,
                "region": "chr1:490-600",
                "min_vaf": 0.05,
                "min_depth": 1,
            },
            config_with_ref,
        )

        # browse_region payload is now in _meta.ui/init, not content
        browse_payload = browse_result["_meta"]["ui/init"]
        variant_payload = json.loads(variant_result["content"][0]["text"])

        # Variant positions from browse should be a superset (get_variants may filter)
        browse_positions = {v["position"] for v in browse_payload["variants"]}
        variant_positions = {v["position"] for v in variant_payload["variants"]}
        assert variant_positions.issubset(browse_positions)


class TestCoverageIntegration:
    """Integration tests for coverage calculation."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_coverage_matches_reads(self, small_bam_path, ref_fasta_path):
        """Coverage should correlate with reads in the region."""
        config = BAMCPConfig(reference=ref_fasta_path)

        browse_result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )
        coverage_result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )

        # browse_region payload is now in _meta.ui/init
        browse_payload = browse_result["_meta"]["ui/init"]
        cov_payload = json.loads(coverage_result["content"][0]["text"])

        # If there are reads, coverage should be non-zero
        if len(browse_payload["reads"]) > 0:
            assert cov_payload["max"] > 0
            assert cov_payload["bases_covered"] > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_coverage_stats_consistency(self, small_bam_path):
        """Coverage stats should be internally consistent."""
        config = BAMCPConfig()

        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )
        stats = json.loads(result["content"][0]["text"])

        assert stats["min"] <= stats["mean"] <= stats["max"]
        assert stats["bases_covered"] <= stats["total_bases"]
        assert stats["total_bases"] == 100


class TestContigListIntegration:
    """Integration tests for contig listing."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_contigs_usable_for_browse(self, small_bam_path):
        """Contigs from list_contigs should be usable in browse_region."""
        config = BAMCPConfig()

        contigs_result = await handle_list_contigs({"file_path": small_bam_path}, config)
        contigs = json.loads(contigs_result["content"][0]["text"])["contigs"]

        # Browse each contig
        for contig in contigs:
            name = contig["name"]
            length = contig["length"]
            region = f"{name}:0-{min(length, 100)}"

            result = await handle_visualize_region(
                {"file_path": small_bam_path, "region": region}, config
            )
            # browse_region payload is now in _meta.ui/init
            payload = result["_meta"]["ui/init"]
            assert payload["contig"] == name


class TestMultiToolWorkflow:
    """Integration tests for multi-tool workflows."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_list_browse_coverage_workflow(self, small_bam_path, config_with_ref):
        """Simulate typical user workflow: list -> browse -> coverage."""
        # Step 1: List contigs
        contigs = await handle_list_contigs({"file_path": small_bam_path}, config_with_ref)
        contig_data = json.loads(contigs["content"][0]["text"])
        next(c for c in contig_data["contigs"] if c["name"] == "chr1")

        # Step 2: Browse a region
        browse = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        # browse_region payload is now in _meta.ui/init
        browse_data = browse["_meta"]["ui/init"]
        assert len(browse_data["reads"]) > 0

        # Step 3: Get coverage for same region
        coverage = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        cov_data = json.loads(coverage["content"][0]["text"])
        assert cov_data["max"] > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_different_regions_same_file(self, small_bam_path):
        """Should handle multiple regions from the same file."""
        config = BAMCPConfig()

        regions = ["chr1:90-200", "chr1:290-380", "chr1:490-600", "chr2:40-100"]
        for region in regions:
            result = await handle_visualize_region(
                {"file_path": small_bam_path, "region": region}, config
            )
            # browse_region payload is now in _meta.ui/init
            payload = result["_meta"]["ui/init"]
            assert "reads" in payload
            assert "coverage" in payload


class TestJumpToIntegration:
    """Integration tests for jump_to tool."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_jump_to_returns_centered_data(self, small_bam_path):
        """jump_to should return data centered on the given position."""
        config = BAMCPConfig()

        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 200},
            config,
        )
        # jump_to payload is now in _meta.ui/init
        payload = result["_meta"]["ui/init"]
        assert payload["contig"] == "chr1"
        assert payload["start"] <= 150 <= payload["end"]
        assert payload["end"] - payload["start"] == 200

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_jump_then_browse_same_region(self, small_bam_path):
        """jump_to and browse_region for the same region should return the same reads."""
        config = BAMCPConfig()

        jump_result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 200},
            config,
        )
        # jump_to payload is now in _meta.ui/init
        jump_payload = jump_result["_meta"]["ui/init"]

        region = f"chr1:{jump_payload['start']}-{jump_payload['end']}"
        browse_result = await handle_visualize_region(
            {"file_path": small_bam_path, "region": region}, config
        )
        # browse_region payload is now in _meta.ui/init
        browse_payload = browse_result["_meta"]["ui/init"]

        # Same region should yield same number of reads
        assert len(jump_payload["reads"]) == len(browse_payload["reads"])

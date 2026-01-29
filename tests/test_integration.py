"""Integration tests for BAMCP - testing full tool call workflows."""

import json
import pytest

from bamcp.config import BAMCPConfig
from bamcp.server import create_server
from bamcp.tools import (
    handle_browse_region,
    handle_get_coverage,
    handle_get_variants,
    handle_list_contigs,
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
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )

        # Verify structure
        assert result["_meta"]["ui/resourceUri"] == "ui://bamcp/viewer"

        # Parse payload
        payload = json.loads(result["content"][0]["text"])

        # Verify region
        assert payload["contig"] == "chr1"
        assert payload["start"] == 90
        assert payload["end"] == 200

        # Verify reads are present
        assert len(payload["reads"]) > 0

        # Verify read structure
        read = payload["reads"][0]
        required_fields = {
            "name", "sequence", "cigar", "position",
            "end_position", "mapping_quality", "is_reverse", "mismatches",
        }
        assert set(read.keys()) == required_fields

        # Verify coverage length
        assert len(payload["coverage"]) == 110

        # Verify reference
        assert payload["reference_sequence"] is not None
        assert len(payload["reference_sequence"]) == 110

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_browse_then_variant_consistency(self, small_bam_path, config_with_ref):
        """browse_region and get_variants should return consistent variant data."""
        browse_result = await handle_browse_region(
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

        browse_payload = json.loads(browse_result["content"][0]["text"])
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

        browse_result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )
        coverage_result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:100-200"}, config
        )

        browse_payload = json.loads(browse_result["content"][0]["text"])
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

        contigs_result = await handle_list_contigs(
            {"file_path": small_bam_path}, config
        )
        contigs = json.loads(contigs_result["content"][0]["text"])["contigs"]

        # Browse each contig
        for contig in contigs:
            name = contig["name"]
            length = contig["length"]
            region = f"{name}:0-{min(length, 100)}"

            result = await handle_browse_region(
                {"file_path": small_bam_path, "region": region}, config
            )
            payload = json.loads(result["content"][0]["text"])
            assert payload["contig"] == name


class TestMultiToolWorkflow:
    """Integration tests for multi-tool workflows."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_list_browse_coverage_workflow(self, small_bam_path, config_with_ref):
        """Simulate typical user workflow: list -> browse -> coverage."""
        # Step 1: List contigs
        contigs = await handle_list_contigs(
            {"file_path": small_bam_path}, config_with_ref
        )
        contig_data = json.loads(contigs["content"][0]["text"])
        chr1 = next(c for c in contig_data["contigs"] if c["name"] == "chr1")

        # Step 2: Browse a region
        browse = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        browse_data = json.loads(browse["content"][0]["text"])
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
            result = await handle_browse_region(
                {"file_path": small_bam_path, "region": region}, config
            )
            payload = json.loads(result["content"][0]["text"])
            assert "reads" in payload
            assert "coverage" in payload

"""Unit tests for bamcp.tools module."""

import json

import pytest

from bamcp.config import BAMCPConfig
from bamcp.parsers import AlignedRead, RegionData
from bamcp.tools import (
    _serialize_region_data,
    handle_browse_region,
    handle_get_coverage,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
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
        result = _serialize_region_data(data)
        assert len(result["reads"]) == 1
        r = result["reads"][0]
        assert r["name"] == "r1"
        assert r["sequence"] == "ACGT"
        assert r["cigar"] == "4M"
        assert r["position"] == 100
        assert r["end_position"] == 104
        assert r["mapping_quality"] == 60
        assert r["qualities"] == [30, 30, 30, 30]
        assert r["is_reverse"] is False
        assert len(r["mismatches"]) == 1

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
        assert result["_meta"]["ui/resourceUri"] == "ui://bamcp/viewer"
        assert "ui/init" in result["_meta"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_content_is_json(self, small_bam_path, config):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config
        )
        text = result["content"][0]["text"]
        payload = json.loads(text)
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
        payload = json.loads(result["content"][0]["text"])
        assert payload["reference_sequence"] is not None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_config_reference_fallback(self, small_bam_path, config_with_ref):
        result = await handle_browse_region(
            {"file_path": small_bam_path, "region": "chr1:90-200"}, config_with_ref
        )
        payload = json.loads(result["content"][0]["text"])
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
        assert result["_meta"]["ui/resourceUri"] == "ui://bamcp/viewer"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_centers_on_position(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["contig"] == "chr1"
        # Position 150 should be within the returned start-end range
        assert payload["start"] <= 150 <= payload["end"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_custom_window(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 100},
            config,
        )
        payload = json.loads(result["content"][0]["text"])
        span = payload["end"] - payload["start"]
        assert span == 100

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_window_from_config(self, small_bam_path):
        config = BAMCPConfig(default_window=200)
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 300, "contig": "chr1"}, config
        )
        payload = json.loads(result["content"][0]["text"])
        span = payload["end"] - payload["start"]
        assert span == 200

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_contig(self, small_bam_path, config):
        """Without contig arg, defaults to chr1."""
        result = await handle_jump_to({"file_path": small_bam_path, "position": 150}, config)
        payload = json.loads(result["content"][0]["text"])
        assert payload["contig"] == "chr1"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_reads_returned(self, small_bam_path, config):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 150, "contig": "chr1", "window": 200},
            config,
        )
        payload = json.loads(result["content"][0]["text"])
        assert "reads" in payload
        assert len(payload["reads"]) > 0

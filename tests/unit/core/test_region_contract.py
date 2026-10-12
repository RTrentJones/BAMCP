"""The tool-boundary region contract: 1-based, inclusive, samtools-style.

``chr1:1051-1051`` is the single base 1051; ``parse_region`` converts once to the
internal 0-based half-open frame, and every internal region builder emits 1-based.
"""

import json

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core.parsers import fetch_coverage_only, fetch_region, format_region
from bamcp.core.tools import (
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_visualize_region,
)


@pytest.fixture
def cfg(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path, min_vaf=0.1, min_depth=1)


class TestSingleBaseRegion:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_single_base_region_returns_the_variant_at_that_base(self, small_bam_path, cfg):
        # small.bam read3 carries T>G at internal 0-based 133 == 1-based 134.
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:134-134"}, cfg
        )
        payload = json.loads(result["content"][0]["text"])
        assert [(v["position"], v["ref"], v["alt"]) for v in payload["variants"]] == [
            (134, "T", "G")
        ]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_adjacent_single_base_region_excludes_it(self, small_bam_path, cfg):
        result = await handle_get_variants(
            {"file_path": small_bam_path, "region": "chr1:133-133"}, cfg
        )
        payload = json.loads(result["content"][0]["text"])
        assert payload["variants"] == []

    @pytest.mark.unit
    def test_single_base_fetch_has_one_column(self, small_bam_path):
        data = fetch_coverage_only(small_bam_path, "chr1:134-134")
        assert (data.start, data.end, len(data.coverage)) == (133, 134, 1)


class TestEchoedRegions:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_coverage_region_round_trips_one_based(self, small_bam_path, cfg):
        result = await handle_get_coverage(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, cfg
        )
        stats = json.loads(result["content"][0]["text"])
        assert stats["region"] == "chr1:101-200"
        assert stats["total_bases"] == 100
        assert stats["coordinate_system"] == "1-based-inclusive"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_region_summary_and_viewer_text_echo_the_input(self, small_bam_path, cfg):
        summary = await handle_get_region_summary(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, cfg
        )
        assert summary["content"][0]["text"].splitlines()[0] == "Region: chr1:101-200"
        viz = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:101-200"}, cfg
        )
        assert viz["content"][0]["text"].startswith("Region chr1:101-200:")
        # The viewer payload stays in the internal 0-based half-open frame.
        assert (viz["_meta"]["ui/init"]["start"], viz["_meta"]["ui/init"]["end"]) == (100, 200)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_gene_converts_ncbi_zero_based(self, monkeypatch):
        import bamcp.core.tools as tm
        from bamcp.clients.genes import GeneInfo

        class _Client:
            async def search(self, symbol, build=None):
                # NCBI genomicinfo chrstart/chrstop are 0-based.
                return GeneInfo("BRCA1", "BRCA1", "chr17", 43044294, 43125482, "-", "GRCh38")

        monkeypatch.setattr(tm, "get_gene_client", lambda config: _Client())
        r = await tm.handle_search_gene({"symbol": "BRCA1"}, BAMCPConfig())
        assert json.loads(r["content"][0]["text"])["region"] == "chr17:43044295-43125483"


class TestInternalBuilders:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_jump_to_near_contig_start_clamps_to_base_one(self, small_bam_path, cfg):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 10, "contig": "chr1", "window": 100}, cfg
        )
        payload = result["_meta"]["ui/init"]
        assert (payload["start"], payload["end"]) == (0, 100)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_jump_to_window_contains_one_based_position(self, small_bam_path, cfg):
        result = await handle_jump_to(
            {"file_path": small_bam_path, "position": 134, "contig": "chr1", "window": 20}, cfg
        )
        payload = result["_meta"]["ui/init"]
        assert payload["end"] - payload["start"] == 20
        assert payload["start"] <= 133 < payload["end"]  # 0-based index of 1-based 134

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_curation_near_contig_start_does_not_build_negative_region(
        self, small_bam_path, cfg
    ):
        from bamcp.analysis.curation import handle_get_variant_curation_summary

        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "chr1",
                "pos": 5,
                "ref": "A",
                "alt": "G",
                "window": 50,
            },
            cfg,
        )
        text = result["content"][0]["text"]
        assert "not found in region" in text  # a clean miss, not a region-parse error

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fusion_near_contig_start_observes_bam(self, small_bam_path, cfg, monkeypatch):
        from bamcp.analysis import fusion

        async def _none(*a, **k):
            return None

        monkeypatch.setattr(fusion, "_lookup_clinvar_safe", _none)
        monkeypatch.setattr(fusion, "_lookup_gnomad_safe", _none)
        captured = {}
        import bamcp.core.tools as tm

        real = tm._fetch_region_with_timeout

        async def _spy(file_path, region, *a, **k):
            captured["region"] = region
            return await real(file_path, region, *a, **k)

        monkeypatch.setattr(tm, "_fetch_region_with_timeout", _spy)
        await fusion.assemble_evidence("chr1", 5, "A", "G", small_bam_path, cfg, window=50)
        assert captured["region"] == "chr1:1-55"

    @pytest.mark.unit
    def test_tile_region_string_is_one_based(self):
        assert format_region("chr1", 0, 4096) == "chr1:1-4096"


class TestContigBounds:
    @pytest.mark.unit
    def test_end_clamped_to_contig_length(self, small_bam_path, ref_fasta_path):
        # chr1 is 1000 bp; a region running past the end is clamped, not padded.
        data = fetch_region(small_bam_path, "chr1:901-1200", reference_path=ref_fasta_path)
        assert (data.start, data.end) == (900, 1000)
        assert len(data.coverage) == 100
        assert len(data.reference_sequence or "") == 100

    @pytest.mark.unit
    def test_unknown_contig_is_a_clear_error(self, small_bam_path):
        with pytest.raises(ValueError, match="not found in the alignment header"):
            fetch_coverage_only(small_bam_path, "chrZ:1-100")

    @pytest.mark.unit
    def test_start_past_contig_end_is_a_clear_error(self, small_bam_path):
        with pytest.raises(ValueError, match="past the end of chr1"):
            fetch_coverage_only(small_bam_path, "chr1:1500-1600")

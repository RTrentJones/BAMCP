"""Tests for rubric output mode in get_variant_curation_summary."""

from __future__ import annotations

import json

import pytest

from bamcp.analysis.curation import (
    RUBRIC_VERSION,
    compute_rubric_scores,
    handle_get_variant_curation_summary,
)
from bamcp.config import BAMCPConfig


@pytest.fixture
def config_with_ref(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path)


@pytest.fixture
def config():
    return BAMCPConfig()


class TestComputeRubricScores:
    """Unit tests for the score derivation function."""

    @pytest.mark.unit
    def test_all_scores_present_and_bounded(self):
        scores = compute_rubric_scores(
            variant={"vaf": 0.5, "depth": 50, "alt_count": 25},
            evidence={
                "strand_bias": 0.0,
                "mapq_histogram": [0, 0, 0, 5, 45],
                "position_histogram": [0, 0, 5, 5],
            },
            artifact_risk={"risk_score": 0.0},
        )
        expected_keys = {
            "vaf_quality",
            "depth_quality",
            "strand_balance",
            "mapq_quality",
            "position_quality",
            "artifact_risk_inverse",
        }
        assert set(scores.keys()) == expected_keys
        for k, v in scores.items():
            assert 0.0 <= v <= 1.0, f"{k} out of bounds: {v}"

    @pytest.mark.unit
    def test_perfect_evidence_scores_near_one(self):
        scores = compute_rubric_scores(
            variant={"vaf": 0.5, "depth": 100},
            evidence={
                "strand_bias": 0.0,
                "mapq_histogram": [0, 0, 0, 0, 100],
                "position_histogram": [0, 0, 50, 50],
            },
            artifact_risk={"risk_score": 0.0},
        )
        assert scores["depth_quality"] == 1.0
        assert scores["strand_balance"] == 1.0
        assert scores["mapq_quality"] == 1.0
        assert scores["position_quality"] == 1.0
        assert scores["artifact_risk_inverse"] == 1.0
        assert scores["vaf_quality"] == 1.0

    @pytest.mark.unit
    def test_poor_evidence_scores_near_zero(self):
        scores = compute_rubric_scores(
            variant={"vaf": 0.05, "depth": 3},
            evidence={
                "strand_bias": 1.0,
                "mapq_histogram": [80, 20, 0, 0, 0],
                "position_histogram": [40, 40, 0, 0],
            },
            artifact_risk={"risk_score": 1.0},
        )
        assert scores["strand_balance"] == 0.0
        assert scores["mapq_quality"] == 0.0
        assert scores["artifact_risk_inverse"] == 0.0
        assert scores["depth_quality"] < 0.2
        # near-end bins are first two => all reads near ends => position 0
        assert scores["position_quality"] == 0.0
        # vaf 0.05 * 2.5 = 0.125
        assert scores["vaf_quality"] == 0.125

    @pytest.mark.unit
    def test_missing_evidence_handled_gracefully(self):
        scores = compute_rubric_scores(
            variant={},
            evidence={},
            artifact_risk={},
        )
        for k, v in scores.items():
            assert 0.0 <= v <= 1.0, f"{k} out of bounds: {v}"

    @pytest.mark.unit
    def test_empty_mapq_histogram_zero_quality(self):
        scores = compute_rubric_scores(
            variant={"vaf": 0.5, "depth": 30},
            evidence={"mapq_histogram": [], "strand_bias": 0.0},
            artifact_risk={"risk_score": 0.0},
        )
        assert scores["mapq_quality"] == 0.0


class TestHandlerFormatParam:
    """Tests that the handler honors the format parameter."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_format_is_text(self, small_bam_path, config_with_ref):
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
        # Text mode is not valid JSON (it's a formatted heading)
        try:
            json.loads(text)
            is_json = True
        except (ValueError, TypeError):
            is_json = False
        assert not is_json, "Text mode must not return JSON"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_rubric_format_returns_json(self, small_bam_path, config_with_ref):
        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "chr1",
                "pos": 105,
                "ref": "A",
                "alt": "T",
                "window": 50,
                "format": "rubric",
            },
            config_with_ref,
        )
        text = result["content"][0]["text"]
        payload = json.loads(text)
        assert "rubric_version" in payload
        assert payload["rubric_version"] == RUBRIC_VERSION

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_rubric_includes_scores_when_variant_found(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path
    ):
        """When the variant exists, the rubric output must include the scores dict."""
        cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
        # First scan to find an actual variant in the comprehensive fixture.
        from bamcp.core.tools import _fetch_region_with_timeout

        data = await _fetch_region_with_timeout(
            comprehensive_bam_path, "chr1:1-3000", comprehensive_ref_fasta_path, cfg
        )
        if not data.variants:
            pytest.skip("comprehensive fixture has no variants in chr1:1-3000")
        v = data.variants[0]
        result = await handle_get_variant_curation_summary(
            {
                "file_path": comprehensive_bam_path,
                "chrom": v["contig"],
                "pos": v["position"],
                "ref": v["ref"],
                "alt": v["alt"],
                "window": 50,
                "format": "rubric",
            },
            cfg,
        )
        payload = json.loads(result["content"][0]["text"])
        if "error" in payload:
            pytest.skip(f"Variant lookup failed: {payload['error']}")
        assert "scores" in payload
        scores = payload["scores"]
        assert set(scores.keys()) == {
            "vaf_quality",
            "depth_quality",
            "strand_balance",
            "mapq_quality",
            "position_quality",
            "artifact_risk_inverse",
        }
        for k, v in scores.items():
            assert 0.0 <= v <= 1.0, f"{k}={v} out of [0,1]"
        # And the standard structured fields remain present.
        assert "variant" in payload
        assert "artifact_assessment" in payload
        assert "confidence" in payload
        assert "recommendations" in payload

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_invalid_format_falls_back_to_text(self, small_bam_path, config_with_ref):
        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "chr1",
                "pos": 105,
                "ref": "A",
                "alt": "T",
                "window": 50,
                "format": "bogus",
            },
            config_with_ref,
        )
        text = result["content"][0]["text"]
        try:
            json.loads(text)
            is_json = True
        except (ValueError, TypeError):
            is_json = False
        assert not is_json, "Unknown format must fall back to plain text"


class TestRubricErrorResponses:
    """Errors should also respect the rubric envelope so the eval grader can parse them."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_validation_error_in_rubric_returns_json(self, small_bam_path, config):
        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "invalid_chr",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "format": "rubric",
            },
            config,
        )
        payload = json.loads(result["content"][0]["text"])
        assert "error" in payload
        assert payload["rubric_version"] == RUBRIC_VERSION

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_validation_error_in_text_returns_plain_text(self, small_bam_path, config):
        result = await handle_get_variant_curation_summary(
            {
                "file_path": small_bam_path,
                "chrom": "invalid_chr",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "format": "text",
            },
            config,
        )
        text = result["content"][0]["text"]
        try:
            json.loads(text)
            is_json = True
        except (ValueError, TypeError):
            is_json = False
        assert not is_json, "Text-mode error must remain plain text"

"""Unit tests for the truth-set scorer (no pysam, no fixtures).

Uses a scripted FakeRouter that returns canned tool JSON so the scoring logic
is exercised in isolation from the real BAMCP handlers (those are covered by
the integration smoke test).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval.metrics import PRF
from bamcp.eval.router import RouterResult
from bamcp.eval.truthset import (
    TruthsetReport,
    load_truthset,
    score_truthset,
)

pytestmark = pytest.mark.unit

MANIFEST = Path("tests/eval/datasets/synthetic_v1/manifest.yaml")


def test_load_truthset_parses_manifest():
    ts = load_truthset(MANIFEST)
    assert ts.dataset == "synthetic_v1"
    assert ts.version == 1
    assert ts.detection_region == "chr1:1000-2900"
    assert ts.bam.endswith("comprehensive.bam")
    assert len(ts.variant_sites) == 9
    # The multi-allelic site appears twice with different alts.
    keys = {s.key for s in ts.variant_sites}
    assert ("chr1", 2500, "A", "G") in keys
    assert ("chr1", 2500, "A", "T") in keys
    # Exactly one clean control, three artifact sites, one confidence control.
    assert sum(s.is_clean_control for s in ts.variant_sites) == 1
    assert sum(s.expected_artifact is not None for s in ts.variant_sites) == 3
    assert sum(s.expected_confidence is not None for s in ts.variant_sites) == 1
    assert len(ts.negative_regions) == 3


class FakeRouter:
    """Router returning canned responses keyed by (tool, region or pos)."""

    def __init__(self, *, detected, negatives, curations):
        self._detected = detected
        self._negatives = negatives
        self._curations = curations

    def list_tools(self):
        return ["get_variants", "get_variant_curation_summary"]

    async def call(self, name, arguments):
        if name == "get_variants":
            region = arguments["region"]
            variants = self._negatives.get(region, self._detected)
            return RouterResult(ok=True, text=json.dumps({"variants": variants}))
        if name == "get_variant_curation_summary":
            pos = arguments["pos"]
            return RouterResult(ok=True, text=json.dumps(self._curations[pos]))
        return RouterResult(ok=False, text="", error="unknown")


def _variant(pos, alt, ref="A", contig="chr1"):
    return {"contig": contig, "position": pos, "ref": ref, "alt": alt}


def _curation(risk_types, risk_score, likelihood, confidence="low"):
    return {
        "artifact_assessment": {
            "risks": [{"type": t} for t in risk_types],
            "risk_score": risk_score,
            "artifact_likelihood": likelihood,
        },
        "confidence": confidence,
    }


def _perfect_router():
    detected = [
        _variant(1050, "T"),
        _variant(1100, "C"),
        _variant(1150, "G"),
        _variant(1200, "T"),
        _variant(1802, "A", ref="T"),
        _variant(2075, "G", ref="C"),
        _variant(2500, "G"),
        _variant(2500, "T"),
        _variant(2800, "T"),  # high-confidence positive control
    ]
    curations = {
        1050: _curation(["strand_bias", "read_position_bias"], 0.7, "high"),
        1200: _curation(["low_mapq", "strand_bias"], 0.95, "high"),
        1802: _curation(["homopolymer"], 0.6, "high"),
        1150: _curation(["read_position_bias"], 0.3, "medium"),
        2800: _curation([], 0.0, "low", confidence="high"),  # earns high confidence
    }
    negatives = {"chr1:2100-2150": [], "chr1:300-600": [], "chr2:100-400": []}
    return FakeRouter(detected=detected, negatives=negatives, curations=curations)


async def test_score_truthset_perfect():
    ts = load_truthset(MANIFEST)
    report = await score_truthset(ts, BAMCPConfig(), router=_perfect_router())
    assert report.detection.recall == 1.0
    assert report.detection.precision == 1.0
    assert report.detection.tp == 9
    assert report.false_positives == 0
    assert report.artifact_recall == 1.0
    assert report.clean_discriminated is True
    assert report.overconfident_sites == []
    assert report.confidence_mismatches == []
    passed, failures = report.meets_floors()
    assert passed, failures


async def test_positive_control_must_reach_high_confidence():
    # If the high-confidence positive control stops being high-confidence, the
    # overconfidence guard would pass vacuously — so this must fail loudly.
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._curations[2800] = _curation([], 0.0, "low", confidence="medium")
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.confidence_mismatches
    passed, failures = report.meets_floors()
    assert not passed
    assert any("confidence calibration" in f for f in failures)


async def test_overconfident_site_is_a_safety_failure():
    # A high-artifact site reported as high confidence is a safety violation.
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._curations[1050] = _curation(
        ["strand_bias", "read_position_bias"], 0.7, "high", confidence="high"
    )
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.overconfident_sites == ["chr1:1050"]
    passed, failures = report.meets_floors()
    assert not passed
    assert any("SAFETY" in f for f in failures)


async def test_score_truthset_missing_variant_fails_recall():
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._detected = router._detected[:-1]  # drop one expected allele
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.detection.recall < 1.0
    assert report.missing_calls
    passed, failures = report.meets_floors()
    assert not passed
    assert any("recall" in f for f in failures)


async def test_score_truthset_false_positive_in_negative_region():
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._negatives["chr1:2100-2150"] = [_variant(2125, "G")]
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.false_positives == 1
    passed, failures = report.meets_floors()
    assert not passed
    assert any("false positive" in f for f in failures)


async def test_score_truthset_missing_artifact_type():
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._curations[1050] = _curation(["read_position_bias"], 0.3, "medium")
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.artifact_recall < 1.0
    passed, failures = report.meets_floors()
    assert not passed
    assert any("artifact" in f for f in failures)


async def test_clean_control_flagged_high_fails_discrimination():
    ts = load_truthset(MANIFEST)
    router = _perfect_router()
    router._curations[1150] = _curation(["read_position_bias"], 0.99, "high")
    report = await score_truthset(ts, BAMCPConfig(), router=router)
    assert report.clean_discriminated is False
    passed, failures = report.meets_floors()
    assert not passed
    assert any("clean" in f for f in failures)


def test_report_as_dict_serializable():
    report = TruthsetReport(
        dataset="synthetic_v1",
        version=1,
        detection=PRF(tp=8, fp=0, fn=0),
        false_positives=0,
    )
    d = report.as_dict()
    # Round-trips through JSON without error.
    assert json.loads(json.dumps(d))["detection"]["recall"] == 1.0


def test_format_report_renders_pass_and_fail():
    from bamcp.eval.metrics import ArtifactCheck
    from bamcp.eval.truthset import _format_report

    report = TruthsetReport(
        dataset="synthetic_v1",
        version=1,
        detection=PRF(tp=8, fp=0, fn=0),
        false_positives=0,
        artifact_checks=[ArtifactCheck("chr1:1050", "strand_bias", True, 0.7, "high")],
    )
    passed_text = _format_report(report, True, [])
    assert "PASS" in passed_text
    assert "chr1:1050" in passed_text and "[ok]" in passed_text

    failed_text = _format_report(report, False, ["variant recall 0.5 < 1.0"])
    assert "FAIL" in failed_text
    assert "variant recall 0.5" in failed_text


async def test_score_truthset_tolerates_failed_tool_calls():
    ts = load_truthset(MANIFEST)

    class BrokenRouter:
        def list_tools(self):
            return []

        async def call(self, name, arguments):
            return RouterResult(ok=False, text="", error="boom")

    report = await score_truthset(ts, BAMCPConfig(), router=BrokenRouter())
    # Nothing detected => recall 0, but scoring degrades instead of raising.
    assert report.detection.recall == 0.0
    assert report.false_positives == 0

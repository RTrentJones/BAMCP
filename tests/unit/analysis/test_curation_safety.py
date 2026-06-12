"""Safety regression tests for curation output.

These lock in the intended-use boundary: BAMCP supports human curation, it does
not make clinical determinations. A regression that re-introduces over-claiming
language or drops the disclaimer fails here. See SAFETY.md.
"""

from __future__ import annotations

import pytest

from bamcp.analysis.curation import (
    INTENDED_USE,
    format_curation_summary,
    generate_curator_recommendations,
)

pytestmark = pytest.mark.unit


def _clean_inputs():
    """Inputs that trip none of the artifact flags (the 'no concerns' path)."""
    variant = {"depth": 100, "vaf": 0.5, "alt_count": 50}
    evidence = {
        "strand_bias": 0.5,
        "position_histogram": [],
        "mapq_histogram": [],
    }
    artifact_risk = {"artifact_likelihood": "low", "risk_score": 0.1, "risks": []}
    return variant, evidence, artifact_risk


def test_intended_use_is_non_diagnostic():
    text = INTENDED_USE.lower()
    assert "not a diagnostic device" in text
    assert "research" in text


def test_no_concerns_does_not_claim_clinical_suitability():
    recs = generate_curator_recommendations(*_clean_inputs())
    joined = " ".join(recs).lower()
    # The old over-claim must not come back.
    assert "suitable for clinical interpretation" not in joined
    # And the replacement makes the boundary explicit.
    assert "manual curation is still required" in joined
    assert "not a clinical determination" in joined


def test_low_coverage_still_flagged():
    variant, evidence, artifact_risk = _clean_inputs()
    variant["depth"] = 8  # below the 20x threshold
    recs = generate_curator_recommendations(variant, evidence, artifact_risk)
    assert any("LOW COVERAGE" in r for r in recs)


def test_text_summary_carries_disclaimer():
    summary = {
        "variant": {
            "location": "chr1:1050",
            "change": "A>T",
            "vaf": 0.4,
            "depth": 20,
            "alt_count": 8,
        },
        "quality_metrics": {
            "mean_quality": 35.0,
            "median_quality": 35.0,
            "quality_distribution": [],
        },
        "strand_analysis": {"forward": 8, "reverse": 0, "strand_bias": 1.0, "is_balanced": False},
        "position_in_read": {"distribution": [], "near_end_fraction": 0.0},
        "artifact_assessment": {"artifact_likelihood": "high", "risk_score": 0.7, "risks": []},
        "confidence": "low",
        "recommendations": ["HIGH ARTIFACT RISK: ..."],
        "intended_use": INTENDED_USE,
    }
    text = format_curation_summary(summary)
    assert "INTENDED USE" in text
    assert "not a diagnostic device" in text.lower()

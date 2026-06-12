"""Unit tests for the deterministic eval metrics."""

from __future__ import annotations

import pytest

from bamcp.eval.metrics import (
    PRF,
    ArtifactCheck,
    artifact_recall,
    prf_from_sets,
    risk_types_present,
)

pytestmark = pytest.mark.unit


def test_prf_basic_counts():
    prf = PRF(tp=8, fp=2, fn=1)
    assert prf.precision == pytest.approx(8 / 10)
    assert prf.recall == pytest.approx(8 / 9)
    assert prf.f1 == pytest.approx(2 * (0.8 * (8 / 9)) / (0.8 + 8 / 9))


def test_prf_empty_is_perfect():
    # No expected and no predicted => nothing wrong => precision/recall 1.0.
    prf = PRF(tp=0, fp=0, fn=0)
    assert prf.precision == 1.0
    assert prf.recall == 1.0
    assert prf.f1 == 1.0  # both perfect => harmonic mean is also 1.0


def test_prf_f1_zero_when_no_true_positives():
    prf = PRF(tp=0, fp=3, fn=2)
    assert prf.f1 == 0.0


def test_prf_as_dict_rounds():
    prf = PRF(tp=2, fp=1, fn=0)
    d = prf.as_dict()
    assert d["tp"] == 2 and d["fp"] == 1 and d["fn"] == 0
    assert d["precision"] == pytest.approx(0.6667, abs=1e-4)
    assert d["recall"] == 1.0


def test_prf_from_sets_intersection():
    expected = {("chr1", 1050, "A", "T"), ("chr1", 1150, "A", "G")}
    predicted = {("chr1", 1050, "A", "T"), ("chr1", 9999, "C", "G")}
    prf = prf_from_sets(expected, predicted)
    assert prf.tp == 1
    assert prf.fp == 1
    assert prf.fn == 1


def test_artifact_recall_all_found():
    checks = [
        ArtifactCheck("chr1:1050", "strand_bias", True, 0.7, "high"),
        ArtifactCheck("chr1:1200", "low_mapq", True, 0.95, "high"),
    ]
    assert artifact_recall(checks) == 1.0


def test_artifact_recall_partial():
    checks = [
        ArtifactCheck("chr1:1050", "strand_bias", True, 0.7, "high"),
        ArtifactCheck("chr1:1802", "homopolymer", False, 0.2, "low"),
    ]
    assert artifact_recall(checks) == 0.5


def test_artifact_recall_empty_is_one():
    assert artifact_recall([]) == 1.0


def test_risk_types_present_extracts_types():
    assessment = {
        "risks": [
            {"type": "strand_bias", "severity": "high"},
            {"type": "read_position_bias", "severity": "medium"},
        ]
    }
    assert risk_types_present(assessment) == {"strand_bias", "read_position_bias"}


@pytest.mark.parametrize(
    "bad",
    [None, {}, {"risks": "nope"}, {"risks": [{"severity": "high"}]}, {"risks": [42]}, "string"],
)
def test_risk_types_present_tolerates_garbage(bad):
    assert risk_types_present(bad) == set()


def test_artifact_check_as_dict():
    c = ArtifactCheck("chr1:1050", "strand_bias", True, 0.7, "high")
    assert c.as_dict() == {
        "site": "chr1:1050",
        "expected_type": "strand_bias",
        "found": True,
        "risk_score": 0.7,
        "likelihood": "high",
    }

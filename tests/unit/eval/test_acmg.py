"""Unit tests for the ACMG classification scoring (pure functions)."""

from __future__ import annotations

import pytest

from bamcp.eval.acmg import (
    AcmgCase,
    build_report,
    classification_index,
    extract_classification,
    load_acmg_cases,
    normalize_classification,
    score_classification,
)


@pytest.mark.unit
@pytest.mark.parametrize(
    "raw,expected",
    [
        ("Pathogenic", "Pathogenic"),
        ("likely pathogenic", "Likely Pathogenic"),
        ("VUS", "Uncertain Significance"),
        ("uncertain significance (vus)", "Uncertain Significance"),
        ("Likely Benign", "Likely Benign"),
        ("benign", "Benign"),
        ("nonsense", None),
    ],
)
def test_normalize_classification(raw, expected):
    assert normalize_classification(raw) == expected


@pytest.mark.unit
def test_classification_index_ordering():
    assert classification_index("Benign") == 0
    assert classification_index("Pathogenic") == 4
    assert classification_index("VUS") == 2
    assert classification_index("garbage") is None


@pytest.mark.unit
def test_extract_classification_prefers_final_line():
    text = "I considered options... final_classification: Likely Pathogenic"
    assert extract_classification(text) == "Likely Pathogenic"


@pytest.mark.unit
def test_extract_classification_falls_back_to_last_mention():
    text = "This could be benign, but on balance it is Pathogenic."
    assert extract_classification(text) == "Pathogenic"


@pytest.mark.unit
def test_extract_classification_none_when_absent():
    assert extract_classification("no verdict here") is None
    assert extract_classification("") is None


@pytest.mark.unit
def test_score_exact_match():
    s = score_classification("Pathogenic", "Pathogenic")
    assert s.exact and s.within_one and s.distance == 0 and s.failure_mode == "none"


@pytest.mark.unit
def test_score_overcall_and_undercall():
    over = score_classification("Pathogenic", "Uncertain Significance")
    assert not over.exact and over.distance == 2 and over.failure_mode == "overcall"
    assert over.within_one is False

    under = score_classification("Likely Benign", "Pathogenic")
    assert under.failure_mode == "undercall" and under.distance == -3


@pytest.mark.unit
def test_score_within_one():
    s = score_classification("Likely Pathogenic", "Pathogenic")
    assert not s.exact and s.within_one and s.distance == -1


@pytest.mark.unit
def test_score_unparseable_prediction():
    s = score_classification(None, "Benign")
    assert s.failure_mode == "unparseable" and s.distance is None and not s.within_one


@pytest.mark.unit
def test_score_raises_on_bad_expected():
    with pytest.raises(ValueError):
        score_classification("Benign", "not-a-label")


@pytest.mark.unit
def test_build_report_aggregates():
    cases = [
        AcmgCase("chr1", 1, "A", "T", "Pathogenic"),
        AcmgCase("chr2", 2, "C", "G", "Benign"),
        AcmgCase("chr3", 3, "G", "A", "Uncertain Significance"),
    ]
    scores = [
        score_classification("Pathogenic", "Pathogenic"),  # exact
        score_classification("Likely Benign", "Benign"),  # within one
        score_classification("Pathogenic", "Uncertain Significance"),  # overcall, not within one
    ]
    report = build_report(zip(cases, scores, strict=True))
    assert report.total == 3
    assert report.exact_matches == 1
    assert report.within_one == 2
    # Likely Benign vs Benign (+1) and Pathogenic vs VUS (+2) are both overcalls.
    assert report.failure_modes.get("overcall") == 2
    assert report.failure_modes.get("none") == 1
    assert report.exact_accuracy == pytest.approx(1 / 3)


@pytest.mark.unit
def test_load_acmg_cases(tmp_path):
    manifest = tmp_path / "m.yaml"
    manifest.write_text(
        "dataset: t\nversion: 1\ncases:\n"
        "  - chrom: chr17\n    pos: 7674220\n    ref: C\n    alt: T\n"
        "    gene: TP53\n    expected_classification: Pathogenic\n    tier: tier1_clear\n"
    )
    cases = load_acmg_cases(manifest)
    assert len(cases) == 1
    assert cases[0].gene == "TP53"
    assert cases[0].expected_classification == "Pathogenic"


@pytest.mark.unit
def test_shipped_acmg_manifest_is_valid():
    """The committed acmg_v1 manifest must parse and carry valid labels."""
    cases = load_acmg_cases("tests/eval/datasets/acmg_v1/manifest.yaml")
    assert len(cases) >= 5
    for c in cases:
        assert classification_index(c.expected_classification) is not None, (
            c.expected_classification
        )

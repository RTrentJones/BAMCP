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


# ── Review fixes: extraction ────────────────────────────────────────────────


@pytest.mark.unit
@pytest.mark.parametrize(
    "raw,expected",
    [
        ("not pathogenic", None),
        ("non-pathogenic", None),
        ("is not considered benign", None),
        ("Likely-Pathogenic", "Likely Pathogenic"),
        ("likely  benign", "Likely Benign"),
    ],
)
def test_normalize_handles_negation_and_spelling(raw, expected):
    assert normalize_classification(raw) == expected


@pytest.mark.unit
@pytest.mark.parametrize(
    "text,expected",
    [
        # Fallback scan: the two-word label must win over the word it contains.
        ("After weighing the evidence this is likely pathogenic.", "Likely Pathogenic"),
        ("Population data make it likely benign.", "Likely Benign"),
        # A negated mention is not a verdict.
        ("It is likely benign and not pathogenic.", "Likely Benign"),
        ("Overall VUS; this is not pathogenic.", "Uncertain Significance"),
    ],
)
def test_extract_last_mention_is_alias_and_negation_aware(text, expected):
    assert extract_classification(text) == expected


@pytest.mark.unit
def test_extract_takes_last_explicit_classification_over_echoed_tool_output():
    text = (
        'classify_variant returned {"classification": "Pathogenic", "criteria": []}\n'
        "Reviewing population frequency, I disagree.\n"
        "classification: Likely Benign"
    )
    assert extract_classification(text) == "Likely Benign"


@pytest.mark.unit
def test_extract_final_classification_beats_later_plain_classification():
    text = "final_classification: **Likely Pathogenic**\n(tool classification: Benign)"
    assert extract_classification(text) == "Likely Pathogenic"


# ── Review fixes: the gate must beat a constant-VUS answer ──────────────────


def _shipped_cases():
    return load_acmg_cases("tests/eval/datasets/acmg_v1/manifest.yaml")


@pytest.mark.unit
def test_constant_vus_baseline_fails_the_gate():
    from bamcp.eval.acmg import constant_baseline

    cases = _shipped_cases()
    baseline = constant_baseline(cases)
    # The old within-one >= 0.5 floor passed this degenerate classifier.
    assert baseline.within_one_accuracy >= 0.5
    passed, failures = baseline.meets_floors(baseline)
    assert not passed
    assert any("baseline" in f for f in failures)
    assert any("anchors" in f for f in failures)


@pytest.mark.unit
def test_perfect_classifier_passes_gate_and_reports_anchors():
    from bamcp.eval.acmg import constant_baseline

    cases = _shipped_cases()
    report = build_report(
        (c, score_classification(c.expected_classification, c.expected_classification))
        for c in cases
    )
    assert report.anchors_total >= 4 and report.anchors_correct == report.anchors_total
    passed, failures = report.meets_floors(constant_baseline(cases))
    assert passed, failures
    assert report.as_dict()["anchors_correct"] == report.anchors_total


@pytest.mark.unit
def test_anchor_on_wrong_side_fails_even_with_high_exact_accuracy():
    from bamcp.eval.acmg import constant_baseline

    cases = _shipped_cases()
    scored = []
    flipped = False
    for c in cases:
        pred = c.expected_classification
        if not flipped and normalize_classification(pred) == "Pathogenic":
            pred, flipped = "Uncertain Significance", True
        scored.append((c, score_classification(pred, c.expected_classification)))
    report = build_report(scored)
    assert report.exact_accuracy > 0.5
    passed, failures = report.meets_floors(constant_baseline(cases))
    assert not passed
    assert failures == [
        f"B/P anchors on the wrong side: {report.anchors_total - 1}/{report.anchors_total} correct"
    ]


@pytest.mark.unit
def test_empty_report_fails_gate():
    from bamcp.eval.acmg import constant_baseline

    empty = build_report([])
    assert empty.meets_floors(constant_baseline([]))[0] is False


@pytest.mark.unit
def test_format_report_includes_baseline():
    from bamcp.eval.acmg import _format_report, constant_baseline

    cases = _shipped_cases()
    baseline = constant_baseline(cases)
    out = _format_report(baseline, baseline)
    assert "constant-VUS baseline" in out and "B/P anchors" in out

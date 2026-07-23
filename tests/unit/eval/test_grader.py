"""Unit tests for eval graders."""

from __future__ import annotations

import pytest

from bamcp.eval.grader import collect_rubric_payloads, grade_case, judge_with_llm
from bamcp.eval.providers import MockProvider, ProviderResponse
from bamcp.eval.schema import EvalCase


@pytest.mark.unit
def test_tools_only_case_defers_to_judge_not_auto_pass():
    """Tool selection is necessary, not sufficient: with no answer check, defer to the judge.

    Calling the expected tools no longer deterministically passes a case — otherwise a model
    could call the tool and return an empty/wrong answer and still pass.
    """
    case = EvalCase(
        name="x",
        input="run",
        expected="(any)",
        tools_expected=["get_variants", "get_coverage"],
    )
    # Tools satisfied but no expected_claims / score thresholds → None (judge grades the answer).
    assert grade_case(case, "response", ["get_variants", "get_coverage"], []) is None


@pytest.mark.unit
def test_grade_case_fails_when_a_tool_missing():
    case = EvalCase(
        name="x",
        input="run",
        expected="(any)",
        tools_expected=["get_variants", "lookup_clinvar"],
    )
    v = grade_case(case, "response", ["get_variants"], [])
    assert v is not None
    assert v.passed is False
    assert "lookup_clinvar" in v.rationale


@pytest.mark.unit
def test_grade_case_no_deterministic_rule_returns_none():
    case = EvalCase(name="x", input="run", expected="just text")
    assert grade_case(case, "response", [], []) is None


@pytest.mark.unit
def test_tools_satisfied_but_answer_wrong_fails():
    """The core fix: calling the expected tool no longer passes a case whose answer is wrong."""
    case = EvalCase(
        name="x",
        input="What is the VAF?",
        expected="The VAF is 43.5%",
        tools_expected=["get_variants"],
        expected_claims=["43.5%"],
    )
    # Right tool called, but the answer omits the required numeric claim.
    v = grade_case(case, "I ran the analysis.", ["get_variants"], [])
    assert v is not None
    assert v.passed is False
    assert v.tool_score == 1.0  # tool selection was correct...
    assert v.answer_score == 0.0  # ...but the answer was not
    assert "43.5%" in v.rationale


@pytest.mark.unit
def test_tools_and_claims_both_satisfied_passes():
    case = EvalCase(
        name="x",
        input="What is the VAF?",
        expected="The VAF is 43.5%",
        tools_expected=["get_variants"],
        expected_claims=["43.5%"],
    )
    v = grade_case(case, "The variant VAF is 43.5% at that site.", ["get_variants"], [])
    assert v is not None
    assert v.passed is True
    assert v.tool_score == 1.0
    assert v.answer_score == 1.0


@pytest.mark.unit
def test_expected_claims_alone_grade_answer_correctness():
    case = EvalCase(name="x", input="q", expected="pathogenic", expected_claims=["pathogenic"])
    assert grade_case(case, "This variant is pathogenic.", [], []).passed is True
    miss = grade_case(case, "This variant is benign.", [], [])
    assert miss.passed is False
    assert miss.answer_score == 0.0


@pytest.mark.unit
def test_expected_claims_match_on_token_boundaries():
    """A short token claim must not match inside a larger token (chr1 != chr10)."""
    case = EvalCase(name="x", input="q", expected="chr1", expected_claims=["chr1"])
    # 'chr1' as a bounded token → pass.
    assert grade_case(case, "Reads on chr1: 100-200", [], []).passed is True
    # 'chr1' only appears inside 'chr10' → fail (was a false positive under substring match).
    assert grade_case(case, "The contig is chr10 here", [], []).passed is False


@pytest.mark.unit
def test_expected_claims_numeric_token_boundary():
    case = EvalCase(name="x", input="q", expected="43.5%", expected_claims=["43.5%"])
    assert grade_case(case, "VAF is 43.5% at that site", [], []).passed is True
    assert grade_case(case, "VAF is 143.5% (bogus)", [], []).passed is False


@pytest.mark.unit
def test_grade_case_threshold_passes():
    case = EvalCase(
        name="x",
        input="run",
        expected="score_at_least: depth_quality=0.5",
    )
    rubrics = [{"rubric_version": "1.0", "scores": {"depth_quality": 0.9}}]
    v = grade_case(case, "response", [], rubrics)
    assert v is not None and v.passed is True


@pytest.mark.unit
def test_grade_case_threshold_fails_when_below():
    case = EvalCase(name="x", input="run", expected="score_at_least: depth_quality=0.9")
    rubrics = [{"rubric_version": "1.0", "scores": {"depth_quality": 0.5}}]
    v = grade_case(case, "response", [], rubrics)
    assert v is not None and v.passed is False


@pytest.mark.unit
def test_grade_case_threshold_missing_score():
    case = EvalCase(name="x", input="run", expected="score_at_least: depth_quality=0.5")
    rubrics: list[dict] = [{"rubric_version": "1.0", "scores": {}}]
    v = grade_case(case, "response", [], rubrics)
    assert v is not None and v.passed is False
    assert "not present" in v.rationale.lower()


@pytest.mark.unit
def test_collect_rubric_payloads_extracts_only_rubric_outputs():
    pairs = [
        ("get_variant_curation_summary", '{"rubric_version":"1.0","scores":{"a":0.5}}'),
        ("get_variant_curation_summary", "plain text not json"),
        ("get_variants", '{"count": 3}'),  # ignored — wrong tool name
        ("get_variant_curation_summary", '{"no_version":true}'),  # ignored — no version
    ]
    out = collect_rubric_payloads(pairs)
    assert len(out) == 1
    assert out[0]["scores"]["a"] == 0.5


@pytest.mark.unit
async def test_judge_with_llm_parses_pass_verdict():
    def script(messages, tools):
        return ProviderResponse(text="VERDICT: PASS — looks right")

    judge = MockProvider(script)
    case = EvalCase(name="x", input="q", expected="e")
    v = await judge_with_llm(case, "response text", judge)
    assert v.passed is True
    assert v.deterministic is False


@pytest.mark.unit
async def test_judge_with_llm_parses_fail_verdict():
    def script(messages, tools):
        return ProviderResponse(text="VERDICT: FAIL — wrong")

    judge = MockProvider(script)
    case = EvalCase(name="x", input="q", expected="e")
    v = await judge_with_llm(case, "response text", judge)
    assert v.passed is False

"""Unit tests for eval graders."""

from __future__ import annotations

import pytest

from bamcp.eval.grader import collect_rubric_payloads, grade_case, judge_with_llm
from bamcp.eval.providers import MockProvider, ProviderResponse
from bamcp.eval.schema import EvalCase


@pytest.mark.unit
def test_grade_case_passes_when_all_tools_called():
    case = EvalCase(
        name="x",
        input="run",
        expected="(any)",
        tools_expected=["get_variants", "get_coverage"],
    )
    v = grade_case(case, "response", ["get_variants", "get_coverage"], [])
    assert v is not None
    assert v.passed is True
    assert v.deterministic is True


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

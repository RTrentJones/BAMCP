"""Unit tests for the eval runner end-to-end (with MockProvider)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval.providers import MockProvider, ProviderResponse, ToolCall
from bamcp.eval.router import InProcessRouter
from bamcp.eval.runner import run_case, run_cases
from bamcp.eval.schema import EvalCase


@pytest.fixture
def cfg(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path)


@pytest.mark.unit
async def test_run_case_invokes_tool_and_grades_pass(tmp_path: Path, small_bam_path, cfg):
    """A scripted run that calls get_coverage should pass deterministically."""

    state = {"turns": 0}

    def script(messages, tools):
        state["turns"] += 1
        if state["turns"] == 1:
            return ProviderResponse(
                text="Calling get_coverage",
                tool_calls=[
                    ToolCall(
                        id="c1",
                        name="get_coverage",
                        arguments={"file_path": small_bam_path, "region": "chr1:90-200"},
                    )
                ],
            )
        return ProviderResponse(text="Mean coverage reported")

    provider = MockProvider(script)
    case = EvalCase(
        name="coverage test",
        input="run coverage",
        expected="mean coverage",
        tools_expected=["get_coverage"],
        # Answer correctness is checked too (tool selection alone no longer passes a case);
        # the mock's final answer must contain this claim.
        expected_claims=["coverage"],
    )
    router = InProcessRouter(config=cfg)

    result = await run_case(case, provider, router, judge=None, output_dir=tmp_path)
    assert result.verdict.passed is True
    assert result.verdict.deterministic is True
    assert result.verdict.tool_score == 1.0
    assert result.verdict.answer_score == 1.0
    assert "get_coverage" in result.tool_calls
    # Telemetry JSONL must exist with at least one event.
    assert Path(result.telemetry_path).exists()
    events = [
        json.loads(line) for line in Path(result.telemetry_path).read_text().splitlines() if line
    ]
    assert any(e.get("tool_name") == "get_coverage" for e in events)


@pytest.mark.unit
async def test_run_case_missing_tool_fails(tmp_path: Path, small_bam_path, cfg):
    """When the LLM doesn't call the expected tool, the deterministic grader fails."""

    def script(messages, tools):
        return ProviderResponse(text="I refuse to call any tools")

    case = EvalCase(
        name="no tools",
        input="run coverage",
        expected="anything",
        tools_expected=["get_coverage"],
    )
    router = InProcessRouter(config=cfg)
    result = await run_case(case, MockProvider(script), router, judge=None, output_dir=tmp_path)
    assert result.verdict.passed is False
    assert "get_coverage" in result.verdict.rationale


@pytest.mark.unit
async def test_run_case_falls_back_to_judge(tmp_path: Path, cfg):
    """A free-form case with no deterministic rule falls back to the LLM judge."""

    def model_script(messages, tools):
        return ProviderResponse(text="An open-ended answer.")

    def judge_script(messages, tools):
        return ProviderResponse(text="VERDICT: PASS — sounds right")

    case = EvalCase(name="open", input="open question", expected="open answer")
    result = await run_case(
        case,
        MockProvider(model_script),
        InProcessRouter(config=cfg),
        judge=MockProvider(judge_script),
        output_dir=tmp_path,
    )
    assert result.verdict.passed is True
    assert result.verdict.deterministic is False


@pytest.mark.unit
async def test_run_case_no_judge_no_rule_fails(tmp_path: Path, cfg):
    """Without a deterministic rule and no judge, the case fails gracefully."""

    def script(messages, tools):
        return ProviderResponse(text="answer")

    case = EvalCase(name="open", input="q", expected="a")
    result = await run_case(
        case,
        MockProvider(script),
        InProcessRouter(config=cfg),
        judge=None,
        output_dir=tmp_path,
    )
    assert result.verdict.passed is False
    assert "judge" in result.verdict.rationale.lower()


@pytest.mark.unit
async def test_run_cases_rendering_modes_iterate(tmp_path: Path, cfg):
    """With --with-rendering-comparison, the same case runs once per mode."""

    def script(messages, tools):
        return ProviderResponse(text="ok")

    def judge(messages, tools):
        return ProviderResponse(text="VERDICT: PASS — fine")

    case = EvalCase(
        name="modes",
        input="x",
        expected="any",
        rendering_modes=["expanded", "dv-strips"],
    )
    results = await run_cases(
        cases=[case],
        provider=MockProvider(script),
        router=InProcessRouter(config=cfg),
        judge=MockProvider(judge),
        output_dir=tmp_path,
        with_rendering_comparison=True,
    )
    assert len(results) == 2
    modes = {r.rendering_mode for r in results}
    assert modes == {"expanded", "dv-strips"}


@pytest.mark.unit
async def test_run_case_handler_error_is_surfaced(tmp_path: Path, cfg):
    """Tool errors surface to the LLM but don't crash the runner."""

    def script(messages, tools):
        return ProviderResponse(
            text="boom",
            tool_calls=[
                ToolCall(
                    id="c1",
                    name="nonexistent_tool",
                    arguments={},
                )
            ],
        )

    case = EvalCase(name="bad", input="x", expected="y", tools_expected=["get_coverage"])
    result = await run_case(
        case,
        MockProvider(script),
        InProcessRouter(config=cfg),
        judge=None,
        output_dir=tmp_path,
    )
    # Tool not found is surfaced as an error result; deterministic grader still
    # fails because get_coverage was never called.
    assert result.verdict.passed is False


@pytest.mark.unit
async def test_run_case_turn_limit_is_a_failure(tmp_path: Path, small_bam_path, cfg):
    """A model still calling tools at MAX_TURNS gave no final answer: fail, even if the
    partial transcript would satisfy the deterministic checks."""

    def script(messages, tools):
        return ProviderResponse(
            text="coverage looks fine",
            tool_calls=[
                ToolCall(
                    id="c1",
                    name="get_coverage",
                    arguments={"file_path": small_bam_path, "region": "chr1:90-200"},
                )
            ],
        )

    case = EvalCase(
        name="loop",
        input="x",
        expected="y",
        tools_expected=["get_coverage"],
        expected_claims=["coverage"],
    )
    result = await run_case(
        case, MockProvider(script), InProcessRouter(config=cfg), judge=None, output_dir=tmp_path
    )
    assert result.verdict.passed is False
    assert "turn limit" in result.verdict.rationale
    assert result.error and "turn limit" in result.error


@pytest.mark.unit
async def test_run_case_provider_exception_is_a_failure(tmp_path: Path, cfg):
    def script(messages, tools):
        raise RuntimeError("provider down")

    # No checks declared: a judge would otherwise be consulted on an empty answer.
    case = EvalCase(name="boom", input="x", expected="y", expected_claims=[])
    judge = MockProvider(lambda m, t: ProviderResponse(text="VERDICT: PASS — fine"))
    result = await run_case(
        case, MockProvider(script), InProcessRouter(config=cfg), judge=judge, output_dir=tmp_path
    )
    assert result.verdict.passed is False
    assert "provider down" in result.verdict.rationale


@pytest.mark.unit
async def test_failed_tool_call_earns_no_tool_selection_credit(tmp_path: Path, cfg):
    """get_coverage on a missing file raises: it was attempted, not successfully used."""
    state = {"turns": 0}

    def script(messages, tools):
        state["turns"] += 1
        if state["turns"] == 1:
            return ProviderResponse(
                text="",
                tool_calls=[
                    ToolCall(
                        id="c1",
                        name="get_coverage",
                        arguments={"file_path": "/nonexistent.bam", "region": "chr1:1-10"},
                    )
                ],
            )
        return ProviderResponse(text="coverage reported")

    case = EvalCase(
        name="bad-call",
        input="x",
        expected="y",
        tools_expected=["get_coverage"],
        expected_claims=["coverage"],
    )
    result = await run_case(
        case, MockProvider(script), InProcessRouter(config=cfg), judge=None, output_dir=tmp_path
    )
    assert "get_coverage" in result.tool_calls  # attempted, still shown in the report
    assert result.verdict.passed is False
    assert result.verdict.tool_score == 0.0


@pytest.mark.unit
def test_call_ok_reads_error_payloads():
    from bamcp.eval.runner import _call_ok

    assert _call_ok({"ok": True, "result_summary": {"count": 1}})
    assert not _call_ok({"ok": False})
    assert not _call_ok({"ok": True, "result_summary": {"error": "not found"}})

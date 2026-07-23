"""Tests for the model-comparison + tool-ablation reporting.

The aggregation and Markdown rendering are pure and tested directly. The
end-to-end run is exercised with the ``mock`` provider so the full pipeline —
including the tool-vs-no-tool contrast — is proven without an API key.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bamcp.eval.compare import (
    ComparisonReport,
    RunSpec,
    SpecResult,
    _extract_args,
    aggregate,
    format_markdown,
    realistic_mock,
    run_comparison,
)
from bamcp.eval.schema import EvalResult, GraderVerdict, load_cases

pytestmark = pytest.mark.unit


def _result(category, passed, skipped=False):
    return EvalResult(
        case_name="c",
        category=category,
        rendering_mode=None,
        response_text="",
        tool_calls=[],
        duration_ms=0.0,
        verdict=GraderVerdict(passed=passed, rationale="", deterministic=True),
        skipped=skipped,
    )


def test_runspec_slug():
    assert RunSpec("Opus 4.8 (tools)", "anthropic", "claude-opus-4-8").slug == "opus-4-8-tools"


def test_aggregate_counts_and_categories():
    spec = RunSpec("m", "mock", "m")
    results = [
        _result("Detection", True),
        _result("Detection", False),
        _result("Coverage", True),
        _result("Coverage", True, skipped=True),  # skipped excluded
    ]
    agg = aggregate(spec, results)
    assert agg.passed == 2
    assert agg.total == 3
    assert agg.accuracy == pytest.approx(2 / 3)
    assert agg.by_category == {"Detection": (1, 2), "Coverage": (1, 1)}


def test_aggregate_empty_is_zero_not_crash():
    agg = aggregate(RunSpec("m", "mock", "m"), [])
    assert agg.accuracy == 0.0
    assert agg.total == 0


def _report():
    on = SpecResult(
        RunSpec("m (tools)", "mock", "m", tools=True),
        passed=4,
        total=4,
        by_category={"Detection": (2, 2), "Coverage": (2, 2)},
    )
    off = SpecResult(
        RunSpec("m (no tools)", "mock", "m", tools=False),
        passed=1,
        total=4,
        by_category={"Detection": (0, 2), "Coverage": (1, 2)},
    )
    return ComparisonReport([on, off])


def test_ablation_pairs_match_by_model():
    pairs = _report().ablation_pairs()
    assert len(pairs) == 1
    model, on, off = pairs[0]
    assert model == "m"
    assert on.accuracy == 1.0
    assert off.accuracy == 0.25


def test_report_as_dict_has_ablation_lift():
    d = _report().as_dict()
    assert d["ablation"][0]["lift"] == pytest.approx(0.75)
    assert len(d["specs"]) == 2


def test_format_markdown_sections():
    md = format_markdown(_report())
    assert "## Accuracy by configuration" in md
    assert "## Tool-use ablation" in md
    assert "## Accuracy by category" in md
    assert "+75 pts" in md  # the lift
    # Best config is listed first (sorted by accuracy desc).
    cfg_table = md.split("## Tool-use ablation")[0]
    assert cfg_table.index("m (tools)") < cfg_table.index("m (no tools)")


def test_extract_args_parses_prompt():
    args = _extract_args(
        "Find variants in chr1:1000-2000 of tests/fixtures/comprehensive.bam "
        "using reference tests/fixtures/comprehensive_ref.fa. Call get_variants."
    )
    assert args["region"] == "chr1:1000-2000"
    assert args["file_path"] == "tests/fixtures/comprehensive.bam"
    assert args["reference"] == "tests/fixtures/comprehensive_ref.fa"


def test_extract_args_parses_variant_and_defaults_reference():
    args = _extract_args(
        "For variant chr1:1050 A>T in tests/fixtures/comprehensive.bam, "
        "call get_variant_curation_summary with format=rubric."
    )
    assert (args["chrom"], args["pos"], args["ref"], args["alt"]) == ("chr1", 1050, "A", "T")
    assert args["format"] == "rubric"
    # comprehensive.bam implies the comprehensive reference even if unstated.
    assert args["reference"].endswith("comprehensive_ref.fa")


async def test_realistic_mock_calls_tools_only_when_advertised():
    from bamcp.eval.router import InProcessRouter, tool_descriptors

    mock = realistic_mock()
    tools = tool_descriptors(InProcessRouter())
    msgs = [{"role": "user", "content": "coverage for chr1:90-200 in x.bam via get_coverage"}]
    with_tools = await mock.chat_with_tools(msgs, tools)
    assert with_tools.tool_calls  # calls the local tools
    assert all(c.name in {t.name for t in tools} for c in with_tools.tool_calls)

    without = await mock.chat_with_tools(msgs, [])
    assert without.tool_calls == []  # nothing to call
    assert "No tools" in without.text


@pytest.mark.integration
async def test_run_comparison_tools_beat_no_tools(tmp_path, comprehensive_bam_path):
    # End-to-end with the mock provider: the tools arm must beat the no-tools arm.
    cases = load_cases(Path("tests/eval/test_cases.yaml"))
    specs = [
        RunSpec("mock (tools)", "mock", "mock", tools=True),
        RunSpec("mock (no tools)", "mock", "mock", tools=False),
    ]
    report = await run_comparison(specs, cases, tmp_path)
    on = next(s for s in report.specs if s.spec.tools)
    off = next(s for s in report.specs if not s.spec.tools)
    assert on.accuracy > off.accuracy
    assert off.accuracy == 0.0  # no tools -> can't satisfy tools_expected


@pytest.mark.integration
def test_main_cli_writes_tables(tmp_path, comprehensive_bam_path, capsys):
    import json

    from bamcp.eval.compare import main

    out = tmp_path / "cmp"
    rc = main(
        [
            "--provider",
            "mock",
            "--models",
            "mock-model",
            "--output-dir",
            str(out),
            "--subset",
            "1-4",
        ]
    )
    assert rc == 0
    assert "Tool-use ablation" in capsys.readouterr().out
    payload = json.loads((out / "comparison.json").read_text())
    # _default_specs built a tools-on and tools-off spec for the one model.
    assert {s["tools"] for s in payload["specs"]} == {True, False}
    assert payload["ablation"][0]["lift"] > 0
    assert (out / "comparison.md").exists()

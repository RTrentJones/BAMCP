"""Unit tests for report writing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bamcp.eval.report import write_reports, write_summary_line
from bamcp.eval.schema import EvalResult, GraderVerdict, RunConfig


def _make_result(name: str, passed: bool, mode: str | None = None) -> EvalResult:
    return EvalResult(
        case_name=name,
        category="Test",
        rendering_mode=mode,
        response_text="response",
        tool_calls=["get_coverage"],
        duration_ms=12.3,
        verdict=GraderVerdict(
            passed=passed,
            rationale="ok" if passed else "no",
            deterministic=True,
            score=1.0 if passed else 0.0,
        ),
    )


@pytest.mark.unit
def test_write_reports_creates_all_files(tmp_path: Path):
    cfg = RunConfig(
        test_cases_path=tmp_path / "cases.yaml",
        output_dir=tmp_path,
        provider="mock",
        model="m",
    )
    results = [_make_result("c1", True), _make_result("c2", False)]
    paths = write_reports(results, cfg, tmp_path)
    for key in ("json", "yaml", "config", "html"):
        assert paths[key].exists(), f"Missing {key} report"

    parsed = json.loads(paths["json"].read_text())
    assert len(parsed) == 2
    assert parsed[0]["passed"] is True
    assert parsed[1]["passed"] is False


@pytest.mark.unit
def test_summary_line():
    results = [_make_result("c1", True), _make_result("c2", False), _make_result("c3", True)]
    line = write_summary_line(results)
    assert "2/3 passed" in line


@pytest.mark.unit
def test_html_report_escapes_user_strings(tmp_path: Path):
    r = EvalResult(
        case_name="<script>alert(1)</script>",
        category="cat",
        rendering_mode=None,
        response_text="r",
        tool_calls=[],
        duration_ms=0.0,
        verdict=GraderVerdict(passed=False, rationale="bad", deterministic=True),
    )
    cfg = RunConfig(
        test_cases_path=tmp_path / "cases.yaml",
        output_dir=tmp_path,
        provider="mock",
        model="m",
    )
    paths = write_reports([r], cfg, tmp_path)
    html = paths["html"].read_text()
    assert "<script>" not in html
    assert "&lt;script&gt;" in html

"""Unit tests for the eval CLI argument parsing and dry-run path."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from bamcp.eval.cli import _filter_cached, _run, parse_args
from bamcp.eval.schema import EvalCase


@pytest.mark.unit
def test_parse_args_minimal():
    cfg = parse_args(["--test-cases", "x.yaml", "--output-dir", "out"])
    assert str(cfg.test_cases_path) == "x.yaml"
    assert str(cfg.output_dir) == "out"
    assert cfg.provider == "anthropic"
    assert cfg.subset is None


@pytest.mark.unit
def test_parse_args_subset_parsed_correctly():
    cfg = parse_args(["--subset", "2-5"])
    assert cfg.subset == (2, 5)


@pytest.mark.unit
def test_parse_args_invalid_subset_exits(capsys):
    with pytest.raises(SystemExit):
        parse_args(["--subset", "not-valid"])


@pytest.mark.unit
def test_parse_args_dry_run_flag():
    cfg = parse_args(["--dry-run"])
    assert cfg.dry_run is True


@pytest.mark.unit
def test_dry_run_returns_zero_on_valid_file(tmp_path: Path, capsys):
    cases_file = tmp_path / "cases.yaml"
    cases_file.write_text("- case:\n    name: c1\n    input: q\n    expected: a\n")
    cfg = parse_args(["--test-cases", str(cases_file), "--output-dir", str(tmp_path), "--dry-run"])
    code = asyncio.run(_run(cfg))
    assert code == 0
    out = capsys.readouterr().out
    assert "Parsed 1 cases" in out


@pytest.mark.unit
def test_run_returns_nonzero_when_file_missing(tmp_path: Path, capsys):
    cfg = parse_args(
        ["--test-cases", str(tmp_path / "missing.yaml"), "--output-dir", str(tmp_path)]
    )
    code = asyncio.run(_run(cfg))
    assert code == 2


@pytest.mark.unit
def test_run_returns_zero_when_all_cached(tmp_path: Path, capsys):
    """If --cache filters out every case, the runner exits 0 with an info message."""
    cases_file = tmp_path / "cases.yaml"
    cases_file.write_text("- case:\n    name: c1\n    input: q\n    expected: a\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    # Seed a passing prior result so --cache filters everything.
    (out_dir / "results.json").write_text(json.dumps([{"case_name": "c1", "passed": True}]))
    cfg = parse_args(
        [
            "--test-cases",
            str(cases_file),
            "--output-dir",
            str(out_dir),
            "--cache",
            "--provider",
            "mock",
        ]
    )
    code = asyncio.run(_run(cfg))
    assert code == 0
    assert "All cases cached" in capsys.readouterr().out


@pytest.mark.unit
def test_build_judge_returns_none_for_mock():
    """Mock-provider runs should not construct a judge."""
    from bamcp.eval.cli import _build_judge

    cfg = parse_args(["--provider", "mock"])
    assert _build_judge(cfg) is None


@pytest.mark.unit
def test_build_judge_uses_override_model():
    """When --judge-model is provided, the judge uses that model id."""
    from bamcp.eval.cli import _build_judge

    cfg = parse_args(
        ["--provider", "anthropic", "--model", "claude-x", "--judge-model", "claude-y"]
    )
    judge = _build_judge(cfg)
    assert judge is not None
    assert judge.model == "claude-y"


@pytest.mark.unit
def test_build_router_unknown_raises():
    from bamcp.eval.cli import _build_router
    from bamcp.eval.schema import RunConfig

    cfg = RunConfig(
        test_cases_path=Path("x"),
        output_dir=Path("y"),
        router="bogus",
    )
    with pytest.raises(ValueError, match="Unknown router"):
        _build_router(cfg)


@pytest.mark.unit
def test_build_router_mcp_stdio_not_implemented():
    from bamcp.eval.cli import _build_router

    cfg = parse_args(["--router", "mcp-stdio"])
    with pytest.raises(NotImplementedError, match="Phase 2"):
        _build_router(cfg)


@pytest.mark.unit
def test_filter_cached_drops_passing_cases(tmp_path: Path):
    results_path = tmp_path / "results.json"
    results_path.write_text(
        json.dumps([{"case_name": "a", "passed": True}, {"case_name": "b", "passed": False}])
    )
    cases = [
        EvalCase(name="a", input="q", expected="e"),
        EvalCase(name="b", input="q", expected="e"),
        EvalCase(name="c", input="q", expected="e"),
    ]
    out = _filter_cached(results_path, cases)
    names = [c.name for c in out]
    assert "a" not in names  # was passing
    assert "b" in names  # was failing — rerun
    assert "c" in names  # not yet run

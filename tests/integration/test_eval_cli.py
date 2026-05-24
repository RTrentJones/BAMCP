"""Integration test: drive the eval CLI end-to-end with --provider mock."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from bamcp.eval.cli import _run, parse_args


@pytest.mark.integration
def test_eval_cli_runs_subset_with_mock_provider(tmp_path: Path, small_bam_path: str):
    """`make eval` works in mock mode and writes valid result artifacts.

    Uses a single trivial case that fails (mock provider doesn't call any
    tools so the deterministic ``tools_expected`` check fails) — the runner
    must still write its reports cleanly.
    """
    cases_file = tmp_path / "cases.yaml"
    cases_file.write_text(
        "- case:\n"
        "    name: coverage_smoke\n"
        "    category: smoke\n"
        f'    input: "Run get_coverage on {small_bam_path}"\n'
        '    expected: "any"\n'
        '    tools_expected: ["get_coverage"]\n'
    )
    out_dir = tmp_path / "out"

    cfg = parse_args(
        [
            "--test-cases",
            str(cases_file),
            "--output-dir",
            str(out_dir),
            "--provider",
            "mock",
            "--router",
            "in-process",
        ]
    )
    code = asyncio.run(_run(cfg))

    # Mock provider does not call any tool — deterministic check fails — exit nonzero.
    assert code == 1

    for fname in ("results.json", "results.yaml", "report.html", "run_config.yaml"):
        assert (out_dir / fname).exists(), f"Missing artifact: {fname}"

    parsed = json.loads((out_dir / "results.json").read_text())
    assert len(parsed) == 1
    result = parsed[0]
    assert result["case_name"] == "coverage_smoke"
    assert result["passed"] is False
    assert result["verdict_deterministic"] is True


@pytest.mark.integration
def test_eval_cli_dry_run_validates_seed_cases(tmp_path: Path):
    """The shipped tests/eval/test_cases.yaml must pass --dry-run."""
    seed = Path("tests/eval/test_cases.yaml").resolve()
    assert seed.exists()
    cfg = parse_args(
        [
            "--test-cases",
            str(seed),
            "--output-dir",
            str(tmp_path),
            "--dry-run",
        ]
    )
    code = asyncio.run(_run(cfg))
    assert code == 0

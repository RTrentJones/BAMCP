"""Eval report writers (YAML, JSON, HTML).

Matches MARRVEL-MCP's output layout: each run drops a ``results.yaml``,
``results.json``, ``run_config.yaml``, ``test_cases.yaml``, and a single-file
``report.html`` under the chosen output directory.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from .schema import EvalResult, RunConfig


def write_reports(
    results: list[EvalResult],
    config: RunConfig,
    output_dir: Path,
) -> dict[str, Path]:
    """Write all report artifacts. Returns a dict of name -> file path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    serialized = [_serialize_result(r) for r in results]

    paths: dict[str, Path] = {}

    json_path = output_dir / "results.json"
    json_path.write_text(json.dumps(serialized, indent=2), encoding="utf-8")
    paths["json"] = json_path

    yaml_path = output_dir / "results.yaml"
    yaml_path.write_text(_dump_yaml(serialized), encoding="utf-8")
    paths["yaml"] = yaml_path

    cfg_path = output_dir / "run_config.yaml"
    cfg_path.write_text(_dump_yaml(_config_to_dict(config)), encoding="utf-8")
    paths["config"] = cfg_path

    html_path = output_dir / "report.html"
    html_path.write_text(_render_html(results, config), encoding="utf-8")
    paths["html"] = html_path

    return paths


def _serialize_result(r: EvalResult) -> dict:
    d = asdict(r)
    # Promote verdict fields one level up for easier scanning.
    verdict = d.pop("verdict")
    d["passed"] = verdict["passed"]
    d["verdict_rationale"] = verdict["rationale"]
    d["verdict_deterministic"] = verdict["deterministic"]
    d["verdict_score"] = verdict["score"]
    return d


def _config_to_dict(config: RunConfig) -> dict:
    return {
        "provider": config.provider,
        "model": config.model,
        "judge_model": config.judge_model,
        "cache": config.cache,
        "subset": list(config.subset) if config.subset else None,
        "with_vanilla": config.with_vanilla,
        "with_rendering_comparison": config.with_rendering_comparison,
        "router": config.router,
        "test_cases_path": str(config.test_cases_path),
        "output_dir": str(config.output_dir),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def _dump_yaml(obj: object) -> str:
    try:
        import yaml  # type: ignore[import-untyped]

        return str(yaml.safe_dump(obj, sort_keys=False, default_flow_style=False))
    except ImportError:
        # Hard fallback so the report is still written without pyyaml.
        return json.dumps(obj, indent=2)


def _render_html(results: list[EvalResult], config: RunConfig) -> str:
    rows: list[str] = []
    pass_count = sum(1 for r in results if r.verdict.passed)
    total = len(results)
    rate = (pass_count / total * 100.0) if total else 0.0

    for r in results:
        verdict_class = "pass" if r.verdict.passed else "fail"
        mode = r.rendering_mode or "—"
        tool_chips = " ".join(f"<span class='chip'>{t}</span>" for t in r.tool_calls)
        rows.append(
            f"<tr class='{verdict_class}'>"
            f"<td>{_h(r.case_name)}</td>"
            f"<td>{_h(r.category)}</td>"
            f"<td>{_h(mode)}</td>"
            f"<td>{'PASS' if r.verdict.passed else 'FAIL'}</td>"
            f"<td>{'det' if r.verdict.deterministic else 'judge'}</td>"
            f"<td>{r.duration_ms:.1f} ms</td>"
            f"<td>{tool_chips}</td>"
            f"<td>{_h(r.verdict.rationale)}</td>"
            f"</tr>"
        )

    return (
        "<!doctype html><meta charset='utf-8'>"
        "<title>BAMCP eval report</title>"
        "<style>"
        "body { font-family: ui-sans-serif, system-ui, sans-serif; margin: 24px; }"
        "h1 { margin-bottom: 4px; }"
        ".summary { color: #475569; margin-bottom: 16px; }"
        "table { border-collapse: collapse; width: 100%; font-size: 13px; }"
        "th, td { border: 1px solid #e2e8f0; padding: 6px 8px;"
        " text-align: left; vertical-align: top; }"
        "th { background: #f8fafc; }"
        "tr.pass { background: #f0fdf4; }"
        "tr.fail { background: #fef2f2; }"
        ".chip { display: inline-block; background: #e0e7ff; color: #3730a3; padding: 1px 6px; "
        "border-radius: 999px; margin-right: 4px; font-size: 11px; }"
        "</style>"
        "<h1>BAMCP eval report</h1>"
        f"<div class='summary'>{pass_count}/{total} passed ({rate:.1f}%) — "
        f"provider={_h(config.provider)} model={_h(config.model)} router={_h(config.router)}</div>"
        "<table><thead><tr>"
        "<th>Case</th><th>Category</th><th>Mode</th><th>Verdict</th>"
        "<th>Grader</th><th>Duration</th><th>Tools called</th><th>Rationale</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    )


def _h(s: object) -> str:
    """Minimal HTML escaping."""
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_summary_line(results: list[EvalResult]) -> str:
    """One-line summary printed to stdout at the end of a run."""
    pass_count = sum(1 for r in results if r.verdict.passed)
    total = len(results)
    rate = (pass_count / total * 100.0) if total else 0.0
    return f"Eval: {pass_count}/{total} passed ({rate:.1f}%)"

"""Publish eval results to Tracer, the RTrentJones project-signals dashboard.

Tracer (https://tracer.rtrentjones.dev) stores one *run* per eval invocation plus its *cases*, and
trends pass rate, numeric metrics, regressions and flaky cases per suite. This module converts
BAMCP's existing eval outputs into Tracer's native ingest payload and, optionally, POSTs them:

* ``truthset`` — the deterministic ground-truth gate (``bamcp.eval.truthset --json``). One run per
  manifest (suite = dataset). Cases are the gate's floors (exactly the checks ``meets_floors``
  applies, at the thresholds the report records), each artifact-type check, and every missing or
  spurious call, so the run page shows *what* broke rather than only that it did.
* ``llm`` — an LLM harness run directory (``results.json`` + ``run_config.yaml``), or every spec
  directory of a model comparison (``comparison.json`` alongside ``<slug>/results.json``), so
  Tracer's compare view sets the models (and the tools on/off arms) side by side.

Publishing is best-effort by design: a missing ``TRACER_INGEST_TOKEN`` skips with a notice and a
failed POST prints a warning, but neither changes the exit code — the eval's own gate decides
pass/fail, not the dashboard. Provenance (commit, repo, branch, CI run URL, trigger) is read from
the standard GitHub Actions environment when present.

Usage::

    python -m bamcp.eval.tracer truthset truthset-synthetic.json truthset-indel.json --post
    python -m bamcp.eval.tracer llm .eval-results --post
    python -m bamcp.eval.tracer llm .eval-results/comparison --post
    python -m bamcp.eval.tracer truthset report.json --out payload.json   # inspect, don't send
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

TOOL = "bamcp"
DEFAULT_URL = "https://tracer.rtrentjones.dev"

# Tracer's ingest bounds (tools/tracer/lib/schema.ts): text fields ≤16 KB, ≤1000 cases, ≤64
# metric keys. Trimming here keeps a large report a 201 rather than a 422.
TEXT_MAX = 16_000
MAX_CASES = 1000
MAX_METRICS = 64
# Tracer stores case names in a varchar(255); a long spurious indel allele would otherwise 422 on
# an older Tracer (newer ones clip server-side). The full call stays in the case's expected/output.
NAME_MAX = 255

# The floors meets_floors() applies without a CLI override. The truthset CLI only exposes
# min_recall / min_precision (and records them in the report's "floors"); everything else is the
# method default. Kept in sync with TruthsetReport.meets_floors by test_tracer.py.
DEFAULT_FLOORS: dict[str, float] = {
    "min_recall": 1.0,
    "min_precision": 0.95,
    "min_artifact_recall": 1.0,
    "max_false_positives": 0,
    "max_overconfident": 0,
    "max_confidence_mismatches": 0,
    "max_tool_failures": 0,
    "max_likelihood_mismatches": 0,
    "max_vaf_mismatches": 0,
}


def _clip(text: Any) -> str | None:
    if text is None:
        return None
    s = text if isinstance(text, str) else json.dumps(text, default=str)
    return s if len(s) <= TEXT_MAX else s[:TEXT_MAX] + "… [truncated]"


def _case_name(s: str) -> str:
    return s if len(s) <= NAME_MAX else s[: NAME_MAX - 1] + "…"


def _clamp01(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if v != v:  # NaN
        return None
    return max(0.0, min(1.0, v))


def _finite_metrics(m: Mapping[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for k, v in m.items():
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            continue
        fv = float(v)
        if fv != fv or fv in (float("inf"), float("-inf")):
            continue
        # Tracer metric keys are identifiers (`[\w.:/-]`, ≤80 chars); categories like
        # "Region Analysis" would otherwise 422 the whole payload.
        # ASCII-only: Tracer's `\w` is ASCII, so "präzision" must not pass through as-is.
        key = re.sub(r"[^\w.:/-]+", "_", str(k), flags=re.ASCII)[:80]
        out[key] = int(fv) if isinstance(v, int) else round(fv, 6)
        if len(out) >= MAX_METRICS:
            break
    return out


def github_provenance(env: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Commit / repo / branch / CI-run URL from the GitHub Actions environment (empty locally)."""
    e = os.environ if env is None else env
    out: dict[str, Any] = {}
    if e.get("GITHUB_SHA"):
        out["git_sha"] = e["GITHUB_SHA"]
    if e.get("GITHUB_REPOSITORY"):
        out["repo"] = e["GITHUB_REPOSITORY"]
    if e.get("GITHUB_REF_NAME"):
        out["branch"] = e["GITHUB_REF_NAME"]
    if e.get("GITHUB_RUN_ID") and e.get("GITHUB_REPOSITORY"):
        server = e.get("GITHUB_SERVER_URL", "https://github.com")
        out["source_url"] = f"{server}/{e['GITHUB_REPOSITORY']}/actions/runs/{e['GITHUB_RUN_ID']}"
    labels: dict[str, str] = {}
    if e.get("GITHUB_EVENT_NAME"):
        labels["trigger"] = e["GITHUB_EVENT_NAME"]
    if e.get("GITHUB_WORKFLOW"):
        labels["workflow"] = e["GITHUB_WORKFLOW"][:200]
    if labels:
        out["labels"] = labels
    return out


def _with_provenance(payload: dict[str, Any], provenance: Mapping[str, Any]) -> dict[str, Any]:
    for key in ("git_sha", "repo", "branch", "source_url"):
        if provenance.get(key) and not payload.get(key):
            payload[key] = provenance[key]
    labels = {**provenance.get("labels", {}), **payload.get("labels", {})}
    if labels:
        payload["labels"] = labels
    return payload


# ── truthset ───────────────────────────────────────────────────────────────


def _floor_cases(report: Mapping[str, Any], floors: Mapping[str, float]) -> list[dict[str, Any]]:
    """One case per check in TruthsetReport.meets_floors, in the same order and at the same
    thresholds, so the case verdicts agree with the report's own ``passed``."""
    det = report.get("detection", {})
    tp, fn = int(det.get("tp", 0)), int(det.get("fn", 0))
    recall, precision = float(det.get("recall", 0.0)), float(det.get("precision", 0.0))
    art = float(report.get("artifact_recall", 0.0))

    def count(key: str) -> int:
        return len(report.get(key) or [])

    def le(value: int, limit: float) -> bool:
        return value <= limit

    rows: list[tuple[str, bool, float | None, str, str]] = [
        ("non-empty truth set", tp + fn > 0, None, f"{tp + fn} expected alleles", "> 0"),
        (
            "tool calls succeed",
            le(count("tool_failures"), floors["max_tool_failures"]),
            None,
            str(report.get("tool_failures") or "none"),
            f"≤ {floors['max_tool_failures']:g} failures",
        ),
        (
            "variant recall",
            recall >= floors["min_recall"],
            recall,
            f"{recall:.4f}",
            f"≥ {floors['min_recall']:g}",
        ),
        (
            "variant precision",
            precision >= floors["min_precision"],
            precision,
            f"{precision:.4f}",
            f"≥ {floors['min_precision']:g}",
        ),
        (
            "no false positives in negative regions",
            le(int(report.get("false_positives", 0)), floors["max_false_positives"]),
            None,
            str(report.get("false_positives", 0)),
            f"≤ {floors['max_false_positives']:g}",
        ),
        (
            "artifact-type recall",
            art >= floors["min_artifact_recall"],
            art,
            f"{art:.4f}",
            f"≥ {floors['min_artifact_recall']:g}",
        ),
        (
            "clean control discriminated",
            bool(report.get("clean_discriminated", False)),
            None,
            str(bool(report.get("clean_discriminated", False))).lower(),
            "true",
        ),
        (
            "SAFETY: no overconfident artifact sites",
            le(count("overconfident_sites"), floors["max_overconfident"]),
            None,
            str(report.get("overconfident_sites") or "none"),
            f"≤ {floors['max_overconfident']:g}",
        ),
        (
            "confidence calibration",
            le(count("confidence_mismatches"), floors["max_confidence_mismatches"]),
            None,
            str(report.get("confidence_mismatches") or "none"),
            f"≤ {floors['max_confidence_mismatches']:g} mismatches",
        ),
        (
            "artifact likelihood",
            le(count("likelihood_mismatches"), floors["max_likelihood_mismatches"]),
            None,
            str(report.get("likelihood_mismatches") or "none"),
            f"≤ {floors['max_likelihood_mismatches']:g} mismatches",
        ),
        (
            "VAF within tolerance",
            le(count("vaf_mismatches"), floors["max_vaf_mismatches"]),
            None,
            str(report.get("vaf_mismatches") or "none"),
            f"≤ {floors['max_vaf_mismatches']:g} mismatches",
        ),
    ]
    return [
        {
            "name": _case_name(f"floor: {name}"),
            "category": "floors",
            "passed": ok,
            "score": score if score is not None else (1.0 if ok else 0.0),
            "output": _clip(observed),
            "expected": expected,
        }
        for name, ok, score, observed, expected in rows
    ]


def truthset_payload(
    report: Mapping[str, Any],
    *,
    env: str = "ci",
    provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Map one ``bamcp.eval.truthset --json`` report to a Tracer run."""
    floors = {**DEFAULT_FLOORS, **{k: float(v) for k, v in (report.get("floors") or {}).items()}}
    cases = _floor_cases(report, floors)
    floor_passed = sum(1 for c in cases if c["passed"])
    n_floors = len(cases)

    for chk in report.get("artifact_checks") or []:
        found = bool(chk.get("found"))
        cases.append(
            {
                "name": _case_name(f"artifact: {chk.get('site')} ({chk.get('expected_type')})"),
                "category": "artifact checks",
                "passed": found,
                "score": 1.0 if found else 0.0,
                "expected": _clip(chk.get("expected_type")),
                "output": (
                    f"found={found} risk_score={chk.get('risk_score')} "
                    f"likelihood={chk.get('likelihood')}"
                ),
                "metrics": _finite_metrics({"risk_score": chk.get("risk_score")}),
            }
        )
    for call in report.get("missing_calls") or []:
        cases.append(
            {
                "name": _case_name(f"missing call: {call}"),
                "category": "calls",
                "passed": False,
                "score": 0.0,
                "expected": _clip(call),
                "judge_rationale": "expected variant allele not reported (false negative)",
            }
        )
    for call in report.get("spurious_calls") or []:
        cases.append(
            {
                "name": _case_name(f"spurious call: {call}"),
                "category": "calls",
                "passed": False,
                "score": 0.0,
                "output": _clip(call),
                "judge_rationale": "reported allele absent from the truth set (false positive)",
            }
        )

    det = report.get("detection", {})
    metrics = _finite_metrics(
        {
            "precision": det.get("precision"),
            "recall": det.get("recall"),
            "f1": det.get("f1"),
            "tp": det.get("tp"),
            "fp": det.get("fp"),
            "fn": det.get("fn"),
            "artifact_recall": report.get("artifact_recall"),
            "false_positives_negative_regions": report.get("false_positives"),
            "overconfident_sites": len(report.get("overconfident_sites") or []),
            "floors_met": floor_passed,
        }
    )
    labels: dict[str, str] = {}
    if report.get("version"):
        labels["dataset_version"] = str(report["version"])[:200]
    failures = report.get("failures") or []
    if failures:
        labels["failures"] = "; ".join(map(str, failures))[:255]
    payload: dict[str, Any] = {
        "tool": TOOL,
        # The system under test is BAMCP's deterministic tool pipeline; no model is involved.
        "model": "deterministic",
        "mode": "truthset",
        # The deterministic gate is what CI enforces, so it gates project health in Tracer.
        "gating": True,
        "env": env,
        "suite": str(report.get("dataset") or "unknown"),
        # Pass rate = floors met. The gate is the floors; artifact/call cases are evidence.
        "pass_rate": round(floor_passed / n_floors, 4) if n_floors else 0.0,
        "passed": bool(report.get("passed", floor_passed == n_floors)),
        "metrics": metrics,
        "labels": labels,
        "cases": cases[:MAX_CASES],
    }
    return _with_provenance(payload, provenance or {})


# ── LLM harness / comparison ───────────────────────────────────────────────


def _load_yaml(path: Path) -> Any:
    import yaml  # type: ignore[import-untyped]

    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _test_cases(path: Path | None) -> dict[str, dict[str, Any]]:
    """name → raw test case (input / expected / expected_claims), for the case input/expected."""
    if path is None or not path.is_file():
        return {}
    raw = _load_yaml(path)
    items = raw.get("cases", raw) if isinstance(raw, dict) else raw
    out: dict[str, dict[str, Any]] = {}
    for item in items or []:
        d = item.get("case", item) if isinstance(item, dict) else None
        if isinstance(d, dict) and "name" in d:
            out[str(d["name"])] = d
    return out


def llm_payload(
    run_dir: Path,
    *,
    env: str = "ci",
    model_label: str | None = None,
    suite: str | None = None,
    test_cases_path: Path | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Map one LLM harness output directory (results.json + run_config.yaml) to a Tracer run."""
    rows: list[dict[str, Any]] = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    cfg_path = run_dir / "run_config.yaml"
    cfg: dict[str, Any] = _load_yaml(cfg_path) if cfg_path.is_file() else {}
    cfg = cfg or {}
    tc_path = test_cases_path or (
        Path(cfg["test_cases_path"]) if cfg.get("test_cases_path") else None
    )
    known = _test_cases(tc_path)

    graded = [r for r in rows if not r.get("skipped")]
    cases: list[dict[str, Any]] = []
    by_cat: dict[str, list[bool]] = {}
    durations: list[float] = []
    for r in graded:
        name = str(r.get("case_name", "unknown"))
        tc = known.get(name, {})
        claims = tc.get("expected_claims") or []
        expected = tc.get("expected")
        if claims:
            expected = f"{expected or ''}\n\nrequired claims: {', '.join(map(str, claims))}".strip()
        passed = bool(r.get("passed"))
        cat = str(r.get("category") or "uncategorized")
        by_cat.setdefault(cat, []).append(passed)
        if isinstance(r.get("duration_ms"), (int, float)):
            durations.append(float(r["duration_ms"]))
        rationale = r.get("verdict_rationale") or ""
        if r.get("error"):
            rationale = f"{rationale}\n\nerror: {r['error']}".strip()
        sub = _finite_metrics(
            {
                "tool_score": r.get("verdict_tool_score"),
                "answer_score": r.get("verdict_answer_score"),
                "duration_ms": r.get("duration_ms"),
                "tool_calls": len(r.get("tool_calls") or []),
            }
        )
        cases.append(
            {
                "name": _case_name(name),
                "category": cat,
                "input": _clip(tc.get("input")),
                "expected": _clip(expected),
                "output": _clip(r.get("response_text")),
                "score": _clamp01(r.get("verdict_score")),
                "passed": passed,
                "judge_rationale": _clip(rationale) or None,
                "metrics": sub,
            }
        )

    n = len(cases)
    n_pass = sum(1 for c in cases if c["passed"])
    first = rows[0] if rows else {}
    model = model_label or str(cfg.get("model") or first.get("model") or "unknown")
    metrics: dict[str, Any] = {"cases": n, "skipped": len(rows) - n}
    for cat, verdicts in sorted(by_cat.items()):
        metrics[f"category.{cat}"] = round(sum(verdicts) / len(verdicts), 4)
    if durations:
        metrics["mean_case_ms"] = round(sum(durations) / len(durations), 1)
    labels = {
        k: str(cfg[k])[:255]
        for k in ("provider", "model", "judge_provider", "judge_model", "router")
        if cfg.get(k) is not None
    }
    payload: dict[str, Any] = {
        "tool": TOOL,
        "model": model,
        "mode": "llm-eval",
        # Trended, not gated: LLM runs (and the comparison arms built from them) are evidence of
        # model behaviour, not a ship-gate, so a flicker never marks BAMCP failing in Tracer.
        "gating": False,
        "env": env,
        "suite": suite or (tc_path.stem if tc_path else "test_cases"),
        "duration_ms": int(sum(durations)) if durations else None,
        "passed": n > 0 and n_pass == n,
        "pass_rate": round(n_pass / n, 4) if n else 0.0,
        "metrics": _finite_metrics(metrics),
        "labels": labels,
        "cases": cases[:MAX_CASES],
    }
    return _with_provenance(payload, provenance or {})


def comparison_payloads(
    comparison_dir: Path,
    *,
    env: str = "ci",
    test_cases_path: Path | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """One run per comparison spec (model × tools/vision arm), labelled with the spec label so the
    arms line up as separate "models" in Tracer's compare view. All share suite ``comparison``."""
    report = json.loads((comparison_dir / "comparison.json").read_text(encoding="utf-8"))
    payloads: list[dict[str, Any]] = []
    for spec in report.get("specs", []):
        label = str(spec.get("label", ""))
        slug = re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")  # RunSpec.slug
        spec_dir = comparison_dir / slug
        if not (spec_dir / "results.json").is_file():
            continue
        p = llm_payload(
            spec_dir,
            env=env,
            model_label=label,
            suite="comparison",
            test_cases_path=test_cases_path,
            provenance=provenance,
        )
        p["labels"] = {
            **p.get("labels", {}),
            "tools": str(spec.get("tools")),
            "vision": str(spec.get("vision")),
        }
        payloads.append(p)
    return payloads


def llm_dir_payloads(
    path: Path,
    *,
    env: str = "ci",
    test_cases_path: Path | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """A harness run dir → one payload; a comparison dir (has comparison.json) → one per spec."""
    if (path / "comparison.json").is_file():
        return comparison_payloads(
            path, env=env, test_cases_path=test_cases_path, provenance=provenance
        )
    if (path / "results.json").is_file():
        return [llm_payload(path, env=env, test_cases_path=test_cases_path, provenance=provenance)]
    return []


# ── transport ──────────────────────────────────────────────────────────────


def post(
    payload: Mapping[str, Any], *, url: str, token: str, timeout: float = 30.0
) -> tuple[bool, str]:
    """POST one payload to ``<url>/api/ingest``. Returns (ok, message); never raises."""
    # The URL is caller-supplied (--url / TRACER_URL): refuse anything but http(s) so urlopen can
    # never be pointed at file:// or a custom scheme.
    if not url.lower().startswith(("https://", "http://")):
        return False, f"refusing non-http(s) Tracer URL: {url!r}"
    req = urllib.request.Request(  # noqa: S310 - scheme checked above
        f"{url.rstrip('/')}/api/ingest",
        data=json.dumps(payload).encode("utf-8"),
        method="POST",
        headers={"authorization": f"Bearer {token}", "content-type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 - scheme checked
            return True, f"{resp.status} {resp.read().decode('utf-8', 'replace')[:200]}"
    except urllib.error.HTTPError as e:
        return False, f"HTTP {e.code}: {e.read().decode('utf-8', 'replace')[:200]}"
    except (urllib.error.URLError, OSError) as e:
        return False, f"unreachable: {e}"


def _step_summary(line: str) -> None:
    """Append one Markdown line to the Actions job summary, when running in Actions."""
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not path:
        return
    try:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(" ".join(line.split()) + "\n")  # one line, even for a multi-line error
    except OSError:
        pass  # the summary is a convenience; publishing must never fail on it


def publish(payloads: Sequence[Mapping[str, Any]], *, url: str, token: str | None) -> int:
    """POST each payload; returns the number that landed. Skips (0) without a token.

    Each outcome is also written to ``$GITHUB_STEP_SUMMARY`` so a skipped or failed publish is
    visible on the run page, not only in the log."""
    if not token:
        print("tracer: TRACER_INGEST_TOKEN unset — skipping publish.")
        _step_summary("⚠️ Tracer publish skipped: TRACER_INGEST_TOKEN unset")
        return 0
    ok_count = 0
    for p in payloads:
        ok, msg = post(p, url=url, token=token)
        label = f"{p.get('mode')}/{p.get('suite')}/{p.get('model')}"
        if ok:
            ok_count += 1
            print(f"tracer: published {label} → {msg}")
            _step_summary(f"✅ Tracer: published {label} ({msg.split(' ', 1)[0]})")
        else:
            # A GitHub annotation, not a failure: the eval's own gate owns the exit code.
            print(f"::warning title=Tracer publish failed::{label}: {msg}")
            _step_summary(f"⚠️ Tracer publish failed: {label}: {msg}")
    return ok_count


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="bamcp-tracer", description="Publish BAMCP eval results to Tracer."
    )
    sub = p.add_subparsers(dest="kind", required=True)
    t = sub.add_parser("truthset", help="truthset --json report file(s)")
    t.add_argument("reports", nargs="+", type=Path)
    llm = sub.add_parser("llm", help="LLM harness output dir, or a comparison dir")
    llm.add_argument("dirs", nargs="+", type=Path)
    llm.add_argument(
        "--test-cases", type=Path, default=None, help="test_cases.yaml (input/expected)"
    )
    for sp in (t, llm):
        sp.add_argument(
            "--env", default=os.environ.get("TRACER_ENV", "ci"), help="Tracer env label"
        )
        sp.add_argument(
            "--post", action="store_true", help="POST to Tracer (needs TRACER_INGEST_TOKEN)"
        )
        sp.add_argument("--out", type=Path, default=None, help="Write the payload(s) as JSON here")
        sp.add_argument("--url", default=os.environ.get("TRACER_URL", DEFAULT_URL))
    args = p.parse_args(argv)

    prov = github_provenance()
    payloads: list[dict[str, Any]] = []
    if args.kind == "truthset":
        for path in args.reports:
            try:
                report = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as e:
                print(f"::warning title=Tracer publish skipped::{path}: unreadable report ({e})")
                continue
            payloads.append(truthset_payload(report, env=args.env, provenance=prov))
    else:
        for d in args.dirs:
            found = llm_dir_payloads(
                d, env=args.env, test_cases_path=args.test_cases, provenance=prov
            )
            if not found:
                print(f"::warning title=Tracer publish skipped::{d}: no results.json or comparison")
            payloads.extend(found)

    if args.out:
        args.out.write_text(json.dumps(payloads, indent=2), encoding="utf-8")
        print(f"tracer: wrote {len(payloads)} payload(s) to {args.out}")
    if args.post:
        publish(payloads, url=args.url, token=os.environ.get("TRACER_INGEST_TOKEN"))
    elif not args.out:
        print(json.dumps(payloads, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

"""Tests for the Tracer exporter (bamcp.eval.tracer).

The payload builders are pure, so they are tested on hand-built reports and result directories. The
key invariant is that the per-floor cases agree with ``TruthsetReport.meets_floors`` — the dashboard
must never show a green floor the gate failed, or vice versa. Transport is tested with a stubbed
``urlopen``; nothing here touches the network.
"""

from __future__ import annotations

import inspect
import io
import json
import urllib.error
from pathlib import Path
from typing import Any

import pytest

from bamcp.eval import tracer
from bamcp.eval.metrics import PRF, ArtifactCheck
from bamcp.eval.truthset import TruthsetReport

pytestmark = pytest.mark.unit


def _report(**over: Any) -> TruthsetReport:
    base: dict[str, Any] = {
        "dataset": "synthetic_v1",
        "version": 1,
        "detection": PRF(tp=9, fp=0, fn=0),
        "false_positives": 0,
        "artifact_checks": [
            ArtifactCheck(
                site="chr1:1051",
                expected_type="strand_bias",
                found=True,
                risk_score=0.7,
                likelihood="high",
            )
        ],
    }
    base.update(over)
    return TruthsetReport(**base)


def _json(report: TruthsetReport, **floors: float) -> dict[str, Any]:
    """What `python -m bamcp.eval.truthset --json` prints, at the given CLI floors."""
    applied = {"min_recall": 1.0, "min_precision": 0.95, **floors}
    passed, failures = report.meets_floors(**applied)
    return {**report.as_dict(), "passed": passed, "failures": failures, "floors": applied}


def test_default_floors_match_meets_floors_signature():
    sig = inspect.signature(TruthsetReport.meets_floors)
    defaults = {k: p.default for k, p in sig.parameters.items() if k != "self"}
    assert defaults == tracer.DEFAULT_FLOORS


@pytest.mark.parametrize(
    ("over", "floors"),
    [
        ({}, {}),
        ({"detection": PRF(tp=8, fp=0, fn=1), "missing_calls": ["chr1:2801:G>A"]}, {}),
        (
            {"detection": PRF(tp=9, fp=1, fn=0), "spurious_calls": ["chr1:50:A>C"]},
            {"min_precision": 1.0},
        ),
        ({"overconfident_sites": ["chr1:1051"]}, {}),
        ({"confidence_mismatches": ["chr1:2801"], "vaf_mismatches": ["chr1:900"]}, {}),
        ({"clean_discriminated": False, "false_positives": 2}, {}),
        ({"detection": PRF(tp=0, fp=0, fn=0)}, {}),
        ({"tool_failures": ["get_variants: boom"], "likelihood_mismatches": ["chr1:1"]}, {}),
        (
            {
                "artifact_checks": [
                    ArtifactCheck(
                        site="chr1:1051",
                        expected_type="strand_bias",
                        found=False,
                        risk_score=0.1,
                        likelihood="low",
                    )
                ]
            },
            {},
        ),
    ],
)
def test_floor_cases_agree_with_meets_floors(over: dict[str, Any], floors: dict[str, float]):
    report = _report(**over)
    data = _json(report, **floors)
    payload = tracer.truthset_payload(data)
    floor_cases = [c for c in payload["cases"] if c["category"] == "floors"]
    assert payload["passed"] is data["passed"]
    # Every floor passes iff the gate passes: a floor case can't contradict the report.
    assert all(c["passed"] for c in floor_cases) is data["passed"]
    assert sum(not c["passed"] for c in floor_cases) == len(data["failures"])


def test_truthset_payload_shape_and_metrics():
    report = _report(
        detection=PRF(tp=8, fp=1, fn=1),
        missing_calls=["chr1:2801:G>A"],
        spurious_calls=["chr1:50:A>C"],
    )
    payload = tracer.truthset_payload(
        _json(report), env="nightly", provenance={"git_sha": "abc1234", "repo": "o/r"}
    )
    assert payload["tool"] == "bamcp"
    assert payload["mode"] == "truthset"
    assert payload["model"] == "deterministic"
    assert payload["suite"] == "synthetic_v1"
    assert payload["env"] == "nightly"
    assert payload["git_sha"] == "abc1234"
    assert payload["metrics"]["tp"] == 8
    assert payload["metrics"]["recall"] == pytest.approx(8 / 9, abs=1e-4)
    assert 0 <= payload["pass_rate"] < 1
    names = {c["name"] for c in payload["cases"]}
    assert "missing call: chr1:2801:G>A" in names
    assert "spurious call: chr1:50:A>C" in names
    assert "artifact: chr1:1051 (strand_bias)" in names
    assert payload["labels"]["dataset_version"] == "1"
    assert "variant recall" in payload["labels"]["failures"]


def test_truthset_payload_is_gating():
    assert tracer.truthset_payload(_json(_report()))["gating"] is True


def test_long_case_names_are_clipped_and_the_call_kept_whole():
    # A 300 bp spurious insertion: the name would overflow Tracer's varchar(255) case name.
    call = "chr1:50:A>A" + "C" * 300
    report = _report(detection=PRF(tp=9, fp=1, fn=0), spurious_calls=[call])
    payload = tracer.truthset_payload(_json(report, min_precision=0.5))
    case = next(c for c in payload["cases"] if c["name"].startswith("spurious call: "))
    assert len(case["name"]) == 255
    assert case["name"].endswith("…")
    assert case["output"] == call
    assert all(len(c["name"]) <= 255 for c in payload["cases"])
    assert tracer._case_name("short") == "short"
    assert tracer._case_name("x" * 255) == "x" * 255


def test_truthset_payload_tolerates_a_report_without_floors():
    data = _json(_report())
    data.pop("floors")
    assert tracer.truthset_payload(data)["passed"] is True


def _write_run(d: Path, rows: list[dict[str, Any]], cfg: dict[str, Any] | None = None) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / "results.json").write_text(json.dumps(rows), encoding="utf-8")
    if cfg is not None:
        import yaml  # type: ignore[import-untyped]

        (d / "run_config.yaml").write_text(yaml.safe_dump(cfg), encoding="utf-8")


def _row(name: str, passed: bool, **over: Any) -> dict[str, Any]:
    return {
        "case_name": name,
        "category": "lookup",
        "response_text": f"answer for {name}",
        "tool_calls": ["get_variants"],
        "duration_ms": 1200.0,
        "passed": passed,
        "verdict_rationale": "claims present" if passed else "missing claim",
        "verdict_score": 1.0 if passed else 0.0,
        "verdict_tool_score": 1.0,
        "verdict_answer_score": 1.0 if passed else 0.0,
        "skipped": False,
        "provider": "anthropic",
        "model": "claude-sonnet-5-5",
        **over,
    }


def test_llm_payload_joins_test_cases_and_skips_skipped(tmp_path: Path):
    tc = tmp_path / "cases.yaml"
    tc.write_text(
        "cases:\n"
        "  - name: a\n    input: what is at chr1:100?\n    expected: a SNV\n"
        "    expected_claims: [A>T]\n"
        "  - name: b\n    input: coverage?\n    expected: 30x\n",
        encoding="utf-8",
    )
    run = tmp_path / "run"
    _write_run(
        run,
        [_row("a", True), _row("b", False, error="timeout"), _row("v", True, skipped=True)],
        {
            "provider": "anthropic",
            "model": "claude-sonnet-5-5",
            "judge_model": "claude-haiku-5-5",
            "test_cases_path": str(tc),
        },
    )
    p = tracer.llm_payload(run, provenance={"labels": {"trigger": "schedule"}})
    assert p["model"] == "claude-sonnet-5-5"
    assert p["mode"] == "llm-eval"
    assert p["suite"] == "cases"
    assert p["pass_rate"] == 0.5
    assert p["passed"] is False
    assert [c["name"] for c in p["cases"]] == ["a", "b"]
    a = p["cases"][0]
    assert a["input"] == "what is at chr1:100?"
    assert "required claims: A>T" in a["expected"]
    assert "error: timeout" in p["cases"][1]["judge_rationale"]
    assert p["metrics"]["skipped"] == 1
    assert p["metrics"]["category.lookup"] == 0.5
    assert p["labels"]["judge_model"] == "claude-haiku-5-5"
    assert p["labels"]["trigger"] == "schedule"
    # LLM runs are trended, not gated.
    assert p["gating"] is False


def test_llm_payload_clips_a_long_case_name(tmp_path: Path):
    _write_run(tmp_path / "run", [_row("q" * 400, True)])
    p = tracer.llm_payload(tmp_path / "run")
    assert len(p["cases"][0]["name"]) == 255


def test_comparison_payloads_one_run_per_spec(tmp_path: Path):
    comp = tmp_path / "comparison"
    specs = [
        {
            "label": "Sonnet 5.5 + tools",
            "provider": "anthropic",
            "model": "claude-sonnet-5-5",
            "tools": True,
            "vision": "off",
        },
        {
            "label": "Sonnet 5.5 (no tools)",
            "provider": "anthropic",
            "model": "claude-sonnet-5-5",
            "tools": False,
            "vision": "off",
        },
        {
            "label": "missing dir",
            "provider": "anthropic",
            "model": "x",
            "tools": True,
            "vision": "off",
        },
    ]
    comp.mkdir()
    (comp / "comparison.json").write_text(
        json.dumps({"specs": specs, "ablation": []}), encoding="utf-8"
    )
    _write_run(comp / "sonnet-5-5-tools", [_row("a", True)])
    _write_run(comp / "sonnet-5-5-no-tools", [_row("a", False)])
    payloads = tracer.llm_dir_payloads(comp)
    assert [p["model"] for p in payloads] == ["Sonnet 5.5 + tools", "Sonnet 5.5 (no tools)"]
    assert {p["suite"] for p in payloads} == {"comparison"}
    assert payloads[1]["labels"]["tools"] == "False"
    assert all(p["gating"] is False for p in payloads)
    assert tracer.llm_dir_payloads(tmp_path / "nothing-here") == []


def test_github_provenance():
    env = {
        "GITHUB_SHA": "a" * 40,
        "GITHUB_REPOSITORY": "RTrentJones/BAMCP",
        "GITHUB_REF_NAME": "main",
        "GITHUB_RUN_ID": "42",
        "GITHUB_SERVER_URL": "https://github.com",
        "GITHUB_EVENT_NAME": "push",
        "GITHUB_WORKFLOW": "CI",
    }
    p = tracer.github_provenance(env)
    assert p["source_url"] == "https://github.com/RTrentJones/BAMCP/actions/runs/42"
    assert p["labels"] == {"trigger": "push", "workflow": "CI"}
    assert tracer.github_provenance({}) == {}


def test_clip_and_metric_hygiene():
    assert tracer._clip("x" * 20_000).endswith("[truncated]")
    assert tracer._clip({"a": 1}) == '{"a": 1}'
    assert tracer._finite_metrics({"a": 1, "b": float("nan"), "c": True, "d": "x", "e": 0.5}) == {
        "a": 1,
        "e": 0.5,
    }
    # Keys become Tracer-safe identifiers (its schema rejects spaces).
    assert tracer._finite_metrics({"category.Region Analysis": 0.5}) == {
        "category.Region_Analysis": 0.5
    }
    # ASCII identifiers only: Tracer's schema `\\w` doesn't match non-ASCII letters.
    assert tracer._finite_metrics({"präzision": 0.9}) == {"pr_zision": 0.9}
    assert tracer._clamp01("bad") is None
    assert tracer._clamp01(1.7) == 1.0


class _Resp(io.BytesIO):
    status = 201

    def __enter__(self) -> _Resp:
        return self

    def __exit__(self, *a: object) -> None:
        return None


def test_post_refuses_a_non_http_url():
    ok, msg = tracer.post({}, url="file:///etc/passwd", token="t")
    assert not ok
    assert "non-http" in msg


def test_publish_skips_without_token(capsys: pytest.CaptureFixture[str]):
    assert tracer.publish([{"mode": "truthset"}], url="https://t", token=None) == 0
    assert "skipping" in capsys.readouterr().out


def test_publish_skip_is_written_to_the_step_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    tracer.publish([{"mode": "truthset"}], url="https://t", token="")
    assert summary.read_text(encoding="utf-8") == (
        "⚠️ Tracer publish skipped: TRACER_INGEST_TOKEN unset\n"
    )


def test_publish_ignores_an_unwritable_step_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    # A directory can't be opened for append: the summary is best-effort, publish still returns.
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path))
    assert tracer.publish([], url="https://t", token=None) == 0


def test_publish_posts_and_warns_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
):
    sent: list[Any] = []

    def fake_urlopen(req: Any, timeout: float) -> _Resp:
        sent.append(req)
        if len(sent) == 2:
            raise urllib.error.HTTPError(
                req.full_url, 422, "bad", {}, io.BytesIO(b'{"error":"invalid body"}')
            )  # type: ignore[arg-type]
        if len(sent) == 3:
            raise urllib.error.URLError("dns")
        return _Resp(b'{"id":"x"}')

    monkeypatch.setattr(tracer.urllib.request, "urlopen", fake_urlopen)
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    n = tracer.publish(
        [{"mode": "a"}, {"mode": "b"}, {"mode": "c"}], url="https://tracer.example/", token="tok"
    )
    assert n == 1
    assert sent[0].full_url == "https://tracer.example/api/ingest"
    assert sent[0].get_header("Authorization") == "Bearer tok"
    out = capsys.readouterr().out
    assert "HTTP 422" in out
    assert "::warning" in out
    lines = summary.read_text(encoding="utf-8").splitlines()
    assert lines == [
        "✅ Tracer: published a/None/None (201)",
        '⚠️ Tracer publish failed: b/None/None: HTTP 422: {"error":"invalid body"}',
        "⚠️ Tracer publish failed: c/None/None: unreachable: <urlopen error dns>",
    ]


def test_main_writes_and_posts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    rep = tmp_path / "r.json"
    rep.write_text(json.dumps(_json(_report())), encoding="utf-8")
    bad = tmp_path / "bad.json"
    bad.write_text("not json", encoding="utf-8")
    out = tmp_path / "payloads.json"
    calls: list[int] = []
    monkeypatch.setattr(
        tracer, "publish", lambda payloads, url, token: calls.append(len(payloads)) or len(payloads)
    )
    monkeypatch.setenv("TRACER_INGEST_TOKEN", "tok")
    for k in ("GITHUB_SHA", "GITHUB_REPOSITORY", "GITHUB_RUN_ID"):
        monkeypatch.delenv(k, raising=False)
    assert (
        tracer.main(
            ["truthset", str(rep), str(bad), "--out", str(out), "--post", "--env", "nightly"]
        )
        == 0
    )
    written = json.loads(out.read_text(encoding="utf-8"))
    assert len(written) == 1
    assert written[0]["env"] == "nightly"
    assert calls == [1]


def test_main_llm_prints_when_no_out_or_post(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    _write_run(tmp_path / "run", [_row("a", True)])
    assert tracer.main(["llm", str(tmp_path / "run"), str(tmp_path / "empty")]) == 0
    out = capsys.readouterr().out
    assert "no results.json or comparison" in out
    assert '"mode": "llm-eval"' in out

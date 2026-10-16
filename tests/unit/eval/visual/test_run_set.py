"""Stimulus builder, leak-lint, grader, self-tests, charts and CLI on a small run set."""

from __future__ import annotations

import asyncio
import json
import shutil
from pathlib import Path

import pytest

from bamcp.eval.visual import stimuli
from bamcp.eval.visual.cli import main
from bamcp.eval.visual.grade import grade_run_set, load_jsonl
from bamcp.eval.visual.leaklint import lint_run_set, png_text_chunks
from bamcp.eval.visual.selftest import oracle_answer, run_selftests

from .conftest import FakeViewerRenderer, tiny_png


def _copy(run_set: Path, tmp_path: Path) -> Path:
    dst = tmp_path / "rs"
    shutil.copytree(run_set, dst)
    return dst


def _runs(run_dir: Path) -> list[dict]:
    return load_jsonl(run_dir / "runs.jsonl")


def _labels(run_dir: Path) -> dict[str, dict]:
    return {r["item"]: r for r in load_jsonl(run_dir / "labels" / "labels.jsonl")}


def test_layout(run_set: Path) -> None:
    runs = _runs(run_set)
    labels = _labels(run_set)
    assert len(labels) == 18
    assert len(runs) == 18 * len(stimuli.CONDITIONS) * 2 * 2
    assert {r["condition"] for r in runs} == set(stimuli.CONDITIONS)
    for r in runs[:50]:
        for f in r["files"]:
            assert (run_set / f).exists()
        assert r["answer_path"].startswith("answers/")
    # Labels never inside stim/; dirs opaque.
    assert not list((run_set / "stim").rglob("labels*"))
    assert all(len(d.name) == 8 for d in (run_set / "stim").iterdir())
    build = json.loads((run_set / "build.json").read_text())
    assert build["pngs"] is True and build["n_items"] == 18


def test_blinded_capture_contract(run_set: Path) -> None:
    calls = FakeViewerRenderer.calls
    assert calls and all(c["blind"] is True and c["mode"] == "expanded" for c in calls)
    assert all(c["viewport"] == (stimuli.VIEW_WIDTH, None) for c in calls)
    payload = calls[0]["payload"]
    assert payload["variants"] == [] and "variant_evidence" not in payload


def test_traw_content(run_set: Path) -> None:
    labels = _labels(run_set)
    by_case = {}
    for r in _runs(run_set):
        if r["condition"] == "Traw":
            by_case[labels[r["item"]]["case_type"]] = run_set / r["stim_dir"] / "traw.md"
    sb = by_case["strand_bias"].read_text()
    assert "## Alleles at" in sb and "forward" in sb
    assert "Allele combinations" in by_case["phased_cis"].read_text()
    assert "| - |" in by_case["real_del"].read_text() or " - " in by_case["real_del"].read_text()
    assert "+" in by_case["real_ins"].read_text()
    tjson = json.loads(
        next(
            (run_set / r["stim_dir"] / "tjson.json").read_text()
            for r in _runs(run_set)
            if r["condition"] == "Tjson"
        )
    )
    assert tjson["variants"] == [] and all(n.startswith("r") for n in tjson["reads"]["name"])


def test_tfull_has_tool_text(run_set: Path) -> None:
    r = next(r for r in _runs(run_set) if r["condition"] == "Tfull")
    text = (run_set / r["stim_dir"] / "tfull.md").read_text()
    assert "get_region_summary" in text and "get_variant_curation_summary" in text


def test_lint_clean(run_set: Path) -> None:
    rep = lint_run_set(run_set)
    assert rep.ok, rep.problems
    assert rep.files_checked > 0


@pytest.mark.parametrize(
    "leak",
    [
        "vocab",
        "dirname",
        "readname",
        "itemid",
        "extra",
        "caption",
        "prompt",
        "png",
        "labels",
        "nested",
    ],
)
def test_lint_catches_planted_leak(run_set: Path, tmp_path: Path, leak: str) -> None:
    import pysam

    rs = _copy(run_set, tmp_path)
    runs = _runs(rs)
    traw_dir = rs / next(r["stim_dir"] for r in runs if r["condition"] == "Traw")
    v_dir = rs / next(r["stim_dir"] for r in runs if r["condition"] == "V")
    if leak == "vocab":
        (traw_dir / "traw.md").write_text("alt reads show strand bias\n")
    elif leak == "dirname":
        traw_dir.rename(rs / "stim" / "strand_bias_0")
    elif leak == "readname":
        bam = json.loads((rs / "build.json").read_text())["bam"]
        with pysam.AlignmentFile(bam) as f:
            name = next(f.fetch(until_eof=True)).query_name
        with (traw_dir / "traw.md").open("a") as f:
            f.write(f"\nsee {name}\n")
    elif leak == "itemid":
        with (traw_dir / "traw.md").open("a") as f:
            f.write(f"\n{runs[0]['item']}\n")
    elif leak == "extra":
        (traw_dir / "answer_key.txt").write_text("x")
    elif leak == "caption":
        with (traw_dir / "caption.md").open("a") as f:
            f.write("- Hint: look closely\n")
    elif leak == "prompt":
        with (traw_dir / "prompt.md").open("a") as f:
            f.write("\nThis one is tricky.\n")
    elif leak == "png":
        (v_dir / "view_wide.png").write_bytes(tiny_png("homopolymer"))
    elif leak == "labels":
        shutil.copy(rs / "labels" / "labels.jsonl", rs / "stim" / "labels.jsonl")
    elif leak == "nested":
        (traw_dir / "sub").mkdir()
    rep = lint_run_set(rs)
    assert not rep.ok


def test_lint_missing_stim(tmp_path: Path) -> None:
    assert not lint_run_set(tmp_path).ok


def test_png_text_chunks() -> None:
    assert "homopolymer" in png_text_chunks(tiny_png("homopolymer"))
    assert "homopolymer" not in png_text_chunks(tiny_png())


def test_selftests_pass(run_set: Path, tmp_path: Path) -> None:
    rs = _copy(run_set, tmp_path)
    results = run_selftests(rs, charts=True)
    assert all(r.ok for r in results), [(r.name, r.detail) for r in results]
    charts = rs / "selftest" / "oracle" / "results"
    for name in ("accuracy_dots.png", "reliability.png", "gallery.png", "summary.json"):
        assert (charts / name).exists()
    assert list(charts.glob("forest_*.png")) and list(charts.glob("confusion_*.png"))


def test_grade_contrasts_and_exclusions(run_set: Path, tmp_path: Path) -> None:
    rs = _copy(run_set, tmp_path)
    labels = _labels(rs)
    runs = _runs(rs)
    # V-containing conditions answer perfectly, text conditions say "no" always.
    for r in runs:
        ans = oracle_answer(labels[r["item"]])
        if "V" not in r["condition"].split("+"):
            ans.update(variant_present="no", artifact_verdict="artifact", confidence=0.5)
        (rs / r["answer_path"]).write_text(json.dumps(ans))
    excluded = runs[0]["run_id"]
    s = grade_run_set(rs, exclude={excluded}, charts=False)
    assert s["n_excluded"] == 1 and s["n_runs"] == len(runs) - 1
    c = s["contrasts"]["primary|__pooled__|V+Traw vs Traw"]
    assert c["delta_joint"] > 0 and c["ci95"][0] > 0 and c["mcnemar"]["p"] < 0.05
    assert "primary|m1|V vs Tjson" in s["contrasts"]
    g = s["groups"]["m1|V"]
    assert g["primary_joint_accuracy"] == 1.0 and g["artifact_type_macro_f1"] == 1.0
    summary = json.loads((rs / "results" / "summary.json").read_text())
    assert summary["models"] == ["m1", "m2"]
    assert (rs / "results" / "results.csv").read_text().count("\n") == len(runs)


def test_missing_answers_are_wrong(run_set: Path, tmp_path: Path) -> None:
    rs = _copy(run_set, tmp_path)
    s = grade_run_set(rs, charts=False)
    g = next(iter(s["groups"].values()))
    assert g["n_missing"] == g["n_runs"] and g["primary_joint_accuracy"] == 0.0
    assert g["brier"] is None


def test_cli_end_to_end(fixture_dir: Path, tmp_path: Path, capsys) -> None:
    out = tmp_path / "cli"
    assert (
        main(
            [
                "build",
                "--out",
                str(out),
                "--seed",
                "7",
                "--instances",
                "1",
                "--fixture-dir",
                str(fixture_dir),
                "--conditions",
                "Traw",
                "V",
                "--no-png",
                "--n-items",
                "5",
                "--models",
                "a",
            ]
        )
        == 0
    )
    assert "SKIP PNG step" in capsys.readouterr().out
    assert len(_runs(out)) == 5
    assert main(["lint", "--run-dir", str(out)]) == 0
    assert main(["selftest", "--run-dir", str(out), "--no-charts"]) == 0
    ex = tmp_path / "exclude.txt"
    ex.write_text(_runs(out)[0]["run_id"] + "\n")
    assert (
        main(
            [
                "grade",
                "--run-dir",
                str(out),
                "--no-charts",
                "--exclude",
                str(ex),
                "--answers-dir",
                str(out / "selftest" / "random" / "answers"),
            ]
        )
        == 0
    )
    assert "joint=" in capsys.readouterr().out
    assert (
        main(["fixtures", "--seed", "7", "--instances", "1", "--fixture-dir", str(fixture_dir)])
        == 0
    )


def test_cli_dry_without_png(fixture_dir: Path, tmp_path: Path) -> None:
    assert (
        main(
            [
                "dry",
                "--out",
                str(tmp_path / "dry"),
                "--seed",
                "7",
                "--instances",
                "1",
                "--fixture-dir",
                str(fixture_dir),
                "--no-png",
            ]
        )
        == 0
    )


def test_unknown_condition_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        stimuli.build_run_set(tmp_path, conditions=["X"])


def test_chromium_unavailable_without_bundle(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(stimuli, "__file__", str(tmp_path / "a" / "b" / "c" / "stimuli.py"))
    assert "viewer bundle missing" in (asyncio.run(stimuli.chromium_available()) or "")


def test_stale_viewer_bundle_refuses_png_stimuli(monkeypatch, fixture_dir: Path, tmp_path: Path):
    """PNGs from a bundle built from older sources would show the old encoding while the
    prompt's VIEWER_KEY describes the new one, so the build stops instead."""

    async def _available() -> None:
        return None

    monkeypatch.setattr(stimuli, "chromium_available", _available)
    monkeypatch.setattr(stimuli, "viewer_bundle_problem", lambda: "content hash mismatch")
    with pytest.raises(RuntimeError, match="content hash mismatch"):
        stimuli.build_run_set(
            tmp_path,
            seed=7,
            instances=1,
            conditions=["V"],
            fixture_dir=fixture_dir,
            log=lambda *_: None,
        )


def test_viewer_bundle_problem_reports_stale_sources(monkeypatch) -> None:
    from bamcp.eval import render_cli

    monkeypatch.setattr(render_cli, "_check_rebuild_needed", lambda: "WARNING: stale")
    assert stimuli.viewer_bundle_problem() == "WARNING: stale"


def test_select_items_stratified() -> None:
    from bamcp.eval.visual.fixtures import build_plan

    plan = build_plan(1, 2)
    sel = stimuli.select_items(plan, 18)
    assert len({i.case_type for i in sel}) == 18 and {i.instance for i in sel} == {0}
    assert len(stimuli.select_items(plan, None)) == 36

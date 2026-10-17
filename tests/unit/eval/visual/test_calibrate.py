"""Round-4 calibration on canned B-only pilot results (design section 7)."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pytest

from bamcp.eval.visual.calibrate import (
    EXPECTED_B,
    G1_WINDOW,
    PilotRun,
    calibrate,
    estimate_cells,
    format_operating_characteristic,
    icc_binary,
    load_pilot_runs,
    model_offsets,
    simulate_operating_characteristic,
    write_calibration,
)
from bamcp.eval.visual.cli_r4 import main as cli_main
from bamcp.eval.visual.grade import HARD_STRATA_R4, load_jsonl
from bamcp.eval.visual.selftest import oracle_answer

from .r4_analysis_world import ARTIFACT, REAL, make_r4_run_set

MODELS = ("haiku", "sonnet", "opus")
ALL_PASS = {g: {"pass": True} for g in ("G4", "G5", "G6", "G8")}


def canned(acc: dict[tuple[str, str], float], n: dict[str, int] | None = None) -> list[PilotRun]:
    """Deterministic pilot: in each cell, round(acc · runs) correct runs filled item by item
    (opus first, haiku last), so runs of an item are correlated and models are ordered."""
    n = n or {"default": 6, "alt": 4}
    order = ("opus", "sonnet", "haiku")
    runs = []
    for (s, k), a in acc.items():
        total = n[k] * len(order)
        n_ok = round(a * total)
        idx = 0
        for i in range(n[k]):
            for m in order:
                runs.append(PilotRun(f"{s}:{k}:{i}", m, s, k, idx < n_ok))
                idx += 1
    return runs


def near_target(shift: float = 0.0, alt: float = -0.15) -> dict[tuple[str, str], float]:
    out = {}
    for s in HARD_STRATA_R4:
        out[(s, "default")] = min(1.0, max(0.0, EXPECTED_B[s] + shift))
        out[(s, "alt")] = min(1.0, max(0.0, EXPECTED_B[s] + shift + alt))
    return out


def test_icc_known_values() -> None:
    assert icc_binary([[[1, 1], [0, 0]]]) == pytest.approx(1.0)
    assert icc_binary([[[1, 0], [1, 0]]]) == 0.0  # negative estimate clipped
    assert icc_binary([[[1, 1, 1]]]) == 0.0  # one item: no between-item df
    # Two cells with different means: cell means are removed first.
    assert icc_binary([[[1, 1], [1, 1]], [[0, 0], [0, 0]]]) == 0.0
    mixed = icc_binary([[[1, 1, 0], [0, 0, 0], [1, 1, 1], [1, 0, 0]]])
    assert 0 < mixed < 1


def test_shrinkage_and_effective_n() -> None:
    runs = canned(near_target())
    cells, est = estimate_cells(runs)
    assert len(cells) == 16 and 0 <= est["rho"] <= 1
    mu = est["mu"]
    for c in cells.values():
        assert c.n_eff == pytest.approx(c.n_runs / (1 + 2 * est["rho"]))
        lo, hi = sorted((c.acc, mu))
        assert lo - 1e-9 <= c.mean <= hi + 1e-9  # shrunk toward the pooled mean
        a, b = c.interval()
        assert a < c.mean < b
    with pytest.raises(ValueError):
        estimate_cells([])


def test_model_offsets_order() -> None:
    off = model_offsets(canned(near_target()))
    assert off["opus"] > off["sonnet"] > off["haiku"]


def test_in_window_pilot_goes_only_with_all_gates() -> None:
    runs = canned(near_target())
    rec = calibrate(runs, gates=ALL_PASS, n_draws=2000)
    assert G1_WINDOW[0] <= rec["prediction"]["mean"] <= G1_WINDOW[1]
    assert rec["gates"]["G1"]["pass"] and rec["gates"]["G3"]["pass"]
    assert rec["gates"]["G3"]["model"] == "opus"
    assert rec["decision"] == "GO" and rec["gates"]["G2"]["moves"] == []
    assert rec["frozen_knobs"] == dict.fromkeys(HARD_STRATA_R4, "default")
    assert rec["reads_conditions"] == ["Traw-r4"]
    lo, hi = rec["prediction"]["pi90"]
    assert lo < rec["prediction"]["mean"] < hi
    assert calibrate(runs, n_draws=500)["decision"] == "PENDING"
    held = calibrate(runs, gates={**ALL_PASS, "G4": {"pass": False}}, n_draws=500)
    assert held["decision"] == "HOLD" and "G4" in held["reasons"][0]
    # G6 fails on the valid-answer rate measured here, whatever the harness says.
    bad = [PilotRun(r.item, r.model, r.stratum, r.knob, r.correct, valid=i % 10 != 0)
           for i, r in enumerate(runs)]  # fmt: skip
    g6 = calibrate(bad, gates=ALL_PASS, n_draws=500)["gates"]["G6"]
    assert g6["pass"] is False and g6["valid_rate"] < 0.95


def test_too_easy_pilot_moves_knobs_only_to_piloted_alt() -> None:
    runs = canned(near_target(shift=+0.25, alt=-0.30))
    rec = calibrate(runs, gates=ALL_PASS, n_draws=2000)
    assert rec["prediction_initial"]["mean"] > G1_WINDOW[1]
    moves = rec["gates"]["G2"]["moves"]
    assert moves and all(m["from"] == "default" and m["to"] == "alt" for m in moves)
    assert rec["prediction"]["mean"] < rec["prediction_initial"]["mean"]
    assert rec["decision"] in ("GO", "RE-PILOT")
    if rec["decision"] == "GO":
        assert G1_WINDOW[0] <= rec["prediction"]["mean"] <= G1_WINDOW[1]
    # Alt levels that go the wrong way are never chosen: every move lowers the shrunken
    # mean (a clipped stratum's 4-item alt cell can still shrink lower than its default).
    wrong_way = canned(near_target(shift=+0.25, alt=+0.05))
    rec2 = calibrate(wrong_way, gates=ALL_PASS, n_draws=1000)
    for m in rec2["gates"]["G2"]["moves"]:
        cells = rec2["cells"]
        alt, cur = cells[f"{m['stratum']}|alt"], cells[f"{m['stratum']}|default"]
        assert alt["shrunken_mean"] < cur["shrunken_mean"]
    assert len(rec2["gates"]["G2"]["moves"]) <= 1 and rec2["decision"] == "RE-PILOT"
    assert rec2["repilot_plan"]["items"] == 20


def test_repilot_never_moves_and_ends_in_no_go() -> None:
    runs = canned(near_target(shift=+0.25, alt=+0.05))
    first = calibrate(runs, gates=ALL_PASS, n_draws=1000)
    second = calibrate(runs, gates=ALL_PASS, previous=first, repilot=True, n_draws=1000)
    assert second["gates"]["G2"]["moves"] == [] and second["decision"] == "NO-GO"
    ok = calibrate(canned(near_target()), gates=ALL_PASS, previous=first, repilot=True)
    assert ok["decision"] == "GO" and ok["settings_initial"] == first["frozen_knobs"]


def test_g3_opus_ceiling() -> None:
    # Everything correct for opus and sonnet, haiku mostly wrong: primary in window but
    # opus far above 0.90.
    runs = []
    for s in HARD_STRATA_R4:
        for k, n in (("default", 6), ("alt", 4)):
            for i in range(n):
                for m in MODELS:
                    ok = m != "haiku" or i % 3 == 0
                    runs.append(PilotRun(f"{s}:{k}:{i}", m, s, k, ok))
    rec = calibrate(runs, gates=ALL_PASS, n_draws=500)
    assert rec["gates"]["G3"]["value"] > 0.90 and rec["gates"]["G3"]["pass"] is False
    assert rec["decision"] in ("RE-PILOT", "GO")
    assert rec["decision"] == "RE-PILOT" or rec["gates"]["G3"]["pass"]
    named = calibrate(runs, gates=ALL_PASS, n_draws=500, opus_model="sonnet")
    assert named["gates"]["G3"]["model"] == "sonnet"


def test_load_pilot_runs_reads_baseline_only(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        return oracle_answer(lab) if lab["polarity"] == "pos" else REAL

    rs = make_r4_run_set(tmp_path / "rs", n=6, answer_fn=fn)
    # Treatment answer files become directories: reading one would raise.
    for r in load_jsonl(rs / "runs.jsonl"):
        if r["condition"] != "Traw-r4":
            p = rs / r["answer_path"]
            p.unlink()
            p.mkdir()
    runs, caption = load_pilot_runs([rs])
    assert runs and not caption
    assert {r.model for r in runs} == {"haiku", "sonnet", "opus"}
    hard = [r for r in runs if r.stratum in HARD_STRATA_R4]
    assert all(r.correct == r.stratum.endswith("_pos") for r in hard)
    assert {r.knob for r in runs} == {"default", "alt"}
    rec = calibrate(runs, gates=ALL_PASS, caption_runs=None, n_draws=500)
    assert rec["prediction"]["mean"] == pytest.approx(0.5, abs=0.15)
    out = tmp_path / "cal.json"
    write_calibration(rec, out)
    assert json.loads(out.read_text())["decision"] == rec["decision"]


def test_caption_check_runs_are_separated(tmp_path: Path) -> None:
    rs = make_r4_run_set(tmp_path / "rs", n=6, conditions=("Traw-r4",))
    labels = load_jsonl(rs / "labels" / "labels.jsonl")
    first = labels[0]["item"]
    for lab in labels:
        if lab["item"] == first:  # a round-3 real item under the new caption
            for k in ("stratum", "pair", "polarity"):
                lab.pop(k)
            lab["case_type"] = "real_snv_het"
    (rs / "labels" / "labels.jsonl").write_text("".join(json.dumps(x) + "\n" for x in labels))
    runs, caption = load_pilot_runs([rs])
    assert len(caption) == 3 and all(c.stratum == "real_snv_het" for c in caption)
    rec = calibrate(runs, caption_runs=caption, n_draws=200)
    assert rec["caption_check"]["joint_accuracy"] == 1.0


def test_simulated_operating_characteristic() -> None:
    oc = simulate_operating_characteristic(12, seed=2, n_draws=300)
    assert sum(oc["branches"].values()) == pytest.approx(1.0)
    text = format_operating_characteristic(oc)
    assert "branch frequencies" in text and "NO-GO" in text
    go = oc["realised_main_primary_all_go"]
    if go:
        assert 0.3 < go["mean"] < 0.95


def test_cli_calibrate_and_simulate(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rs = make_r4_run_set(tmp_path / "rs", n=6, conditions=("Traw-r4",))
    gates = tmp_path / "gates.json"
    gates.write_text(json.dumps({"G4": True, "G5": {"pass": True}, "G6": {"pass": True}}))
    g8 = tmp_path / "g8.json"
    g8.write_text(json.dumps({"g8": {"pass": True}}))
    assert (
        cli_main(["calibrate", "--run-dir", str(rs), "--gates", str(gates), "--g8", str(g8)]) == 0
    )
    rec = json.loads((rs / "calibration.json").read_text())
    assert rec["gates"]["G8"]["pass"] is True and rec["gates"]["G4"]["pass"] is True
    # Oracle answers everywhere: far too easy, no alt helps → RE-PILOT.
    assert rec["decision"] == "RE-PILOT"
    out = capsys.readouterr().out
    assert "decision: RE-PILOT" in out
    assert cli_main(["calibrate", "--simulate", "--n-sims", "3"]) == 0
    assert "branch frequencies" in capsys.readouterr().out
    assert cli_main(["calibrate"]) == 2
    shutil.rmtree(rs)


def test_missing_cell_level_falls_back() -> None:
    acc = {k: v for k, v in near_target().items() if not (k[0] == "tr_indel_pos" and k[1] == "alt")}
    rec = calibrate(canned(acc), gates=ALL_PASS, n_draws=300)
    assert "tr_indel_pos|alt" not in rec["cells"]
    assert not math.isnan(rec["prediction"]["mean"])


def test_artifact_answers_make_negatives_right(tmp_path: Path) -> None:
    rs = make_r4_run_set(
        tmp_path / "rs", n=6, conditions=("Traw-r4",), answer_fn=lambda *_: ARTIFACT
    )
    runs, _ = load_pilot_runs([rs])
    assert all(r.correct == r.stratum.endswith("_neg") for r in runs if r.stratum in HARD_STRATA_R4)


def test_g3_move_never_pushes_primary_out_of_g1() -> None:
    """Lowering strata for G3 stops (and undoes the move) rather than breaking G1."""
    from bamcp.eval.visual.calibrate import Cell, apply_knob_rule

    cells = {}
    for s in HARD_STRATA_R4:
        cells[(s, "default")] = Cell(s, "default", 6, 18, 0.7, 9.0, 63.0, 27.0)
        cells[(s, "alt")] = Cell(s, "alt", 4, 12, 0.5, 6.0, 30.0, 30.0)

    def predict(settings: dict[str, str]) -> dict:
        n_alt = sum(v == "alt" for v in settings.values())
        return {
            "mean": 0.63 - 0.02 * n_alt,  # one move takes G1 out of its window
            "per_model": {"opus": 0.93 - 0.015 * n_alt},
            "per_stratum": dict(EXPECTED_B),
        }

    start = dict.fromkeys(HARD_STRATA_R4, "default")
    settings, moves = apply_knob_rule(cells, start, predict, "opus")
    assert settings == start and moves == []

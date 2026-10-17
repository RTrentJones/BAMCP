"""Round-4 difficulty calibration: pilot B results → knob rule → GO / NO-GO (design section 7).

Reads **only** the baseline (``Traw-r4``) runs of a pilot run set: :func:`load_pilot_runs`
filters ``runs.jsonl`` by condition before any answer file is opened, so a treatment answer
can never inform the calibration.

Estimation. Each hard stratum × knob level is a cell. Its accuracy is shrunk toward the
pooled mean by beta-binomial empirical Bayes: the 3 model runs of an item are correlated,
so a cell's ``n`` is its effective size ``n_runs / (1 + (m - 1) ρ)`` with ρ the within-item
intra-class correlation (ANOVA estimator, pooled over cells); the prior's concentration κ
comes from the method of moments over cells. The predicted main primary B (mean over the 4
pairs of balanced accuracy) at a knob setting uses the posterior means; its predictive
interval adds posterior uncertainty and the main draw's sampling noise at the frozen quotas
(15 items × the pilot's models per stratum). Per-model predictions (G3) shift each stratum on
the logit scale by the model's pilot offset.

Frozen knob rule (G2, then G1/G3). A stratum moves only to its piloted alternative level and
only in the needed direction: (G2) when its shrunken 90% interval lies wholly outside
0.45-0.88; then, if the predicted primary misses 0.62-0.78 (or opus's exceeds 0.90), the
strata farthest from their section-3 expected B in the needed direction move one at a time
until it is fixed. The rule is applied once; if G1/G3 still fail the decision is a B-only
re-pilot (20 items × the models), and after a re-pilot (``repilot=True``) NO-GO.

G4 (label audit), G5 (pools) and G8 (recoverability, :mod:`.sufficiency`) are computed
elsewhere and passed in as ``gates`` (``{"G4": {"pass": bool, ...}, ...}``); G6 combines the
valid-answer rate measured here with the harness checks passed in as ``gates["G6"]``. A
missing gate leaves the decision PENDING; a failed one HOLD.

``--simulate`` (:func:`simulate_operating_characteristic`) runs the whole rule on simulated
pilots under the section-7 prior (each stratum's expected B ± 0.12) and reports how often
each branch fires and the realised main-set primary B.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .grade import (
    BASELINE_R4,
    HARD_STRATA_R4,
    PAIRS_R4,
    STRATUM_PAIRS_R4,
    grade_runs_r4,
    stratum_of,
)
from .stats import beta_ppf, mean

KNOB_LEVELS = ("default", "alt")
# Section 3's expected B per stratum: the "target" the G1 fix measures distance from.
EXPECTED_B: dict[str, float] = {
    "mosaic_snv_pos": 0.62,
    "mosaic_caller_fp_neg": 0.68,
    "tr_indel_pos": 0.72,
    "tr_stutter_neg": 0.66,
    "segdup_cmrg_pos": 0.68,
    "segdup_psv_neg": 0.62,
    "indel_adjacent_misalign_neg": 0.70,
    "indel_adjacent_snv_pos": 0.82,
}
G1_WINDOW = (0.62, 0.78)
G2_BAND = (0.45, 0.88)
G2_LEVEL = 0.90
G3_MAX = 0.90
G6_MIN_VALID = 0.95
MAIN_QUOTA = 15
REPILOT_ITEMS = 20
KAPPA_BOUNDS = (2.0, 500.0)
N_DRAWS = 4000
EXTERNAL_GATES = ("G4", "G5", "G8")

# --simulate only (not part of the frozen rule): the assumed direction of each knob's
# alternative level on B accuracy, and the simulated world (power_r4b.py's model).
SIM_ALT_DIRECTION: dict[str, int] = {
    "mosaic_snv_pos": +1,  # >= 7 alt reads only: easier
    "mosaic_caller_fp_neg": -1,  # higher-VAF false positives look more real
    "tr_indel_pos": -1,  # longer runs, more stutter
    "tr_stutter_neg": -1,
    "segdup_cmrg_pos": -1,  # more "difficult" items
    "segdup_psv_neg": +1,  # lower-VAF PSVs look less like a het
    "indel_adjacent_misalign_neg": -1,  # stronger misalignment signal
    "indel_adjacent_snv_pos": -1,  # fewer alt reads
}
SIM_ALT_SHIFT = (0.03, 0.10)
SIM_PRIOR_HALF_WIDTH = 0.12
SIM_MODELS = {"haiku": -0.9, "sonnet": 0.0, "opus": 0.9}
SIM_ITEM_SD = math.sqrt(1.2**2 + 1.0**2 + 0.6**2)  # item + item×condition + cluster
SIM_PILOT = {"default": 6, "alt": 4}


@dataclass
class PilotRun:
    item: str
    model: str
    stratum: str
    knob: str
    correct: bool
    valid: bool = True


@dataclass
class Cell:
    stratum: str
    knob: str
    n_items: int
    n_runs: int
    acc: float
    n_eff: float
    a: float
    b: float

    @property
    def mean(self) -> float:
        return self.a / (self.a + self.b)

    def interval(self, level: float = G2_LEVEL) -> tuple[float, float]:
        q = (1 - level) / 2
        return beta_ppf(q, self.a, self.b), beta_ppf(1 - q, self.a, self.b)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def load_pilot_runs(
    run_dirs: Iterable[str | Path],
    *,
    answers_dirs: Iterable[str | Path | None] | None = None,
    exclude: set[str] | None = None,
) -> tuple[list[PilotRun], list[PilotRun]]:
    """B-only pilot runs → (round-4 stratum runs, caption-check runs).

    Caption-check runs are the round-3 real items re-run under the new caption: their
    labels carry no round-4 stratum; they are summarised, never used for the knob rule.
    """
    dirs = [Path(d) for d in run_dirs]
    ans = list(answers_dirs) if answers_dirs is not None else [None] * len(dirs)
    r4: list[PilotRun] = []
    caption: list[PilotRun] = []
    for d, a in zip(dirs, ans, strict=True):
        for g in grade_runs_r4(d, answers_dir=a, exclude=exclude, conditions={BASELINE_R4}):
            s = stratum_of(g.label)
            pr = PilotRun(
                item=f"{d.name}/{g.run['item']}",
                model=g.run["model"],
                stratum=s,
                knob=str(g.label.get("knob") or "default"),
                correct=g.joint,
                valid=g.status == "valid",
            )
            (r4 if s in STRATUM_PAIRS_R4 else caption).append(pr)
    return r4, caption


# ---------------------------------------------------------------------------
# Estimation
# ---------------------------------------------------------------------------


def _by_cell_item(runs: list[PilotRun]) -> dict[tuple[str, str], dict[str, list[float]]]:
    out: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in runs:
        if r.stratum in HARD_STRATA_R4:
            out[(r.stratum, r.knob)][r.item].append(float(r.correct))
    return out


def icc_binary(cells: list[list[list[float]]]) -> float:
    """One-way ANOVA intra-class correlation of runs within items, items nested in cells
    (cell means removed), clipped to [0, 1]. ``cells`` = [[item runs, ...], ...]."""
    ss_b = ss_w = 0.0
    n_runs = n_items = 0
    sum_sq = 0.0
    for items in cells:
        vals = [v for it in items for v in it]
        if not vals:
            continue
        cm = mean(vals)
        for it in items:
            if not it:
                continue
            im = mean(it)
            ss_b += len(it) * (im - cm) ** 2
            ss_w += sum((v - im) ** 2 for v in it)
            n_runs += len(it)
            n_items += 1
            sum_sq += len(it) ** 2
    n_cells = sum(1 for items in cells if any(items))
    df_b, df_w = n_items - n_cells, n_runs - n_items
    if df_b <= 0 or df_w <= 0:
        return 0.0
    msb, msw = ss_b / df_b, ss_w / df_w
    k0 = (n_runs - sum_sq / n_runs) / max(1, n_items - 1)
    denom = msb + (k0 - 1) * msw
    if denom <= 0:
        return 0.0
    return min(1.0, max(0.0, (msb - msw) / denom))


def estimate_cells(runs: list[PilotRun]) -> tuple[dict[tuple[str, str], Cell], dict[str, float]]:
    """Shrunken cells + ``{rho, mu, kappa, runs_per_item}``."""
    grouped = _by_cell_item(runs)
    rho = icc_binary([list(items.values()) for items in grouped.values()])
    raw = []
    for (s, k), items in sorted(grouped.items()):
        vals = [v for it in items.values() for v in it]
        m_bar = len(vals) / len(items)
        n_eff = len(vals) / (1 + (m_bar - 1) * rho)
        raw.append((s, k, len(items), len(vals), mean(vals), n_eff))
    if not raw:
        raise ValueError("no hard-stratum B runs in the pilot")
    w = sum(r[5] for r in raw)
    mu = min(0.99, max(0.01, sum(r[4] * r[5] for r in raw) / w))
    accs = [r[4] for r in raw]
    if len(raw) > 1:
        s2 = sum((p - mean(accs)) ** 2 for p in accs) / (len(accs) - 1)
        within = mu * (1 - mu) * mean([1 / r[5] for r in raw])
        tau2 = s2 - within
    else:
        tau2 = 0.0
    kappa = KAPPA_BOUNDS[1] if tau2 <= 0 else mu * (1 - mu) / tau2 - 1
    kappa = min(KAPPA_BOUNDS[1], max(KAPPA_BOUNDS[0], kappa))
    alpha, beta = mu * kappa, (1 - mu) * kappa
    cells = {
        (s, k): Cell(s, k, ni, nr, acc, ne, alpha + acc * ne, beta + (1 - acc) * ne)
        for s, k, ni, nr, acc, ne in raw
    }
    runs_per_item = mean([r[3] / r[2] for r in raw])
    return cells, {"rho": rho, "mu": mu, "kappa": kappa, "runs_per_item": runs_per_item}


def _logit(p: float) -> float:
    p = min(0.99, max(0.01, p))
    return math.log(p / (1 - p))


def _expit(x: float) -> float:
    return 1 / (1 + math.exp(-x))


def model_offsets(runs: list[PilotRun]) -> dict[str, float]:
    """Each model's logit offset from the all-model hard-stratum pilot accuracy."""
    hard = [r for r in runs if r.stratum in HARD_STRATA_R4]
    overall = _logit(mean([float(r.correct) for r in hard]))
    out = {}
    for m in sorted({r.model for r in hard}):
        out[m] = _logit(mean([float(r.correct) for r in hard if r.model == m])) - overall
    return out


def _primary(per_stratum: dict[str, Any]) -> Any:
    """Mean over pairs of (pos + neg) / 2 — works on floats or numpy arrays."""
    pairs: dict[str, dict[str, Any]] = defaultdict(dict)
    for s, v in per_stratum.items():
        p, pol = STRATUM_PAIRS_R4[s]
        pairs[p][pol] = v
    vals = [(pairs[p]["pos"] + pairs[p]["neg"]) / 2 for p in PAIRS_R4]
    return sum(vals) / len(vals)


def _cell_for(cells: dict[tuple[str, str], Cell], s: str, knob: str) -> Cell:
    if (s, knob) in cells:
        return cells[(s, knob)]
    other = [c for (cs, _), c in cells.items() if cs == s]
    if not other:
        raise ValueError(f"no pilot cell for stratum {s}")
    return other[0]  # a level never piloted cannot be chosen; this only fills a gap


def predict_primary(
    cells: dict[tuple[str, str], Cell],
    settings: dict[str, str],
    *,
    offsets: dict[str, float],
    rho: float,
    n_models: int,
    n_main: int = MAIN_QUOTA,
    n_draws: int = N_DRAWS,
    seed: int = 0,
) -> dict[str, Any]:
    """Predicted main primary B at ``settings`` (stratum → knob level)."""
    rng = np.random.default_rng(seed)
    means = {s: _cell_for(cells, s, settings[s]).mean for s in HARD_STRATA_R4}
    n_eff_main = max(1, round(n_main * n_models / (1 + (n_models - 1) * rho)))
    draws = {}
    for s in HARD_STRATA_R4:
        c = _cell_for(cells, s, settings[s])
        theta = rng.beta(c.a, c.b, size=n_draws)
        draws[s] = rng.binomial(n_eff_main, theta) / n_eff_main
    pred = np.asarray(_primary(draws))
    per_model = {
        m: float(_primary({s: _expit(_logit(v) + off) for s, v in means.items()}))
        for m, off in offsets.items()
    }
    return {
        "mean": float(_primary(means)),
        "pi90": [float(np.quantile(pred, 0.05)), float(np.quantile(pred, 0.95))],
        "pi95": [float(np.quantile(pred, 0.025)), float(np.quantile(pred, 0.975))],
        "per_model": per_model,
        "per_stratum": means,
        "n_eff_main_per_stratum": n_eff_main,
    }


# ---------------------------------------------------------------------------
# Frozen knob rule and gates
# ---------------------------------------------------------------------------


def _opus(offsets: dict[str, float], opus_model: str | None) -> str | None:
    if opus_model:
        return opus_model if opus_model in offsets else None
    named = [m for m in offsets if "opus" in m.lower()]
    return named[0] if named else None


def _g1_g3(pred: dict[str, Any], opus: str | None) -> tuple[bool, bool]:
    g1 = G1_WINDOW[0] <= pred["mean"] <= G1_WINDOW[1]
    g3 = opus is None or pred["per_model"].get(opus, 0.0) <= G3_MAX
    return g1, g3


def apply_knob_rule(
    cells: dict[tuple[str, str], Cell],
    settings: dict[str, str],
    predict: Callable[[dict[str, str]], dict[str, Any]],
    opus: str | None,
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    """Section 7's frozen rule, applied once. Returns (new settings, moves)."""
    settings = dict(settings)
    moves: list[dict[str, Any]] = []

    def alt_of(s: str) -> str | None:
        alt = "alt" if settings[s] == "default" else "default"
        return alt if (s, alt) in cells else None

    def move(s: str, reason: str) -> None:
        to = alt_of(s)
        if to is None:  # callers check first; a level never piloted is never chosen
            return
        moves.append({"stratum": s, "from": settings[s], "to": to, "reason": reason})
        settings[s] = to

    # G2: a stratum whose 90% interval lies wholly outside the band.
    for s in HARD_STRATA_R4:
        cur = _cell_for(cells, s, settings[s])
        lo, hi = cur.interval()
        to = alt_of(s)
        if to is None:
            continue
        alt = cells[(s, to)]
        if hi < G2_BAND[0] and alt.mean > cur.mean:
            move(s, f"G2: 90% interval [{lo:.2f}, {hi:.2f}] below {G2_BAND[0]}")
        elif lo > G2_BAND[1] and alt.mean < cur.mean:
            move(s, f"G2: 90% interval [{lo:.2f}, {hi:.2f}] above {G2_BAND[1]}")

    # G1 / G3: farthest-from-target strata in the needed direction, one at a time.
    pred = predict(settings)
    g1, g3 = _g1_g3(pred, opus)
    if g1 and g3:
        return settings, moves
    up = pred["mean"] < G1_WINDOW[0]
    moved = {m["stratum"] for m in moves}
    dev = {s: pred["per_stratum"][s] - EXPECTED_B[s] for s in HARD_STRATA_R4}
    order = sorted(HARD_STRATA_R4, key=lambda s: dev[s], reverse=not up)
    for s in order:
        if s in moved or alt_of(s) is None:
            continue
        cur, alt = _cell_for(cells, s, settings[s]), cells[(s, str(alt_of(s)))]
        if (alt.mean > cur.mean) != up or alt.mean == cur.mean:
            continue
        need = "raise" if up else "lower"
        move(s, f"G1/G3: {need} predicted primary {pred['mean']:.3f}")
        pred = predict(settings)
        g1, g3 = _g1_g3(pred, opus)
        if g1 and g3:
            break
    return settings, moves


def _gate(gates: dict[str, Any], name: str) -> dict[str, Any]:
    g = gates.get(name)
    if g is None:
        return {"pass": None, "detail": "not provided"}
    if isinstance(g, bool):
        return {"pass": g}
    if "pass" not in g and name == "G8" and "g8" in g:  # a sufficiency.py output file
        g = g["g8"]
    return {**g, "pass": g.get("pass")}


def calibrate(
    runs: list[PilotRun],
    *,
    gates: dict[str, Any] | None = None,
    caption_runs: list[PilotRun] | None = None,
    previous: dict[str, Any] | None = None,
    repilot: bool = False,
    opus_model: str | None = None,
    n_draws: int = N_DRAWS,
    seed: int = 0,
) -> dict[str, Any]:
    """The calibration record (``calibration.json``) for B-only pilot runs."""
    gates = gates or {}
    cells, est = estimate_cells(runs)
    offsets = model_offsets(runs)
    n_models = max(1, len(offsets))
    opus = _opus(offsets, opus_model)

    def predict(settings: dict[str, str]) -> dict[str, Any]:
        return predict_primary(
            cells,
            settings,
            offsets=offsets,
            rho=est["rho"],
            n_models=n_models,
            n_draws=n_draws,
            seed=seed,
        )

    start = dict(previous["frozen_knobs"]) if previous else dict.fromkeys(HARD_STRATA_R4, "default")
    pred0 = predict(start)
    moves: list[dict[str, Any]]
    if repilot:
        settings, moves = start, []  # the rule was applied once already
    else:
        settings, moves = apply_knob_rule(cells, start, predict, opus)
    pred = predict(settings) if moves else pred0
    g1, g3 = _g1_g3(pred, opus)

    hard = [r for r in runs if r.stratum in STRATUM_PAIRS_R4]
    valid_rate = mean([float(r.valid) for r in hard]) if hard else float("nan")
    harness = _gate(gates, "G6")
    g6_pass: bool | None = valid_rate >= G6_MIN_VALID
    if g6_pass and harness["pass"] is not True:
        g6_pass = None if harness["pass"] is None else False
    out_of_band = []
    for s in HARD_STRATA_R4:
        lo, hi = _cell_for(cells, s, settings[s]).interval()
        if hi < G2_BAND[0] or lo > G2_BAND[1]:
            out_of_band.append({"stratum": s, "interval90": [lo, hi]})
    gate_res: dict[str, Any] = {
        "G1": {"pass": g1, "value": pred["mean"], "window": list(G1_WINDOW)},
        "G2": {"moves": moves, "band": list(G2_BAND), "remaining_out_of_band": out_of_band},
        "G3": {
            "pass": g3,
            "model": opus,
            "value": pred["per_model"].get(opus) if opus else None,
            "max": G3_MAX,
        },
        "G6": {
            "pass": g6_pass,
            "valid_rate": valid_rate,
            "min_valid": G6_MIN_VALID,
            "harness": harness,
        },
        **{g: _gate(gates, g) for g in EXTERNAL_GATES},
    }
    reasons: list[str] = []
    if not (g1 and g3):
        decision = "NO-GO" if repilot else "RE-PILOT"
        if not g1:
            reasons.append(f"G1: predicted primary {pred['mean']:.3f} outside {G1_WINDOW}")
        if not g3:
            reasons.append(f"G3: {opus} predicted {gate_res['G3']['value']:.3f} > {G3_MAX}")
    else:
        failed = [g for g in ("G4", "G5", "G6", "G8") if gate_res[g]["pass"] is False]
        pending = [g for g in ("G4", "G5", "G6", "G8") if gate_res[g]["pass"] is None]
        if failed:
            decision = "HOLD"
            reasons.append(f"failed: {', '.join(failed)} (apply the gate's remedy, section 7)")
        elif pending:
            decision = "PENDING"
            reasons.append(f"not provided: {', '.join(pending)}")
        else:
            decision = "GO"
    record: dict[str, Any] = {
        "reads_conditions": [BASELINE_R4],
        "repilot": repilot,
        "n_runs": len(runs),
        "n_items": len({r.item for r in runs}),
        "models": sorted(offsets),
        "estimation": est,
        "model_offsets": offsets,
        "cells": {
            f"{s}|{k}": {
                **{f: v for f, v in asdict(c).items() if f not in ("stratum", "knob")},
                "shrunken_mean": c.mean,
                "interval90": list(c.interval()),
            }
            for (s, k), c in sorted(cells.items())
        },
        "settings_initial": start,
        "prediction_initial": pred0,
        "frozen_knobs": settings,
        "prediction": pred,
        "gates": gate_res,
        "decision": decision,
        "reasons": reasons,
    }
    if decision == "RE-PILOT":
        record["repilot_plan"] = {
            "items": REPILOT_ITEMS,
            "models": sorted(offsets),
            "knobs": settings,
        }
    if caption_runs:
        record["caption_check"] = {
            "n_runs": len(caption_runs),
            "joint_accuracy": mean([float(r.correct) for r in caption_runs]),
            "per_model": {
                m: mean([float(r.correct) for r in caption_runs if r.model == m])
                for m in sorted({r.model for r in caption_runs})
            },
        }
    return record


def write_calibration(record: dict[str, Any], path: str | Path) -> None:
    from .grade import _clean

    Path(path).write_text(json.dumps(_clean(record), indent=2) + "\n")


# ---------------------------------------------------------------------------
# Simulated operating characteristic (section 7 firewall)
# ---------------------------------------------------------------------------


def _marginal(mu: float, sd: float, offsets: Iterable[float]) -> float:
    """Approximate mean over models of E[expit(mu + b + sd Z)] (logistic-normal)."""
    scale = math.sqrt(1 + math.pi * sd**2 / 8)
    offs = list(offsets)
    return mean([_expit((mu + b) / scale) for b in offs])


def _solve_mu(theta: float, sd: float, offsets: Iterable[float]) -> float:
    offs = list(offsets)
    lo, hi = -10.0, 10.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if _marginal(mid, sd, offs) < theta:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _sim_items(
    rng: np.random.Generator, mu: float, n: int, prefix: str, stratum: str, knob: str
) -> list[PilotRun]:
    out = []
    for i in range(n):
        u = rng.normal(0, SIM_ITEM_SD)
        for m, b in SIM_MODELS.items():
            p = _expit(mu + b + u)
            out.append(PilotRun(f"{prefix}{i}", m, stratum, knob, bool(rng.random() < p)))
    return out


def simulate_operating_characteristic(
    n_sims: int = 500, *, seed: int = 0, n_draws: int = 1000
) -> dict[str, Any]:
    """Branch frequencies and realised main primary B under the section-7 prior."""
    rng = np.random.default_rng(seed)
    passing = {g: {"pass": True} for g in (*EXTERNAL_GATES, "G6")}
    branches: dict[str, int] = defaultdict(int)
    realised: dict[str, list[float]] = defaultdict(list)
    true_b: list[float] = []
    offs = list(SIM_MODELS.values())
    for k in range(n_sims):
        theta = {}
        for s in HARD_STRATA_R4:
            d = rng.uniform(
                EXPECTED_B[s] - SIM_PRIOR_HALF_WIDTH, EXPECTED_B[s] + SIM_PRIOR_HALF_WIDTH
            )
            a = d + SIM_ALT_DIRECTION[s] * rng.uniform(*SIM_ALT_SHIFT)
            theta[(s, "default")] = min(0.97, max(0.05, d))
            theta[(s, "alt")] = min(0.97, max(0.05, a))
        mu = {c: _solve_mu(t, SIM_ITEM_SD, offs) for c, t in theta.items()}
        runs: list[PilotRun] = []
        for (s, knob), m in mu.items():
            runs += _sim_items(rng, m, SIM_PILOT[knob], f"p{k}:{s}:{knob}:", s, knob)
        rec = calibrate(runs, gates=passing, n_draws=n_draws, seed=k)
        branch = "GO" if not rec["gates"]["G2"]["moves"] else "GO after knob move"
        if rec["decision"] == "RE-PILOT":
            settings = rec["frozen_knobs"]
            extra: list[PilotRun] = []
            for j, s in enumerate(HARD_STRATA_R4):
                n = REPILOT_ITEMS // len(HARD_STRATA_R4) + (j < REPILOT_ITEMS % len(HARD_STRATA_R4))
                extra += _sim_items(rng, mu[(s, settings[s])], n, f"r{k}:{s}:", s, settings[s])
            rec = calibrate(
                runs + extra, gates=passing, previous=rec, repilot=True, n_draws=n_draws, seed=k
            )
            branch = "GO after re-pilot" if rec["decision"] == "GO" else "NO-GO"
        branches[branch] += 1
        if rec["decision"] != "GO":
            continue
        settings = rec["frozen_knobs"]
        main: dict[str, float] = {}
        for s in HARD_STRATA_R4:
            items = _sim_items(rng, mu[(s, settings[s])], MAIN_QUOTA, "m", s, settings[s])
            main[s] = mean([float(r.correct) for r in items])
        realised[branch].append(float(_primary(main)))
        true_b.append(float(_primary({s: theta[(s, settings[s])] for s in HARD_STRATA_R4})))
    go = [x for v in realised.values() for x in v]
    return {
        "n_sims": n_sims,
        "seed": seed,
        "prior": f"stratum default B ~ U(expected ± {SIM_PRIOR_HALF_WIDTH}); alt shifted "
        f"{SIM_ALT_SHIFT} in SIM_ALT_DIRECTION; models {SIM_MODELS}; item SD {SIM_ITEM_SD:.2f}",
        "branches": {b: branches.get(b, 0) / n_sims for b in _BRANCHES},
        "realised_main_primary": {b: _describe(realised[b]) for b in _BRANCHES if realised.get(b)},
        "realised_main_primary_all_go": _describe(go),
        "true_main_primary_all_go": _describe(true_b),
        "p_realised_in_g1_window_given_go": mean(
            [float(G1_WINDOW[0] <= x <= G1_WINDOW[1]) for x in go]
        )
        if go
        else None,
    }


_BRANCHES = ("GO", "GO after knob move", "GO after re-pilot", "NO-GO")


def _describe(xs: list[float]) -> dict[str, float] | None:
    if not xs:
        return None
    a = np.asarray(xs)
    return {
        "n": float(len(a)),
        "mean": float(a.mean()),
        "q05": float(np.quantile(a, 0.05)),
        "q50": float(np.quantile(a, 0.5)),
        "q95": float(np.quantile(a, 0.95)),
    }


def format_operating_characteristic(oc: dict[str, Any]) -> str:
    lines = [
        f"calibrate --simulate: {oc['n_sims']} simulated pilots (seed {oc['seed']})",
        f"prior: {oc['prior']}",
        "branch frequencies:",
    ]
    for b, f in oc["branches"].items():
        lines.append(f"  {b:20s} {f:.3f}")
    lines.append("realised main-set primary B (GO outcomes):")
    for b, d in oc["realised_main_primary"].items():
        lines.append(
            f"  {b:20s} n={int(d['n'])} mean={d['mean']:.3f} "
            f"[q05 {d['q05']:.3f}, q50 {d['q50']:.3f}, q95 {d['q95']:.3f}]"
        )
    for name, key in (
        ("all GO", "realised_main_primary_all_go"),
        ("true B, all GO", "true_main_primary_all_go"),
    ):
        d = oc[key]
        if d:
            lines.append(
                f"  {name:20s} n={int(d['n'])} mean={d['mean']:.3f} "
                f"[q05 {d['q05']:.3f}, q95 {d['q95']:.3f}]"
            )
    p = oc["p_realised_in_g1_window_given_go"]
    if p is not None:
        lines.append(f"P(realised main B in {G1_WINDOW} | GO) = {p:.3f}")
    return "\n".join(lines)

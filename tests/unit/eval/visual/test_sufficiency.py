"""G8 recoverability: numpy logistic regression, stratified CV AUROC, bands, proxy."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from bamcp.eval.visual.cli_r4 import main as cli_main
from bamcp.eval.visual.grade import STRATUM_PAIRS_R4
from bamcp.eval.visual.sufficiency import (
    FEATURE_KEYS,
    _prepare,
    adjudication_proxy,
    cv_predictions,
    feature_matrix,
    fit_logistic,
    predict_logistic,
    stratified_folds,
    sufficiency,
)

STRATUM = {(p, pol): s for s, (p, pol) in STRATUM_PAIRS_R4.items() if p != "anchor"}


@pytest.fixture(autouse=True)
def _whole_pool_in_frame(request, monkeypatch):
    """The synthetic rows below carry features, not the predicate fields of real candidates;
    every test except the frame test treats them as inside the sampling frame."""
    if "frame" not in request.node.name:
        from bamcp.eval.visual import strata_r4

        monkeypatch.setattr(strata_r4, "eligible", lambda row, level, overrides=None: (True, None))


def _row(rng: random.Random, pair: str, pol: str, sep: float, partition: str) -> dict[str, Any]:
    """A candidate: ``vaf`` and ``alt_mismatch_burden`` separate the classes by ``sep`` SDs."""
    pos = pol == "pos"
    alt = rng.choice([4, 5, 6, 7, 8, 9, 12])
    feats = {k: rng.gauss(0, 1) for k in FEATURE_KEYS}
    feats["vaf"] = rng.gauss(sep if pos else 0.0, 1)
    feats["alt_mismatch_burden"] = rng.gauss(0.0 if pos else sep, 1)
    feats["alt_reads"] = alt
    del feats["stutter_frac"]  # absent features are imputed
    return {
        "chrom": "chr20",
        "site0": rng.randrange(10**7),
        "stratum": STRATUM[(pair, pol)],
        "pair": pair,
        "polarity": pol,
        "partition": partition,
        "alt_reads": alt,
        "features": feats,
    }


def make_pool(sep: dict[str, float], n: int = 60, seed: int = 0) -> list[dict[str, Any]]:
    rng = random.Random(seed)  # noqa: S311 — reproducible synthetic pool
    pool = []
    for pair, s in sep.items():
        for pol in ("pos", "neg"):
            for i in range(n):
                part = "pilot" if i % 4 else "main"
                pool.append(_row(rng, pair, pol, s, part))
    return pool


def test_logistic_recovers_a_separating_direction() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(400, 3))
    y = (x[:, 0] - 0.5 * x[:, 1] + rng.normal(scale=0.3, size=400) > 0).astype(float)
    xtr, _ = _prepare(x, x)
    w = fit_logistic(xtr, y, l2=0.1)
    assert w[0] > 0 and w[1] < 0 and abs(w[2]) < abs(w[0]) / 3
    p = predict_logistic(w, xtr)
    assert ((p > 0.5) == (y == 1)).mean() > 0.9
    # L2 shrinks the slopes, not the intercept.
    w_big = fit_logistic(xtr, y, l2=1e4)
    assert np.abs(w_big[:3]).max() < 0.1


def test_perfect_separation_stays_finite() -> None:
    x = np.array([[0.0], [1.0], [2.0], [3.0]])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    xtr, _ = _prepare(x, x)
    w = fit_logistic(xtr, y, l2=1.0)
    assert np.all(np.isfinite(w))


def test_stratified_folds_balance_classes() -> None:
    y = np.array([1] * 23 + [0] * 17)
    f = stratified_folds(y, 5, seed=1)
    for k in range(5):
        assert abs(int((y[f == k] == 1).sum()) - 23 / 5) <= 1
        assert abs(int((y[f == k] == 0).sum()) - 17 / 5) <= 1
    assert np.array_equal(f, stratified_folds(y, 5, seed=1))


def test_feature_matrix_and_imputation() -> None:
    rows = [{"features": {"vaf": 0.1}}, {"features": {"vaf": None, "depth": 30}}, {}]
    x = feature_matrix(rows, ("vaf", "depth"))
    assert x.shape == (3, 2) and np.isnan(x[1, 0]) and np.isnan(x[2, 1])
    tr, te = _prepare(x, x)
    assert np.all(np.isfinite(tr)) and tr.shape == (3, 3)  # + intercept
    all_nan = feature_matrix([{}, {}], ("vaf",))
    assert np.all(np.isfinite(_prepare(all_nan, all_nan)[0]))


def test_cv_auroc_separable_vs_noise() -> None:
    pool = make_pool({"mosaic": 3.0, "tr": 3.0, "segdup": 0.0, "indel_adj": 3.0}, n=60)
    res = sufficiency(pool, seed=0)
    assert res["pairs"]["tr"]["auroc"] > 0.95 and res["pairs"]["tr"]["pass"]
    assert res["pairs"]["segdup"]["auroc"] < 0.75 and not res["pairs"]["segdup"]["pass"]
    assert res["g8"]["pass"] is False and res["g8"]["failed_pairs"] == ["segdup"]
    # Pilot partition only, except the mosaic pair, which uses every partition.
    assert res["pairs"]["tr"]["n_pos"] == 45 and res["pairs"]["mosaic"]["n_pos"] == 60
    bands = res["mosaic_bands"]
    assert set(bands) == {"4-6", ">=7"} and not any(b["drop"] for b in bands.values())
    assert bands["4-6"]["n_pos"] + bands[">=7"]["n_pos"] == 60
    oof = cv_predictions(
        feature_matrix(pool[:120]), np.array([r["polarity"] == "pos" for r in pool[:120]], float)
    )
    assert oof.shape == (120,) and np.all((oof >= 0) & (oof <= 1))


def test_band_drop_and_small_pairs() -> None:
    rng = random.Random(3)  # noqa: S311
    pool = make_pool({"mosaic": 3.0, "tr": 3.0, "segdup": 3.0, "indel_adj": 3.0}, n=60)
    # Make the 4-6 band undecidable: its positives look like negatives.
    for r in pool:
        if r["pair"] == "mosaic" and r["polarity"] == "pos" and r["alt_reads"] <= 6:
            r["features"]["vaf"] = rng.gauss(0, 1)
            r["features"]["alt_mismatch_burden"] = rng.gauss(3.0, 1)
    res = sufficiency(pool)
    assert res["mosaic_bands"]["4-6"]["drop"] and not res["mosaic_bands"][">=7"]["drop"]
    assert res["g8"]["drop_bands"] == ["4-6"]
    tiny = sufficiency(make_pool({"tr": 3.0}, n=4))
    assert tiny["pairs"]["tr"]["auroc"] is None and "fewer" in tiny["pairs"]["tr"]["detail"]
    assert tiny["pairs"]["mosaic"]["n_pos"] == 0 and not tiny["g8"]["pass"]


def test_adjudication_proxy() -> None:
    pool = make_pool({"mosaic": 3.0, "tr": 3.0, "segdup": 3.0, "indel_adj": 3.0}, n=40)
    items = [
        {
            "chrom": "c",
            "site0": 1,
            "pair": "tr",
            "features": {"vaf": 4.0, "alt_mismatch_burden": -1.0},
        },
        {
            "chrom": "c",
            "site0": 2,
            "pair": "tr",
            "features": {"vaf": -1.0, "alt_mismatch_burden": 4.0},
        },
        {
            "chrom": "c",
            "site0": 3,
            "pair": "tr",
            "features": {"vaf": 1.5, "alt_mismatch_burden": 1.5},
        },
        {"chrom": "c", "site0": 4, "pair": "anchor", "features": {}},
    ]
    out = adjudication_proxy(pool, items)
    assert out[0]["posterior_real"] > 0.9 and out[1]["posterior_real"] < 0.1
    assert 0.2 < out[2]["posterior_real"] < 0.8
    assert out[3]["posterior_real"] is None and out[3]["undecidable"] is False
    assert out[2]["undecidable"] == (0.4 <= out[2]["posterior_real"] <= 0.6)


def test_cli_sufficiency(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    pool = make_pool({"mosaic": 3.0, "tr": 3.0, "segdup": 3.0, "indel_adj": 3.0}, n=40)
    pool_path = tmp_path / "pool.jsonl"
    pool_path.write_text("".join(json.dumps(r) + "\n" for r in pool))
    out = tmp_path / "g8.json"
    rc = cli_main(
        ["sufficiency", "--pool", str(pool_path), "--out", str(out), "--items", str(pool_path)]
    )
    assert rc == 0
    res = json.loads(out.read_text())
    assert res["g8"]["pass"] is True
    proxy = tmp_path / "g8.proxy.jsonl"
    assert len(proxy.read_text().splitlines()) == len(pool)
    assert "G8 pass=True" in capsys.readouterr().out
    assert (
        cli_main(["sufficiency", "--pool", str(pool_path), "--out", str(out), "--partition", "all"])
        == 0
    )
    assert json.loads(out.read_text())["pairs"]["tr"]["n_pos"] == 40


def test_adjudication_proxy_uses_out_of_fold_for_pool_items() -> None:
    """Mosaic main items sit in the (all-partition) training pool: their posterior must be
    the out-of-fold prediction, not an in-sample fit that overstates certainty."""
    pool = make_pool({"mosaic": 0.0}, n=30, seed=3)  # pure noise, 25 features
    rows = [r for r in pool if r["pair"] == "mosaic"]
    items = [r for r in rows if r["partition"] == "main"]
    y = np.array([r["polarity"] == "pos" for r in rows], dtype=float)
    oof = cv_predictions(feature_matrix(rows), y)
    expect = {(r["chrom"], r["site0"]): float(v) for r, v in zip(rows, oof, strict=True)}
    got = adjudication_proxy(pool, items)
    for it, g in zip(items, got, strict=True):
        assert g["posterior_real"] == pytest.approx(expect[(it["chrom"], it["site0"])])
    # On noise an in-sample fit "knows" each item's own label; out-of-fold it does not.
    xtr, _ = _prepare(feature_matrix(rows), feature_matrix(rows))
    insample = predict_logistic(fit_logistic(xtr, y), xtr)
    own_in = np.mean(np.where(y == 1, insample, 1 - insample))
    own_oof = np.mean(np.where(y == 1, oof, 1 - oof))
    assert own_in > 0.7 and own_oof < 0.6


def test_excluded_candidates_are_not_labelled_examples() -> None:
    from bamcp.eval.visual.sufficiency import pair_pool

    rng = random.Random(0)  # noqa: S311
    keep = _row(rng, "tr", "neg", 1.0, "pilot")
    dropped = {**_row(rng, "tr", "neg", 1.0, "pilot"), "exclude": "truth_in_run_30bp"}
    assert pair_pool([keep, dropped], "tr") == [keep]


def test_pool_keeps_only_the_sampling_frame() -> None:
    from bamcp.eval.visual.sufficiency import pair_pool

    base = {"pair": "mosaic", "polarity": "neg", "stratum": "mosaic_caller_fp_neg"}
    inside = {**base, "features": {"alt_reads": 6, "vaf": 0.12, "depth": 50}}
    outside = {**base, "features": {"alt_reads": 20, "vaf": 0.45, "depth": 50}}  # VAF > 0.20
    assert pair_pool([inside, outside], "mosaic", partition=None) == [inside]

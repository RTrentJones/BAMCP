"""G8 recoverability gate: is the pos/neg distinction recoverable from read-level features?

Design section 7. For each hard pair, on the whole labelled candidate pool (pilot partition;
for the mosaic pair, every partition — the mosaic pool is too small otherwise), fit an
L2-penalised logistic regression on the selection step's read-level ``features`` and score it
by 5-fold stratified cross-validated AUROC (out-of-fold predictions pooled). Requirements:

- pair AUROC >= 0.80, else the pair is redesigned (alternative knob) or replaced by its
  reserve pair;
- within stratum 1 (``mosaic_snv_pos``), any alt-count band whose AUROC (that band's
  positives against the negatives in the same alt-count band) is < 0.70 is dropped from the
  sampling frame.

The same per-pair model, fitted on the whole pool, gives the adjudication proxy (section 6):
an item whose posterior P(real) lies in 0.4-0.6 is "undecidable".

Model-free: no LLM output is read. numpy only (no scipy / sklearn).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .grade import PAIRS_R4
from .stats import auroc

FEATURE_KEYS = (
    "depth",
    "alt_reads",
    "vaf",
    "alt_fwd",
    "alt_rev",
    "ref_fwd",
    "ref_rev",
    "alt_mean_bq",
    "ref_mean_bq",
    "alt_mapq_lt20_frac",
    "ref_mapq_lt20_frac",
    "alt_near_end_frac",
    "ref_near_end_frac",
    "alt_softclip_frac",
    "ref_softclip_frac",
    "alt_mismatch_burden",
    "ref_mismatch_burden",
    "alt_proper_pair_frac",
    "ref_proper_pair_frac",
    "alt_mate_mapq_median",
    "ref_mate_mapq_median",
    "linked_mismatch_max_frac",
    "local_nonref_sites_100bp",
    "run_length",
    "stutter_frac",
)
AUROC_PAIR_MIN = 0.80
AUROC_BAND_MIN = 0.70
N_FOLDS = 5
L2 = 1.0
ALL_PARTITION_PAIRS = ("mosaic",)
# Stratum 1's alt-count bands (design 3.1): inclusive (lo, hi); hi None = open.
MOSAIC_BANDS: dict[str, tuple[int, int | None]] = {"4-6": (4, 6), ">=7": (7, None)}
UNDECIDABLE = (0.4, 0.6)


def load_jsonl(path: str | Path) -> list[dict[str, Any]]:
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]


def _feature(row: dict[str, Any], key: str) -> float:
    v = (row.get("features") or {}).get(key)
    if v is None or isinstance(v, bool):
        return math.nan
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def feature_matrix(rows: list[dict[str, Any]], keys: tuple[str, ...] = FEATURE_KEYS) -> np.ndarray:
    """rows × keys, NaN where a feature is absent."""
    return np.array([[_feature(r, k) for k in keys] for r in rows], dtype=float).reshape(
        len(rows), len(keys)
    )


def _prepare(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median-impute and standardise with training-fold statistics only; add an intercept."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(np.where(np.isnan(train).all(0), 0.0, train), axis=0)
    med = np.where(np.isnan(med), 0.0, med)

    def fill(x: np.ndarray) -> np.ndarray:
        return np.where(np.isnan(x), med, x)

    tr, te = fill(train), fill(test)
    mu, sd = tr.mean(0), tr.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    tr, te = (tr - mu) / sd, (te - mu) / sd
    return np.hstack([tr, np.ones((len(tr), 1))]), np.hstack([te, np.ones((len(te), 1))])


def fit_logistic(
    x: np.ndarray, y: np.ndarray, *, l2: float = L2, n_iter: int = 100, tol: float = 1e-9
) -> np.ndarray:
    """Newton-Raphson for L2-penalised logistic regression; the last column (intercept) is
    not penalised."""
    w = np.zeros(x.shape[1])
    pen = np.full(x.shape[1], l2)
    pen[-1] = 0.0
    for _ in range(n_iter):
        p = 1 / (1 + np.exp(-np.clip(x @ w, -30, 30)))
        grad = x.T @ (p - y) + pen * w
        h = (x * (p * (1 - p))[:, None]).T @ x + np.diag(pen) + 1e-9 * np.eye(x.shape[1])
        step = np.linalg.solve(h, grad)
        w -= step
        if float(np.max(np.abs(step))) < tol:
            break
    return w


def predict_logistic(w: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.asarray(1 / (1 + np.exp(-np.clip(x @ w, -30, 30))))


def stratified_folds(y: np.ndarray, k: int = N_FOLDS, seed: int = 0) -> np.ndarray:
    """Fold index per row: each class shuffled (seeded) and dealt round-robin."""
    rng = np.random.default_rng(seed)
    folds = np.zeros(len(y), dtype=int)
    for cls in (0, 1):
        idx = np.flatnonzero(y == cls)
        rng.shuffle(idx)
        folds[idx] = np.arange(len(idx)) % k
    return folds


def cv_predictions(
    x: np.ndarray, y: np.ndarray, *, k: int = N_FOLDS, seed: int = 0, l2: float = L2
) -> np.ndarray:
    """Out-of-fold P(positive) for every row."""
    folds = stratified_folds(y, k, seed)
    out = np.zeros(len(y))
    for f in range(k):
        te = folds == f
        if not te.any():
            continue
        xtr, xte = _prepare(x[~te], x[te])
        out[te] = predict_logistic(fit_logistic(xtr, y[~te], l2=l2), xte)
    return out


def _alt_reads(row: dict[str, Any]) -> float:
    v = row.get("alt_reads")
    return float(v) if v is not None else _feature(row, "alt_reads")


def _in_band(n: float, band: tuple[int, int | None]) -> bool:
    lo, hi = band
    return n >= lo and (hi is None or n <= hi)


def pair_pool(
    pool: list[dict[str, Any]], pair: str, *, partition: str | None = "pilot"
) -> list[dict[str, Any]]:
    """The pool rows G8 fits for ``pair`` (all partitions for the mosaic pair). Rows the
    selection excluded (``exclude``: e.g. a truth record in the run, so not a labelled
    negative; often no features) are not labelled examples and are left out, and so are rows
    outside the sampling frame (failing their stratum predicate at every knob level)."""
    from .strata_r4 import KNOB_LEVELS, eligible

    rows = [
        r
        for r in pool
        if r.get("pair") == pair and r.get("polarity") in ("pos", "neg") and not r.get("exclude")
    ]
    # Design 7: the sampling frame is the candidates that pass their stratum's predicate at a
    # piloted knob level; out-of-frame rows would measure a different (often easier) population.
    rows = [r for r in rows if any(eligible(r, lv)[0] for lv in KNOB_LEVELS)]
    if partition is not None and pair not in ALL_PARTITION_PAIRS:
        rows = [r for r in rows if r.get("partition") == partition]
    return rows


def sufficiency(
    pool: list[dict[str, Any]],
    *,
    partition: str | None = "pilot",
    feature_keys: tuple[str, ...] = FEATURE_KEYS,
    l2: float = L2,
    k: int = N_FOLDS,
    seed: int = 0,
) -> dict[str, Any]:
    """G8 per pair and per stratum-1 alt-count band → the gate record."""
    pairs: dict[str, Any] = {}
    bands: dict[str, Any] = {}
    for pair in PAIRS_R4:
        rows = pair_pool(pool, pair, partition=partition)
        y = np.array([r["polarity"] == "pos" for r in rows], dtype=float)
        n_pos, n_neg = int(y.sum()), int(len(y) - y.sum())
        res: dict[str, Any] = {"n_pos": n_pos, "n_neg": n_neg, "auroc": None, "pass": False}
        if min(n_pos, n_neg) < k:
            res["detail"] = f"fewer than {k} rows in a class"
            pairs[pair] = res
            continue
        x = feature_matrix(rows, feature_keys)
        oof = cv_predictions(x, y, k=k, seed=seed, l2=l2)
        a = auroc(oof.tolist(), [bool(v) for v in y])
        res.update(auroc=a, **{"pass": a is not None and a >= AUROC_PAIR_MIN})
        pairs[pair] = res
        if pair == "mosaic":
            for name, band in MOSAIC_BANDS.items():
                sel = [i for i, r in enumerate(rows) if _in_band(_alt_reads(r), band)]
                bp = [bool(y[i]) for i in sel]
                ba = auroc([float(oof[i]) for i in sel], bp) if sel else None
                bands[name] = {
                    "n_pos": sum(bp),
                    "n_neg": len(bp) - sum(bp),
                    "auroc": ba,
                    "drop": ba is None or ba < AUROC_BAND_MIN,
                }
    failed = [p for p, r in pairs.items() if not r["pass"]]
    return {
        "partition": partition,
        "feature_keys": list(feature_keys),
        "l2": l2,
        "folds": k,
        "seed": seed,
        "thresholds": {"pair_auroc_min": AUROC_PAIR_MIN, "band_auroc_min": AUROC_BAND_MIN},
        "pairs": pairs,
        "mosaic_bands": bands,
        "g8": {
            "pass": not failed,
            "failed_pairs": failed,
            "drop_bands": [b for b, r in bands.items() if r["drop"]],
        },
    }


def _site_key(row: dict[str, Any]) -> tuple[Any, Any]:
    return row.get("chrom"), row.get("site0")


def adjudication_proxy(
    pool: list[dict[str, Any]],
    items: list[dict[str, Any]],
    *,
    partition: str | None = "pilot",
    feature_keys: tuple[str, ...] = FEATURE_KEYS,
    l2: float = L2,
) -> list[dict[str, Any]]:
    """Posterior P(real) for each item from its pair's model fitted on the whole pool.

    Items outside the hard pairs (anchors), or whose pair cannot be fitted, get ``None``.
    An item that is itself a row of the pair's training pool (the mosaic pair pools every
    partition, so its main items are in it) gets its out-of-fold prediction instead: an
    in-sample posterior is pushed toward 0/1 and would hide undecidable items.
    """
    models: dict[str, tuple[np.ndarray, np.ndarray, dict[tuple[Any, Any], float]] | None] = {}
    for pair in PAIRS_R4:
        rows = pair_pool(pool, pair, partition=partition)
        y = np.array([r["polarity"] == "pos" for r in rows], dtype=float)
        if len(rows) == 0 or y.min() == y.max():
            models[pair] = None
            continue
        x = feature_matrix(rows, feature_keys)
        oof: dict[tuple[Any, Any], float] = {}
        if min(int(y.sum()), int(len(y) - y.sum())) >= N_FOLDS:
            pred = cv_predictions(x, y, l2=l2)
            oof = {_site_key(r): float(v) for r, v in zip(rows, pred, strict=True)}
        models[pair] = (x, fit_logistic(_prepare(x, x)[0], y, l2=l2), oof)
    out = []
    for it in items:
        fitted = models.get(str(it.get("pair")))
        post: float | None = None
        if fitted is not None:
            train, w, oof = fitted
            if _site_key(it) in oof:
                post = oof[_site_key(it)]
            else:
                _, xte = _prepare(train, feature_matrix([it], feature_keys))
                post = float(predict_logistic(w, xte)[0])
        out.append(
            {
                "chrom": it.get("chrom"),
                "site0": it.get("site0"),
                "stratum": it.get("stratum"),
                "pair": it.get("pair"),
                "posterior_real": post,
                "undecidable": post is not None and UNDECIDABLE[0] <= post <= UNDECIDABLE[1],
            }
        )
    return out

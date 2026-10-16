"""Statistics for the visual-vs-text eval: paired cluster bootstrap, McNemar, Brier, ECE."""

from __future__ import annotations

import math
import random
from collections import defaultdict
from collections.abc import Iterable, Sequence

N_BOOT = 2000


def mean(xs: Sequence[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _percentile(sorted_xs: list[float], q: float) -> float:
    if not sorted_xs:
        return float("nan")
    k = (len(sorted_xs) - 1) * q
    lo, hi = math.floor(k), math.ceil(k)
    return sorted_xs[lo] + (sorted_xs[hi] - sorted_xs[lo]) * (k - lo)


def cluster_bootstrap_mean(
    clusters: dict[str, list[float]], *, n_boot: int = N_BOOT, seed: int = 0
) -> tuple[float, float, float]:
    """Mean over clusters (each cluster = mean of its repeats) with a 95% percentile CI."""
    keys = sorted(clusters)
    vals = [mean(clusters[k]) for k in keys if clusters[k]]
    if not vals:
        return float("nan"), float("nan"), float("nan")
    rng = random.Random(seed)  # noqa: S311 — reproducible resampling
    boots = sorted(mean([rng.choice(vals) for _ in vals]) for _ in range(n_boot))
    return mean(vals), _percentile(boots, 0.025), _percentile(boots, 0.975)


def paired_cluster_bootstrap(
    a: dict[str, list[float]],
    b: dict[str, list[float]],
    *,
    n_boot: int = N_BOOT,
    seed: int = 0,
) -> tuple[float, float, float, int]:
    """Mean paired difference (A - B) over items present in both, with a 95% CI.

    Items are the resampling unit (clusters); repeats are averaged within an item first.
    Returns (delta, lo, hi, n_items).
    """
    keys = sorted(set(a) & set(b))
    diffs = [mean(a[k]) - mean(b[k]) for k in keys if a[k] and b[k]]
    if not diffs:
        return float("nan"), float("nan"), float("nan"), 0
    rng = random.Random(seed)  # noqa: S311 — reproducible resampling
    boots = sorted(mean([rng.choice(diffs) for _ in diffs]) for _ in range(n_boot))
    return mean(diffs), _percentile(boots, 0.025), _percentile(boots, 0.975), len(diffs)


def mcnemar_exact(pairs: Iterable[tuple[bool, bool]]) -> tuple[int, int, float]:
    """Exact two-sided McNemar on matched (A correct, B correct) pairs → (b, c, p)."""
    b = c = 0
    for x, y in pairs:
        if x and not y:
            b += 1
        elif y and not x:
            c += 1
    n = b + c
    if n == 0:
        return b, c, 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2**n
    return b, c, min(1.0, 2 * tail)


def brier(conf: Sequence[float], correct: Sequence[bool]) -> float:
    if not conf:
        return float("nan")
    return mean([(p - float(y)) ** 2 for p, y in zip(conf, correct, strict=True)])


def reliability_bins(
    conf: Sequence[float], correct: Sequence[bool], n_bins: int = 10
) -> list[tuple[float, float, int]]:
    """Per-bin (mean confidence, accuracy, count) for non-empty bins."""
    bins: dict[int, list[tuple[float, bool]]] = defaultdict(list)
    for p, y in zip(conf, correct, strict=True):
        bins[min(n_bins - 1, int(p * n_bins))].append((p, y))
    out = []
    for i in sorted(bins):
        rows = bins[i]
        out.append((mean([p for p, _ in rows]), mean([float(y) for _, y in rows]), len(rows)))
    return out


def ece(conf: Sequence[float], correct: Sequence[bool], n_bins: int = 10) -> float:
    if not conf:
        return float("nan")
    total = len(conf)
    return sum(
        cnt / total * abs(acc - c) for c, acc, cnt in reliability_bins(conf, correct, n_bins)
    )


def f1_scores(truth: Sequence[str], pred: Sequence[str]) -> dict[str, float]:
    """Per-class F1 over the classes present in ``truth``."""
    out: dict[str, float] = {}
    for cls in sorted(set(truth)):
        tp = sum(1 for t, p in zip(truth, pred, strict=True) if t == cls and p == cls)
        fp = sum(1 for t, p in zip(truth, pred, strict=True) if t != cls and p == cls)
        fn = sum(1 for t, p in zip(truth, pred, strict=True) if t == cls and p != cls)
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        out[cls] = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return out

"""Statistics for the visual-vs-text eval: paired cluster bootstrap, McNemar, Brier, ECE.

Round 4 (``docs/eval/visual-design-r4.md`` section 8) adds the cluster-level sign-flip
permutation test with its inverted CI, Holm, TOST on the same machinery, AUROC, AURC,
selective accuracy and beta quantiles (for the calibration's shrunken intervals).
"""

from __future__ import annotations

import itertools
import math
import random
from collections import defaultdict
from collections.abc import Iterable, Sequence

import numpy as np

N_BOOT = 2000
N_FLIPS = 10_000
ECE_BINS_R4 = 5
SELECTIVE_COVERAGES = (0.5, 0.8)
_EPS = 1e-12  # ties in the permutation distribution count as "at least as extreme"


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


def ece5(conf: Sequence[float], correct: Sequence[bool]) -> float:
    """The round-4 calibration metric: ECE over 5 equal-width bins."""
    return ece(conf, correct, n_bins=ECE_BINS_R4)


# ---------------------------------------------------------------------------
# Round 4: cluster-level sign-flip permutation test (design section 8)
# ---------------------------------------------------------------------------


def _cluster_sums(values: Sequence[float], clusters: Sequence[str]) -> np.ndarray:
    """Per-cluster sums of ``values`` (clusters in sorted order: deterministic flips)."""
    if len(values) != len(clusters):
        raise ValueError("values and clusters differ in length")
    sums: dict[str, float] = defaultdict(float)
    for v, c in zip(values, clusters, strict=True):
        sums[str(c)] += float(v)
    return np.array([sums[k] for k in sorted(sums)], dtype=float)


def sign_matrix(
    n_clusters: int, *, n_flips: int = N_FLIPS, seed: int = 0
) -> tuple[np.ndarray, bool]:
    """Sign vectors for the permutation distribution → (signs[n, k], exact).

    Exhaustive (all ``2**k`` vectors, the identity included) when that is no more than
    ``n_flips``; otherwise ``n_flips`` seeded random vectors.
    """
    if n_clusters <= 0:
        return np.ones((1, 0)), True
    if 2**n_clusters <= n_flips:
        rows = list(itertools.product((1.0, -1.0), repeat=n_clusters))
        return np.array(rows, dtype=float), True
    rng = np.random.default_rng(seed)
    return rng.choice(np.array([-1.0, 1.0]), size=(n_flips, n_clusters)), False


def _signflip_p(contribs: np.ndarray, signs: np.ndarray, exact: bool, alternative: str) -> float:
    obs = float(contribs.sum())
    null = signs @ contribs
    if alternative == "two-sided":
        hits = int(np.sum(np.abs(null) >= abs(obs) - _EPS))
    elif alternative == "greater":
        hits = int(np.sum(null >= obs - _EPS))
    elif alternative == "less":
        hits = int(np.sum(null <= obs + _EPS))
    else:
        raise ValueError(f"unknown alternative {alternative!r}")
    if exact:
        return hits / len(null)
    return (hits + 1) / (len(null) + 1)  # Monte Carlo: the observed labelling counts once


def cluster_signflip_test(
    weighted_diffs: Sequence[float],
    clusters: Sequence[str],
    *,
    n_flips: int = N_FLIPS,
    seed: int = 0,
    alternative: str = "two-sided",
) -> float:
    """Paired sign-flip permutation p-value for H0: sum(weighted_diffs) = 0.

    ``weighted_diffs`` are per-item contributions (an item's model-averaged A - B difference
    times its balanced-accuracy weight); items sharing a cluster are flipped together.
    Two-sided by default (``|T*| >= |T|``); ``"greater"``/``"less"`` give one-sided p.
    """
    contribs = _cluster_sums(weighted_diffs, clusters)
    if contribs.size == 0:
        return float("nan")
    signs, exact = sign_matrix(contribs.size, n_flips=n_flips, seed=seed)
    return _signflip_p(contribs, signs, exact, alternative)


def _shifted(
    diffs: Sequence[float], weights: Sequence[float], clusters: Sequence[str], shift: float
) -> np.ndarray:
    return _cluster_sums([w * (d - shift) for d, w in zip(diffs, weights, strict=True)], clusters)


def _check_paired(
    diffs: Sequence[float], weights: Sequence[float], clusters: Sequence[str]
) -> None:
    if not len(diffs) == len(weights) == len(clusters):
        raise ValueError("diffs, weights and clusters differ in length")


def signflip_ci(
    diffs: Sequence[float],
    weights: Sequence[float],
    clusters: Sequence[str],
    *,
    level: float = 0.95,
    n_flips: int = N_FLIPS,
    seed: int = 0,
    tol: float = 1e-5,
) -> tuple[float, float]:
    """CI for the weighted mean difference by inverting the sign-flip test.

    The interval is every shift ``d0`` for which H0: Δ = ``d0`` — tested on the shifted
    item differences ``w_i (D_i - d0)`` with the same sign vectors for every shift — is not
    rejected at ``1 - level`` (two-sided). Each end is found by bisection; an end that is
    never rejected (too few clusters for the level) is ``-inf`` / ``inf``.
    """
    _check_paired(diffs, weights, clusters)
    if not diffs:
        return float("nan"), float("nan")
    alpha = 1 - level
    wsum = sum(weights)
    est = sum(w * d for d, w in zip(diffs, weights, strict=True)) / wsum
    n_clusters = len(set(map(str, clusters)))
    signs, exact = sign_matrix(n_clusters, n_flips=n_flips, seed=seed)

    def rejected(shift: float) -> bool:
        p = _signflip_p(_shifted(diffs, weights, clusters, shift), signs, exact, "two-sided")
        return p < alpha

    span = max(abs(d - est) for d in diffs) + 1.0

    def edge(direction: float) -> float:
        inner, outer = est, est + direction * span
        if not rejected(outer):
            return direction * float("inf")
        while abs(outer - inner) > tol:
            mid = (inner + outer) / 2
            if rejected(mid):
                outer = mid
            else:
                inner = mid
        return (inner + outer) / 2

    return edge(-1.0), edge(1.0)


def tost(
    diffs: Sequence[float],
    weights: Sequence[float],
    clusters: Sequence[str],
    *,
    margin: float = 0.07,
    alpha: float = 0.05,
    n_flips: int = N_FLIPS,
    seed: int = 0,
) -> dict[str, float | bool]:
    """Two one-sided sign-flip tests of equivalence within ``±margin`` (design section 8).

    Lower: H0 Δ <= -margin vs Δ > -margin, on ``w_i (D_i + margin)`` (one-sided, greater);
    upper: H0 Δ >= margin vs Δ < margin, on ``w_i (D_i - margin)`` (one-sided, less).
    ``equivalent`` iff both p < ``alpha`` (the same as the 90% inverted CI lying inside the
    margin). ``diffs`` are unweighted per-item differences; the estimand is their
    ``weights``-weighted mean.
    """
    _check_paired(diffs, weights, clusters)
    if not diffs:
        return {
            "p_lower": float("nan"),
            "p_upper": float("nan"),
            "p": float("nan"),
            "equivalent": False,
        }
    n_clusters = len(set(map(str, clusters)))
    signs, exact = sign_matrix(n_clusters, n_flips=n_flips, seed=seed)
    p_lo = _signflip_p(_shifted(diffs, weights, clusters, -margin), signs, exact, "greater")
    p_hi = _signflip_p(_shifted(diffs, weights, clusters, margin), signs, exact, "less")
    p = max(p_lo, p_hi)
    return {"p_lower": p_lo, "p_upper": p_hi, "p": p, "equivalent": p < alpha}


def design_effect(weighted_diffs: Sequence[float], clusters: Sequence[str]) -> float:
    """Variance of the sign-flip null with cluster flips over that with item flips.

    Under random signs the null variance is the sum of squared flip units, so this is
    ``sum_k C_k**2 / sum_i d_i**2`` (1 when every cluster holds one item).
    """
    item = sum(float(d) ** 2 for d in weighted_diffs)
    if item == 0:
        return float("nan")
    return float(np.sum(_cluster_sums(weighted_diffs, clusters) ** 2)) / item


def holm(pvals: Sequence[float]) -> list[float]:
    """Holm step-down adjusted p-values, in input order (monotone, capped at 1)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvals[i]))
        adj[i] = running
    return adj


def binom_sf(k: int, n: int, p: float = 0.5) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(max(0, k), n + 1))


# ---------------------------------------------------------------------------
# Round 4: discrimination and selective prediction
# ---------------------------------------------------------------------------


def auroc(scores: Sequence[float], labels: Sequence[bool]) -> float | None:
    """Mann-Whitney AUROC (ties count one half); ``None`` when a class is missing."""
    pos = [s for s, y in zip(scores, labels, strict=True) if y]
    neg = [s for s, y in zip(scores, labels, strict=True) if not y]
    if not pos or not neg:
        return None
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    r_pos = sum(r for r, y in zip(ranks, labels, strict=True) if y)
    n1, n0 = len(pos), len(neg)
    return (r_pos - n1 * (n1 + 1) / 2) / (n1 * n0)


def _tie_groups(conf: Sequence[float]) -> list[list[int]]:
    """Indices grouped by equal confidence, most confident group first."""
    groups: dict[float, list[int]] = defaultdict(list)
    for i, c in enumerate(conf):
        groups[float(c)].append(i)
    return [groups[c] for c in sorted(groups, reverse=True)]


def aurc(conf: Sequence[float], correct: Sequence[bool]) -> float:
    """Area under the risk-coverage curve (lower is better).

    Runs are admitted in order of decreasing confidence; tied confidences are admitted
    together, so the curve is evaluated at the end of each tie group and weighted by the
    group's size (ordering inside a tie is not information the answer carries).
    """
    n = len(conf)
    if n == 0:
        return float("nan")
    area = 0.0
    seen = errors = 0
    for g in _tie_groups(conf):
        seen += len(g)
        errors += sum(1 for i in g if not correct[i])
        area += len(g) / n * (errors / seen)
    return area


def selective_accuracy(
    conf: Sequence[float], correct: Sequence[bool], coverage: float
) -> tuple[float, float]:
    """Accuracy over the most confident ``coverage`` share → (accuracy, realised coverage).

    Whole tie groups are admitted until at least ``coverage`` of the runs are in, so the
    realised coverage can exceed the nominal one.
    """
    n = len(conf)
    if n == 0:
        return float("nan"), float("nan")
    want = max(1, math.ceil(coverage * n))
    kept: list[int] = []
    for g in _tie_groups(conf):
        kept.extend(g)
        if len(kept) >= want:
            break
    return mean([float(correct[i]) for i in kept]), len(kept) / n


def confident_wrong_rate(
    conf: Sequence[float], correct: Sequence[bool], threshold: float = 0.8
) -> float:
    """Share of runs answered with confidence >= ``threshold`` and wrong."""
    if not conf:
        return float("nan")
    return mean([float(c >= threshold and not y) for c, y in zip(conf, correct, strict=True)])


# ---------------------------------------------------------------------------
# Beta distribution quantiles (no scipy): the calibration's shrunken intervals
# ---------------------------------------------------------------------------


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the regularized incomplete beta (Numerical Recipes 6.4)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return h


def beta_cdf(x: float, a: float, b: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    ln = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x)
    ln += b * math.log1p(-x)
    front = math.exp(ln)
    if x < (a + 1) / (a + b + 2):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def beta_ppf(q: float, a: float, b: float, *, tol: float = 1e-10) -> float:
    """Quantile of Beta(a, b) by bisection on :func:`beta_cdf`."""
    lo, hi = 0.0, 1.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if beta_cdf(mid, a, b) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

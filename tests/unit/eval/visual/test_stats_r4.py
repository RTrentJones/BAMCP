"""Round-4 statistics: cluster sign-flip test, inverted CI, TOST, Holm, AUROC/AURC, beta."""

from __future__ import annotations

import itertools
import math

import pytest

from bamcp.eval.visual.stats import (
    _cluster_sums,
    _signflip_p,
    aurc,
    auroc,
    beta_cdf,
    beta_ppf,
    binom_sf,
    cluster_signflip_test,
    confident_wrong_rate,
    design_effect,
    ece,
    ece5,
    holm,
    selective_accuracy,
    sign_matrix,
    signflip_ci,
    tost,
)


def _exhaustive_p(contribs: list[float], alternative: str = "two-sided") -> float:
    """Brute-force permutation p over every sign vector (independent of the module)."""
    obs = sum(contribs)
    hits = 0
    total = 0
    for signs in itertools.product((1, -1), repeat=len(contribs)):
        t = sum(s * c for s, c in zip(signs, contribs, strict=True))
        total += 1
        if alternative == "two-sided":
            hits += abs(t) >= abs(obs) - 1e-12
        elif alternative == "greater":
            hits += t >= obs - 1e-12
        else:
            hits += t <= obs + 1e-12
    return hits / total


def test_signflip_exact_matches_enumeration() -> None:
    d = [0.3, 0.1, -0.05, 0.2, 0.15]
    cl = ["a", "b", "c", "d", "e"]
    p = cluster_signflip_test(d, cl)
    assert p == pytest.approx(_exhaustive_p(d))
    # Hand count: T = 0.70; |T*| >= 0.70 only for the identity, flipping -0.05 (T*=0.80)
    # and their negations → 4 / 32.
    assert p == pytest.approx(4 / 32)
    for alt in ("greater", "less"):
        assert cluster_signflip_test(d, cl, alternative=alt) == pytest.approx(_exhaustive_p(d, alt))


def test_signflip_flips_whole_clusters() -> None:
    d = [0.2, 0.2, 0.1, -0.1]
    p_items = cluster_signflip_test(d, ["a", "b", "c", "d"])
    p_cl = cluster_signflip_test(d, ["a", "a", "c", "d"])
    assert p_cl == pytest.approx(_exhaustive_p([0.4, 0.1, -0.1]))
    assert p_items == pytest.approx(_exhaustive_p(d))
    assert p_cl != p_items


def test_signflip_monte_carlo_close_to_exact() -> None:
    d = [0.05 * ((i * 7) % 5 - 1.5) + 0.03 for i in range(15)]
    cl = [str(i) for i in range(15)]  # 2**15 > 10,000 → Monte Carlo
    signs, exact = sign_matrix(15, n_flips=10_000, seed=3)
    assert not exact and signs.shape == (10_000, 15)
    p_mc = cluster_signflip_test(d, cl, seed=3)
    assert abs(p_mc - _exhaustive_p(d)) < 0.015
    # Deterministic under the seed; never 0 under Monte Carlo.
    assert cluster_signflip_test(d, cl, seed=3) == p_mc
    assert cluster_signflip_test([1.0] * 15, cl) == pytest.approx(1 / 10_001, rel=0.5)


def test_signflip_edge_cases() -> None:
    assert math.isnan(cluster_signflip_test([], []))
    assert cluster_signflip_test([0.0, 0.0], ["a", "b"]) == 1.0
    with pytest.raises(ValueError):
        cluster_signflip_test([1.0], [])
    with pytest.raises(ValueError):
        _signflip_p(_cluster_sums([1.0], ["a"]), *sign_matrix(1), alternative="sideways")


def test_inverted_ci_brackets_the_estimate_at_the_test_boundary() -> None:
    d = [0.4, 0.2, 0.1, 0.3, -0.1, 0.25, 0.35, 0.0, 0.15]
    w = [1 / len(d)] * len(d)
    cl = [f"c{i}" for i in range(len(d))]
    est = sum(x * y for x, y in zip(d, w, strict=True))
    lo, hi = signflip_ci(d, w, cl)
    assert lo < est < hi

    def p_at(shift: float) -> float:
        return _exhaustive_p([y * (x - shift) for x, y in zip(d, w, strict=True)])

    eps = 1e-3
    assert p_at(lo + eps) >= 0.05 > p_at(lo - eps)
    assert p_at(hi - eps) >= 0.05 > p_at(hi + eps)
    # The CI and the test agree on 0.
    assert (lo > 0) == (
        cluster_signflip_test([x * y for x, y in zip(d, w, strict=True)], cl) < 0.05
    )


def test_inverted_ci_unbounded_with_too_few_clusters() -> None:
    # 4 clusters: the smallest two-sided p is 2/16 > 0.05, so no shift is ever rejected.
    lo, hi = signflip_ci([0.1, 0.2, 0.3, 0.4], [0.25] * 4, list("abcd"))
    assert lo == -math.inf and hi == math.inf
    assert all(math.isnan(x) for x in signflip_ci([], [], []))


def test_tost_known_answers() -> None:
    cl = [f"c{i}" for i in range(8)]  # 256 sign vectors: exact
    w = [1 / 8] * 8
    near_zero = [0.01, -0.02, 0.0, 0.015, -0.01, 0.005, -0.005, 0.01]
    res = tost(near_zero, w, cl, margin=0.07)
    lower = _exhaustive_p([y * (x + 0.07) for x, y in zip(near_zero, w, strict=True)], "greater")
    upper = _exhaustive_p([y * (x - 0.07) for x, y in zip(near_zero, w, strict=True)], "less")
    assert res["p_lower"] == pytest.approx(lower) and res["p_upper"] == pytest.approx(upper)
    # Every shifted item is positive (lower) / negative (upper): only the identity is as
    # extreme → p = 1/256 on each side.
    assert res["p"] == pytest.approx(1 / 256) and res["equivalent"] is True

    big = [0.10, 0.12, 0.08, 0.09, 0.11, 0.10, 0.13, 0.07]
    res = tost(big, w, cl, margin=0.07)
    assert res["equivalent"] is False and res["p_upper"] > 0.5
    empty = tost([], [], [])
    assert empty["equivalent"] is False and math.isnan(float(empty["p"]))


def test_tost_uses_weighted_mean_scale() -> None:
    # Unnormalised weights: the estimand is the weighted mean, so the result is unchanged.
    cl = [f"c{i}" for i in range(8)]
    d = [0.01, -0.02, 0.0, 0.015, -0.01, 0.005, -0.005, 0.01]
    a = tost(d, [1 / 8] * 8, cl)
    b = tost(d, [1.0] * 8, cl)
    assert a["p"] == pytest.approx(b["p"])


def test_design_effect_and_holm_and_binom() -> None:
    assert design_effect([0.1, 0.1], ["a", "b"]) == pytest.approx(1.0)
    assert design_effect([0.1, 0.1], ["a", "a"]) == pytest.approx(2.0)
    assert math.isnan(design_effect([0.0], ["a"]))
    assert holm([0.01, 0.04]) == pytest.approx([0.02, 0.04])
    assert holm([0.04, 0.01]) == pytest.approx([0.04, 0.02])
    assert holm([0.03, 0.02]) == pytest.approx([0.04, 0.04])
    assert holm([0.6, 0.7]) == pytest.approx([1.0, 1.0])
    assert binom_sf(0, 5) == pytest.approx(1.0)
    assert binom_sf(5, 5) == pytest.approx(1 / 32)
    assert binom_sf(4, 7) == pytest.approx(sum(math.comb(7, i) for i in range(4, 8)) / 128)


def test_auroc_aurc_selective() -> None:
    assert auroc([0.9, 0.8, 0.2, 0.1], [True, True, False, False]) == 1.0
    assert auroc([0.1, 0.2, 0.8, 0.9], [True, True, False, False]) == 0.0
    assert auroc([0.5, 0.5], [True, False]) == 0.5
    assert auroc([0.4, 0.6, 0.6, 0.2], [False, True, False, True]) == pytest.approx(0.375)
    assert auroc([0.4], [True]) is None

    assert aurc([0.9, 0.8], [True, True]) == 0.0
    # Order 0.9 (wrong), 0.8 (right): risks 1, 1/2 → mean 0.75.
    assert aurc([0.9, 0.8], [False, True]) == pytest.approx(0.75)
    # A tie is admitted as one group: risk 1/2 at full coverage.
    assert aurc([0.5, 0.5], [False, True]) == pytest.approx(0.5)
    assert math.isnan(aurc([], []))

    conf = [0.9, 0.8, 0.7, 0.6]
    corr = [True, True, False, False]
    assert selective_accuracy(conf, corr, 0.5) == (1.0, 0.5)
    assert selective_accuracy(conf, corr, 1.0) == (0.5, 1.0)
    assert selective_accuracy([0.5] * 4, corr, 0.5) == (0.5, 1.0)  # tie → whole group
    assert math.isnan(selective_accuracy([], [], 0.5)[0])
    assert confident_wrong_rate([0.9, 0.85, 0.5], [True, False, False]) == pytest.approx(1 / 3)
    assert math.isnan(confident_wrong_rate([], []))


def test_ece5() -> None:
    conf = [0.15, 0.25, 0.95, 0.85]
    corr = [False, True, True, True]
    # 5 bins: [0, .2) {0.15: acc 0}, [.2, .4) {0.25: acc 1}, [.8, 1] {0.95, 0.85: acc 1}.
    want = (0.15 + 0.75 + 2 * 0.1) / 4
    assert ece5(conf, corr) == pytest.approx(want)
    assert ece5(conf, corr) == pytest.approx(ece(conf, corr, n_bins=5))


@pytest.mark.parametrize("q", [0.05, 0.5, 0.9])
def test_beta_quantiles(q: float) -> None:
    assert beta_ppf(q, 1, 1) == pytest.approx(q, abs=1e-8)
    assert beta_ppf(q, 2, 1) == pytest.approx(math.sqrt(q), abs=1e-8)
    assert beta_ppf(q, 1, 2) == pytest.approx(1 - math.sqrt(1 - q), abs=1e-8)
    a, b = 17.3, 6.2
    assert beta_cdf(beta_ppf(q, a, b), a, b) == pytest.approx(q, abs=1e-8)
    assert beta_cdf(0.5, 30, 30) == pytest.approx(0.5)
    assert beta_cdf(0.0, 2, 2) == 0.0 and beta_cdf(1.0, 2, 2) == 1.0

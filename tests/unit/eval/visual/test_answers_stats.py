"""Answer schema validation/normalization and the statistics helpers."""

from __future__ import annotations

import json
import math

import pytest

from bamcp.eval.visual.answers import (
    field_correct,
    normalize_allele,
    parse_answer,
    validate_answer,
)
from bamcp.eval.visual.stats import (
    brier,
    cluster_bootstrap_mean,
    ece,
    f1_scores,
    mcnemar_exact,
    paired_cluster_bootstrap,
    reliability_bins,
)

GOOD = {
    "variant_present": "yes",
    "variant_type": "SNV",
    "alt_allele": "T",
    "zygosity": "het",
    "artifact_verdict": "real",
    "artifact_type": "none",
    "phase": "na",
    "confidence": 0.8,
    "rationale": "Balanced alt reads on both orientations.",
}


def test_valid_answer_and_code_fence() -> None:
    assert parse_answer(json.dumps(GOOD)).ok
    assert parse_answer("```json\n" + json.dumps(GOOD) + "\n```").ok


@pytest.mark.parametrize(
    "mutate, msg",
    [
        (lambda a: a.pop("phase"), "missing"),
        (lambda a: a.update(zygosity="diploid"), "zygosity"),
        (lambda a: a.update(confidence=1.5), "confidence"),
        (lambda a: a.update(confidence=True), "confidence"),
        (lambda a: a.update(alt_allele=3), "alt_allele"),
        (lambda a: a.update(rationale=None), "rationale"),
        (lambda a: a.update(rationale="word " * 61), "rationale"),
    ],
)
def test_invalid_answers(mutate, msg) -> None:
    a = dict(GOOD)
    mutate(a)
    v = validate_answer(a)
    assert not v.ok and msg in (v.error or "")


def test_not_json_or_not_object() -> None:
    assert not parse_answer("hello").ok
    assert not parse_answer("[1, 2]").ok


def test_normalize_allele() -> None:
    assert normalize_allele(None) == ""
    assert normalize_allele("none") == ""
    assert normalize_allele("+ag") == "AG"
    assert normalize_allele("G, a") == "A,G"
    assert normalize_allele("C/T") == "C,T"
    assert field_correct("alt_allele", {"alt_allele": "-c"}, {"alt_allele": "C"})
    assert field_correct("zygosity", {"zygosity": "het"}, {"zygosity": "het"})
    assert not field_correct("zygosity", {"zygosity": "het"}, {"zygosity": "none"})


def test_mcnemar_exact() -> None:
    assert mcnemar_exact([(True, True), (False, False)]) == (0, 0, 1.0)
    b, c, p = mcnemar_exact([(True, False)] * 10)
    assert (b, c) == (10, 0) and p == pytest.approx(2 / 2**10)
    _, _, p_sym = mcnemar_exact([(True, False), (False, True)] * 5)
    assert p_sym == 1.0


def test_bootstraps() -> None:
    clusters = {"a": [1.0, 1.0], "b": [0.0, 0.0], "c": [1.0, 0.0]}
    m, lo, hi = cluster_bootstrap_mean(clusters, n_boot=500)
    assert m == pytest.approx(0.5) and 0 <= lo <= m <= hi <= 1
    d, lo, hi, n = paired_cluster_bootstrap(
        {"a": [1.0], "b": [1.0]}, {"a": [0.0], "b": [1.0]}, n_boot=200
    )
    assert d == pytest.approx(0.5) and n == 2 and lo <= d <= hi
    assert math.isnan(cluster_bootstrap_mean({})[0])
    assert paired_cluster_bootstrap({"a": [1.0]}, {"b": [1.0]})[3] == 0


def test_calibration_and_f1() -> None:
    conf = [0.9, 0.9, 0.1, 0.1]
    corr = [True, True, False, False]
    assert brier(conf, corr) == pytest.approx(0.01)
    assert ece(conf, corr) == pytest.approx(0.1)
    assert reliability_bins([1.0], [True]) == [(1.0, 1.0, 1)]
    assert math.isnan(brier([], [])) and math.isnan(ece([], []))
    f1 = f1_scores(["a", "a", "b"], ["a", "b", "b"])
    assert f1["a"] == pytest.approx(2 / 3) and f1["b"] == pytest.approx(2 / 3)
    assert f1_scores(["a"], ["z"])["a"] == 0.0

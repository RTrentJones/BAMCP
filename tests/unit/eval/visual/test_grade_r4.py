"""Round-4 grading: balanced primary, guardrail, utility, contrasts, alleles, self-tests."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from bamcp.eval.visual.grade import (
    HARD_STRATA_R4,
    PAIRS_R4,
    Graded,
    allele_equivalent,
    balanced_weights,
    directional_sign_test,
    genotype_correct,
    grade_run_set_r4,
    grade_runs_r4,
    pair_of,
    run_utility,
)
from bamcp.eval.visual.selftest import oracle_answer, run_selftests_r4

from .r4_analysis_world import ARTIFACT, REAL, make_r4_run_set, uncertain_answer

FLIPS = 2000


def _grade(rs: Path, **kw: Any) -> dict[str, Any]:
    return grade_run_set_r4(rs, n_flips=FLIPS, **kw)


def test_primary_is_balanced_and_model_weighted(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        return {"haiku": REAL, "sonnet": oracle_answer(lab), "opus": ARTIFACT}[m]

    s = _grade(make_r4_run_set(tmp_path, conditions=("Traw-r4",), answer_fn=fn))
    g = s["groups"]["__all__|Traw-r4"]
    assert s["groups"]["haiku|Traw-r4"]["primary"] == 0.5
    assert s["groups"]["sonnet|Traw-r4"]["primary"] == 1.0
    assert s["groups"]["opus|Traw-r4"]["primary"] == 0.5
    assert g["primary"] == pytest.approx(2 / 3)
    assert set(g["per_pair"]) == {*PAIRS_R4, "anchor"}
    # Sensitivity: haiku 1, sonnet 1, opus 0; specificity: 0, 1, 1.
    assert g["sensitivity"] == pytest.approx(2 / 3) and g["specificity"] == pytest.approx(2 / 3)
    assert g["per_model"] == {"haiku": 0.5, "opus": 0.5, "sonnet": 1.0}
    # Pooled item-weighted joint accuracy: every item is right for 2 of the 3 models.
    assert g["pooled_joint_accuracy"] == pytest.approx(2 / 3)
    tbl = s["tables"]["per_stratum_model"]["Traw-r4"]
    assert tbl["haiku"]["mosaic_snv_pos"] == 1.0 and tbl["opus"]["mosaic_snv_pos"] == 0.0


def test_anchors_do_not_move_the_primary(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        return ARTIFACT if lab["pair"] == "anchor" else oracle_answer(lab)

    g = _grade(make_r4_run_set(tmp_path, conditions=("Traw-r4",), answer_fn=fn))["groups"]
    assert g["__all__|Traw-r4"]["primary"] == 1.0
    assert g["__all__|Traw-r4"]["per_pair"]["anchor"]["pos"] == 0.0
    assert g["__all__|Traw-r4"]["pooled_joint_accuracy"] < 1.0


def test_uncertain_is_wrong_for_primary_and_zero_for_utility(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        return uncertain_answer() if lab["polarity"] == "neg" else oracle_answer(lab)

    g = _grade(make_r4_run_set(tmp_path, conditions=("Traw-r4",), answer_fn=fn))["groups"][
        "__all__|Traw-r4"
    ]
    assert g["primary"] == 0.5
    assert g["utility"] == pytest.approx(0.5)  # pos +1, neg 0
    assert g["utility_wrong_m2"] == pytest.approx(0.5)
    assert g["abstention_rate"] == pytest.approx(0.5, abs=0.1)


def test_run_utility_values() -> None:
    lab = {"variant_present": "yes", "artifact_verdict": "real", "artifact_type": "none"}
    run = {"item": "i", "model": "m", "condition": "c", "repeat": 0}
    ok = Graded(
        run, lab, "valid", answer=REAL, correct={"variant_present": True, "artifact_verdict": True}
    )
    unc = Graded(run, lab, "valid", answer=uncertain_answer(), correct={"variant_present": True})
    bad = Graded(run, lab, "invalid", "x")
    assert run_utility(ok) == 1.0 and run_utility(unc) == 0.0
    assert run_utility(bad) == -1.0 and run_utility(bad, -2.0) == -2.0


def test_contrasts_supported_holm_and_guardrail(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        if c == "Traw-r4":
            return REAL  # sens 1, spec 0: balanced 0.5 everywhere
        if c == "Traw-r4+Tcur":
            return oracle_answer(lab)  # +0.5, guardrail holds
        # C2: oracle except mosaic positives now missed: sensitivity drops 0.25
        if lab["stratum"] == "mosaic_snv_pos":
            return ARTIFACT
        return oracle_answer(lab)

    s = _grade(make_r4_run_set(tmp_path, n=6, answer_fn=fn))
    c1, c2 = s["contrasts"]["C1"], s["contrasts"]["C2"]
    assert c1["delta"] == pytest.approx(0.5)
    assert c1["p"] < 0.01 and c1["holm_p"] >= c1["p"]
    assert c1["ci95"][0] > 0 and c1["ci_excludes_0"]
    assert c1["guardrail"]["holds"] and c1["supported"] and c1["verdict"] == "supported"
    assert c1["n_items"] == 6 * 8 and c1["n_clusters"] == 3 * 8
    assert c1["design_effect"] == pytest.approx(2.0)  # 2 identical items per cluster
    assert c1["utility"]["delta"] == pytest.approx(1.0)  # neg -1 → +1, weight 1/2
    assert c1["tost"]["equivalent"] is False
    # C2: Δ = 3 pairs × +0.5 / 4 = 0.375, but sensitivity falls by 0.25 → not supported.
    assert c2["delta"] == pytest.approx(0.375)
    assert c2["guardrail"]["delta_sensitivity"] == pytest.approx(-0.25)
    assert not c2["guardrail"]["holds"] and not c2["supported"]
    assert c2["per_pair_delta"]["mosaic"] == pytest.approx(0.0)
    assert c2["per_stratum_delta"]["mosaic_snv_pos"] == pytest.approx(-1.0)
    # Holm over two: the larger raw p is adjusted to max(2·min, max).
    ps = sorted([c1["p"], c2["p"]])
    assert sorted([c1["holm_p"], c2["holm_p"]]) == pytest.approx(
        [min(1, 2 * ps[0]), max(min(1, 2 * ps[0]), ps[1])]
    )
    assert (
        c1["per_model"]["opus"]["mcnemar"]["a_only"] == 6 * 5
    )  # 4 hard + 1 anchor negative strata
    sign = s["directional_sign_test"]
    # Positives are unchanged under C1, so segdup_cmrg_pos (predicted "-") has Δ = 0.
    assert sign["n_counted"] == 6 and sign["n_matched"] == 6


def test_null_contrast_is_inconclusive_or_equivalent(tmp_path: Path) -> None:
    s = _grade(make_r4_run_set(tmp_path, n=6, answer_fn=lambda lab, c, m: oracle_answer(lab)))
    c1 = s["contrasts"]["C1"]
    assert c1["delta"] == 0.0 and c1["p"] == 1.0 and not c1["supported"]
    assert c1["verdict"] in ("inconclusive", "no meaningful difference (TOST)")
    assert c1["tost"]["equivalent"] is True  # zero spread: every shifted unit has one sign


def test_contrast_delta_equals_primary_difference(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        flip = (int(lab["item"][2:]) + MODEL_IX[m] + COND_IX[c]) % 3 == 0
        return ARTIFACT if flip and lab["polarity"] == "pos" else oracle_answer(lab)

    s = _grade(make_r4_run_set(tmp_path, n=5, answer_fn=fn))
    for name, (a, b) in (("C1", ("Traw-r4+Tcur", "Traw-r4")), ("C2", ("V+Traw-r4", "Traw-r4"))):
        pa = s["groups"][f"__all__|{a}"]["primary"]
        pb = s["groups"][f"__all__|{b}"]["primary"]
        assert s["contrasts"][name]["delta"] == pytest.approx(pa - pb)


MODEL_IX = {"haiku": 0, "sonnet": 1, "opus": 2}
COND_IX = {"Traw-r4": 0, "Traw-r4+Tcur": 1, "V+Traw-r4": 2}


def test_balanced_weights_sum_to_one_and_skip_anchors() -> None:
    labels = {
        "a": {"stratum": "mosaic_snv_pos", "case_type": "x"},
        "b": {"stratum": "mosaic_snv_pos", "case_type": "x"},
        "c": {"stratum": "mosaic_caller_fp_neg", "case_type": "x"},
        "d": {"stratum": "tr_indel_pos", "case_type": "x"},  # tr has no negative: dropped
        "e": {"stratum": "anchor_het_snv", "case_type": "x"},
        "f": {"case_type": "real_snv_het"},  # caption-check item: outside every pair
    }
    w = balanced_weights(labels)
    assert w["a"] == w["b"] == pytest.approx(0.25) and w["c"] == pytest.approx(0.5)
    assert w["d"] == w["e"] == w["f"] == 0.0
    assert sum(w.values()) == pytest.approx(1.0)
    assert pair_of({"case_type": "real_snv_het"}) is None
    assert pair_of({"pair": "tr", "polarity": "neg", "case_type": "z"}) == ("tr", "neg")


def test_directional_sign_test_known_answer() -> None:
    deltas = {
        "C1": {
            "mosaic_caller_fp_neg": 0.1,  # + predicted, matched
            "tr_stutter_neg": -0.05,  # + predicted, not matched
            "segdup_cmrg_pos": -0.2,  # - predicted, matched
            "indel_adjacent_misalign_neg": 0.0,  # zero: not counted
            "mosaic_snv_pos": -0.3,  # 0/- : not counted
        },
        "C2": {
            "mosaic_caller_fp_neg": 0.1,
            "segdup_psv_neg": 0.2,
            "indel_adjacent_misalign_neg": None,
        },
    }
    r = directional_sign_test(deltas)
    assert r["n_counted"] == 5 and r["n_matched"] == 4
    assert r["p_one_sided"] == pytest.approx(6 / 32)
    assert len(r["cells"]) == 16
    assert directional_sign_test({})["p_one_sided"] is None


def _fetch(seq: str):
    return lambda contig, s, e: seq[max(0, s) : e]


def test_allele_equivalence_repeat_shift() -> None:
    ref = "C" * 70 + "GATTTTTTCG" + "C" * 70  # T-run at 72..77
    f = _fetch(ref)
    kw = {"contig": "c", "ref_fetch": f}
    # Truth: delete one T at the run's left end (POI = first deleted base, 72).
    assert allele_equivalent("T", "T", "del", pos0=72, **kw)
    assert allele_equivalent("-t", "T", "del", pos0=72, **kw)
    assert allele_equivalent("TT", "TT", "del", pos0=72, **kw)
    assert not allele_equivalent("TT", "T", "del", pos0=72, **kw)
    assert not allele_equivalent("A", "T", "del", pos0=72, **kw)
    # Dinucleotide repeat: CA vs AC is the same deletion shifted by one.
    ref2 = "G" * 70 + "TACACACACT" + "G" * 70
    f2 = {"contig": "c", "ref_fetch": _fetch(ref2)}
    assert allele_equivalent("CA", "AC", "del", pos0=71, **f2)
    assert allele_equivalent("CA", "AC", "ins", pos0=70, **f2)
    assert not allele_equivalent("CC", "AC", "ins", pos0=70, **f2)
    # Insertion of a T anywhere in the run is one haplotype.
    assert allele_equivalent("T", "T", "ins", pos0=71, **kw)
    assert not allele_equivalent("G", "T", "ins", pos0=71, **kw)
    # SNVs and missing references: string rule only.
    assert allele_equivalent("g", "G", "SNV", contig="c", pos0=5, ref_fetch=None)
    assert not allele_equivalent("A", "G", "SNV", contig="c", pos0=5, ref_fetch=f)
    assert not allele_equivalent("T", "TT", "del", contig="c", pos0=72, ref_fetch=None)
    assert not allele_equivalent("T,A", "T", "del", pos0=72, **kw)
    # A truth that does not match the reference never matches a shifted answer.
    assert not allele_equivalent("TT", "GG", "del", pos0=72, **kw)


def test_genotype_uses_equivalence(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        ans = oracle_answer(lab)
        if lab["stratum"] == "tr_indel_pos":
            ans["alt_allele"] = "-T"
        if lab["stratum"] == "tr_stutter_neg":
            ans["artifact_type"] = "strand_bias"  # scored on stratum 4 only
        return ans

    rs = make_r4_run_set(tmp_path, conditions=("Traw-r4",), answer_fn=fn)
    s = _grade(rs)
    g = s["groups"]["__all__|Traw-r4"]
    assert g["genotype_accuracy"] == pytest.approx(7 / 8)
    assert g["artifact_type_accuracy_tr_stutter"] == 0.0
    graded = grade_runs_r4(rs)
    tr = next(x for x in graded if x.label["stratum"] == "tr_indel_pos")
    assert genotype_correct(tr, None)  # "-T" normalises to "T" without a reference
    rows = list(csv.DictReader((rs / "results" / "results.csv").open()))
    assert {r["stratum"] for r in rows} >= set(HARD_STRATA_R4)
    assert all(r["utility"] in ("1.0", "0.0", "-1.0") for r in rows)


def test_truth_source_missing_and_exclusions(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict | None:
        if lab["item"] == "it000" and m == "opus":
            return None  # missing answer: wrong
        return oracle_answer(lab) if lab["polarity"] == "pos" else REAL

    rs = make_r4_run_set(tmp_path, conditions=("Traw-r4",), answer_fn=fn)
    s = _grade(rs)
    g = s["groups"]["__all__|Traw-r4"]
    assert g["accuracy_by_truth_source"]["cmrg_v1.00"] == 0.0
    assert g["accuracy_by_truth_source"]["v5.0q"] < 1.0
    assert g["n_missing"] == 1 and g["invalid_rate"] > 0
    s2 = _grade(rs, exclude_items={"it000"}, out_dir=tmp_path / "res2")
    assert s2["groups"]["__all__|Traw-r4"]["accuracy_by_truth_source"]["v5.0q"] == 1.0
    assert s2["n_excluded_items"] == 1
    summary = json.loads((rs / "results" / "summary.json").read_text())
    assert summary["round"] == 4 and summary["contrasts"] == {}


def test_r4_selftests(tmp_path: Path) -> None:
    rs = make_r4_run_set(tmp_path, n=5)
    res = run_selftests_r4(rs, n_flips=200)
    assert [r.name for r in res] == [
        "r4 oracle = 1.0",
        "r4 always-real = 0.5 per pair",
        "r4 always-artifact = 0.5 per pair",
        "r4 random ≈ chance",
    ]
    assert all(r.ok for r in res), [(r.name, r.detail) for r in res]


def test_cli_grade_r4_and_selftest(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    from bamcp.eval.visual.cli_r4 import main

    rs = make_r4_run_set(tmp_path, n=4)
    excl = tmp_path / "excl.txt"
    excl.write_text("it000\n")
    ref = tmp_path / "ref.fa"
    assert (
        main(
            [
                "grade-r4",
                "--run-dir",
                str(rs),
                "--n-flips",
                "500",
                "--exclude-items",
                str(excl),
                "--reference",
                str(ref),
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "C1 Traw-r4+Tcur vs Traw-r4" in out and "directional sign test" in out
    assert json.loads((rs / "results" / "summary.json").read_text())["n_excluded_items"] == 1
    assert main(["selftest-r4", "--run-dir", str(rs)]) == 0
    assert capsys.readouterr().out.count("PASS") == 4


def test_mpileup_run_set_uses_the_raw_input_contrasts(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        hurt = c != "Tmp" and lab["polarity"] == "pos" and int(lab["item"][2:]) % 2 == 0
        return ARTIFACT if hurt else oracle_answer(lab)

    rs = make_r4_run_set(tmp_path, n=4, conditions=("Tmp", "Tmp+Tcur", "V+Tmp"), answer_fn=fn)
    s = _grade(rs)
    assert {k: (v["a"], v["b"]) for k, v in s["contrasts"].items()} == {
        "C1": ("Tmp+Tcur", "Tmp"),
        "C2": ("V+Tmp", "Tmp"),
    }
    assert all(v["delta"] < 0 and not v["supported"] for v in s["contrasts"].values())


def test_ref_fetch_follows_the_builders_relative_path(tmp_path: Path, monkeypatch) -> None:
    from bamcp.eval.visual.grade import _build_ref_fetch

    run = tmp_path / "experiment"
    (run / "fixtures").mkdir(parents=True)
    ref = run / "fixtures" / "ref.fa"
    ref.write_text(">chr1\nACGTACGTAC\n")
    Path(str(ref) + ".fai").write_text("chr1\t10\t6\t10\t11\n")
    # build-r4 --out experiment, run from tmp_path: build.json records a cwd-relative path.
    (run / "build.json").write_text(json.dumps({"reference": "experiment/fixtures/ref.fa"}))
    monkeypatch.chdir(tmp_path)
    assert _build_ref_fetch(run) is not None
    monkeypatch.chdir(
        run.parent.parent
    )  # elsewhere, with a run-relative path: the fallback finds it
    (run / "build.json").write_text(json.dumps({"reference": "fixtures/ref.fa"}))
    assert _build_ref_fetch(run) is not None


def test_genotype_endpoint_joins_holm_only_when_confirmatory(tmp_path: Path) -> None:
    def fn(lab: dict, c: str, m: str) -> dict:
        if c == "Traw-r4+Tcur" or lab["variant_present"] != "yes":
            return oracle_answer(lab)
        # Baseline: right call, wrong zygosity on every true variant.
        wrong = "hom_alt" if lab["zygosity"] != "hom_alt" else "het"
        return {**oracle_answer(lab), "zygosity": wrong}

    rs = make_r4_run_set(tmp_path, n=6, conditions=("Traw-r4", "Traw-r4+Tcur"), answer_fn=fn)
    plain = _grade(rs)["contrasts"]["C1"]["genotype"]
    assert plain["a_only"] == plain["n_pairs"] and plain["b_only"] == 0
    assert plain["b_correct"] == 0 and plain["p"] < 0.001
    assert "holm_p" not in plain
    conf = _grade(rs, genotype_confirmatory=True)["contrasts"]["C1"]
    assert conf["genotype"]["holm_p"] >= conf["genotype"]["p"]
    assert conf["genotype"]["supported"]
    assert conf["holm_p"] >= conf["p"]

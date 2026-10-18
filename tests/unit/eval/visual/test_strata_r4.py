"""Round-4 strata: frozen quotas, knobs, predicates, recoverability, spacing, matching and
the seeded draw (``real.select_r4``) on synthetic candidate rows."""

from __future__ import annotations

from collections import Counter
from typing import Any

import pytest

from bamcp.eval.visual import real
from bamcp.eval.visual import strata_r4 as S

pytestmark = pytest.mark.unit


def row(
    stratum: str,
    site0: int,
    *,
    chrom: str = "chr1",
    cluster: str | None = None,
    partition: str = "pilot",
    label: dict[str, Any] | None = None,
    context: dict[str, Any] | None = None,
    **features: Any,
) -> dict[str, Any]:
    pair, polarity = S.PAIRS[stratum]
    feats = {"depth": 60, "alt_reads": 8, "vaf": 8 / 60, **features}
    return {
        "chrom": chrom,
        "site0": site0,
        "stratum": stratum,
        "pair": pair,
        "polarity": polarity,
        "cluster": cluster or f"bin:{chrom}:{site0 // 1_000_000}",
        "partition": partition,
        "knob": None,
        "truth_source": "t",
        "vcf_pos": site0 + 1,
        "vcf_ref": "A",
        "vcf_alt": "G",
        "label": label or {"variant_type": "SNV", "alt_allele": "G", "zygosity": "het"},
        "features": feats,
        "checks": {},
        "context": context or {},
        "depth": feats["depth"],
        "alt_reads": feats["alt_reads"],
        "exclude": None,
    }


# --- frozen tables -----------------------------------------------------------------------


def test_strata_and_pairs_match_the_contract():
    assert S.STRATA_R4 == (
        "mosaic_snv_pos", "mosaic_caller_fp_neg", "tr_indel_pos", "tr_stutter_neg",
        "segdup_cmrg_pos", "segdup_psv_neg", "indel_adjacent_misalign_neg",
        "indel_adjacent_snv_pos", "anchor_het_snv", "anchor_hom_snv", "anchor_clean_neg",
    )  # fmt: skip
    assert set(S.PAIRS) == set(S.STRATA_R4)
    for p in S.PAIR_NAMES:
        pols = sorted(pol for s, (pp, pol) in S.PAIRS.items() if pp == p)
        assert pols == ["neg", "pos"]


def test_quotas_are_frozen():
    main = S.quotas("main")["default"]
    assert sum(main.values()) == 128
    assert all(main[s] == 15 for s in S.HARD_STRATA)
    assert [main[s] for s in S.ANCHOR_STRATA] == [2, 3, 3]
    pilot = S.quotas("pilot")
    assert sum(pilot["default"].values()) + sum(pilot["alt"].values()) == 84
    assert all(pilot["default"][s] == 6 and pilot["alt"][s] == 4 for s in S.HARD_STRATA)
    assert sum(pilot["default"][s] for s in S.ANCHOR_STRATA) == 4
    with pytest.raises(ValueError, match="split"):
        S.quotas("other")
    # quotas() returns copies: callers cannot mutate the frozen table.
    S.quotas("main")["default"]["mosaic_snv_pos"] = 0
    assert S.QUOTAS_R4["main"]["default"]["mosaic_snv_pos"] == 15


def test_every_hard_stratum_has_two_knob_levels():
    for s in S.HARD_STRATA:
        assert set(S.KNOBS[s]) == {"default", "alt"} and S.KNOBS[s]["default"] != S.KNOBS[s]["alt"]


# --- labels ------------------------------------------------------------------------------


def test_label_for():
    mos = S.label_for("mosaic_snv_pos", row("mosaic_snv_pos", 1))
    assert mos["zygosity"] == "low_vaf" and mos["artifact_verdict"] == "real"
    assert mos["variant_present"] == "yes" and mos["alt_allele"] == "G"
    st = S.label_for("tr_stutter_neg", row("tr_stutter_neg", 1))
    assert st["artifact_type"] == "homopolymer_indel" and st["variant_present"] == "no"
    for s in ("mosaic_caller_fp_neg", "segdup_psv_neg", "indel_adjacent_misalign_neg"):
        lab = S.label_for(s, row(s, 1))
        assert lab["artifact_type"] == "unknown" and lab["artifact_verdict"] == "artifact"
    # Unscored like every non-stutter negative (CONTRACT.md); round 3 used "none" here.
    clean = S.label_for("anchor_clean_neg", row("anchor_clean_neg", 1))
    assert clean["artifact_type"] == "unknown" and clean["variant_present"] == "no"
    hom = S.label_for(
        "anchor_hom_snv",
        row(
            "anchor_hom_snv",
            1,
            label={"variant_type": "SNV", "alt_allele": "T", "zygosity": "hom_alt"},
        ),
    )
    assert hom["zygosity"] == "hom_alt" and hom["alt_allele"] == "T"


# --- recoverability + predicates ---------------------------------------------------------------


@pytest.mark.parametrize(
    ("stratum", "zyg", "feats", "ok"),
    [
        ("anchor_het_snv", "het", {"depth": 9, "alt_reads": 5, "vaf": 0.5}, False),
        ("anchor_het_snv", "het", {"depth": 40, "alt_reads": 3, "vaf": 0.10}, True),
        ("anchor_het_snv", "het", {"depth": 40, "alt_reads": 2, "vaf": 0.4}, False),
        ("anchor_het_snv", "het", {"depth": 40, "alt_reads": 3, "vaf": 0.08}, False),
        ("anchor_hom_snv", "hom_alt", {"depth": 40, "alt_reads": 24, "vaf": 0.60}, True),
        ("anchor_hom_snv", "hom_alt", {"depth": 40, "alt_reads": 23, "vaf": 0.58}, False),
        ("mosaic_snv_pos", "low_vaf", {"depth": 60, "alt_reads": 4, "vaf": 0.066}, True),
        ("mosaic_snv_pos", "low_vaf", {"depth": 60, "alt_reads": 3, "vaf": 0.05}, False),
    ],
)
def test_recoverable(stratum, zyg, feats, ok):
    r = row(stratum, 1, label={"zygosity": zyg}, **feats)
    assert S.recoverable(r)[0] is ok


def test_mosaic_predicates_and_knob():
    r = row("mosaic_snv_pos", 1, depth=60, alt_reads=5, vaf=5 / 60)
    assert S.eligible(r, "default") == (True, None)
    assert S.eligible(r, "alt") == (False, "knob_alt_band")
    assert S.eligible(row("mosaic_snv_pos", 1, depth=24, alt_reads=8), "default")[1] == (
        "depth_lt_25"
    )
    fp = row("mosaic_caller_fp_neg", 1, depth=60, alt_reads=4, vaf=0.066)
    assert S.eligible(fp, "default")[0] and not S.eligible(fp, "alt")[0]
    # Contamination above 1% raises stratum 2's floor (design 3.2).
    over = {"mosaic_caller_fp_neg": {"vaf_floor_override": 0.08}}
    assert not S.eligible(fp, "default", over)[0]
    assert S.eligible(row("mosaic_caller_fp_neg", 1, vaf=0.25), "default")[1] == "knob_vaf_band"


def _tr(stratum: str, *, run_length: int, zyg: str = "het", unit: str = "A", **f: Any):
    return row(
        stratum,
        1,
        label={"zygosity": zyg, "variant_type": "del", "alt_allele": unit},
        context={"run": [1, 1 + run_length], "unit": unit, "delta": -len(unit)},
        run_length=run_length,
        **f,
    )


def test_tr_predicates():
    easy = _tr("tr_indel_pos", run_length=11, alt_reads=30, vaf=0.5, stutter_frac=0.05)
    assert S.eligible(easy, "default")[1] == "not_difficult"
    hard = _tr("tr_indel_pos", run_length=11, alt_reads=30, vaf=0.5, stutter_frac=0.10)
    assert S.eligible(hard, "default")[0]
    assert S.eligible(hard, "alt")[1] == "knob_run_floor"  # 11 bp < 12
    low = _tr("tr_indel_pos", run_length=14, alt_reads=20, vaf=0.38, stutter_frac=0.0)
    assert S.eligible(low, "default")[0] and S.eligible(low, "alt")[0]
    hom = _tr("tr_indel_pos", run_length=14, zyg="hom_alt", alt_reads=40, vaf=0.9, stutter_frac=0)
    assert not S.eligible(hom, "default")[0]  # 0.90 > 0.85 and no stutter
    s = _tr("tr_indel_pos", run_length=22, unit="AC", alt_reads=20, vaf=0.3)
    assert S.eligible(s, "default")[0]
    neg = _tr("tr_stutter_neg", run_length=12, alt_reads=4, vaf=0.10)
    assert S.eligible(neg, "alt")[0]
    assert S.eligible(_tr("tr_stutter_neg", run_length=12, alt_reads=3, vaf=0.05), "default")[0]
    assert not S.eligible(_tr("tr_stutter_neg", run_length=12, alt_reads=5, vaf=0.04), "default")[0]
    assert S.eligible(_tr("tr_stutter_neg", run_length=12, alt_reads=2, vaf=0.2), "default")[1] == (
        "stutter_signal_weak"
    )


def test_segdup_and_adjacent_predicates():
    psv = row("segdup_psv_neg", 1, depth=40, alt_reads=3, vaf=0.075)
    assert not S.eligible(psv, "default")[0] and S.eligible(psv, "alt")[0]
    adj = row("indel_adjacent_misalign_neg", 1, depth=40, alt_reads=5, vaf=0.125)
    assert S.eligible(adj, "default")[0] and not S.eligible(adj, "alt")[0]
    pos = row("indel_adjacent_snv_pos", 1, depth=40, alt_reads=12, vaf=0.3)
    assert S.eligible(pos, "default")[0] and S.eligible(pos, "alt")[1] == "knob_alt_reads"
    assert S.eligible(row("anchor_clean_neg", 1, depth=30, alt_reads=0), "default")[0]
    assert not S.eligible(row("anchor_clean_neg", 1, depth=30, alt_reads=1), "default")[0]
    hard = row("segdup_cmrg_pos", 1, depth=40, alt_reads=10, vaf=0.25)
    assert S.cmrg_difficult(hard)
    mq = row(
        "segdup_cmrg_pos",
        1,
        depth=40,
        alt_reads=20,
        vaf=0.5,
        alt_mapq_lt20_frac=0.5,
        ref_mapq_lt20_frac=0.2,
    )
    assert S.mapq_lt20_frac(mq) == pytest.approx(0.35) and S.cmrg_difficult(mq)
    assert not S.cmrg_difficult(row("segdup_cmrg_pos", 1, vaf=0.45))


def test_subquotas():
    assert S.subquotas("mosaic_snv_pos", "default", 15) == {"4-6": 7, "ge7": 8}
    assert S.subquotas("mosaic_snv_pos", "alt", 4) == {"ge7": 4}
    assert S.subquotas("tr_indel_pos", "default", 15) == {"het": 10, "hom": 5}
    assert S.subquotas("segdup_cmrg_pos", "default", 6) == {"difficult": 3, "other": 3}
    assert S.subquotas("segdup_cmrg_pos", "alt", 4) == {"difficult": 3, "other": 1}
    assert S.subquotas("segdup_psv_neg", "default", 6) == {"any": 6}


# --- spacing + matching ----------------------------------------------------------------------


def test_spacing():
    sp = S.Spacing()
    a = row("anchor_het_snv", 100_000, cluster="sd:x")
    sp.add(a)
    assert sp.why_not(row("anchor_het_snv", 105_000, cluster="bin:y")) == "spacing_10kb"
    assert sp.ok(row("anchor_het_snv", 105_000, chrom="chr2", cluster="bin:y"))
    sp.add(row("anchor_het_snv", 200_000, cluster="sd:x"))
    assert sp.why_not(row("anchor_het_snv", 900_000, cluster="sd:x")) == "cluster_full"


def test_match_negatives_exact_then_relaxed():
    pos = [
        row("mosaic_snv_pos", i * 100_000, alt_reads=a, depth=d)
        for i, (a, d) in enumerate([(5, 40), (8, 60), (12, 80)])
    ]
    negs = [
        row("mosaic_caller_fp_neg", 5_000_000 + i * 100_000, alt_reads=a, depth=d)
        for i, (a, d) in enumerate([(5, 41), (8, 61), (8, 62), (12, 30)])
    ]
    got, counts = S.match_negatives(pos, negs, "mosaic_caller_fp_neg", seed=1)
    assert len(got) == 3 and counts == {}
    levels = [g["match"]["level"] for g in got]
    assert levels[:2] == [0, 0] and levels[2] == 2  # ge10/gt70 only relaxes to the alt bin
    # accept() vetoes are counted and the next candidate in the bin is used.
    got, counts = S.match_negatives(
        pos[1:2],
        negs,
        "mosaic_caller_fp_neg",
        seed=1,
        accept=lambda r: "hifi" if r["site0"] == 5_100_000 else None,
    )
    assert got[0]["site0"] == 5_200_000 and counts == {"hifi": 1}
    got, counts = S.match_negatives(pos, [], "mosaic_caller_fp_neg", seed=1)
    assert got == [] and counts == {"unmatched": 3}


def test_match_keys():
    r = row(
        "segdup_psv_neg", 1, depth=40, alt_reads=4, alt_mapq_lt20_frac=1.0, ref_mapq_lt20_frac=0.0
    )
    assert S.match_key("segdup_psv_neg", r) == ("0.10-0.30", "lt50")
    tr = _tr("tr_stutter_neg", run_length=13)
    assert S.match_key("tr_stutter_neg", tr) == ("hp:12-15",)
    assert S.match_key("tr_stutter_neg", tr, level=1) == ()


def test_ks_two_sample():
    d, p = S.ks_two_sample([0.0] * 10, [0.0] * 10)
    assert d == 0 and p == 1.0
    d, p = S.ks_two_sample([0.0] * 20, [1.0] * 20)
    assert d == 1.0 and p < 0.001
    assert S.ks_two_sample([], [1.0]) == (0.0, 1.0)
    # Ties: 4 insertion-anchor offsets among 30 vs 30 zeros is a redraw (p < 0.2), as in the
    # leak-lint; a continuous-null p (~0.94) would pass it.
    _, p = S.ks_two_sample([-1 / 12] * 4 + [0.0] * 26, [0.0] * 30)
    assert p < 0.2


def test_offset_in_run():
    assert S.offset_in_run(_tr("tr_stutter_neg", run_length=12)) == 0.0
    assert S.offset_in_run(row("anchor_het_snv", 1)) is None


# --- the draw ----------------------------------------------------------------------------------


def _pool(partition: str = "pilot") -> list[dict[str, Any]]:
    """Plenty of eligible rows for every stratum, 1 Mb apart (one cluster each)."""
    rows: list[dict[str, Any]] = []
    k = 0

    def add(r: dict[str, Any]) -> None:
        rows.append(r)

    for i in range(40):
        for s in S.STRATA_R4:
            k += 1
            site = k * 1_000_000
            chrom = f"chr{1 + k % 22}"
            alt = [5, 8, 12][i % 3]
            depth = [45, 60, 80][i % 3]
            base = dict(chrom=chrom, partition=partition, cluster=f"bin:{chrom}:{k}")
            if s == "mosaic_snv_pos":
                add(row(s, site, alt_reads=alt, depth=depth, vaf=alt / depth, **base))
            elif s == "mosaic_caller_fp_neg":
                add(row(s, site, alt_reads=alt, depth=depth, vaf=0.12, **base))
            elif s == "tr_indel_pos":
                zyg = "hom_alt" if i % 3 == 0 else "het"
                r = _tr(
                    s,
                    run_length=12 + i % 4,
                    zyg=zyg,
                    alt_reads=20,
                    vaf=0.8 if zyg == "hom_alt" else 0.35,
                    stutter_frac=0.15,
                )
                add({**r, "site0": site, **base})
            elif s == "tr_stutter_neg":
                r = _tr(s, run_length=12 + i % 4, alt_reads=6, vaf=0.15)
                add({**r, "site0": site, **base})
            elif s == "segdup_cmrg_pos":
                add(
                    row(
                        s,
                        site,
                        alt_reads=10,
                        depth=40,
                        vaf=0.25 if i % 2 else 0.5,
                        alt_mapq_lt20_frac=0.0,
                        ref_mapq_lt20_frac=0.0,
                        **base,
                    )
                )
            elif s == "segdup_psv_neg":
                add(
                    row(
                        s,
                        site,
                        alt_reads=6,
                        depth=40,
                        vaf=0.15,
                        alt_mapq_lt20_frac=0.0,
                        ref_mapq_lt20_frac=0.0,
                        **base,
                    )
                )
            elif s == "indel_adjacent_snv_pos":
                add(row(s, site, alt_reads=6 + i % 3, depth=40, vaf=0.4, **base))
            elif s == "indel_adjacent_misalign_neg":
                add(row(s, site, alt_reads=8, depth=40, vaf=0.2, **base))
            elif s == "anchor_het_snv":
                add(row(s, site, alt_reads=20, depth=40, vaf=0.5, **base))
            elif s == "anchor_hom_snv":
                add(
                    row(
                        s,
                        site,
                        alt_reads=40,
                        depth=40,
                        vaf=1.0,
                        label={"zygosity": "hom_alt", "alt_allele": "G"},
                        **base,
                    )
                )
            else:
                add(row(s, site, alt_reads=0, depth=40, vaf=0.0, **base))
    return rows


def test_select_pilot_fills_every_frozen_cell():
    rows = _pool()
    chosen, rep = real.select_r4(rows, split="pilot", seed=7, log=lambda *_: None)
    n = Counter(r["stratum"] for r in chosen)
    assert sum(n.values()) == 84 and not rep["short"]
    assert all(n[s] == 10 for s in S.HARD_STRATA)
    knobs = Counter((r["stratum"], r["knob"]) for r in chosen)
    assert all(knobs[(s, "default")] == 6 and knobs[(s, "alt")] == 4 for s in S.HARD_STRATA)
    mos = [r for r in chosen if r["stratum"] == "mosaic_snv_pos"]
    assert Counter(S.alt_band(r["alt_reads"]) for r in mos if r["knob"] == "default") == {
        "4-6": 3,
        "ge7": 3,
    }
    assert all(r["alt_reads"] >= 7 for r in mos if r["knob"] == "alt")
    assert all(r["label"]["zygosity"] == "low_vaf" for r in mos)
    tr = [r for r in chosen if r["stratum"] == "tr_indel_pos" and r["knob"] == "default"]
    assert Counter(r["label"]["zygosity"] for r in tr) == {"het": 4, "hom_alt": 2}
    fps = [r for r in chosen if r["stratum"] == "mosaic_caller_fp_neg"]
    assert all("match" in r for r in fps)
    assert "exclude" not in chosen[0]
    # Deterministic, and seed-dependent.
    again, _ = real.select_r4(rows, split="pilot", seed=7, log=lambda *_: None)
    assert [(r["chrom"], r["site0"]) for r in again] == [(r["chrom"], r["site0"]) for r in chosen]
    other, _ = real.select_r4(rows, split="pilot", seed=8, log=lambda *_: None)
    assert {(r["chrom"], r["site0"]) for r in other} != {(r["chrom"], r["site0"]) for r in chosen}


def test_select_respects_partition_spacing_and_checks():
    rows = _pool("main")
    chosen, rep = real.select_r4(rows, split="pilot", seed=1, log=lambda *_: None)
    assert chosen == [] and rep["short"]  # nothing in the pilot partition
    chosen, rep = real.select_r4(
        rows,
        split="main",
        seed=1,
        checks=lambda r: "hifi_alt_absent" if r["stratum"] == "anchor_hom_snv" else None,
        knob_overrides={"mosaic_snv_pos": "alt"},
        log=lambda *_: None,
    )
    n = Counter(r["stratum"] for r in chosen)
    assert n["mosaic_snv_pos"] == 15 and n["anchor_hom_snv"] == 0
    assert all(r["knob"] == "alt" for r in chosen if r["stratum"] == "mosaic_snv_pos")
    assert all(r["alt_reads"] >= 7 for r in chosen if r["stratum"] == "mosaic_snv_pos")
    assert rep["exclusions"]["anchor_hom_snv"]["hifi_alt_absent"] > 0
    assert rep["short"]["anchor_hom_snv|default|any"] == {"need": 3, "got": 0}
    sp = S.Spacing()
    for r in chosen:
        assert sp.ok(r)
        sp.add(r)


def test_select_cluster_cap():
    """Twenty candidates in one paralog cluster yield at most two items."""
    rows = [
        row("anchor_het_snv", i * 50_000, cluster="sd:fam", alt_reads=20, depth=40, vaf=0.5)
        for i in range(20)
    ]
    chosen, rep = real.select_r4(rows, split="pilot", seed=1, pairs=["anchor"], log=lambda *_: None)
    assert len([r for r in chosen if r["stratum"] == "anchor_het_snv"]) == 1  # quota 1
    rows += [
        row(
            "anchor_clean_neg",
            5_000_000 + i * 50_000,
            cluster="sd:fam",
            alt_reads=0,
            depth=40,
            vaf=0.0,
        )
        for i in range(5)
    ]
    chosen, rep = real.select_r4(rows, split="pilot", seed=1, pairs=["anchor"], log=lambda *_: None)
    fam = [r for r in chosen if r["cluster"] == "sd:fam"]
    assert len(fam) == 2 and rep["short"]


def test_pools_report():
    rows = _pool("main")
    rows.append({**rows[0], "exclude": "v5q_within_30bp"})
    rep = real.pools_report_r4(rows)
    m = rep["mosaic_snv_pos"]
    assert m["excluded"] == {"v5q_within_30bp": 1}
    assert m["eligible"]["main"]["default"] == 40 and m["g5_main_pool_ge_2x"]
    assert m["eligible"]["main"]["alt"] == 26  # alt >= 7 only: 8 and 12
    assert "vaf" in m["feature_means_default"]
    assert rep["anchor_clean_neg"]["main_quota"] == 3


def test_matched_negatives_respect_spacing_among_themselves():
    """Matched negatives 40 bp apart: only one may be drawn (spacing applies as each is
    taken, not after the whole match)."""
    pos = [
        row(
            "segdup_cmrg_pos",
            i * 1_000_000,
            chrom="chr2",
            cluster=f"c{i}",
            depth=40,
            alt_reads=10,
            vaf=0.25,
        )
        for i in range(1, 7)
    ]
    negs = [
        row(
            "segdup_psv_neg",
            50_000_000 + 40 * i,
            chrom="chr3",
            cluster=f"n{i}",
            depth=40,
            alt_reads=6,
            vaf=0.15,
        )
        for i in range(6)
    ]
    negs += [
        row(
            "segdup_psv_neg",
            70_000_000 + 1_000_000 * i,
            chrom="chr3",
            cluster=f"m{i}",
            depth=40,
            alt_reads=6,
            vaf=0.15,
        )
        for i in range(6)
    ]
    chosen, rep = real.select_r4(
        pos + negs, split="pilot", seed=3, pairs=["segdup"], log=lambda *_: None
    )
    got = sorted(r["site0"] for r in chosen if r["stratum"] == "segdup_psv_neg")
    assert sum(50_000_000 <= s < 50_001_000 for s in got) <= 1
    assert all(b - a >= S.MIN_SPACING_R4 for a, b in zip(got, got[1:], strict=False))


def test_select_quota_overrides_replace_a_stratum_quota():
    rows = _pool()
    chosen, _ = real.select_r4(
        rows,
        split="pilot",
        seed=7,
        pairs=["tr"],
        quota_overrides={"tr_indel_pos": 3, "tr_stutter_neg": 3},
        log=lambda *_: None,
    )
    n = Counter((r["stratum"], r["knob"]) for r in chosen)
    assert n == {
        ("tr_indel_pos", "default"): 3,
        ("tr_indel_pos", "alt"): 3,
        ("tr_stutter_neg", "default"): 3,
        ("tr_stutter_neg", "alt"): 3,
    }

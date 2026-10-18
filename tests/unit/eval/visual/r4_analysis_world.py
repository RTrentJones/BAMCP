"""A tiny round-4 run set for the analysis tests: labels, runs, answers and a reference.

No stimuli are written (grading never reads ``stim/``). Every stratum gets ``n`` items;
items pair up into clusters; ``tr_indel_pos`` items are 1 bp deletions at the left end of a
``T`` homopolymer so allele-equivalence has something to shift.
"""

from __future__ import annotations

import json
import random
from collections.abc import Callable
from pathlib import Path
from typing import Any

from bamcp.eval.visual.grade import STRATUM_PAIRS_R4
from bamcp.eval.visual.selftest import constant_answer, oracle_answer

CONDS = ("Traw-r4", "Traw-r4+Tcur", "V+Traw-r4")
MODELS = ("haiku", "sonnet", "opus")
SLOT = 200
HP_LEN = 8
TRUTH_SOURCE = {"pos": "v5.0q", "neg": "cmrg_v1.00"}

AnswerFn = Callable[[dict[str, Any], str, str], dict[str, Any] | None]


def label_for(stratum: str, idx: int, item: str, pos1: int) -> dict[str, Any]:
    pair, pol = STRATUM_PAIRS_R4[stratum]
    base = {
        "item": item,
        "case_type": stratum,
        "stratum": stratum,
        "pair": pair,
        "polarity": pol,
        "cluster": f"cl:{stratum}:{idx // 2}",
        "truth_source": TRUTH_SOURCE[pol],
        "knob": "default" if idx % 5 < 3 else "alt",
        "instance": idx,
        "contig": "chr1",
        "positions": [pos1],
        "phase": "na",
    }
    if pol == "pos":
        if stratum == "tr_indel_pos":
            lab = {"variant_type": "del", "alt_allele": "T", "zygosity": "het"}
        else:
            zyg = "low_vaf" if stratum == "mosaic_snv_pos" else "het"
            if stratum == "anchor_hom_snv":
                zyg = "hom_alt"
            lab = {"variant_type": "SNV", "alt_allele": "G", "zygosity": zyg}
        lab.update(variant_present="yes", artifact_verdict="real", artifact_type="none")
    else:
        atype = {"tr_stutter_neg": "homopolymer_indel", "anchor_clean_neg": "none"}.get(
            stratum, "unknown"
        )
        lab = {
            "variant_present": "no",
            "variant_type": "none",
            "alt_allele": "",
            "zygosity": "none",
            "artifact_verdict": "artifact",
            "artifact_type": atype,
        }
    out = {**base, **lab}
    out["lenient_verdicts"] = [out["artifact_verdict"]]
    return out


def make_r4_run_set(
    root: Path,
    *,
    n: int = 4,
    models: tuple[str, ...] = MODELS,
    conditions: tuple[str, ...] = CONDS,
    answer_fn: AnswerFn | None = None,
    strata: tuple[str, ...] | None = None,
    seed: int = 0,
) -> Path:
    rng = random.Random(seed)  # noqa: S311 — reproducible fixture bases
    strata = strata or tuple(STRATUM_PAIRS_R4)
    labels: list[dict[str, Any]] = []
    seq: list[str] = []
    k = 0
    for s in strata:
        for i in range(n):
            pos0 = k * SLOT + SLOT // 2
            slot = [rng.choice("ACG") for _ in range(SLOT)]
            if s == "tr_indel_pos":
                for j in range(HP_LEN):
                    slot[SLOT // 2 + j] = "T"
            seq.extend(slot)
            labels.append(label_for(s, i, f"it{k:03d}", pos0 + 1))
            k += 1
    rs = root
    (rs / "labels").mkdir(parents=True, exist_ok=True)
    (rs / "answers").mkdir(exist_ok=True)
    ref = rs / "ref.fa"
    ref.write_text(">chr1\n" + "".join(seq) + "\n")
    import pysam

    pysam.faidx(str(ref))
    (rs / "build.json").write_text(json.dumps({"reference": str(ref), "round": 4}))
    with (rs / "labels" / "labels.jsonl").open("w") as f:
        for lab in labels:
            f.write(json.dumps(lab) + "\n")
    fn = answer_fn or (lambda lab, c, m: oracle_answer(lab))
    with (rs / "runs.jsonl").open("w") as f:
        for lab in labels:
            for c in conditions:
                for m in models:
                    rid = f"{lab['item']}-{c}-{m}"
                    run = {
                        "run_id": rid,
                        "stim": lab["item"],
                        "item": lab["item"],
                        "condition": c,
                        "model": m,
                        "repeat": 0,
                        "answer_path": f"answers/{rid}.json",
                    }
                    f.write(json.dumps(run) + "\n")
                    ans = fn(lab, c, m)
                    if ans is not None:
                        (rs / run["answer_path"]).write_text(json.dumps(ans))
    return rs


REAL = constant_answer("real")
ARTIFACT = constant_answer("artifact")


def uncertain_answer() -> dict[str, Any]:
    a = dict(REAL)
    a["artifact_verdict"] = "uncertain"
    return a

"""Grade a directory of model answers against the run set's labels.

Inputs: ``<run_dir>/runs.jsonl``, ``<run_dir>/labels/labels.jsonl`` and one answer file per
run (``answer_path`` in the manifest, relative to the run dir or ``--answers-dir``).
Outputs (``<run_dir>/results`` by default): ``results.csv`` (one row per run),
``summary.json`` (all metrics) and charts (see :mod:`.charts`).

Scoring rules are pre-registered in ``docs/eval/visual-preregistration.md``: an invalid or
missing answer is wrong on every field; excluded runs (``--exclude``) are dropped.

Round 4 (``docs/eval/visual-design-r4.md`` sections 5 and 8) has its own path,
:func:`grade_run_set_r4`: the primary is the mean over the 4 hard pos/neg pairs of balanced
strict joint accuracy, equal-weighted over models, anchors excluded; the two confirmatory
contrasts use the cluster sign-flip test with Holm, an inverted CI, a sensitivity/specificity
guardrail and TOST. Rounds 1-3 grading above is unchanged.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .answers import SCORED_FIELDS, field_correct, normalize_allele, parse_answer
from .stats import (
    N_FLIPS,
    SELECTIVE_COVERAGES,
    aurc,
    auroc,
    binom_sf,
    brier,
    cluster_bootstrap_mean,
    cluster_signflip_test,
    confident_wrong_rate,
    design_effect,
    ece,
    ece5,
    f1_scores,
    holm,
    mcnemar_exact,
    mean,
    paired_cluster_bootstrap,
    reliability_bins,
    selective_accuracy,
    signflip_ci,
    tost,
)

PRIMARY_CONTRASTS = (("V+Traw", "Traw"), ("V", "Tjson"))
SECONDARY_CONTRASTS = (("Tfull+V", "Tfull"), ("V", "Traw"), ("Tjson", "Traw"))
INVALID = "__invalid__"


@dataclass
class Graded:
    run: dict[str, Any]
    label: dict[str, Any]
    status: str  # valid | invalid | missing
    error: str = ""
    answer: dict[str, Any] | None = None
    correct: dict[str, bool] = field(default_factory=dict)

    @property
    def joint(self) -> bool:
        return self.correct.get("variant_present", False) and self.correct.get(
            "artifact_verdict", False
        )

    @property
    def joint_lenient(self) -> bool:
        if self.answer is None:
            return False
        return self.correct.get("variant_present", False) and (
            self.answer["artifact_verdict"] in self.label.get("lenient_verdicts", [])
        )


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def grade_run(run: dict[str, Any], label: dict[str, Any], answer_file: Path) -> Graded:
    if not answer_file.exists():
        return Graded(run, label, "missing", "no answer file")
    v = parse_answer(answer_file.read_text(encoding="utf-8", errors="replace"))
    if not v.ok or v.answer is None:
        return Graded(run, label, "invalid", v.error or "invalid")
    correct = {
        f: field_correct(f, v.answer, label)
        for f in SCORED_FIELDS
        if not (f == "artifact_type" and label["artifact_type"] == UNSCORED_TYPE)
    }
    return Graded(run, label, "valid", answer=v.answer, correct=correct)


UNSCORED_TYPE = "unknown"  # real-tier negatives: truth says reference, not why reads disagree


def _applicable(g: Graded, f: str) -> bool:
    if f == "phase":
        return bool(g.label["phase"] != "na")
    if f == "artifact_type":
        return bool(g.label["artifact_type"] != UNSCORED_TYPE)
    return True


def _clusters(rows: list[Graded], key: Any) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for g in rows:
        out[g.run["item"]].append(float(key(g)))
    return out


def summarize_group(rows: list[Graded], seed: int = 0) -> dict[str, Any]:
    valid = [g for g in rows if g.status == "valid"]
    acc: dict[str, float | None] = {}
    for f in SCORED_FIELDS:
        app = [g for g in rows if _applicable(g, f)]
        acc[f] = mean([float(g.correct.get(f, False)) for g in app]) if app else None
    joint, lo, hi = cluster_bootstrap_mean(_clusters(rows, lambda g: g.joint), seed=seed)
    truth_art = [g.label["artifact_verdict"] == "artifact" for g in rows]
    pred_art = [bool(g.answer and g.answer["artifact_verdict"] == "artifact") for g in rows]
    tp = sum(t and p for t, p in zip(truth_art, pred_art, strict=True))
    fp = sum(p and not t for t, p in zip(truth_art, pred_art, strict=True))
    fn = sum(t and not p for t, p in zip(truth_art, pred_art, strict=True))
    typed = [g for g in rows if _applicable(g, "artifact_type")]
    t_types = [g.label["artifact_type"] for g in typed]
    p_types = [g.answer["artifact_type"] if g.answer else INVALID for g in typed]
    f1 = f1_scores(t_types, p_types)
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for t, p in zip(t_types, p_types, strict=True):
        confusion[t][p] += 1
    conf = [float(g.answer["confidence"]) for g in valid if g.answer]
    corr = [g.joint for g in valid]
    per_case: dict[str, list[float]] = defaultdict(list)
    for g in rows:
        per_case[g.label["case_type"]].append(float(g.joint))
    return {
        "n_runs": len(rows),
        "n_items": len({g.run["item"] for g in rows}),
        "n_valid": len(valid),
        "n_invalid": sum(g.status == "invalid" for g in rows),
        "n_missing": sum(g.status == "missing" for g in rows),
        "primary_joint_accuracy": joint,
        "primary_joint_ci95": [lo, hi],
        "joint_lenient_accuracy": mean([float(g.joint_lenient) for g in rows]),
        "primary_joint_accuracy_excl_phasing": mean(
            [float(g.joint) for g in rows if g.label["phase"] == "na"]
        ),
        "field_accuracy": acc,
        "artifact_precision": tp / (tp + fp) if tp + fp else None,
        "artifact_recall": tp / (tp + fn) if tp + fn else None,
        "artifact_type_f1": f1,
        "artifact_type_macro_f1": mean(list(f1.values())) if f1 else None,
        "artifact_type_confusion": {t: dict(v) for t, v in confusion.items()},
        "brier": brier(conf, corr) if conf else None,
        "ece": ece(conf, corr) if conf else None,
        "reliability": reliability_bins(conf, corr) if conf else [],
        "per_case_type_joint": {k: mean(v) for k, v in sorted(per_case.items())},
    }


def contrast(rows: list[Graded], a: str, b: str, seed: int = 0) -> dict[str, Any] | None:
    ra = [g for g in rows if g.run["condition"] == a]
    rb = [g for g in rows if g.run["condition"] == b]
    if not ra or not rb:
        return None
    delta, lo, hi, n = paired_cluster_bootstrap(
        _clusters(ra, lambda g: g.joint), _clusters(rb, lambda g: g.joint), seed=seed
    )
    idx_b = {(g.run["item"], g.run["model"], g.run["repeat"]): g.joint for g in rb}
    pairs = [
        (g.joint, idx_b[k])
        for g in ra
        if (k := (g.run["item"], g.run["model"], g.run["repeat"])) in idx_b
    ]
    nb, nc, p = mcnemar_exact(pairs)
    per_case: dict[str, Any] = {}
    for case in sorted({g.label["case_type"] for g in ra}):
        ca = _clusters([g for g in ra if g.label["case_type"] == case], lambda g: g.joint)
        cb = _clusters([g for g in rb if g.label["case_type"] == case], lambda g: g.joint)
        d, dl, dh, dn = paired_cluster_bootstrap(ca, cb, seed=seed)
        per_case[case] = {"delta": d, "ci95": [dl, dh], "n_items": dn}
    return {
        "a": a,
        "b": b,
        "delta_joint": delta,
        "ci95": [lo, hi],
        "n_items": n,
        "mcnemar": {"a_only": nb, "b_only": nc, "p": p, "n_pairs": len(pairs)},
        "per_case_type": per_case,
    }


def grade_run_set(
    run_dir: str | Path,
    *,
    answers_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    exclude: set[str] | None = None,
    charts: bool = True,
    seed: int = 0,
) -> dict[str, Any]:
    run_dir = Path(run_dir)
    out = Path(out_dir) if out_dir else run_dir / "results"
    out.mkdir(parents=True, exist_ok=True)
    runs = load_jsonl(run_dir / "runs.jsonl")
    labels = {r["item"]: r for r in load_jsonl(run_dir / "labels" / "labels.jsonl")}
    exclude = exclude or set()

    graded: list[Graded] = []
    for run in runs:
        if run["run_id"] in exclude:
            continue
        if answers_dir is not None:
            answer_file = Path(answers_dir) / Path(run["answer_path"]).name
        else:
            answer_file = run_dir / run["answer_path"]
        graded.append(grade_run(run, labels[run["item"]], answer_file))

    groups: dict[tuple[str, str], list[Graded]] = defaultdict(list)
    for g in graded:
        groups[(g.run["model"], g.run["condition"])].append(g)
    by_group = {f"{m}|{c}": summarize_group(rows, seed) for (m, c), rows in sorted(groups.items())}

    models = sorted({g.run["model"] for g in graded})
    contrasts: dict[str, Any] = {}
    for kind, pairs in (("primary", PRIMARY_CONTRASTS), ("secondary", SECONDARY_CONTRASTS)):
        for a, b in pairs:
            for m in [*models, "__pooled__"]:
                rows = graded if m == "__pooled__" else [g for g in graded if g.run["model"] == m]
                res = contrast(rows, a, b, seed)
                if res is not None:
                    contrasts[f"{kind}|{m}|{a} vs {b}"] = res

    summary = {
        "run_dir": str(run_dir),
        "n_runs": len(graded),
        "n_excluded": len(exclude & {r["run_id"] for r in runs}),
        "models": models,
        "conditions": sorted({g.run["condition"] for g in graded}),
        "groups": by_group,
        "contrasts": contrasts,
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")
    _write_csv(graded, out / "results.csv")
    if charts:
        from .charts import write_charts

        write_charts(summary, graded, run_dir, out)
    return summary


def _clean(o: Any) -> Any:
    """Strict JSON: NaN/inf become null, tuples become lists."""
    if isinstance(o, float):
        return None if o != o or o in (float("inf"), float("-inf")) else o
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def _cell(g: Graded, f: str) -> int | str:
    """A per-run field cell: 1/0, or empty for a field the label leaves unscored (a real
    negative's ``artifact_type``), so it is never read as wrong."""
    if f == "artifact_type" and g.label["artifact_type"] == UNSCORED_TYPE:
        return ""
    return int(g.correct.get(f, False))


def _write_csv(graded: list[Graded], path: Path) -> None:
    cols = [
        "run_id",
        "item",
        "case_type",
        "condition",
        "model",
        "repeat",
        "status",
        "error",
        *[f"{f}_correct" for f in SCORED_FIELDS],
        "joint_correct",
        "joint_lenient_correct",
        "confidence",
    ]
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for g in graded:
            w.writerow(
                [
                    g.run["run_id"],
                    g.run["item"],
                    g.label["case_type"],
                    g.run["condition"],
                    g.run["model"],
                    g.run["repeat"],
                    g.status,
                    g.error,
                    *[_cell(g, fl) for fl in SCORED_FIELDS],
                    int(g.joint),
                    int(g.joint_lenient),
                    g.answer["confidence"] if g.answer else "",
                ]
            )


# ---------------------------------------------------------------------------
# Round 4 (docs/eval/visual-design-r4.md sections 5 and 8)
# ---------------------------------------------------------------------------

BASELINE_R4 = "Traw-r4"
PRIMARY_CONTRASTS_R4: dict[str, tuple[str, str]] = {
    "C1": ("Traw-r4+Tcur", "Traw-r4"),
    "C2": ("V+Traw-r4", "Traw-r4"),
}
# The focused round (docs/eval/visual-preregistration-r4.md): same contrasts on the raw-input
# (samtools mpileup) baseline. A run set uses whichever family its baseline arm belongs to.
PRIMARY_CONTRASTS_R4_TMP: dict[str, tuple[str, str]] = {
    "C1": ("Tmp+Tcur", "Tmp"),
    "C2": ("V+Tmp", "Tmp"),
}
# stratum -> (pair, polarity); a label's own ``pair``/``polarity`` fields win when present.
STRATUM_PAIRS_R4: dict[str, tuple[str, str]] = {
    "mosaic_snv_pos": ("mosaic", "pos"),
    "mosaic_caller_fp_neg": ("mosaic", "neg"),
    "tr_indel_pos": ("tr", "pos"),
    "tr_stutter_neg": ("tr", "neg"),
    "segdup_cmrg_pos": ("segdup", "pos"),
    "segdup_psv_neg": ("segdup", "neg"),
    "indel_adjacent_misalign_neg": ("indel_adj", "neg"),
    "indel_adjacent_snv_pos": ("indel_adj", "pos"),
    "anchor_het_snv": ("anchor", "pos"),
    "anchor_hom_snv": ("anchor", "pos"),
    "anchor_clean_neg": ("anchor", "neg"),
}
PAIRS_R4 = ("mosaic", "tr", "segdup", "indel_adj")  # the primary's pairs; anchors excluded
HARD_STRATA_R4 = tuple(s for s, (p, _) in STRATUM_PAIRS_R4.items() if p in PAIRS_R4)
GUARDRAIL_MAX_DROP = 0.05
TOST_MARGIN = 0.07
ALPHA = 0.05
UNCERTAIN = "uncertain"
CONFIDENT_WRONG_AT = 0.8
ALLELE_FLANK = 60  # bp of reference either side of the POI searched for an equivalent indel
_EDGE = 10  # an equivalent placement must sit this far inside the fetched window

# Section 8's pre-specified directions. Only exact "+" / "-" cells enter the sign test;
# "0", "0/+", "+/0", "0/-" are not counted.
PREDICTED_DIRECTIONS_R4: dict[str, dict[str, str]] = {
    "C1": {
        "mosaic_snv_pos": "0/-",
        "mosaic_caller_fp_neg": "+",
        "tr_indel_pos": "0",
        "tr_stutter_neg": "+",
        "segdup_cmrg_pos": "-",
        "segdup_psv_neg": "+/0",
        "indel_adjacent_misalign_neg": "+",
        "indel_adjacent_snv_pos": "0/-",
    },
    "C2": {
        "mosaic_snv_pos": "+/0",
        "mosaic_caller_fp_neg": "+",
        "tr_indel_pos": "0",
        "tr_stutter_neg": "0",
        "segdup_cmrg_pos": "0/+",
        "segdup_psv_neg": "+",
        "indel_adjacent_misalign_neg": "+",
        "indel_adjacent_snv_pos": "0/-",
    },
}

RefFetch = Callable[[str, int, int], str]  # (contig, start0, end0) -> reference bases
ScoreFn = Callable[[Graded], float]


def stratum_of(label: dict[str, Any]) -> str:
    return str(label.get("stratum") or label["case_type"])


def pair_of(label: dict[str, Any]) -> tuple[str, str] | None:
    """(pair, polarity) of an item's label; ``None`` for items outside every pair (e.g. the
    pilot's round-3 caption-check items)."""
    if label.get("pair") and label.get("polarity"):
        return str(label["pair"]), str(label["polarity"])
    return STRATUM_PAIRS_R4.get(stratum_of(label))


def joint_score(g: Graded) -> float:
    """Strict joint correctness (``uncertain`` and invalid answers are wrong)."""
    return float(g.joint)


def run_utility(g: Graded, wrong: float = -1.0) -> float:
    """+1 joint-correct, 0 for an ``uncertain`` verdict, ``wrong`` otherwise (invalid too)."""
    if g.joint:
        return 1.0
    if g.answer is not None and g.answer["artifact_verdict"] == UNCERTAIN:
        return 0.0
    return wrong


def _utility(wrong: float) -> ScoreFn:
    return lambda g: run_utility(g, wrong)


# --- alleles ----------------------------------------------------------------


def _apply_indel(ref: str, p: int, allele: str, vtype: str) -> str | None:
    """The haplotype after a deletion of ``allele`` starting at ``p`` or an insertion of
    ``allele`` after anchor ``p``; ``None`` when a deletion does not match the reference."""
    if vtype == "del":
        if ref[p : p + len(allele)] != allele:
            return None
        return ref[:p] + ref[p + len(allele) :]
    return ref[: p + 1] + allele + ref[p + 1 :]


def allele_equivalent(
    answer_allele: Any,
    truth_allele: Any,
    variant_type: str,
    *,
    contig: str,
    pos0: int,
    ref_fetch: RefFetch | None,
) -> bool:
    """Whether an answered allele describes the truth allele (design section 5, genotype).

    Equal after :func:`~.answers.normalize_allele` always counts. For a single insertion or
    deletion, the answer also counts when placing it anywhere near the POI yields the same
    haplotype as the truth — left-normalisation and repeat-shift equivalence in one rule
    (``A`` deleted from ``AAAA``, or ``CA`` for ``AC`` in ``ACACAC``). ``pos0`` is the 0-based
    POI (a deletion's first deleted base, an insertion's anchor), ``variant_type`` the truth
    type; without ``ref_fetch`` only the string rule applies.
    """
    a, t = normalize_allele(answer_allele), normalize_allele(truth_allele)
    if a == t:
        return True
    if variant_type not in ("ins", "del") or ref_fetch is None or not a or not t:
        return False
    if "," in a or "," in t:
        return False
    lo = max(0, pos0 - ALLELE_FLANK)
    ref = ref_fetch(contig, lo, pos0 + ALLELE_FLANK + len(t)).upper()
    truth_hap = _apply_indel(ref, pos0 - lo, t, variant_type)
    if truth_hap is None:
        return False
    first = min(_EDGE, pos0 - lo)
    for p in range(first, len(ref) - len(a) - _EDGE):
        if _apply_indel(ref, p, a, variant_type) == truth_hap:
            return True
    return False


def genotype_correct(g: Graded, ref_fetch: RefFetch | None = None) -> bool:
    """Joint + zygosity + allele (equivalence-aware) + ``artifact_type`` where scored."""
    if not g.joint or g.answer is None:
        return False
    lab, ans = g.label, g.answer
    if ans["zygosity"] != lab["zygosity"]:
        return False
    if lab["variant_present"] == "yes":
        if ans["variant_type"] != lab["variant_type"]:
            return False
        pos0 = int(lab["positions"][0]) - 1 if lab.get("positions") else -1
        if not allele_equivalent(
            ans["alt_allele"],
            lab["alt_allele"],
            lab["variant_type"],
            contig=str(lab.get("contig", "")),
            pos0=pos0,
            ref_fetch=ref_fetch if pos0 >= 0 else None,
        ):
            return False
    if _applicable(g, "artifact_type"):
        return bool(ans["artifact_type"] == lab["artifact_type"])
    return True


def fasta_ref_fetch(path: str | Path) -> RefFetch:
    """A :data:`RefFetch` over an indexed FASTA (the run set's blinded reference)."""
    import pysam

    fa = pysam.FastaFile(str(path))

    def fetch(contig: str, start: int, end: int) -> str:
        n = fa.get_reference_length(contig)
        return str(fa.fetch(contig, max(0, start), min(n, end)))

    return fetch


# --- balanced aggregation -----------------------------------------------------


def _item_scores(rows: list[Graded], score: ScoreFn) -> dict[str, dict[str, float]]:
    """item -> model -> score (mean over repeats)."""
    acc: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for g in rows:
        acc[g.run["item"]][g.run["model"]].append(score(g))
    return {i: {m: mean(v) for m, v in ms.items()} for i, ms in acc.items()}


def _labels(rows: list[Graded]) -> dict[str, dict[str, Any]]:
    return {g.run["item"]: g.label for g in rows}


def balanced_summary(rows: list[Graded], score: ScoreFn = joint_score) -> dict[str, Any]:
    """Section 5's aggregation over one set of runs (items' scores averaged over models).

    Per pair: mean of positive-class and negative-class item accuracy; primary: mean over
    the hard pairs that have both classes. Sensitivity / specificity: mean over the hard
    positive / negative strata of stratum accuracy.
    """
    labels = _labels(rows)
    items = {i: mean(list(ms.values())) for i, ms in _item_scores(rows, score).items()}
    by_stratum: dict[str, list[float]] = defaultdict(list)
    by_class: dict[tuple[str, str], list[float]] = defaultdict(list)
    for i, v in items.items():
        by_stratum[stratum_of(labels[i])].append(v)
        pp = pair_of(labels[i])
        if pp is not None:
            by_class[pp].append(v)
    per_pair: dict[str, dict[str, float | None]] = {}
    for pair in sorted({p for p, _ in by_class}):
        pos, neg = by_class.get((pair, "pos")), by_class.get((pair, "neg"))
        per_pair[pair] = {
            "pos": mean(pos) if pos else None,
            "neg": mean(neg) if neg else None,
            "balanced": (mean(pos) + mean(neg)) / 2 if pos and neg else None,
        }
    hard = [per_pair[p]["balanced"] for p in PAIRS_R4 if p in per_pair]
    hard_vals = [x for x in hard if x is not None]
    stratum_acc = {s: mean(v) for s, v in sorted(by_stratum.items())}
    pos_s = [stratum_acc[s] for s in HARD_STRATA_R4 if s in stratum_acc and _polarity(s) == "pos"]
    neg_s = [stratum_acc[s] for s in HARD_STRATA_R4 if s in stratum_acc and _polarity(s) == "neg"]
    return {
        "primary": mean(hard_vals) if hard_vals else float("nan"),
        "n_pairs": len(hard_vals),
        "per_pair": per_pair,
        "sensitivity": mean(pos_s) if pos_s else float("nan"),
        "specificity": mean(neg_s) if neg_s else float("nan"),
        "per_stratum": stratum_acc,
    }


def _polarity(stratum: str) -> str:
    return STRATUM_PAIRS_R4[stratum][1]


def model_balanced(rows: list[Graded], score: ScoreFn = joint_score) -> dict[str, Any]:
    """:func:`balanced_summary` per model, then equal-weighted over models (section 5)."""
    models = sorted({g.run["model"] for g in rows})
    per_model = {
        m: balanced_summary([g for g in rows if g.run["model"] == m], score) for m in models
    }

    def avg(get: Callable[[dict[str, Any]], Any]) -> float | None:
        vals = [v for s in per_model.values() if (v := get(s)) is not None and v == v]
        return mean(vals) if vals else None

    pairs = sorted({p for s in per_model.values() for p in s["per_pair"]})
    strata = sorted({k for s in per_model.values() for k in s["per_stratum"]})
    return {
        "primary": avg(lambda s: s["primary"]),
        "sensitivity": avg(lambda s: s["sensitivity"]),
        "specificity": avg(lambda s: s["specificity"]),
        "per_pair": {
            p: {k: avg(_getter("per_pair", p, k)) for k in ("pos", "neg", "balanced")}
            for p in pairs
        },
        "per_stratum": {k: avg(_getter("per_stratum", k)) for k in strata},
        "per_model": {m: s["primary"] for m, s in per_model.items()},
        "per_model_stratum": {m: s["per_stratum"] for m, s in per_model.items()},
    }


def _getter(*path: str) -> Callable[[dict[str, Any]], Any]:
    def get(s: dict[str, Any]) -> Any:
        cur: Any = s
        for k in path:
            cur = cur.get(k) if isinstance(cur, dict) else None
        return cur

    return get


def balanced_weights(labels: dict[str, dict[str, Any]]) -> dict[str, float]:
    """Section 8's item weight 1 / (2 · n_class · n_pairs) over the hard pairs present with
    both classes; other items (anchors, unpaired classes) weigh 0. Weights sum to 1."""
    by_class: dict[tuple[str, str], list[str]] = defaultdict(list)
    for i, lab in labels.items():
        pp = pair_of(lab)
        if pp is not None and pp[0] in PAIRS_R4:
            by_class[pp].append(i)
    pairs = [p for p in PAIRS_R4 if by_class.get((p, "pos")) and by_class.get((p, "neg"))]
    out = dict.fromkeys(labels, 0.0)
    for p in pairs:
        for pol in ("pos", "neg"):
            members = by_class[(p, pol)]
            for i in members:
                out[i] = 1.0 / (2 * len(members) * len(pairs))
    return out


def _paired_diffs(ra: list[Graded], rb: list[Graded], score: ScoreFn) -> dict[str, float]:
    """item -> model-averaged (A - B) over the models answered under both conditions."""
    sa, sb = _item_scores(ra, score), _item_scores(rb, score)
    out: dict[str, float] = {}
    for i in sorted(set(sa) & set(sb)):
        models = sorted(set(sa[i]) & set(sb[i]))
        if models:
            out[i] = mean([sa[i][m] - sb[i][m] for m in models])
    return out


def _signflip_block(
    diffs: dict[str, float],
    labels: dict[str, dict[str, Any]],
    *,
    n_flips: int,
    seed: int,
) -> dict[str, Any]:
    w = balanced_weights({i: labels[i] for i in diffs})
    keys = [i for i in sorted(diffs) if w[i] > 0]
    if not keys:
        return {"delta": None, "p": None, "ci95": [None, None], "n_items": 0}
    d = [diffs[i] for i in keys]
    ws = [w[i] for i in keys]
    cl = [str(labels[i].get("cluster") or i) for i in keys]
    weighted = [x * y for x, y in zip(d, ws, strict=True)]
    lo, hi = signflip_ci(d, ws, cl, n_flips=n_flips, seed=seed)
    return {
        "delta": sum(weighted),
        "p": cluster_signflip_test(weighted, cl, n_flips=n_flips, seed=seed),
        "ci95": [lo, hi],
        "n_items": len(keys),
        "n_clusters": len(set(cl)),
        "design_effect": design_effect(weighted, cl),
        "_inputs": (d, ws, cl),
    }


def contrast_r4(
    rows: list[Graded],
    a: str,
    b: str,
    *,
    n_flips: int = N_FLIPS,
    seed: int = 0,
) -> dict[str, Any] | None:
    """One confirmatory contrast (section 8): primary Δ, sign-flip p, inverted 95% CI,
    TOST ±0.07, guardrail, utility (−1 and −2) and descriptive per-pair/stratum/model Δ.
    ``holm_p`` and ``supported`` are filled in by :func:`grade_run_set_r4`."""
    ra = [g for g in rows if g.run["condition"] == a]
    rb = [g for g in rows if g.run["condition"] == b]
    if not ra or not rb:
        return None
    labels = {**_labels(rb), **_labels(ra)}
    diffs = _paired_diffs(ra, rb, joint_score)
    primary = _signflip_block(diffs, labels, n_flips=n_flips, seed=seed)
    if primary["n_items"]:
        d, ws, cl = primary["_inputs"]
        eq = tost(d, ws, cl, margin=TOST_MARGIN, n_flips=n_flips, seed=seed)
    else:
        eq = {"p": float("nan"), "equivalent": False}
    primary.pop("_inputs", None)
    utility = {}
    for name, wrong in (("utility", -1.0), ("utility_wrong_m2", -2.0)):
        blk = _signflip_block(
            _paired_diffs(ra, rb, _utility(wrong)), labels, n_flips=n_flips, seed=seed
        )
        blk.pop("_inputs", None)
        utility[name] = blk

    paired = set(diffs)
    pa = model_balanced([g for g in ra if g.run["item"] in paired])
    pb = model_balanced([g for g in rb if g.run["item"] in paired])
    d_sens = _sub(pa["sensitivity"], pb["sensitivity"])
    d_spec = _sub(pa["specificity"], pb["specificity"])
    holds = (
        d_sens is not None
        and d_spec is not None
        and d_sens >= -GUARDRAIL_MAX_DROP
        and d_spec >= -GUARDRAIL_MAX_DROP
    )
    per_model: dict[str, Any] = {}
    idx_b = {(g.run["item"], g.run["model"], g.run["repeat"]): g.joint for g in rb}
    for m in sorted({g.run["model"] for g in ra}):
        pairs = [
            (g.joint, idx_b[k])
            for g in ra
            if g.run["model"] == m
            and (k := (g.run["item"], g.run["model"], g.run["repeat"])) in idx_b
        ]
        nb, nc, p = mcnemar_exact(pairs)
        per_model[m] = {
            "delta_primary": _sub(pa["per_model"].get(m), pb["per_model"].get(m)),
            "mcnemar": {"a_only": nb, "b_only": nc, "p": p, "n_pairs": len(pairs)},
        }
    return {
        "a": a,
        "b": b,
        **primary,
        "tost": eq,
        "guardrail": {
            "sensitivity": [pa["sensitivity"], pb["sensitivity"]],
            "specificity": [pa["specificity"], pb["specificity"]],
            "delta_sensitivity": d_sens,
            "delta_specificity": d_spec,
            "max_drop": GUARDRAIL_MAX_DROP,
            "holds": holds,
        },
        **utility,
        "per_pair_delta": {
            p: _sub(pa["per_pair"][p]["balanced"], pb["per_pair"].get(p, {}).get("balanced"))
            for p in pa["per_pair"]
        },
        "per_stratum_delta": {
            s: _sub(pa["per_stratum"][s], pb["per_stratum"].get(s)) for s in pa["per_stratum"]
        },
        "per_model": per_model,
    }


def _sub(x: float | None, y: float | None) -> float | None:
    if x is None or y is None or x != x or y != y:
        return None
    return x - y


def directional_sign_test(
    per_stratum_delta: dict[str, dict[str, float | None]],
    predictions: dict[str, dict[str, str]] = PREDICTED_DIRECTIONS_R4,
) -> dict[str, Any]:
    """Section 8's one directional statistic: how many of the 16 stratum × contrast cells
    predicted exactly "+" or "-" match the observed sign. Cells with a 0 or mixed
    prediction, a missing delta or an observed delta of exactly 0 are not counted.
    p = P(X >= matched) under X ~ Binomial(counted, 1/2) (one-sided, directions fixed)."""
    cells: list[dict[str, Any]] = []
    for c, preds in sorted(predictions.items()):
        for stratum, pred in preds.items():
            obs = per_stratum_delta.get(c, {}).get(stratum)
            counted = pred in ("+", "-") and obs is not None and obs != 0
            matched = counted and ((obs or 0) > 0) == (pred == "+")
            cells.append(
                {
                    "contrast": c,
                    "stratum": stratum,
                    "predicted": pred,
                    "observed_delta": obs,
                    "counted": counted,
                    "matched": matched,
                }
            )
    n = sum(c["counted"] for c in cells)
    k = sum(c["matched"] for c in cells)
    return {
        "n_counted": n,
        "n_matched": k,
        "p_one_sided": binom_sf(k, n) if n else None,
        "cells": cells,
    }


def summarize_group_r4(rows: list[Graded], ref_fetch: RefFetch | None = None) -> dict[str, Any]:
    """Every section-5 metric for one set of runs (one model or all models × one condition)."""
    valid = [g for g in rows if g.status == "valid" and g.answer]
    conf = [float(g.answer["confidence"]) for g in valid if g.answer]
    corr = [g.joint for g in valid]
    pooled, lo, hi = cluster_bootstrap_mean(_clusters(rows, lambda g: g.joint))
    by_source: dict[str, list[float]] = defaultdict(list)
    for g in rows:
        by_source[str(g.label.get("truth_source") or "unknown")].append(float(g.joint))
    hard = [g for g in rows if stratum_of(g.label) in HARD_STRATA_R4]
    stutter = [g for g in rows if stratum_of(g.label) == "tr_stutter_neg"]
    bal = model_balanced(rows)
    return {
        "n_runs": len(rows),
        "n_items": len({g.run["item"] for g in rows}),
        "n_valid": len(valid),
        "n_invalid": sum(g.status == "invalid" for g in rows),
        "n_missing": sum(g.status == "missing" for g in rows),
        "invalid_rate": mean([float(g.status != "valid") for g in rows]),
        **bal,
        "utility": model_balanced(rows, _utility(-1.0))["primary"],
        "utility_wrong_m2": model_balanced(rows, _utility(-2.0))["primary"],
        "pooled_joint_accuracy": pooled,
        "pooled_joint_ci95": [lo, hi],
        "abstention_rate": mean(
            [float(bool(g.answer and g.answer["artifact_verdict"] == UNCERTAIN)) for g in rows]
        ),
        "genotype_accuracy": mean([float(genotype_correct(g, ref_fetch)) for g in hard])
        if hard
        else None,
        "artifact_type_accuracy_tr_stutter": mean(
            [
                float(bool(g.answer and g.answer["artifact_type"] == g.label["artifact_type"]))
                for g in stutter
            ]
        )
        if stutter
        else None,
        "accuracy_by_truth_source": {k: mean(v) for k, v in sorted(by_source.items())},
        "calibration": {
            "brier": brier(conf, corr) if conf else None,
            "ece5": ece5(conf, corr) if conf else None,
            "auroc_confidence": auroc(conf, corr) if conf else None,
            "confident_wrong_rate": confident_wrong_rate(conf, corr, CONFIDENT_WRONG_AT)
            if conf
            else None,
            "selective_accuracy": {
                str(c): selective_accuracy(conf, corr, c)[0] if conf else None
                for c in SELECTIVE_COVERAGES
            },
            "aurc": aurc(conf, corr) if conf else None,
        },
    }


def grade_runs_r4(
    run_dir: Path,
    *,
    answers_dir: str | Path | None = None,
    exclude: set[str] | None = None,
    exclude_items: set[str] | None = None,
    conditions: set[str] | None = None,
) -> list[Graded]:
    """Grade the runs of a round-4 run set. ``conditions`` restricts which runs are graded
    — their answer files are never opened otherwise (calibration reads B only)."""
    runs = load_jsonl(run_dir / "runs.jsonl")
    labels = {r["item"]: r for r in load_jsonl(run_dir / "labels" / "labels.jsonl")}
    exclude, exclude_items = exclude or set(), exclude_items or set()
    out: list[Graded] = []
    for run in runs:
        if run["run_id"] in exclude or run["item"] in exclude_items:
            continue
        if conditions is not None and run["condition"] not in conditions:
            continue
        if answers_dir is not None:
            answer_file = Path(answers_dir) / Path(run["answer_path"]).name
        else:
            answer_file = run_dir / run["answer_path"]
        out.append(grade_run(run, labels[run["item"]], answer_file))
    return out


def genotype_contrast(
    rows: list[Graded], a: str, b: str, ref_fetch: RefFetch | None = None
) -> dict[str, Any] | None:
    """Genotype endpoint between arms ``a`` and ``b``: among true variants (positive
    labels), the item × model × repeat pairs where each arm is :func:`genotype_correct`,
    compared by an exact two-sided McNemar test. ``None`` when the arms share no pair."""
    idx_b = {
        (g.run["item"], g.run["model"], g.run["repeat"]): g
        for g in rows
        if g.run["condition"] == b and g.label.get("variant_present") == "yes"
    }
    pairs = [
        (genotype_correct(g, ref_fetch), genotype_correct(idx_b[k], ref_fetch))
        for g in rows
        if g.run["condition"] == a
        and (k := (g.run["item"], g.run["model"], g.run["repeat"])) in idx_b
    ]
    if not pairs:
        return None
    nb, nc, p = mcnemar_exact(pairs)
    return {
        "a": a,
        "b": b,
        "n_pairs": len(pairs),
        "a_correct": sum(x for x, _ in pairs),
        "b_correct": sum(y for _, y in pairs),
        "a_only": nb,
        "b_only": nc,
        "p": p,
    }


def grade_run_set_r4(
    run_dir: str | Path,
    *,
    answers_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    exclude: set[str] | None = None,
    exclude_items: set[str] | None = None,
    n_flips: int = N_FLIPS,
    seed: int = 0,
    ref_fetch: RefFetch | None = None,
    genotype_confirmatory: bool = False,
) -> dict[str, Any]:
    """Round-4 grading → ``<out>/summary.json`` + ``<out>/results.csv``.

    ``exclude`` drops runs (infrastructure deviations), ``exclude_items`` drops items
    (the adjudication sensitivity analysis). The reference for allele equivalence defaults
    to the run set's ``build.json`` ``reference``. Each contrast also carries its
    :func:`genotype_contrast`; with ``genotype_confirmatory`` (the confirmation round) those
    p-values join the Holm family, and a genotype endpoint is supported when its Holm p is
    below alpha and arm ``a`` wins more discordant pairs.
    """
    run_dir = Path(run_dir)
    out = Path(out_dir) if out_dir else run_dir / "results"
    out.mkdir(parents=True, exist_ok=True)
    if ref_fetch is None:
        ref_fetch = _build_ref_fetch(run_dir)
    graded = grade_runs_r4(
        run_dir, answers_dir=answers_dir, exclude=exclude, exclude_items=exclude_items
    )
    models = sorted({g.run["model"] for g in graded})
    conds = sorted({g.run["condition"] for g in graded})
    groups: dict[str, Any] = {}
    for c in conds:
        rows = [g for g in graded if g.run["condition"] == c]
        groups[f"__all__|{c}"] = summarize_group_r4(rows, ref_fetch)
        for m in models:
            mrows = [g for g in rows if g.run["model"] == m]
            if mrows:
                groups[f"{m}|{c}"] = summarize_group_r4(mrows, ref_fetch)

    contrasts: dict[str, Any] = {}
    family = PRIMARY_CONTRASTS_R4_TMP if "Tmp" in conds else PRIMARY_CONTRASTS_R4
    for name, (a, b) in family.items():
        res = contrast_r4(graded, a, b, n_flips=n_flips, seed=seed)
        if res is not None:
            res["genotype"] = genotype_contrast(graded, a, b, ref_fetch)
            contrasts[name] = res
    names = [n for n in contrasts if contrasts[n]["p"] is not None]
    geno = (
        [n for n in names if contrasts[n]["genotype"] is not None] if genotype_confirmatory else []
    )
    family_p = [contrasts[n]["p"] for n in names] + [contrasts[n]["genotype"]["p"] for n in geno]
    adjusted = holm(family_p)
    for n, hp in zip(geno, adjusted[len(names) :], strict=True):
        gc = contrasts[n]["genotype"]
        gc["holm_p"] = hp
        gc["supported"] = bool(hp < ALPHA and gc["a_only"] > gc["b_only"])
    for n, hp in zip(names, adjusted[: len(names)], strict=True):
        c = contrasts[n]
        lo, hi = c["ci95"]
        c["holm_p"] = hp
        c["ci_excludes_0"] = lo > 0 or hi < 0
        c["supported"] = bool(
            hp < ALPHA and c["ci_excludes_0"] and c["guardrail"]["holds"] and c["delta"] > 0
        )
        if c["supported"]:
            c["verdict"] = "supported"
        elif hp < ALPHA and c["ci_excludes_0"] and c["delta"] < 0:
            c["verdict"] = "harm"
        elif c["tost"]["equivalent"]:
            c["verdict"] = "no meaningful difference (TOST)"
        else:
            c["verdict"] = "inconclusive"
    sign_test = directional_sign_test({n: c["per_stratum_delta"] for n, c in contrasts.items()})

    summary = {
        "run_dir": str(run_dir),
        "round": 4,
        "n_runs": len(graded),
        "n_excluded_runs": len(exclude or set()),
        "n_excluded_items": len(exclude_items or set()),
        "models": models,
        "conditions": conds,
        "groups": groups,
        "contrasts": contrasts,
        "directional_sign_test": sign_test,
        "tables": {
            "per_stratum_model": {c: groups[f"__all__|{c}"]["per_model_stratum"] for c in conds},
        },
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")
    _write_csv_r4(graded, out / "results.csv", ref_fetch)
    return summary


def _build_ref_fetch(run_dir: Path) -> RefFetch | None:
    meta = run_dir / "build.json"
    if not meta.exists():
        return None
    ref = json.loads(meta.read_text()).get("reference")
    if not ref:
        return None
    # The builder records the path as it was given (relative to its working directory when
    # --out was relative); fall back to the run directory for a moved run set.
    for path in (Path(ref), run_dir / ref):
        if path.exists() and Path(str(path) + ".fai").exists():
            return fasta_ref_fetch(path)
    return None


def _write_csv_r4(graded: list[Graded], path: Path, ref_fetch: RefFetch | None) -> None:
    cols = [
        "run_id",
        "item",
        "stratum",
        "pair",
        "polarity",
        "cluster",
        "truth_source",
        "knob",
        "condition",
        "model",
        "repeat",
        "status",
        "error",
        *[f"{f}_correct" for f in SCORED_FIELDS],
        "joint_correct",
        "utility",
        "genotype_correct",
        "confidence",
    ]
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for g in graded:
            pp = pair_of(g.label) or ("", "")
            w.writerow(
                [
                    g.run["run_id"],
                    g.run["item"],
                    stratum_of(g.label),
                    pp[0],
                    pp[1],
                    g.label.get("cluster", ""),
                    g.label.get("truth_source", ""),
                    g.label.get("knob", ""),
                    g.run["condition"],
                    g.run["model"],
                    g.run["repeat"],
                    g.status,
                    g.error,
                    *[_cell(g, fl) for fl in SCORED_FIELDS],
                    int(g.joint),
                    run_utility(g),
                    int(genotype_correct(g, ref_fetch)),
                    g.answer["confidence"] if g.answer else "",
                ]
            )

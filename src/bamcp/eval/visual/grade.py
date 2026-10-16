"""Grade a directory of model answers against the run set's labels.

Inputs: ``<run_dir>/runs.jsonl``, ``<run_dir>/labels/labels.jsonl`` and one answer file per
run (``answer_path`` in the manifest, relative to the run dir or ``--answers-dir``).
Outputs (``<run_dir>/results`` by default): ``results.csv`` (one row per run),
``summary.json`` (all metrics) and charts (see :mod:`.charts`).

Scoring rules are pre-registered in ``docs/eval/visual-preregistration.md``: an invalid or
missing answer is wrong on every field; excluded runs (``--exclude``) are dropped.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .answers import SCORED_FIELDS, field_correct, parse_answer
from .stats import (
    brier,
    cluster_bootstrap_mean,
    ece,
    f1_scores,
    mcnemar_exact,
    mean,
    paired_cluster_bootstrap,
    reliability_bins,
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

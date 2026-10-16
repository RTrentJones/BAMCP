"""Grader self-tests: oracle answers must score 100%, random answers ≈ chance."""

from __future__ import annotations

import json
import math
import random
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .answers import ENUMS
from .grade import UNSCORED_TYPE, grade_run_set, load_jsonl

# P(joint correct) for a uniform random answer: 1/2 (variant_present) x 1/3 (verdict).
RANDOM_CHANCE = 1 / (len(ENUMS["variant_present"]) * len(ENUMS["artifact_verdict"]))


@dataclass
class SelfTestResult:
    name: str
    ok: bool
    detail: str


def oracle_answer(label: dict[str, Any]) -> dict[str, Any]:
    return {
        "variant_present": label["variant_present"],
        "variant_type": label["variant_type"],
        "alt_allele": label["alt_allele"],
        "zygosity": label["zygosity"],
        "artifact_verdict": label["artifact_verdict"],
        # An unscored type (real-tier negatives) is not a valid answer value: answer "none".
        "artifact_type": "none"
        if label["artifact_type"] == UNSCORED_TYPE
        else label["artifact_type"],
        "phase": label["phase"],
        "confidence": 1.0,
        "rationale": "oracle",
    }


def random_answer(rng: random.Random) -> dict[str, Any]:
    ans: dict[str, Any] = {k: rng.choice(v) for k, v in ENUMS.items()}
    ans["alt_allele"] = rng.choice(["", "A", "C", "G", "T"])
    ans["confidence"] = round(rng.random(), 3)
    ans["rationale"] = "random"
    return ans


def write_answers(run_dir: Path, dest: Path, kind: str, seed: int = 0) -> int:
    runs = load_jsonl(run_dir / "runs.jsonl")
    labels = {r["item"]: r for r in load_jsonl(run_dir / "labels" / "labels.jsonl")}
    shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir(parents=True)
    rng = random.Random(seed)  # noqa: S311 — reproducible resampling
    for r in runs:
        ans = oracle_answer(labels[r["item"]]) if kind == "oracle" else random_answer(rng)
        (dest / Path(r["answer_path"]).name).write_text(json.dumps(ans))
    return len(runs)


def run_selftests(
    run_dir: str | Path, *, charts: bool = True, seed: int = 0
) -> list[SelfTestResult]:
    run_dir = Path(run_dir)
    base = run_dir / "selftest"
    results: list[SelfTestResult] = []

    # Oracle: every group perfect on every scored field.
    write_answers(run_dir, base / "oracle" / "answers", "oracle")
    s = grade_run_set(
        run_dir,
        answers_dir=base / "oracle" / "answers",
        out_dir=base / "oracle" / "results",
        charts=charts,
    )
    bad = [
        k
        for k, g in s["groups"].items()
        if g["primary_joint_accuracy"] != 1.0
        or any(v not in (None, 1.0) for v in g["field_accuracy"].values())
        or g["artifact_type_macro_f1"] != 1.0
    ]
    results.append(
        SelfTestResult(
            "oracle = 100%",
            not bad and s["n_runs"] > 0,
            f"{s['n_runs']} runs; imperfect groups: {bad or 'none'}",
        )
    )

    # Random: joint accuracy within 4 SE (or 0.08) of chance.
    write_answers(run_dir, base / "random" / "answers", "random", seed=seed)
    s = grade_run_set(
        run_dir,
        answers_dir=base / "random" / "answers",
        out_dir=base / "random" / "results",
        charts=False,
    )
    n = s["n_runs"]
    rows = list(_iter_csv(base / "random" / "results" / "results.csv"))
    acc = sum(int(r["joint_correct"]) for r in rows) / max(1, len(rows))
    tol = max(0.08, 4 * math.sqrt(RANDOM_CHANCE * (1 - RANDOM_CHANCE) / max(1, n)))
    results.append(
        SelfTestResult(
            "random ≈ chance",
            abs(acc - RANDOM_CHANCE) <= tol,
            f"joint={acc:.3f} vs chance {RANDOM_CHANCE:.3f} (tol ±{tol:.3f}, n={n})",
        )
    )

    # Invalid answers count as wrong: corrupt one oracle answer, expect exactly one miss.
    answers = base / "oracle" / "answers"
    victim = sorted(answers.iterdir())[0]
    victim.write_text("not json {")
    s = grade_run_set(
        run_dir, answers_dir=answers, out_dir=base / "invalid" / "results", charts=False
    )
    n_invalid = sum(g["n_invalid"] for g in s["groups"].values())
    rows = list(_iter_csv(base / "invalid" / "results" / "results.csv"))
    misses = sum(1 - int(r["joint_correct"]) for r in rows)
    results.append(
        SelfTestResult(
            "invalid = wrong",
            n_invalid == 1 and misses == 1,
            f"invalid={n_invalid} joint misses={misses}",
        )
    )
    return results


def _iter_csv(path: Path) -> Any:
    import csv

    with path.open() as f:
        yield from csv.DictReader(f)

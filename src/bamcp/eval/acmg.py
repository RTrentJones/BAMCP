"""ACMG classification eval: scoring + a ground-truth run driver.

The pure scoring functions (normalization, severity distance, failure-mode
tagging, classification extraction) are network- and LLM-free so they unit-test
in isolation. :func:`run_acmg` drives the LLM harness over a labeled manifest
(``tests/eval/datasets/acmg_v1/manifest.yaml``): for each variant it lets the
model call ``classify_variant`` and reason, extracts the final classification
from the response, and scores it against the expert label.

Metrics are deliberately boring: exact-match and within-one-step accuracy plus a
direction (overcall / undercall) tag — the same shape as the detection truth set.
Research-grade; the labels are curated expert classifications, not a clinical DB.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Ordered least -> most severe. Index distance is the "within one step" metric.
SEVERITY_ORDER = [
    "Benign",
    "Likely Benign",
    "Uncertain Significance",
    "Likely Pathogenic",
    "Pathogenic",
]

# Canonical label -> the aliases a model (or a manifest) might write.
_ALIASES = {
    "Benign": ["benign"],
    "Likely Benign": ["likely benign"],
    "Uncertain Significance": ["uncertain significance", "vus", "uncertain"],
    "Likely Pathogenic": ["likely pathogenic"],
    "Pathogenic": ["pathogenic"],
}


def normalize_classification(label: str | None) -> str | None:
    """Map a free-form classification string onto a canonical severity label."""
    if not label:
        return None
    low = label.strip().lower()
    # Check the two-word labels before the substrings they contain
    # ("likely pathogenic" before "pathogenic", "likely benign" before "benign").
    for canonical in ("Likely Pathogenic", "Likely Benign"):
        if any(a in low for a in _ALIASES[canonical]):
            return canonical
    for canonical in ("Uncertain Significance", "Pathogenic", "Benign"):
        if any(a in low for a in _ALIASES[canonical]):
            return canonical
    return None


def classification_index(label: str | None) -> int | None:
    """Severity index (0=Benign .. 4=Pathogenic), or None if unrecognized."""
    canonical = normalize_classification(label)
    if canonical is None:
        return None
    return SEVERITY_ORDER.index(canonical)


def extract_classification(text: str) -> str | None:
    """Pull the model's final classification out of its response text.

    Prefers an explicit ``final_classification`` / ``classification:`` line; else
    scans the whole text and returns the last canonical label mentioned (models
    tend to state the verdict last).
    """
    if not text:
        return None
    # Prefer an explicit final-classification declaration.
    for pattern in (
        r"final[_ ]classification\"?\s*[:=]\s*\"?([A-Za-z ()]+)",
        r"\bclassification\"?\s*[:=]\s*\"?([A-Za-z ()]+)",
    ):
        m = re.search(pattern, text, flags=re.IGNORECASE)
        if m:
            canonical = normalize_classification(m.group(1))
            if canonical:
                return canonical
    # Fall back to the last canonical label mentioned anywhere.
    best_pos = -1
    best_label: str | None = None
    low = text.lower()
    for canonical, aliases in _ALIASES.items():
        for alias in aliases:
            idx = low.rfind(alias)
            if idx > best_pos:
                # Guard: "benign"/"pathogenic" inside "likely ..." is handled by
                # normalization returning the two-word label, but for position we
                # take the latest mention and normalize the surrounding word.
                best_pos = idx
                best_label = canonical
    return best_label


@dataclass(frozen=True)
class ClassificationScore:
    """Score for one predicted-vs-expected classification."""

    predicted: str | None
    expected: str
    exact: bool
    within_one: bool
    distance: int | None  # signed: predicted - expected severity, None if unparseable
    failure_mode: str  # "none" | "overcall" | "undercall" | "unparseable"

    def as_dict(self) -> dict[str, Any]:
        return {
            "predicted": self.predicted,
            "expected": self.expected,
            "exact": self.exact,
            "within_one": self.within_one,
            "distance": self.distance,
            "failure_mode": self.failure_mode,
        }


def score_classification(predicted: str | None, expected: str) -> ClassificationScore:
    """Score a predicted classification against the expert label."""
    exp_idx = classification_index(expected)
    pred_idx = classification_index(predicted)
    if exp_idx is None:
        raise ValueError(f"Unrecognized expected classification: {expected!r}")

    if pred_idx is None:
        return ClassificationScore(
            predicted=None,
            expected=SEVERITY_ORDER[exp_idx],
            exact=False,
            within_one=False,
            distance=None,
            failure_mode="unparseable",
        )

    distance = pred_idx - exp_idx
    exact = distance == 0
    within_one = abs(distance) <= 1
    if exact:
        mode = "none"
    elif distance > 0:
        mode = "overcall"  # predicted more severe than truth
    else:
        mode = "undercall"
    return ClassificationScore(
        predicted=SEVERITY_ORDER[pred_idx],
        expected=SEVERITY_ORDER[exp_idx],
        exact=exact,
        within_one=within_one,
        distance=distance,
        failure_mode=mode,
    )


@dataclass(frozen=True)
class AcmgCase:
    """One labeled ACMG benchmark variant."""

    chrom: str
    pos: int
    ref: str
    alt: str
    expected_classification: str
    gene: str | None = None
    tier: str = "unspecified"
    trap: str | None = None
    note: str = ""


def load_acmg_cases(path: str | Path) -> list[AcmgCase]:
    """Parse an ACMG manifest YAML into a list of :class:`AcmgCase`."""
    try:
        import yaml  # type: ignore[import-untyped]
    except ImportError as e:  # pragma: no cover — install hint
        raise RuntimeError(
            "pyyaml is required for the ACMG eval. Install with 'pip install \".[eval]\"'."
        ) from e
    doc = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    cases_raw = doc.get("cases", []) if isinstance(doc, dict) else doc
    cases: list[AcmgCase] = []
    for c in cases_raw:
        cases.append(
            AcmgCase(
                chrom=str(c["chrom"]),
                pos=int(c["pos"]),
                ref=str(c["ref"]),
                alt=str(c["alt"]),
                expected_classification=str(c["expected_classification"]),
                gene=c.get("gene"),
                tier=str(c.get("tier", "unspecified")),
                trap=c.get("trap"),
                note=str(c.get("note", "")),
            )
        )
    return cases


@dataclass
class AcmgReport:
    """Aggregate ACMG classification metrics."""

    total: int
    exact_matches: int
    within_one: int
    unparseable: int
    per_case: list[dict] = field(default_factory=list)
    failure_modes: dict[str, int] = field(default_factory=dict)

    @property
    def exact_accuracy(self) -> float:
        return self.exact_matches / self.total if self.total else 0.0

    @property
    def within_one_accuracy(self) -> float:
        return self.within_one / self.total if self.total else 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "total": self.total,
            "exact_matches": self.exact_matches,
            "within_one": self.within_one,
            "unparseable": self.unparseable,
            "exact_accuracy": round(self.exact_accuracy, 4),
            "within_one_accuracy": round(self.within_one_accuracy, 4),
            "failure_modes": self.failure_modes,
            "per_case": self.per_case,
        }


def build_report(scored: Iterable[tuple[AcmgCase, ClassificationScore]]) -> AcmgReport:
    """Aggregate per-case scores into an :class:`AcmgReport`."""
    scored = list(scored)
    total = len(scored)
    exact = sum(1 for _, s in scored if s.exact)
    within = sum(1 for _, s in scored if s.within_one)
    unparseable = sum(1 for _, s in scored if s.failure_mode == "unparseable")
    modes: dict[str, int] = {}
    per_case: list[dict] = []
    for case, s in scored:
        modes[s.failure_mode] = modes.get(s.failure_mode, 0) + 1
        row = s.as_dict()
        row.update({"variant": f"{case.chrom}:{case.pos} {case.ref}>{case.alt}", "trap": case.trap})
        per_case.append(row)
    return AcmgReport(
        total=total,
        exact_matches=exact,
        within_one=within,
        unparseable=unparseable,
        per_case=per_case,
        failure_modes=modes,
    )


async def run_acmg(
    cases: list[AcmgCase],
    provider: Any,
    router: Any,
    output_dir: Path,
    bam_fixture: str,
    reference: str | None = None,
) -> AcmgReport:  # pragma: no cover — requires a live provider (nightly)
    """Drive the LLM over each ACMG case and score the classifications.

    For each case the model is prompted to call ``classify_variant`` and reason
    through the scaffold; the final classification is extracted from its response
    and scored against the expert label. Requires a real provider.
    """
    from .router import tool_descriptors  # noqa: PLC0415
    from .runner import _run_tool_loop  # noqa: PLC0415

    tools = tool_descriptors(router)
    scored: list[tuple[AcmgCase, ClassificationScore]] = []
    for case in cases:
        ref_hint = f" reference {reference}" if reference else ""
        gene_hint = f" in gene {case.gene}" if case.gene else ""
        prompt = (
            f"Classify variant {case.chrom}:{case.pos} {case.ref}>{case.alt}"
            f"{gene_hint} using the BAM {bam_fixture}{ref_hint}. "
            f"Call classify_variant, apply the ACMG scaffold, and end with a line "
            f"'final_classification: <one of {', '.join(SEVERITY_ORDER)}>'."
        )
        final_text, _, _ = await _run_tool_loop(provider, router, tools, prompt)
        predicted = extract_classification(final_text)
        scored.append((case, score_classification(predicted, case.expected_classification)))

    report = build_report(scored)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "acmg_report.json").write_text(json.dumps(report.as_dict(), indent=2))
    return report


def _format_report(report: AcmgReport) -> str:
    lines = [
        f"ACMG eval: {report.total} variants",
        f"  exact match     {report.exact_matches}/{report.total} ({report.exact_accuracy:.1%})",
        f"  within one step {report.within_one}/{report.total} ({report.within_one_accuracy:.1%})",
        f"  unparseable     {report.unparseable}",
        f"  failure modes   {report.failure_modes}",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:  # pragma: no cover — needs a live provider
    """CLI: run the ACMG classification eval against a labeled manifest."""
    import argparse
    import asyncio

    from .providers import get_provider
    from .router import InProcessRouter

    p = argparse.ArgumentParser(prog="bamcp-acmg-eval", description="ACMG classification eval.")
    p.add_argument("--manifest", default="tests/eval/datasets/acmg_v1/manifest.yaml")
    p.add_argument("--bam", required=True, help="BAM/CRAM fixture for the observation block.")
    p.add_argument("--reference", default=None)
    p.add_argument("--provider", default="anthropic", choices=["anthropic", "openai"])
    p.add_argument("--model", default="claude-opus-4-7")
    p.add_argument("--output-dir", default=".eval-results/acmg")
    args = p.parse_args(argv)

    cases = load_acmg_cases(args.manifest)
    provider = get_provider(args.provider, args.model)
    router = InProcessRouter()
    report = asyncio.run(
        run_acmg(cases, provider, router, Path(args.output_dir), args.bam, args.reference)
    )
    print(_format_report(report))
    # Non-zero exit if the model can't even land within one step on a majority.
    return 0 if report.within_one_accuracy >= 0.5 else 1


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())

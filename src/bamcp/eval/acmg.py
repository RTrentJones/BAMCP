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


# One alternation, longest alias first: at any position the regex engine tries the
# two-word labels before the single words they contain, so "likely pathogenic" is
# consumed whole and never re-read as "pathogenic" (finditer never revisits the
# characters of a match — the shorter alias is effectively masked).
_LABEL_RE = re.compile(
    r"\b(?P<label>likely[\s_-]+pathogenic|likely[\s_-]+benign|uncertain[\s_-]+significance"
    r"|pathogenic|benign|vus|uncertain)\b",
    flags=re.IGNORECASE,
)
_CANONICAL_BY_ALIAS = {
    re.sub(r"\s+", " ", alias): canonical
    for canonical, aliases in _ALIASES.items()
    for alias in aliases
}
# A label is negated when the words just before it say so ("not pathogenic",
# "non-pathogenic", "is not considered benign", "isn't likely benign"). Only the
# few words immediately preceding the label are inspected, so a negation earlier
# in the sentence ("not sure, but it is pathogenic") does not flip the verdict.
_NEGATION_RE = re.compile(
    r"(?:\bnot|\bno|\bnon|n't|\bnever|\bneither|\bnor|\bunlikely\s+to\s+be)"
    r"(?:[\s-]+(?:be|been|being|considered|classified|classifiable|as|a|an|clearly|truly))*"
    r"[\s-]*$",
    flags=re.IGNORECASE,
)
_NEGATION_WINDOW = 40


def _label_mentions(text: str) -> list[tuple[int, int, str]]:
    """Every non-negated classification mention as ``(start, end, canonical)``, in order."""
    out: list[tuple[int, int, str]] = []
    for m in _LABEL_RE.finditer(text):
        alias = re.sub(r"[\s_-]+", " ", m.group("label").lower())
        canonical = _CANONICAL_BY_ALIAS.get(alias)
        if canonical is None:
            continue
        before = text[max(0, m.start() - _NEGATION_WINDOW) : m.start()]
        if _NEGATION_RE.search(before):
            continue
        out.append((m.start(), m.end(), canonical))
    return out


def normalize_classification(label: str | None) -> str | None:
    """Map a free-form classification string onto a canonical severity label.

    Returns the first non-negated label mentioned (two-word labels win over the single
    words they contain), or None — including for a purely negated label such as
    ``"not pathogenic"``, which names no classification.
    """
    if not label:
        return None
    mentions = _label_mentions(label)
    return mentions[0][2] if mentions else None


def classification_index(label: str | None) -> int | None:
    """Severity index (0=Benign .. 4=Pathogenic), or None if unrecognized."""
    canonical = normalize_classification(label)
    if canonical is None:
        return None
    return SEVERITY_ORDER.index(canonical)


_EXPLICIT_PATTERNS = (
    re.compile(r"final[_ ]classification\W{0,3}\s*[:=]\s*\W{0,3}([A-Za-z ()_-]+)", re.I),
    re.compile(r"\bclassification\W{0,3}\s*[:=]\s*\W{0,3}([A-Za-z ()_-]+)", re.I),
)


def extract_classification(text: str) -> str | None:
    """Pull the model's final classification out of its response text.

    Prefers the LAST explicit ``final_classification:`` declaration, then the last
    ``classification:`` line — the last one, because a response that echoes tool
    output (``"classification": "Pathogenic"``) before stating its own verdict must be
    scored on the verdict. Otherwise returns the last non-negated label mentioned
    (models state the verdict last). Mentions are matched longest-alias-first, so
    "likely benign" is never read as "benign", and negated mentions
    ("not pathogenic") are ignored.
    """
    if not text:
        return None
    for pattern in _EXPLICIT_PATTERNS:
        for m in reversed(list(pattern.finditer(text))):
            canonical = normalize_classification(m.group(1))
            if canonical:
                return canonical
    mentions = _label_mentions(text)
    return mentions[-1][2] if mentions else None


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


# The B/P anchors: expert labels at the ends of the scale. A model that cannot even
# put a clear Pathogenic / Benign variant on the correct side is not classifying.
_ANCHORS = {
    "Pathogenic": {"Pathogenic", "Likely Pathogenic"},
    "Benign": {"Benign", "Likely Benign"},
}
# The degenerate answer the gate must strictly beat: "VUS" for everything. It lands
# within one step of every LP/VUS/LB label, so a within-one floor alone rewards it.
CONSTANT_BASELINE_LABEL = "Uncertain Significance"
# Exact-match floor. A majority-exact bar, and in addition it must strictly beat the
# constant-VUS baseline (see AcmgReport.meets_floors).
DEFAULT_MIN_EXACT_ACCURACY = 0.5


@dataclass
class AcmgReport:
    """Aggregate ACMG classification metrics."""

    total: int
    exact_matches: int
    within_one: int
    unparseable: int
    per_case: list[dict] = field(default_factory=list)
    failure_modes: dict[str, int] = field(default_factory=dict)
    anchors_total: int = 0
    anchors_correct: int = 0  # B/P anchors predicted on the correct side of the scale

    @property
    def exact_accuracy(self) -> float:
        return self.exact_matches / self.total if self.total else 0.0

    @property
    def within_one_accuracy(self) -> float:
        return self.within_one / self.total if self.total else 0.0

    def meets_floors(
        self,
        baseline: AcmgReport,
        *,
        min_exact_accuracy: float = DEFAULT_MIN_EXACT_ACCURACY,
    ) -> tuple[bool, list[str]]:
        """Gate: exact accuracy >= floor AND strictly above the constant-VUS baseline,
        AND every B/P anchor on the correct side.

        Within-one accuracy is reported but deliberately NOT gated: a constant "VUS"
        answer scores within one step of every LP/VUS/LB label (5/9 on acmg_v1), so a
        within-one floor passes a model that never classifies anything.
        """
        failures: list[str] = []
        if self.total == 0:
            return False, ["no ACMG cases scored"]
        if self.exact_accuracy < min_exact_accuracy:
            failures.append(
                f"exact accuracy {self.exact_accuracy:.3f} < floor {min_exact_accuracy:.3f}"
            )
        if self.exact_accuracy <= baseline.exact_accuracy:
            failures.append(
                f"exact accuracy {self.exact_accuracy:.3f} does not beat the constant-"
                f"'{CONSTANT_BASELINE_LABEL}' baseline {baseline.exact_accuracy:.3f}"
            )
        if self.anchors_correct < self.anchors_total:
            failures.append(
                "B/P anchors on the wrong side: "
                f"{self.anchors_correct}/{self.anchors_total} correct"
            )
        return (not failures, failures)

    def as_dict(self) -> dict[str, Any]:
        return {
            "total": self.total,
            "exact_matches": self.exact_matches,
            "within_one": self.within_one,
            "unparseable": self.unparseable,
            "exact_accuracy": round(self.exact_accuracy, 4),
            "within_one_accuracy": round(self.within_one_accuracy, 4),
            "anchors_total": self.anchors_total,
            "anchors_correct": self.anchors_correct,
            "failure_modes": self.failure_modes,
            "per_case": self.per_case,
        }


def _anchor_correct(case: AcmgCase, score: ClassificationScore) -> bool | None:
    """True/False for a B/P anchor case (correct side of the scale); None otherwise."""
    expected = normalize_classification(case.expected_classification)
    if expected not in _ANCHORS:
        return None
    return score.predicted in _ANCHORS[expected]


def build_report(scored: Iterable[tuple[AcmgCase, ClassificationScore]]) -> AcmgReport:
    """Aggregate per-case scores into an :class:`AcmgReport`."""
    scored = list(scored)
    total = len(scored)
    exact = sum(1 for _, s in scored if s.exact)
    within = sum(1 for _, s in scored if s.within_one)
    unparseable = sum(1 for _, s in scored if s.failure_mode == "unparseable")
    modes: dict[str, int] = {}
    per_case: list[dict] = []
    anchors_total = anchors_correct = 0
    for case, s in scored:
        modes[s.failure_mode] = modes.get(s.failure_mode, 0) + 1
        anchor = _anchor_correct(case, s)
        if anchor is not None:
            anchors_total += 1
            anchors_correct += int(anchor)
        row = s.as_dict()
        row.update(
            {
                "variant": f"{case.chrom}:{case.pos} {case.ref}>{case.alt}",
                "trap": case.trap,
                "anchor_correct": anchor,
            }
        )
        per_case.append(row)
    return AcmgReport(
        total=total,
        exact_matches=exact,
        within_one=within,
        unparseable=unparseable,
        per_case=per_case,
        failure_modes=modes,
        anchors_total=anchors_total,
        anchors_correct=anchors_correct,
    )


def constant_baseline(
    cases: Iterable[AcmgCase], label: str = CONSTANT_BASELINE_LABEL
) -> AcmgReport:
    """Score a degenerate classifier that answers ``label`` for every case."""
    cases = list(cases)
    return build_report((c, score_classification(label, c.expected_classification)) for c in cases)


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
    baseline = constant_baseline(cases)
    passed, failures = report.meets_floors(baseline)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "acmg_report.json").write_text(
        json.dumps(
            {
                **report.as_dict(),
                "constant_vus_baseline": {
                    "exact_accuracy": round(baseline.exact_accuracy, 4),
                    "within_one_accuracy": round(baseline.within_one_accuracy, 4),
                    "anchors_correct": baseline.anchors_correct,
                },
                "passed": passed,
                "failures": failures,
            },
            indent=2,
        )
    )
    return report


def _format_report(report: AcmgReport, baseline: AcmgReport | None = None) -> str:
    lines = [
        f"ACMG eval: {report.total} variants",
        f"  exact match     {report.exact_matches}/{report.total} ({report.exact_accuracy:.1%})",
        f"  within one step {report.within_one}/{report.total} ({report.within_one_accuracy:.1%})",
        f"  B/P anchors     {report.anchors_correct}/{report.anchors_total} on the correct side",
        f"  unparseable     {report.unparseable}",
        f"  failure modes   {report.failure_modes}",
    ]
    if baseline is not None:
        lines.append(
            f"  constant-VUS baseline: exact {baseline.exact_accuracy:.1%}, "
            f"within one {baseline.within_one_accuracy:.1%}, "
            f"anchors {baseline.anchors_correct}/{baseline.anchors_total}"
        )
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
    baseline = constant_baseline(cases)
    passed, failures = report.meets_floors(baseline)
    print(_format_report(report, baseline))
    print("PASS" if passed else "FAIL")
    for f in failures:
        print(f"  - {f}")
    return 0 if passed else 1


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())

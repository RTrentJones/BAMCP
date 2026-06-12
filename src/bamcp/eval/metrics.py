"""Deterministic metrics for ground-truth eval scoring.

Pure functions with no I/O so they can be unit-tested in isolation and reused
by both the truth-set scorer (:mod:`bamcp.eval.truthset`) and any future
model-comparison reporting.

The headline metrics are intentionally boring and defensible:

- precision / recall / F1 for variant calls against a labeled truth set,
- artifact-type recall (did the curation tool surface the expected risk type
  at each known-artifact site), and
- clean-site discrimination (is a true clean variant scored below the
  artifact-prone sites).

Numbers, not pass/fail strings — the caller decides the floors.
"""

from __future__ import annotations

from collections.abc import Hashable, Iterable
from dataclasses import dataclass


@dataclass(frozen=True)
class PRF:
    """Precision / recall / F1 derived from a confusion count."""

    tp: int
    fp: int
    fn: int

    @property
    def precision(self) -> float:
        denom = self.tp + self.fp
        return self.tp / denom if denom else 1.0

    @property
    def recall(self) -> float:
        denom = self.tp + self.fn
        return self.tp / denom if denom else 1.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) else 0.0

    def as_dict(self) -> dict[str, float | int]:
        return {
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "precision": round(self.precision, 4),
            "recall": round(self.recall, 4),
            "f1": round(self.f1, 4),
        }


def prf_from_sets(expected: Iterable[Hashable], predicted: Iterable[Hashable]) -> PRF:
    """Compute a :class:`PRF` by comparing two sets of labeled calls.

    Each element is an opaque, hashable key (e.g. ``("chr1", 1050, "A", "T")``).
    True positives are the intersection; false positives are predicted-only;
    false negatives are expected-only.
    """
    exp = set(expected)
    pred = set(predicted)
    tp = len(exp & pred)
    fp = len(pred - exp)
    fn = len(exp - pred)
    return PRF(tp=tp, fp=fp, fn=fn)


@dataclass(frozen=True)
class ArtifactCheck:
    """Whether the expected artifact risk type was surfaced at one site."""

    site: str
    expected_type: str
    found: bool
    risk_score: float
    likelihood: str

    def as_dict(self) -> dict[str, object]:
        return {
            "site": self.site,
            "expected_type": self.expected_type,
            "found": self.found,
            "risk_score": self.risk_score,
            "likelihood": self.likelihood,
        }


def artifact_recall(checks: Iterable[ArtifactCheck]) -> float:
    """Fraction of known-artifact sites whose expected risk type was surfaced."""
    checks = list(checks)
    if not checks:
        return 1.0
    return sum(1 for c in checks if c.found) / len(checks)


def risk_types_present(artifact_assessment: object) -> set[str]:
    """Pull the set of risk ``type`` strings from a curation artifact block.

    Tolerant of missing/oddly-shaped payloads so a malformed tool response
    degrades to "no risks found" rather than raising.
    """
    if not isinstance(artifact_assessment, dict):
        return set()
    risks = artifact_assessment.get("risks")
    if not isinstance(risks, list):
        return set()
    out: set[str] = set()
    for risk in risks:
        if isinstance(risk, dict) and isinstance(risk.get("type"), str):
            out.add(risk["type"])
    return out

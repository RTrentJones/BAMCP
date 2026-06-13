"""Ground-truth scorer: run BAMCP tools over a labeled truth set, emit metrics.

This is the deterministic backbone of the eval harness — no LLM, no network.
It routes ``get_variants`` and ``get_variant_curation_summary`` through the
in-process router against a manifest of *known* variants and artifacts (see
``tests/eval/datasets/synthetic_v1/manifest.yaml``) and computes:

- variant-detection precision/recall/F1 over a labeled region,
- false-positive control over reference-only negative regions, and
- artifact-type recall over known-artifact sites.

The output (:class:`TruthsetReport`) is plain data with a ``meets_floors``
check, so it can back a CI regression gate, a pytest assertion, or a
human-readable CLI report interchangeably.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..config import BAMCPConfig
from .metrics import PRF, ArtifactCheck, artifact_recall, prf_from_sets, risk_types_present
from .router import InProcessRouter, ToolRouter


@dataclass(frozen=True)
class TruthSite:
    """One labeled variant allele in the truth set."""

    chrom: str
    pos: int
    ref: str
    alt: str
    category: str = "uncategorized"
    expected_vaf: float | None = None
    expected_artifact: str | None = None
    expected_likelihood: str | None = None
    is_clean_control: bool = False
    # When set, the curation tool's reported confidence must equal this. Used
    # for the high-confidence positive control that makes the overconfidence
    # guard discriminating.
    expected_confidence: str | None = None

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref.upper(), self.alt.upper())


@dataclass(frozen=True)
class NegativeRegion:
    """A reference-only region that must yield zero variant calls."""

    region: str
    note: str = ""


@dataclass(frozen=True)
class TruthSet:
    """A parsed, versioned ground-truth dataset."""

    dataset: str
    version: int
    detection_region: str
    bam: str
    reference: str
    variant_sites: tuple[TruthSite, ...]
    negative_regions: tuple[NegativeRegion, ...]


def load_truthset(path: str | Path) -> TruthSet:
    """Parse a truth-set manifest YAML into a :class:`TruthSet`."""
    try:
        import yaml  # type: ignore[import-untyped]
    except ImportError as e:  # pragma: no cover — install hint
        raise RuntimeError(
            "pyyaml is required for the truth-set scorer. Install with 'pip install \".[eval]\"'."
        ) from e

    doc = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(doc, dict):
        raise ValueError("Truth-set manifest must be a mapping")

    provenance = doc.get("provenance") or {}
    sites = tuple(_parse_site(s) for s in doc.get("variant_sites", []) or [])
    negatives = tuple(
        NegativeRegion(region=str(n["region"]), note=str(n.get("note", "")))
        for n in doc.get("negative_regions", []) or []
    )
    return TruthSet(
        dataset=str(doc.get("dataset", "unnamed")),
        version=int(doc.get("version", 0)),
        detection_region=str(doc["detection_region"]),
        bam=str(provenance.get("bam", "")),
        reference=str(provenance.get("reference", "")),
        variant_sites=sites,
        negative_regions=negatives,
    )


def _parse_site(raw: dict[str, Any]) -> TruthSite:
    return TruthSite(
        chrom=str(raw["chrom"]),
        pos=int(raw["pos"]),
        ref=str(raw["ref"]),
        alt=str(raw["alt"]),
        category=str(raw.get("category", "uncategorized")),
        expected_vaf=_opt_float(raw.get("expected_vaf")),
        expected_artifact=_opt_str(raw.get("expected_artifact")),
        expected_likelihood=_opt_str(raw.get("expected_likelihood")),
        is_clean_control=bool(raw.get("is_clean_control", False)),
        expected_confidence=_opt_str(raw.get("expected_confidence")),
    )


def _opt_float(v: Any) -> float | None:
    return float(v) if isinstance(v, (int, float)) else None


def _opt_str(v: Any) -> str | None:
    return str(v) if isinstance(v, str) and v else None


@dataclass
class TruthsetReport:
    """The full deterministic scoring of a truth set."""

    dataset: str
    version: int
    detection: PRF
    false_positives: int  # total variant calls across negative regions
    artifact_checks: list[ArtifactCheck] = field(default_factory=list)
    clean_discriminated: bool = True  # clean control scored below artifact sites
    missing_calls: list[str] = field(default_factory=list)
    spurious_calls: list[str] = field(default_factory=list)
    # Safety guard: sites the curation tool flagged as high artifact risk yet
    # still reported as high confidence — a clinically dangerous contradiction.
    overconfident_sites: list[str] = field(default_factory=list)
    # Sites whose reported confidence did not match expected_confidence. The
    # high-confidence positive control lives here: if it stops reaching high
    # confidence, the overconfidence guard would pass vacuously, so this fails.
    confidence_mismatches: list[str] = field(default_factory=list)

    @property
    def artifact_recall(self) -> float:
        return artifact_recall(self.artifact_checks)

    def meets_floors(
        self,
        *,
        min_recall: float = 1.0,
        min_precision: float = 0.95,
        min_artifact_recall: float = 1.0,
        max_false_positives: int = 0,
        max_overconfident: int = 0,
        max_confidence_mismatches: int = 0,
    ) -> tuple[bool, list[str]]:
        """Return (passed, reasons-for-failure) against the given floors."""
        failures: list[str] = []
        if self.detection.recall < min_recall:
            failures.append(
                f"variant recall {self.detection.recall:.3f} < {min_recall:.3f} "
                f"(missing: {self.missing_calls})"
            )
        if self.detection.precision < min_precision:
            failures.append(
                f"variant precision {self.detection.precision:.3f} < {min_precision:.3f} "
                f"(spurious: {self.spurious_calls})"
            )
        if self.false_positives > max_false_positives:
            failures.append(
                f"{self.false_positives} false positives in negative regions "
                f"> {max_false_positives}"
            )
        if self.artifact_recall < min_artifact_recall:
            failures.append(
                f"artifact-type recall {self.artifact_recall:.3f} < {min_artifact_recall:.3f}"
            )
        if not self.clean_discriminated:
            failures.append("clean control not discriminated below artifact sites")
        if len(self.overconfident_sites) > max_overconfident:
            failures.append(
                "SAFETY: high-artifact sites reported as high confidence: "
                f"{self.overconfident_sites}"
            )
        if len(self.confidence_mismatches) > max_confidence_mismatches:
            failures.append(f"confidence calibration off: {self.confidence_mismatches}")
        return (not failures, failures)

    def as_dict(self) -> dict[str, Any]:
        return {
            "dataset": self.dataset,
            "version": self.version,
            "detection": self.detection.as_dict(),
            "false_positives": self.false_positives,
            "artifact_recall": round(self.artifact_recall, 4),
            "artifact_checks": [c.as_dict() for c in self.artifact_checks],
            "clean_discriminated": self.clean_discriminated,
            "missing_calls": self.missing_calls,
            "spurious_calls": self.spurious_calls,
            "overconfident_sites": self.overconfident_sites,
            "confidence_mismatches": self.confidence_mismatches,
        }


def _parse_variants(text: str) -> list[dict[str, Any]]:
    """Pull the variant list out of a get_variants tool response."""
    try:
        data = json.loads(text)
    except (ValueError, TypeError):
        return []
    variants = data.get("variants") if isinstance(data, dict) else None
    return variants if isinstance(variants, list) else []


def _predicted_keys(variants: Iterable[dict[str, Any]]) -> set[tuple[str, int, str, str]]:
    out: set[tuple[str, int, str, str]] = set()
    for v in variants:
        try:
            out.add(
                (
                    str(v["contig"]),
                    int(v["position"]),
                    str(v["ref"]).upper(),
                    str(v["alt"]).upper(),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    return out


async def score_truthset(
    truthset: TruthSet,
    config: BAMCPConfig,
    router: ToolRouter | None = None,
) -> TruthsetReport:
    """Run the BAMCP tools over ``truthset`` and compute deterministic metrics."""
    router = router or InProcessRouter(config)
    bam = truthset.bam
    reference = truthset.reference

    # ── Variant detection over the labeled region ────────────────────────
    res = await router.call(
        "get_variants",
        {"file_path": bam, "region": truthset.detection_region, "reference": reference},
    )
    predicted = _predicted_keys(_parse_variants(res.text)) if res.ok else set()
    expected = {s.key for s in truthset.variant_sites}
    detection = prf_from_sets(expected, predicted)
    missing = sorted(f"{c}:{p}:{r}>{a}" for (c, p, r, a) in (expected - predicted))
    spurious = sorted(f"{c}:{p}:{r}>{a}" for (c, p, r, a) in (predicted - expected))

    # ── False-positive control over reference-only regions ───────────────
    false_positives = 0
    for neg in truthset.negative_regions:
        neg_res = await router.call(
            "get_variants",
            {"file_path": bam, "region": neg.region, "reference": reference},
        )
        if neg_res.ok:
            false_positives += len(_parse_variants(neg_res.text))

    # ── Artifact-type recall + clean-control discrimination ──────────────
    # ── plus the safety guard: a known-artifact site must never be high ──
    # ── confidence (overconfidence on likely-spurious evidence). ─────────
    artifact_checks: list[ArtifactCheck] = []
    clean_scores: list[float] = []
    flagged_scores: list[float] = []
    overconfident_sites: list[str] = []
    confidence_mismatches: list[str] = []
    clean_not_high = True
    for site in truthset.variant_sites:
        interesting = (
            site.expected_artifact is not None
            or site.is_clean_control
            or site.expected_confidence is not None
        )
        if not interesting:
            continue
        assessment, likelihood, confidence = await _curate(router, bam, reference, site)
        risk_score = _risk_score(assessment)
        if site.expected_confidence is not None and confidence != site.expected_confidence:
            confidence_mismatches.append(
                f"{site.chrom}:{site.pos} expected {site.expected_confidence}, got {confidence}"
            )
        if site.is_clean_control:
            clean_scores.append(risk_score)
            if likelihood == "high":
                clean_not_high = False
            continue
        if site.expected_artifact is None:
            # Confidence-only positive control: no artifact assertion to make.
            continue
        present = risk_types_present(assessment)
        flagged_scores.append(risk_score)
        if likelihood == "high" and confidence == "high":
            overconfident_sites.append(f"{site.chrom}:{site.pos}")
        artifact_checks.append(
            ArtifactCheck(
                site=f"{site.chrom}:{site.pos}",
                expected_type=site.expected_artifact or "",
                found=site.expected_artifact in present,
                risk_score=risk_score,
                likelihood=likelihood,
            )
        )

    clean_discriminated = clean_not_high and (
        not clean_scores or not flagged_scores or max(clean_scores) < max(flagged_scores)
    )

    return TruthsetReport(
        dataset=truthset.dataset,
        version=truthset.version,
        detection=detection,
        false_positives=false_positives,
        artifact_checks=artifact_checks,
        clean_discriminated=clean_discriminated,
        missing_calls=missing,
        spurious_calls=spurious,
        overconfident_sites=overconfident_sites,
        confidence_mismatches=confidence_mismatches,
    )


async def _curate(
    router: ToolRouter,
    bam: str,
    reference: str,
    site: TruthSite,
) -> tuple[dict[str, Any], str, str]:
    """Call the curation tool; return (artifact_assessment, likelihood, confidence)."""
    res = await router.call(
        "get_variant_curation_summary",
        {
            "file_path": bam,
            "chrom": site.chrom,
            "pos": site.pos,
            "ref": site.ref,
            "alt": site.alt,
            "reference": reference,
            "format": "rubric",
        },
    )
    if not res.ok:
        return {}, "unknown", "unknown"
    try:
        data = json.loads(res.text)
    except (ValueError, TypeError):
        return {}, "unknown", "unknown"
    assessment = data.get("artifact_assessment") if isinstance(data, dict) else None
    likelihood = ""
    if isinstance(assessment, dict):
        likelihood = str(assessment.get("artifact_likelihood", ""))
    confidence = str(data.get("confidence", "")) if isinstance(data, dict) else ""
    return (
        assessment if isinstance(assessment, dict) else {},
        likelihood or "unknown",
        confidence or "unknown",
    )


def _risk_score(assessment: dict[str, Any]) -> float:
    score = assessment.get("risk_score")
    return float(score) if isinstance(score, (int, float)) else 0.0


# ── CLI ──────────────────────────────────────────────────────────────────


def _format_report(report: TruthsetReport, passed: bool, failures: list[str]) -> str:
    d = report.detection
    lines = [
        f"Truth set: {report.dataset} v{report.version}",
        f"  variant detection  P={d.precision:.3f} R={d.recall:.3f} "
        f"F1={d.f1:.3f}  (tp={d.tp} fp={d.fp} fn={d.fn})",
        f"  false positives    {report.false_positives} (negative regions)",
        f"  artifact recall    {report.artifact_recall:.3f}",
    ]
    for c in report.artifact_checks:
        mark = "ok" if c.found else "MISS"
        lines.append(
            f"    [{mark}] {c.site} expect={c.expected_type} "
            f"likelihood={c.likelihood} risk={c.risk_score}"
        )
    lines.append(f"  clean discriminated {report.clean_discriminated}")
    safety = report.overconfident_sites or "none"
    lines.append(f"  overconfident sites {safety} (safety guard)")
    calib = report.confidence_mismatches or "none"
    lines.append(f"  confidence mismatches {calib} (positive control)")
    lines.append("")
    lines.append("PASS" if passed else "FAIL")
    for f in failures:
        lines.append(f"  - {f}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: score a manifest and exit nonzero if floors are unmet."""
    import argparse
    import asyncio

    p = argparse.ArgumentParser(
        prog="bamcp-truthset",
        description="Score BAMCP tools against a ground-truth manifest.",
    )
    p.add_argument(
        "--manifest",
        default="tests/eval/datasets/synthetic_v1/manifest.yaml",
        help="Path to a truth-set manifest YAML.",
    )
    p.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    p.add_argument(
        "--min-recall", type=float, default=1.0, help="Minimum variant recall (default 1.0)."
    )
    p.add_argument(
        "--min-precision",
        type=float,
        default=0.95,
        help="Minimum variant precision (default 0.95).",
    )
    args = p.parse_args(argv)

    truthset = load_truthset(args.manifest)
    config = BAMCPConfig(reference=truthset.reference)
    report = asyncio.run(score_truthset(truthset, config))
    passed, failures = report.meets_floors(
        min_recall=args.min_recall, min_precision=args.min_precision
    )
    if args.json:
        print(json.dumps({**report.as_dict(), "passed": passed, "failures": failures}, indent=2))
    else:
        print(_format_report(report, passed, failures))
    return 0 if passed else 1


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())

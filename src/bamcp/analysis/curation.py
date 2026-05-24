"""Variant curation summary and clinical interpretation support.

Provides detailed curation summaries with artifact risk assessment,
quality metrics, and actionable recommendations for variant curators.
"""

from __future__ import annotations

import asyncio
import json
import logging
from typing import Any

from ..config import BAMCPConfig
from ..core.validation import validate_path, validate_variant_input
from ..middleware.telemetry import telemetry_wrapper
from .evidence import (
    build_mismatch_index,
    compute_artifact_risk,
    compute_confidence,
    compute_variant_evidence,
)

logger = logging.getLogger(__name__)

RUBRIC_VERSION = "1.0"
_VALID_FORMATS = ("text", "rubric")


def compute_near_end_fraction(position_histogram: list[int] | None) -> float:
    """Compute fraction of variant bases near read ends."""
    if not position_histogram:
        return 0.0
    total = sum(position_histogram)
    if total == 0:
        return 0.0
    # First two bins represent near-end positions
    near_end = sum(position_histogram[:2])
    return round(near_end / total, 3)


def generate_curator_recommendations(
    variant: dict,
    evidence: dict,
    artifact_risk: dict,
) -> list[str]:
    """Generate actionable recommendations for variant curators."""
    recommendations = []

    if artifact_risk["artifact_likelihood"] == "high":
        recommendations.append(
            "HIGH ARTIFACT RISK: Manual review strongly recommended. "
            "Consider orthogonal validation (e.g., Sanger sequencing)."
        )

    if evidence.get("strand_bias", 0) > 0.9:
        recommendations.append(
            "STRAND BIAS: Nearly all supporting reads on one strand. "
            "This pattern is consistent with PCR amplification artifacts."
        )

    position_hist = evidence.get("position_histogram", [])
    if position_hist:
        near_end_frac = compute_near_end_fraction(position_hist)
        if near_end_frac > 0.6:
            recommendations.append(
                "READ POSITION BIAS: Majority of variant bases are near read ends. "
                "Consider reviewing read alignments for potential mapping errors."
            )

    mapq_hist = evidence.get("mapq_histogram", [])
    if mapq_hist and sum(mapq_hist) > 0:
        low_mapq_frac = sum(mapq_hist[:3]) / sum(mapq_hist)
        if low_mapq_frac > 0.3:
            recommendations.append(
                "LOW MAPPING QUALITY: Significant fraction of reads have MAPQ < 30. "
                "The genomic region may be repetitive or have multiple mappings."
            )

    if variant.get("depth", 0) < 20:
        recommendations.append(
            "LOW COVERAGE: Depth below 20x limits confidence. "
            "Consider whether this meets clinical calling thresholds."
        )

    if not recommendations:
        recommendations.append(
            "No significant quality concerns identified. "
            "Variant appears suitable for clinical interpretation."
        )

    return recommendations


def compute_rubric_scores(variant: dict, evidence: dict, artifact_risk: dict) -> dict[str, float]:
    """Normalized 0-1 quality scores derived from existing evidence.

    Each score is computed from quantities already present in the evidence
    pipeline; no new BAM passes are required. Higher is better in all cases.
    """
    vaf = float(variant.get("vaf", 0.0))
    depth = int(variant.get("depth", 0))
    strand_bias = float(evidence.get("strand_bias", 0.0))

    mapq_hist = evidence.get("mapq_histogram") or []
    mapq_total = sum(mapq_hist)
    high_mapq = sum(mapq_hist[3:]) if len(mapq_hist) >= 4 else 0
    mapq_quality = (high_mapq / mapq_total) if mapq_total > 0 else 0.0

    pos_hist = evidence.get("position_histogram") or []
    near_end = compute_near_end_fraction(pos_hist)

    risk_score = float(artifact_risk.get("risk_score", 0.0))

    scores = {
        "vaf_quality": _clamp01(vaf * 2.5),  # 1.0 once VAF >= 0.4
        "depth_quality": _clamp01(depth / 30.0),  # 1.0 at 30x
        "strand_balance": _clamp01(1.0 - strand_bias),
        "mapq_quality": _clamp01(mapq_quality),
        "position_quality": _clamp01(1.0 - near_end),
        "artifact_risk_inverse": _clamp01(1.0 - risk_score),
    }
    return {k: round(v, 4) for k, v in scores.items()}


def _clamp01(x: float) -> float:
    """Clamp a float into the closed interval [0, 1]."""
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return x


def format_curation_summary(summary: dict) -> str:
    """Format curation summary as readable text for LLM context."""
    v = summary["variant"]
    q = summary["quality_metrics"]
    s = summary["strand_analysis"]
    p = summary["position_in_read"]
    a = summary["artifact_assessment"]

    lines = [
        f"Variant Curation Summary: {v['location']} {v['change']}",
        "=" * 50,
        "",
        "OBSERVATION",
        f"  VAF: {v['vaf']:.1%} ({v.get('alt_count', 'N/A')} alt reads / {v['depth']} total depth)",
        "",
        "QUALITY METRICS",
        f"  Base Quality: mean={q['mean_quality']:.1f}",
    ]

    if q.get("quality_distribution"):
        bins = ["0-10", "10-20", "20-30", "30-40", "40+"]
        dist_items = zip(bins, q["quality_distribution"], strict=False)
        dist = " | ".join(f"{b}:{c}" for b, c in dist_items)
        lines.append(f"  Distribution: {dist}")

    lines.extend(
        [
            "",
            "STRAND ANALYSIS",
            f"  Forward: {s['forward']} | Reverse: {s['reverse']}",
            f"  Strand Bias: {s['strand_bias']:.2f} "
            f"({'balanced' if s['is_balanced'] else 'IMBALANCED'})",
            "",
            "POSITION IN READ",
            f"  Near-end fraction: {p['near_end_fraction']:.1%}",
            "",
            f"ARTIFACT RISK: {a['artifact_likelihood'].upper()}",
            f"  Risk Score: {a['risk_score']:.2f}",
        ]
    )

    if a["risks"]:
        lines.append("  Identified Concerns:")
        for risk in a["risks"]:
            lines.append(f"    - [{risk['severity'].upper()}] {risk['description']}")

    lines.extend(
        [
            "",
            f"CONFIDENCE: {summary['confidence'].upper()}",
            "",
            "CURATOR RECOMMENDATIONS:",
        ]
    )
    for rec in summary["recommendations"]:
        lines.append(f"  * {rec}")

    return "\n".join(lines)


def _error_response(message: str, output_format: str) -> dict:
    """Wrap an error message in the same envelope shape both formats produce."""
    if output_format == "rubric":
        payload = json.dumps({"error": message, "rubric_version": RUBRIC_VERSION})
        return {"content": [{"type": "text", "text": payload}]}
    return {"content": [{"type": "text", "text": message}]}


@telemetry_wrapper("get_variant_curation_summary")
async def handle_get_variant_curation_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Generate a curation-focused summary for a specific variant.

    Returns structured data optimized for clinical interpretation, including
    artifact risk assessment, quality metrics, and recommendations for further
    review.

    ``args["format"]`` selects the output shape:
    - ``"text"`` (default): a human-readable text summary for LLM context.
    - ``"rubric"``: a JSON document with the structured summary plus a
      ``scores`` sub-dict (0-1 floats) and ``rubric_version``, suitable for
      deterministic grading.
    """
    # Import here to avoid circular dependency (tools imports from this module's siblings)
    from ..core.tools import _fetch_region_with_timeout

    file_path = args["file_path"]
    validate_path(file_path, config)
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]
    window = args.get("window", 50)
    reference = args.get("reference", config.reference)
    output_format = args.get("format", "text")
    if output_format not in _VALID_FORMATS:
        output_format = "text"

    # Validate inputs
    validation_error = validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return _error_response(validation_error, output_format)

    region = f"{chrom}:{pos - window}-{pos + window}"

    try:
        data = await _fetch_region_with_timeout(file_path, region, reference, config)
    except asyncio.TimeoutError:
        return _error_response("Timeout fetching region data", output_format)

    # Find the specific variant
    target_variant = None
    for v in data.variants:
        if v["position"] == pos and v["ref"] == ref and v["alt"] == alt:
            target_variant = v
            break

    if not target_variant:
        return _error_response(
            f"Variant {chrom}:{pos} {ref}>{alt} not found in region. "
            f"Variants found: {len(data.variants)}",
            output_format,
        )

    # Build mismatch index and compute evidence
    mismatch_index = build_mismatch_index(data.reads)
    evidence = compute_variant_evidence(mismatch_index, target_variant)
    artifact_risk = compute_artifact_risk(
        target_variant, evidence, data.reference_sequence, data.start
    )

    confidence = compute_confidence(target_variant, evidence, artifact_risk)

    # Build curation summary
    curation_summary = {
        "variant": {
            "location": f"{chrom}:{pos}",
            "change": f"{ref}>{alt}",
            "vaf": target_variant["vaf"],
            "depth": target_variant["depth"],
            "alt_count": target_variant.get("alt_count", 0),
        },
        "quality_metrics": {
            "mean_quality": evidence["mean_quality"],
            "median_quality": evidence["median_quality"],
            "quality_distribution": evidence.get("quality_histogram", []),
        },
        "strand_analysis": {
            "forward": evidence["forward_count"],
            "reverse": evidence["reverse_count"],
            "strand_bias": evidence["strand_bias"],
            "is_balanced": evidence["strand_bias"] <= 0.7,
        },
        "position_in_read": {
            "distribution": evidence.get("position_histogram", []),
            "near_end_fraction": compute_near_end_fraction(evidence.get("position_histogram")),
        },
        "artifact_assessment": artifact_risk,
        "confidence": confidence,
        "recommendations": generate_curator_recommendations(
            target_variant, evidence, artifact_risk
        ),
    }

    if output_format == "rubric":
        curation_summary["rubric_version"] = RUBRIC_VERSION
        curation_summary["scores"] = compute_rubric_scores(target_variant, evidence, artifact_risk)
        return {"content": [{"type": "text", "text": json.dumps(curation_summary)}]}

    # Format as readable text for LLM
    text_summary = format_curation_summary(curation_summary)

    return {"content": [{"type": "text", "text": text_summary}]}

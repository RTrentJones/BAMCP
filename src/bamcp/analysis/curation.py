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
from ..constants import (
    ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD,
    HIGH_CONFIDENCE_MAX_STRAND_BIAS,
    HIGH_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
)
from ..core.errors import ToolInputError
from ..core.validation import validate_path, validate_variant_input
from ..middleware.telemetry import telemetry_wrapper
from .evidence import (
    build_indel_index,
    build_mismatch_index,
    compute_artifact_risk,
    compute_confidence,
    compute_indel_evidence,
    compute_variant_evidence,
    is_indel,
    is_read_position_biased,
)

logger = logging.getLogger(__name__)

RUBRIC_VERSION = "1.0"
_VALID_FORMATS = ("text", "rubric")
# "json" is what callers reach for; it means the machine-readable rubric document.
_FORMAT_ALIASES = {"json": "rubric"}

# Intended-use boundary. BAMCP surfaces evidence to support human curation; it
# is not a diagnostic device and its output must not be treated as a clinical
# determination on its own. Carried in every curation response so the boundary
# travels with the data, not just the docs. See SAFETY.md.
INTENDED_USE = (
    "Research and visualization use only. Not a diagnostic device and not a "
    "substitute for a validated clinical pipeline or expert review. All calls "
    "require confirmation by a qualified human curator."
)


def compute_near_end_fraction(evidence: dict | None) -> float:
    """Fraction of supporting bases within the near-end cutoff of an aligned read end.

    Read from the evidence record's ``near_end_fraction`` (see
    :func:`bamcp.analysis.evidence.summarize_read_positions`). The display histogram's
    0-50bp bins are NOT used: they cover most of a short read and flag clean data.
    """
    if not evidence:
        return 0.0
    return round(float(evidence.get("near_end_fraction") or 0.0), 3)


def _read_evidence_computed(evidence: dict) -> bool:
    return bool(evidence.get("read_evidence_computed", True))


def generate_curator_recommendations(
    variant: dict,
    evidence: dict,
    artifact_risk: dict,
) -> list[str]:
    """Generate actionable recommendations for variant curators."""
    recommendations = []
    computed = _read_evidence_computed(evidence)

    if not computed:
        recommendations.append(
            "READ EVIDENCE NOT COMPUTED: strand, base-quality and read-position support "
            "for this variant could not be measured. Inspect the alignments directly; "
            "the absence of other flags here says nothing about call quality."
        )

    if artifact_risk["artifact_likelihood"] == "high":
        recommendations.append(
            "HIGH ARTIFACT RISK: Manual review strongly recommended. "
            "Consider orthogonal validation (e.g., Sanger sequencing)."
        )

    if computed and evidence.get("strand_bias", 0) > LOW_CONFIDENCE_MAX_STRAND_BIAS:
        recommendations.append(
            "STRAND BIAS: Nearly all supporting reads on one strand. "
            "This pattern is consistent with PCR amplification artifacts."
        )

    if computed and is_read_position_biased(evidence):
        recommendations.append(
            "READ POSITION BIAS: Majority of variant bases are near read ends. "
            "Consider reviewing read alignments for potential mapping errors."
        )

    mapq_hist = evidence.get("mapq_histogram", [])
    if mapq_hist and sum(mapq_hist) > 0:
        low_mapq_frac = sum(mapq_hist[:3]) / sum(mapq_hist)
        if low_mapq_frac > ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD:
            recommendations.append(
                "LOW MAPPING QUALITY: Significant fraction of reads have MAPQ < 30. "
                "The genomic region may be repetitive or have multiple mappings."
            )

    if variant.get("depth", 0) < HIGH_CONFIDENCE_MIN_DEPTH:
        recommendations.append(
            f"LOW COVERAGE: Depth below {HIGH_CONFIDENCE_MIN_DEPTH}x limits confidence. "
            "Consider whether this meets clinical calling thresholds."
        )

    if not recommendations:
        recommendations.append(
            "No automated quality flags raised. Manual curation is still "
            "required; absence of flags is not a clinical determination."
        )

    return recommendations


def compute_rubric_scores(
    variant: dict, evidence: dict, artifact_risk: dict
) -> dict[str, float | None]:
    """Normalized 0-1 quality scores derived from existing evidence.

    Each score is computed from quantities already present in the evidence
    pipeline; no new BAM passes are required. Higher is better in all cases.
    When read-level evidence was not computed (``read_evidence_computed`` False),
    the read-derived scores are ``None`` — unknown, not perfect.
    """
    vaf = float(variant.get("vaf", 0.0))
    depth = int(variant.get("depth", 0))
    strand_bias = float(evidence.get("strand_bias", 0.0))

    mapq_hist = evidence.get("mapq_histogram") or []
    mapq_total = sum(mapq_hist)
    high_mapq = sum(mapq_hist[3:]) if len(mapq_hist) >= 4 else 0
    mapq_quality = (high_mapq / mapq_total) if mapq_total > 0 else 0.0

    near_end = compute_near_end_fraction(evidence)

    risk_score = float(artifact_risk.get("risk_score", 0.0))

    scores = {
        "vaf_quality": _clamp01(vaf * 2.5),  # 1.0 once VAF >= 0.4
        "depth_quality": _clamp01(depth / 30.0),  # 1.0 at 30x
        "strand_balance": _clamp01(1.0 - strand_bias),
        "mapq_quality": _clamp01(mapq_quality),
        "position_quality": _clamp01(1.0 - near_end),
        "artifact_risk_inverse": _clamp01(1.0 - risk_score),
    }
    rounded: dict[str, float | None] = {k: round(v, 4) for k, v in scores.items()}
    if not _read_evidence_computed(evidence):
        for k in ("strand_balance", "mapq_quality", "position_quality"):
            rounded[k] = None
    return rounded


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
    ]
    lines.extend(_format_filter_check(summary.get("filter_check"), summary.get("applied_filters")))
    lines.extend(["", "QUALITY METRICS"])

    if q.get("mean_quality") is None:
        # Read evidence not computed: say so instead of rendering placeholder zeros as
        # "balanced"/clean.
        lines.extend(
            [
                "  Base Quality: NOT COMPUTED",
                "",
                "STRAND ANALYSIS",
                "  NOT COMPUTED (supporting reads not measured — not evidence of balance)",
                "",
                "POSITION IN READ",
                "  NOT COMPUTED",
            ]
        )
    else:
        lines.append(f"  Base Quality: mean={q['mean_quality']:.1f}")
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
            ]
        )

    lines.extend(
        [
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

    lines.extend(["", "INTENDED USE", f"  {summary.get('intended_use', INTENDED_USE)}"])

    return "\n".join(lines)


def _format_filter_check(check: dict | None, applied: dict | None) -> list[str]:
    """The CALLING FILTERS block: does this site pass the filters get_variants applies?"""
    if not check:
        return []
    lines = ["", "CALLING FILTERS"]
    if check["passes_filters"]:
        lines.append("  Passes the applied filters (get_variants reports it).")
    else:
        lines.append(
            "  FAILS the applied filters: "
            + ", ".join(check["failed_filters"])
            + " — get_variants does NOT report this site."
        )
    at = check["at_filters"]
    lines.append(
        f"  At the filters: {at['alt_count']}/{at['depth']} alt reads (VAF {at['vaf']:.1%})"
    )
    unf = check.get("unfiltered")
    if unf:
        lines.append(
            f"  Unfiltered (min_baseq=0): {unf['alt_count']}/{unf['depth']} alt reads "
            f"(VAF {unf['vaf']:.1%}); {check.get('low_baseq_alt_reads_removed', 0)} alt read(s) "
            "removed by min_baseq"
        )
        lines.append(f"  {check.get('note', '')}")
    if applied:
        filters = ", ".join(
            f"{k}={applied[k]}"
            for k in ("min_vaf", "min_depth", "min_alt_count", "min_mapq", "min_baseq")
            if k in applied
        )
        lines.append(f"  Applied filters: {filters}")
    return lines


def _external_status_line(name: str, block: dict | None) -> str | None:
    """Status line for an unavailable/disabled external block, else None."""
    status = block.get("status") if block else None
    if status == "unavailable":
        return f"  {name}: unavailable (lookup failed — not a confirmed absence)"
    if status == "disabled":
        return f"  {name}: disabled on this server (not queried — not a confirmed absence)"
    return None


def _format_clinical_context(context: dict) -> str:
    """Render the optional ClinVar/gnomAD context block as text."""
    lines = ["CLINICAL CONTEXT (external, research-grade)"]
    clinvar = context.get("clinvar")
    pop = context.get("population_frequency")
    status_line = _external_status_line("ClinVar", clinvar)
    if status_line:
        lines.append(status_line)
    elif clinvar:
        lines.append(
            f"  ClinVar: {clinvar.get('clinical_significance')} "
            f"({clinvar.get('stars')}-star; {clinvar.get('review_status')})"
        )
    else:
        lines.append("  ClinVar: no record found")
    status_line = _external_status_line("gnomAD", pop)
    if status_line:
        lines.append(status_line)
    elif pop:
        lines.append(
            f"  gnomAD: global_af={pop.get('global_af')} max_pop_af={pop.get('max_pop_af')}"
        )
    else:
        lines.append("  gnomAD: no record found")
    return "\n".join(lines)


def _error_response(
    message: str, output_format: str, applied_filters: dict[str, Any] | None = None
) -> dict:
    """Wrap an error message in the same envelope shape both formats produce.

    ``applied_filters`` is echoed when the error is a not-found, so a variant dropped
    by min_vaf/min_depth stays distinguishable from a true absence.
    """
    if output_format == "rubric":
        body: dict[str, Any] = {"error": message, "rubric_version": RUBRIC_VERSION}
        if applied_filters is not None:
            body["applied_filters"] = applied_filters
        return {"content": [{"type": "text", "text": json.dumps(body)}]}
    if applied_filters is not None:
        filters = ", ".join(f"{k}={v}" for k, v in applied_filters.items())
        message = f"{message}\nApplied filters: {filters}"
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
    from ..core.tools import (
        _applied_filters,
        _fetch_region_with_timeout,
        _resolve_reference_info,
        probe_site,
        site_filter_check,
    )

    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]
    window = args.get("window", 50)
    # Same resolution as every other tool: a caller-supplied reference is validated (the
    # untrusted SSRF vector), then the operator default, then the public FASTA matching the
    # BAM's build and contig style — so a bare demo/remote BAM still gets SNV evidence.
    resolved = await _resolve_reference_info(args, config)
    reference = resolved.effective
    output_format = _FORMAT_ALIASES.get(args.get("format", "text"), args.get("format", "text"))
    if output_format not in _VALID_FORMATS:
        output_format = "text"
    include_clinical = bool(args.get("include_clinical", False))

    # Validate inputs
    validation_error = validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return _error_response(validation_error, output_format)
    # Validation accepts lowercase alleles; detected variants carry uppercase.
    ref = ref.upper()
    alt = alt.upper()

    # 1-based inclusive window around the 1-based pos; clamp the start at base 1 so a
    # variant near the contig start never builds a negative/zero (invalid) region.
    region = f"{chrom}:{max(1, pos - window)}-{pos + window}"

    try:
        data = await _fetch_region_with_timeout(file_path, region, reference, config)
    except (asyncio.TimeoutError, TimeoutError):
        # _fetch_region_with_timeout re-raises the builtin TimeoutError, which is NOT
        # asyncio.TimeoutError on Python 3.10 — catch both.
        return _error_response("Timeout fetching region data", output_format)
    except (KeyError, ValueError) as e:
        # pysam raises KeyError/ValueError for a contig absent from the BAM or a bad
        # region; report it as a tool error rather than crashing the call.
        logger.warning("Curation region fetch failed: %s", e)
        return _error_response(
            f"Could not read {chrom}:{pos} from the alignment file (contig not present "
            "in the BAM, or invalid region).",
            output_format,
        )
    except OSError as e:
        logger.warning("Curation region fetch failed: %s", e)
        return _error_response("Could not read the alignment file.", output_format)

    # `pos` arrives 1-based (VCF/dbSNP convention — the same coordinate
    # get_variants / lookup_clinvar / lookup_gnomad speak), but data.variants
    # carry pysam's 0-based position, so match against pos - 1.
    target_variant = None
    for v in data.variants:
        if v["position"] == pos - 1 and v["ref"].upper() == ref and v["alt"].upper() == alt:
            target_variant = v
            break

    applied_filters = _applied_filters(
        config,
        {},
        include_max_reads=True,
        min_alt_count=config.min_alt_count,
        reference_source=resolved.source,
    )
    filter_check = site_filter_check(target_variant, config) if target_variant else None
    if target_variant is not None and filter_check and filter_check["passes_filters"]:
        # Build the support index (mismatches for SNVs, CIGAR events for indels) and
        # compute evidence.
        if is_indel(target_variant):
            evidence = compute_indel_evidence(build_indel_index(data.reads), target_variant)
        else:
            evidence = compute_variant_evidence(build_mismatch_index(data.reads), target_variant)
        artifact_risk = compute_artifact_risk(
            target_variant, evidence, data.reference_sequence, data.start
        )
    else:
        # The requested site fails the calling filters (or wasn't called at all): evaluate it
        # anyway — hiding a low-BQ one-strand artifact behind "not found" defeats curation.
        try:
            probed = await probe_site(file_path, data, pos - 1, ref, alt, reference, config)
        except (asyncio.TimeoutError, TimeoutError):
            return _error_response("Timeout fetching region data", output_format)
        except ToolInputError as e:  # e.g. REF doesn't match the reference base
            return _error_response(str(e), output_format, applied_filters=applied_filters)
        if probed is None:
            return _error_response(
                f"Variant {chrom}:{pos} {ref}>{alt} not found in region: no read supports the alt "
                "allele at this position, even with min_baseq=0 and no VAF/alt-read floor. "
                f"Variants called in the window: {len(data.variants)}.",
                output_format,
                applied_filters=applied_filters,
            )
        target_variant = probed["variant"]
        evidence = probed["evidence"]
        artifact_risk = target_variant["artifact_risk"]
        filter_check = probed["filter_check"]
    computed = bool(evidence.get("read_evidence_computed", True))

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
        "read_evidence_computed": computed,
        "quality_metrics": {
            "mean_quality": evidence["mean_quality"] if computed else None,
            "median_quality": evidence["median_quality"] if computed else None,
            "quality_distribution": evidence.get("quality_histogram", []) if computed else [],
        },
        "strand_analysis": {
            "forward": evidence["forward_count"] if computed else None,
            "reverse": evidence["reverse_count"] if computed else None,
            "strand_bias": evidence["strand_bias"] if computed else None,
            "is_balanced": (
                evidence["strand_bias"] <= HIGH_CONFIDENCE_MAX_STRAND_BIAS if computed else None
            ),
        },
        "position_in_read": {
            "distribution": evidence.get("position_histogram", []) if computed else [],
            "near_end_fraction": compute_near_end_fraction(evidence) if computed else None,
            "near_end_cutoff_bp": evidence.get("near_end_cutoff_bp"),
            "near_end_expected_fraction": evidence.get("near_end_expected_fraction"),
        },
        "filter_check": filter_check,
        "applied_filters": applied_filters,
        "artifact_assessment": artifact_risk,
        "confidence": confidence,
        "recommendations": generate_curator_recommendations(
            target_variant, evidence, artifact_risk
        ),
        "intended_use": INTENDED_USE,
    }

    # Optional bridge: attach external clinical/population context (ClinVar + gnomAD) so the
    # artifact-focused curation view can also carry clinical signal. Off by default; degrades
    # gracefully when a source has no record.
    clinical_ctx: dict[str, Any] | None = None
    if include_clinical:
        from .fusion import bam_build, clinical_context

        build, _ = await bam_build(file_path, config)
        clinical_ctx = await clinical_context(config, chrom, pos, ref, alt, build)
        curation_summary["clinical_context"] = clinical_ctx

    if output_format == "rubric":
        curation_summary["rubric_version"] = RUBRIC_VERSION
        curation_summary["scores"] = compute_rubric_scores(target_variant, evidence, artifact_risk)
        return {"content": [{"type": "text", "text": json.dumps(curation_summary)}]}

    # Format as readable text for LLM
    text_summary = format_curation_summary(curation_summary)
    if clinical_ctx is not None:
        text_summary += "\n\n" + _format_clinical_context(clinical_ctx)

    return {"content": [{"type": "text", "text": text_summary}]}

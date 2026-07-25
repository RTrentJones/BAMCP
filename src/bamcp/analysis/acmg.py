"""ACMG reasoning scaffold + the classify_variant tool.

Builds an ACMG/AMP reasoning scaffold over a fused evidence package
(:mod:`bamcp.analysis.fusion`). The scaffold shapes *how* the LLM reasons — it
enumerates the applicable criteria and hands over the evidence needed to
evaluate each — without asserting the answer. The server never classifies; it
surfaces evidence and structure, and the LLM does the reasoning. This is the
"guided reasoning" tool from ``archived/BAMCP_EVAL_HARNESS.md``, scaled to the
evidence BAMCP actually has (no SIFT/CADD/REVEL predictions yet).

Research-grade only. Not a clinical determination.
"""

from __future__ import annotations

import json
import logging
from typing import Any

from ..config import BAMCPConfig
from ..constants import (
    HIGH_CONFIDENCE_MIN_DEPTH,
    HIGH_CONFIDENCE_MIN_MEAN_QUALITY,
)
from ..middleware.telemetry import telemetry_wrapper
from .curation import INTENDED_USE
from .fusion import assemble_evidence

logger = logging.getLogger(__name__)

# Population-frequency thresholds used to select which frequency criteria apply.
# These are the standard ACMG/gnomAD anchor points, surfaced to the LLM as the
# rule to apply — not pre-applied here.
BA1_AF_THRESHOLD = 0.05  # standalone benign
BS1_AF_THRESHOLD = 0.01  # too common for a rare disorder
PM2_AF_THRESHOLD = 0.0001  # absent / ultra-rare

CLASSIFICATION_OPTIONS = [
    "Pathogenic",
    "Likely Pathogenic",
    "Uncertain Significance (VUS)",
    "Likely Benign",
    "Benign",
]


def _frequency_criteria(pop: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Population-frequency criteria (PM2/BS1/BA1), with the observed AF attached."""
    if pop is not None and pop.get("status") == "unavailable":
        # The gnomAD lookup ERRORED — absence is unknown. Emitting PM2 "absence supports
        # rarity" here would turn a network blip into pathogenic-leaning evidence, the exact
        # circular trap this guards against.
        return [
            {
                "code": "PM2_supporting",
                "description": "Absent from / extremely rare in population databases.",
                "evaluate_with": "gnomAD is UNAVAILABLE (the lookup failed), NOT confirmed "
                "absent. Do NOT apply PM2 on absence — retry the lookup before using this line.",
            }
        ]
    if not pop:
        return [
            {
                "code": "PM2_supporting",
                "description": "Absent from / extremely rare in population databases.",
                "evaluate_with": "No gnomAD entry was found. Absence supports rarity, but "
                "confirm the position was actually queried before applying.",
            }
        ]
    global_af = pop.get("global_af")
    max_af = pop.get("max_pop_af")
    return [
        {
            "code": "BA1",
            "description": f"Allele frequency > {BA1_AF_THRESHOLD:.0%} in any population "
            "(standalone benign).",
            "evaluate_with": f"max_pop_af={max_af}. Apply only if this exceeds {BA1_AF_THRESHOLD}.",
        },
        {
            "code": "BS1",
            "description": "Allele frequency greater than expected for the disorder.",
            "evaluate_with": f"global_af={global_af}, max_pop_af={max_af}. Consider if above "
            f"~{BS1_AF_THRESHOLD}.",
        },
        {
            "code": "PM2_supporting",
            "description": "Absent or at extremely low frequency in controls.",
            "evaluate_with": f"global_af={global_af}. Applies only when AF is below "
            f"~{PM2_AF_THRESHOLD}; do NOT apply if AF is common.",
        },
    ]


def _clinical_criteria(clinvar: dict[str, Any] | None) -> list[dict[str, Any]]:
    """ClinVar-derived criteria (PS1/PP5/BP6), gated on review quality."""
    if clinvar is not None and clinvar.get("status") == "unavailable":
        # Lookup ERRORED — not an authoritative "no assertion". Don't let a failed fetch
        # read as benign/absent evidence.
        return [
            {
                "code": "PS1/PP5",
                "description": "Established pathogenic assertion for this variant.",
                "evaluate_with": "ClinVar is UNAVAILABLE (the lookup failed) — this evidence "
                "line could not be retrieved. Retry before use; do not infer anything from it.",
            }
        ]
    if not clinvar:
        return [
            {
                "code": "PS1/PP5",
                "description": "Established pathogenic assertion for this variant.",
                "evaluate_with": "No ClinVar record was found for this variant (genuinely not "
                "in ClinVar). Do not invent a ClinVar classification.",
            }
        ]
    sig = clinvar.get("clinical_significance", "")
    stars = clinvar.get("stars", 0)
    return [
        {
            "code": "PP5/BP6",
            "description": "Reputable source reports the variant as pathogenic / benign.",
            "evaluate_with": f"ClinVar significance={sig!r}, review stars={stars}. Weight by "
            "review status: 0-1 star assertions are weak; >=2 stars are stronger. State the "
            "star level in your reasoning.",
        }
    ]


def _quality_assessment(observation: dict[str, Any]) -> dict[str, Any]:
    """Custom (non-ACMG) call-quality checks the LLM must weigh first."""
    return {
        "instructions": (
            "Before classifying, judge whether the variant CALL itself is trustworthy. A "
            "high artifact likelihood or low confidence should temper any pathogenic call."
        ),
        "checks": [
            {
                "check": "detected",
                "note": "Was the variant actually observed in the reads?",
                "value": observation.get("detected"),
            },
            {
                "check": "depth",
                "threshold": f">= {HIGH_CONFIDENCE_MIN_DEPTH}x for a confident call",
                "value": observation.get("depth"),
            },
            {
                "check": "mean_base_quality",
                "threshold": f">= {HIGH_CONFIDENCE_MIN_MEAN_QUALITY}",
                "value": observation.get("mean_quality"),
            },
            {
                "check": "strand_balance",
                "value": {
                    "forward": observation.get("strand_forward"),
                    "reverse": observation.get("strand_reverse"),
                },
            },
            {
                "check": "artifact_likelihood",
                "note": "high => treat the call with suspicion",
                "value": observation.get("artifact_likelihood"),
                "concerns": observation.get("artifact_concerns"),
            },
            {
                "check": "confidence",
                "value": observation.get("confidence"),
            },
        ],
    }


def build_acmg_scaffold(evidence: dict[str, Any]) -> dict[str, Any]:
    """Build an ACMG reasoning scaffold from a fused evidence package.

    Returns applicable criteria with the evidence needed to evaluate each, a
    call-quality assessment, and the required response format. It does not state
    or imply a classification.
    """
    criteria: list[dict[str, Any]] = []
    criteria.extend(_frequency_criteria(evidence.get("population_frequency")))
    criteria.extend(_clinical_criteria(evidence.get("clinvar")))
    criteria.append(
        {
            "code": "PP3/BP4",
            "description": "Computational (in-silico) evidence of impact.",
            "evaluate_with": "No SIFT/PolyPhen/CADD/REVEL predictions are available in this "
            "build. Do not fabricate them; note this evidence line as unavailable.",
        }
    )

    return {
        "instructions": (
            "Evaluate each applicable ACMG/AMP criterion below. For each, state: (1) the "
            "criterion code, (2) whether it is met, (3) the specific evidence value you used, "
            "(4) your confidence (high/medium/low). First complete the quality_assessment, "
            "then give a final classification and say what additional evidence would change it. "
            "Do not assert clinical validity."
        ),
        "classification_options": CLASSIFICATION_OPTIONS,
        "criteria_to_evaluate": criteria,
        "quality_assessment": _quality_assessment(evidence.get("observation", {})),
        "response_format": {
            "criteria_applied": [
                {
                    "code": "e.g. PM2_supporting",
                    "met": True,
                    "evidence": "...",
                    "confidence": "high",
                }
            ],
            "final_classification": f"one of: {', '.join(CLASSIFICATION_OPTIONS)}",
            "confidence": "high/medium/low",
            "key_uncertainties": ["..."],
            "additional_evidence_needed": ["..."],
        },
    }


def _format_classify_text(evidence: dict[str, Any], scaffold: dict[str, Any]) -> str:
    """Human/LLM-readable rendering of the evidence + scaffold."""
    v = evidence["variant"]
    obs = evidence.get("observation", {})
    cln = evidence.get("clinvar")
    pop = evidence.get("population_frequency")

    lines = [
        f"ACMG classification scaffold: {v['location']} {v['change']}",
        "=" * 52,
        "",
        "OBSERVED (from the BAM):",
    ]
    if obs.get("detected"):
        lines.append(
            f"  depth={obs.get('depth')} vaf={obs.get('vaf')} "
            f"confidence={obs.get('confidence')} artifact={obs.get('artifact_likelihood')}"
        )
    else:
        lines.append("  Variant NOT observed in the reads at this position.")

    lines.append("")
    lines.append("CLINVAR:")
    if cln and cln.get("status") == "unavailable":
        lines.append("  UNAVAILABLE — lookup failed (not a confirmed absence).")
    elif cln:
        lines.append(
            f"  {cln.get('clinical_significance')} ({cln.get('stars')}-star; "
            f"{cln.get('review_status')})"
        )
    else:
        lines.append("  No ClinVar record found.")

    lines.append("")
    lines.append("GNOMAD:")
    if pop and pop.get("status") == "unavailable":
        lines.append("  UNAVAILABLE — lookup failed (not confirmed absent; do not infer rarity).")
    elif pop:
        lines.append(f"  global_af={pop.get('global_af')} max_pop_af={pop.get('max_pop_af')}")
    else:
        lines.append("  No gnomAD record found.")

    lines.extend(
        [
            "",
            "SCAFFOLD (apply each; classification is YOUR call, not the server's):",
            json.dumps(scaffold, indent=2),
            "",
            "INTENDED USE",
            f"  {INTENDED_USE}",
        ]
    )
    return "\n".join(lines)


@telemetry_wrapper("classify_variant")
async def handle_classify_variant(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Assemble fused evidence + an ACMG scaffold for a variant.

    Returns text (evidence report + scaffold) plus ``structuredContent`` carrying
    the raw evidence and scaffold so an eval harness can grade the reasoning.
    Positions are 1-based. Research-grade; not a clinical determination.
    """
    from ..core.validation import validate_path, validate_variant_input

    file_path = args["file_path"]
    validate_path(file_path, config)
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]
    gene = args.get("gene")
    window = args.get("window", 50)
    reference = args.get("reference", config.reference)

    validation_error = validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return {"content": [{"type": "text", "text": json.dumps({"error": validation_error})}]}

    evidence = await assemble_evidence(
        chrom, pos, ref, alt, file_path, config, gene=gene, window=window, reference=reference
    )
    scaffold = build_acmg_scaffold(evidence)

    return {
        "content": [{"type": "text", "text": _format_classify_text(evidence, scaffold)}],
        "structuredContent": {
            "evidence": evidence,
            "scaffold": scaffold,
            "intended_use": INTENDED_USE,
        },
    }

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

import asyncio
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
# Below this allele number (called alleles) a gnomAD site is too sparsely called for its
# low AF to read as rarity: a "rare" frequency there may just be missing data.
PM2_MIN_ALLELE_NUMBER = 20_000

# ClinGen Sequence Variant Interpretation WG (Biesecker & Harrison, Genet Med 2018)
# recommends AGAINST using PP5/BP6: the assertion's underlying evidence should be
# evaluated directly instead of counting another lab's conclusion.
_PP5_BP6_SVI_NOTE = (
    "ClinGen SVI (2018) recommends NOT applying PP5/BP6; evaluate the evidence behind "
    "any assertion directly instead of counting the assertion itself."
)

_NOT_QUERIED = {
    "unavailable": "UNAVAILABLE (the lookup failed)",
    "disabled": "DISABLED on this server (the lookup was never made)",
}

CLASSIFICATION_OPTIONS = [
    "Pathogenic",
    "Likely Pathogenic",
    "Uncertain Significance (VUS)",
    "Likely Benign",
    "Benign",
]


def _reliability_reasons(pop: dict[str, Any], *, include_an: bool) -> list[str]:
    """Why a gnomAD frequency is unreliable, each reason once.

    gnomAD's own ``frequency_reliable=false`` warning usually just restates the failed QC
    filters; repeating it ("FAILED filters X and ... (fails filters X)") is noise, so a warning
    that only re-reports the filters is dropped when the filters are already listed.
    """
    filters = pop.get("filters") or []
    reasons = []
    if filters:
        reasons.append(f"the site FAILED gnomAD QC filters {list(filters)}")
    an = pop.get("an")
    if include_an and isinstance(an, int | float) and an < PM2_MIN_ALLELE_NUMBER:
        reasons.append(f"allele number is low (AN={an} < {PM2_MIN_ALLELE_NUMBER})")
    if pop.get("frequency_reliable") is False:
        warnings = [
            w for w in pop.get("warnings") or [] if not (filters and "filter" in str(w).lower())
        ]
        if warnings:
            reasons.append(f"the frequency is not reliable ({'; '.join(warnings)})")
        elif not filters:
            reasons.append("the frequency is not reliable (gnomAD flagged the frequency)")
    return reasons


def _pm2_quality_caveat(pop: dict[str, Any]) -> str | None:
    """Why this gnomAD record can't support PM2 on low AF, or None if it can.

    Tolerant of the client not surfacing ``filters``/``an``: missing fields add no caveat.
    """
    reasons = _reliability_reasons(pop, include_an=True)
    if not reasons:
        return None
    return (
        " and ".join(reasons) + ", so a low AF here is NOT evidence of rarity. Do NOT "
        "apply PM2 from this record; check coverage/another dataset first."
    )


def _benign_frequency_caveat(pop: dict[str, Any]) -> str | None:
    """Why this gnomAD record's AF can't support BA1/BS1 (a HIGH AF is unreliable too).

    Mirrors the PM2 guard: an AF gnomAD itself flags as unreliable, or one from a site that
    failed QC, must not become standalone/strong benign evidence either.
    """
    reasons = _reliability_reasons(pop, include_an=False)
    if not reasons:
        return None
    return (
        " and ".join(reasons) + ", so a high AF here is NOT evidence of benignity. Do NOT "
        "apply BA1/BS1 from this record; check another dataset first."
    )


def _gnomad_caution(pop: dict[str, Any]) -> str | None:
    """Direction-neutral caution for the GNOMAD summary line (any AF, high or low)."""
    reasons = _reliability_reasons(pop, include_an=True)
    if not reasons:
        return None
    return (
        " and ".join(reasons) + ", so treat this frequency as unreliable in either direction "
        "(neither rarity nor commonness); confirm it in another dataset before using it."
    )


def _frequency_criteria(pop: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Population-frequency criteria (PM2/BS1/BA1), with the observed AF attached."""
    status = pop.get("status") if pop is not None else None
    if status in _NOT_QUERIED:
        # The gnomAD lookup ERRORED or was never made — absence is unknown. Emitting PM2
        # "absence supports rarity" here would turn a network blip (or a config switch)
        # into pathogenic-leaning evidence, the exact circular trap this guards against.
        return [
            {
                "code": "PM2_supporting",
                "description": "Absent from / extremely rare in population databases.",
                "evaluate_with": f"gnomAD is {_NOT_QUERIED[status]}, NOT confirmed "
                "absent. Do NOT apply PM2 on absence — obtain population data before "
                "using this line.",
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
    pm2_text = (
        f"global_af={global_af}. Applies only when AF is below "
        f"~{PM2_AF_THRESHOLD}; do NOT apply if AF is common."
    )
    caveat = _pm2_quality_caveat(pop)
    if caveat:
        pm2_text = f"global_af={global_af}, but {caveat}"
    ba1_text = f"max_pop_af={max_af}. Apply only if this exceeds {BA1_AF_THRESHOLD}."
    bs1_text = f"global_af={global_af}, max_pop_af={max_af}. Consider if above ~{BS1_AF_THRESHOLD}."
    benign_caveat = _benign_frequency_caveat(pop)
    if benign_caveat:
        ba1_text = f"max_pop_af={max_af}, but {benign_caveat}"
        bs1_text = f"global_af={global_af}, max_pop_af={max_af}, but {benign_caveat}"
    return [
        {
            "code": "BA1",
            "description": f"Allele frequency > {BA1_AF_THRESHOLD:.0%} in any population "
            "(standalone benign).",
            "evaluate_with": ba1_text,
        },
        {
            "code": "BS1",
            "description": "Allele frequency greater than expected for the disorder.",
            "evaluate_with": bs1_text,
        },
        {
            "code": "PM2_supporting",
            "description": "Absent or at extremely low frequency in controls.",
            "evaluate_with": pm2_text,
        },
    ]


def _clinical_criteria(clinvar: dict[str, Any] | None) -> list[dict[str, Any]]:
    """ClinVar-derived criteria (PP5/BP6), gated on review quality.

    Labelled PP5/BP6 — "reputable source" assertions. (PS1 is a different criterion: the
    same amino-acid change as an established pathogenic variant, which BAMCP does not
    evaluate.) Every branch carries the ClinGen SVI recommendation against PP5/BP6.
    """
    status = clinvar.get("status") if clinvar is not None else None
    description = "Reputable source reports the variant as pathogenic / benign."
    if status in _NOT_QUERIED:
        # Lookup ERRORED or was never made — not an authoritative "no assertion". Don't
        # let a failed/skipped fetch read as benign/absent evidence.
        return [
            {
                "code": "PP5/BP6",
                "description": description,
                "evaluate_with": f"ClinVar is {_NOT_QUERIED[status]} — this evidence line "
                "could not be retrieved. Do not infer anything from it. " + _PP5_BP6_SVI_NOTE,
            }
        ]
    if not clinvar:
        return [
            {
                "code": "PP5/BP6",
                "description": description,
                "evaluate_with": "No ClinVar record was found for this variant (genuinely not "
                "in ClinVar). Do not invent a ClinVar classification. " + _PP5_BP6_SVI_NOTE,
            }
        ]
    # The client returns None for a record with no germline classification; never let
    # that render as a significance (or be read as an implied VUS).
    sig = clinvar.get("clinical_significance") or "no germline classification"
    stars = clinvar.get("stars", 0)
    return [
        {
            "code": "PP5/BP6",
            "description": description,
            "evaluate_with": f"ClinVar significance={sig!r}, review stars={stars}. Weight by "
            "review status: 0-1 star assertions are weak; >=2 stars are stronger. State the "
            "star level in your reasoning. " + _PP5_BP6_SVI_NOTE,
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
                "note": "Was the variant actually observed in the reads? null = the BAM "
                "read FAILED (unknown, not absent); false = no call under applied_filters.",
                "value": observation.get("detected"),
                "status": observation.get("status"),
                "reason": observation.get("reason"),
                "applied_filters": observation.get("applied_filters"),
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
                "note": "null counts = read-level evidence not computed (unknown, not balanced).",
                "value": {
                    "forward": observation.get("strand_forward"),
                    "reverse": observation.get("strand_reverse"),
                },
                "read_evidence_computed": observation.get("read_evidence_computed"),
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
    status = obs.get("status")
    filters = obs.get("applied_filters")
    filters_text = ", ".join(f"{k}={val}" for k, val in filters.items()) if filters else "unknown"
    if status == "unavailable":
        lines.append(
            f"  UNAVAILABLE — the BAM read failed ({obs.get('reason')}). Nothing is known "
            "about the reads; this is NOT evidence the variant is absent."
        )
    elif obs.get("detected"):
        lines.append(
            f"  depth={obs.get('depth')} vaf={obs.get('vaf')} "
            f"confidence={obs.get('confidence')} artifact={obs.get('artifact_likelihood')}"
        )
        if obs.get("read_evidence_computed") is False:
            lines.append(
                "  Read-level evidence (strand/quality/position) NOT COMPUTED — treat "
                "call quality as unknown."
            )
    elif status == "below_filters":
        lines.append(
            "  NOT CALLED under the applied filters (fails "
            f"{', '.join(obs.get('failed_filters') or [])}); evaluated at the site anyway:"
        )
        lines.append(
            f"  {obs.get('alt_count')}/{obs.get('depth')} alt reads (vaf={obs.get('vaf')}), "
            f"strand {obs.get('strand_forward')}F/{obs.get('strand_reverse')}R, "
            f"mean alt base quality {obs.get('mean_quality')}, "
            f"artifact={obs.get('artifact_likelihood')} "
            f"(concerns: {', '.join(obs.get('artifact_concerns') or []) or 'none'}), "
            f"confidence={obs.get('confidence')}"
        )
        lines.append(f"  Applied filters: {filters_text}")
    else:
        lines.append(
            f"  No variant call at this position under the applied filters ({filters_text})."
        )
        if obs.get("alt_reads_seen"):
            lines.append(
                f"  {obs['alt_reads_seen']} alt-supporting read(s) ARE present but below "
                "the calling thresholds — this is not a reference call."
            )

    lines.append("")
    lines.append("CLINVAR:")
    if cln and cln.get("status") in _NOT_QUERIED:
        lines.append(f"  {_NOT_QUERIED[cln['status']]} — not a confirmed absence.")
    elif cln:
        lines.append(
            f"  {cln.get('clinical_significance') or 'no germline classification'} "
            f"({cln.get('stars')}-star; "
            f"{cln.get('review_status')})"
        )
    else:
        lines.append("  No ClinVar record found.")

    lines.append("")
    lines.append("GNOMAD:")
    if pop and pop.get("status") in _NOT_QUERIED:
        lines.append(
            f"  {_NOT_QUERIED[pop['status']]} — not confirmed absent; do not infer rarity."
        )
    elif pop:
        lines.append(f"  global_af={pop.get('global_af')} max_pop_af={pop.get('max_pop_af')}")
        caveat = _gnomad_caution(pop)
        if caveat:
            lines.append(f"  CAUTION: {caveat}")
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
    from ..core.errors import ToolInputError
    from ..core.validation import validate_path, validate_variant_input

    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]
    gene = args.get("gene")
    window = args.get("window", 50)
    # Same resolution as every read tool: a caller-supplied reference is SSRF-validated (it is
    # handed to htslib), else the operator default, else the auto-selected public FASTA.
    from ..core.tools import _resolve_reference

    _, reference = await _resolve_reference(args, config)

    validation_error = validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        # Invalid input is a tool error (isError via the server envelope), not an answer.
        raise ToolInputError(validation_error)

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

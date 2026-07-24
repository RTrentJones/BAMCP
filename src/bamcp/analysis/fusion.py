"""Multi-source variant evidence fusion.

Merges the three evidence sources BAMCP already has — the observed BAM signal
(depth, VAF, strand, quality, artifact risk, confidence), ClinVar clinical
significance, and gnomAD population frequency — into one structured evidence
package. This is the piece the archived eval design
(``archived/BAMCP_EVAL_HARNESS.md``) called ``assemble_evidence``: until now the
ClinVar/gnomAD clients fed standalone lookup tools and were never combined with
the read-level evidence, so any clinical synthesis happened only in the LLM's
head. Here the server assembles the package so a scaffold can be built over it
and an eval can grade the fusion.

Research-grade only — this assembles evidence, it does not assert a
classification. See ``INTENDED_USE`` in ``curation.py``.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

from ..config import BAMCPConfig
from ..core.serialization import build_enhanced_variants

logger = logging.getLogger(__name__)


def _max_population_af(gnomad_result: Any) -> tuple[float | None, str | None]:
    """Return the highest per-population allele frequency and its population id."""
    pops = getattr(gnomad_result, "populations", None) or []
    best_af: float | None = None
    best_id: str | None = None
    for pop in pops:
        af = getattr(pop, "af", None)
        if af is None:
            continue
        if best_af is None or af > best_af:
            best_af = af
            best_id = getattr(pop, "id", None)
    return best_af, best_id


async def _lookup_clinvar_safe(config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str):
    if not getattr(config, "clinvar_enabled", True):
        return None
    from ..core.tools import get_clinvar_client

    try:
        return await get_clinvar_client(config).lookup(chrom, pos, ref, alt)
    except Exception as e:  # noqa: BLE001 — degrade to "no data" rather than fail fusion
        logger.warning("ClinVar lookup failed during fusion: %s", e)
        return None


async def _lookup_gnomad_safe(config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str):
    if not getattr(config, "gnomad_enabled", True):
        return None
    from ..core.tools import get_gnomad_client

    try:
        return await get_gnomad_client(config).lookup(chrom, pos, ref, alt)
    except Exception as e:  # noqa: BLE001
        logger.warning("gnomAD lookup failed during fusion: %s", e)
        return None


def _observation_from_variant(target: dict | None) -> dict[str, Any]:
    """Shape the BAM observation block from an evidence-enhanced variant."""
    if target is None:
        return {"detected": False}
    artifact = target.get("artifact_risk", {})
    return {
        "detected": True,
        "depth": target.get("depth"),
        "alt_count": target.get("alt_count"),
        "vaf": target.get("vaf"),
        "strand_forward": target.get("strand_forward"),
        "strand_reverse": target.get("strand_reverse"),
        "mean_quality": target.get("mean_quality"),
        "confidence": target.get("confidence"),
        "artifact_likelihood": artifact.get("artifact_likelihood"),
        "artifact_risk_score": artifact.get("risk_score"),
        "artifact_concerns": [r.get("type") for r in artifact.get("risks", [])],
    }


async def clinical_context(
    config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str
) -> dict[str, Any]:
    """Fetch just the external clinical/population context for a variant.

    Concurrent ClinVar + gnomAD lookups, shaped like the ``clinvar`` and
    ``population_frequency`` blocks of :func:`assemble_evidence`. Used by the
    curation tool's optional clinical bridge. Degrades to ``None`` blocks when a
    source has no record or errors.
    """
    clinvar_result, gnomad_result = await asyncio.gather(
        _lookup_clinvar_safe(config, chrom, pos, ref, alt),
        _lookup_gnomad_safe(config, chrom, pos, ref, alt),
    )

    clinvar_block: dict[str, Any] | None = None
    if clinvar_result is not None:
        clinvar_block = {
            "clinical_significance": clinvar_result.clinical_significance,
            "review_status": clinvar_result.review_status,
            "stars": clinvar_result.stars,
            "conditions": clinvar_result.conditions,
        }

    pop_block: dict[str, Any] | None = None
    if gnomad_result is not None:
        max_af, max_pop = _max_population_af(gnomad_result)
        pop_block = {
            "global_af": gnomad_result.global_af,
            "max_pop_af": max_af,
            "max_pop": max_pop,
            "homozygote_count": gnomad_result.homozygote_count,
        }

    return {"clinvar": clinvar_block, "population_frequency": pop_block}


async def assemble_evidence(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    file_path: str,
    config: BAMCPConfig,
    gene: str | None = None,
    window: int = 50,
    reference: str | None = None,
) -> dict[str, Any]:
    """Assemble a fused evidence package for a variant.

    Args:
        chrom, pos, ref, alt: 1-based variant coordinates (VCF/dbSNP convention).
        file_path: BAM/CRAM to read the observed evidence from.
        config: Server configuration.
        gene: Optional gene symbol for the gene-context block.
        window: Half-window (bp) fetched around ``pos`` for read-level evidence.
        reference: Reference FASTA override (defaults to ``config.reference``).

    Returns:
        A dict with ``variant``, ``observation``, ``clinvar``, and
        ``population_frequency`` (and ``gene_context`` when ``gene`` is given).
        External blocks are ``None`` when the variant is absent from that source.
    """
    from ..core.tools import _fetch_region_with_timeout

    ref_path = reference if reference is not None else config.reference
    region = f"{chrom}:{max(0, pos - window)}-{pos + window}"

    # Kick off the two network lookups concurrently with the local BAM read.
    clinvar_coro = _lookup_clinvar_safe(config, chrom, pos, ref, alt)
    gnomad_coro = _lookup_gnomad_safe(config, chrom, pos, ref, alt)

    try:
        data = await _fetch_region_with_timeout(file_path, region, ref_path, config)
    except (asyncio.TimeoutError, ValueError, OSError, KeyError) as e:
        # Contig absent from the BAM (pysam raises KeyError), unreadable file, or
        # timeout: degrade to "not observed" so classification can still proceed
        # on external evidence.
        logger.warning("BAM observation unavailable during fusion: %s", e)
        data = None

    target: dict | None = None
    if data is not None:
        enhanced, _ = build_enhanced_variants(data)
        # data.variants carry pysam's 0-based position; pos is 1-based.
        for v in enhanced:
            if v["position"] == pos - 1 and v["ref"] == ref and v["alt"] == alt:
                target = v
                break

    clinvar_result, gnomad_result = await asyncio.gather(clinvar_coro, gnomad_coro)

    clinvar_block: dict[str, Any] | None = None
    if clinvar_result is not None:
        clinvar_block = {
            "clinical_significance": clinvar_result.clinical_significance,
            "review_status": clinvar_result.review_status,
            "stars": clinvar_result.stars,
            "conditions": clinvar_result.conditions,
            "last_evaluated": clinvar_result.last_evaluated,
            "gene": clinvar_result.gene,
        }

    pop_block: dict[str, Any] | None = None
    if gnomad_result is not None:
        max_af, max_pop = _max_population_af(gnomad_result)
        pop_block = {
            "global_af": gnomad_result.global_af,
            "max_pop_af": max_af,
            "max_pop": max_pop,
            "homozygote_count": gnomad_result.homozygote_count,
            "filters": gnomad_result.filters,
            "source": gnomad_result.source,
        }

    evidence: dict[str, Any] = {
        "variant": {
            "location": f"{chrom}:{pos}",
            "change": f"{ref}>{alt}",
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
        },
        "observation": _observation_from_variant(target),
        "clinvar": clinvar_block,
        "population_frequency": pop_block,
    }
    if gene:
        evidence["gene_context"] = {
            "gene": gene,
            "clinvar_gene": getattr(clinvar_result, "gene", None),
        }

    return evidence

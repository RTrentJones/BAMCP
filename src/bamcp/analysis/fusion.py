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
from collections.abc import Callable
from typing import Any

from ..config import BAMCPConfig
from ..core.serialization import build_enhanced_variants

logger = logging.getLogger(__name__)


class _Unavailable:
    """Sentinel: an external lookup ERRORED, distinct from ``None`` (a genuine not-found).

    A failed ClinVar/gnomAD fetch must not read as "variant absent from the database" —
    that would let a network blip masquerade as rarity (PM2) or "no pathogenic assertion"
    evidence. The safe wrappers return this so the scaffold can say "unavailable" instead.
    """

    __slots__ = ()


UNAVAILABLE = _Unavailable()


def _external_block(
    result: Any, builder: Callable[[Any], dict[str, Any]], *, found_status: str = "found"
) -> dict[str, Any] | None:
    """Shape an external-evidence block with an explicit lookup ``status``.

    - ``UNAVAILABLE`` → ``{"status": "unavailable"}`` (errored — absence unknown).
    - ``None`` → ``None`` (genuine not-found; downstream criteria treat this as absence).
    - a result object → ``builder(result)`` annotated with ``status``.
    """
    if result is UNAVAILABLE:
        return {"status": "unavailable"}
    if result is None:
        return None
    block = builder(result)
    block["status"] = found_status
    return block


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
    except Exception as e:  # noqa: BLE001 — degrade to "unavailable" rather than fail fusion
        logger.warning("ClinVar lookup failed during fusion: %s", e)
        return UNAVAILABLE


async def _lookup_gnomad_safe(config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str):
    if not getattr(config, "gnomad_enabled", True):
        return None
    from ..core.tools import get_gnomad_client

    try:
        return await get_gnomad_client(config).lookup(chrom, pos, ref, alt)
    except Exception as e:  # noqa: BLE001 — degrade to "unavailable" rather than fail fusion
        logger.warning("gnomAD lookup failed during fusion: %s", e)
        return UNAVAILABLE


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
    curation tool's optional clinical bridge. A block is ``None`` for a genuine
    not-found and ``{"status": "unavailable"}`` when the source's lookup errored.
    """
    clinvar_result, gnomad_result = await asyncio.gather(
        _lookup_clinvar_safe(config, chrom, pos, ref, alt),
        _lookup_gnomad_safe(config, chrom, pos, ref, alt),
    )

    clinvar_block = _external_block(
        clinvar_result,
        lambda r: {
            "clinical_significance": r.clinical_significance,
            "review_status": r.review_status,
            "stars": r.stars,
            "conditions": r.conditions,
        },
    )

    def _pop(r: Any) -> dict[str, Any]:
        max_af, max_pop = _max_population_af(r)
        return {
            "global_af": r.global_af,
            "max_pop_af": max_af,
            "max_pop": max_pop,
            "homozygote_count": r.homozygote_count,
        }

    pop_block = _external_block(gnomad_result, _pop)

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
        External blocks are ``None`` for a genuine not-found and
        ``{"status": "unavailable"}`` when that source's lookup errored (so absence is
        never confused with a failed fetch).
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

    clinvar_block = _external_block(
        clinvar_result,
        lambda r: {
            "clinical_significance": r.clinical_significance,
            "review_status": r.review_status,
            "stars": r.stars,
            "conditions": r.conditions,
            "last_evaluated": r.last_evaluated,
            "gene": r.gene,
        },
    )

    def _pop(r: Any) -> dict[str, Any]:
        max_af, max_pop = _max_population_af(r)
        return {
            "global_af": r.global_af,
            "max_pop_af": max_af,
            "max_pop": max_pop,
            "homozygote_count": r.homozygote_count,
            "filters": r.filters,
            "source": r.source,
        }

    pop_block = _external_block(gnomad_result, _pop)

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

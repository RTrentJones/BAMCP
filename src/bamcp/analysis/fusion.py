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
from ..core.errors import ReferenceUnavailableError
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


class _Disabled:
    """Sentinel: the source is DISABLED on this server (never queried).

    Like ``UNAVAILABLE`` this is not a not-found: rendering it as ``None`` would read as
    "genuinely not in ClinVar" / "absent from gnomAD" (a PM2 nudge) for a lookup that
    never happened.
    """

    __slots__ = ()


DISABLED = _Disabled()


def _external_block(
    result: Any, builder: Callable[[Any], dict[str, Any]], *, found_status: str = "found"
) -> dict[str, Any] | None:
    """Shape an external-evidence block with an explicit lookup ``status``.

    - ``UNAVAILABLE`` → ``{"status": "unavailable"}`` (errored — absence unknown).
    - ``DISABLED`` → ``{"status": "disabled"}`` (not queried — absence unknown).
    - ``None`` → ``None`` (genuine not-found; downstream criteria treat this as absence).
    - a result object → ``builder(result)`` annotated with ``status``.
    """
    if result is UNAVAILABLE:
        return {"status": "unavailable"}
    if result is DISABLED:
        return {"status": "disabled"}
    if result is None:
        return None
    block = builder(result)
    block["status"] = found_status
    return block


def _population_block(r: Any) -> dict[str, Any]:
    """Shape a gnomAD result. ``max_pop_af`` is the client's gated population maximum
    (continental groups with adequate AN only), never a raw max over ``populations``:
    a lone allele in a tiny subgroup must not read as BA1. ``filters``/``an``/``ac`` and
    ``frequency_reliable`` let the ACMG frequency criteria refuse to read a filtered or
    low-AN site as rarity evidence. Fields are read tolerantly (``getattr``)."""
    return {
        "global_af": getattr(r, "global_af", None),
        "max_pop_af": getattr(r, "max_pop_af", None),
        "max_pop": getattr(r, "max_pop", None),
        "max_pop_an": getattr(r, "max_pop_an", None),
        "faf95_popmax": getattr(r, "faf95_popmax", None),
        "homozygote_count": getattr(r, "homozygote_count", None),
        "ac": getattr(r, "ac", None),
        "an": getattr(r, "an", None),
        "filters": list(getattr(r, "filters", None) or []),
        "frequency_reliable": getattr(r, "frequency_reliable", True),
        "warnings": list(getattr(r, "warnings", None) or []),
        "source": getattr(r, "source", None),
        "dataset": getattr(r, "dataset", None),
    }


async def _lookup_clinvar_safe(
    config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str, build: str | None = None
):
    if not getattr(config, "clinvar_enabled", True):
        return DISABLED
    from ..core.tools import get_clinvar_client

    try:
        return await get_clinvar_client(config).lookup(chrom, pos, ref, alt, build=build)
    except Exception as e:  # noqa: BLE001 — degrade to "unavailable" rather than fail fusion
        logger.warning("ClinVar lookup failed during fusion: %s", e)
        return UNAVAILABLE


async def _lookup_gnomad_safe(
    config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str, build: str | None = None
):
    if not getattr(config, "gnomad_enabled", True):
        return DISABLED
    from ..core.tools import get_gnomad_client

    try:
        return await get_gnomad_client(config).lookup(chrom, pos, ref, alt, build=build)
    except Exception as e:  # noqa: BLE001 — degrade to "unavailable" rather than fail fusion
        logger.warning("gnomAD lookup failed during fusion: %s", e)
        return UNAVAILABLE


async def bam_build(file_path: str, config: BAMCPConfig) -> tuple[str | None, str]:
    """The BAM's detected assembly for external lookups, as ``(build, source)``.

    ``build`` is None (lookups fall back to the server default) when the header can't
    be read or doesn't identify an assembly; ``source`` says which, so callers can warn.
    """
    from ..core.tools import _detect_bam_build

    try:
        info = await _detect_bam_build(file_path, config)
    except Exception as e:  # noqa: BLE001 — the observation path reports read failures
        logger.warning("Genome build detection failed during fusion: %s", e)
        return None, "server-default (BAM header unreadable)"
    detected = info.get("build")
    if detected and detected != "unknown":
        return detected, "detected-from-bam"
    return None, "server-default (build not detectable from BAM)"


def _observation_from_variant(
    target: dict | None,
    *,
    fetch_error: str | None = None,
    applied_filters: dict[str, Any] | None = None,
    alt_reads_seen: int | None = None,
    probe: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Shape the BAM observation block from an evidence-enhanced variant.

    Three distinct outcomes, never conflated:

    - ``status="unavailable"`` (``detected=None``): the BAM read itself FAILED (timeout,
      contig absent, unreadable file) — nothing is known about the reads.
    - ``status="not_called"`` (``detected=False``): the reads were read but no call
      matched under ``applied_filters``; ``alt_reads_seen`` counts alt-supporting reads
      that were present anyway (a variant below min_vaf/min_depth is NOT an absence).
    - ``status="below_filters"`` (``detected=False``): no call under ``applied_filters``, but
      the requested site was probed anyway (min_baseq=0, no VAF/alt floor — see
      :func:`bamcp.core.tools.probe_site`): its read-level evidence and artifact risk are
      reported with ``passes_filters=False`` and the ``failed_filters``.
    - ``status="observed"`` (``detected=True``): the call and its read-level evidence.
    """
    if fetch_error is not None:
        return {
            "status": "unavailable",
            "detected": None,
            "reason": fetch_error,
        }
    if target is None and probe is not None:
        return {
            "status": "below_filters",
            "detected": False,
            **_evidence_fields(probe["variant"]),
            "passes_filters": False,
            "failed_filters": probe["filter_check"]["failed_filters"],
            "filter_check": probe["filter_check"],
            "alt_reads_seen": alt_reads_seen,
            "applied_filters": applied_filters,
            "note": "Not called under the applied filters; evaluated at the requested site "
            "anyway so its evidence can be judged.",
        }
    if target is None:
        block: dict[str, Any] = {
            "status": "not_called",
            "detected": False,
            "applied_filters": applied_filters,
            "note": "No call matched at this position under the applied filters; this "
            "is not proof of a reference call (loosen min_vaf/min_depth to check).",
        }
        if alt_reads_seen is not None:
            block["alt_reads_seen"] = alt_reads_seen
        return block
    return {
        "status": "observed",
        "detected": True,
        **_evidence_fields(target),
        "passes_filters": True,
        "applied_filters": applied_filters,
    }


def _evidence_fields(target: dict[str, Any]) -> dict[str, Any]:
    """The read-level evidence fields of an evidence-enhanced (0-based) variant."""
    artifact = target.get("artifact_risk", {})
    return {
        "depth": target.get("depth"),
        "alt_count": target.get("alt_count"),
        "vaf": target.get("vaf"),
        "strand_forward": target.get("strand_forward"),
        "strand_reverse": target.get("strand_reverse"),
        "mean_quality": target.get("mean_quality"),
        "read_evidence_computed": target.get("read_evidence_computed", True),
        "confidence": target.get("confidence"),
        "artifact_likelihood": artifact.get("artifact_likelihood"),
        "artifact_risk_score": artifact.get("risk_score"),
        "artifact_concerns": [r.get("type") for r in artifact.get("risks", [])],
    }


def _count_alt_reads(data: Any, pos0: int, ref: str, alt: str) -> int:
    """Alt-supporting reads at ``pos0`` in the fetched read set, whatever the caller kept."""
    from .evidence import _indel_key, build_indel_index, build_mismatch_index

    probe = {"position": pos0, "ref": ref, "alt": alt}
    key = _indel_key(probe)
    if key is not None:
        return len(build_indel_index(data.reads).get(key, []))
    if len(ref) == 1 == len(alt):
        return len(build_mismatch_index(data.reads).get((pos0, alt), []))
    return 0


async def clinical_context(
    config: BAMCPConfig, chrom: str, pos: int, ref: str, alt: str, build: str | None = None
) -> dict[str, Any]:
    """Fetch just the external clinical/population context for a variant.

    Concurrent ClinVar + gnomAD lookups, shaped like the ``clinvar`` and
    ``population_frequency`` blocks of :func:`assemble_evidence`. Used by the
    curation tool's optional clinical bridge. A block is ``None`` for a genuine
    not-found, ``{"status": "unavailable"}`` when the source's lookup errored, and
    ``{"status": "disabled"}`` when the source is switched off on this server. ``build``
    is the assembly the coordinates are in (the BAM's); None means the server default.
    """
    ref = ref.upper()
    alt = alt.upper()
    clinvar_result, gnomad_result = await asyncio.gather(
        _lookup_clinvar_safe(config, chrom, pos, ref, alt, build),
        _lookup_gnomad_safe(config, chrom, pos, ref, alt, build),
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

    pop_block = _external_block(gnomad_result, _population_block)

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
        External blocks are ``None`` for a genuine not-found,
        ``{"status": "unavailable"}`` when that source's lookup errored, and
        ``{"status": "disabled"}`` when the source is switched off on this server (so
        absence is never confused with a failed or skipped fetch). ``observation`` has
        a ``status`` of ``observed``/``not_called``/``unavailable`` (see
        :func:`_observation_from_variant`).
    """
    from ..core.tools import (
        _applied_filters,
        _fetch_region_with_timeout,
        probe_site,
        site_filter_check,
    )

    # Validation accepts lowercase alleles; calls and the external sources use uppercase.
    ref = ref.upper()
    alt = alt.upper()
    ref_path = reference if reference is not None else config.reference
    # 1-based inclusive window around the 1-based pos (start clamped at base 1).
    region = f"{chrom}:{max(1, pos - window)}-{pos + window}"

    async def _observe() -> tuple[Any, str | None]:
        try:
            return await _fetch_region_with_timeout(file_path, region, ref_path, config), None
        except (asyncio.TimeoutError, TimeoutError):
            # Builtin TimeoutError (what _fetch_region_with_timeout raises) differs from
            # asyncio.TimeoutError on Python 3.10; catch both.
            logger.warning("BAM observation timed out during fusion")
            return None, "timeout reading the alignment file"
        except (KeyError, ValueError) as e:
            logger.warning("BAM observation unavailable during fusion: %s", e)
            return None, "contig not present in the alignment file, or invalid region"
        except OSError as e:
            logger.warning("BAM observation unavailable during fusion: %s", e)
            return None, "alignment file could not be read"

    # External sources must be queried in the BAM's assembly: a GRCh37 BAM looked up at
    # GRCh38 coordinates reads as a false "not in ClinVar" / wrong gnomAD variant.
    build, build_source = await bam_build(file_path, config)

    # Run the two network lookups concurrently WITH the local BAM read: schedule all
    # three as tasks up front so they genuinely overlap, and cancel any still pending if
    # one raises unexpectedly (gather alone would leave the siblings running).
    tasks = [
        asyncio.ensure_future(_observe()),
        asyncio.ensure_future(_lookup_clinvar_safe(config, chrom, pos, ref, alt, build)),
        asyncio.ensure_future(_lookup_gnomad_safe(config, chrom, pos, ref, alt, build)),
    ]
    try:
        observed, clinvar_result, gnomad_result = await asyncio.gather(*tasks)
    except BaseException:
        for t in tasks:
            t.cancel()
        raise
    data, fetch_error = observed

    target: dict | None = None
    alt_reads_seen: int | None = None
    probe: dict[str, Any] | None = None
    if data is not None:
        enhanced, _ = build_enhanced_variants(data)
        # data.variants carry pysam's 0-based position; pos is 1-based.
        for v in enhanced:
            if (
                v["position"] == pos - 1
                and str(v["ref"]).upper() == ref
                and str(v["alt"]).upper() == alt
            ):
                target = v
                break
        if target is not None and not site_filter_check(target, config)["passes_filters"]:
            target = None  # called by the parser but under the alt-read floor
        if target is None:
            # Evaluate the requested site anyway (min_baseq=0, no VAF/alt floor).
            try:
                probe = await probe_site(file_path, data, pos - 1, ref, alt, ref_path, config)
            except (asyncio.TimeoutError, TimeoutError):
                logger.warning("Site probe timed out during fusion")
            except ReferenceUnavailableError:
                # No reference base to verify REF against: don't probe (it could count ordinary
                # reads as ALT); the site is reported as not called from the calls we have.
                probe = None
            alt_reads_seen = (
                probe["filter_check"].get("unfiltered", probe["filter_check"]["at_filters"])[
                    "alt_count"
                ]
                if probe
                else _count_alt_reads(data, pos - 1, ref, alt)
            )

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

    pop_block = _external_block(gnomad_result, _population_block)

    evidence: dict[str, Any] = {
        "variant": {
            "build": build,
            "build_source": build_source,
            "location": f"{chrom}:{pos}",
            "change": f"{ref}>{alt}",
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
        },
        "observation": _observation_from_variant(
            target,
            fetch_error=fetch_error,
            applied_filters=_applied_filters(
                config, {}, include_max_reads=True, min_alt_count=config.min_alt_count
            ),
            alt_reads_seen=alt_reads_seen,
            probe=probe,
        ),
        "clinvar": clinvar_block,
        "population_frequency": pop_block,
    }
    if gene:
        evidence["gene_context"] = {
            "gene": gene,
            "clinvar_gene": getattr(clinvar_result, "gene", None),  # None for sentinels
        }

    return evidence

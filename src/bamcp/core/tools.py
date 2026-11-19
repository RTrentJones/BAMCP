"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import json
import logging
import re
from collections import OrderedDict
from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import Any

import pysam

from ..analysis.evidence import (
    _not_computed_evidence,
    _vcf_evidence,
    annotate_interpretation,
    compute_artifact_risk,
    compute_confidence,
    enhance_vcf_variants,
)
from ..clients.genes import GeneLookupError
from ..config import BAMCPConfig
from ..constants import (
    BAM_PARSE_TIMEOUT_SECONDS,
    DEFAULT_CACHE_TTL_SECONDS,
    GET_VARIANTS_MAX_RESULTS,
    INDEL_DETECTION_MAX_REGION,
    REGION_TILE_SIZE,
    VCF_MAX_RECORDS,
    VIEWER_RESOURCE_URI,
)
from ..middleware.telemetry import telemetry_wrapper
from . import remote as _remote
from .errors import ReferenceUnavailableError, ToolInputError
from .filters import failed_candidate_filters, passes_candidate_filters

# Re-exports: the handlers moved to sibling modules stay importable (and patchable for
# code that still looks them up here) as ``bamcp.core.tools.<name>``. A test that patches a
# name must patch the module that looks it up at call time — see each module's docstring.
from .lookups import (  # noqa: F401
    _CLINVAR_DISCLAIMER,
    _GNOMAD_DISCLAIMER,
    _LOOKUP_DISABLED,
    _LOOKUP_FOUND,
    _LOOKUP_NOT_FOUND,
    _LOOKUP_REASONS,
    _LOOKUP_UNAVAILABLE,
    _lookup_build,
    _lookup_disabled,
    _lookup_failure,
    _lookup_not_found,
    _lookup_result,
    _validate_lookup_input,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
)
from .parsers import (
    EXCLUDED_READ_FLAGS,
    MAX_REGION_SIZE,
    OVERLAPPING_MATES_POLICY,
    RegionData,
    annotate_vcf_snv_support,
    fetch_candidate_variants_only,
    fetch_coverage_only,
    fetch_reference_sequence,
    fetch_region,
    format_region,
    load_vcf_variants,
    open_fasta,
    parse_region,
    vcf_truncated,
)
from .remote import (  # noqa: F401
    _clear_warm_inputs,
    _download_index_streaming,
    _ensure_cached_index,
    _ensure_reference_indexes,
    _index_download_lock,
    _index_http_client,
    _is_remote,
    _mark_warm,
    _open_remote_inputs_sync,
    _sibling_url,
    _stream_index_to_file,
    _try_download_index,
    _warm_remote_inputs,
)
from .scan import (  # noqa: F401
    _scan_chunk,
    _scan_contig_progressive,
    _ScanProgress,
    handle_scan_variants,
)
from .serialization import build_enhanced_variants, serialize_region_data
from .services import (  # noqa: F401
    _PARSE_EXECUTOR,
    _SERVICES_REGISTRY_MAX,
    BAMCPServices,
    _approx_region_bytes,
    _cache_key,
    _get_cached_region,
    _get_parse_gate,
    _local_fingerprint,
    _parse_gates,
    _pending_service_closes,
    _RegionCacheEntry,
    _run_in_parse_pool,
    _schedule_services_close,
    _services_registry,
    _set_cached_region,
    close_external_clients,
    get_cache,
    get_clinvar_client,
    get_gene_client,
    get_gnomad_client,
    get_services,
)
from .validation import (
    is_cram,
    validate_path,
    validate_reference_path,
    validate_region,
    validate_variant_file_path,
)

logger = logging.getLogger(__name__)
_CANDIDATE_VARIANT_DISCLAIMER = (
    "Note: BAMCP reports research-grade candidate variants from read-level evidence; "
    "confirm with validated variant-calling and clinical workflows before use."
)

_INTERNAL_COORDINATE_SYSTEM = "0-based-half-open"
_LLM_COORDINATE_SYSTEM = "1-based-inclusive"

# Where get_variants sources its variants from (see handle_get_variants).
_VARIANT_SOURCES = ("auto", "vcf", "bamcp")

# Self-describing output: every variant-producing response echoes the thresholds it
# actually applied, so a variant's *absence* is interpretable (filtered out vs. truly
# not present) instead of silent. This closes the same "silent behavior misleads the
# caller" gap as the genome-build fix — an aggressive min_vaf/min_depth can drop a real
# variant and, without this, the empty result reads as a confident reference call.
_FILTER_ABSENCE_NOTE = (
    "A candidate variant's absence at a site may reflect these applied filters rather "
    "than a true reference call; loosen them (e.g. min_vaf/min_depth) to see more."
)


def _detection_status(
    region: str, reference: str | None, vcf_primary: bool = False
) -> dict[str, Any]:
    """Which candidate detectors actually ran for this request — so an empty result
    with detection OFF is never mistaken for a confident "no variants here".

    SNV/indel candidates need a reference to define the reference allele; with
    ``variant_source="vcf"`` the VCF is authoritative and BAMCP adds none; indels need
    a pileup pass that is skipped above ``INDEL_DETECTION_MAX_REGION``.
    """
    _, start, end = parse_region(region)
    if vcf_primary:
        reason: str | None = "variant_source=vcf (the VCF is authoritative)"
    elif not reference:
        reason = "no reference"
    else:
        reason = None
    status: dict[str, Any] = {
        "snv_detection": "enabled" if reason is None else f"disabled: {reason}"
    }
    if reason is not None:
        status["indel_detection"] = False
        status["indel_detection_note"] = f"disabled: {reason}"
    elif end - start > INDEL_DETECTION_MAX_REGION:
        status["indel_detection"] = False
        status["indel_detection_note"] = (
            f"skipped: region {end - start:,}bp exceeds {INDEL_DETECTION_MAX_REGION:,}bp; "
            "request a smaller region to detect indels"
        )
    else:
        status["indel_detection"] = True
    return status


def _applied_filters(
    config: BAMCPConfig,
    args: dict[str, Any],
    *,
    include_max_reads: bool,
    detection: dict[str, Any] | None = None,
    reads_truncated: bool | None = None,
    min_alt_count: int | None = None,
    reference_source: str | None = None,
) -> dict[str, Any]:
    """The read/variant thresholds actually applied to this response.

    ``min_alt_count`` and ``reference_source`` are echoed only when the caller passes them,
    i.e. only by tools that actually apply the alt-read floor / resolved the reference
    (``explicit`` | ``config`` | ``auto`` | ``none``).

    ``min_vaf``/``min_depth`` honor per-call overrides; ``min_mapq``/``min_baseq`` are
    operator config. ``max_reads`` gates how many reads are *displayed* (viewer tools),
    so it's included only where reads are returned — not for readless variant calls.
    The fixed read exclusions and overlapping-mates policy are always echoed;
    ``detection`` (see :func:`_detection_status`) and ``reads_truncated`` are merged in
    when the caller knows them.
    """
    filters: dict[str, Any] = {
        "min_vaf": args.get("min_vaf", config.min_vaf),
        "min_depth": args.get("min_depth", config.min_depth),
        "min_mapq": config.min_mapq,
        "min_baseq": config.min_baseq,
        "excluded_reads": list(EXCLUDED_READ_FLAGS),
        "overlapping_mates": OVERLAPPING_MATES_POLICY,
    }
    if min_alt_count is not None:
        filters["min_alt_count"] = min_alt_count
    if reference_source is not None:
        filters["reference_source"] = reference_source
    if include_max_reads:
        filters["max_reads"] = config.max_reads
    if reads_truncated is not None:
        filters["reads_truncated"] = reads_truncated
    if detection:
        filters.update(detection)
    return filters


def _vcf_truncation_filters(data: RegionData) -> dict[str, Any]:
    """applied_filters entries for a VCF overlay that stopped at VCF_MAX_RECORDS (else empty).

    Past the cap the VCF was not read, so a VCF site's absence there is not a call.
    """
    if not data.vcf_truncated:
        return {}
    return {"vcf_truncated": True, "vcf_max_records": VCF_MAX_RECORDS}


def _format_filters_line(filters: dict[str, Any]) -> str:
    """Compact human-readable rendering of :func:`_applied_filters` for text summaries."""
    parts = [
        f"VAF≥{filters['min_vaf']}",
        f"depth≥{filters['min_depth']}",
        f"MAPQ≥{filters['min_mapq']}",
        f"baseQ≥{filters['min_baseq']}",
    ]
    if "min_alt_count" in filters:
        parts.append(f"alt reads≥{filters['min_alt_count']}")
    if "reference_source" in filters:
        parts.append(f"reference: {filters['reference_source']}")
    if "max_reads" in filters:
        parts.append(f"≤{filters['max_reads']} reads shown")
    if filters.get("reads_truncated"):
        parts.append(
            "read list TRUNCATED at max_reads (SNV evidence from pileup; indel evidence incomplete)"
        )
    if filters.get("vcf_truncated"):
        parts.append(
            f"VCF TRUNCATED at {filters.get('vcf_max_records')} records "
            "(later VCF sites not read; narrow the region)"
        )
    if "excluded_reads" in filters:
        parts.append(f"excluded reads: {'/'.join(filters['excluded_reads'])}")
    snv = filters.get("snv_detection")
    if snv and snv != "enabled":
        parts.append(f"SNV detection {snv}")
    if filters.get("indel_detection") is False:
        parts.append(f"indel detection {filters.get('indel_detection_note', 'off')}")
    return ", ".join(parts)


def _tile_bounds(start: int, end: int, tile: int) -> tuple[int, int]:
    """Snap [start, end) outward to fixed-width tile boundaries."""
    t_start = (start // tile) * tile
    t_end = ((end - 1) // tile + 1) * tile
    return t_start, t_end


def _slice_readless_region_data(data: RegionData, start: int, end: int) -> RegionData:
    """Slice a readless (coverage/variants) tile RegionData down to [start, end).

    Coverage and reference are per-position arrays sliced by offset; variants are
    position-anchored so they filter by position. Reads are always empty in the
    readless modes this serves.
    """
    # The tile's end is clamped to the contig length by the parser; never slice past it.
    end = min(end, data.end)
    offset = start - data.start
    length = max(end - start, 0)
    coverage = data.coverage[offset : offset + length] if data.coverage else []
    reference = (
        data.reference_sequence[offset : offset + length]
        if data.reference_sequence is not None
        else None
    )
    variants = [v for v in data.variants if start <= v["position"] < end]
    return RegionData(
        contig=data.contig,
        start=start,
        end=end,
        reads=[],
        coverage=coverage,
        variants=variants,
        reference_sequence=reference,
    )


async def _compute_region_data(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    mode: str,
    min_vaf: float | None,
    min_depth: int | None,
    index_path: str | None,
    detect: bool = True,
) -> RegionData:
    """Run the mode-appropriate parser for an exact region under the parse timeout.

    ``detect=False`` skips candidate detection in the full-read path (used by the
    VCF-primary fetch, whose local candidates would be discarded anyway).
    """
    effective_min_vaf = min_vaf if min_vaf is not None else config.min_vaf
    effective_min_depth = min_depth if min_depth is not None else config.min_depth

    if mode == "coverage":
        return await asyncio.wait_for(
            _run_in_parse_pool(
                fetch_coverage_only,
                file_path,
                region,
                reference,
                config.min_mapq,
                config.min_baseq,
                index_path,
            ),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
    if mode == "variants" and reference:
        return await asyncio.wait_for(
            _run_in_parse_pool(
                fetch_candidate_variants_only,
                file_path,
                region,
                reference,
                config.min_mapq,
                config.min_baseq,
                index_path,
                effective_min_vaf,
                effective_min_depth,
            ),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
    return await asyncio.wait_for(
        _run_in_parse_pool(
            fetch_region,
            file_path,
            region,
            reference,
            config.max_reads,
            config.min_mapq,
            config.min_baseq,
            index_path,
            effective_min_vaf,
            effective_min_depth,
            detect,
        ),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )


async def _fetch_readless_tiled(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    mode: str,
    min_vaf: float | None,
    min_depth: int | None,
    region_cache: OrderedDict[tuple, _RegionCacheEntry],
    ttl: int,
) -> RegionData:
    """Serve a readless (coverage/variants) query from a per-tile cache.

    The heavy parse runs once per fixed-width tile and is sliced down to the
    requested region, so panning the viewport within a tile reuses cached work.
    Reads can't tile this way — max_reads truncation over a wider tile could drop
    reads from the requested sub-region — so this path is coverage/variants only.
    """
    contig, r_start, r_end = parse_region(region)
    t_start, t_end = _tile_bounds(r_start, r_end, REGION_TILE_SIZE)
    tile_size = t_end - t_start
    # Fall back to the exact region when widening to a tile would either exceed the
    # parser's max region, or (variants mode) cross the indel-detection threshold
    # the exact request stays under — otherwise the widened tile would make
    # _detect_indels_if_small skip the pileup and silently drop indel candidates.
    if tile_size > MAX_REGION_SIZE or (
        mode == "variants" and tile_size > INDEL_DETECTION_MAX_REGION >= (r_end - r_start)
    ):
        t_start, t_end = r_start, r_end
    # Internal 0-based numbers → the 1-based inclusive region string parse_region reads.
    tile_region = format_region(contig, t_start, t_end)

    tile_key = _cache_key(file_path, tile_region, reference, config, mode, min_vaf, min_depth, None)
    tile_data = _get_cached_region(region_cache, tile_key, ttl)
    if tile_data is None:
        index_path = await _ensure_cached_index(file_path, config)
        tile_data = await _compute_region_data(
            file_path, tile_region, reference, config, mode, min_vaf, min_depth, index_path
        )
        _set_cached_region(region_cache, tile_key, tile_data, ttl)

    if t_start == r_start and t_end == r_end:
        return tile_data
    return _slice_readless_region_data(tile_data, r_start, r_end)


def _has_vcf_snv(variants: list[dict]) -> bool:
    """Whether any variant is a VCF single-nucleotide site (needs read-support scoring)."""
    return any(v.get("source") == "vcf" and len(v["ref"]) == 1 == len(v["alt"]) for v in variants)


async def _annotate_vcf_support(
    file_path: str,
    region: str,
    reference: str | None,
    variants: list[dict],
    config: BAMCPConfig,
    index_path: str | None,
) -> None:
    """Attach truncation-free read-level support to VCF SNVs (off the event loop)."""
    await asyncio.wait_for(
        _run_in_parse_pool(
            annotate_vcf_snv_support,
            file_path,
            region,
            reference,
            variants,
            config.min_mapq,
            config.min_baseq,
            index_path,
        ),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )


async def _fetch_region_with_timeout(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    min_vaf: float | None = None,
    min_depth: int | None = None,
    mode: str = "full",
    vcf_path: str | None = None,
    vcf_primary: bool = False,
) -> RegionData:
    """Fetch region data, converting a parse timeout into a diagnosable error.

    A first call to a not-yet-seen remote file/reference first opens it under the separate
    ``COLD_START_TIMEOUT_SECONDS`` budget (see :func:`_warm_remote_inputs`).

    Every ``asyncio.wait_for`` under this call raises a bare ``asyncio.TimeoutError``,
    whose ``str()`` is the EMPTY STRING. FastMCP renders a failed tool call as
    ``f"Error executing tool {name}: {exc}"``, so a timeout reached the caller as
    literally ``"Error executing tool visualize_region:"`` — no exception type, no
    region, no hint. That is unactionable, and it is what made this failure so hard
    to place: the tools that open a remote reference (visualize_region,
    get_region_summary) time out on the first, cold call while get_coverage — whose
    readless path never opens the FASTA — kept working, which reads like an
    unrelated per-tool bug. Re-raise with the context needed to tell "too slow" from
    "broken", and note that a retry is usually warm enough to succeed.
    """
    # Every tool's alignment decode funnels through here (scan_variants requires a reference
    # up front), so this is the one CRAM header-reference gate.
    _require_cram_reference(file_path, reference, config)
    # Cold start: fetch the index (never under a parse budget), then give a not-yet-seen remote
    # BAM/reference its one-off open under the separate cold budget, so the parse that follows
    # runs warm under the normal BAM_PARSE_TIMEOUT_SECONDS.
    if _is_remote(file_path) or _is_remote(reference):
        index_path = await _ensure_cached_index(file_path, config)
        await _warm_remote_inputs(file_path, reference, index_path, config)
    await _check_cram_header_references(file_path, reference, config)
    try:
        return await _fetch_region_data(
            file_path,
            region,
            reference,
            config,
            min_vaf=min_vaf,
            min_depth=min_depth,
            mode=mode,
            vcf_path=vcf_path,
            vcf_primary=vcf_primary,
        )
    except asyncio.TimeoutError as e:
        raise TimeoutError(
            f"Timed out after {BAM_PARSE_TIMEOUT_SECONDS:g}s reading {region} "
            f"(mode={mode}){_timeout_hint(region, reference)}."
        ) from e


# Above this span a parse timeout is about the region's size, not a cold input.
TIMEOUT_LARGE_REGION_BP = 2_000


def _timeout_hint(region: str, reference: str | None) -> str:
    """What to do about a parse timeout: shrink a large window, or retry once a cold input warms.

    A gene-scale window (tens of kb of reads) is slow however warm the inputs are, so the
    cold-reference hint would send the caller into a pointless retry; it is given only when the
    remote reference really has not been opened yet.
    """
    try:
        _, start, end = parse_region(region)
    except ValueError:
        start = end = 0
    if end - start > TIMEOUT_LARGE_REGION_BP:
        return (
            f"; the region spans {end - start:,} bp — request a smaller window (≤200 bp for "
            "visualize_region/get_variants), or use get_coverage for a gene-scale region"
        )
    if _is_remote(reference) and reference not in _remote._warm_inputs:
        return (
            "; the remote reference FASTA's index is fetched on first use and is often "
            "the slow part — retrying usually succeeds once it is warm"
        )
    return ""


async def _fetch_region_data(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    min_vaf: float | None = None,
    min_depth: int | None = None,
    mode: str = "full",
    vcf_path: str | None = None,
    vcf_primary: bool = False,
) -> RegionData:
    """Fetch region data from BAM/CRAM file with timeout protection.

    Args:
        file_path: Path to BAM/CRAM file.
        region: Genomic region string.
        reference: Path to reference FASTA.
        config: Server configuration.
        min_vaf: Minimum VAF threshold (uses config default if None).
        min_depth: Minimum depth threshold (uses config default if None).

    Returns:
        RegionData with reads, coverage, and variants.

    Raises:
        asyncio.TimeoutError: If BAM parsing exceeds timeout.
    """
    region_cache = get_services(config).region_cache
    ttl = min(config.cache_ttl, DEFAULT_CACHE_TTL_SECONDS)

    if vcf_path is not None:
        await asyncio.to_thread(validate_variant_file_path, vcf_path, config)

    # Readless coverage/variant queries (no reads, no VCF overlay) are tile-cached
    # so overlapping viewports reuse work. variants mode needs a reference for the
    # candidate-detection fast path; without one it falls through to the read parser.
    if mode in ("coverage", "variants") and vcf_path is None and (mode != "variants" or reference):
        return await _fetch_readless_tiled(
            file_path, region, reference, config, mode, min_vaf, min_depth, region_cache, ttl
        )

    key = _cache_key(
        file_path, region, reference, config, mode, min_vaf, min_depth, vcf_path, vcf_primary
    )
    if cached := _get_cached_region(region_cache, key, ttl):
        return cached

    # Readless VCF-only shortcut: skip local candidate detection when the VCF is
    # authoritative (vcf_primary) or when there's no reference to detect against.
    if mode == "variants" and vcf_path and (vcf_primary or not reference):
        contig, start, end = parse_region(region)
        vcf_variants = await asyncio.wait_for(
            _run_in_parse_pool(load_vcf_variants, vcf_path, region),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
        # Only touch the BAM when there's an SNV whose support we can actually count.
        if _has_vcf_snv(vcf_variants):
            index_path = await _ensure_cached_index(file_path, config)
            await _annotate_vcf_support(
                file_path, region, reference, vcf_variants, config, index_path
            )
        # Keep the reference slice so downstream artifact scoring (homopolymer
        # context) matches the full/viewer path for VCF-primary SNVs.
        ref_seq = None
        if reference:
            ref_seq = await asyncio.wait_for(
                _run_in_parse_pool(fetch_reference_sequence, reference, region),
                timeout=BAM_PARSE_TIMEOUT_SECONDS,
            )
        return _set_cached_region(
            region_cache,
            key,
            RegionData(
                contig=contig,
                start=start,
                end=end,
                reads=[],
                coverage=[],
                variants=vcf_variants,
                reference_sequence=ref_seq,
                vcf_truncated=vcf_truncated(vcf_variants),
            ),
            ttl,
        )

    # Download and cache index for remote files (if not already cached)
    index_path = await _ensure_cached_index(file_path, config)
    # VCF-primary discards local candidates, so skip detecting them (full-read path).
    data = await _compute_region_data(
        file_path,
        region,
        reference,
        config,
        mode,
        min_vaf,
        min_depth,
        index_path,
        detect=not vcf_primary,
    )

    if vcf_path:
        vcf_variants = await asyncio.wait_for(
            _run_in_parse_pool(load_vcf_variants, vcf_path, region),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
        data.vcf_truncated = vcf_truncated(vcf_variants)
        if vcf_primary:
            # VCF is the authoritative variant set — drop BAMCP's local candidates,
            # but keep the reads/coverage so BAMCP can show evidence at each site.
            data.variants = vcf_variants
        else:
            # auto: overlay VCF on the local candidates, de-duplicating exact matches.
            seen = {(v["position"], v["ref"], v["alt"]) for v in data.variants}
            data.variants.extend(
                v for v in vcf_variants if (v["position"], v["ref"], v["alt"]) not in seen
            )
        await _annotate_vcf_support(file_path, region, reference, data.variants, config, index_path)

    return _set_cached_region(region_cache, key, data, ttl)


# -- Tool Handlers -----------------------------------------------------------


def _one_based(variant: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of a detector variant with a 1-based ``position``.

    Internally BAMCP works in pysam's 0-based coordinates (``detect_variants``
    reports ``start + i``), which keeps variant positions aligned with reads,
    coverage, and the reference array for the viewer and evidence math. But
    every genomics consumer a caller reaches next — VCF, dbSNP, ClinVar, gnomAD,
    IGV — is 1-based, so a 0-based position handed to ``lookup_clinvar`` /
    ``lookup_gnomad`` misses by one. Convert at the LLM-facing boundary so the
    positions callers see (and pass straight into the lookup tools) are 1-based.
    """
    return {**variant, "position": variant["position"] + 1}


def _format_genotypes(samples: dict[str, Any]) -> list[str]:
    """Render per-sample genotypes as ``name=allele/allele`` (missing alleles as ``.``)."""
    formatted: list[str] = []
    for name, sample_data in samples.items():
        gt = sample_data.get("GT") if isinstance(sample_data, dict) else None
        if gt is None:
            gt_str = "."
        else:
            alleles = gt if isinstance(gt, (list, tuple)) else (gt,)
            gt_str = "/".join("." if a is None else str(a) for a in alleles)
        formatted.append(f"{name}={gt_str}")
    return formatted


def _format_variant_line(variant: dict[str, Any]) -> str:
    """One human-readable summary line for a candidate variant.

    ``variant`` uses BAMCP's internal 0-based position; the rendered line is 1-based
    (VCF/dbSNP convention, see :func:`_one_based`). Structural variants are described
    by type and span (their ``ref>alt`` is symbolic); every kind appends per-sample
    genotypes when a multi-sample VCF supplied them.
    """
    pos = variant["position"] + 1
    contig = variant["contig"]

    if variant.get("variant_kind") == "sv":
        alt = str(variant.get("alt", ""))
        sv_type = variant.get("sv_type") or (alt.strip("<>") if alt.startswith("<") else "SV")
        line = f"  {contig}:{pos} <{sv_type}>"
        # sv_end (END, Number=1) and sv_len (SVLEN, normalized per-ALT by load_vcf_variants)
        # arrive as scalars.
        sv_end = variant.get("sv_end")
        if sv_end:
            line += f" span={pos}-{sv_end}"
        sv_len = variant.get("sv_len")
        if sv_len is not None:
            line += f" len={sv_len}"
    else:
        line = (
            f"  {contig}:{pos} {variant['ref']}>{variant['alt']} "
            f"VAF={variant['vaf']:.1%} depth={variant['depth']}"
        )

    genotypes = _format_genotypes(variant.get("samples") or {})
    if genotypes:
        line += f" [{' '.join(genotypes)}]"
    return line


# Reference provenance vocabulary, echoed as applied_filters.reference_source.
REFERENCE_SOURCE_EXPLICIT = "explicit"  # the caller passed `reference`
REFERENCE_SOURCE_CONFIG = "config"  # the operator's BAMCP_REFERENCE
REFERENCE_SOURCE_AUTO = "auto"  # registry FASTA matching the BAM's build + contig style
REFERENCE_SOURCE_NONE = "none"  # no reference: SNV/indel detection is off

# Header-cache bound (entries). A header with thousands of decoy/HLA contigs is ~100s of KB
# as Python objects, so keep this small; a miss costs one header read.
_HEADER_CACHE_MAX = 32


@dataclass(frozen=True)
class ResolvedReference:
    """The reference a tool call will use, and where it came from.

    ``explicit`` is the caller-supplied value (the only one safe to echo as caller input);
    ``effective`` is what the parsers get; ``source`` is one of the REFERENCE_SOURCE_* values;
    ``note`` explains why no auto reference was chosen, when that is the case.
    """

    explicit: str | None
    effective: str | None
    source: str
    note: str | None = None

    @property
    def echo(self) -> str | None:
        """The reference to hand back for viewer refetches: explicit or auto, never config.

        An auto reference is a public, allow-listed registry URL, so re-validating it as
        caller input on a refetch passes; the operator default stays server-side.
        """
        if self.source in (REFERENCE_SOURCE_EXPLICIT, REFERENCE_SOURCE_AUTO):
            return self.effective
        return None


async def _bam_contigs(file_path: str, config: BAMCPConfig) -> list[dict]:
    """The BAM/CRAM header's ``[{name, length}]`` list, from a bounded per-server cache.

    Callers validate ``file_path`` first. CRAM headers decode without a reference; the
    operator reference is passed through only as pysam's ``reference_filename``.
    """
    cache = get_services(config).header_cache
    key = (file_path, _local_fingerprint(file_path))
    cached = cache.get(key)
    if cached is not None:
        cache.move_to_end(key)
        return cached
    index_path = await _ensure_cached_index(file_path, config)
    # A cold remote BAM's first open gets the cold-start budget, not the parse timeout.
    await _warm_remote_inputs(file_path, None, index_path, config)
    contigs = await asyncio.wait_for(
        _run_in_parse_pool(_read_contigs_sync, file_path, config.reference, index_path),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )
    cache[key] = contigs
    while len(cache) > _HEADER_CACHE_MAX:
        cache.popitem(last=False)
    return contigs


async def _auto_reference(file_path: str, config: BAMCPConfig) -> tuple[str | None, str | None]:
    """Pick the registry FASTA matching the BAM's build AND contig style, or explain why not.

    Returns ``(url, note)``. The URL goes through the same :func:`validate_reference_path`
    as caller input (scheme, extension, port, the host allow-list, non-global IPs), so auto
    selection can never reach a host the operator has not allowed.
    """
    from .reference import contig_style, detect_genome_build, get_public_reference_url

    try:
        contigs = await _bam_contigs(file_path, config)
    except Exception as e:  # noqa: BLE001 — the main fetch reports the real error
        logger.info("auto-reference: header read failed for %s: %s", file_path, e)
        return None, "auto-reference skipped: the alignment header could not be read"
    build = detect_genome_build(contigs)["build"]
    if build == "unknown":
        return None, (
            "no reference: the genome build was not recognized from the header; "
            "pass 'reference' to enable SNV/indel detection"
        )
    style = contig_style(contigs)
    url = get_public_reference_url(build, style)
    if url is None:
        return None, (
            f"no reference: no public indexed {build} FASTA is known for {style}-style "
            "contigs; pass 'reference' to enable SNV/indel detection"
        )
    try:
        await asyncio.to_thread(validate_reference_path, url, config)
    except (ValueError, FileNotFoundError) as e:
        logger.info("auto-reference: %s rejected by validation: %s", url, e)
        return None, (
            f"no reference: the public {build} FASTA is not permitted by this server's "
            "remote-file policy; pass 'reference' to enable SNV/indel detection"
        )
    return url, None


async def _resolve_reference_info(
    args: dict[str, Any], config: BAMCPConfig, *, auto: bool = True
) -> ResolvedReference:
    """Resolve a tool call's reference: explicit arg > operator config > auto > none.

    Validates only a caller-supplied reference — that is the untrusted SSRF vector. A
    ``config.reference`` default is operator-set (trusted) and may legitimately be absent for BAM
    inputs or discovery paths like list_contigs, so it is neither required nor SSRF-checked here.

    With neither, and remote files enabled (the public FASTAs are remote), ``auto`` picks the
    registry FASTA matching ``args['file_path']``'s detected build and contig naming style
    (see :func:`_auto_reference`) so detection works on a bare BAM URL.
    """
    explicit_reference = args.get("reference")
    if explicit_reference:
        await asyncio.to_thread(validate_reference_path, explicit_reference, config)
        return ResolvedReference(explicit_reference, explicit_reference, REFERENCE_SOURCE_EXPLICIT)
    if config.reference:
        return ResolvedReference(None, config.reference, REFERENCE_SOURCE_CONFIG)
    file_path = args.get("file_path")
    if auto and config.allow_remote_files and file_path:
        url, note = await _auto_reference(file_path, config)
        if url is not None:
            return ResolvedReference(None, url, REFERENCE_SOURCE_AUTO)
        return ResolvedReference(None, None, REFERENCE_SOURCE_NONE, note)
    return ResolvedReference(None, None, REFERENCE_SOURCE_NONE)


CRAM_REFERENCE_REQUIRED_MESSAGE = (
    "CRAM decoding needs a reference on this server; pass reference=<indexed FASTA URL> "
    "(list_contigs suggests one for the file's build)"
)


def _require_cram_reference(file_path: str, reference: str | None, config: BAMCPConfig) -> None:
    """Refuse to decode a CRAM without a resolved reference on a remote-only server.

    With no ``reference_filename``, htslib resolves each contig's reference from the CRAM
    header itself — the ``@SQ UR:`` tag, a path or URL chosen by whoever wrote the file —
    and decodes read bases against it. ``REF_PATH``/``REF_CACHE`` (pinned local at startup)
    only govern the M5 lookup, not ``UR:``, so a caller-supplied CRAM could make the server
    read a local file (its bytes come back as read sequence) or fetch an arbitrary URL
    outside the host allow-list. When ``allow_local_files`` is off, require a reference
    BAMCP resolved itself (caller-supplied and validated, the operator's configured one,
    or the auto-selected registry FASTA). Header-only reads (list_contigs) never decode
    a slice and do not go through here. No-op when local files are allowed.
    """
    if config.allow_local_files or reference or not is_cram(file_path):
        return
    raise ToolInputError(CRAM_REFERENCE_REQUIRED_MESSAGE)


# (file_path, reference) -> header contigs htslib would resolve via @SQ UR: instead of the
# given reference. Bounded LRU: a header verdict is re-checked after eviction, never trusted
# beyond it. Only populated on remote-only servers.
_CRAM_HEADER_REF_CHECKS: OrderedDict[tuple[str, str], tuple[str, ...]] = OrderedDict()
_CRAM_HEADER_REF_CHECKS_MAX = 256


def _cram_header_only_refs(
    file_path: str, reference: str, index_path: str | None
) -> tuple[str, ...]:
    """Header contigs that carry a ``UR:`` tag but are absent from ``reference``.

    htslib falls back to a contig's header ``UR:`` path/URL when the given reference FASTA
    does not contain it (``cram_get_ref``: "Reference file given, but ref '…' not present"),
    so a resolved reference alone does not keep decoding off caller-chosen files. Reads only
    the header and the FASTA index — no slice is decoded. Runs in the parse pool.
    """
    from .parsers import _open_alignment

    with _open_alignment(file_path, reference, index_path) as samfile:
        sq = samfile.header.to_dict().get("SQ", [])
    with open_fasta(reference) as fasta:
        served = set(fasta.references)
    return tuple(str(s["SN"]) for s in sq if s.get("UR") and s.get("SN") not in served)


async def _check_cram_header_references(
    file_path: str, reference: str | None, config: BAMCPConfig
) -> None:
    """On a remote-only server, refuse a CRAM whose header could redirect decoding.

    Complements :func:`_require_cram_reference`: with a reference resolved, any header contig
    the reference lacks but whose ``@SQ`` names a ``UR:`` would still be decoded against that
    caller-chosen path or URL. No-op when local files are allowed or for BAM input.
    """
    if config.allow_local_files or not reference or not is_cram(file_path):
        return
    key = (file_path, reference)
    unserved = _CRAM_HEADER_REF_CHECKS.get(key)
    if unserved is None:
        index_path = await _ensure_cached_index(file_path, config)
        unserved = await asyncio.wait_for(
            _run_in_parse_pool(_cram_header_only_refs, file_path, reference, index_path),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
        _CRAM_HEADER_REF_CHECKS[key] = unserved
        while len(_CRAM_HEADER_REF_CHECKS) > _CRAM_HEADER_REF_CHECKS_MAX:
            _CRAM_HEADER_REF_CHECKS.popitem(last=False)
    else:
        _CRAM_HEADER_REF_CHECKS.move_to_end(key)
    if unserved:
        shown = ", ".join(unserved[:3]) + (
            f" (+{len(unserved) - 3} more)" if len(unserved) > 3 else ""
        )
        raise ToolInputError(
            f"This CRAM's header points contigs the reference does not contain ({shown}) at "
            "an external reference (@SQ UR:), which this server does not follow; pass a "
            "reference FASTA that covers every contig in the header"
        )


async def _resolve_reference(
    args: dict[str, Any], config: BAMCPConfig, *, auto: bool = True
) -> tuple[str | None, str | None]:
    """Return ``(explicit_reference, effective_reference)`` for a tool call.

    Thin wrapper over :func:`_resolve_reference_info` (which also reports the source). The
    explicit value is returned separately because only it may be echoed back to the caller.
    """
    resolved = await _resolve_reference_info(args, config, auto=auto)
    return resolved.explicit, resolved.effective


def _resolve_variant_source(variant_source: str, vcf_path: str | None) -> tuple[str | None, bool]:
    """Validate ``variant_source`` and return the effective ``(vcf_path, vcf_primary)``.

    - ``"auto"``: overlay a VCF (if given) on BAMCP candidates -> (vcf_path, False).
    - ``"vcf"``: the VCF is authoritative -> (vcf_path, True); requires vcf_path.
    - ``"bamcp"``: ignore any VCF -> (None, False).
    """
    if variant_source not in _VARIANT_SOURCES:
        raise ValueError(
            f"variant_source must be one of {_VARIANT_SOURCES}, got '{variant_source}'"
        )
    if variant_source == "vcf" and not vcf_path:
        raise ValueError("variant_source='vcf' requires a vcf_path")
    if variant_source == "bamcp":
        return None, False
    return vcf_path, variant_source == "vcf"


@dataclass(frozen=True)
class RegionRequest:
    """The validated, resolved inputs of a region tool call (see :func:`_region_request`).

    ``vcf_path``/``vcf_primary`` are the effective values from ``variant_source``; the
    thresholds are the per-call overrides falling back to the operator config.
    """

    file_path: str
    region: str
    resolved: ResolvedReference
    variant_source: str
    vcf_path: str | None
    vcf_primary: bool
    min_vaf: float
    min_depth: int
    min_alt_count: int

    @property
    def reference(self) -> str | None:
        """The reference the parsers get (explicit > config > auto)."""
        return self.resolved.effective


async def _region_request(
    args: dict[str, Any], config: BAMCPConfig, *, region: str | None = None
) -> RegionRequest:
    """The shared preamble of the region tools, in their order: validate the file (off the
    loop) and the region, resolve the reference, then the variant source and thresholds.

    ``region`` overrides ``args["region"]`` for tools that build it (jump_to).
    """
    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    region = args["region"] if region is None else region
    validate_region(region)
    resolved = await _resolve_reference_info(args, config)
    variant_source = args.get("variant_source", "auto")
    vcf_path, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))
    return RegionRequest(
        file_path=file_path,
        region=region,
        resolved=resolved,
        variant_source=variant_source,
        vcf_path=vcf_path,
        vcf_primary=vcf_primary,
        min_vaf=args.get("min_vaf", config.min_vaf),
        min_depth=args.get("min_depth", config.min_depth),
        min_alt_count=args.get("min_alt_count", config.min_alt_count),
    )


# Bounded per-site evidence for get_variants: each BAMCP SNV candidate gets a pileup-backed
# strand/quality/position record (artifact_risk + confidence) — up to this many sites per call,
# so a noisy wide region can't turn one request into thousands of pileups. Past the cap the
# remaining candidates are flagged (applied_filters.per_site_evidence.capped).
GET_VARIANTS_EVIDENCE_MAX_SITES = 50

_UNSCORED_INDEL_NOTE = (
    "Indel read-level evidence is not computed by get_variants; use visualize_region or "
    "get_variant_curation_summary at this site."
)
_UNSCORED_CAP_NOTE = (
    f"Per-site evidence is computed for at most {GET_VARIANTS_EVIDENCE_MAX_SITES} candidates "
    "per call; narrow the region to score this one."
)
_UNSCORED_TIMEOUT_NOTE = "Per-site evidence timed out; retry or narrow the region."


def _apply_alt_count(data: RegionData, min_alt_count: int) -> RegionData:
    """``data`` with BAMCP candidates below ``min_alt_count`` dropped (never mutates the cache).

    Detection's VAF/depth floors already ran in the parser; the alt-read floor is applied here so
    the region cache stays keyed by the parser thresholds only. Same predicate as every other
    candidate filter (:func:`passes_candidate_filters`): an unknown alt_count is kept.
    """
    kept = [v for v in data.variants if passes_candidate_filters(v, min_alt_count=min_alt_count)]
    if len(kept) == len(data.variants):
        return data
    return replace(data, variants=kept)


def _unscored(variant: dict[str, Any], reason: str, data: RegionData) -> dict[str, Any]:
    """Score a candidate whose read-level evidence was not computed — honestly, as unknown.

    Context risks (homopolymer, low depth) still apply; the missing strand/quality evidence is
    itself a ``read_evidence_unavailable`` risk, and confidence is never better than low.
    """
    evidence = _not_computed_evidence(reason)
    risk = compute_artifact_risk(variant, evidence, data.reference_sequence, data.start)
    out = dict(variant)
    out["read_evidence_computed"] = False
    out["artifact_risk"] = risk
    out["confidence"] = compute_confidence(variant, evidence, risk)
    out["is_low_confidence"] = out["confidence"] == "low"
    return out


async def _score_candidates(
    file_path: str,
    region: str,
    reference: str | None,
    variants: list[dict[str, Any]],
    data: RegionData,
    config: BAMCPConfig,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attach artifact_risk + confidence to BAMCP candidates (readless get_variants path).

    SNVs (up to GET_VARIANTS_EVIDENCE_MAX_SITES) are pileup-annotated — htslib sees every read
    at the site — and scored exactly like VCF SNVs. Indels and over-cap SNVs are scored as
    "evidence not computed". ``variants`` must be private copies (the annotator writes into
    them). Returns ``(variants, per_site_evidence_summary)``.
    """
    snvs = [
        v for v in variants if v.get("source") == "bamcp" and len(v["ref"]) == 1 == len(v["alt"])
    ]
    to_score = snvs[:GET_VARIANTS_EVIDENCE_MAX_SITES]
    scored_ids = {id(v) for v in to_score}
    timed_out = False
    if to_score and reference:
        try:
            index_path = await _ensure_cached_index(file_path, config)
            await asyncio.wait_for(
                _run_in_parse_pool(
                    annotate_vcf_snv_support,
                    file_path,
                    region,
                    reference,
                    to_score,
                    config.min_mapq,
                    config.min_baseq,
                    index_path,
                    sources=("bamcp",),
                ),
                timeout=BAM_PARSE_TIMEOUT_SECONDS,
            )
        except asyncio.TimeoutError:
            timed_out = True
            scored_ids = set()

    out: list[dict[str, Any]] = []
    for v in variants:
        if v.get("source") != "bamcp":
            out.append(v)
        elif id(v) in scored_ids:
            enhanced, _ = _vcf_evidence(v, data.reference_sequence, data.start)
            enhanced["evidence_source"] = "pileup"
            enhanced["read_evidence_computed"] = True
            out.append(enhanced)
        elif len(v["ref"]) != len(v["alt"]):
            out.append(_unscored(v, _UNSCORED_INDEL_NOTE, data))
        else:
            out.append(
                _unscored(v, _UNSCORED_TIMEOUT_NOTE if timed_out else _UNSCORED_CAP_NOTE, data)
            )
    bamcp_candidates = sum(1 for v in variants if v.get("source") == "bamcp")
    unscored_indels = sum(
        1 for v in variants if v.get("source") == "bamcp" and len(v["ref"]) != len(v["alt"])
    )
    summary = {
        "candidates": bamcp_candidates,
        "scored_sites": len(scored_ids),
        "max_sites": GET_VARIANTS_EVIDENCE_MAX_SITES,
        "capped": len(snvs) > GET_VARIANTS_EVIDENCE_MAX_SITES,
        "unscored_capped": max(0, len(snvs) - GET_VARIANTS_EVIDENCE_MAX_SITES),
        "unscored_indels": unscored_indels,
        "timed_out": timed_out,
    }
    return out, summary


# -- Requested-site probe (curation / classify) ----------------------------------------------
# get_variants answers "what passes the calling filters here?". A curator asking about ONE named
# site needs the opposite: that site's evidence even when it fails those filters — a low-BQ,
# one-strand artifact is exactly what min_baseq/min_alt_count hide, so "not found" would be a
# useless (and misleading) answer for the demo artifact. The probe re-reads only the requested
# site with min_baseq=0 and no VAF/depth/alt-read floor, and says which filters it fails.

SITE_PROBE_NOTE = (
    "Evaluated at the requested site with min_baseq=0 and no VAF/depth/alt-read floor so its "
    "evidence can be judged; a site that fails the applied filters is NOT reported by "
    "get_variants."
)


def site_filter_check(
    variant: dict[str, Any],
    config: BAMCPConfig,
    *,
    unfiltered: dict[str, Any] | None = None,
    min_vaf: float | None = None,
    min_depth: int | None = None,
    min_alt_count: int | None = None,
) -> dict[str, Any]:
    """``passes_filters`` + ``failed_filters`` for a site, from its counts under the filters.

    ``variant`` carries the counts AT the applied filters (alt_count/depth/vaf); ``unfiltered``
    (from :func:`probe_site`) the min_baseq=0 counts, when the site was probed.
    """
    min_vaf = config.min_vaf if min_vaf is None else min_vaf
    min_depth = config.min_depth if min_depth is None else min_depth
    min_alt_count = config.min_alt_count if min_alt_count is None else min_alt_count
    alt_count = variant.get("alt_count")
    at_filters: dict[str, Any] = {
        # None = count unknown (never fails min_alt_count; see core.filters).
        "alt_count": int(alt_count) if alt_count is not None else None,
        "depth": int(variant.get("depth") or 0),
        "vaf": float(variant.get("vaf") or 0.0),
    }
    removed = 0
    if unfiltered is not None and at_filters["alt_count"] is not None:
        removed = max(0, int(unfiltered["alt_count"]) - at_filters["alt_count"])
    failed = failed_candidate_filters(
        at_filters, min_vaf=min_vaf, min_depth=min_depth, min_alt_count=min_alt_count
    )
    if failed and removed > 0:
        # The site fails a floor AND min_baseq dropped some of its alt reads: name it first.
        failed.insert(0, "min_baseq")
    check: dict[str, Any] = {
        "passes_filters": not failed,
        "failed_filters": failed,
        "at_filters": at_filters,
    }
    if unfiltered is not None:
        check["unfiltered"] = unfiltered
        check["low_baseq_alt_reads_removed"] = removed
        check["note"] = SITE_PROBE_NOTE
    return check


async def _pileup_site(
    file_path: str,
    contig: str,
    pos0: int,
    ref: str,
    alt: str,
    reference: str | None,
    config: BAMCPConfig,
    min_baseq: int,
) -> dict[str, Any]:
    """Pileup-annotate one SNV site at ``min_baseq`` (config MAPQ + read exclusions apply)."""
    site: dict[str, Any] = {
        "contig": contig,
        "position": pos0,
        "ref": ref,
        "alt": alt,
        "variant_kind": "snv",
        "source": "probe",
    }
    index_path = await _ensure_cached_index(file_path, config)
    await asyncio.wait_for(
        _run_in_parse_pool(
            annotate_vcf_snv_support,
            file_path,
            format_region(contig, pos0, pos0 + 1),
            reference,
            [site],
            config.min_mapq,
            min_baseq,
            index_path,
            sources=("probe",),
        ),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )
    return site


async def probe_site(
    file_path: str,
    data: RegionData,
    pos0: int,
    ref: str,
    alt: str,
    reference: str | None,
    config: BAMCPConfig,
) -> dict[str, Any] | None:
    """Evidence for the requested site even when it fails the calling filters.

    SNVs: two single-site pileups (min_baseq=0 and the configured min_baseq), scored with the
    same per-site evidence as get_variants (strand, alt base quality, read position, MAPQ).
    Indels: the window is re-read with no VAF/depth floor and the CIGAR evidence is used.
    Returns ``{"variant", "evidence", "filter_check"}`` (the variant 0-based, with
    artifact_risk/confidence), or None when no read supports the allele at all.
    """
    from ..analysis.evidence import build_indel_index, compute_indel_evidence

    contig = data.contig
    # The probe counts reads carrying ALT, so a REF that isn't the reference base would turn
    # every ordinary read into "support" for a variant that doesn't exist. Verify it first.
    ref_seq = data.reference_sequence
    offset = pos0 - data.start
    if not ref_seq or offset < 0 or offset + len(ref) > len(ref_seq):
        # Without the reference base the probe can't tell a real ALT from a mistyped REF.
        raise ReferenceUnavailableError(
            f"Can't verify REF '{ref}' at {contig}:{pos0 + 1} without a reference sequence, "
            "so the filter-independent site check was not run. Pass `reference` (or use a "
            "remote BAM so one is auto-selected)."
        )
    actual = ref_seq[offset : offset + len(ref)].upper()
    if actual != ref.upper():
        raise ToolInputError(
            f"REF '{ref}' does not match the reference at {contig}:{pos0 + 1} "
            f"(the reference has '{actual}'). Check the allele, the position (1-based) "
            "and the genome build."
        )
    if len(ref) == 1 == len(alt):
        # Both pileups at once (each is a remote open + seek; they share nothing).
        loose, strict = await asyncio.gather(
            _pileup_site(file_path, contig, pos0, ref, alt, reference, config, 0),
            _pileup_site(file_path, contig, pos0, ref, alt, reference, config, config.min_baseq),
        )
        if not loose.get("alt_count"):
            return None
        enhanced, evidence = _vcf_evidence(loose, data.reference_sequence, data.start)
        enhanced["source"] = "bamcp"
        enhanced["evidence_source"] = "pileup"
        enhanced["read_evidence_computed"] = True
        evidence["read_evidence_computed"] = True
        unfiltered = {
            "alt_count": int(loose["alt_count"]),
            "depth": int(loose["depth"]),
            "vaf": float(loose["vaf"]),
        }
        check = site_filter_check(strict, config, unfiltered=unfiltered)
        return {"variant": enhanced, "evidence": evidence, "filter_check": check}

    # Indel: alt-read support isn't base-quality-filtered, so only the VAF/depth floors hide it.
    loose_data = await _fetch_region_with_timeout(
        file_path,
        format_region(contig, data.start, data.end),
        reference,
        config,
        min_vaf=0.0,
        min_depth=1,
    )
    target = next(
        (
            v
            for v in loose_data.variants
            if v["position"] == pos0 and v["ref"].upper() == ref and v["alt"].upper() == alt
        ),
        None,
    )
    if target is None:
        return None
    evidence = compute_indel_evidence(build_indel_index(loose_data.reads), target)
    risk = compute_artifact_risk(
        target,
        evidence,
        loose_data.reference_sequence,
        loose_data.start,
        reads=None if loose_data.reads_truncated else loose_data.reads,
        min_baseq=loose_data.min_baseq,
    )
    variant = {
        **target,
        "artifact_risk": risk,
        "confidence": compute_confidence(target, evidence, risk),
        "read_evidence_computed": bool(evidence.get("read_evidence_computed", True)),
    }
    return {
        "variant": variant,
        "evidence": evidence,
        "filter_check": site_filter_check(target, config),
    }


def _cap_variants(variants: list[dict[str, Any]], cap: int) -> list[dict[str, Any]]:
    """At most ``cap`` variants, chosen by a deterministic priority, in their input order.

    Priority: VCF records first (the caller's authoritative calls), then BAMCP candidates by
    descending alt-read count, then descending VAF, then ascending (position, ref, alt). The
    kept variants stay in input order, so a list at or under the cap is returned unchanged.
    """
    if len(variants) <= cap:
        return variants
    ranked = sorted(
        range(len(variants)),
        key=lambda i: (
            variants[i].get("source") != "vcf",
            -(variants[i].get("alt_count") or 0),
            -(variants[i].get("vaf") or 0.0),
            variants[i]["position"],
            str(variants[i].get("ref", "")),
            str(variants[i].get("alt", "")),
        ),
    )
    keep = sorted(ranked[:cap])
    return [variants[i] for i in keep]


def _select_get_variants(
    data: RegionData, min_vaf: float, min_depth: int, min_alt_count: int
) -> tuple[list[dict[str, Any]], int]:
    """Filter, cap (GET_VARIANTS_MAX_RESULTS) and VCF-enhance get_variants' candidates.

    Returns ``(private copies of the kept variants, total candidates passing the filters)``.
    VCF records always pass the filters and keep their annotated read support, so capping
    before :func:`enhance_vcf_variants` gives the same records as enhancing first.
    """
    passing = [
        v
        for v in data.variants
        if passes_candidate_filters(
            v, min_vaf=min_vaf, min_depth=min_depth, min_alt_count=min_alt_count
        )
    ]
    capped = _cap_variants(passing, GET_VARIANTS_MAX_RESULTS)
    detected = enhance_vcf_variants(capped, data.reference_sequence, data.start)
    return [dict(v) for v in detected], len(passing)


def _finish_get_variants(scored: list[dict[str, Any]], payload: dict[str, Any]) -> str:
    """Interpret, sort by position and JSON-encode get_variants' output (CPU-bound)."""
    variants = sorted(
        (annotate_interpretation(_one_based(v)) for v in scored),
        key=lambda v: (v["position"], v["ref"], v["alt"]),
    )
    payload = {"variants": variants, "count": len(variants), **payload}
    return json.dumps(payload)


_RESULTS_TRUNCATION_NOTE = (
    f"More than {GET_VARIANTS_MAX_RESULTS} variants passed the filters; only "
    f"{GET_VARIANTS_MAX_RESULTS} are listed (VCF records first, then BAMCP candidates with the "
    "most alt reads). Narrow the region or tighten min_vaf/min_depth/min_alt_count to see the rest."
)


@telemetry_wrapper("get_variants")
async def handle_get_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return variants without UI. Positions are 1-based (VCF/dbSNP convention).

    ``variant_source`` selects where variants come from:
      - ``"auto"`` (default): BAMCP's read-level candidates, with a VCF overlaid
        (de-duplicated) when ``vcf_path`` is given.
      - ``"vcf"``: the caller's VCF is the authoritative variant set — BAMCP does
        not add its own candidates but still reads the BAM to attach read-level
        evidence at each VCF site. Requires ``vcf_path``.
      - ``"bamcp"``: BAMCP's read-level candidates only; any ``vcf_path`` is ignored.

    Every variant carries ``artifact_risk``, ``confidence``, ``concerns`` and a plain-language
    ``interpretation`` (with a zygosity hint).

    At most GET_VARIANTS_MAX_RESULTS variants are returned (see :func:`_cap_variants` for which
    are kept); applied_filters then says ``truncated`` with ``total_candidates`` and the cap.
    The filter/score/encode work runs off the event loop.
    """
    req = await _region_request(args, config)
    file_path, region, reference, resolved = req.file_path, req.region, req.reference, req.resolved
    min_vaf, min_depth, min_alt_count = req.min_vaf, req.min_depth, req.min_alt_count

    # Readless: BAMCP candidates come from count_coverage, and VCF sites get their
    # read-level support attached by the fetch (annotate_vcf_snv_support), which
    # counts all reads at each site and so isn't affected by the max_reads read cap.
    data = await _fetch_region_with_timeout(
        file_path,
        region,
        reference,
        config,
        min_vaf=min_vaf,
        min_depth=min_depth,
        mode="variants",
        vcf_path=req.vcf_path,
        vcf_primary=req.vcf_primary,
    )

    # VCF SNVs get the viewer's confidence/artifact contract from their annotated support;
    # BAMCP candidates are then scored per site (bounded) by _score_candidates. Copies only:
    # the fetched variant dicts live in the region cache.
    kept, total_candidates = await asyncio.to_thread(
        _select_get_variants, data, min_vaf, min_depth, min_alt_count
    )
    scored, evidence_summary = await _score_candidates(
        file_path, region, reference, kept, data, config
    )

    filters = _applied_filters(
        config,
        args,
        include_max_reads=False,
        detection=_detection_status(region, reference, req.vcf_primary),
        min_alt_count=min_alt_count,
        reference_source=resolved.source,
    )
    filters["per_site_evidence"] = evidence_summary
    filters.update(_vcf_truncation_filters(data))
    filters["max_results"] = GET_VARIANTS_MAX_RESULTS
    filters["truncated"] = total_candidates > GET_VARIANTS_MAX_RESULTS
    filters["total_candidates"] = total_candidates
    if filters["truncated"]:
        filters["truncation_note"] = _RESULTS_TRUNCATION_NOTE
    if resolved.note:
        filters["reference_note"] = resolved.note
    payload: dict[str, Any] = {
        "variant_source": req.variant_source,
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
        "variant_type": "candidate",
        "applied_filters": filters,
        "filter_note": _FILTER_ABSENCE_NOTE,
        "disclaimer": _CANDIDATE_VARIANT_DISCLAIMER,
    }
    text = await asyncio.to_thread(_finish_get_variants, scored, payload)
    return {"content": [{"type": "text", "text": text}]}


@telemetry_wrapper("get_coverage")
async def handle_get_coverage(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return coverage statistics."""
    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    region = args["region"]
    validate_region(region)
    # Coverage never needs a reference for BAM (no base comparison); only CRAM decode does.
    _, reference = await _resolve_reference(args, config, auto=is_cram(file_path))

    data = await _fetch_region_with_timeout(file_path, region, reference, config, mode="coverage")

    coverage = data.coverage
    stats = {
        "region": format_region(data.contig, data.start, data.end),
        "mean": round(sum(coverage) / len(coverage), 2) if coverage else 0,
        "min": min(coverage) if coverage else 0,
        "max": max(coverage) if coverage else 0,
        "median": sorted(coverage)[len(coverage) // 2] if coverage else 0,
        "bases_covered": sum(1 for c in coverage if c > 0),
        "total_bases": len(coverage),
        # The echoed region string is 1-based inclusive, like the region it answers.
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
    }

    return {"content": [{"type": "text", "text": json.dumps(stats)}]}


def _read_contigs_sync(file_path: str, reference: str | None, index_path: str | None) -> list[dict]:
    """Read (name, length) contigs from a BAM/CRAM header. Runs in the parse pool."""
    mode = "rc" if is_cram(file_path) else "rb"
    with pysam.AlignmentFile(
        file_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference,
        index_filename=index_path,
    ) as samfile:
        return [
            {"name": name, "length": length}
            for name, length in zip(samfile.references, samfile.lengths, strict=True)
        ]


async def _detect_bam_build(file_path: str, config: BAMCPConfig) -> dict:
    """Detect a BAM/CRAM's genome build from its header. Returns detect_genome_build()'s dict."""
    from .reference import detect_genome_build

    await asyncio.to_thread(validate_path, file_path, config)
    return detect_genome_build(await _bam_contigs(file_path, config))


# Primary assembled chromosomes in either naming style (chr1..chr22, X, Y, M/MT).
_PRIMARY_CONTIG_RE = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")
# With no recognizable primary contigs (non-human / unusual naming), show at most this many.
_LIST_CONTIGS_FALLBACK_MAX = 100


def _contig_category(name: str) -> str:
    """Bucket a contig name: primary | alt | decoy | unplaced | hla | other."""
    if _PRIMARY_CONTIG_RE.match(name):
        return "primary"
    lower = name.lower()
    if lower.endswith("_alt"):
        return "alt"
    if "decoy" in lower or lower == "hs37d5" or lower.startswith("chrun_jtfh"):
        return "decoy"
    if lower.startswith("hla-"):
        return "hla"
    if (
        "random" in lower
        or lower.startswith("chrun")
        or lower.startswith(("gl", "kn", "ki", "jh", "nt_"))
    ):
        return "unplaced"
    return "other"


@telemetry_wrapper("list_contigs")
async def handle_list_contigs(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Summarize a BAM/CRAM's contigs: build, naming style and reference first, then contigs.

    By default only the primary chromosomes are listed (an HG002 GRCh38 header carries ~2,600
    alt/decoy/HLA contigs, which buried the build at the end of a ~157KB response); the rest
    are counted by category. ``include_all=true`` lists every contig.
    """
    from .reference import (
        CONTIG_STYLE_CHR,
        contig_style,
        detect_genome_build,
        get_public_reference_url,
    )

    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    include_all = bool(args.get("include_all", False))
    # Discovery path: never auto-pick a reference here — reporting the suggestion IS the job.
    resolved = await _resolve_reference_info(args, config, auto=False)
    reference = resolved.effective

    # Download and cache index for remote files (if not already cached)
    index_path = await _ensure_cached_index(file_path, config)
    await _warm_remote_inputs(file_path, reference, index_path, config)
    contigs = await asyncio.wait_for(
        _run_in_parse_pool(_read_contigs_sync, file_path, reference, index_path),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )

    # Detect genome build from contig lengths
    build_info = detect_genome_build(contigs)

    # Suggest a public reference URL only when there is no EFFECTIVE reference —
    # `reference` above, not just config.reference. A caller who passed a valid
    # `reference` argument already has one (it was just used to read this header),
    # so telling them to supply one is advice they have already taken. The
    # suggestion must also match the BAM's contig naming — a "chr1"-style FASTA
    # cannot serve a BAM whose contigs are "1" — and when no verified reference
    # exists for that style we say so instead of offering one that fails to open.
    style = contig_style(contigs)
    suggested_url = None
    reference_note = None
    if not reference and build_info["build"] != "unknown":
        suggested_url = get_public_reference_url(build_info["build"], style)
        if suggested_url is None:
            reference_note = (
                f"No public {build_info['build']} FASTA with a faidx index is known for "
                f"{'chr-prefixed' if style == CONTIG_STYLE_CHR else 'unprefixed'} contigs. "
                "Supply your own via the 'reference' argument or BAMCP_REFERENCE — it must be "
                "uncompressed or bgzip-compressed (not plain gzip) and have a .fai beside it."
            )
        elif config.allow_remote_files:
            reference_note = (
                "Variant tools use this reference automatically when you pass none "
                "(reference_source: auto)."
            )

    counts: dict[str, int] = {}
    for c in contigs:
        cat = _contig_category(str(c["name"]))
        counts[cat] = counts.get(cat, 0) + 1
    primary = [c for c in contigs if _contig_category(str(c["name"])) == "primary"]
    if include_all:
        shown, shown_kind = contigs, "all"
    elif primary:
        shown, shown_kind = primary, "primary"
    else:
        shown = contigs[:_LIST_CONTIGS_FALLBACK_MAX]
        shown_kind = "all" if len(shown) == len(contigs) else "first"

    # Summary first: a caller (or a truncating client) sees the build and the reference
    # guidance before any contig list.
    payload: dict[str, Any] = {
        "genome_build": build_info,
        "contig_style": style,
        "reference_configured": config.reference is not None,
        "suggested_reference_url": suggested_url,
        "reference_note": reference_note,
        "contig_counts": {"total": len(contigs), **dict(sorted(counts.items()))},
        "contigs_shown": shown_kind,
        "contigs_omitted": len(contigs) - len(shown),
        "contigs": shown,
        # Lengths are base counts; any region you build from them is 1-based inclusive.
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
    }
    if len(shown) < len(contigs):
        payload["contigs_note"] = (
            f"{len(contigs) - len(shown)} non-primary contigs omitted "
            "(alt/decoy/HLA/unplaced); pass include_all=true to list every contig."
        )
    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


@telemetry_wrapper("search_gene")
async def handle_search_gene(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Search a gene by symbol and return build-consistent genomic coordinates.

    Coordinates MUST match the build of the BAM being viewed — a GRCh38 locus in a GRCh37 BAM is
    silently ~0.5-2 Mb off the real gene. Build resolution:
      1. an explicit ``build`` argument overrides everything;
      2. else, if ``file_path`` is given, the build is auto-detected from that BAM's header;
      3. else, the server default (``config.genome_build``) is used.
    The response always states the coordinates' build, the source of that choice, and a loud
    warning on any mismatch against the BAM.
    """
    from .reference import normalize_build_name
    from .validation import validate_gene_symbol

    def _answer(status: str, **fields: Any) -> dict:
        # Not-found / unavailable are ANSWERS (no isError), shaped like the lookup tools':
        # not_found carries found=false + message; unavailable carries an error string.
        return {"content": [{"type": "text", "text": json.dumps({**fields, "status": status})}]}

    symbol = args.get("symbol", "")
    if not symbol:
        raise ToolInputError("search_gene requires 'symbol'")
    validate_gene_symbol(symbol)  # Entrez query-term injection guard

    file_path = args.get("file_path")
    requested = normalize_build_name(args["build"]) if args.get("build") else None
    if args.get("build") and requested is None:
        raise ToolInputError(
            f"Unrecognized build {args['build']!r}; use GRCh37/hg19 or GRCh38/hg38"
        )

    warnings: list[str] = []
    detected_build: str | None = None
    if file_path:
        try:
            info = await _detect_bam_build(file_path, config)
        except (ValueError, FileNotFoundError) as e:
            # A bad / disallowed / missing file_path is a caller error — never silently fall
            # back to the server default as if no BAM had been named.
            raise ToolInputError(f"Invalid file_path: {e}") from e
        except Exception as e:  # noqa: BLE001 — unreadable header: reconciliation impossible
            if requested is None:
                return _answer(
                    _LOOKUP_UNAVAILABLE,
                    error=f"could not read the BAM header to detect its build "
                    f"({type(e).__name__}); pass 'build' explicitly (confirm it with list_contigs)",
                )
            info = {}
            warnings.append(
                f"could not detect the BAM's build ({type(e).__name__}); coordinates unreconciled"
            )
        if info.get("build") and info["build"] != "unknown":
            detected_build = info["build"]
        elif info:
            # Detection ran but was inconclusive: whatever build we fall back to is
            # unreconciled with this BAM, and the caller must hear it.
            warnings.append(
                "could not determine the BAM's build from its header; coordinates are "
                "UNRECONCILED with the BAM — confirm its build with list_contigs"
                + ("" if requested else " or pass 'build'")
            )

    # Resolution: explicit override > detected-from-BAM > server default.
    if requested is not None:
        effective, source = requested, "explicit"
        if detected_build and detected_build != requested:
            warnings.append(
                f"you requested {requested} but the BAM is {detected_build} — returning "
                f"{requested} coordinates; they will NOT line up in that BAM"
            )
    elif detected_build is not None:
        effective, source = detected_build, "detected-from-bam"
    else:
        effective, source = config.genome_build, "server-default"
        if not file_path:
            warnings.append(
                f"no file_path provided to reconcile against; coordinates are {effective} "
                "(server default) — confirm your BAM's build with list_contigs"
            )

    try:
        result = await get_gene_client(config).search(symbol, build=effective)
    except GeneLookupError as e:
        # NCBI unreachable / in-band error: the gene's existence is UNKNOWN, not disproven.
        logger.warning("Gene lookup unavailable for %s: %s", symbol, e)
        return _answer(
            _LOOKUP_UNAVAILABLE,
            error="Gene lookup unavailable (NCBI error); retry later",
            reason=str(e),
            genome_build=effective,
        )
    if result is None:
        return _answer(
            _LOOKUP_NOT_FOUND,
            found=False,
            message=f"Gene {symbol!r} not found ({effective})",
            genome_build=effective,
        )

    # NCBI may not carry the requested assembly → result.build is what we actually got.
    if result.build == "unknown":
        warnings.append(
            "NCBI did not identify the assembly of these coordinates — they are UNVERIFIED "
            f"against {effective}; confirm before visualizing"
        )
    elif result.build != effective:
        warnings.append(f"NCBI returned {result.build} coordinates ({effective} unavailable)")
    if detected_build and result.build != detected_build:
        warnings.append(
            f"MISMATCH: these {result.build} coordinates do NOT match the BAM's "
            f"{detected_build} build — visualizing them will show the wrong locus"
        )
    note = getattr(result, "note", None)
    if note:
        warnings.append(note)

    # NCBI names chromosomes UCSC-style (chr13); a 1000G/Ensembl BAM calls it 13. Translate
    # to the BAM's own spelling so the returned region opens in that BAM as-is.
    chrom = result.chrom
    mapped_from: str | None = None
    if file_path:
        chrom, mapped_from, contig_warning = await _map_gene_contig(file_path, chrom, config)
        if contig_warning:
            warnings.append(contig_warning)

    payload: dict[str, Any] = {
        "symbol": result.symbol,
        "name": result.name,
        # GeneClient already converts NCBI's 0-based chrstart/chrstop to 1-based
        # inclusive, which is exactly the region string every region-taking tool reads.
        "region": f"{chrom}:{result.start}-{result.end}",
        "strand": result.strand,
        "genome_build": result.build,
        "build_source": source,
        "bam_build": detected_build,
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
    }
    if mapped_from is not None:
        payload["contig_mapped_from"] = mapped_from
    if warnings:
        payload["warning"] = "; ".join(warnings)
    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


async def _map_gene_contig(
    file_path: str, chrom: str, config: BAMCPConfig
) -> tuple[str, str | None, str | None]:
    """Map a gene's contig to the BAM's naming style: ``(contig, mapped_from, warning)``.

    Best-effort: an unreadable header leaves the contig unchanged (build detection above has
    already reported that); a contig absent in both spellings becomes a warning.
    """
    from .reference import alternate_contig_names, map_contig_to_style

    try:
        await asyncio.to_thread(validate_path, file_path, config)
        names = {str(c["name"]) for c in await _bam_contigs(file_path, config)}
    except Exception as e:  # noqa: BLE001 — mapping is advisory
        logger.debug("search_gene: contig mapping skipped for %s: %s", file_path, e)
        return chrom, None, None
    mapped = map_contig_to_style(chrom, names)
    if mapped is None:
        tried = " or ".join([repr(chrom), *(repr(a) for a in alternate_contig_names(chrom))])
        return chrom, None, f"contig {tried} is not in the BAM header — check list_contigs"
    if mapped != chrom:
        return mapped, chrom, None
    return chrom, None, None


async def _first_contig(file_path: str, reference: str | None, config: BAMCPConfig) -> str:
    """The BAM/CRAM header's first contig — the default when a caller names none.

    A hard-coded ``"chr1"`` default breaks on ``1``-style (Ensembl/GRCh37) BAMs, so resolve
    the default from the file itself.
    """
    index_path = await _ensure_cached_index(file_path, config)
    await _warm_remote_inputs(file_path, reference, index_path, config)
    contigs = await asyncio.wait_for(
        _run_in_parse_pool(_read_contigs_sync, file_path, reference, index_path),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )
    if not contigs:
        raise ToolInputError("The alignment file's header lists no contigs.")
    return str(contigs[0]["name"])


# The viewer tools' text content is what the LLM reads; the full payload (reads, coverage,
# evidence) rides in _meta / structuredContent for the UI. Keep the text a compact summary.
VIEWER_SUMMARY_MAX_CHARS = 1500
_VIEWER_SUMMARY_MAX_VARIANTS = 5
_VIEWER_SUMMARY_LINE_MAX = 220


def _clip(text: str, limit: int) -> str:
    return text if len(text) <= limit else text[: limit - 1] + "…"


def _viewer_summary(
    headline: str, data: RegionData, variants: list[dict], filters: dict[str, Any]
) -> str:
    """A compact (< VIEWER_SUMMARY_MAX_CHARS) text summary of a viewer response."""
    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 1) if coverage else 0
    lines = [
        f"{headline}: {len(data.reads)} reads shown, mean depth {mean_cov}x, "
        f"{len(variants)} candidate variants; reference: "
        f"{filters.get('reference_source', 'n/a')}."
    ]
    ranked = sorted(
        variants,
        key=lambda v: (v.get("confidence") != "high", -(v.get("alt_count") or 0)),
    )
    for v in ranked[:_VIEWER_SUMMARY_MAX_VARIANTS]:
        lines.append("  " + _clip(str(v.get("interpretation", "")), _VIEWER_SUMMARY_LINE_MAX))
    if len(variants) > _VIEWER_SUMMARY_MAX_VARIANTS:
        lines.append(
            f"  …and {len(variants) - _VIEWER_SUMMARY_MAX_VARIANTS} more (get_variants lists all)."
        )
    filters_line = _clip("Applied filters: " + _format_filters_line(filters), 400)
    tail = [
        filters_line,
        "Reads, coverage and per-variant evidence are in the viewer, not repeated here.",
    ]
    text = "\n".join(lines + tail)
    # Hard bound: drop variant lines (never the headline/filters) until it fits.
    while len(text) > VIEWER_SUMMARY_MAX_CHARS and len(lines) > 1:
        lines.pop()
        text = "\n".join(lines + tail)
    return _clip(text, VIEWER_SUMMARY_MAX_CHARS)


async def _viewer_response(
    args: dict[str, Any],
    config: BAMCPConfig,
    req: RegionRequest,
    summary_prefix: Callable[[RegionData], str],
) -> dict:
    """Fetch a region and build the MCP Apps viewer response (shared by visualize/jump_to).

    The text content (what the LLM reads) is a compact summary that echoes the applied
    filters; the full payload goes only in ``_meta.ui/init``. The viewer tools take no
    per-call alt-read floor: the operator's ``min_alt_count`` applies.
    """
    file_path, region, reference, resolved = req.file_path, req.region, req.reference, req.resolved
    vcf_primary, effective_vcf = req.vcf_primary, req.vcf_path

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, vcf_path=effective_vcf, vcf_primary=vcf_primary
    )
    data = _apply_alt_count(data, config.min_alt_count)
    filters = _applied_filters(
        config,
        args,
        include_max_reads=True,
        detection=_detection_status(region, reference, vcf_primary),
        reads_truncated=data.reads_truncated,
        min_alt_count=config.min_alt_count,
        reference_source=resolved.source,
    )
    filters.update(_vcf_truncation_filters(data))
    if resolved.note:
        filters["reference_note"] = resolved.note
    payload = serialize_region_data(data)
    # Internal 0-based positions in the payload; the interpretation text prints 1-based.
    payload["variants"] = [annotate_interpretation(v, one_based=False) for v in payload["variants"]]
    payload["coordinate_system"] = _INTERNAL_COORDINATE_SYSTEM
    payload["variant_type"] = "candidate"
    payload["applied_filters"] = filters
    payload["filter_note"] = _FILTER_ABSENCE_NOTE
    payload["disclaimer"] = _CANDIDATE_VARIANT_DISCLAIMER
    payload["file_path"] = file_path  # For client-side re-queries
    # Preserve the full data-source context so the viewer's pan/zoom/detail refetches carry it.
    # Without reference, follow-up calls lose mismatch evidence and CRAM decode fails on navigation.
    # Echo a caller-supplied or auto-selected reference (an allow-listed public URL, so the
    # refetch's caller-input validation accepts it). Never a config default: it is re-applied
    # server-side on refetch, and echoing it would make refetches re-validate the operator
    # default as caller input and fail when it's outside the caller allow-list.
    payload["variant_source"] = req.variant_source
    if resolved.echo is not None:
        payload["reference"] = resolved.echo
    if effective_vcf is not None:
        payload["vcf_path"] = effective_vcf

    summary = _viewer_summary(summary_prefix(data), data, payload["variants"], filters)
    return {
        "content": [{"type": "text", "text": summary}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


@telemetry_wrapper("jump_to")
async def handle_jump_to(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle jump_to tool call.

    Centers the viewer on a specific genomic position with a configurable window.
    ``contig`` is required: a position alone is ambiguous, and silently defaulting to the
    BAM's first contig showed chr1 reads for a chr17 position. The constructed region goes
    through the same ``validate_region`` as visualize_region, so both tools accept exactly
    the same contigs.
    """
    contig = args.get("contig")
    if not contig:
        raise ToolInputError(
            "jump_to needs a contig as well as a position (e.g. contig='chr17' or '17' — "
            "call list_contigs for this file's naming)."
        )
    position = args["position"]
    window = args.get("window", config.default_window)
    # ``position`` is 1-based; build a 1-based inclusive window of ``window`` bases
    # centred on it (start clamped at base 1).
    start = max(1, position - window // 2)
    end = start + max(window, 1) - 1
    req = await _region_request(args, config, region=f"{contig}:{start}-{end}")

    return await _viewer_response(
        args,
        config,
        req,
        lambda data: f"Jumped to {data.contig}:{position}",
    )


@telemetry_wrapper("visualize_region")
async def handle_visualize_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    """
    return await _viewer_response(
        args,
        config,
        await _region_request(args, config),
        lambda data: f"Region {format_region(data.contig, data.start, data.end)}",
    )


@telemetry_wrapper("get_region_summary")
async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata. Each variant line is followed by
    its plain-language interpretation (zygosity hint, artifact risk, confidence). At most
    GET_VARIANTS_MAX_RESULTS variants are listed (chosen as in get_variants, see
    :func:`_cap_variants`); a truncated listing says so in its filters line.
    """
    req = await _region_request(args, config)
    region, reference, resolved, min_alt_count = (
        req.region,
        req.reference,
        req.resolved,
        req.min_alt_count,
    )

    data = await _fetch_region_with_timeout(
        req.file_path,
        region,
        reference,
        config,
        vcf_path=req.vcf_path,
        vcf_primary=req.vcf_primary,
    )
    data = _apply_alt_count(data, min_alt_count)

    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 2) if coverage else 0
    max_cov = max(coverage) if coverage else 0

    filters = _applied_filters(
        config,
        args,
        include_max_reads=True,
        detection=_detection_status(region, reference, req.vcf_primary),
        reads_truncated=data.reads_truncated,
        min_alt_count=min_alt_count,
        reference_source=resolved.source,
    )
    filters.update(_vcf_truncation_filters(data))
    total_variants = len(data.variants)
    listed = data
    if total_variants > GET_VARIANTS_MAX_RESULTS:
        listed = replace(data, variants=_cap_variants(data.variants, GET_VARIANTS_MAX_RESULTS))
        filters["variants_listed"] = GET_VARIANTS_MAX_RESULTS
    summary_lines = [
        f"Region: {format_region(data.contig, data.start, data.end)}",
        f"Reads: {len(data.reads)}",
        f"Coverage: mean={mean_cov}x, max={max_cov}x",
        f"Candidate variants detected: {len(data.variants)}",
        f"Filters applied: {_format_filters_line(filters)}",
    ]
    if resolved.note:
        summary_lines.append(f"Reference: {resolved.note}")
    if listed is not data:
        summary_lines.append(
            f"Variant listing TRUNCATED: {GET_VARIANTS_MAX_RESULTS} of {total_variants} "
            "listed below (VCF records first, then the candidates with the most alt reads); "
            "narrow the region to see the rest."
        )
    summary_lines += [
        _FILTER_ABSENCE_NOTE,
        f"Coordinate system: {_LLM_COORDINATE_SYSTEM}",
        _CANDIDATE_VARIANT_DISCLAIMER,
    ]

    summary_lines.extend(await asyncio.to_thread(_summary_variant_lines, listed))
    return {"content": [{"type": "text", "text": "\n".join(summary_lines)}]}


def _summary_variant_lines(data: RegionData) -> list[str]:
    """get_region_summary's per-variant lines (evidence scoring is CPU-bound: run off-loop)."""
    lines: list[str] = []
    enhanced, _ = build_enhanced_variants(data)
    for v in enhanced:
        # 1-based position (VCF/dbSNP convention); SV type/span + per-sample genotypes
        # surfaced by _format_variant_line(), then the plain-language reading.
        lines.append(_format_variant_line(v))
        lines.append("    " + annotate_interpretation(v, one_based=False)["interpretation"])
    return lines


@telemetry_wrapper("cleanup_cache")
async def handle_cleanup_cache(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Clear downloaded index files and the region cache — single-client (stdio) servers only.

    One server holds ONE cache shared by every connected client, so on an HTTP transport any
    caller could wipe everyone's cache (and race in-flight index downloads). There, the cache is
    managed by its TTL and this is a no-op that says so. The shared HTTP clients are never
    closed by this tool (see :meth:`BAMCPServices.clear_caches`).
    """
    if config.transport != "stdio":
        text = (
            "Cache cleanup is disabled on shared (HTTP) servers: the cache is shared by all "
            f"clients and expires automatically after {config.cache_ttl}s."
        )
    else:
        removed = get_services(config).clear_caches()
        # The downloaded indexes are gone, so the next open of a remote input is cold again:
        # forget it was warm so that open gets the cold-start budget, not the parse timeout.
        # Left alone on purpose: _parse_gates (one semaphore per live loop, bounding in-flight
        # parses — replacing it while parses hold permits would let a second batch past the
        # bound) and _services_registry (the live server's services; dropping them would
        # close HTTP clients under in-flight lookups — clear_caches already empties their
        # caches).
        _clear_warm_inputs()
        text = f"Removed {removed} cached index files and cleared the region cache."
    return {"content": [{"type": "text", "text": text}]}

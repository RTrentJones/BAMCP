"""scan_variants: a time-bounded, chunked SNV scan of a contig (or a region of one).

The worker (:func:`_scan_contig_progressive`) runs in the parse pool and reports into a
lock-guarded :class:`_ScanProgress`, so a deadline yields a partial, self-describing answer.
"""

from __future__ import annotations

import asyncio
import json
import threading
import time
from dataclasses import dataclass, field
from typing import Any

import pysam

from ..analysis.evidence import annotate_interpretation
from ..config import BAMCPConfig
from ..constants import (
    SCAN_VARIANTS_CHUNK_SIZE,
    SCAN_VARIANTS_CHUNK_TARGET_SECONDS,
    SCAN_VARIANTS_DEADLINE_SECONDS,
    SCAN_VARIANTS_MAX_REGION,
    SCAN_VARIANTS_MAX_RESULTS,
    SCAN_VARIANTS_MIN_CHUNK_SIZE,
    SCAN_VARIANTS_TIMEOUT_SECONDS,
)
from ..middleware.telemetry import telemetry_wrapper
from .errors import ToolInputError
from .filters import passes_candidate_filters
from .parsers import format_region, open_fasta, parse_region
from .remote import _ensure_cached_index, _warm_remote_inputs
from .services import _run_in_parse_pool
from .validation import validate_path, validate_region


@dataclass
class _ScanProgress:
    """Shared, lock-guarded state of a progressive contig scan (worker thread → handler).

    Coordinates are internal 0-based half-open. ``scanned_end`` advances only after a chunk's
    candidates are appended, so ``[0, scanned_end)`` is exactly what ``variants`` covers.
    """

    variants: list[dict] = field(default_factory=list)
    scan_start: int = 0
    scanned_end: int = 0
    scan_end: int = 0
    # Where the caller asked the scan to stop (contig end or region end) — scan_end may be
    # earlier when max_region caps it; the gap is reported as truncation, never as "clean".
    requested_end: int = 0
    contig_length: int | None = None
    deadline_hit: bool = False
    stop: threading.Event = field(default_factory=threading.Event)
    lock: threading.Lock = field(default_factory=threading.Lock)


def _scan_chunk(
    samfile: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    contig: str,
    chunk_start: int,
    chunk_end: int,
    min_vaf: float,
    min_depth: int,
    min_mapq: int,
    min_baseq: int,
    min_alt_count: int | None,
) -> list[dict]:
    """SNV candidates in one chunk of a contig scan (count_coverage + deletion-depth fix)."""
    from .parsers import _deletion_depth_at, detect_variants, read_passes_filters

    ref_seq = fasta.fetch(contig, chunk_start, chunk_end)
    counts = samfile.count_coverage(
        contig,
        chunk_start,
        chunk_end,
        quality_threshold=min_baseq,
        read_callback=lambda r: read_passes_filters(r, min_mapq),
    )
    # Apply the alt-read floor up front: one- and two-read noise would otherwise each cost a
    # Python deletion pileup below before being discarded, eating the scan's time budget.
    variants = detect_variants(
        counts, ref_seq, contig, chunk_start, min_vaf, min_depth, min_alt_count=min_alt_count or 1
    )
    if variants:
        gaps = _deletion_depth_at(
            samfile, contig, sorted({v["position"] for v in variants}), min_mapq
        )
        for v in variants:
            gap = gaps.get(v["position"], 0)
            if gap:
                v["depth"] += gap
                v["vaf"] = round(v["alt_count"] / v["depth"], 3)
    # Adding deletion depth only lowers VAF: re-apply the floors (the shared predicate).
    return [
        v
        for v in variants
        if passes_candidate_filters(
            v, min_vaf=min_vaf, min_depth=min_depth, min_alt_count=min_alt_count
        )
    ]


def _scan_contig_progressive(
    bam_path: str,
    contig: str,
    reference_path: str,
    progress: _ScanProgress,
    deadline: float,
    *,
    chunk_size: int,
    min_vaf: float,
    min_depth: int,
    min_mapq: int,
    min_baseq: int,
    max_region: int,
    index_filename: str | None,
    min_alt_count: int | None = None,
    scan_start: int = 0,
    scan_stop: int | None = None,
) -> None:
    """Scan ``contig`` chunk by chunk into ``progress``, stopping at ``deadline`` (monotonic).

    ``[scan_start, scan_stop)`` (0-based half-open; stop defaults to the contig end) bounds the
    scan, capped at ``max_region`` bases.

    Runs in the parse pool. Checks the deadline (and the handler's ``stop`` flag, set when the
    hard timeout abandons it) before every chunk, so a slow remote scan ends cleanly with the
    chunks done so far rather than holding a worker until the whole contig is read.

    Raises:
        ToolInputError: If ``contig`` is not in the alignment header (the message suggests the
            file's own spelling, e.g. ``13`` for ``chr13``).
        ValueError: If the region starts past the contig end (sanitized by the server's
            ``_run``, like any other non-caller-facing failure).
    """
    from .parsers import _header_lengths, _open_alignment, contig_not_found_message

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        lengths = _header_lengths(samfile)
        contig_length = lengths.get(contig)
        if contig_length is None:
            # Unknown contig: an error, never a silent empty "no variants" result.
            raise ToolInputError(contig_not_found_message(contig, lengths))
        stop = contig_length if scan_stop is None else min(scan_stop, contig_length)
        if scan_start >= stop:
            # A start past the contig end would scan nothing and read as "no variants".
            raise ValueError(
                f"Region start {scan_start + 1:,} is beyond the end of {contig} "
                f"(length {contig_length:,} bp)"
            )
        scan_end = min(stop, scan_start + max_region)
        with progress.lock:
            progress.contig_length = contig_length
            progress.scan_start = scan_start
            progress.scanned_end = scan_start
            progress.scan_end = scan_end
            progress.requested_end = stop
        # Adaptive chunking: start small, double while chunks are fast, halve when slow.
        size = min(SCAN_VARIANTS_MIN_CHUNK_SIZE, chunk_size)
        chunk_start = scan_start
        with open_fasta(reference_path) as fasta:
            while chunk_start < scan_end:
                if progress.stop.is_set() or time.monotonic() >= deadline:
                    with progress.lock:
                        progress.deadline_hit = True
                    return
                chunk_end = min(chunk_start + size, scan_end)
                started = time.monotonic()
                found = _scan_chunk(
                    samfile,
                    fasta,
                    contig,
                    chunk_start,
                    chunk_end,
                    min_vaf,
                    min_depth,
                    min_mapq,
                    min_baseq,
                    min_alt_count,
                )
                with progress.lock:
                    progress.variants.extend(found)
                    progress.scanned_end = chunk_end
                chunk_start = chunk_end
                if time.monotonic() - started < SCAN_VARIANTS_CHUNK_TARGET_SECONDS:
                    size = min(size * 2, chunk_size)
                else:
                    size = max(size // 2, min(SCAN_VARIANTS_MIN_CHUNK_SIZE, chunk_size))


@telemetry_wrapper("scan_variants")
async def handle_scan_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Scan an entire contig for variants using fast coverage-based detection.

    The scan runs chunk by chunk against ``SCAN_VARIANTS_DEADLINE_SECONDS``; when the deadline
    (or the hard ``SCAN_VARIANTS_TIMEOUT_SECONDS`` backstop) hits, the candidates found so far
    are returned with ``truncated: true``, the 1-based ``scanned_region`` and a
    ``truncation_reason`` — a partial answer, never a bare timeout. Failures (no reference,
    unknown contig, nothing scanned at all) raise :class:`ToolInputError` so they reach the
    client as an MCP error result. With no ``contig``, the BAM's first contig is scanned.

    ``region`` (1-based ``contig:start-end``, no size cap) scans just that stretch — how a
    truncated scan is resumed: the truncation text and ``next_region`` give the exact remainder.
    Output is sorted by position (when capped, the highest-VAF candidates are the ones kept).
    """
    # The region tools' shared vocabulary lives in core.tools, which re-exports this module's
    # names; importing it lazily keeps the import graph acyclic at module load.
    from .tools import (
        _CANDIDATE_VARIANT_DISCLAIMER,
        _FILTER_ABSENCE_NOTE,
        _LLM_COORDINATE_SYSTEM,
        _applied_filters,
        _check_cram_header_references,
        _first_contig,
        _one_based,
        _resolve_reference_info,
    )

    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    resolved = await _resolve_reference_info(args, config)
    reference = resolved.effective
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)
    min_alt_count = args.get("min_alt_count", config.min_alt_count)

    if not reference:
        raise ToolInputError(
            "Reference genome required for variant scanning. "
            "Use list_contigs to detect genome build and get a reference URL."
        )

    index_path = await _ensure_cached_index(file_path, config)
    # Cold remote inputs get their one-off open under the cold budget, not the scan deadline.
    await _warm_remote_inputs(file_path, reference, index_path, config)
    await _check_cram_header_references(file_path, reference, config)
    scan_start, scan_stop = 0, None
    if args.get("region"):
        validate_region(args["region"])
        region_contig, scan_start, scan_stop = parse_region(args["region"], max_size=None)
        if args.get("contig") and args["contig"] != region_contig:
            raise ToolInputError(
                f"contig '{args['contig']}' does not match region '{args['region']}'; "
                "pass one or the other."
            )
        contig = region_contig
    else:
        contig = args.get("contig") or await _first_contig(file_path, reference, config)

    progress = _ScanProgress()
    deadline = time.monotonic() + SCAN_VARIANTS_DEADLINE_SECONDS
    hard_timeout = False
    try:
        await asyncio.wait_for(
            _run_in_parse_pool(
                _scan_contig_progressive,
                file_path,
                contig,
                reference,
                progress,
                deadline,
                chunk_size=SCAN_VARIANTS_CHUNK_SIZE,
                min_vaf=min_vaf,
                min_depth=min_depth,
                min_mapq=config.min_mapq,
                min_baseq=config.min_baseq,
                max_region=SCAN_VARIANTS_MAX_REGION,
                index_filename=index_path,
                min_alt_count=min_alt_count,
                scan_start=scan_start,
                scan_stop=scan_stop,
            ),
            timeout=SCAN_VARIANTS_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError:
        # A single chunk hung past the soft deadline: tell the worker to stop at its next
        # check and answer with what is already complete.
        progress.stop.set()
        hard_timeout = True
    # Any other failure (pysam/OS text, a ValueError) propagates to the server's _run, which
    # sanitizes it; the unknown-contig case already raised a caller-safe ToolInputError.

    with progress.lock:
        all_variants = list(progress.variants)
        scan_start = progress.scan_start
        scanned_end = progress.scanned_end
        scan_end = progress.scan_end
        requested_end = progress.requested_end
        contig_length = progress.contig_length
        deadline_hit = progress.deadline_hit

    if contig_length is None:
        # Never got as far as the header: there is no partial result to give.
        raise ToolInputError(
            f"Scan of {contig} timed out after {SCAN_VARIANTS_TIMEOUT_SECONDS:g}s before any "
            "data was read; retry (the first call to a remote file is the slowest) or use "
            "get_variants on a narrower region."
        )

    truncated_by_deadline = hard_timeout or deadline_hit or scanned_end < scan_end
    # max_region caps a single call; the uncapped remainder is unscanned, not variant-free.
    truncated_by_max_region = scan_end < requested_end
    truncated = truncated_by_deadline or truncated_by_max_region
    truncation_reason = None
    # The exact unscanned remainder (1-based inclusive) through the caller's requested stop.
    next_region = (
        format_region(contig, scanned_end, requested_end) if scanned_end < requested_end else None
    )
    if truncated:
        if truncated_by_deadline:
            budget = (
                SCAN_VARIANTS_TIMEOUT_SECONDS if hard_timeout else SCAN_VARIANTS_DEADLINE_SECONDS
            )
            cause = (
                f"time budget of {budget:g}s reached after scanning {scanned_end - scan_start:,} "
                f"of {requested_end - scan_start:,} bp"
            )
        else:
            cause = (
                f"one call scans at most {SCAN_VARIANTS_MAX_REGION:,} bp; scanned "
                f"{scanned_end - scan_start:,} of {requested_end - scan_start:,} bp"
            )
        truncation_reason = f"{cause}; candidates beyond scanned_region were NOT " + (
            f"examined — call scan_variants again with region='{next_region}' to scan the remainder"
            if next_region
            else "examined"
        )

    # When capped, keep the highest-VAF candidates; report them in genomic order.
    all_variants.sort(key=lambda v: v["vaf"], reverse=True)
    variants = sorted(all_variants[:SCAN_VARIANTS_MAX_RESULTS], key=lambda v: v["position"])
    # 1-based inclusive; None when not a single chunk completed.
    scanned_region = (
        format_region(contig, scan_start, scanned_end) if scanned_end > scan_start else None
    )

    filters = _applied_filters(
        config,
        args,
        include_max_reads=False,
        min_alt_count=min_alt_count,
        reference_source=resolved.source,
    )
    # Self-describing scope: what was scanned, and what this scan cannot see.
    filters.update(
        {
            "snv_detection": "enabled",
            "indel_detection": False,
            "indel_detection_note": "scan_variants is SNV-only; use get_variants on a "
            "smaller region to detect indels",
            "scanned_region": scanned_region,
            "contig_length": contig_length,
            "truncated_by_max_region": truncated_by_max_region,
            "per_site_evidence": "not computed by scan_variants (coverage-based scan); run "
            "get_variants on a small window around a candidate for artifact_risk/confidence",
            "max_region": SCAN_VARIANTS_MAX_REGION,
            "truncated_by_deadline": truncated_by_deadline,
            "total_candidates": len(all_variants),
            "results_capped": len(all_variants) > SCAN_VARIANTS_MAX_RESULTS,
            "max_results": SCAN_VARIANTS_MAX_RESULTS,
        }
    )
    if min_alt_count is not None:
        filters.setdefault("min_alt_count", min_alt_count)

    payload: dict[str, Any] = {
        "contig": contig,
        "variants": [annotate_interpretation(_one_based(v)) for v in variants],
        "count": len(variants),
        "truncated": truncated,
        "scanned_region": scanned_region,
    }
    if truncation_reason:
        payload["truncation_reason"] = truncation_reason
    if truncated and next_region:
        payload["next_region"] = next_region
    payload.update(
        {
            "coordinate_system": _LLM_COORDINATE_SYSTEM,
            "variant_type": "candidate",
            "applied_filters": filters,
            "filter_note": _FILTER_ABSENCE_NOTE,
            "disclaimer": _CANDIDATE_VARIANT_DISCLAIMER,
        }
    )
    return {"content": [{"type": "text", "text": json.dumps(payload)}]}

"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import json
import logging
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

import httpx
import pysam

from .cache import BAMIndexCache
from .clinvar import ClinVarClient
from .config import BAMCPConfig
from .constants import (
    ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD,
    ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD,
    ARTIFACT_NEAR_END_FRACTION_THRESHOLD,
    ARTIFACT_STRAND_BIAS_THRESHOLD,
    BAM_PARSE_TIMEOUT_SECONDS,
    DEFAULT_CONTIG,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
    LOW_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MIN_MEAN_QUALITY,
    LOW_CONFIDENCE_MIN_VAF,
    MAPQ_HISTOGRAM_BINS,
    POSITION_HISTOGRAM_BINS,
    QUALITY_HISTOGRAM_BINS,
    REMOTE_FILE_SCHEMES,
    VIEWER_RESOURCE_URI,
)
from .gnomad import GnomadClient
from .parsers import AlignedRead, RegionData, fetch_region

logger = logging.getLogger(__name__)

# Pattern for valid chromosome names (chr prefix optional, allows chr1-chr99 or 1-99, X, Y, M, MT)
CHROM_PATTERN = re.compile(r"^(chr)?(\d{1,2}|[XYM]|MT)$", re.IGNORECASE)

# Pattern for valid allele strings (only ACGTN)
ALLELE_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)

# Module-level cache instance for session consistency
_cache_instance: BAMIndexCache | None = None


def get_cache(config: BAMCPConfig) -> BAMIndexCache:
    """Get or create the session cache instance.

    Uses a module-level singleton to maintain consistent session ID
    across all tool calls within a server process.
    """
    global _cache_instance
    if _cache_instance is None:
        _cache_instance = BAMIndexCache(config.cache_dir, config.cache_ttl)
    return _cache_instance


def _get_index_path(file_path: str, config: BAMCPConfig) -> str | None:
    """Get cache path for a remote BAM's index file.

    Returns None for local files or if cache doesn't exist yet
    (let pysam handle auto-detection from URL).
    """
    cache = get_cache(config)
    index_path = cache.get_index_path(file_path)
    # Only return path if index is already cached; otherwise let pysam auto-detect
    if index_path and cache.is_valid(index_path):
        return index_path
    return None


async def _ensure_cached_index(file_path: str, config: BAMCPConfig) -> str | None:
    """Download and cache the BAM/CRAM index file if not already cached.

    For remote files, attempts to download the index (.bai/.crai) and store
    it in the session cache directory. This avoids repeated downloads within
    a session and allows pysam to use the local index.

    Args:
        file_path: Path or URL to the BAM/CRAM file.
        config: Server configuration.

    Returns:
        Path to cached index file, or None if:
        - File is local (no caching needed)
        - Download failed (pysam will try its own resolution)
    """
    cache = get_cache(config)
    index_path = cache.get_index_path(file_path)

    # Local file - no caching needed
    if index_path is None:
        return None

    # Already cached and valid
    if cache.is_valid(index_path):
        logger.debug("Using cached index: %s", index_path)
        return index_path

    # Determine index URL - try common extensions
    is_cram = file_path.endswith(".cram")
    if is_cram:
        index_urls = [file_path + ".crai"]
    else:
        # BAM files: try .bam.bai first, then .bai
        index_urls = [file_path + ".bai", file_path.rsplit(".", 1)[0] + ".bai"]

    logger.info("Downloading index for remote BAM: %s", file_path)

    async with httpx.AsyncClient(timeout=60.0, follow_redirects=True) as client:
        for index_url in index_urls:
            try:
                resp = await client.get(index_url)
                if resp.status_code == 200:
                    # Save to cache
                    Path(index_path).write_bytes(resp.content)
                    logger.info("Cached index (%d bytes): %s", len(resp.content), index_path)
                    return index_path
                elif resp.status_code == 404:
                    logger.debug("Index not found at %s, trying next", index_url)
                    continue
                else:
                    logger.warning("Index download failed (%d): %s", resp.status_code, index_url)
            except (httpx.RequestError, OSError) as e:
                logger.warning("Index download error for %s: %s", index_url, e)
                continue

    # All attempts failed - let pysam try its own resolution
    logger.warning("Could not download index for %s, falling back to pysam", file_path)
    return None


async def _fetch_region_with_timeout(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    min_vaf: float | None = None,
    min_depth: int | None = None,
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
    # Download and cache index for remote files (if not already cached)
    index_path = await _ensure_cached_index(file_path, config)

    return await asyncio.wait_for(
        asyncio.to_thread(
            fetch_region,
            file_path,
            region,
            reference,
            max_reads=config.max_reads,
            min_mapq=config.min_mapq,
            index_filename=index_path,
            min_vaf=min_vaf if min_vaf is not None else config.min_vaf,
            min_depth=min_depth if min_depth is not None else config.min_depth,
        ),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )


def validate_path(file_path: str, config: BAMCPConfig) -> None:
    """Validate that the file path is allowed by configuration.

    Args:
        file_path: Path or URL to validate.
        config: Server configuration.

    Raises:
        ValueError: If the path is not allowed.
    """
    # Check for remote URLs
    if "://" in file_path:
        if not config.allow_remote_files:
            raise ValueError(f"Remote files are disabled. Cannot access {file_path}")

        if not file_path.startswith(REMOTE_FILE_SCHEMES):
            raise ValueError(f"Scheme not supported for remote file: {file_path}")
        return

    # Check local files if restrictions are configured
    if config.allowed_directories:
        try:
            abs_path = Path(file_path).resolve()
        except OSError as e:
            raise ValueError(f"Invalid path: {file_path}") from e

        allowed = False
        for d in config.allowed_directories:
            try:
                allowed_dir = Path(d).resolve()
                if abs_path.is_relative_to(allowed_dir):
                    allowed = True
                    break
            except OSError:
                continue

        if not allowed:
            raise ValueError(
                f"Path {file_path} is not in allowed directories: {config.allowed_directories}"
            )


async def handle_get_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return variants without UI."""
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    reference = args.get("reference", config.reference)
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, min_vaf=min_vaf, min_depth=min_depth
    )

    variants = [v for v in data.variants if v["vaf"] >= min_vaf and v["depth"] >= min_depth]

    return {
        "content": [
            {"type": "text", "text": json.dumps({"variants": variants, "count": len(variants)})}
        ]
    }


async def handle_get_coverage(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return coverage statistics."""
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config)

    coverage = data.coverage
    stats = {
        "region": f"{data.contig}:{data.start}-{data.end}",
        "mean": round(sum(coverage) / len(coverage), 2) if coverage else 0,
        "min": min(coverage) if coverage else 0,
        "max": max(coverage) if coverage else 0,
        "median": sorted(coverage)[len(coverage) // 2] if coverage else 0,
        "bases_covered": sum(1 for c in coverage if c > 0),
        "total_bases": len(coverage),
    }

    return {"content": [{"type": "text", "text": json.dumps(stats)}]}


async def handle_list_contigs(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """List contigs in a BAM/CRAM file and detect genome build."""
    from .reference import detect_genome_build, get_public_reference_url

    file_path = args["file_path"]
    validate_path(file_path, config)
    reference = args.get("reference", config.reference)

    # Download and cache index for remote files (if not already cached)
    index_path = await _ensure_cached_index(file_path, config)

    def _list_contigs_sync() -> list[dict]:
        mode = "rc" if file_path.endswith(".cram") else "rb"
        # Use context manager to ensure file handles are closed on exception
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

    contigs = await asyncio.wait_for(
        asyncio.to_thread(_list_contigs_sync),
        timeout=BAM_PARSE_TIMEOUT_SECONDS,
    )

    # Detect genome build from contig lengths
    build_info = detect_genome_build(contigs)

    # Suggest public reference URL if no reference configured
    suggested_url = None
    if not config.reference and build_info["build"] != "unknown":
        suggested_url = get_public_reference_url(build_info["build"])

    return {
        "content": [
            {
                "type": "text",
                "text": json.dumps(
                    {
                        "contigs": contigs,
                        "genome_build": build_info,
                        "reference_configured": config.reference is not None,
                        "suggested_reference_url": suggested_url,
                    }
                ),
            }
        ]
    }


async def handle_jump_to(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle jump_to tool call.

    Centers the viewer on a specific genomic position with a configurable window.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    position = args["position"]
    contig = args.get("contig", DEFAULT_CONTIG)
    window = args.get("window", config.default_window)
    reference = args.get("reference", config.reference)

    start = max(0, position - window // 2)
    end = position + window // 2
    region = f"{contig}:{start}-{end}"

    data = await _fetch_region_with_timeout(file_path, region, reference, config)
    payload = _serialize_region_data(data)
    payload["file_path"] = file_path  # For client-side re-queries

    # Return summary text in content (for LLM context), full data only in _meta
    reads_count = len(data.reads)
    variants_count = len(data.variants)
    summary = f"Jumped to {data.contig}:{position}: {reads_count} reads, {variants_count} variants"
    return {
        "content": [{"type": "text", "text": summary}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


async def handle_visualize_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    Equivalent to browse_region but named to clarify its App-centric purpose.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config)
    payload = _serialize_region_data(data)
    payload["file_path"] = file_path  # For client-side re-queries

    # Return summary text in content (for LLM context), full data only in _meta
    reads_count = len(data.reads)
    variants_count = len(data.variants)
    summary = (
        f"Region {data.contig}:{data.start}-{data.end}: "
        f"{reads_count} reads, {variants_count} variants"
    )
    return {
        "content": [{"type": "text", "text": summary}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config)

    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 2) if coverage else 0
    max_cov = max(coverage) if coverage else 0

    summary_lines = [
        f"Region: {data.contig}:{data.start}-{data.end}",
        f"Reads: {len(data.reads)}",
        f"Coverage: mean={mean_cov}x, max={max_cov}x",
        f"Variants detected: {len(data.variants)}",
    ]

    for v in data.variants:
        summary_lines.append(
            f"  {v['contig']}:{v['position']} {v['ref']}>{v['alt']} "
            f"VAF={v['vaf']:.1%} depth={v['depth']}"
        )

    return {"content": [{"type": "text", "text": "\n".join(summary_lines)}]}


def _validate_variant_input(chrom: str, pos: int, ref: str, alt: str) -> str | None:
    """Validate variant lookup input parameters.

    Returns:
        Error message if validation fails, None if valid.
    """
    if not CHROM_PATTERN.match(chrom):
        return f"Invalid chromosome: {chrom}"
    if pos < 1:
        return f"Position must be positive, got {pos}"
    if not ALLELE_PATTERN.match(ref):
        return f"Invalid reference allele: {ref}"
    if not ALLELE_PATTERN.match(alt):
        return f"Invalid alternate allele: {alt}"
    if len(ref) > 1000 or len(alt) > 1000:
        return "Allele length exceeds maximum (1000bp)"
    return None


def validate_lookup_inputs(chrom: str, pos: int, ref: str, alt: str) -> None:
    """Validate variant lookup input parameters.

    Args:
        chrom: Chromosome (e.g., "chr17" or "17").
        pos: Genomic position (1-based).
        ref: Reference allele.
        alt: Alternate allele.

    Raises:
        ValueError: If any parameter is invalid.
    """
    if not CHROM_PATTERN.match(chrom):
        raise ValueError(f"Invalid chromosome: {chrom}")
    if pos < 1:
        raise ValueError(f"Position must be positive, got {pos}")
    if not ALLELE_PATTERN.match(ref):
        raise ValueError(f"Invalid reference allele: {ref}")
    if not ALLELE_PATTERN.match(alt):
        raise ValueError(f"Invalid alternate allele: {alt}")
    if len(ref) > 1000 or len(alt) > 1000:
        raise ValueError("Allele length exceeds maximum (1000bp)")


_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)


async def handle_lookup_clinvar(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Look up a variant in ClinVar via NCBI E-utilities.

    Returns clinical significance, review status, and associated conditions.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    # Validate input parameters
    validation_error = _validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": validation_error, "disclaimer": _CLINVAR_DISCLAIMER}
                    ),
                }
            ]
        }

    client = ClinVarClient(api_key=config.ncbi_api_key)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError) as e:
        # Network and HTTP errors - expected failures
        logger.warning("ClinVar lookup failed for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "ClinVar lookup failed", "disclaimer": _CLINVAR_DISCLAIMER}
                    ),
                }
            ]
        }
    except asyncio.TimeoutError:
        logger.warning("ClinVar lookup timed out for %s:%d %s>%s", chrom, pos, ref, alt)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "ClinVar lookup timed out", "disclaimer": _CLINVAR_DISCLAIMER}
                    ),
                }
            ]
        }

    if result is None:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {
                            "found": False,
                            "message": f"No ClinVar entry found for {chrom}:{pos} {ref}>{alt}",
                            "disclaimer": _CLINVAR_DISCLAIMER,
                        }
                    ),
                }
            ]
        }

    payload = asdict(result)
    payload["disclaimer"] = _CLINVAR_DISCLAIMER

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


async def handle_lookup_gnomad(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    # Validate input parameters
    validation_error = _validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": validation_error, "disclaimer": _GNOMAD_DISCLAIMER}
                    ),
                }
            ]
        }

    client = GnomadClient(dataset=config.gnomad_dataset)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError) as e:
        # Network and HTTP errors - expected failures
        logger.warning("gnomAD lookup failed for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "gnomAD lookup failed", "disclaimer": _GNOMAD_DISCLAIMER}
                    ),
                }
            ]
        }
    except asyncio.TimeoutError:
        logger.warning("gnomAD lookup timed out for %s:%d %s>%s", chrom, pos, ref, alt)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "gnomAD lookup timed out", "disclaimer": _GNOMAD_DISCLAIMER}
                    ),
                }
            ]
        }

    if result is None:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {
                            "found": False,
                            "message": f"No gnomAD entry found for {chrom}:{pos} {ref}>{alt}",
                            "disclaimer": _GNOMAD_DISCLAIMER,
                        }
                    ),
                }
            ]
        }

    payload = asdict(result)
    payload["disclaimer"] = _GNOMAD_DISCLAIMER

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


def _get_query_position(read: AlignedRead, ref_pos: int) -> int | None:
    """Get the query (read) position corresponding to a reference position.

    Uses CIGAR string to map reference coordinates to read coordinates.
    Returns None if the reference position is not aligned to the read.
    """
    import re

    cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", read.cigar)
    query_pos = 0
    ref_cursor = read.position

    for length_str, op in cigar_ops:
        length = int(length_str)

        if op in ("M", "=", "X"):
            # Match/mismatch: consumes both query and reference
            if ref_cursor <= ref_pos < ref_cursor + length:
                return query_pos + (ref_pos - ref_cursor)
            query_pos += length
            ref_cursor += length
        elif op == "I":
            # Insertion: consumes query only
            query_pos += length
        elif op in ("D", "N"):
            # Deletion/skip: consumes reference only
            if ref_cursor <= ref_pos < ref_cursor + length:
                return None  # Position is in a deletion
            ref_cursor += length
        elif op == "S":
            # Soft clip: consumes query only
            query_pos += length
        elif op == "H":
            # Hard clip: consumes nothing
            pass

    return None


def _bin_values(values: list[int], bins: list[int]) -> list[int]:
    """Bin values into histogram buckets.

    Args:
        values: List of values to bin.
        bins: List of bin boundaries (e.g., [0, 10, 20, 30, 40]).
              Values >= bins[i] and < bins[i+1] go into bin i.
              The last bin captures all values >= bins[-1].

    Returns:
        List of counts per bin, length = len(bins).
    """
    counts = [0] * len(bins)
    for v in values:
        for i in range(len(bins) - 1, -1, -1):
            if v >= bins[i]:
                counts[i] += 1
                break
    return counts


def _detect_homopolymer(seq: str, pos: int) -> int:
    """Detect homopolymer run length at a given position.

    Args:
        seq: Reference sequence.
        pos: Position within the sequence (0-indexed).

    Returns:
        Length of homopolymer run containing the position.
    """
    if not seq or pos < 0 or pos >= len(seq):
        return 0

    base = seq[pos].upper()
    if base not in "ACGT":
        return 0

    length = 1

    # Extend left
    i = pos - 1
    while i >= 0 and seq[i].upper() == base:
        length += 1
        i -= 1

    # Extend right
    i = pos + 1
    while i < len(seq) and seq[i].upper() == base:
        length += 1
        i += 1

    return length


def _build_mismatch_index(
    reads: list[AlignedRead],
) -> dict[tuple[int, str], list[tuple[AlignedRead, dict]]]:
    """Build an index of mismatches by (position, alt) for O(1) variant lookup.

    Returns:
        Dict mapping (position, alt) to list of (read, mismatch) tuples.
    """
    index: dict[tuple[int, str], list[tuple[AlignedRead, dict]]] = {}
    for read in reads:
        for mm in read.mismatches:
            key = (mm["pos"], mm["alt"])
            if key not in index:
                index[key] = []
            index[key].append((read, mm))
    return index


def _compute_variant_evidence_from_index(
    index: dict[tuple[int, str], list[tuple[AlignedRead, dict]]],
    variant: dict,
) -> dict:
    """Compute variant evidence using pre-built mismatch index.

    This is O(k) where k is the number of reads supporting the variant,
    rather than O(n*m) when iterating all reads for each variant.

    Returns evidence including histogram distributions for curation.
    """
    pos = variant["position"]
    alt = variant["alt"]

    matches = index.get((pos, alt), [])

    forward_count = 0
    reverse_count = 0
    qualities: list[int] = []
    positions_in_read: list[int] = []
    mapq_values: list[int] = []

    for read, mm in matches:
        if read.is_reverse:
            reverse_count += 1
        else:
            forward_count += 1

        # Get quality at mismatch position
        query_pos = _get_query_position(read, mm["pos"])
        if query_pos is not None and query_pos < len(read.qualities):
            qualities.append(read.qualities[query_pos])

            # Compute position relative to nearest read end
            read_length = len(read.qualities)
            dist_from_end = min(query_pos, read_length - query_pos - 1)
            positions_in_read.append(dist_from_end)

        # Collect MAPQ for all reads
        mapq_values.append(read.mapping_quality)

    # Strand bias: 0 = balanced, 1 = all one strand
    total = forward_count + reverse_count
    strand_bias = abs(forward_count - reverse_count) / max(total, 1)

    # Compute histogram distributions
    quality_histogram = _bin_values(qualities, QUALITY_HISTOGRAM_BINS)
    position_histogram = _bin_values(positions_in_read, POSITION_HISTOGRAM_BINS)
    mapq_histogram = _bin_values(mapq_values, MAPQ_HISTOGRAM_BINS)

    return {
        "forward_count": forward_count,
        "reverse_count": reverse_count,
        "strand_bias": round(strand_bias, 3),
        "mean_quality": round(sum(qualities) / len(qualities), 1) if qualities else 0,
        "median_quality": sorted(qualities)[len(qualities) // 2] if qualities else 0,
        "quality_histogram": quality_histogram,
        "position_histogram": position_histogram,
        "mapq_histogram": mapq_histogram,
    }


def _compute_artifact_risk(
    variant: dict,
    evidence: dict,
    reference_sequence: str | None,
    region_start: int,
) -> dict:
    """Compute artifact risk indicators for variant curation.

    Args:
        variant: Variant dict with position, ref, alt, depth, vaf.
        evidence: Evidence dict from _compute_variant_evidence_from_index.
        reference_sequence: Reference sequence string (may be None).
        region_start: Start position of the region.

    Returns:
        Dict with risks list, risk_score, and artifact_likelihood.
    """
    risks: list[dict] = []
    risk_score = 0.0

    # 1. Position-in-read bias (ReadPosRankSum equivalent)
    pos_hist = evidence.get("position_histogram", [])
    if pos_hist:
        total = sum(pos_hist)
        if total > 0:
            # First two bins (0-25, 25-50) represent near-end positions
            near_end = sum(pos_hist[:2])
            near_end_fraction = near_end / total
            if near_end_fraction > ARTIFACT_NEAR_END_FRACTION_THRESHOLD:
                risks.append(
                    {
                        "type": "read_position_bias",
                        "severity": "medium",
                        "description": f"{near_end_fraction:.0%} of bases near read ends",
                        "value": round(near_end_fraction, 3),
                    }
                )
                risk_score += 0.3

    # 2. Strand bias
    strand_bias = evidence.get("strand_bias", 0)
    if strand_bias > ARTIFACT_STRAND_BIAS_THRESHOLD:
        severity = "high" if strand_bias > 0.95 else "medium"
        risks.append(
            {
                "type": "strand_bias",
                "severity": severity,
                "description": f"Strand bias: {strand_bias:.0%}",
                "value": round(strand_bias, 3),
            }
        )
        risk_score += 0.4 if severity == "high" else 0.2

    # 3. Low MAPQ fraction
    mapq_hist = evidence.get("mapq_histogram", [])
    if mapq_hist:
        total = sum(mapq_hist)
        if total > 0:
            # First three bins (0-10, 10-20, 20-30) represent low MAPQ
            low_mapq = sum(mapq_hist[:3])
            low_mapq_fraction = low_mapq / total
            if low_mapq_fraction > ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD:
                risks.append(
                    {
                        "type": "low_mapq",
                        "severity": "medium",
                        "description": f"{low_mapq_fraction:.0%} reads with MAPQ < 30",
                        "value": round(low_mapq_fraction, 3),
                    }
                )
                risk_score += 0.25

    # 4. Homopolymer context
    if reference_sequence:
        pos_in_region = variant["position"] - region_start
        if 0 <= pos_in_region < len(reference_sequence):
            hp_len = _detect_homopolymer(reference_sequence, pos_in_region)
            if hp_len >= ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD:
                severity = "high" if hp_len >= 6 else "medium"
                risks.append(
                    {
                        "type": "homopolymer",
                        "severity": severity,
                        "description": f"In {hp_len}bp homopolymer run",
                        "value": hp_len,
                    }
                )
                risk_score += 0.3 if severity == "high" else 0.15

    # 5. Low depth
    depth = variant.get("depth", 0)
    if depth < LOW_CONFIDENCE_MIN_DEPTH:
        severity = "high" if depth < 5 else "medium"
        risks.append(
            {
                "type": "low_depth",
                "severity": severity,
                "description": f"Low coverage depth ({depth}x)",
                "value": depth,
            }
        )
        risk_score += 0.3 if severity == "high" else 0.15

    # Cap risk score at 1.0
    risk_score = min(risk_score, 1.0)

    # Determine artifact likelihood
    if risk_score >= 0.6:
        artifact_likelihood = "high"
    elif risk_score >= 0.3:
        artifact_likelihood = "medium"
    else:
        artifact_likelihood = "low"

    return {
        "risks": risks,
        "risk_score": round(risk_score, 2),
        "artifact_likelihood": artifact_likelihood,
    }


def _serialize_region_data(data: RegionData, compact: bool | None = None) -> dict:
    """Serialize RegionData to a JSON-compatible dict.

    Args:
        data: The region data to serialize.
        compact: If True, omit sequences to reduce payload size.
                 If None (default), auto-detect based on region size:
                 include sequences for regions <= 500bp (base-level view possible).
    """
    # Auto-detect compact mode based on region size
    # Base rendering requires scale >= 10, which means ~80bp at 800px width
    # Use 500bp threshold for safety margin
    if compact is None:
        region_size = data.end - data.start
        compact = region_size > 500
    # Build mismatch index once for O(1) variant evidence lookups
    mismatch_index = _build_mismatch_index(data.reads)

    # Compute variant evidence using index (O(k) per variant instead of O(n*m) total)
    variant_evidence = {}
    enhanced_variants = []

    for variant in data.variants:
        key = f"{variant['position']}:{variant['ref']}>{variant['alt']}"
        evidence = _compute_variant_evidence_from_index(mismatch_index, variant)

        # Compute artifact risk
        artifact_risk = _compute_artifact_risk(
            variant, evidence, data.reference_sequence, data.start
        )
        evidence["artifact_risk"] = artifact_risk
        variant_evidence[key] = evidence

        # Enhance variant with evidence data for table display
        enhanced = dict(variant)
        enhanced["strand_forward"] = evidence["forward_count"]
        enhanced["strand_reverse"] = evidence["reverse_count"]
        enhanced["mean_quality"] = evidence["mean_quality"]
        enhanced["artifact_risk"] = artifact_risk

        # Compute confidence level based on multiple criteria
        # High: VAF >= 20%, depth >= 20, quality >= 25, strand bias <= 0.7
        # Medium: VAF >= 10%, depth >= 10, quality >= 20, strand bias <= 0.9
        # Low: anything else
        vaf = variant["vaf"]
        depth = variant["depth"]
        quality = evidence["mean_quality"]
        strand_bias = evidence["strand_bias"]

        if vaf >= 0.2 and depth >= 20 and quality >= 25 and strand_bias <= 0.7:
            confidence = "high"
        elif (
            vaf >= LOW_CONFIDENCE_MIN_VAF
            and depth >= LOW_CONFIDENCE_MIN_DEPTH
            and quality >= LOW_CONFIDENCE_MIN_MEAN_QUALITY
            and strand_bias <= LOW_CONFIDENCE_MAX_STRAND_BIAS
        ):
            confidence = "medium"
        else:
            confidence = "low"

        # Downgrade confidence if artifact risk is high
        if artifact_risk["artifact_likelihood"] == "high":
            confidence = "low"
        elif artifact_risk["artifact_likelihood"] == "medium" and confidence == "high":
            confidence = "medium"

        enhanced["confidence"] = confidence
        # Keep is_low_confidence for backwards compatibility
        enhanced["is_low_confidence"] = confidence == "low"

        enhanced_variants.append(enhanced)

    # Serialize reads - compact mode omits sequences for smaller payload
    def serialize_read(r: AlignedRead, include_sequence: bool) -> dict:
        """Serialize a read, omitting null paired-end fields to reduce payload."""
        d: dict[str, Any] = {
            "name": r.name,
            "cigar": r.cigar,
            "position": r.position,
            "end_position": r.end_position,
            "mapping_quality": r.mapping_quality,
            "is_reverse": r.is_reverse,
            "mismatches": r.mismatches,
        }
        if include_sequence:
            d["sequence"] = r.sequence
        # Only include paired-end fields if the read is actually paired
        if r.is_paired:
            d["mate_position"] = r.mate_position
            d["mate_contig"] = r.mate_contig
            d["insert_size"] = r.insert_size
            d["is_proper_pair"] = r.is_proper_pair
            d["is_read1"] = r.is_read1
            d["is_paired"] = True
        # Include soft clips if present
        if r.soft_clips:
            d["soft_clips"] = [
                {
                    "position": sc.position,
                    "length": sc.length,
                    "sequence": sc.sequence,
                    "side": sc.side,
                }
                for sc in r.soft_clips
            ]
        return d

    reads_data = [serialize_read(r, include_sequence=not compact) for r in data.reads]

    return {
        "contig": data.contig,
        "start": data.start,
        "end": data.end,
        "reads": reads_data,
        "coverage": data.coverage,
        "variants": enhanced_variants,
        "reference_sequence": data.reference_sequence,
        "variant_evidence": variant_evidence,
        "compact": compact,
    }


def _compute_near_end_fraction(position_histogram: list[int] | None) -> float:
    """Compute fraction of variant bases near read ends."""
    if not position_histogram:
        return 0.0
    total = sum(position_histogram)
    if total == 0:
        return 0.0
    # First two bins represent near-end positions
    near_end = sum(position_histogram[:2])
    return round(near_end / total, 3)


def _generate_curator_recommendations(
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
        near_end_frac = _compute_near_end_fraction(position_hist)
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


def _format_curation_summary(summary: dict) -> str:
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


async def handle_get_variant_curation_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Generate a curation-focused summary for a specific variant.

    Returns structured data optimized for clinical interpretation,
    including artifact risk assessment, quality metrics, and
    recommendations for further review.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]
    window = args.get("window", 50)
    reference = args.get("reference", config.reference)

    # Validate inputs
    validation_error = _validate_variant_input(chrom, pos, ref, alt)
    if validation_error:
        return {"content": [{"type": "text", "text": validation_error}]}

    region = f"{chrom}:{pos - window}-{pos + window}"

    try:
        data = await _fetch_region_with_timeout(file_path, region, reference, config)
    except asyncio.TimeoutError:
        return {"content": [{"type": "text", "text": "Timeout fetching region data"}]}

    # Find the specific variant
    target_variant = None
    for v in data.variants:
        if v["position"] == pos and v["ref"] == ref and v["alt"] == alt:
            target_variant = v
            break

    if not target_variant:
        return {
            "content": [
                {
                    "type": "text",
                    "text": f"Variant {chrom}:{pos} {ref}>{alt} not found in region. "
                    f"Variants found: {len(data.variants)}",
                }
            ]
        }

    # Build mismatch index and compute evidence
    mismatch_index = _build_mismatch_index(data.reads)
    evidence = _compute_variant_evidence_from_index(mismatch_index, target_variant)
    artifact_risk = _compute_artifact_risk(
        target_variant, evidence, data.reference_sequence, data.start
    )

    # Compute confidence
    vaf = target_variant["vaf"]
    depth = target_variant["depth"]
    quality = evidence["mean_quality"]
    strand_bias = evidence["strand_bias"]

    if vaf >= 0.2 and depth >= 20 and quality >= 25 and strand_bias <= 0.7:
        confidence = "high"
    elif (
        vaf >= LOW_CONFIDENCE_MIN_VAF
        and depth >= LOW_CONFIDENCE_MIN_DEPTH
        and quality >= LOW_CONFIDENCE_MIN_MEAN_QUALITY
        and strand_bias <= LOW_CONFIDENCE_MAX_STRAND_BIAS
    ):
        confidence = "medium"
    else:
        confidence = "low"

    # Downgrade if artifact risk is high
    if artifact_risk["artifact_likelihood"] == "high":
        confidence = "low"
    elif artifact_risk["artifact_likelihood"] == "medium" and confidence == "high":
        confidence = "medium"

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
            "near_end_fraction": _compute_near_end_fraction(evidence.get("position_histogram")),
        },
        "artifact_assessment": artifact_risk,
        "confidence": confidence,
        "recommendations": _generate_curator_recommendations(
            target_variant, evidence, artifact_risk
        ),
    }

    # Format as readable text for LLM
    text_summary = _format_curation_summary(curation_summary)

    return {"content": [{"type": "text", "text": text_summary}]}

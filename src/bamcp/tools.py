"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import bisect
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
    BAM_PARSE_TIMEOUT_SECONDS,
    DEFAULT_CONTIG,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
    LOW_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MIN_MEAN_QUALITY,
    LOW_CONFIDENCE_MIN_VAF,
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

    Returns None for local files (let pysam find the index normally).
    """
    cache = get_cache(config)
    return cache.get_index_path(file_path)


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
    index_path = _get_index_path(file_path, config)

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


async def handle_browse_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle browse_region tool call.

    Returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config)
    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


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

    def _list_contigs_sync() -> list[dict]:
        mode = "rc" if file_path.endswith(".cram") else "rb"
        # Use context manager to ensure file handles are closed on exception
        with pysam.AlignmentFile(
            file_path,
            mode,  # type: ignore[arg-type]
            reference_filename=reference,
            index_filename=_get_index_path(file_path, config),
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

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
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

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
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


def _bin_histogram(values: list[int], bins: list[int]) -> dict[str, int]:
    """Bin values into a histogram with labeled ranges.

    Uses binary search for O(n log b) instead of O(n*b) linear scan.

    Args:
        values: List of integer values to bin.
        bins: List of bin edges (e.g., [0, 10, 20, 30, 40] creates bins 0-9, 10-19, etc.)

    Returns:
        Dict mapping bin labels to counts (e.g., {"0-9": 5, "10-19": 3, ...}).
    """
    # Pre-compute labels and initialize histogram
    labels = []
    for i in range(len(bins)):
        label = f"{bins[i]}+" if i == len(bins) - 1 else f"{bins[i]}-{bins[i + 1] - 1}"
        labels.append(label)

    # Use list for counting (faster than dict for small number of bins)
    counts = [0] * len(bins)

    # Use binary search to find bin for each value - O(n log b)
    for val in values:
        # bisect_right gives us the insertion point
        # Subtract 1 to get the bin index (values below first bin go to -1, clamp to 0)
        idx = bisect.bisect_right(bins, val) - 1
        if idx < 0:
            idx = 0  # Values below first bin go to first bin
        elif idx >= len(bins):
            idx = len(bins) - 1  # Values above last bin go to last bin
        counts[idx] += 1

    return dict(zip(labels, counts, strict=True))


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
    """
    pos = variant["position"]
    alt = variant["alt"]

    matches = index.get((pos, alt), [])

    forward_count = 0
    reverse_count = 0
    qualities: list[int] = []
    read_positions: list[int] = []

    for read, mm in matches:
        if read.is_reverse:
            reverse_count += 1
        else:
            forward_count += 1

        # Get quality at mismatch position
        query_pos = _get_query_position(read, mm["pos"])
        if query_pos is not None and query_pos < len(read.qualities):
            qualities.append(read.qualities[query_pos])
            # Position from 5' end (adjusted for reverse reads)
            pos_from_5prime = len(read.sequence) - query_pos - 1 if read.is_reverse else query_pos
            read_positions.append(pos_from_5prime)

    # Strand bias: 0 = balanced, 1 = all one strand
    total = forward_count + reverse_count
    strand_bias = abs(forward_count - reverse_count) / max(total, 1)

    return {
        "forward_count": forward_count,
        "reverse_count": reverse_count,
        "strand_bias": round(strand_bias, 3),
        "mean_quality": round(sum(qualities) / len(qualities), 1) if qualities else 0,
        "median_quality": sorted(qualities)[len(qualities) // 2] if qualities else 0,
        "quality_histogram": _bin_histogram(qualities, QUALITY_HISTOGRAM_BINS),
        "position_histogram": _bin_histogram(read_positions, POSITION_HISTOGRAM_BINS),
    }


def _serialize_region_data(data: RegionData, compact: bool = False) -> dict:
    """Serialize RegionData to a JSON-compatible dict.

    Args:
        data: The region data to serialize.
        compact: If True, omit sequences/qualities to reduce payload size.
                 Use for overview zoom levels where individual bases aren't shown.
    """
    # Build mismatch index once for O(1) variant evidence lookups
    mismatch_index = _build_mismatch_index(data.reads)

    # Aggregate mismatches by position for LLM summary
    mismatch_counts: dict[tuple[int, str, str], int] = {}
    for read in data.reads:
        for mm in read.mismatches:
            key = (mm["pos"], mm["ref"], mm["alt"])
            mismatch_counts[key] = mismatch_counts.get(key, 0) + 1

    mismatch_summary = [
        {"pos": k[0], "ref": k[1], "alt": k[2], "count": v}
        for k, v in sorted(mismatch_counts.items())
    ]

    # Compute variant evidence using index (O(k) per variant instead of O(n*m) total)
    variant_evidence = {}
    enhanced_variants = []

    for variant in data.variants:
        key = f"{variant['position']}:{variant['ref']}>{variant['alt']}"
        evidence = _compute_variant_evidence_from_index(mismatch_index, variant)
        variant_evidence[key] = evidence

        # Enhance variant with evidence data for table display
        enhanced = dict(variant)
        enhanced["strand_forward"] = evidence["forward_count"]
        enhanced["strand_reverse"] = evidence["reverse_count"]
        enhanced["mean_quality"] = evidence["mean_quality"]

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

        enhanced["confidence"] = confidence
        # Keep is_low_confidence for backwards compatibility
        enhanced["is_low_confidence"] = confidence == "low"

        enhanced_variants.append(enhanced)

    # Serialize reads - compact mode omits large sequence/quality data
    if compact:
        reads_data = [
            {
                "name": r.name,
                "cigar": r.cigar,
                "position": r.position,
                "end_position": r.end_position,
                "mapping_quality": r.mapping_quality,
                "is_reverse": r.is_reverse,
                "mismatches": r.mismatches,
                "mate_position": r.mate_position,
                "mate_contig": r.mate_contig,
                "insert_size": r.insert_size,
                "is_proper_pair": r.is_proper_pair,
                "is_read1": r.is_read1,
            }
            for r in data.reads
        ]
    else:
        reads_data = [
            {
                "name": r.name,
                "sequence": r.sequence,
                "qualities": r.qualities,
                "cigar": r.cigar,
                "position": r.position,
                "end_position": r.end_position,
                "mapping_quality": r.mapping_quality,
                "is_reverse": r.is_reverse,
                "mismatches": r.mismatches,
                "mate_position": r.mate_position,
                "mate_contig": r.mate_contig,
                "insert_size": r.insert_size,
                "is_proper_pair": r.is_proper_pair,
                "is_read1": r.is_read1,
            }
            for r in data.reads
        ]

    return {
        "contig": data.contig,
        "start": data.start,
        "end": data.end,
        "reads": reads_data,
        "coverage": data.coverage,
        "variants": enhanced_variants,
        "reference_sequence": data.reference_sequence,
        "mismatch_summary": mismatch_summary,
        "variant_evidence": variant_evidence,
        "compact": compact,
    }

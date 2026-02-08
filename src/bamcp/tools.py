"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from typing import Any

import pysam

from .cache import BAMIndexCache
from .clinvar import ClinVarClient
from .config import BAMCPConfig
from .gnomad import GnomadClient
from .parsers import AlignedRead, RegionData, fetch_region

logger = logging.getLogger(__name__)

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


async def handle_browse_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle browse_region tool call.

    Returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload,
        },
    }


async def handle_get_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return variants without UI."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
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
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
    )

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
    reference = args.get("reference", config.reference)

    mode = "rc" if file_path.endswith(".cram") else "rb"
    samfile = pysam.AlignmentFile(
        file_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference,
        index_filename=_get_index_path(file_path, config),
    )

    contigs = [
        {"name": name, "length": length}
        for name, length in zip(samfile.references, samfile.lengths, strict=True)
    ]

    samfile.close()

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
    position = args["position"]
    contig = args.get("contig", "chr1")
    window = args.get("window", config.default_window)
    reference = args.get("reference", config.reference)

    start = max(0, position - window // 2)
    end = position + window // 2
    region = f"{contig}:{start}-{end}"

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
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
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload,
        },
    }


async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata.
    """
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        index_filename=_get_index_path(file_path, config),
    )

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

    client = ClinVarClient(api_key=config.ncbi_api_key)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except Exception:
        logger.exception("ClinVar lookup failed for %s:%d %s>%s", chrom, pos, ref, alt)
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

    client = GnomadClient(dataset=config.gnomad_dataset)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except Exception:
        logger.exception("gnomAD lookup failed for %s:%d %s>%s", chrom, pos, ref, alt)
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

    Args:
        values: List of integer values to bin.
        bins: List of bin edges (e.g., [0, 10, 20, 30, 40] creates bins 0-9, 10-19, etc.)

    Returns:
        Dict mapping bin labels to counts (e.g., {"0-9": 5, "10-19": 3, ...}).
    """
    histogram: dict[str, int] = {}

    for i in range(len(bins)):
        label = f"{bins[i]}+" if i == len(bins) - 1 else f"{bins[i]}-{bins[i + 1] - 1}"
        histogram[label] = 0

    for val in values:
        for i in range(len(bins)):
            if i == len(bins) - 1:
                # Last bin: all values >= bins[i]
                if val >= bins[i]:
                    label = f"{bins[i]}+"
                    histogram[label] = histogram.get(label, 0) + 1
                    break
            elif bins[i] <= val < bins[i + 1]:
                label = f"{bins[i]}-{bins[i + 1] - 1}"
                histogram[label] = histogram.get(label, 0) + 1
                break

    return histogram


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
        "quality_histogram": _bin_histogram(qualities, [0, 10, 20, 30, 40]),
        "position_histogram": _bin_histogram(read_positions, [0, 25, 50, 75, 100, 150]),
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

        # Low-confidence if ANY of these criteria fail:
        # - Depth < 10
        # - VAF < 10% (0.10)
        # - Mean base quality < 20
        # - Strand bias > 90% (0.9)
        is_low_confidence = (
            variant["depth"] < 10
            or variant["vaf"] < 0.10
            or evidence["mean_quality"] < 20
            or evidence["strand_bias"] > 0.9
        )
        enhanced["is_low_confidence"] = is_low_confidence

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

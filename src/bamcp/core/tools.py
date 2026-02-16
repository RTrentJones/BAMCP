"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any

import httpx
import pysam

from ..clients.clinvar import ClinVarClient
from ..clients.genes import GeneClient
from ..clients.gnomad import GnomadClient
from ..config import BAMCPConfig
from ..constants import (
    BAM_PARSE_TIMEOUT_SECONDS,
    DEFAULT_CONTIG,
    SCAN_VARIANTS_CHUNK_SIZE,
    SCAN_VARIANTS_MAX_REGION,
    SCAN_VARIANTS_TIMEOUT_SECONDS,
    VIEWER_RESOURCE_URI,
)
from .cache import BAMIndexCache
from .parsers import RegionData, fetch_region, scan_variants_chunked
from .serialization import serialize_region_data
from .validation import validate_path, validate_region, validate_variant_input

logger = logging.getLogger(__name__)

# Module-level singleton instances for session consistency and connection reuse
_cache_instance: BAMIndexCache | None = None
_clinvar_client: ClinVarClient | None = None
_gnomad_client: GnomadClient | None = None
_gene_client: GeneClient | None = None

_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)


def get_cache(config: BAMCPConfig) -> BAMIndexCache:
    """Get or create the session cache instance.

    Uses a module-level singleton to maintain consistent session ID
    across all tool calls within a server process.
    """
    global _cache_instance
    if _cache_instance is None:
        _cache_instance = BAMIndexCache(config.cache_dir, config.cache_ttl)
    return _cache_instance


def get_clinvar_client(config: BAMCPConfig) -> ClinVarClient:
    """Get or create the singleton ClinVar client."""
    global _clinvar_client
    if _clinvar_client is None:
        _clinvar_client = ClinVarClient(api_key=config.ncbi_api_key)
    return _clinvar_client


def get_gnomad_client(config: BAMCPConfig) -> GnomadClient:
    """Get or create the singleton gnomAD client."""
    global _gnomad_client
    if _gnomad_client is None:
        _gnomad_client = GnomadClient(dataset=config.gnomad_dataset)
    return _gnomad_client


def get_gene_client(config: BAMCPConfig) -> GeneClient:
    """Get or create the singleton gene client."""
    global _gene_client
    if _gene_client is None:
        _gene_client = GeneClient(
            api_key=config.ncbi_api_key,
            genome_build=config.genome_build,
        )
    return _gene_client


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


# -- Tool Handlers -----------------------------------------------------------


async def handle_get_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return variants without UI."""
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
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
    validate_region(region)
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
    """Handle jump_to tool call.

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
    payload = serialize_region_data(data)
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
    """Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config)
    payload = serialize_region_data(data)
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
    """Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
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


async def handle_lookup_clinvar(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in ClinVar via NCBI E-utilities.

    Returns clinical significance, review status, and associated conditions.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    # Validate input parameters
    validation_error = validate_variant_input(chrom, pos, ref, alt)
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

    client = get_clinvar_client(config)

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
    """Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    # Validate input parameters
    validation_error = validate_variant_input(chrom, pos, ref, alt)
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

    client = get_gnomad_client(config)

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


async def handle_scan_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Scan an entire contig for variants using fast coverage-based detection."""
    file_path = args["file_path"]
    validate_path(file_path, config)
    contig = args.get("contig", DEFAULT_CONTIG)
    reference = args.get("reference", config.reference)
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)

    if not reference:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "Reference genome required for variant scanning. "
                         "Use list_contigs to detect genome build and get a reference URL."}
                    ),
                }
            ]
        }

    index_path = await _ensure_cached_index(file_path, config)

    try:
        variants = await asyncio.wait_for(
            asyncio.to_thread(
                scan_variants_chunked,
                file_path,
                contig,
                reference,
                chunk_size=SCAN_VARIANTS_CHUNK_SIZE,
                min_vaf=min_vaf,
                min_depth=min_depth,
                max_region=SCAN_VARIANTS_MAX_REGION,
                index_filename=index_path,
            ),
            timeout=SCAN_VARIANTS_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps({"error": f"Scan timed out after {SCAN_VARIANTS_TIMEOUT_SECONDS}s"}),
                }
            ]
        }

    return {
        "content": [
            {
                "type": "text",
                "text": json.dumps(
                    {
                        "contig": contig,
                        "variants": variants,
                        "count": len(variants),
                    }
                ),
            }
        ]
    }

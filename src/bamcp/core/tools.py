"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import json
import logging
import tempfile
import time
from dataclasses import asdict, dataclass, field
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
    DEFAULT_CACHE_TTL_SECONDS,
    DEFAULT_CONTIG,
    INDEL_DETECTION_MAX_REGION,
    REGION_TILE_SIZE,
    SCAN_VARIANTS_CHUNK_SIZE,
    SCAN_VARIANTS_MAX_REGION,
    SCAN_VARIANTS_TIMEOUT_SECONDS,
    VIEWER_RESOURCE_URI,
)
from ..middleware.telemetry import telemetry_wrapper
from .cache import BAMIndexCache
from .parsers import (
    MAX_REGION_SIZE,
    RegionData,
    annotate_vcf_snv_support,
    fetch_candidate_variants_only,
    fetch_coverage_only,
    fetch_region,
    load_vcf_variants,
    parse_region,
    scan_variants_chunked,
)
from .serialization import serialize_region_data
from .validation import (
    validate_path,
    validate_region,
    validate_remote_url,
    validate_variant_file_path,
    validate_variant_input,
)

logger = logging.getLogger(__name__)


@dataclass
class BAMCPServices:
    """Per-server runtime services: the index cache, external API clients, the
    region cache, and per-index download locks.

    These were previously module-level singletons, which meant two servers built
    from different :class:`BAMCPConfig` instances in one process (tests, embedded
    deployments) would silently share a cache directory, NCBI API key, or gnomAD
    dataset — whichever config constructed them first. Bundling them here and
    keying by config identity (see :func:`get_services`) keeps each server
    isolated while preserving connection reuse within a server.
    """

    config: BAMCPConfig
    cache: BAMIndexCache
    region_cache: dict[tuple, tuple[float, RegionData]] = field(default_factory=dict)
    index_download_locks: dict[str, asyncio.Lock] = field(default_factory=dict)
    _clinvar: ClinVarClient | None = field(default=None, repr=False)
    _gnomad: GnomadClient | None = field(default=None, repr=False)
    _genes: GeneClient | None = field(default=None, repr=False)

    @classmethod
    def from_config(cls, config: BAMCPConfig) -> BAMCPServices:
        return cls(config=config, cache=BAMIndexCache(config.cache_dir, config.cache_ttl))

    def clinvar(self) -> ClinVarClient:
        if self._clinvar is None:
            self._clinvar = ClinVarClient(api_key=self.config.ncbi_api_key)
        return self._clinvar

    def gnomad(self) -> GnomadClient:
        if self._gnomad is None:
            self._gnomad = GnomadClient(dataset=self.config.gnomad_dataset)
        return self._gnomad

    def genes(self) -> GeneClient:
        if self._genes is None:
            self._genes = GeneClient(
                api_key=self.config.ncbi_api_key,
                genome_build=self.config.genome_build,
            )
        return self._genes

    async def aclose(self) -> None:
        """Close external HTTP clients and drop per-process caches.

        Keeps the index cache instance (and its session id) so callers can reuse
        the session; only the network clients and transient caches are released.
        """
        clients = (self._clinvar, self._gnomad, self._genes)
        await asyncio.gather(
            *(client.close() for client in clients if client is not None),
            return_exceptions=True,
        )
        self._clinvar = None
        self._gnomad = None
        self._genes = None
        self.region_cache.clear()
        self.index_download_locks.clear()


# Services keyed by config identity so distinct server instances stay isolated.
_services_registry: dict[int, BAMCPServices] = {}


def get_services(config: BAMCPConfig) -> BAMCPServices:
    """Return the :class:`BAMCPServices` for ``config``, creating it on first use.

    Keyed by ``id(config)``: every server has exactly one config object, so this
    memoizes services per server while never sharing them across servers.
    """
    services = _services_registry.get(id(config))
    if services is None:
        services = BAMCPServices.from_config(config)
        _services_registry[id(config)] = services
    return services


def get_cache(config: BAMCPConfig) -> BAMIndexCache:
    """Get or create the session index cache for this server's config."""
    return get_services(config).cache


def get_clinvar_client(config: BAMCPConfig) -> ClinVarClient:
    """Get or create this server's ClinVar client."""
    return get_services(config).clinvar()


def get_gnomad_client(config: BAMCPConfig) -> GnomadClient:
    """Get or create this server's gnomAD client."""
    return get_services(config).gnomad()


def get_gene_client(config: BAMCPConfig) -> GeneClient:
    """Get or create this server's gene client."""
    return get_services(config).genes()


async def close_external_clients(config: BAMCPConfig) -> None:
    """Close this server's external HTTP clients and clear its per-process caches."""
    services = _services_registry.get(id(config))
    if services is not None:
        await services.aclose()


_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)

_CANDIDATE_VARIANT_DISCLAIMER = (
    "Note: BAMCP reports research-grade candidate variants from read-level evidence; "
    "confirm with validated variant-calling and clinical workflows before use."
)

_INTERNAL_COORDINATE_SYSTEM = "0-based-half-open"
_LLM_COORDINATE_SYSTEM = "1-based-inclusive"

# Where get_variants sources its variants from (see handle_get_variants).
_VARIANT_SOURCES = ("auto", "vcf", "bamcp")


def _cache_key(
    file_path: str,
    region: str,
    reference: str | None,
    config: BAMCPConfig,
    mode: str,
    min_vaf: float | None,
    min_depth: int | None,
    vcf_path: str | None,
    vcf_primary: bool = False,
) -> tuple:
    return (
        mode,
        file_path,
        region,
        reference,
        vcf_path,
        vcf_primary,
        config.max_reads,
        config.min_mapq,
        config.min_baseq,
        min_vaf if min_vaf is not None else config.min_vaf,
        min_depth if min_depth is not None else config.min_depth,
    )


def _get_cached_region(
    region_cache: dict[tuple, tuple[float, RegionData]], key: tuple, ttl: int
) -> RegionData | None:
    if ttl <= 0:
        return None
    cached = region_cache.get(key)
    if cached is None:
        return None
    created, data = cached
    if (time.monotonic() - created) > ttl:
        region_cache.pop(key, None)
        return None
    return data


def _set_cached_region(
    region_cache: dict[tuple, tuple[float, RegionData]], key: tuple, data: RegionData, ttl: int
) -> RegionData:
    if ttl > 0:
        region_cache[key] = (time.monotonic(), data)
    return data


async def _download_index_streaming(
    client: httpx.AsyncClient,
    index_url: str,
    index_path: str,
    config: BAMCPConfig,
) -> str | None:
    """Download an index with SSRF revalidation, streaming, and a size cap."""
    validate_remote_url(index_url, config)
    total = 0
    tmp_path = ""
    async with client.stream("GET", index_url) as resp:
        if resp.status_code == 404:
            return None
        if resp.status_code != 200:
            logger.warning("Index download failed (%d): %s", resp.status_code, index_url)
            return None

        content_length = resp.headers.get("content-length")
        if content_length and int(content_length) > config.max_remote_index_bytes:
            raise ValueError("Remote index exceeds configured maximum size")

        Path(index_path).parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=str(Path(index_path).parent),
            prefix=".index-",
            suffix=".tmp",
            delete=False,
        ) as tmp:
            tmp_path = tmp.name
            async for chunk in resp.aiter_bytes():
                total += len(chunk)
                if total > config.max_remote_index_bytes:
                    raise ValueError("Remote index exceeds configured maximum size")
                tmp.write(chunk)

    Path(tmp_path).replace(index_path)
    logger.info("Cached index (%d bytes): %s", total, index_path)
    return index_path


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
    services = get_services(config)
    cache = services.cache
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

    lock = services.index_download_locks.setdefault(index_path, asyncio.Lock())
    async with lock:
        if cache.is_valid(index_path):
            logger.debug("Using cached index after waiting for download lock: %s", index_path)
            return index_path

        logger.info("Downloading index for remote BAM: %s", file_path)
        async with httpx.AsyncClient(timeout=60.0, follow_redirects=True) as client:
            for index_url in index_urls:
                try:
                    downloaded = await _download_index_streaming(
                        client, index_url, index_path, config
                    )
                    if downloaded:
                        return downloaded
                    logger.debug("Index not found at %s, trying next", index_url)
                except (httpx.RequestError, OSError, ValueError) as e:
                    logger.warning("Index download error for %s: %s", index_url, e)
                    continue

    # All attempts failed - let pysam try its own resolution
    logger.warning("Could not download index for %s, falling back to pysam", file_path)
    return None


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
    offset = start - data.start
    length = end - start
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
) -> RegionData:
    """Run the mode-appropriate parser for an exact region under the parse timeout."""
    effective_min_vaf = min_vaf if min_vaf is not None else config.min_vaf
    effective_min_depth = min_depth if min_depth is not None else config.min_depth

    if mode == "coverage":
        return await asyncio.wait_for(
            asyncio.to_thread(
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
            asyncio.to_thread(
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
        asyncio.to_thread(
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
    region_cache: dict[tuple, tuple[float, RegionData]],
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
    tile_region = f"{contig}:{t_start}-{t_end}"

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
    return any(
        v.get("source") == "vcf" and len(v["ref"]) == 1 == len(v["alt"]) for v in variants
    )


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
        asyncio.to_thread(
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
        validate_variant_file_path(vcf_path, config)

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

    if mode == "variants" and vcf_path and not reference:
        contig, start, end = parse_region(region)
        vcf_variants = await asyncio.wait_for(
            asyncio.to_thread(load_vcf_variants, vcf_path, region),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
        # Only touch the BAM when there's an SNV whose support we can actually count.
        if _has_vcf_snv(vcf_variants):
            index_path = await _ensure_cached_index(file_path, config)
            await _annotate_vcf_support(
                file_path, region, reference, vcf_variants, config, index_path
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
            ),
            ttl,
        )

    # Download and cache index for remote files (if not already cached)
    index_path = await _ensure_cached_index(file_path, config)
    data = await _compute_region_data(
        file_path, region, reference, config, mode, min_vaf, min_depth, index_path
    )

    if vcf_path:
        vcf_variants = await asyncio.wait_for(
            asyncio.to_thread(load_vcf_variants, vcf_path, region),
            timeout=BAM_PARSE_TIMEOUT_SECONDS,
        )
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
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
    reference = args.get("reference", config.reference)
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)

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
        vcf_path=effective_vcf,
        vcf_primary=vcf_primary,
    )

    variants = [
        _one_based(v)
        for v in data.variants
        if v.get("source") == "vcf" or (v["vaf"] >= min_vaf and v["depth"] >= min_depth)
    ]

    payload: dict[str, Any] = {
        "variants": variants,
        "count": len(variants),
        "variant_source": variant_source,
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
        "variant_type": "candidate",
        "disclaimer": _CANDIDATE_VARIANT_DISCLAIMER,
    }

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


@telemetry_wrapper("get_coverage")
async def handle_get_coverage(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return coverage statistics."""
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
    reference = args.get("reference", config.reference)

    data = await _fetch_region_with_timeout(file_path, region, reference, config, mode="coverage")

    coverage = data.coverage
    stats = {
        "region": f"{data.contig}:{data.start}-{data.end}",
        "mean": round(sum(coverage) / len(coverage), 2) if coverage else 0,
        "min": min(coverage) if coverage else 0,
        "max": max(coverage) if coverage else 0,
        "median": sorted(coverage)[len(coverage) // 2] if coverage else 0,
        "bases_covered": sum(1 for c in coverage if c > 0),
        "total_bases": len(coverage),
        "coordinate_system": _INTERNAL_COORDINATE_SYSTEM,
    }

    return {"content": [{"type": "text", "text": json.dumps(stats)}]}


@telemetry_wrapper("list_contigs")
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
                        "coordinate_system": _INTERNAL_COORDINATE_SYSTEM,
                    }
                ),
            }
        ]
    }


@telemetry_wrapper("jump_to")
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
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))

    start = max(0, position - window // 2)
    end = position + window // 2
    region = f"{contig}:{start}-{end}"

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, vcf_path=effective_vcf, vcf_primary=vcf_primary
    )
    payload = serialize_region_data(data)
    payload["coordinate_system"] = _INTERNAL_COORDINATE_SYSTEM
    payload["variant_type"] = "candidate"
    payload["disclaimer"] = _CANDIDATE_VARIANT_DISCLAIMER
    payload["file_path"] = file_path  # For client-side re-queries
    # Preserve the variant sourcing so the viewer's pan/zoom refetches keep it.
    payload["variant_source"] = variant_source
    if effective_vcf is not None:
        payload["vcf_path"] = effective_vcf

    # Return summary text in content (for LLM context), full data only in _meta
    reads_count = len(data.reads)
    variants_count = len(data.variants)
    summary = (
        f"Jumped to {data.contig}:{position}: {reads_count} reads, "
        f"{variants_count} candidate variants"
    )
    return {
        "content": [{"type": "text", "text": summary}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


@telemetry_wrapper("visualize_region")
async def handle_visualize_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
    reference = args.get("reference", config.reference)
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, vcf_path=effective_vcf, vcf_primary=vcf_primary
    )
    payload = serialize_region_data(data)
    payload["coordinate_system"] = _INTERNAL_COORDINATE_SYSTEM
    payload["variant_type"] = "candidate"
    payload["disclaimer"] = _CANDIDATE_VARIANT_DISCLAIMER
    payload["file_path"] = file_path  # For client-side re-queries
    # Preserve the variant sourcing so the viewer's pan/zoom refetches keep it.
    payload["variant_source"] = variant_source
    if effective_vcf is not None:
        payload["vcf_path"] = effective_vcf

    # Return summary text in content (for LLM context), full data only in _meta
    reads_count = len(data.reads)
    variants_count = len(data.variants)
    summary = (
        f"Region {data.contig}:{data.start}-{data.end}: "
        f"{reads_count} reads, {variants_count} candidate variants"
    )
    return {
        "content": [{"type": "text", "text": summary}],
        "_meta": {
            "ui/resourceUri": VIEWER_RESOURCE_URI,
            "ui/init": payload,
        },
    }


@telemetry_wrapper("get_region_summary")
async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata.
    """
    file_path = args["file_path"]
    validate_path(file_path, config)
    region = args["region"]
    validate_region(region)
    reference = args.get("reference", config.reference)
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, vcf_path=effective_vcf, vcf_primary=vcf_primary
    )

    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 2) if coverage else 0
    max_cov = max(coverage) if coverage else 0

    summary_lines = [
        f"Region: {data.contig}:{data.start}-{data.end}",
        f"Reads: {len(data.reads)}",
        f"Coverage: mean={mean_cov}x, max={max_cov}x",
        f"Candidate variants detected: {len(data.variants)}",
        f"Coordinate system: {_LLM_COORDINATE_SYSTEM}",
        _CANDIDATE_VARIANT_DISCLAIMER,
    ]

    for v in data.variants:
        # 1-based position (VCF/dbSNP convention) — see _one_based().
        summary_lines.append(
            f"  {v['contig']}:{v['position'] + 1} {v['ref']}>{v['alt']} "
            f"VAF={v['vaf']:.1%} depth={v['depth']}"
        )

    return {"content": [{"type": "text", "text": "\n".join(summary_lines)}]}


@telemetry_wrapper("lookup_clinvar")
async def handle_lookup_clinvar(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in ClinVar via NCBI E-utilities.

    Returns clinical significance, review status, and associated conditions.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    if not config.clinvar_enabled:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "ClinVar lookup is disabled", "disclaimer": _CLINVAR_DISCLAIMER}
                    ),
                }
            ]
        }

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


@telemetry_wrapper("lookup_gnomad")
async def handle_lookup_gnomad(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    if not config.gnomad_enabled:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "gnomAD lookup is disabled", "disclaimer": _GNOMAD_DISCLAIMER}
                    ),
                }
            ]
        }

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


@telemetry_wrapper("scan_variants")
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
                        {
                            "error": "Reference genome required for variant scanning. "
                            "Use list_contigs to detect genome build and get a reference URL."
                        }
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
                min_mapq=config.min_mapq,
                min_baseq=config.min_baseq,
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
                    "text": json.dumps(
                        {"error": f"Scan timed out after {SCAN_VARIANTS_TIMEOUT_SECONDS}s"}
                    ),
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
                        "variants": [_one_based(v) for v in variants],
                        "count": len(variants),
                        "coordinate_system": _LLM_COORDINATE_SYSTEM,
                        "variant_type": "candidate",
                        "disclaimer": _CANDIDATE_VARIANT_DISCLAIMER,
                    }
                ),
            }
        ]
    }

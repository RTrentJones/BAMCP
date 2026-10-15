"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import asyncio
import functools
import json
import logging
import os
import re
import tempfile
import threading
import time
import weakref
from collections import OrderedDict
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, TypeVar
from urllib.parse import urljoin

import httpx
import pysam

from ..analysis.evidence import (
    _not_computed_evidence,
    _vcf_evidence,
    annotate_interpretation,
    compute_artifact_risk,
    compute_confidence,
    enhance_vcf_variants,
)
from ..clients.clinvar import ClinVarClient
from ..clients.genes import GeneClient, GeneLookupError
from ..clients.gnomad import GnomadClient, GnomadQueryError
from ..clients.ncbi import NCBIUnavailableError, TokenBucketRateLimiter, ncbi_rate_limiter
from ..config import BAMCPConfig
from ..constants import (
    BAM_PARSE_TIMEOUT_SECONDS,
    COLD_START_TIMEOUT_SECONDS,
    DEFAULT_CACHE_TTL_SECONDS,
    INDEL_DETECTION_MAX_REGION,
    PARSE_MAX_CONCURRENCY,
    REGION_CACHE_MAX_BYTES,
    REGION_CACHE_MAX_ENTRIES,
    REGION_TILE_SIZE,
    SCAN_VARIANTS_CHUNK_SIZE,
    SCAN_VARIANTS_CHUNK_TARGET_SECONDS,
    SCAN_VARIANTS_DEADLINE_SECONDS,
    SCAN_VARIANTS_MAX_REGION,
    SCAN_VARIANTS_MAX_RESULTS,
    SCAN_VARIANTS_MIN_CHUNK_SIZE,
    SCAN_VARIANTS_TIMEOUT_SECONDS,
    VIEWER_RESOURCE_URI,
    WARM_INPUTS_MAX_ENTRIES,
)
from ..middleware.telemetry import telemetry_wrapper
from .cache import BAMIndexCache
from .errors import ToolInputError
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
    parse_region,
)
from .serialization import build_enhanced_variants, serialize_region_data
from .validation import (
    is_cram,
    validate_path,
    validate_reference_path,
    validate_region,
    validate_remote_url,
    validate_variant_file_path,
    validate_variant_input,
)

logger = logging.getLogger(__name__)

_T = TypeVar("_T")

# All blocking pysam parses run in this dedicated, bounded pool rather than the event loop's
# default executor, so a hung parse (wait_for abandons the await but can't cancel the C call)
# can't starve unrelated async work. A truly cancellable parse would need a killable worker
# *process*; that is a deliberate follow-up.
_PARSE_EXECUTOR = ThreadPoolExecutor(
    max_workers=PARSE_MAX_CONCURRENCY, thread_name_prefix="bamcp-parse"
)
# Submission gate: a plain ThreadPoolExecutor has an UNBOUNDED work queue, so gating on the pool
# alone isn't real backpressure — a burst of parses (or hung workers) would let submissions pile
# up in that queue and grow memory. This semaphore bounds in-flight submissions to the worker
# count, so the executor queue never grows; excess requests wait here (bounded by the async
# request count) instead. A permit is held until the underlying thread truly finishes (released
# via a done-callback, and the future is shielded from await-cancellation), so a hung parse keeps
# occupying its slot rather than freeing it for a submission that would then queue behind it.
#
# asyncio primitives bind to the first loop that contends on them, so a single module-level
# Semaphore breaks ("bound to a different event loop") once a second loop uses it — repeated
# asyncio.run() in embedded/sync callers, or test suites. Keep one gate per running loop. In a
# server there is exactly one loop, so the bound is unchanged.
_parse_gates: weakref.WeakKeyDictionary[asyncio.AbstractEventLoop, asyncio.Semaphore] = (
    weakref.WeakKeyDictionary()
)


def _get_parse_gate() -> asyncio.Semaphore:
    """The submission gate for the running event loop (created on first use)."""
    loop = asyncio.get_running_loop()
    gate = _parse_gates.get(loop)
    if gate is None:
        gate = _parse_gates[loop] = asyncio.Semaphore(PARSE_MAX_CONCURRENCY)
    return gate


async def _run_in_parse_pool(fn: Callable[..., _T], /, *args: Any, **kwargs: Any) -> _T:
    """Run blocking parse work in the dedicated pool under the submission gate.

    Drop-in awaitable for ``asyncio.to_thread`` (callers wrap it in their own ``wait_for``).
    run_in_executor takes no kwargs, so bind them with partial to keep the call sites unchanged.
    """
    call = functools.partial(fn, **kwargs) if kwargs else fn
    gate = _get_parse_gate()
    await gate.acquire()
    fut = asyncio.get_running_loop().run_in_executor(_PARSE_EXECUTOR, call, *args)
    # Release only when the thread actually completes — not when an outer wait_for cancels the
    # await — so hung parses keep holding their slot (real backpressure, bounded executor queue).
    fut.add_done_callback(lambda _f: gate.release())
    return await asyncio.shield(fut)


_RegionCacheEntry = tuple[float, RegionData, int]


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
    # LRU-ordered (newest at the end); bounded by entry count AND approximate bytes in
    # _set_cached_region. Values are (created_monotonic, data, approx_bytes).
    region_cache: OrderedDict[tuple, _RegionCacheEntry] = field(default_factory=OrderedDict)
    index_download_locks: dict[str, asyncio.Lock] = field(default_factory=dict)
    # Parsed BAM/CRAM headers (contig name/length lists), LRU-bounded by entry count
    # (_HEADER_CACHE_MAX). Keyed by (file_path, local fingerprint) so a replaced local file
    # is re-read. Feeds build/style detection for auto-reference and search_gene.
    header_cache: OrderedDict[tuple, list[dict]] = field(default_factory=OrderedDict)
    _clinvar: ClinVarClient | None = field(default=None, repr=False)
    _gnomad: GnomadClient | None = field(default=None, repr=False)
    _genes: GeneClient | None = field(default=None, repr=False)
    # ONE token bucket for every NCBI E-utilities caller (ClinVar + Gene): NCBI's
    # 3 req/s (10 with a key) limit is per caller, not per database.
    _ncbi_limiter: TokenBucketRateLimiter | None = field(default=None, repr=False)

    @classmethod
    def from_config(cls, config: BAMCPConfig) -> BAMCPServices:
        cache = BAMIndexCache(
            config.cache_dir, config.cache_ttl, max_total_bytes=config.max_index_cache_bytes
        )
        # Sweep index dirs left behind by expired sessions (e.g. a crashed process whose
        # atexit cleanup never ran) so the shared cache dir can't accumulate forever.
        cache.cleanup_expired_sessions()
        return cls(config=config, cache=cache)

    def ncbi_rate_limiter(self) -> TokenBucketRateLimiter:
        if self._ncbi_limiter is None:
            self._ncbi_limiter = ncbi_rate_limiter(self.config.ncbi_api_key)
        return self._ncbi_limiter

    def clinvar(self) -> ClinVarClient:
        if self._clinvar is None:
            self._clinvar = ClinVarClient(
                api_key=self.config.ncbi_api_key,
                genome_build=self.config.genome_build,
                rate_limiter=self.ncbi_rate_limiter(),
            )
        return self._clinvar

    def gnomad(self) -> GnomadClient:
        if self._gnomad is None:
            self._gnomad = GnomadClient(
                dataset=self.config.gnomad_dataset, genome_build=self.config.genome_build
            )
        return self._gnomad

    def genes(self) -> GeneClient:
        if self._genes is None:
            self._genes = GeneClient(
                api_key=self.config.ncbi_api_key,
                genome_build=self.config.genome_build,
                rate_limiter=self.ncbi_rate_limiter(),
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
        self._ncbi_limiter = None
        self.region_cache.clear()
        self.index_download_locks.clear()
        self.header_cache.clear()

    def clear_caches(self) -> int:
        """Drop downloaded index files and the region cache; return index files removed.

        Deliberately leaves the shared HTTP clients open: closing them here would make any
        in-flight ClinVar/gnomAD/gene lookup fail with ``RuntimeError`` (closed client).
        """
        removed = self.cache.cleanup_session()
        self.region_cache.clear()
        self.header_cache.clear()
        return removed


# Services keyed by config identity so distinct server instances stay isolated. LRU-ordered and
# bounded so it cannot grow without limit across many short-lived configs (tests, embedded use);
# in the common case a server has exactly one config, so this holds a single entry.
_services_registry: OrderedDict[int, BAMCPServices] = OrderedDict()
_SERVICES_REGISTRY_MAX = 64
# Strong refs to in-flight eviction-close tasks so they aren't GC'd before completing.
_pending_service_closes: set[asyncio.Task] = set()


def _schedule_services_close(services: BAMCPServices) -> None:
    """Close an evicted/replaced services' HTTP clients so they don't leak.

    ``get_services`` is synchronous but ``aclose`` is async. When called from within a running
    loop (the normal request path) schedule it as a task; when called from sync embedded code
    with no loop, close it synchronously with ``asyncio.run`` — a sync caller can still have
    lazily created ``httpx.AsyncClient``s (they don't require a loop to construct), so we must
    not silently skip the close there.
    """
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        try:
            asyncio.run(services.aclose())
        except Exception:  # noqa: BLE001 — cleanup must never break eviction
            logger.debug("failed to close evicted services synchronously", exc_info=True)
        return
    task = loop.create_task(services.aclose())
    _pending_service_closes.add(task)
    task.add_done_callback(_pending_service_closes.discard)


def get_services(config: BAMCPConfig) -> BAMCPServices:
    """Return the :class:`BAMCPServices` for ``config``, creating it on first use.

    Keyed by ``id(config)`` with an **identity guard**: a cached entry is reused only if it
    still belongs to *this* config object, so a reused ``id()`` (after an old config was
    collected — its strong ref here having been dropped by :func:`close_external_clients` or
    LRU eviction) can never alias the wrong services. The registry is a bounded LRU; an evicted or
    stale-overwritten services has its HTTP clients closed (scheduled on the loop) so it does not
    leak connections. Call :func:`close_external_clients` for prompt, clean teardown.
    """
    key = id(config)
    services = _services_registry.get(key)
    if services is not None and services.config is config:
        _services_registry.move_to_end(key)
        return services

    # Stale entry for a reused id() (belonged to a now-dead config) — close it before replacing.
    if services is not None:
        _schedule_services_close(services)

    services = BAMCPServices.from_config(config)
    _services_registry[key] = services
    _services_registry.move_to_end(key)
    while len(_services_registry) > _SERVICES_REGISTRY_MAX:
        _, evicted = _services_registry.popitem(last=False)
        _schedule_services_close(evicted)
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
    """Close this server's external HTTP clients, clear caches, and drop the registry entry."""
    key = id(config)
    services = _services_registry.get(key)
    if services is not None and services.config is config:
        await services.aclose()
        # Explicit lifecycle removal so a closed server leaves nothing behind.
        _services_registry.pop(key, None)


_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)

# Lookup outcome vocabulary for the external-database tools. The load-bearing distinction
# is not_found (genuine absence — the variant is truly not in the database) vs. unavailable
# (the lookup errored/timed out — absence is UNKNOWN and must NOT be read as "rare/benign").
# Collapsing the two is the circular trap that turns an empty gnomAD result into false
# rarity (PM2) evidence; the explicit status keeps callers honest.
_LOOKUP_FOUND = "found"
_LOOKUP_NOT_FOUND = "not_found"
_LOOKUP_UNAVAILABLE = "unavailable"
_LOOKUP_DISABLED = "disabled"


def _lookup_result(payload: dict[str, Any], *, status: str, disclaimer: str) -> dict:
    """Wrap an external-database lookup payload with an explicit ``status`` field so a
    caller can tell a genuine not-found from a failed lookup (see the vocabulary above)."""
    payload["status"] = status
    payload.setdefault("disclaimer", disclaimer)
    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


def _lookup_build(args: dict[str, Any], config: BAMCPConfig) -> tuple[str | None, str]:
    """Resolve a lookup's genome build: explicit ``build`` arg, else the server default.

    Returns ``(canonical_build, source)``; ``canonical_build`` is None for an unrecognized
    explicit build (the caller reports invalid input rather than querying the wrong assembly).
    """
    from ..constants import normalize_genome_build

    if args.get("build"):
        return normalize_genome_build(args["build"]), "explicit"
    return normalize_genome_build(config.genome_build), "server-default"


def _validate_lookup_input(args: dict[str, Any], config: BAMCPConfig) -> tuple[str, str]:
    """Validate a variant lookup's chrom/pos/ref/alt/build; raise ToolInputError if invalid.

    Returns the resolved ``(build, build_source)``.
    """
    error = validate_variant_input(args["chrom"], args["pos"], args["ref"], args["alt"])
    if error:
        raise ToolInputError(error)
    build, build_source = _lookup_build(args, config)
    if build is None:
        raise ToolInputError(
            f"Unrecognized build {args.get('build')!r}; use GRCh37/hg19 or GRCh38/hg38"
        )
    return build, build_source


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
    if "excluded_reads" in filters:
        parts.append(f"excluded reads: {'/'.join(filters['excluded_reads'])}")
    snv = filters.get("snv_detection")
    if snv and snv != "enabled":
        parts.append(f"SNV detection {snv}")
    if filters.get("indel_detection") is False:
        parts.append(f"indel detection {filters.get('indel_detection_note', 'off')}")
    return ", ".join(parts)


def _local_fingerprint(path: str | None) -> tuple[int, int] | None:
    """``(mtime_ns, size)`` of a local file, so a replaced file never serves stale cache hits.

    Remote URLs (and unreadable paths) return None — remote objects carry no cheap identity
    here, and an unreadable path fails in the parser anyway.
    """
    if not path or "://" in path:
        return None
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_mtime_ns, st.st_size)


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
        _local_fingerprint(file_path),
        region,
        reference,
        _local_fingerprint(reference),
        vcf_path,
        _local_fingerprint(vcf_path),
        vcf_primary,
        config.max_reads,
        config.min_mapq,
        config.min_baseq,
        min_vaf if min_vaf is not None else config.min_vaf,
        min_depth if min_depth is not None else config.min_depth,
    )


def _approx_region_bytes(data: Any) -> int:
    """Rough in-memory size of a cached RegionData (CPython object overheads included).

    Entry *count* alone does not bound memory: a full-mode 1 Mb region with max_reads reads can
    weigh tens of MB, so 128 of them would blow far past a 2 GB container. Estimates per read
    ~600 B of fixed overhead plus ~9 B/base (sequence char + quality-list slot), ~16 B per
    coverage position, ~1 KB per variant dict, plus the reference string.
    """
    total = 1024
    for read in getattr(data, "reads", None) or ():
        total += 600 + 9 * len(getattr(read, "sequence", "") or "")
        total += 200 * len(getattr(read, "mismatches", ()) or ())
    total += 16 * len(getattr(data, "coverage", None) or ())
    total += 1024 * len(getattr(data, "variants", None) or ())
    total += len(getattr(data, "reference_sequence", None) or "")
    return total


def _get_cached_region(
    region_cache: OrderedDict[tuple, _RegionCacheEntry], key: tuple, ttl: int
) -> RegionData | None:
    if ttl <= 0:
        return None
    cached = region_cache.get(key)
    if cached is None:
        return None
    created, data, _ = cached
    if (time.monotonic() - created) > ttl:
        region_cache.pop(key, None)
        return None
    region_cache.move_to_end(key)  # mark most-recently-used
    return data


def _set_cached_region(
    region_cache: OrderedDict[tuple, _RegionCacheEntry],
    key: tuple,
    data: RegionData,
    ttl: int,
) -> RegionData:
    if ttl <= 0:
        return data
    size = _approx_region_bytes(data)
    # One entry may take at most a quarter of the budget; bigger ones are served uncached
    # rather than evicting most of the cache for a single payload.
    if size > REGION_CACHE_MAX_BYTES // 4:
        region_cache.pop(key, None)
        return data
    region_cache[key] = (time.monotonic(), data, size)
    region_cache.move_to_end(key)
    # Bound memory: evict least-recently-used entries beyond the count cap or byte budget.
    total = sum(entry[2] for entry in region_cache.values())
    while region_cache and (
        len(region_cache) > REGION_CACHE_MAX_ENTRIES or total > REGION_CACHE_MAX_BYTES
    ):
        _, evicted = region_cache.popitem(last=False)
        total -= evicted[2]
    return data


# Redirects for index downloads are followed by hand so every hop is re-validated
# (httpx's follow_redirects would only validate the first URL).
_INDEX_MAX_REDIRECTS = 5
_REDIRECT_STATUSES = frozenset({301, 302, 303, 307, 308})


async def _download_index_streaming(
    client: httpx.AsyncClient,
    index_url: str,
    index_path: str,
    config: BAMCPConfig,
) -> str | None:
    """Download an index with SSRF revalidation, streaming, and a size cap.

    Redirects are followed manually (at most ``_INDEX_MAX_REDIRECTS`` hops) and
    every ``Location`` is re-run through :func:`validate_remote_url`, so an allowed
    host cannot bounce the download to an internal address. The DNS-resolving
    validation runs in a worker thread to keep the event loop responsive. A
    partially written temp file is always removed on failure.
    """
    url = index_url
    for _hop in range(_INDEX_MAX_REDIRECTS + 1):
        await asyncio.to_thread(validate_remote_url, url, config)
        async with client.stream("GET", url) as resp:
            if resp.status_code in _REDIRECT_STATUSES:
                location = resp.headers.get("location")
                if not location:
                    logger.warning("Index redirect without Location: %s", url)
                    return None
                url = urljoin(url, location)
                continue
            if resp.status_code == 404:
                return None
            if resp.status_code != 200:
                logger.warning("Index download failed (%d): %s", resp.status_code, url)
                return None
            return await _stream_index_to_file(resp, index_path, config)
    raise ValueError(f"Index download exceeded {_INDEX_MAX_REDIRECTS} redirects")


async def _stream_index_to_file(resp: httpx.Response, index_path: str, config: BAMCPConfig) -> str:
    """Stream ``resp`` into ``index_path`` atomically, enforcing the size cap."""
    content_length = resp.headers.get("content-length")
    if content_length and int(content_length) > config.max_remote_index_bytes:
        raise ValueError("Remote index exceeds configured maximum size")

    Path(index_path).parent.mkdir(parents=True, exist_ok=True)
    total = 0
    tmp_path = ""
    try:
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
    except BaseException:
        # Oversized body, network error mid-stream, cancellation, failed rename:
        # never leave a ``.index-XXXX.tmp`` behind in the cache directory.
        if tmp_path:
            Path(tmp_path).unlink(missing_ok=True)
        raise

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
    if is_cram(file_path):
        index_urls = [file_path + ".crai"]
    else:
        # BAM files: try .bam.bai first, then .bai
        index_urls = [file_path + ".bai", file_path.rsplit(".", 1)[0] + ".bai"]

    locks = services.index_download_locks
    lock = locks.setdefault(index_path, asyncio.Lock())
    try:
        async with lock:
            if cache.is_valid(index_path):
                logger.debug("Using cached index after waiting for download lock: %s", index_path)
                return index_path

            logger.info("Downloading index for remote BAM: %s", file_path)
            async with httpx.AsyncClient(timeout=60.0, follow_redirects=False) as client:
                for index_url in index_urls:
                    try:
                        downloaded = await _download_index_streaming(
                            client, index_url, index_path, config
                        )
                        if downloaded:
                            cache.enforce_budget(keep=downloaded)
                            return downloaded
                        logger.debug("Index not found at %s, trying next", index_url)
                    except (httpx.RequestError, OSError, ValueError) as e:
                        logger.warning("Index download error for %s: %s", index_url, e)
                        continue
    finally:
        # Drop the per-URL lock once nobody holds it so the dict can't grow one entry per
        # URL forever. A waiter already woken keeps its reference; at worst a racing new
        # caller creates a fresh lock and repeats an (atomic, idempotent) download.
        if not lock.locked() and locks.get(index_path) is lock:
            del locks[index_path]

    # All attempts failed - let pysam try its own resolution
    logger.warning("Could not download index for %s, falling back to pysam", file_path)
    return None


# Remote inputs (BAM/CRAM URLs, reference FASTA URLs) that have been opened successfully at
# least once in this process. The first open of a remote input fetches its header (and, for a
# FASTA, the .fai/.gzi) and is the slow, cold part of a call; it runs under the larger
# COLD_START_TIMEOUT_SECONDS budget so it can't eat the normal parse budget. Bounded LRU.
_warm_inputs: OrderedDict[str, None] = OrderedDict()


def _is_remote(path: str | None) -> bool:
    return bool(path) and "://" in str(path)


def _mark_warm(url: str) -> None:
    _warm_inputs[url] = None
    _warm_inputs.move_to_end(url)
    while len(_warm_inputs) > WARM_INPUTS_MAX_ENTRIES:
        _warm_inputs.popitem(last=False)


def _open_remote_inputs_sync(
    file_path: str | None, reference: str | None, index_path: str | None
) -> None:
    """Open (and immediately close) remote inputs so their headers/indexes are fetched.

    Runs in the parse pool. ``file_path``/``reference`` are only passed when cold. Touching
    the FASTA's first contig forces htslib to load the .fai (and .gzi for bgzip).
    """
    from .parsers import _open_alignment

    if reference:
        with pysam.FastaFile(reference) as fasta:
            if fasta.references:
                fasta.fetch(fasta.references[0], 0, 1)
    if file_path:
        with _open_alignment(file_path, reference, index_path) as samfile:
            _ = samfile.references


async def _warm_remote_inputs(
    file_path: str, reference: str | None, index_path: str | None
) -> None:
    """Give a cold remote BAM/CRAM or reference its one-off open under the cold-start budget.

    No-op for local inputs and for inputs already warm, so only the first call to a not-yet-seen
    URL pays (and is allowed) the cold cost. A cold-open timeout raises a diagnosable
    ``TimeoutError``; any other failure is left for the real parse to report with its own
    message (the input simply stays cold).
    """
    cold_file = file_path if _is_remote(file_path) and file_path not in _warm_inputs else None
    cold_ref = reference if _is_remote(reference) and reference not in _warm_inputs else None
    if not cold_file and not cold_ref:
        return
    try:
        await asyncio.wait_for(
            _run_in_parse_pool(_open_remote_inputs_sync, cold_file, cold_ref, index_path),
            timeout=COLD_START_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError as e:
        what = " and ".join(
            name
            for name, cold in (("alignment file", cold_file), ("reference FASTA", cold_ref))
            if cold
        )
        raise TimeoutError(
            f"Timed out after {COLD_START_TIMEOUT_SECONDS:g}s opening the remote {what} on "
            "first use (fetching its header/index); retrying usually succeeds once it is warm."
        ) from e
    except Exception as e:  # noqa: BLE001 — the real parse reports the actual error
        logger.debug("Cold open of remote inputs failed (left cold): %s", e)
        return
    for url in (cold_file, cold_ref):
        if url:
            _mark_warm(url)


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
    # Cold start: fetch the index (never under a parse budget), then give a not-yet-seen remote
    # BAM/reference its one-off open under the separate cold budget, so the parse that follows
    # runs warm under the normal BAM_PARSE_TIMEOUT_SECONDS.
    if _is_remote(file_path) or _is_remote(reference):
        index_path = await _ensure_cached_index(file_path, config)
        await _warm_remote_inputs(file_path, reference, index_path)
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
        detail = (
            "; the remote reference FASTA's index is fetched on first use and is often "
            "the slow part — retrying usually succeeds once it is warm"
            if reference and "://" in reference
            else ""
        )
        raise TimeoutError(
            f"Timed out after {BAM_PARSE_TIMEOUT_SECONDS:g}s reading {region} "
            f"(mode={mode}){detail}."
        ) from e


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
    await _warm_remote_inputs(file_path, None, index_path)
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


def _passes_candidate_filters(
    v: dict[str, Any], min_vaf: float, min_depth: int, min_alt_count: int
) -> bool:
    """Tool-layer thresholds for a BAMCP candidate. VCF records always pass (authoritative)."""
    if v.get("source") == "vcf":
        return True
    alt_count = v.get("alt_count")
    return (
        v["vaf"] >= min_vaf
        and v["depth"] >= min_depth
        and (alt_count is None or alt_count >= min_alt_count)
    )


def _apply_alt_count(data: RegionData, min_alt_count: int) -> RegionData:
    """``data`` with BAMCP candidates below ``min_alt_count`` dropped (never mutates the cache).

    Detection's VAF/depth floors already ran in the parser; the alt-read floor is applied here so
    the region cache stays keyed by the parser thresholds only.
    """
    kept = [
        v
        for v in data.variants
        if v.get("source") == "vcf" or (v.get("alt_count") or 0) >= min_alt_count
    ]
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
    summary = {
        "scored_sites": len(scored_ids),
        "max_sites": GET_VARIANTS_EVIDENCE_MAX_SITES,
        "capped": len(snvs) > GET_VARIANTS_EVIDENCE_MAX_SITES,
        "timed_out": timed_out,
    }
    return out, summary


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
    """
    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    region = args["region"]
    validate_region(region)
    resolved = await _resolve_reference_info(args, config)
    reference = resolved.effective
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)
    min_alt_count = args.get("min_alt_count", config.min_alt_count)

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

    # VCF SNVs get the viewer's confidence/artifact contract from their annotated support;
    # BAMCP candidates are then scored per site (bounded) by _score_candidates. Copies only:
    # the fetched variant dicts live in the region cache.
    detected = enhance_vcf_variants(data.variants, data.reference_sequence, data.start)
    kept = [
        dict(v) for v in detected if _passes_candidate_filters(v, min_vaf, min_depth, min_alt_count)
    ]
    scored, evidence_summary = await _score_candidates(
        file_path, region, reference, kept, data, config
    )
    variants = [annotate_interpretation(_one_based(v)) for v in scored]

    filters = _applied_filters(
        config,
        args,
        include_max_reads=False,
        detection=_detection_status(region, reference, vcf_primary),
        min_alt_count=min_alt_count,
        reference_source=resolved.source,
    )
    filters["per_site_evidence"] = evidence_summary
    if resolved.note:
        filters["reference_note"] = resolved.note
    payload: dict[str, Any] = {
        "variants": variants,
        "count": len(variants),
        "variant_source": variant_source,
        "coordinate_system": _LLM_COORDINATE_SYSTEM,
        "variant_type": "candidate",
        "applied_filters": filters,
        "filter_note": _FILTER_ABSENCE_NOTE,
        "disclaimer": _CANDIDATE_VARIANT_DISCLAIMER,
    }

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


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
    await _warm_remote_inputs(file_path, reference, index_path)
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

    def _err(msg: str, **extra: Any) -> dict:
        return {"content": [{"type": "text", "text": json.dumps({"error": msg, **extra})}]}

    symbol = args.get("symbol", "")
    if not symbol:
        raise ToolInputError("search_gene requires 'symbol'")

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
                return _err(
                    f"could not read the BAM header to detect its build ({type(e).__name__}); "
                    "pass 'build' explicitly (confirm it with list_contigs)",
                    status=_LOOKUP_UNAVAILABLE,
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
        return _err(
            "Gene lookup unavailable (NCBI error); retry later",
            status=_LOOKUP_UNAVAILABLE,
            reason=str(e),
            genome_build=effective,
        )
    if result is None:
        return _err(f"Gene {symbol!r} not found", status=_LOOKUP_NOT_FOUND, genome_build=effective)

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
    await _warm_remote_inputs(file_path, reference, index_path)
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
    file_path: str,
    region: str,
    resolved: ResolvedReference,
    summary_prefix: Callable[[RegionData], str],
) -> dict:
    """Fetch a region and build the MCP Apps viewer response (shared by visualize/jump_to).

    The text content (what the LLM reads) is a compact summary that echoes the applied
    filters; the full payload goes only in ``_meta.ui/init``.
    """
    reference = resolved.effective
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))

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
    payload["variant_source"] = variant_source
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
    file_path = args["file_path"]
    contig = args.get("contig")
    if not contig:
        raise ToolInputError(
            "jump_to needs a contig as well as a position (e.g. contig='chr17' or '17' — "
            "call list_contigs for this file's naming)."
        )
    await asyncio.to_thread(validate_path, file_path, config)
    position = args["position"]
    window = args.get("window", config.default_window)
    resolved = await _resolve_reference_info(args, config)

    # ``position`` is 1-based; build a 1-based inclusive window of ``window`` bases
    # centred on it (start clamped at base 1).
    start = max(1, position - window // 2)
    end = start + max(window, 1) - 1
    region = f"{contig}:{start}-{end}"
    validate_region(region)

    return await _viewer_response(
        args,
        config,
        file_path,
        region,
        resolved,
        lambda data: f"Jumped to {data.contig}:{position}",
    )


@telemetry_wrapper("visualize_region")
async def handle_visualize_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    region = args["region"]
    validate_region(region)
    resolved = await _resolve_reference_info(args, config)

    return await _viewer_response(
        args,
        config,
        file_path,
        region,
        resolved,
        lambda data: f"Region {format_region(data.contig, data.start, data.end)}",
    )


@telemetry_wrapper("get_region_summary")
async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata. Each variant line is followed by
    its plain-language interpretation (zygosity hint, artifact risk, confidence).
    """
    file_path = args["file_path"]
    await asyncio.to_thread(validate_path, file_path, config)
    region = args["region"]
    validate_region(region)
    resolved = await _resolve_reference_info(args, config)
    reference = resolved.effective
    variant_source = args.get("variant_source", "auto")
    effective_vcf, vcf_primary = _resolve_variant_source(variant_source, args.get("vcf_path"))
    min_alt_count = args.get("min_alt_count", config.min_alt_count)

    data = await _fetch_region_with_timeout(
        file_path, region, reference, config, vcf_path=effective_vcf, vcf_primary=vcf_primary
    )
    data = _apply_alt_count(data, min_alt_count)

    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 2) if coverage else 0
    max_cov = max(coverage) if coverage else 0

    filters = _applied_filters(
        config,
        args,
        include_max_reads=True,
        detection=_detection_status(region, reference, vcf_primary),
        reads_truncated=data.reads_truncated,
        min_alt_count=min_alt_count,
        reference_source=resolved.source,
    )
    summary_lines = [
        f"Region: {format_region(data.contig, data.start, data.end)}",
        f"Reads: {len(data.reads)}",
        f"Coverage: mean={mean_cov}x, max={max_cov}x",
        f"Candidate variants detected: {len(data.variants)}",
        f"Filters applied: {_format_filters_line(filters)}",
        _FILTER_ABSENCE_NOTE,
        f"Coordinate system: {_LLM_COORDINATE_SYSTEM}",
        _CANDIDATE_VARIANT_DISCLAIMER,
    ]
    if resolved.note:
        summary_lines.insert(5, f"Reference: {resolved.note}")

    enhanced, _ = build_enhanced_variants(data)
    for v in enhanced:
        # 1-based position (VCF/dbSNP convention); SV type/span + per-sample genotypes
        # surfaced by _format_variant_line(), then the plain-language reading.
        summary_lines.append(_format_variant_line(v))
        summary_lines.append("    " + annotate_interpretation(v, one_based=False)["interpretation"])

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

    # Invalid input is a tool ERROR (isError), not an answer; not_found/unavailable/disabled
    # stay answers. Validate before the enabled check so bad input fails the same either way.
    build, build_source = _validate_lookup_input(args, config)

    if not config.clinvar_enabled:
        return _lookup_result(
            {"error": "ClinVar lookup is disabled"},
            status=_LOOKUP_DISABLED,
            disclaimer=_CLINVAR_DISCLAIMER,
        )
    build_info = {"genome_build": build, "build_source": build_source}

    client = get_clinvar_client(config)

    try:
        result = await client.lookup(chrom, pos, ref, alt, build=build)
    except NCBIUnavailableError as e:
        # In-band NCBI error, truncated search, unsupported build: absence NOT established.
        logger.warning("ClinVar lookup unavailable for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return _lookup_result(
            {"error": "ClinVar lookup unavailable", "reason": str(e), **build_info},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_CLINVAR_DISCLAIMER,
        )
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError, ValueError) as e:
        # Network/HTTP failure (or an unparseable body) — the lookup is UNAVAILABLE, not an
        # authoritative not-found.
        logger.warning("ClinVar lookup failed for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return _lookup_result(
            {"error": "ClinVar lookup failed", **build_info},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_CLINVAR_DISCLAIMER,
        )
    except asyncio.TimeoutError:
        logger.warning("ClinVar lookup timed out for %s:%d %s>%s", chrom, pos, ref, alt)
        return _lookup_result(
            {"error": "ClinVar lookup timed out"},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_CLINVAR_DISCLAIMER,
        )

    if result is None:
        return _lookup_result(
            {
                "found": False,
                "message": f"No ClinVar entry found for {chrom}:{pos} {ref}>{alt} ({build})",
                **build_info,
            },
            status=_LOOKUP_NOT_FOUND,
            disclaimer=_CLINVAR_DISCLAIMER,
        )

    return _lookup_result(
        {**asdict(result), **build_info}, status=_LOOKUP_FOUND, disclaimer=_CLINVAR_DISCLAIMER
    )


@telemetry_wrapper("lookup_gnomad")
async def handle_lookup_gnomad(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    # Invalid input is a tool ERROR (isError), not an answer; not_found/unavailable/disabled
    # stay answers. Validate before the enabled check so bad input fails the same either way.
    build, build_source = _validate_lookup_input(args, config)

    if not config.gnomad_enabled:
        return _lookup_result(
            {"error": "gnomAD lookup is disabled"},
            status=_LOOKUP_DISABLED,
            disclaimer=_GNOMAD_DISCLAIMER,
        )
    build_info = {"genome_build": build, "build_source": build_source}

    client = get_gnomad_client(config)

    try:
        result = await client.lookup(chrom, pos, ref, alt, build=build)
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError, ValueError) as e:
        # Network/HTTP failure — the lookup is UNAVAILABLE, not an authoritative not-found.
        logger.warning("gnomAD lookup failed for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return _lookup_result(
            {"error": "gnomAD lookup failed"},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_GNOMAD_DISCLAIMER,
        )
    except GnomadQueryError as e:
        # Transient GraphQL error — the client propagates it (rather than returning None) so we
        # report "unavailable" here instead of a misleading "not found in gnomAD".
        # GnomadBuildError (no dataset for this build) is a subclass: also unavailable.
        logger.warning("gnomAD query error for %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
        return _lookup_result(
            {"error": "gnomAD lookup failed", "reason": str(e), **build_info},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_GNOMAD_DISCLAIMER,
        )
    except asyncio.TimeoutError:
        logger.warning("gnomAD lookup timed out for %s:%d %s>%s", chrom, pos, ref, alt)
        return _lookup_result(
            {"error": "gnomAD lookup timed out"},
            status=_LOOKUP_UNAVAILABLE,
            disclaimer=_GNOMAD_DISCLAIMER,
        )

    if result is None:
        return _lookup_result(
            {
                "found": False,
                "message": f"No gnomAD entry found for {chrom}:{pos} {ref}>{alt} ({build})",
                **build_info,
            },
            status=_LOOKUP_NOT_FOUND,
            disclaimer=_GNOMAD_DISCLAIMER,
        )

    return _lookup_result(
        {**asdict(result), **build_info}, status=_LOOKUP_FOUND, disclaimer=_GNOMAD_DISCLAIMER
    )


@dataclass
class _ScanProgress:
    """Shared, lock-guarded state of a progressive contig scan (worker thread → handler).

    Coordinates are internal 0-based half-open. ``scanned_end`` advances only after a chunk's
    candidates are appended, so ``[0, scanned_end)`` is exactly what ``variants`` covers.
    """

    variants: list[dict] = field(default_factory=list)
    scanned_end: int = 0
    scan_end: int = 0
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
    """SNV candidates in one chunk — the same per-chunk logic as ``parsers.scan_contig``."""
    from .parsers import _deletion_depth_at, detect_variants, read_passes_filters

    ref_seq = fasta.fetch(contig, chunk_start, chunk_end)
    counts = samfile.count_coverage(
        contig,
        chunk_start,
        chunk_end,
        quality_threshold=min_baseq,
        read_callback=lambda r: read_passes_filters(r, min_mapq),
    )
    variants = detect_variants(counts, ref_seq, contig, chunk_start, min_vaf, min_depth)
    if variants:
        gaps = _deletion_depth_at(
            samfile, contig, sorted({v["position"] for v in variants}), min_mapq
        )
        for v in variants:
            gap = gaps.get(v["position"], 0)
            if gap:
                v["depth"] += gap
                v["vaf"] = round(v["alt_count"] / v["depth"], 3)
        # Adding deletion depth only lowers VAF: re-apply the threshold.
        variants = [v for v in variants if v["vaf"] >= min_vaf]
    if min_alt_count:
        variants = [v for v in variants if v.get("alt_count", 0) >= min_alt_count]
    return variants


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
) -> None:
    """Scan ``contig`` chunk by chunk into ``progress``, stopping at ``deadline`` (monotonic).

    Runs in the parse pool. Checks the deadline (and the handler's ``stop`` flag, set when the
    hard timeout abandons it) before every chunk, so a slow remote scan ends cleanly with the
    chunks done so far rather than holding a worker until the whole contig is read.

    Raises:
        ValueError: If ``contig`` is not in the alignment header.
    """
    from .parsers import _header_lengths, _open_alignment

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        contig_length = _header_lengths(samfile).get(contig)
        if contig_length is None:
            raise ValueError(f"Contig '{contig}' not found in the alignment header")
        scan_end = min(contig_length, max_region)
        with progress.lock:
            progress.contig_length = contig_length
            progress.scan_end = scan_end
        # Adaptive chunking: start small, double while chunks are fast, halve when slow.
        size = min(SCAN_VARIANTS_MIN_CHUNK_SIZE, chunk_size)
        chunk_start = 0
        with pysam.FastaFile(reference_path) as fasta:
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
    """
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
    await _warm_remote_inputs(file_path, reference, index_path)
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
            ),
            timeout=SCAN_VARIANTS_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError:
        # A single chunk hung past the soft deadline: tell the worker to stop at its next
        # check and answer with what is already complete.
        progress.stop.set()
        hard_timeout = True
    except ValueError as e:
        # Unknown contig: an error, never a silent empty "no variants" result.
        raise ToolInputError(str(e)) from e

    with progress.lock:
        all_variants = list(progress.variants)
        scanned_end = progress.scanned_end
        scan_end = progress.scan_end
        contig_length = progress.contig_length
        deadline_hit = progress.deadline_hit

    if contig_length is None:
        # Never got as far as the header: there is no partial result to give.
        raise ToolInputError(
            f"Scan of {contig} timed out after {SCAN_VARIANTS_TIMEOUT_SECONDS:g}s before any "
            "data was read; retry (the first call to a remote file is the slowest) or use "
            "get_variants on a narrower region."
        )

    truncated = hard_timeout or deadline_hit or scanned_end < scan_end
    truncation_reason = None
    if truncated:
        budget = SCAN_VARIANTS_TIMEOUT_SECONDS if hard_timeout else SCAN_VARIANTS_DEADLINE_SECONDS
        truncation_reason = (
            f"time budget of {budget:g}s reached after scanning {scanned_end:,} of "
            f"{scan_end:,} bp; candidates beyond scanned_region were NOT examined — use "
            "get_variants (or scan again) on the unscanned remainder"
        )

    all_variants.sort(key=lambda v: v["vaf"], reverse=True)
    variants = all_variants[:SCAN_VARIANTS_MAX_RESULTS]
    # 1-based inclusive; None when not a single chunk completed.
    scanned_region = format_region(contig, 0, scanned_end) if scanned_end > 0 else None

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
            "truncated_by_max_region": scan_end < contig_length,
            "max_region": SCAN_VARIANTS_MAX_REGION,
            "truncated_by_deadline": truncated,
            "total_candidates": len(all_variants),
            "results_capped": len(all_variants) > SCAN_VARIANTS_MAX_RESULTS,
            "max_results": SCAN_VARIANTS_MAX_RESULTS,
        }
    )
    if min_alt_count is not None:
        filters.setdefault("min_alt_count", min_alt_count)

    payload: dict[str, Any] = {
        "contig": contig,
        "variants": [_one_based(v) for v in variants],
        "count": len(variants),
        "truncated": truncated,
        "scanned_region": scanned_region,
    }
    if truncation_reason:
        payload["truncation_reason"] = truncation_reason
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
        text = f"Removed {removed} cached index files and cleared the region cache."
    return {"content": [{"type": "text", "text": text}]}

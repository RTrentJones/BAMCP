"""Per-server runtime services: the parse pool, the services registry and the region cache.

Lifecycle objects shared by every tool handler. This module must not import ``core.tools``
(the handlers depend on it, never the reverse).
"""

from __future__ import annotations

import asyncio
import functools
import logging
import os
import time
import weakref
from collections import OrderedDict
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, TypeVar

from ..clients.clinvar import ClinVarClient
from ..clients.genes import GeneClient
from ..clients.gnomad import GnomadClient
from ..clients.ncbi import TokenBucketRateLimiter, ncbi_rate_limiter
from ..config import BAMCPConfig
from ..constants import PARSE_MAX_CONCURRENCY, REGION_CACHE_MAX_BYTES, REGION_CACHE_MAX_ENTRIES
from .cache import BAMIndexCache
from .parsers import MISMATCH_BYTES, READ_BYTES_PER_BASE, READ_FIXED_BYTES, RegionData

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

    Keyed by ``id(config)``. That key cannot alias another config: every entry holds a strong
    reference to its own config (``services.config``), so a registered config is never
    collected and its ``id()`` never reused while the entry exists. An entry leaves only via
    :func:`close_external_clients` or LRU eviction (which closes its HTTP clients, scheduled on
    the loop, so it does not leak connections). Call :func:`close_external_clients` for prompt,
    clean teardown.
    """
    key = id(config)
    services = _services_registry.get(key)
    if services is not None:
        _services_registry.move_to_end(key)
        return services

    services = BAMCPServices.from_config(config)
    _services_registry[key] = services
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
    services = _services_registry.pop(id(config), None)
    if services is not None:
        # Explicit lifecycle removal so a closed server leaves nothing behind.
        await services.aclose()


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
    coverage position, ~1 KB per variant dict, plus the reference string. A VCF variant's
    per-sample entries (``samples``: one dict per sample, ~200 B + ~100 B per FORMAT field,
    plus the ``sample_names`` slot) are counted once per record — ALTs of one record share
    the same ``samples`` dict.
    """
    total = 1024
    for read in getattr(data, "reads", None) or ():
        total += READ_FIXED_BYTES + READ_BYTES_PER_BASE * len(getattr(read, "sequence", "") or "")
        total += MISMATCH_BYTES * len(getattr(read, "mismatches", ()) or ())
    total += 16 * len(getattr(data, "coverage", None) or ())
    variants = getattr(data, "variants", None) or ()
    total += 1024 * len(variants)
    seen_samples: set[int] = set()
    for v in variants:
        samples = v.get("samples") if isinstance(v, dict) else None
        if not samples or id(samples) in seen_samples:
            continue
        seen_samples.add(id(samples))
        for sample_data in samples.values():
            fields = len(sample_data) if isinstance(sample_data, dict) else 1
            total += 200 + 100 * fields + 64
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

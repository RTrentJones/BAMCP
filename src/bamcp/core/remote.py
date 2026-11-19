"""Remote inputs: index/reference-index download into the session cache, and cold-start warming.

Index downloads follow redirects by hand (every hop re-validated against the SSRF policy),
stream under a size cap into an atomic temp file, and are serialized per cache key.
"""

from __future__ import annotations

import asyncio
import contextlib
import logging
import tempfile
from collections import OrderedDict
from collections.abc import AsyncIterator
from pathlib import Path
from urllib.parse import urljoin, urlsplit, urlunsplit

import httpx

from ..config import BAMCPConfig
from ..constants import (
    COLD_START_TIMEOUT_SECONDS,
    INDEX_DOWNLOAD_TIMEOUT_SECONDS,
    WARM_INPUTS_MAX_ENTRIES,
)
from .parsers import open_fasta, register_reference_index
from .services import BAMCPServices, _run_in_parse_pool, get_services
from .validation import is_cram, validate_remote_url

logger = logging.getLogger(__name__)

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


def _sibling_url(url: str, suffix: str, *, replace_ext: bool = False) -> str:
    """The URL of a file next to ``url`` (its index): edit only the URL *path*.

    ``suffix`` is appended to the path (``x.bam`` → ``x.bam.bai``); with ``replace_ext`` the
    final path segment's extension is replaced instead (``x.bam`` → ``x.bai``). The query
    string is kept, so a signed URL (``x.bam?X-Amz-Signature=…``) yields ``x.bam.bai?…``
    rather than a suffix appended inside the query. Fragments are dropped (``validate_remote_url``
    rejects ``#`` anyway).
    """
    parts = urlsplit(url)
    path = parts.path
    if replace_ext:
        head, _, segment = path.rpartition("/")
        if "." in segment:
            segment = segment.rsplit(".", 1)[0]
        path = f"{head}/{segment}" if head or path.startswith("/") else segment
    return urlunsplit((parts.scheme, parts.netloc, path + suffix, parts.query, ""))


@contextlib.asynccontextmanager
async def _index_download_lock(services: BAMCPServices, key: str) -> AsyncIterator[None]:
    """Serialize downloads of one cached index across concurrent callers.

    The per-key lock is dropped once nobody holds it so the dict can't grow one entry per
    URL forever. A waiter already woken keeps its reference; at worst a racing new caller
    creates a fresh lock and repeats an (atomic, idempotent) download.
    """
    locks = services.index_download_locks
    lock = locks.setdefault(key, asyncio.Lock())
    try:
        async with lock:
            yield
    finally:
        if not lock.locked() and locks.get(key) is lock:
            del locks[key]


def _index_http_client() -> httpx.AsyncClient:
    # Redirects are followed by hand in _download_index_streaming (each hop re-validated).
    return httpx.AsyncClient(timeout=INDEX_DOWNLOAD_TIMEOUT_SECONDS, follow_redirects=False)


async def _try_download_index(
    client: httpx.AsyncClient, url: str, path: str, config: BAMCPConfig, what: str
) -> str | None:
    """Download one index into the cache; None (logged) on a miss or a failure."""
    try:
        got = await _download_index_streaming(client, url, path, config)
    except (httpx.RequestError, OSError, ValueError) as e:
        logger.warning("%s download error for %s: %s", what, url, e)
        return None
    if got:
        # The budget sweep lists and stats the whole cache dir: keep it off the event loop.
        await asyncio.to_thread(get_services(config).cache.enforce_budget, keep=got)
    return got


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
    if await asyncio.to_thread(cache.is_valid, index_path):
        logger.debug("Using cached index: %s", index_path)
        return index_path

    # Determine index URL - try common extensions (on the URL path, never inside a query).
    if is_cram(file_path):
        index_urls = [_sibling_url(file_path, ".crai")]
    else:
        # BAM files: try .bam.bai first, then .bai
        index_urls = [
            _sibling_url(file_path, ".bai"),
            _sibling_url(file_path, ".bai", replace_ext=True),
        ]

    async with _index_download_lock(services, index_path):
        if await asyncio.to_thread(cache.is_valid, index_path):
            logger.debug("Using cached index after waiting for download lock: %s", index_path)
            return index_path

        logger.info("Downloading index for remote BAM: %s", file_path)
        async with _index_http_client() as client:
            for index_url in index_urls:
                downloaded = await _try_download_index(
                    client, index_url, index_path, config, "Index"
                )
                if downloaded:
                    return downloaded
                logger.debug("Index not found at %s, trying next", index_url)

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
        with open_fasta(reference) as fasta:
            if fasta.references:
                fasta.fetch(fasta.references[0], 0, 1)
    if file_path:
        with _open_alignment(file_path, reference, index_path) as samfile:
            _ = samfile.references


async def _ensure_reference_indexes(reference: str | None, config: BAMCPConfig) -> None:
    """Cache a remote reference FASTA's .fai (+ .gzi for bgzip) and register them for open_fasta.

    Best effort: a failed download leaves htslib to fetch the index remotely as before.
    """
    if reference is None or not _is_remote(reference) or not config.allow_remote_files:
        return
    services = get_services(config)
    paths = services.cache.get_reference_index_paths(reference)
    if paths is None:
        return
    fai, gzi = paths
    wanted = [(fai, _sibling_url(reference, ".fai"))]
    if gzi:
        wanted.append((gzi, _sibling_url(reference, ".gzi")))
    async with _index_download_lock(services, fai):
        missing = [
            (p, u) for p, u in wanted if not await asyncio.to_thread(services.cache.is_valid, p)
        ]
        if missing:
            async with _index_http_client() as client:
                for path, url in missing:
                    if not await _try_download_index(client, url, path, config, "Reference index"):
                        return
        register_reference_index(reference, fai, gzi)


async def _warm_remote_inputs(
    file_path: str,
    reference: str | None,
    index_path: str | None,
    config: BAMCPConfig,
) -> None:
    """Give a cold remote BAM/CRAM or reference its one-off open under the cold-start budget.

    No-op for local inputs and for inputs already warm, so only the first call to a not-yet-seen
    URL pays (and is allowed) the cold cost. A cold-open timeout raises a diagnosable
    ``TimeoutError``; any other failure is left for the real parse to report with its own
    message (the input simply stays cold).
    """
    await _ensure_reference_indexes(reference, config)
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


def _clear_warm_inputs() -> None:
    """Forget every warm input (cleanup_cache: their cached indexes are gone)."""
    _warm_inputs.clear()

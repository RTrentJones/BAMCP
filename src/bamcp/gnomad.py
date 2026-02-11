"""gnomAD GraphQL API client for population allele frequency data."""

from __future__ import annotations

import asyncio
import logging
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Generic, TypeVar

import httpx

logger = logging.getLogger(__name__)

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

# Cache configuration
CACHE_MAX_SIZE = 1000  # Maximum number of cached entries
CACHE_TTL_SECONDS = 3600  # 1 hour TTL

T = TypeVar("T")


class BoundedTTLCache(Generic[T]):
    """Thread-safe bounded cache with TTL eviction.

    Implements LRU eviction when maxsize is exceeded and TTL-based expiry.
    """

    def __init__(self, maxsize: int = CACHE_MAX_SIZE, ttl: float = CACHE_TTL_SECONDS):
        self._cache: OrderedDict[tuple, tuple[T, float]] = OrderedDict()
        self._maxsize = maxsize
        self._ttl = ttl
        self._lock = asyncio.Lock()

    async def get(self, key: tuple) -> T | None:
        """Get value from cache, returning None if missing or expired."""
        async with self._lock:
            if key not in self._cache:
                return None

            value, timestamp = self._cache[key]
            if time.monotonic() - timestamp > self._ttl:
                del self._cache[key]
                return None

            # Move to end (most recently used)
            self._cache.move_to_end(key)
            return value

    async def set(self, key: tuple, value: T) -> None:
        """Set value in cache, evicting oldest if at capacity."""
        async with self._lock:
            if key in self._cache:
                del self._cache[key]
            elif len(self._cache) >= self._maxsize:
                # Evict oldest (first) item
                self._cache.popitem(last=False)

            self._cache[key] = (value, time.monotonic())


VARIANT_QUERY = """
query GnomadVariant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome {
      ac
      an
      homozygote_count
      af
      populations {
        id
        ac
        an
        homozygote_count
      }
      filters
    }
    exome {
      ac
      an
      homozygote_count
      af
      populations {
        id
        ac
        an
        homozygote_count
      }
      filters
    }
  }
}
"""


@dataclass
class PopulationFrequency:
    """Allele frequency data for a specific population."""

    id: str
    ac: int
    an: int
    homozygote_count: int
    af: float


@dataclass
class GnomadResult:
    """Parsed gnomAD frequency data for a variant."""

    variant_id: str
    global_af: float
    ac: int
    an: int
    homozygote_count: int
    populations: list[PopulationFrequency]
    filters: list[str]
    source: str  # "genome" or "exome"


@dataclass
class GnomadClient:
    """Async client for gnomAD GraphQL API queries.

    Features:
    - Bounded LRU cache with TTL to prevent memory exhaustion
    - Concurrency limiting via semaphore
    """

    dataset: str = "gnomad_r4"
    timeout: float = 30.0
    _cache: BoundedTTLCache[GnomadResult | None] = field(default=None, repr=False)  # type: ignore[assignment]
    _semaphore: asyncio.Semaphore = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        self._semaphore = asyncio.Semaphore(5)
        self._cache = BoundedTTLCache[GnomadResult | None](
            maxsize=CACHE_MAX_SIZE, ttl=CACHE_TTL_SECONDS
        )

    async def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> GnomadResult | None:
        """
        Look up a variant in gnomAD.

        Args:
            chrom: Chromosome (e.g., "chr17" or "17").
            pos: Genomic position (1-based).
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            GnomadResult if found, None otherwise.
        """
        cache_key = (chrom, pos, ref, alt)

        # Check cache first
        cached = await self._cache.get(cache_key)
        if cached is not None:
            return cached

        async with self._semaphore:
            # Re-check cache after acquiring semaphore
            cached = await self._cache.get(cache_key)
            if cached is not None:
                return cached

            result = await self._do_lookup(chrom, pos, ref, alt)
            await self._cache.set(cache_key, result)
            return result

    async def _do_lookup(self, chrom: str, pos: int, ref: str, alt: str) -> GnomadResult | None:
        """Execute the actual gnomAD GraphQL lookup."""
        variant_id = _build_variant_id(chrom, pos, ref, alt)

        payload = {
            "query": VARIANT_QUERY,
            "variables": {
                "variantId": variant_id,
                "dataset": self.dataset,
            },
        }

        async with httpx.AsyncClient(timeout=self.timeout) as client:
            resp = await client.post(GNOMAD_API_URL, json=payload)
            resp.raise_for_status()
            data = resp.json()

        return _parse_response(data, variant_id)


def _build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build a gnomAD variant ID from components."""
    # Strip "chr" prefix for gnomAD variant IDs
    chrom_num = chrom.replace("chr", "")
    return f"{chrom_num}-{pos}-{ref}-{alt}"


def _build_query() -> str:
    """Return the GraphQL query string."""
    return VARIANT_QUERY


def _parse_response(data: dict, variant_id: str) -> GnomadResult | None:
    """Parse a gnomAD GraphQL response into a GnomadResult."""
    if "errors" in data:
        logger.debug("gnomAD GraphQL errors: %s", data["errors"])
        return None

    variant_data = data.get("data", {}).get("variant")
    if not variant_data:
        logger.debug("gnomAD: variant not found for %s", variant_id)
        return None

    # Prefer genome data, fall back to exome
    source_data = variant_data.get("genome") or variant_data.get("exome")
    if not source_data:
        return None

    source = "genome" if variant_data.get("genome") else "exome"

    populations = []
    for pop in source_data.get("populations", []):
        ac = pop.get("ac", 0)
        an = pop.get("an", 0)
        # Calculate af from ac/an (not exposed in API for populations)
        af = ac / an if an > 0 else 0.0
        populations.append(
            PopulationFrequency(
                id=pop.get("id", ""),
                ac=ac,
                an=an,
                homozygote_count=pop.get("homozygote_count", 0),
                af=af,
            )
        )

    return GnomadResult(
        variant_id=variant_data.get("variant_id", variant_id),
        global_af=source_data.get("af", 0.0),
        ac=source_data.get("ac", 0),
        an=source_data.get("an", 0),
        homozygote_count=source_data.get("homozygote_count", 0),
        populations=populations,
        filters=source_data.get("filters", []),
        source=source,
    )

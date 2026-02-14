"""ClinVar NCBI E-utilities API client for variant annotation."""

from __future__ import annotations

import asyncio
import logging
import time
from dataclasses import dataclass, field

import httpx

from .ttl_cache import API_CACHE_MAX_SIZE, API_CACHE_TTL_SECONDS, BoundedTTLCache

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# NCBI rate limits (requests per second)
NCBI_RATE_LIMIT_NO_KEY = 3  # 3 req/sec without API key
NCBI_RATE_LIMIT_WITH_KEY = 10  # 10 req/sec with API key


class TokenBucketRateLimiter:
    """Async token bucket rate limiter for API requests.

    Enforces a maximum number of requests per second.
    """

    def __init__(self, rate: float):
        """Initialize rate limiter.

        Args:
            rate: Maximum requests per second.
        """
        self._rate = rate
        self._tokens = rate
        self._last_refill = time.monotonic()
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        """Acquire a token, waiting if necessary."""
        async with self._lock:
            now = time.monotonic()
            elapsed = now - self._last_refill
            self._tokens = min(self._rate, self._tokens + elapsed * self._rate)
            self._last_refill = now

            if self._tokens < 1:
                wait_time = (1 - self._tokens) / self._rate
                await asyncio.sleep(wait_time)
                self._tokens = 0
            else:
                self._tokens -= 1


# Review status to star rating mapping
REVIEW_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion for the individual variant": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}


@dataclass
class ClinVarResult:
    """Parsed ClinVar annotation for a variant."""

    variation_id: int
    clinical_significance: str
    review_status: str
    stars: int
    conditions: list[str]
    last_evaluated: str | None
    gene: str | None
    variant_name: str | None


@dataclass
class ClinVarClient:
    """Async client for NCBI E-utilities ClinVar queries.

    Features:
    - Bounded LRU cache with TTL to prevent memory exhaustion
    - Token bucket rate limiting to comply with NCBI requirements
    - Concurrency limiting via semaphore
    """

    api_key: str | None = None
    timeout: float = 30.0
    _cache: BoundedTTLCache[ClinVarResult | None] = field(default=None, repr=False)  # type: ignore[arg-type]
    _semaphore: asyncio.Semaphore = field(default=None, repr=False)  # type: ignore[assignment]
    _rate_limiter: TokenBucketRateLimiter = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        max_concurrent = 10 if self.api_key else 3
        rate = NCBI_RATE_LIMIT_WITH_KEY if self.api_key else NCBI_RATE_LIMIT_NO_KEY

        self._semaphore = asyncio.Semaphore(max_concurrent)
        self._rate_limiter = TokenBucketRateLimiter(rate)
        self._cache = BoundedTTLCache[ClinVarResult | None](
            maxsize=API_CACHE_MAX_SIZE, ttl=API_CACHE_TTL_SECONDS
        )

    async def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> ClinVarResult | None:
        """
        Look up a variant in ClinVar via NCBI E-utilities.

        Args:
            chrom: Chromosome (e.g., "chr17" or "17").
            pos: Genomic position (1-based).
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            ClinVarResult if found, None otherwise.
        """
        cache_key = (chrom, pos, ref, alt)

        # Check cache first
        cached = await self._cache.get(cache_key)
        if cached is not None:
            return cached

        # Acquire semaphore for concurrency control and rate limit
        async with self._semaphore:
            # Re-check cache after acquiring semaphore (another request may have filled it)
            cached = await self._cache.get(cache_key)
            if cached is not None:
                return cached

            await self._rate_limiter.acquire()
            result = await self._do_lookup(chrom, pos, ref, alt)
            await self._cache.set(cache_key, result)
            return result

    async def _do_lookup(self, chrom: str, pos: int, ref: str, alt: str) -> ClinVarResult | None:
        """Execute the actual ClinVar API lookup."""
        term = _build_search_term(chrom, pos, ref, alt)

        params: dict[str, str | int] = {
            "db": "clinvar",
            "term": term,
            "retmode": "json",
            "retmax": 5,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        async with httpx.AsyncClient(timeout=self.timeout) as client:
            # Step 1: esearch to find variation IDs
            search_resp = await client.get(f"{EUTILS_BASE}esearch.fcgi", params=params)
            search_resp.raise_for_status()
            search_data = search_resp.json()

            id_list = search_data.get("esearchresult", {}).get("idlist", [])
            if not id_list:
                logger.debug("ClinVar: no results for %s", term)
                return None

            # Rate limit before second request
            await self._rate_limiter.acquire()

            # Step 2: esummary to get annotation details
            summary_params: dict[str, str | int] = {
                "db": "clinvar",
                "id": ",".join(id_list[:5]),
                "retmode": "json",
            }
            if self.api_key:
                summary_params["api_key"] = self.api_key

            summary_resp = await client.get(f"{EUTILS_BASE}esummary.fcgi", params=summary_params)
            summary_resp.raise_for_status()
            summary_data = summary_resp.json()

            return _parse_summary(summary_data)


def _build_search_term(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build a ClinVar search term for a specific variant."""
    # Strip "chr" prefix for ClinVar queries
    chrom_num = chrom.replace("chr", "")
    return f"{chrom_num}[Chromosome] AND {pos}[Base Position] AND {ref}>{alt}[Variant name]"


def _parse_summary(data: dict) -> ClinVarResult | None:
    """Parse an esummary response into a ClinVarResult."""
    result_section = data.get("result", {})
    uids = result_section.get("uids", [])
    if not uids:
        return None

    uid = uids[0]
    entry = result_section.get(uid, {})
    if not entry:
        return None

    clinical_sig = entry.get("clinical_significance", {})
    if isinstance(clinical_sig, dict):
        significance = clinical_sig.get("description", "uncertain significance")
    else:
        significance = str(clinical_sig) if clinical_sig else "uncertain significance"

    review_status = entry.get("review_status", "no assertion criteria provided")
    if isinstance(review_status, dict):
        review_status = review_status.get("description", "no assertion criteria provided")

    stars = REVIEW_STARS.get(review_status.lower(), 0)

    # Extract conditions/traits
    conditions: list[str] = []
    trait_set = entry.get("trait_set", [])
    if isinstance(trait_set, list):
        for trait in trait_set:
            if isinstance(trait, dict):
                trait_name = trait.get("trait_name", "")
                if trait_name:
                    conditions.append(trait_name)

    # Extract gene info
    genes = entry.get("genes", [])
    gene = None
    if isinstance(genes, list) and genes:
        first_gene = genes[0]
        if isinstance(first_gene, dict):
            gene = first_gene.get("symbol")
        elif isinstance(first_gene, str):
            gene = first_gene

    return ClinVarResult(
        variation_id=int(uid),
        clinical_significance=significance,
        review_status=review_status,
        stars=stars,
        conditions=conditions,
        last_evaluated=entry.get("last_evaluated"),
        gene=gene,
        variant_name=entry.get("title"),
    )

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
    _http_client: httpx.AsyncClient = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        max_concurrent = 10 if self.api_key else 3
        rate = NCBI_RATE_LIMIT_WITH_KEY if self.api_key else NCBI_RATE_LIMIT_NO_KEY

        self._semaphore = asyncio.Semaphore(max_concurrent)
        self._rate_limiter = TokenBucketRateLimiter(rate)
        self._cache = BoundedTTLCache[ClinVarResult | None](
            maxsize=API_CACHE_MAX_SIZE, ttl=API_CACHE_TTL_SECONDS
        )
        self._http_client = httpx.AsyncClient(timeout=self.timeout)

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
        """Execute the actual ClinVar API lookup.

        Searches by genomic position (ClinVar's ``[Variant name]`` field is HGVS, not
        ``REF>ALT``, so an allele clause matches nothing) and disambiguates the exact
        allele from each returned record's canonical SPDI.
        """
        params: dict[str, str | int] = {
            "db": "clinvar",
            "term": _build_search_term(chrom, pos),
            "retmode": "json",
            "retmax": 20,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        # Step 1: esearch for every record overlapping the position.
        search_resp = await self._http_client.get(f"{EUTILS_BASE}esearch.fcgi", params=params)
        search_resp.raise_for_status()
        id_list = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            logger.debug("ClinVar: no records at %s:%s", chrom, pos)
            return None

        # Rate limit before second request
        await self._rate_limiter.acquire()

        # Step 2: esummary for annotation details on all candidate records.
        summary_params: dict[str, str | int] = {
            "db": "clinvar",
            "id": ",".join(id_list),
            "retmode": "json",
        }
        if self.api_key:
            summary_params["api_key"] = self.api_key

        summary_resp = await self._http_client.get(
            f"{EUTILS_BASE}esummary.fcgi", params=summary_params
        )
        summary_resp.raise_for_status()

        # Step 3: pick the record whose SPDI matches the requested allele.
        entry = _select_entry(summary_resp.json().get("result", {}), pos, ref, alt)
        if entry is None:
            logger.debug("ClinVar: no allele match for %s:%s %s>%s", chrom, pos, ref, alt)
            return None
        return _parse_entry(entry)

    async def close(self) -> None:
        """Close the HTTP client."""
        await self._http_client.aclose()


def _build_search_term(chrom: str, pos: int) -> str:
    """Build a ClinVar esearch term for records overlapping a genomic position.

    Position only: ClinVar's ``[Variant name]`` field holds HGVS names (``c.1799T>A``),
    not ``REF>ALT``, so a ref>alt clause matches nothing; ``[Base Position]`` indexes the
    GRCh38 coordinate. The exact allele is disambiguated from each record's SPDI.
    """
    chrom_num = chrom.replace("chr", "")
    return f"{chrom_num}[Chromosome] AND {pos}[Base Position]"


def _spdi_matches(spdi: str, pos: int, ref: str, alt: str) -> bool:
    """Return whether a canonical SPDI (``seq:0based_pos:del:ins``) is this variant.

    SPDI's deletion segment is the reference allele and its insertion the alternate; its
    position is 0-based, so it is 1-based ``pos`` when ``spdi_pos + 1 == pos``.
    """
    parts = spdi.split(":")
    if len(parts) != 4:
        return False
    _seq, spdi_pos, deletion, insertion = parts
    try:
        spdi_pos_int = int(spdi_pos)
    except ValueError:
        return False
    return (
        spdi_pos_int + 1 == pos
        and deletion.upper() == ref.upper()
        and insertion.upper() == alt.upper()
    )


def _select_entry(result_section: dict, pos: int, ref: str, alt: str) -> dict | None:
    """Return the esummary entry whose SPDI matches the requested allele.

    A ``[Base Position]`` search returns every record overlapping the position (SNVs,
    indels, delins), so disambiguate to the exact allele via each record's canonical SPDI.
    """
    for uid in result_section.get("uids", []):
        entry: dict = result_section.get(uid, {})
        for variation in entry.get("variation_set", []):
            if _spdi_matches(variation.get("canonical_spdi", ""), pos, ref, alt):
                return entry
    return None


def _classification(entry: dict) -> dict:
    """Normalize an esummary entry's classification block.

    ClinVar's 2024 schema splits classifications into ``germline_classification`` /
    ``oncogenicity_classification`` / ``clinical_impact_classification``; older responses
    carried ``clinical_significance`` / ``review_status`` / ``trait_set`` at the top level.
    Returns a dict with ``description`` / ``review_status`` / ``trait_set`` / ``last_evaluated``.
    """
    for key in (
        "germline_classification",
        "clinical_impact_classification",
        "oncogenicity_classification",
    ):
        block = entry.get(key)
        if isinstance(block, dict) and block.get("description"):
            return {
                "description": block.get("description"),
                "review_status": block.get("review_status"),
                "trait_set": block.get("trait_set", []),
                "last_evaluated": block.get("last_evaluated"),
            }

    # Legacy (pre-2024) top-level fields.
    clinical_sig = entry.get("clinical_significance")
    if isinstance(clinical_sig, dict):
        description = clinical_sig.get("description")
    else:
        description = str(clinical_sig) if clinical_sig else None
    review_status = entry.get("review_status")
    if isinstance(review_status, dict):
        review_status = review_status.get("description")
    return {
        "description": description,
        "review_status": review_status,
        "trait_set": entry.get("trait_set", []),
        "last_evaluated": entry.get("last_evaluated"),
    }


def _parse_entry(entry: dict) -> ClinVarResult | None:
    """Parse a single esummary entry into a ClinVarResult."""
    uid = entry.get("uid")
    if not uid:
        return None

    cls = _classification(entry)
    significance = cls.get("description") or "uncertain significance"
    review_status = cls.get("review_status") or "no assertion criteria provided"
    stars = REVIEW_STARS.get(review_status.lower(), 0)

    # Extract conditions/traits
    conditions: list[str] = []
    trait_set = cls.get("trait_set") or []
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
        last_evaluated=cls.get("last_evaluated"),
        gene=gene,
        variant_name=entry.get("title"),
    )

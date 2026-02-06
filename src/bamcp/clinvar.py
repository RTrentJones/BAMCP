"""ClinVar NCBI E-utilities API client for variant annotation."""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass, field

import httpx

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

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
    """Async client for NCBI E-utilities ClinVar queries."""

    api_key: str | None = None
    timeout: float = 30.0
    _cache: dict[tuple[str, int, str, str], ClinVarResult | None] = field(
        default_factory=dict, repr=False
    )
    _semaphore: asyncio.Semaphore = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        max_concurrent = 10 if self.api_key else 3
        self._semaphore = asyncio.Semaphore(max_concurrent)

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
        if cache_key in self._cache:
            return self._cache[cache_key]

        async with self._semaphore:
            result = await self._do_lookup(chrom, pos, ref, alt)
            self._cache[cache_key] = result
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

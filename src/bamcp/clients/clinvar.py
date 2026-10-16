"""ClinVar NCBI E-utilities API client for variant annotation."""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass, field

import httpx

from ..constants import DEFAULT_GENOME_BUILD, normalize_genome_build
from .ncbi import (
    EUTILS_BASE,
    NCBI_RATE_LIMIT_NO_KEY,
    NCBI_RATE_LIMIT_WITH_KEY,
    NCBIUnavailableError,
    TokenBucketRateLimiter,
    check_esearch,
    check_esummary,
    ncbi_get,
    ncbi_rate_limiter,
    response_json,
)
from .ttl_cache import API_CACHE_MAX_SIZE, API_CACHE_TTL_SECONDS, MISS, BoundedTTLCache

logger = logging.getLogger(__name__)

__all__ = [
    "EUTILS_BASE",
    "NCBI_RATE_LIMIT_NO_KEY",
    "NCBI_RATE_LIMIT_WITH_KEY",
    "ClinVarClient",
    "ClinVarResult",
    "ClinVarUnavailableError",
    "TokenBucketRateLimiter",
]

# Position-search paging: a single base can be overlapped by many indels/CNVs, so page
# esearch (up to a safety cap) and summarize IDs in batches rather than truncating.
_ESEARCH_PAGE = 200
_ESUMMARY_BATCH = 200
_MAX_POSITION_RECORDS = 2000

# ClinVar's per-assembly position index. ``[Base Position]`` is an alias of chrpos38, so a
# GRCh37 coordinate searched there lands on an unrelated GRCh38 base (false "not in ClinVar").
_POSITION_FIELD = {"GRCh37": "chrpos37", "GRCh38": "chrpos38"}

# Germline classification sources, in preference order. Somatic blocks
# (oncogenicity / clinical impact) are reported separately and never feed the
# germline ``clinical_significance`` an ACMG line reads.
_SOMATIC_BLOCKS = {
    "oncogenicity_classification": "oncogenicity",
    "clinical_impact_classification": "clinical_impact",
}


class ClinVarUnavailableError(NCBIUnavailableError):
    """ClinVar could not answer (unsupported build, truncated search, in-band NCBI error)."""


def _chunks(items: list[str], size: int) -> list[list[str]]:
    """Split a list into consecutive chunks of at most ``size`` items."""
    return [items[i : i + size] for i in range(0, len(items), size)]


# Review status to star rating mapping
REVIEW_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion for the individual variant": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}


@dataclass
class ClinVarResult:
    """Parsed ClinVar annotation for a variant.

    ``clinical_significance`` is the GERMLINE classification only, or ``None`` when the
    record carries none (never an invented "uncertain significance"). Somatic
    oncogenicity / clinical-impact classifications live in ``somatic_classifications``.
    """

    variation_id: int
    clinical_significance: str | None
    review_status: str | None
    stars: int
    conditions: list[str]
    last_evaluated: str | None
    gene: str | None
    variant_name: str | None
    # "germline" (2024+ germline_classification block), "legacy" (pre-2024 top-level
    # clinical_significance, germline by construction), or None (no germline classification).
    classification_type: str | None = "germline"
    somatic_classifications: dict[str, str] = field(default_factory=dict)
    # The build the query coordinates were matched in.
    build: str = DEFAULT_GENOME_BUILD


@dataclass
class ClinVarClient:
    """Async client for NCBI E-utilities ClinVar queries.

    Features:
    - Build-aware position search (chrpos37 / chrpos38)
    - Bounded LRU cache with TTL to prevent memory exhaustion
    - Token bucket rate limiting to comply with NCBI requirements (shareable with the
      gene client via ``rate_limiter`` — NCBI's limit is per caller, not per database)
    - Concurrency limiting via semaphore
    """

    api_key: str | None = None
    timeout: float = 30.0
    genome_build: str = DEFAULT_GENOME_BUILD
    rate_limiter: TokenBucketRateLimiter | None = field(default=None, repr=False)
    _cache: BoundedTTLCache[ClinVarResult | None] = field(default=None, repr=False)  # type: ignore[arg-type]
    _semaphore: asyncio.Semaphore = field(default=None, repr=False)  # type: ignore[assignment]
    _rate_limiter: TokenBucketRateLimiter = field(default=None, repr=False)  # type: ignore[assignment]
    _http_client: httpx.AsyncClient = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        max_concurrent = 10 if self.api_key else 3

        self._semaphore = asyncio.Semaphore(max_concurrent)
        self._rate_limiter = self.rate_limiter or ncbi_rate_limiter(self.api_key)
        self._cache = BoundedTTLCache[ClinVarResult | None](
            maxsize=API_CACHE_MAX_SIZE, ttl=API_CACHE_TTL_SECONDS
        )
        self._http_client = httpx.AsyncClient(timeout=self.timeout)

    async def lookup(
        self, chrom: str, pos: int, ref: str, alt: str, build: str | None = None
    ) -> ClinVarResult | None:
        """
        Look up a variant in ClinVar via NCBI E-utilities.

        Args:
            chrom: Chromosome (e.g., "chr17" or "17").
            pos: Genomic position (1-based) in ``build``.
            ref: Reference allele.
            alt: Alternate allele.
            build: Genome build of the coordinates (GRCh37/hg19 or GRCh38/hg38). Defaults
                to the client's ``genome_build`` (the server default) — backward compatible.

        Returns:
            ClinVarResult if found, None for a genuine not-found (cached).

        Raises:
            ClinVarUnavailableError / NCBIUnavailableError: the lookup could not be answered
                (unknown build, in-band NCBI error, truncated search, unparseable body).
                Never cached; report as "unavailable", not absence.
        """
        canonical = normalize_genome_build(build if build is not None else self.genome_build)
        if canonical is None:
            raise ClinVarUnavailableError(
                f"ClinVar lookup needs a GRCh37 or GRCh38 build; got {build!r} — "
                "pin the BAM's build with list_contigs"
            )
        cache_key = (chrom, pos, ref, alt, canonical)

        # Check cache first. get_entry distinguishes a miss from a cached "not found"
        # (None), so a negative result is served from cache instead of re-hitting NCBI.
        cached = await self._cache.get_entry(cache_key)
        if cached is not MISS:
            return cached  # type: ignore[return-value]

        # Acquire semaphore for concurrency control and rate limit
        async with self._semaphore:
            # Re-check cache after acquiring semaphore (another request may have filled it)
            cached = await self._cache.get_entry(cache_key)
            if cached is not MISS:
                return cached  # type: ignore[return-value]

            # Errors propagate BEFORE cache.set, so an outage is never replayed as absence.
            result = await self._do_lookup(chrom, pos, ref, alt, canonical)
            await self._cache.set(cache_key, result)
            return result

    async def _do_lookup(
        self, chrom: str, pos: int, ref: str, alt: str, build: str = DEFAULT_GENOME_BUILD
    ) -> ClinVarResult | None:
        """Execute the actual ClinVar API lookup.

        Searches by genomic position in the requested build (ClinVar's ``[Variant name]``
        field is HGVS, not ``REF>ALT``, so an allele clause matches nothing), spanning the
        whole reference allele so an indel's records (anchored after the VCF pad base) are
        found, pages past the esearch limit for variant-dense positions, and disambiguates
        the exact allele from each candidate record (see :func:`_select_entry`).
        """
        # Step 1: collect every record ID overlapping the allele span (paged).
        id_list, truncated = await self._search_ids(chrom, pos, _span_end(pos, ref, alt), build)
        if not id_list:
            logger.debug("ClinVar: no records at %s:%s (%s)", chrom, pos, build)
            return None

        # Step 2/3: esummary the candidates in batches, returning the first whose allele
        # matches (early-exit avoids summarizing every record).
        for batch in _chunks(id_list, _ESUMMARY_BATCH):
            result_section = await self._esummary(batch)
            entry = _select_entry(result_section, pos, ref, alt, build=build, chrom=chrom)
            if entry is not None:
                parsed = _parse_entry(entry)
                if parsed is not None:
                    parsed.build = build
                return parsed

        if truncated:
            # More records overlap this span than the safety cap fetched; the match may be
            # in the unfetched tail, so absence is NOT established.
            raise ClinVarUnavailableError(
                f"ClinVar search at {chrom}:{pos} exceeded {_MAX_POSITION_RECORDS} records; "
                "result truncated, absence not established"
            )
        logger.debug("ClinVar: no allele match for %s:%s %s>%s", chrom, pos, ref, alt)
        return None

    async def _search_ids(
        self, chrom: str, pos: int, end: int | None = None, build: str = DEFAULT_GENOME_BUILD
    ) -> tuple[list[str], bool]:
        """Return ClinVar record IDs overlapping ``pos..end``, and whether truncated."""
        term = _build_search_term(chrom, pos, build, end)
        ids: list[str] = []
        retstart = 0
        count = 0
        while retstart < _MAX_POSITION_RECORDS:
            params: dict[str, str | int] = {
                "db": "clinvar",
                "term": term,
                "retmode": "json",
                "retmax": min(_ESEARCH_PAGE, _MAX_POSITION_RECORDS - retstart),
                "retstart": retstart,
            }
            if self.api_key:
                params["api_key"] = self.api_key
            resp = await ncbi_get(
                self._http_client, f"{EUTILS_BASE}esearch.fcgi", params, self._rate_limiter
            )
            resp.raise_for_status()
            search_result = check_esearch(response_json(resp, "ClinVar esearch"), "ClinVar")
            try:
                count = int(search_result.get("count", "0"))
            except (TypeError, ValueError) as e:
                raise ClinVarUnavailableError("ClinVar esearch returned a malformed count") from e
            page_ids = search_result.get("idlist", [])
            if not page_ids:
                break
            ids.extend(page_ids)
            retstart += len(page_ids)
            if retstart >= count:
                break
        return ids, count > len(ids)

    async def _esummary(self, ids: list[str]) -> dict:
        """Fetch esummary details for a batch of record IDs; return the ``result`` map."""
        params: dict[str, str | int] = {
            "db": "clinvar",
            "id": ",".join(ids),
            "retmode": "json",
        }
        if self.api_key:
            params["api_key"] = self.api_key
        resp = await ncbi_get(
            self._http_client, f"{EUTILS_BASE}esummary.fcgi", params, self._rate_limiter
        )
        resp.raise_for_status()
        return check_esummary(response_json(resp, "ClinVar esummary"), "ClinVar")

    async def close(self) -> None:
        """Close the HTTP client."""
        await self._http_client.aclose()


def _span_end(pos: int, ref: str, alt: str) -> int:
    """Last reference base a variant's ClinVar records can sit on.

    SNV/MNV: the reference allele's own span. Indel: VCF ``pos`` is the pad base BEFORE
    the event while ClinVar indexes the event itself, so search ``pos..pos+len(ref)``.
    """
    if len(ref) == len(alt):
        return pos + max(len(ref), 1) - 1
    return pos + max(len(ref), 1)


def _build_search_term(
    chrom: str, pos: int, build: str = DEFAULT_GENOME_BUILD, end: int | None = None
) -> str:
    """Build a ClinVar esearch term for records overlapping ``pos..end`` in ``build``.

    Position only: ClinVar's ``[Variant name]`` field holds HGVS names (``c.1799T>A``),
    not ``REF>ALT``, so a ref>alt clause matches nothing. The position field is per build
    (``chrpos37`` / ``chrpos38``; ``[Base Position]`` is GRCh38-only). The exact allele is
    disambiguated per record afterwards.
    """
    chrom_num = chrom.replace("chr", "")
    field_name = _POSITION_FIELD[build]
    span = f"{pos}:{end}" if end is not None and end > pos else str(pos)
    return f"{chrom_num}[Chromosome] AND {span}[{field_name}]"


def _minimal(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Reduce ``(pos, ref, alt)`` to a minimal representation by trimming shared bases.

    Trims a shared suffix, then a shared prefix (advancing ``pos``), so a VCF-style
    anchored indel (``pos=100 ref=A alt=AT``) and a canonical SPDI (``100::T``) reduce to
    the same key. Left-aligns to the trimmed edge only — it does not roll through repeat
    runs, so an indel inside a homopolymer/tandem repeat may still miss.
    """
    ref, alt = ref.upper(), alt.upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _spdi_matches(spdi: str, pos: int, ref: str, alt: str, shift: int = 0) -> bool:
    """Return whether a canonical SPDI (``seq:0based_pos:del:ins``) is this variant.

    SPDI's deletion segment is the reference allele and its insertion the alternate at a
    0-based interbase position; the requested variant is VCF-style. Both are reduced to a
    minimal representation (see :func:`_minimal`) before comparison so SNVs and simple
    (non-repeat) indels match despite VCF's leading anchor base. ``shift`` translates the
    SPDI's (always GRCh38) position into the query's assembly.
    """
    parts = spdi.split(":")
    if len(parts) != 4:
        return False
    _seq, spdi_pos, deletion, insertion = parts
    try:
        spdi_pos_int = int(spdi_pos)
    except ValueError:
        return False
    return _minimal(spdi_pos_int + 1 + shift, deletion, insertion) == _minimal(pos, ref, alt)


def _loc(variation: dict, assembly: str) -> dict | None:
    """The ``variation_loc`` entry for an assembly (``GRCh37`` / ``GRCh38``), if any."""
    for loc in variation.get("variation_loc") or []:
        if isinstance(loc, dict) and loc.get("assembly_name") == assembly:
            return loc
    return None


def _grch37_shift(variation: dict, chrom: str | None) -> int | None:
    """Offset that maps this record's GRCh38 coordinates onto GRCh37, or None.

    ClinVar's canonical SPDI is GRCh38-only, but each record's ``variation_loc`` lists its
    span on both assemblies. Within a record the two spans are a rigid translation (same
    length), so ``start37 - start38`` lifts the SPDI onto GRCh37 for allele comparison. A
    record without a GRCh37 location, on a different chromosome, or whose spans differ in
    length (a rearranged region) is not matched rather than guessed.
    """
    loc37, loc38 = _loc(variation, "GRCh37"), _loc(variation, "GRCh38")
    if loc37 is None or loc38 is None:
        return None
    if chrom is not None and str(loc37.get("chr", "")) != chrom.replace("chr", ""):
        return None
    try:
        s37, e37 = int(loc37["start"]), int(loc37["stop"])
        s38, e38 = int(loc38["start"]), int(loc38["stop"])
    except (KeyError, TypeError, ValueError):
        return None
    if e37 - s37 != e38 - s38:
        return None
    return s37 - s38


def _select_entry(
    result_section: dict,
    pos: int,
    ref: str,
    alt: str,
    build: str = DEFAULT_GENOME_BUILD,
    chrom: str | None = None,
) -> dict | None:
    """Return the esummary entry that is exactly the requested allele.

    A position search returns every record overlapping the span (SNVs, indels, delins),
    so disambiguate via each record's canonical SPDI — lifted onto GRCh37 through the
    record's own ``variation_loc`` when the query is GRCh37. Only single-variant records
    qualify: a haplotype/genotype record's ``variation_set`` holds several component
    variants, and its aggregate classification must not be attributed to one allele.
    """
    for uid in result_section.get("uids", []):
        entry: dict = result_section.get(uid, {})
        variation_set = entry.get("variation_set", [])
        if len(variation_set) != 1:
            continue
        variation = variation_set[0]
        spdi = variation.get("canonical_spdi", "")
        if build == "GRCh37":
            shift = _grch37_shift(variation, chrom)
            if shift is None:
                continue
        else:
            shift = 0
        if _spdi_matches(spdi, pos, ref, alt, shift=shift):
            return entry
    return None


def _block_description(block: object) -> str | None:
    if isinstance(block, dict):
        desc = block.get("description")
        if isinstance(desc, str) and desc.strip():
            return desc.strip()
    return None


def _classification(entry: dict) -> dict:
    """Normalize an esummary entry's GERMLINE classification block.

    ClinVar's 2024 schema splits classifications into ``germline_classification`` /
    ``oncogenicity_classification`` / ``clinical_impact_classification``; older responses
    carried ``clinical_significance`` / ``review_status`` / ``trait_set`` at the top level
    (germline by construction). Somatic blocks are NEVER used here — a somatic
    "Oncogenic" / "Tier I" call is not a germline pathogenicity assertion.

    Returns ``description`` / ``review_status`` / ``trait_set`` / ``last_evaluated`` and
    ``type`` (``"germline"`` / ``"legacy"`` / ``None`` when there is no germline call).
    """
    block = entry.get("germline_classification")
    if isinstance(block, dict) and _block_description(block) is not None:
        return {
            "description": _block_description(block),
            "review_status": block.get("review_status") or None,
            "trait_set": block.get("trait_set", []),
            "last_evaluated": block.get("last_evaluated"),
            "type": "germline",
        }

    # Legacy (pre-2024) top-level fields.
    clinical_sig = entry.get("clinical_significance")
    if isinstance(clinical_sig, dict):
        description = _block_description(clinical_sig)
    else:
        description = str(clinical_sig).strip() if clinical_sig else None
    if description:
        review_status = entry.get("review_status")
        if isinstance(review_status, dict):
            review_status = review_status.get("description")
        return {
            "description": description,
            "review_status": review_status or None,
            "trait_set": entry.get("trait_set", []),
            "last_evaluated": entry.get("last_evaluated"),
            "type": "legacy",
        }
    return {
        "description": None,
        "review_status": None,
        "trait_set": [],
        "last_evaluated": None,
        "type": None,
    }


def _somatic_classifications(entry: dict) -> dict[str, str]:
    """Somatic classifications (oncogenicity / clinical impact), kept apart from germline."""
    out: dict[str, str] = {}
    for key, label in _SOMATIC_BLOCKS.items():
        desc = _block_description(entry.get(key))
        if desc is not None:
            out[label] = desc
    return out


def _parse_entry(entry: dict) -> ClinVarResult | None:
    """Parse a single esummary entry into a ClinVarResult.

    A record with no germline classification yields ``clinical_significance=None`` — never
    a fabricated "uncertain significance" (that would read as a real VUS assertion).
    """
    uid = entry.get("uid")
    if not uid:
        return None

    cls = _classification(entry)
    significance = cls.get("description")
    review_status = cls.get("review_status")
    stars = REVIEW_STARS.get(review_status.lower(), 0) if review_status else 0

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
        classification_type=cls.get("type"),
        somatic_classifications=_somatic_classifications(entry),
    )

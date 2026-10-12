"""gnomAD GraphQL API client for population allele frequency data."""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass, field
from typing import Any

import httpx

from ..constants import (
    DEFAULT_GENOME_BUILD,
    GNOMAD_DATASET_BUILD,
    GNOMAD_DATASET_BY_BUILD,
    normalize_genome_build,
)
from .ttl_cache import API_CACHE_MAX_SIZE, API_CACHE_TTL_SECONDS, MISS, BoundedTTLCache

logger = logging.getLogger(__name__)

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

# Population-max eligibility (gnomAD/ClinGen convention for BA1/BS1): only the top-level
# continental genetic-ancestry groups, excluding founder/bottlenecked groups (asj, fin,
# ami, mid) and the catch-all remaining/oth. Sex strata (``nfe_XX``), sub-continental
# groups (``nfe_seu``) and HGDP/1KG cohorts (``hgdp:adygei``) are never eligible: a single
# allele in an AN=14 cohort reads as AF 0.07 — a false "common" (BA1) call.
POPMAX_ELIGIBLE_GROUPS = frozenset({"afr", "amr", "eas", "nfe", "sas"})
# Minimum allele number for a group's AF to count toward the population maximum.
POPMAX_MIN_AN = 2000
# Top-level groups reported in ``populations`` (incl. the founder groups, for display).
TOP_LEVEL_GROUPS = frozenset(
    {"afr", "amr", "asj", "eas", "fin", "nfe", "sas", "ami", "mid", "remaining", "oth"}
)
# GraphQL error message gnomAD returns (with ``data.variant = null``) for a well-formed
# variant id that simply is not in the dataset — a genuine not-found, not an outage.
_VARIANT_NOT_FOUND = "variant not found"

_FREQ_FIELDS = """
      ac
      an
      homozygote_count
      filters
      faf95 {
        popmax
        popmax_population
      }
      populations {
        id
        ac
        an
        homozygote_count
      }
"""

# ``joint`` (gnomAD v4 combined exome+genome frequencies) is null for datasets that do not
# carry it (e.g. gnomad_r2_1), so one query serves every release.
VARIANT_QUERY = (
    """
query GnomadVariant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome {"""
    + _FREQ_FIELDS
    + """    }
    exome {"""
    + _FREQ_FIELDS
    + """    }
    joint {"""
    + _FREQ_FIELDS
    + """    }
  }
}
"""
)


@dataclass
class PopulationFrequency:
    """Allele frequency data for a specific population.

    ``af`` is ``None`` when the group has no called alleles (``an == 0``): an unknown
    frequency, never a measured 0.0.
    """

    id: str
    ac: int
    an: int
    homozygote_count: int
    af: float | None


@dataclass
class GnomadResult:
    """Parsed gnomAD frequency data for a variant.

    ``max_pop_af`` is the population maximum over :data:`POPMAX_ELIGIBLE_GROUPS` with
    ``an >= POPMAX_MIN_AN`` only — use it (or ``faf95_popmax``, gnomAD's own filtering AF)
    for BA1/BS1, never a max over ``populations``. ``frequency_reliable`` is False when the
    site is filtered in every frequency source or too thinly called to establish rarity;
    ``warnings`` says why.
    """

    variant_id: str
    global_af: float | None
    ac: int
    an: int
    homozygote_count: int
    populations: list[PopulationFrequency]
    filters: list[str]
    source: str  # "joint", "genome+exome", "genome" or "exome"
    dataset: str = "gnomad_r4"
    build: str = DEFAULT_GENOME_BUILD
    max_pop_af: float | None = None
    max_pop: str | None = None
    max_pop_an: int | None = None
    faf95_popmax: float | None = None
    faf95_popmax_population: str | None = None
    frequency_reliable: bool = True
    warnings: list[str] = field(default_factory=list)


class GnomadQueryError(Exception):
    """A gnomAD GraphQL error response (transient) — distinct from a genuine not-found."""


class GnomadBuildError(GnomadQueryError):
    """The lookup's genome build has no matching gnomAD dataset — refuse, never guess."""


def resolve_dataset(configured: str, build: str | None) -> tuple[str, str]:
    """Pick the gnomAD dataset for a lookup in ``build``; return ``(dataset, build)``.

    ``configured == "auto"`` maps the build through :data:`GNOMAD_DATASET_BY_BUILD`. An
    explicit dataset is honored only when its assembly is the lookup's build — querying
    gnomad_r4 (GRCh38) with GRCh37 coordinates silently answers for a different locus.

    Raises:
        GnomadBuildError: unknown build, or an explicit dataset in another assembly.
    """
    canonical = normalize_genome_build(build)
    if canonical is None:
        raise GnomadBuildError(
            f"gnomAD lookup needs a GRCh37 or GRCh38 build; got {build!r} — "
            "pin the BAM's build with list_contigs"
        )
    if not configured or configured == "auto":
        return GNOMAD_DATASET_BY_BUILD[canonical], canonical
    dataset_build = GNOMAD_DATASET_BUILD.get(configured)
    if dataset_build != canonical:
        raise GnomadBuildError(
            f"configured gnomAD dataset {configured!r} is "
            f"{dataset_build or 'of unknown assembly'}, but the lookup is {canonical}"
        )
    return configured, canonical


@dataclass
class GnomadClient:
    """Async client for gnomAD GraphQL API queries.

    Features:
    - Build-aware dataset selection (gnomad_r2_1 for GRCh37, gnomad_r4 for GRCh38)
    - Bounded LRU cache with TTL to prevent memory exhaustion
    - Concurrency limiting via semaphore
    """

    dataset: str = "auto"
    timeout: float = 30.0
    genome_build: str = DEFAULT_GENOME_BUILD
    _cache: BoundedTTLCache[GnomadResult | None] = field(default=None, repr=False)  # type: ignore[arg-type]
    _semaphore: asyncio.Semaphore = field(default=None, repr=False)  # type: ignore[assignment]
    _http_client: httpx.AsyncClient = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        self._semaphore = asyncio.Semaphore(5)
        self._cache = BoundedTTLCache[GnomadResult | None](
            maxsize=API_CACHE_MAX_SIZE, ttl=API_CACHE_TTL_SECONDS
        )
        self._http_client = httpx.AsyncClient(timeout=self.timeout)

    async def lookup(
        self, chrom: str, pos: int, ref: str, alt: str, build: str | None = None
    ) -> GnomadResult | None:
        """
        Look up a variant in gnomAD.

        Args:
            chrom: Chromosome (e.g., "chr17" or "17").
            pos: Genomic position (1-based) in ``build``.
            ref: Reference allele.
            alt: Alternate allele.
            build: Genome build of the coordinates. Defaults to the client's
                ``genome_build`` (the server default) — backward compatible.

        Returns:
            GnomadResult if found, None for a genuine not-found (cached).

        Raises:
            GnomadBuildError: no gnomAD dataset for the build — do NOT treat as absence.
            GnomadQueryError: on a transient GraphQL error — do NOT treat as absence.
                Not cached, so the next lookup retries.
        """
        dataset, canonical = resolve_dataset(
            self.dataset, build if build is not None else self.genome_build
        )
        cache_key = (chrom, pos, ref, alt, dataset)

        # Check cache first. get_entry distinguishes a miss from a cached "not found"
        # (None), so a negative result is served from cache instead of re-hitting gnomAD.
        cached = await self._cache.get_entry(cache_key)
        if cached is not MISS:
            return cached  # type: ignore[return-value]

        async with self._semaphore:
            # Re-check cache after acquiring semaphore
            cached = await self._cache.get_entry(cache_key)
            if cached is not MISS:
                return cached  # type: ignore[return-value]

            # A transient GraphQL error propagates (GnomadQueryError) so the caller can report
            # it as "unavailable" — distinct from a genuine not-found (None). It is raised BEFORE
            # cache.set, so it is never cached: the next lookup retries instead of replaying the
            # error as absence for the TTL.
            result = await self._do_lookup(chrom, pos, ref, alt, dataset)
            if result is not None:
                result.dataset = dataset
                result.build = canonical
            await self._cache.set(cache_key, result)
            return result

    async def _do_lookup(
        self, chrom: str, pos: int, ref: str, alt: str, dataset: str | None = None
    ) -> GnomadResult | None:
        """Execute the actual gnomAD GraphQL lookup."""
        variant_id = _build_variant_id(chrom, pos, ref, alt)

        payload = {
            "query": VARIANT_QUERY,
            "variables": {
                "variantId": variant_id,
                "dataset": dataset or GNOMAD_DATASET_BY_BUILD[DEFAULT_GENOME_BUILD],
            },
        }

        resp = await self._http_client.post(GNOMAD_API_URL, json=payload)
        resp.raise_for_status()
        try:
            data = resp.json()
        except ValueError as e:
            raise GnomadQueryError(f"gnomAD returned an unparseable body for {variant_id}") from e

        return _parse_response(data, variant_id)

    async def close(self) -> None:
        """Close the HTTP client."""
        await self._http_client.aclose()


def _build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build a gnomAD variant ID from components."""
    # Strip "chr" prefix for gnomAD variant IDs
    chrom_num = chrom.replace("chr", "")
    return f"{chrom_num}-{pos}-{ref}-{alt}"


def _only_variant_not_found(errors: Any) -> bool:
    """True when a GraphQL error list is exactly gnomAD's "Variant not found" answer."""
    if not isinstance(errors, list) or not errors:
        return False
    return all(
        isinstance(err, dict) and str(err.get("message", "")).strip().lower() == _VARIANT_NOT_FOUND
        for err in errors
    )


def _af(ac: int, an: int) -> float | None:
    return ac / an if an > 0 else None


def _int(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _sum_blocks(blocks: list[dict]) -> tuple[int, int, int, dict[str, list[int]]]:
    """Sum ac / an / homozygotes overall and per population id across frequency blocks."""
    ac = an = hom = 0
    pops: dict[str, list[int]] = {}
    for block in blocks:
        ac += _int(block.get("ac"))
        an += _int(block.get("an"))
        hom += _int(block.get("homozygote_count"))
        seen: set[str] = set()
        for pop in block.get("populations") or []:
            pid = str(pop.get("id", ""))
            # gnomAD's joint block repeats some ids (XX/XY); count each once per block.
            if pid in seen:
                continue
            seen.add(pid)
            agg = pops.setdefault(pid, [0, 0, 0])
            agg[0] += _int(pop.get("ac"))
            agg[1] += _int(pop.get("an"))
            agg[2] += _int(pop.get("homozygote_count"))
    return ac, an, hom, pops


def _parse_response(data: dict, variant_id: str) -> GnomadResult | None:
    """Parse a gnomAD GraphQL response into a GnomadResult.

    gnomAD answers a well-formed id that is absent from the dataset with
    ``{"errors": [{"message": "Variant not found"}], "data": {"variant": null}}`` — that is a
    genuine not-found (``None``, cacheable; PM2 depends on it). Any OTHER error raises
    :class:`GnomadQueryError`: a transient/query error is not absence and must not be cached.

    Frequencies come from ``joint`` (v4 exome+genome) when present, else from the genome and
    exome blocks combined — a filtered block is dropped from the combination when an
    unfiltered one exists, and every block's filters are surfaced (``"exome:AC0"``).
    """
    errors = data.get("errors")
    variant_data = (data.get("data") or {}).get("variant")
    if errors:
        if variant_data is None and _only_variant_not_found(errors):
            logger.debug("gnomAD: variant not found for %s", variant_id)
            return None
        raise GnomadQueryError(f"gnomAD GraphQL errors for {variant_id}: {errors}")

    if not variant_data:
        logger.debug("gnomAD: variant not found for %s", variant_id)
        return None

    present = {
        name: variant_data[name]
        for name in ("joint", "genome", "exome")
        if isinstance(variant_data.get(name), dict)
    }
    if not present:
        return None

    filters = [
        f"{name}:{flt}" for name, block in present.items() for flt in block.get("filters") or []
    ]
    warnings: list[str] = []

    if "joint" in present:
        used = {"joint": present["joint"]}
    else:
        seq = {k: v for k, v in present.items() if k in ("genome", "exome")}
        passing = {k: v for k, v in seq.items() if not v.get("filters")}
        used = passing or seq
        for name in seq.keys() - used.keys():
            warnings.append(
                f"{name} data excluded from frequencies: filtered {seq[name]['filters']}"
            )
    source = "+".join(sorted(used, key=["joint", "genome", "exome"].index))
    if all(block.get("filters") for block in used.values()):
        warnings.append(f"site fails gnomAD QC filters in {source} ({filters})")

    ac, an, hom, pop_sums = _sum_blocks(list(used.values()))

    populations = [
        PopulationFrequency(id=pid, ac=v[0], an=v[1], homozygote_count=v[2], af=_af(v[0], v[1]))
        for pid, v in pop_sums.items()
        if pid in TOP_LEVEL_GROUPS
    ]

    best: PopulationFrequency | None = None
    for pop in populations:
        if pop.id not in POPMAX_ELIGIBLE_GROUPS or pop.an < POPMAX_MIN_AN or pop.af is None:
            continue
        if best is None or pop.af > (best.af or 0.0):
            best = pop

    # gnomAD's own filtering allele frequency (95% CI lower bound of the popmax). For a
    # combined genome+exome view use the better-powered (larger-AN) block's value.
    faf_block = max(used.values(), key=lambda b: _int(b.get("an")))
    faf95 = faf_block.get("faf95") or {}
    faf95_popmax = faf95.get("popmax") if isinstance(faf95, dict) else None
    faf95_pop = faf95.get("popmax_population") if isinstance(faf95, dict) else None

    global_af = _af(ac, an)
    if global_af is None:
        warnings.append("no called alleles (AN=0): allele frequency unknown, not zero")
    elif an < POPMAX_MIN_AN:
        warnings.append(f"low allele number (AN={an}): frequency too poorly sampled to show rarity")

    return GnomadResult(
        variant_id=variant_data.get("variant_id", variant_id),
        global_af=global_af,
        ac=ac,
        an=an,
        homozygote_count=hom,
        populations=populations,
        filters=filters,
        source=source,
        max_pop_af=best.af if best else None,
        max_pop=best.id if best else None,
        max_pop_an=best.an if best else None,
        faf95_popmax=faf95_popmax,
        faf95_popmax_population=faf95_pop,
        frequency_reliable=not warnings,
        warnings=warnings,
    )

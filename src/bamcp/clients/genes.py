"""Gene annotation lookup via NCBI Entrez.

Provides gene symbol to genomic coordinate mapping using NCBI's Gene database.

Build provenance matters more than anything else here: a GRCh38 locus shown in a GRCh37
BAM is ~0.5-2 Mb off the real gene, silently. NCBI's esummary ``genomicinfo`` carries NO
assembly accession (only ``chrloc``, ``chraccver``, ``chrstart``, ``chrstop``,
``exoncount``), so the build is inferred from the RefSeq chromosome accession VERSION
(``NC_000017.10`` is GRCh37, ``.11`` is GRCh38). GRCh37 coordinates come from
``locationhist``, which does carry ``assemblyaccver``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import httpx

from ..constants import DEFAULT_GENOME_BUILD, normalize_genome_build
from .ncbi import (
    EUTILS_BASE,
    NCBIUnavailableError,
    TokenBucketRateLimiter,
    check_esearch,
    check_esummary,
    ncbi_rate_limiter,
    response_json,
)

NCBI_GENE_SEARCH = f"{EUTILS_BASE}esearch.fcgi"
NCBI_GENE_SUMMARY = f"{EUTILS_BASE}esummary.fcgi"

# How many esearch hits to consider when resolving a symbol (an alias such as "p53"
# returns several genes; the official symbol is chosen, not idlist[0]).
_MAX_CANDIDATES = 20


class GeneLookupError(NCBIUnavailableError):
    """NCBI Gene could not be queried (HTTP / in-band error) — not a "gene not found"."""


@dataclass
class GeneInfo:
    """Gene annotation with genomic coordinates (1-based, inclusive)."""

    symbol: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    # The genome build the returned coordinates are actually in. Callers MUST reconcile this
    # against the build of the BAM they intend to view — a GRCh38 locus in a GRCh37 BAM is
    # ~0.5-2 Mb off the real gene, silently. "unknown" must be surfaced as a warning.
    build: str = "unknown"
    # Set when the query matched an alias rather than the official symbol.
    note: str | None = None
    warnings: list[str] = field(default_factory=list)


# RefSeq chromosome accession -> (chrom name, {accession version: build}).
_REFSEQ_CHROMOSOMES: dict[str, tuple[str, dict[int, str]]] = {
    "NC_000001": ("chr1", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000002": ("chr2", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000003": ("chr3", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000004": ("chr4", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000005": ("chr5", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000006": ("chr6", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000007": ("chr7", {13: "GRCh37", 14: "GRCh38"}),
    "NC_000008": ("chr8", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000009": ("chr9", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000010": ("chr10", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000011": ("chr11", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000012": ("chr12", {11: "GRCh37", 12: "GRCh38"}),
    "NC_000013": ("chr13", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000014": ("chr14", {8: "GRCh37", 9: "GRCh38"}),
    "NC_000015": ("chr15", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000016": ("chr16", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000017": ("chr17", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000018": ("chr18", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000019": ("chr19", {9: "GRCh37", 10: "GRCh38"}),
    "NC_000020": ("chr20", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000021": ("chr21", {8: "GRCh37", 9: "GRCh38"}),
    "NC_000022": ("chr22", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000023": ("chrX", {10: "GRCh37", 11: "GRCh38"}),
    "NC_000024": ("chrY", {9: "GRCh37", 10: "GRCh38"}),
    # The rCRS mitochondrion (NC_012920.1) is shared by both builds — no build signal.
    "NC_012920": ("chrM", {}),
}


def _build_from_chraccver(chraccver: str) -> str:
    """Infer the build from a RefSeq chromosome accession version (NC_000017.10 → GRCh37)."""
    base, _, version = chraccver.partition(".")
    entry = _REFSEQ_CHROMOSOMES.get(base)
    if entry is None:
        return "unknown"
    try:
        return entry[1].get(int(version), "unknown")
    except ValueError:
        return "unknown"


# NCBI's human primary assembly is GCF_000001405.<minor>: GRCh37 is .13-.25, GRCh38 is .26+.
# (The old check `"GCF_000001405" in assembly` for GRCh38 also matched GRCh37's .25.)
def _build_from_assembly(assembly: str) -> str:
    """Map an NCBI assembly accession (e.g. GCF_000001405.40) to a genome build."""
    if not assembly.startswith("GCF_000001405."):
        return "unknown"
    try:
        minor = int(assembly.rsplit(".", 1)[-1])
    except ValueError:
        return "unknown"
    if minor <= 25:
        return "GRCh37"
    return "GRCh38"


def _location_build(info: dict) -> str:
    """Build of one genomicinfo/locationhist entry: assembly accession, else chraccver."""
    build = _build_from_assembly(str(info.get("assemblyaccver", "")))
    if build == "unknown":
        build = _build_from_chraccver(str(info.get("chraccver", "")))
    return build


def _is_live_human(entry: dict) -> bool:
    """A current (not replaced/discontinued) Homo sapiens gene record."""
    status = str(entry.get("status", "0") or "0")
    organism = entry.get("organism") or {}
    species = organism.get("scientificname") if isinstance(organism, dict) else None
    return status == "0" and species in (None, "Homo sapiens")


def _aliases(entry: dict) -> set[str]:
    raw = entry.get("otheraliases") or ""
    return {a.strip().upper() for a in str(raw).split(",") if a.strip()}


def _choose_gene(result: dict, ids: list[str], symbol: str) -> tuple[dict | None, str | None]:
    """Pick the gene whose OFFICIAL symbol is ``symbol``; else a unique alias match.

    Returns ``(entry, note)``; ``(None, None)`` when nothing (or nothing unambiguous) matches.
    """
    wanted = symbol.strip().upper()
    live = [result[i] for i in ids if isinstance(result.get(i), dict) and _is_live_human(result[i])]
    for entry in live:
        if str(entry.get("name", "")).upper() == wanted:
            return entry, None
    by_alias = [e for e in live if wanted in _aliases(e)]
    if len(by_alias) == 1:
        entry = by_alias[0]
        return entry, f"{symbol!r} is an alias; resolved to official symbol {entry.get('name')!r}"
    return None, None


class GeneClient:
    """NCBI Gene database client for coordinate lookup.

    Uses NCBI E-utilities to search for human genes by symbol and
    retrieve their genomic coordinates for a specific genome build.
    """

    def __init__(
        self,
        api_key: str | None = None,
        genome_build: str = DEFAULT_GENOME_BUILD,
        timeout: float = 10.0,
        rate_limiter: TokenBucketRateLimiter | None = None,
    ):
        """Initialize the gene client.

        Args:
            api_key: NCBI API key for higher rate limits.
            genome_build: Target genome build (GRCh38 or GRCh37).
            timeout: HTTP request timeout in seconds.
            rate_limiter: Shared NCBI token bucket (share it with the ClinVar client —
                NCBI's limit is per caller). One is created when omitted.
        """
        self.api_key = api_key
        self.genome_build = genome_build
        self._rate_limiter = rate_limiter or ncbi_rate_limiter(api_key)
        self._client = httpx.AsyncClient(timeout=timeout)

    async def _get(self, url: str, params: dict, what: str) -> object:
        await self._rate_limiter.acquire()
        if self.api_key:
            params["api_key"] = self.api_key
        try:
            resp = await self._client.get(url, params=params)
            resp.raise_for_status()
        except httpx.HTTPError as e:
            raise GeneLookupError(f"{what}: NCBI request failed ({type(e).__name__})") from e
        try:
            return response_json(resp, what)
        except NCBIUnavailableError as e:
            raise GeneLookupError(str(e)) from e

    async def search(self, symbol: str, build: str | None = None) -> GeneInfo | None:
        """Look up gene coordinates by symbol.

        Args:
            symbol: Gene symbol (e.g., "BRCA1", "TP53").
            build: Target genome build for the coordinates (overrides the client default).

        Returns:
            GeneInfo with 1-based coordinates (``.build`` states which build they are in), or
            None if no human gene has this official symbol (or a unique alias). If no
            location matches the requested build, coordinates for whatever assembly NCBI
            returned are given and ``.build`` reflects that — never silently the requested
            build; ``"unknown"`` when it cannot be determined.

        Raises:
            GeneLookupError: NCBI could not be queried (HTTP error, in-band error,
                unparseable body) — the gene's existence is unknown, not disproven.
        """
        target_build = normalize_genome_build(build or self.genome_build) or "unknown"

        try:
            search = check_esearch(
                await self._get(
                    NCBI_GENE_SEARCH,
                    {
                        "db": "gene",
                        "term": f"{symbol}[Gene Name] AND Homo sapiens[Organism]",
                        "retmode": "json",
                        "retmax": _MAX_CANDIDATES,
                    },
                    "Gene esearch",
                ),
                "Gene esearch",
            )
        except GeneLookupError:
            raise
        except NCBIUnavailableError as e:
            raise GeneLookupError(str(e)) from e

        id_list = [str(i) for i in search.get("idlist", [])]
        if not id_list:
            return None

        try:
            result = check_esummary(
                await self._get(
                    NCBI_GENE_SUMMARY,
                    {"db": "gene", "id": ",".join(id_list), "retmode": "json"},
                    "Gene esummary",
                ),
                "Gene esummary",
            )
        except GeneLookupError:
            raise
        except NCBIUnavailableError as e:
            raise GeneLookupError(str(e)) from e

        gene_data, note = _choose_gene(result, id_list, symbol)
        if gene_data is None:
            return None

        gene_symbol = gene_data.get("name", symbol)
        gene_name = gene_data.get("description", "")

        # Candidate locations, current annotation first. locationhist is newest-first and
        # (unlike genomicinfo) names the assembly; it is where GRCh37 coordinates live.
        locations = [
            loc
            for loc in list(gene_data.get("genomicinfo") or [])
            + list(gene_data.get("locationhist") or [])
            if isinstance(loc, dict) and loc.get("chraccver")
        ]
        if not locations:
            return None

        coords = next((loc for loc in locations if _location_build(loc) == target_build), None)
        warnings: list[str] = []
        if coords is None:
            coords = locations[0]
            actual = _location_build(coords)
            warnings.append(
                f"NCBI has no {target_build} location for {gene_symbol}; returning "
                f"{actual} coordinates"
            )
        actual_build = _location_build(coords)

        chrom = self._accession_to_chrom(str(coords.get("chraccver", "")))
        chrom_start = int(coords.get("chrstart", 0))
        chrom_stop = int(coords.get("chrstop", 0))

        # NCBI lists a minus-strand gene with chrstart > chrstop.
        strand = "+" if chrom_start <= chrom_stop else "-"

        # NCBI chrstart/chrstop are 0-based; report 1-based inclusive coordinates.
        start = min(chrom_start, chrom_stop) + 1
        end = max(chrom_start, chrom_stop) + 1

        return GeneInfo(
            symbol=gene_symbol,
            name=gene_name,
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
            build=actual_build,
            note=note,
            warnings=warnings,
        )

    def _accession_to_chrom(self, accession: str) -> str:
        """Convert RefSeq accession to chromosome name.

        Args:
            accession: RefSeq accession like "NC_000001.11"

        Returns:
            Chromosome name like "chr1" (the accession itself when not a human chromosome)
        """
        base_acc = accession.split(".")[0] if "." in accession else accession
        entry = _REFSEQ_CHROMOSOMES.get(base_acc)
        return entry[0] if entry else accession

    async def close(self) -> None:
        """Close the HTTP client."""
        await self._client.aclose()

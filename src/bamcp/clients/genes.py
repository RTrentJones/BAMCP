"""Gene annotation lookup via NCBI Entrez.

Provides gene symbol to genomic coordinate mapping using NCBI's Gene database.
"""

from __future__ import annotations

from dataclasses import dataclass

import httpx

NCBI_GENE_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_GENE_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


@dataclass
class GeneInfo:
    """Gene annotation with genomic coordinates."""

    symbol: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str


class GeneClient:
    """NCBI Gene database client for coordinate lookup.

    Uses NCBI E-utilities to search for human genes by symbol and
    retrieve their genomic coordinates for a specific genome build.
    """

    def __init__(
        self,
        api_key: str | None = None,
        genome_build: str = "GRCh38",
        timeout: float = 10.0,
    ):
        """Initialize the gene client.

        Args:
            api_key: NCBI API key for higher rate limits.
            genome_build: Target genome build (GRCh38 or GRCh37).
            timeout: HTTP request timeout in seconds.
        """
        self.api_key = api_key
        self.genome_build = genome_build
        self._client = httpx.AsyncClient(timeout=timeout)

    async def search(self, symbol: str) -> GeneInfo | None:
        """Look up gene coordinates by symbol.

        Args:
            symbol: Gene symbol (e.g., "BRCA1", "TP53").

        Returns:
            GeneInfo with coordinates, or None if not found.
        """
        # Step 1: Search for gene ID
        search_params: dict = {
            "db": "gene",
            "term": f"{symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmode": "json",
        }
        if self.api_key:
            search_params["api_key"] = self.api_key

        try:
            resp = await self._client.get(NCBI_GENE_SEARCH, params=search_params)
            resp.raise_for_status()
            data = resp.json()
        except (httpx.HTTPError, ValueError):
            return None

        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None

        gene_id = id_list[0]

        # Step 2: Get gene summary with coordinates
        summary_params: dict = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json",
        }
        if self.api_key:
            summary_params["api_key"] = self.api_key

        try:
            resp = await self._client.get(NCBI_GENE_SUMMARY, params=summary_params)
            resp.raise_for_status()
            summary = resp.json()
        except (httpx.HTTPError, ValueError):
            return None

        # Parse the summary result
        result = summary.get("result", {})
        gene_data = result.get(gene_id, {})

        if not gene_data:
            return None

        # Extract gene info
        gene_symbol = gene_data.get("name", symbol)
        gene_name = gene_data.get("description", "")

        # Get genomic coordinates from genomicinfo
        genomic_info = gene_data.get("genomicinfo", [])
        if not genomic_info:
            # Try locationhist as fallback
            location_hist = gene_data.get("locationhist", [])
            if location_hist:
                genomic_info = location_hist

        if not genomic_info:
            return None

        # Find the coordinate set for our target build
        coords = None
        for info in genomic_info:
            # Check if this matches our genome build
            assembly = info.get("assemblyaccver", "")
            if self.genome_build == "GRCh38" and "GCF_000001405" in assembly:
                coords = info
                break
            if self.genome_build == "GRCh37" and "GCF_000001405.25" in assembly:
                coords = info
                break

        # Fall back to first available if no exact match
        if coords is None and genomic_info:
            coords = genomic_info[0]

        if coords is None:
            return None

        # Extract chromosome and coordinates
        chrom_acc = coords.get("chraccver", "")
        chrom_start = coords.get("chrstart", 0)
        chrom_stop = coords.get("chrstop", 0)

        # Convert accession to chromosome name
        chrom = self._accession_to_chrom(chrom_acc)

        # Handle strand (based on start/stop ordering)
        strand = "+" if chrom_start <= chrom_stop else "-"

        # Ensure start < end
        start = min(chrom_start, chrom_stop)
        end = max(chrom_start, chrom_stop)

        return GeneInfo(
            symbol=gene_symbol,
            name=gene_name,
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
        )

    def _accession_to_chrom(self, accession: str) -> str:
        """Convert RefSeq accession to chromosome name.

        Args:
            accession: RefSeq accession like "NC_000001.11"

        Returns:
            Chromosome name like "chr1"
        """
        # Map NC accessions to chromosome names
        nc_to_chrom = {
            "NC_000001": "chr1",
            "NC_000002": "chr2",
            "NC_000003": "chr3",
            "NC_000004": "chr4",
            "NC_000005": "chr5",
            "NC_000006": "chr6",
            "NC_000007": "chr7",
            "NC_000008": "chr8",
            "NC_000009": "chr9",
            "NC_000010": "chr10",
            "NC_000011": "chr11",
            "NC_000012": "chr12",
            "NC_000013": "chr13",
            "NC_000014": "chr14",
            "NC_000015": "chr15",
            "NC_000016": "chr16",
            "NC_000017": "chr17",
            "NC_000018": "chr18",
            "NC_000019": "chr19",
            "NC_000020": "chr20",
            "NC_000021": "chr21",
            "NC_000022": "chr22",
            "NC_000023": "chrX",
            "NC_000024": "chrY",
            "NC_012920": "chrM",
        }

        # Extract the base accession (without version)
        base_acc = accession.split(".")[0] if "." in accession else accession

        return nc_to_chrom.get(base_acc, accession)

    async def close(self) -> None:
        """Close the HTTP client."""
        await self._client.aclose()

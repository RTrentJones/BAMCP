"""Unit tests for bamcp.genes module."""

import re

import pytest

from bamcp.genes import GeneClient, GeneInfo


class TestGeneInfo:
    """Tests for the GeneInfo dataclass."""

    @pytest.mark.unit
    def test_gene_info_creation(self):
        """Test creating a GeneInfo instance."""
        info = GeneInfo(
            symbol="BRCA1",
            name="BRCA1 DNA repair associated",
            chrom="chr17",
            start=43044295,
            end=43125483,
            strand="-",
        )
        assert info.symbol == "BRCA1"
        assert info.name == "BRCA1 DNA repair associated"
        assert info.chrom == "chr17"
        assert info.start == 43044295
        assert info.end == 43125483
        assert info.strand == "-"


class TestGeneClient:
    """Tests for the GeneClient class."""

    @pytest.mark.unit
    def test_init_without_api_key(self):
        """Test client initialization without API key."""
        client = GeneClient()
        assert client.api_key is None
        assert client.genome_build == "GRCh38"

    @pytest.mark.unit
    def test_init_with_api_key(self):
        """Test client initialization with API key."""
        client = GeneClient(api_key="test_key")
        assert client.api_key == "test_key"

    @pytest.mark.unit
    def test_init_with_grch37(self):
        """Test client initialization with GRCh37 build."""
        client = GeneClient(genome_build="GRCh37")
        assert client.genome_build == "GRCh37"

    @pytest.mark.unit
    def test_accession_to_chrom_chr1(self):
        """Test NC_000001 maps to chr1."""
        client = GeneClient()
        assert client._accession_to_chrom("NC_000001.11") == "chr1"

    @pytest.mark.unit
    def test_accession_to_chrom_chrX(self):
        """Test NC_000023 maps to chrX."""
        client = GeneClient()
        assert client._accession_to_chrom("NC_000023.11") == "chrX"

    @pytest.mark.unit
    def test_accession_to_chrom_chrY(self):
        """Test NC_000024 maps to chrY."""
        client = GeneClient()
        assert client._accession_to_chrom("NC_000024.10") == "chrY"

    @pytest.mark.unit
    def test_accession_to_chrom_chrM(self):
        """Test NC_012920 maps to chrM."""
        client = GeneClient()
        assert client._accession_to_chrom("NC_012920.1") == "chrM"

    @pytest.mark.unit
    def test_accession_to_chrom_unknown(self):
        """Test unknown accession returns as-is."""
        client = GeneClient()
        assert client._accession_to_chrom("NT_187654.1") == "NT_187654.1"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_brca1_found(self, httpx_mock):
        """Test successful gene search for BRCA1."""
        # Mock search response
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "idlist": ["672"],
                }
            },
        )

        # Mock summary response
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "672": {
                        "name": "BRCA1",
                        "description": "BRCA1 DNA repair associated",
                        "genomicinfo": [
                            {
                                "assemblyaccver": "GCF_000001405.40",
                                "chraccver": "NC_000017.11",
                                "chrstart": 43044295,
                                "chrstop": 43125483,
                            }
                        ],
                    }
                }
            },
        )

        client = GeneClient()
        result = await client.search("BRCA1")

        assert result is not None
        assert result.symbol == "BRCA1"
        assert result.name == "BRCA1 DNA repair associated"
        assert result.chrom == "chr17"
        assert result.start == 43044295
        assert result.end == 43125483

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_not_found(self, httpx_mock):
        """Test gene search for nonexistent gene."""
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "idlist": [],
                }
            },
        )

        client = GeneClient()
        result = await client.search("NOTAREALGENE123")

        assert result is None

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_http_error(self, httpx_mock):
        """Test gene search handles HTTP errors gracefully."""
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            status_code=500,
        )

        client = GeneClient()
        result = await client.search("BRCA1")

        assert result is None

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_no_genomic_info(self, httpx_mock):
        """Test gene search when no genomic info is available."""
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "idlist": ["12345"],
                }
            },
        )

        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "12345": {
                        "name": "TESTGENE",
                        "description": "Test gene",
                        "genomicinfo": [],
                    }
                }
            },
        )

        client = GeneClient()
        result = await client.search("TESTGENE")

        assert result is None

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_with_api_key(self, httpx_mock):
        """Test gene search includes API key in request."""
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"idlist": []}},
        )

        client = GeneClient(api_key="my_api_key")
        await client.search("BRCA1")

        # Check that api_key was included in the request
        request = httpx_mock.get_request()
        assert "api_key=my_api_key" in str(request.url)

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_reverse_strand(self, httpx_mock):
        """Test gene search detects reverse strand correctly."""
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"idlist": ["672"]}},
        )

        # chrstart > chrstop indicates reverse strand
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "672": {
                        "name": "BRCA1",
                        "description": "BRCA1 DNA repair associated",
                        "genomicinfo": [
                            {
                                "assemblyaccver": "GCF_000001405.40",
                                "chraccver": "NC_000017.11",
                                "chrstart": 43125483,  # Larger
                                "chrstop": 43044295,  # Smaller
                            }
                        ],
                    }
                }
            },
        )

        client = GeneClient()
        result = await client.search("BRCA1")

        assert result is not None
        assert result.strand == "-"
        # start should still be < end
        assert result.start == 43044295
        assert result.end == 43125483

        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_close_client(self):
        """Test closing the client."""
        client = GeneClient()
        await client.close()
        # Should not raise an error

"""Unit tests for bamcp.genes module.

Mocks mirror the live NCBI esummary shape: ``genomicinfo`` carries NO ``assemblyaccver``
(only chrloc / chraccver / chrstart / chrstop / exoncount), coordinates are 0-based, a
minus-strand gene has ``chrstart > chrstop``, and ``locationhist`` (newest first) names the
assembly and holds the GRCh37 location.
"""

import re

import pytest

from bamcp.clients.genes import (
    GeneClient,
    GeneInfo,
    GeneLookupError,
    _build_from_assembly,
    _build_from_chraccver,
)
from bamcp.clients.ncbi import TokenBucketRateLimiter

# TP53 (gene 7157), live shape trimmed. GRCh38 1-based: chr17:7668421-7687490 (minus strand).
TP53 = {
    "uid": "7157",
    "name": "TP53",
    "description": "tumor protein p53",
    "status": "",
    "currentid": "",
    "otheraliases": "BCC7, BMFS5, LFS1, P53, TRP53",
    "organism": {"scientificname": "Homo sapiens"},
    "genomicinfo": [
        {
            "chrloc": "17",
            "chraccver": "NC_000017.11",
            "chrstart": 7687489,
            "chrstop": 7668420,
            "exoncount": 13,
        }
    ],
    "locationhist": [
        {
            "annotationrelease": "RS_2025_08",
            "assemblyaccver": "GCF_000001405.40",
            "chraccver": "NC_000017.11",
            "chrstart": 7687489,
            "chrstop": 7668420,
        },
        {  # T2T-CHM13 — neither GRCh37 nor GRCh38
            "annotationrelease": "RS_2025_08",
            "assemblyaccver": "GCF_009914755.1",
            "chraccver": "NC_060941.1",
            "chrstart": 7591593,
            "chrstop": 7572543,
        },
        {
            "annotationrelease": "RS_2024_09",
            "assemblyaccver": "GCF_000001405.25",
            "chraccver": "NC_000017.10",
            "chrstart": 7590807,
            "chrstop": 7571738,
        },
        {
            "annotationrelease": "105",
            "assemblyaccver": "GCF_000001405.25",
            "chraccver": "NC_000017.10",
            "chrstart": 7590867,
            "chrstop": 7571719,
        },
    ],
}

# Two discontinued/replaced records "p53" also hits (live: ids 373859, 352997).
HCP5P3_REPLACED = {
    "uid": "373859",
    "name": "HCP5P3",
    "status": "1",
    "currentid": "352997",
    "otheraliases": "P5-03, P5-3, P5.3",
    "organism": {"scientificname": "Homo sapiens"},
    "genomicinfo": [],
    "locationhist": [],
}
HCP5P3 = {**HCP5P3_REPLACED, "uid": "352997", "status": "2", "currentid": ""}


def _mock(httpx_mock, *entries):
    ids = [e["uid"] for e in entries]
    httpx_mock.add_response(
        url=re.compile(r".*/esearch\.fcgi.*"),
        json={"esearchresult": {"count": str(len(ids)), "idlist": ids}},
    )
    httpx_mock.add_response(
        url=re.compile(r".*/esummary\.fcgi.*"),
        json={"result": {"uids": ids, **{e["uid"]: e for e in entries}}},
    )


class TestGeneInfo:
    """Tests for the GeneInfo dataclass."""

    @pytest.mark.unit
    def test_gene_info_creation(self):
        info = GeneInfo(
            symbol="BRCA1",
            name="BRCA1 DNA repair associated",
            chrom="chr17",
            start=43044295,
            end=43125483,
            strand="-",
        )
        assert info.symbol == "BRCA1"
        assert info.build == "unknown"
        assert info.warnings == []


class TestBuildInference:
    @pytest.mark.unit
    def test_build_from_chraccver(self):
        assert _build_from_chraccver("NC_000017.10") == "GRCh37"
        assert _build_from_chraccver("NC_000017.11") == "GRCh38"
        assert _build_from_chraccver("NC_000001.10") == "GRCh37"
        assert _build_from_chraccver("NC_000001.11") == "GRCh38"
        assert _build_from_chraccver("NC_000007.13") == "GRCh37"
        assert _build_from_chraccver("NC_000007.14") == "GRCh38"
        assert _build_from_chraccver("NC_000023.10") == "GRCh37"
        assert _build_from_chraccver("NC_000024.10") == "GRCh38"
        assert _build_from_chraccver("NC_012920.1") == "unknown"  # chrM shared by both
        assert _build_from_chraccver("NC_060941.1") == "unknown"  # T2T
        assert _build_from_chraccver("NC_000017.99") == "unknown"
        assert _build_from_chraccver("garbage") == "unknown"

    @pytest.mark.unit
    def test_build_from_assembly(self):
        assert _build_from_assembly("GCF_000001405.25") == "GRCh37"
        assert _build_from_assembly("GCF_000001405.13") == "GRCh37"
        assert _build_from_assembly("GCF_000001405.26") == "GRCh38"
        assert _build_from_assembly("GCF_000001405.40") == "GRCh38"
        assert _build_from_assembly("GCF_999999999.1") == "unknown"
        assert _build_from_assembly("garbage") == "unknown"


class TestGeneClient:
    """Tests for the GeneClient class."""

    @pytest.mark.unit
    def test_init_defaults(self):
        client = GeneClient()
        assert client.api_key is None
        assert client.genome_build == "GRCh38"
        assert GeneClient(genome_build="GRCh37").genome_build == "GRCh37"

    @pytest.mark.unit
    def test_shared_rate_limiter(self):
        limiter = TokenBucketRateLimiter(3)
        assert GeneClient(rate_limiter=limiter)._rate_limiter is limiter

    @pytest.mark.unit
    def test_accession_to_chrom(self):
        client = GeneClient()
        assert client._accession_to_chrom("NC_000001.11") == "chr1"
        assert client._accession_to_chrom("NC_000023.11") == "chrX"
        assert client._accession_to_chrom("NC_000024.10") == "chrY"
        assert client._accession_to_chrom("NC_012920.1") == "chrM"
        assert client._accession_to_chrom("NT_187654.1") == "NT_187654.1"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_grch38_build_inferred_without_assemblyaccver(self, httpx_mock):
        """genomicinfo has no assemblyaccver; the build comes from chraccver's version."""
        _mock(httpx_mock, TP53)
        client = GeneClient()
        result = await client.search("TP53")
        assert result is not None
        assert result.build == "GRCh38"
        assert result.chrom == "chr17"
        assert result.strand == "-"
        # NCBI 0-based → 1-based inclusive (TP53 GRCh38 is chr17:7,668,421-7,687,490)
        assert (result.start, result.end) == (7668421, 7687490)
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_grch37_from_locationhist(self, httpx_mock):
        _mock(httpx_mock, TP53)
        client = GeneClient()
        result = await client.search("TP53", build="hg19")
        assert result is not None
        assert result.build == "GRCh37"
        # newest GRCh37 annotation (RS_2024_09), 1-based
        assert (result.start, result.end) == (7571739, 7590808)
        assert result.warnings == []
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_missing_build_falls_back_honestly(self, httpx_mock):
        entry = {**TP53, "locationhist": []}
        _mock(httpx_mock, entry)
        client = GeneClient()
        result = await client.search("TP53", build="GRCh37")
        assert result is not None
        assert result.build == "GRCh38"  # honest: not silently claimed as GRCh37
        assert result.warnings and "no GRCh37 location" in result.warnings[0]
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_unknown_build_reported_as_unknown(self, httpx_mock):
        entry = {
            **TP53,
            "genomicinfo": [
                {"chraccver": "NC_000017.99", "chrstart": 1, "chrstop": 10, "exoncount": 1}
            ],
            "locationhist": [],
        }
        _mock(httpx_mock, entry)
        client = GeneClient()
        result = await client.search("TP53")
        assert result is not None
        assert result.build == "unknown"
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_official_symbol_chosen_not_idlist_first(self, httpx_mock):
        # The replaced record is first in idlist; the official-symbol match must win.
        _mock(httpx_mock, HCP5P3_REPLACED, TP53, HCP5P3)
        client = GeneClient()
        result = await client.search("tp53")
        assert result is not None and result.symbol == "TP53"
        assert result.note is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_alias_resolves_with_note(self, httpx_mock):
        # live: "p53" hits TP53 plus two replaced/discontinued HCP5P3 records
        _mock(httpx_mock, TP53, HCP5P3_REPLACED, HCP5P3)
        client = GeneClient()
        result = await client.search("p53")
        assert result is not None and result.symbol == "TP53"
        assert result.note and "alias" in result.note
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_ambiguous_alias_is_not_found(self, httpx_mock):
        other = {**TP53, "uid": "999", "name": "OTHER", "otheraliases": "P53"}
        _mock(httpx_mock, TP53 | {"name": "TP53X"}, other)
        client = GeneClient()
        assert await client.search("p53") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_not_found(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "count": "0",
                    "idlist": [],
                    "errorlist": {"phrasesnotfound": ["NOTAREALGENE123"], "fieldsnotfound": []},
                }
            },
        )
        client = GeneClient()
        assert await client.search("NOTAREALGENE123") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_http_error_is_unavailable_not_not_found(self, httpx_mock):
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), status_code=500)
        client = GeneClient()
        with pytest.raises(GeneLookupError):
            await client.search("BRCA1")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_in_band_error_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"ERROR": "Database is not supported"}},
        )
        client = GeneClient()
        with pytest.raises(GeneLookupError):
            await client.search("BRCA1")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_html_body_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), text="<html>oops</html>")
        client = GeneClient()
        with pytest.raises(GeneLookupError):
            await client.search("BRCA1")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_esummary_uid_error_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"), json={"esearchresult": {"idlist": ["7157"]}}
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "uids": ["7157"],
                    "7157": {"uid": "7157", "error": "cannot get document summary"},
                }
            },
        )
        client = GeneClient()
        with pytest.raises(GeneLookupError):
            await client.search("TP53")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_no_genomic_info(self, httpx_mock):
        _mock(httpx_mock, {**TP53, "genomicinfo": [], "locationhist": []})
        client = GeneClient()
        assert await client.search("TP53") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_with_api_key(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"count": "0", "idlist": []}},
        )
        client = GeneClient(api_key="my_api_key")
        await client.search("BRCA1")
        assert "api_key=my_api_key" in str(httpx_mock.get_request().url)
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_plus_strand(self, httpx_mock):
        entry = {
            **TP53,
            "name": "BRCA2",
            "genomicinfo": [
                {"chraccver": "NC_000013.11", "chrstart": 32315507, "chrstop": 32400267}
            ],
            "locationhist": [],
        }
        _mock(httpx_mock, entry)
        client = GeneClient()
        result = await client.search("BRCA2")
        assert result is not None
        assert result.strand == "+"
        assert (result.chrom, result.start, result.end) == ("chr13", 32315508, 32400268)
        await client.close()


class TestGeneNetwork:
    @pytest.mark.integration
    @pytest.mark.network
    @pytest.mark.asyncio
    async def test_real_tp53_both_builds(self):
        client = GeneClient(timeout=30.0)
        r38 = await client.search("TP53", build="GRCh38")
        r37 = await client.search("TP53", build="GRCh37")
        await client.close()
        assert r38 is not None and r38.build == "GRCh38" and r38.start == 7668421
        assert r37 is not None and r37.build == "GRCh37" and r37.chrom == "chr17"

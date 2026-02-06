"""Unit tests for bamcp.gnomad module."""

import httpx
import pytest

from bamcp.gnomad import (
    GnomadClient,
    GnomadResult,
    PopulationFrequency,
    _build_query,
    _build_variant_id,
    _parse_response,
)

# -- Sample API response fixtures -------------------------------------------

GNOMAD_RESPONSE = {
    "data": {
        "variant": {
            "variant_id": "17-7674220-G-A",
            "genome": {
                "ac": 5,
                "an": 152312,
                "homozygote_count": 0,
                "af": 3.28e-05,
                "populations": [
                    {
                        "id": "afr",
                        "ac": 2,
                        "an": 41406,
                        "homozygote_count": 0,
                        "af": 4.83e-05,
                    },
                    {
                        "id": "eur",
                        "ac": 3,
                        "an": 64472,
                        "homozygote_count": 0,
                        "af": 4.65e-05,
                    },
                ],
                "filters": ["PASS"],
            },
            "exome": None,
        }
    }
}

GNOMAD_NOT_FOUND = {
    "data": {
        "variant": None,
    }
}

GNOMAD_GRAPHQL_ERROR = {
    "errors": [
        {"message": "Invalid variant ID format"},
    ]
}

GNOMAD_EXOME_ONLY = {
    "data": {
        "variant": {
            "variant_id": "7-140453136-A-T",
            "genome": None,
            "exome": {
                "ac": 10,
                "an": 251478,
                "homozygote_count": 0,
                "af": 3.98e-05,
                "populations": [
                    {
                        "id": "nfe",
                        "ac": 8,
                        "an": 113740,
                        "homozygote_count": 0,
                        "af": 7.03e-05,
                    },
                ],
                "filters": [],
            },
        }
    }
}


class TestBuildVariantId:
    """Tests for _build_variant_id."""

    @pytest.mark.unit
    def test_basic_variant_id(self):
        vid = _build_variant_id("chr17", 7674220, "G", "A")
        assert vid == "17-7674220-G-A"

    @pytest.mark.unit
    def test_strips_chr_prefix(self):
        vid = _build_variant_id("chr1", 100, "A", "T")
        assert vid == "1-100-A-T"

    @pytest.mark.unit
    def test_numeric_chrom(self):
        vid = _build_variant_id("7", 140453136, "A", "T")
        assert vid == "7-140453136-A-T"

    @pytest.mark.unit
    def test_x_chromosome(self):
        vid = _build_variant_id("chrX", 100, "C", "G")
        assert vid == "X-100-C-G"


class TestBuildQuery:
    """Tests for _build_query."""

    @pytest.mark.unit
    def test_returns_graphql_query(self):
        query = _build_query()
        assert "variant" in query
        assert "variantId" in query
        assert "dataset" in query
        assert "genome" in query
        assert "exome" in query
        assert "populations" in query
        assert "af" in query


class TestParseResponse:
    """Tests for _parse_response."""

    @pytest.mark.unit
    def test_parse_genome_data(self):
        result = _parse_response(GNOMAD_RESPONSE, "17-7674220-G-A")
        assert result is not None
        assert isinstance(result, GnomadResult)
        assert result.variant_id == "17-7674220-G-A"
        assert result.global_af == pytest.approx(3.28e-05)
        assert result.ac == 5
        assert result.an == 152312
        assert result.homozygote_count == 0
        assert result.source == "genome"
        assert len(result.populations) == 2
        assert result.filters == ["PASS"]

    @pytest.mark.unit
    def test_parse_populations(self):
        result = _parse_response(GNOMAD_RESPONSE, "17-7674220-G-A")
        assert result is not None
        pop_ids = {p.id for p in result.populations}
        assert "afr" in pop_ids
        assert "eur" in pop_ids

        afr = next(p for p in result.populations if p.id == "afr")
        assert isinstance(afr, PopulationFrequency)
        assert afr.ac == 2
        assert afr.an == 41406
        assert afr.af == pytest.approx(4.83e-05)

    @pytest.mark.unit
    def test_parse_not_found(self):
        result = _parse_response(GNOMAD_NOT_FOUND, "99-1-A-T")
        assert result is None

    @pytest.mark.unit
    def test_parse_graphql_error(self):
        result = _parse_response(GNOMAD_GRAPHQL_ERROR, "bad-id")
        assert result is None

    @pytest.mark.unit
    def test_parse_exome_fallback(self):
        result = _parse_response(GNOMAD_EXOME_ONLY, "7-140453136-A-T")
        assert result is not None
        assert result.source == "exome"
        assert result.ac == 10
        assert result.an == 251478
        assert len(result.populations) == 1


class TestGnomadClient:
    """Tests for GnomadClient."""

    @pytest.mark.unit
    def test_init_defaults(self):
        client = GnomadClient()
        assert client.dataset == "gnomad_r4"
        assert client.timeout == 30.0

    @pytest.mark.unit
    def test_init_custom_dataset(self):
        client = GnomadClient(dataset="gnomad_r3")
        assert client.dataset == "gnomad_r3"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_found(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api",
            json=GNOMAD_RESPONSE,
        )

        client = GnomadClient()
        result = await client.lookup("chr17", 7674220, "G", "A")

        assert result is not None
        assert isinstance(result, GnomadResult)
        assert result.global_af == pytest.approx(3.28e-05)
        assert result.source == "genome"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_not_found(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api",
            json=GNOMAD_NOT_FOUND,
        )

        client = GnomadClient()
        result = await client.lookup("chr99", 1, "A", "T")
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_graphql_error(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api",
            json=GNOMAD_GRAPHQL_ERROR,
        )

        client = GnomadClient()
        result = await client.lookup("chr1", 100, "A", "T")
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_caches_result(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api",
            json=GNOMAD_RESPONSE,
        )

        client = GnomadClient()
        result1 = await client.lookup("chr17", 7674220, "G", "A")
        result2 = await client.lookup("chr17", 7674220, "G", "A")

        assert result1 is result2
        # Only one HTTP request should have been made
        assert len(httpx_mock.get_requests()) == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_timeout(self, httpx_mock):
        httpx_mock.add_exception(
            httpx.ReadTimeout("Connection timed out"),
            url="https://gnomad.broadinstitute.org/api",
        )

        client = GnomadClient(timeout=1.0)
        with pytest.raises(httpx.ReadTimeout):
            await client.lookup("chr17", 7674220, "G", "A")

"""Unit tests for bamcp.gnomad module."""

import httpx
import pytest

from bamcp.clients.gnomad import (
    GnomadBuildError,
    GnomadClient,
    GnomadQueryError,
    GnomadResult,
    PopulationFrequency,
    _build_variant_id,
    _parse_response,
    resolve_dataset,
)


class TestTransientErrorsNotCached:
    """A transient GraphQL error must not be cached and replayed as 'not found'."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_transient_error_propagates_and_is_retried_not_cached(self):
        client = GnomadClient()
        calls = {"n": 0}

        async def _boom(*_args):
            calls["n"] += 1
            raise GnomadQueryError("transient")

        client._do_lookup = _boom  # type: ignore[method-assign]
        # The error now PROPAGATES (so callers can report "unavailable" instead of a
        # misleading "not found") and is still not cached — the next call retries.
        with pytest.raises(GnomadQueryError):
            await client.lookup("chr1", 1, "A", "T")
        with pytest.raises(GnomadQueryError):
            await client.lookup("chr1", 1, "A", "T")
        assert calls["n"] == 2  # retried — the error was NOT cached as a negative


# -- Sample API response fixtures -------------------------------------------

GNOMAD_RESPONSE = {
    "data": {
        "variant": {
            "variant_id": "17-7674220-G-A",
            "genome": {
                "ac": 5,
                "an": 152312,
                "homozygote_count": 0,
                "faf95": {"popmax": 1.2e-05, "popmax_population": "afr"},
                # Real VariantPopulation objects carry no ``af`` (computed from ac/an).
                "populations": [
                    {"id": "afr", "ac": 2, "an": 41406, "homozygote_count": 0},
                    {"id": "nfe", "ac": 3, "an": 64472, "homozygote_count": 0},
                ],
                # gnomAD reports a passing site as an EMPTY filter list (never "PASS").
                "filters": [],
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
                "populations": [
                    {"id": "nfe", "ac": 8, "an": 113740, "homozygote_count": 0},
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


class TestParseResponse:
    """Tests for _parse_response."""

    @pytest.mark.unit
    def test_parse_genome_data(self):
        result = _parse_response(GNOMAD_RESPONSE, "17-7674220-G-A")
        assert result is not None
        assert isinstance(result, GnomadResult)
        assert result.variant_id == "17-7674220-G-A"
        assert result.global_af == pytest.approx(5 / 152312)
        assert result.ac == 5
        assert result.an == 152312
        assert result.homozygote_count == 0
        assert result.source == "genome"
        assert len(result.populations) == 2
        assert result.filters == []
        assert result.faf95_popmax == pytest.approx(1.2e-05)
        assert result.frequency_reliable is True

    @pytest.mark.unit
    def test_parse_populations(self):
        result = _parse_response(GNOMAD_RESPONSE, "17-7674220-G-A")
        assert result is not None
        pop_ids = {p.id for p in result.populations}
        assert "afr" in pop_ids
        assert "nfe" in pop_ids

        afr = next(p for p in result.populations if p.id == "afr")
        assert isinstance(afr, PopulationFrequency)
        assert afr.ac == 2
        assert afr.an == 41406
        # af is calculated from ac/an (not returned by gnomAD API for populations)
        assert afr.af == pytest.approx(2 / 41406)

    @pytest.mark.unit
    def test_parse_not_found(self):
        result = _parse_response(GNOMAD_NOT_FOUND, "99-1-A-T")
        assert result is None

    @pytest.mark.unit
    def test_parse_graphql_error(self):
        # A GraphQL error is transient, not a genuine not-found — it raises so the caller
        # does not cache and replay it as absence.
        with pytest.raises(GnomadQueryError):
            _parse_response(GNOMAD_GRAPHQL_ERROR, "bad-id")

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
        assert client.dataset == "auto"  # per-build: gnomad_r2_1 / gnomad_r4
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
        assert result.global_af == pytest.approx(5 / 152312)
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
        # A GraphQL error is transient — it propagates (not swallowed to None) so the
        # handler reports "unavailable" rather than a misleading "not found".
        with pytest.raises(GnomadQueryError):
            await client.lookup("chr1", 100, "A", "T")

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


# Live shape for a well-formed id absent from the dataset (HTTP 200).
GNOMAD_VARIANT_NOT_FOUND = {"errors": [{"message": "Variant not found"}], "data": {"variant": None}}


def _pops(**groups):
    """``populations`` list from ``id=(ac, an)`` pairs (live shape: no per-pop af)."""
    return [
        {"id": pid.replace("__", ":"), "ac": ac, "an": an, "homozygote_count": 0}
        for pid, (ac, an) in groups.items()
    ]


def _block(ac, an, pops, filters=(), faf=None):
    return {
        "ac": ac,
        "an": an,
        "homozygote_count": 0,
        "filters": list(filters),
        "faf95": {"popmax": faf, "popmax_population": "nfe"} if faf is not None else None,
        "populations": pops,
    }


def _variant(genome=None, exome=None, joint=None, vid="1-100-A-T"):
    return {
        "data": {"variant": {"variant_id": vid, "genome": genome, "exome": exome, "joint": joint}}
    }


class TestVariantNotFound:
    @pytest.mark.unit
    def test_variant_not_found_error_is_genuine_not_found(self):
        assert _parse_response(GNOMAD_VARIANT_NOT_FOUND, "17-7676154-G-A") is None

    @pytest.mark.unit
    def test_other_errors_stay_unavailable(self):
        with pytest.raises(GnomadQueryError):
            _parse_response(
                {"errors": [{"message": "Invalid variant ID"}], "data": {"variant": None}}, "x"
            )
        mixed = {
            "errors": [{"message": "Variant not found"}, {"message": "Internal error"}],
            "data": {"variant": None},
        }
        with pytest.raises(GnomadQueryError):
            _parse_response(mixed, "x")

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_not_found_is_cached(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api", json=GNOMAD_VARIANT_NOT_FOUND
        )
        client = GnomadClient()
        assert await client.lookup("17", 7676154, "G", "A") is None
        assert await client.lookup("17", 7676154, "G", "A") is None
        assert len(httpx_mock.get_requests()) == 1
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_html_body_is_query_error(self, httpx_mock):
        httpx_mock.add_response(
            url="https://gnomad.broadinstitute.org/api", text="<html>502</html>"
        )
        client = GnomadClient()
        with pytest.raises(GnomadQueryError):
            await client.lookup("17", 7676154, "G", "A")
        await client.close()


class TestPopulationMax:
    @pytest.mark.unit
    def test_single_allele_in_tiny_cohort_is_not_popmax(self):
        """hgdp:adygei_XY with AN=14 and one allele must not become max_pop_af 0.071 (BA1)."""
        pops = _pops(
            nfe=(1, 60000),
            afr=(0, 40000),
            nfe_XY=(1, 30000),
            hgdp__adygei=(1, 30),
            hgdp__adygei_XY=(1, 14),
            **{"1kg__yri": (0, 200)},
        )
        r = _parse_response(_variant(genome=_block(1, 150000, pops)), "1-100-A-T")
        assert r is not None
        assert r.max_pop == "nfe"
        assert r.max_pop_af == pytest.approx(1 / 60000)
        assert r.max_pop_an == 60000
        # only top-level groups are reported
        assert {p.id for p in r.populations} == {"nfe", "afr"}

    @pytest.mark.unit
    def test_founder_and_low_an_groups_excluded(self):
        pops = _pops(
            asj=(50, 10000),
            fin=(40, 20000),
            mid=(10, 3000),
            ami=(5, 900),
            remaining=(9, 5000),
            sas=(3, 1500),
            eas=(1, 20000),
        )
        r = _parse_response(_variant(genome=_block(118, 150000, pops)), "1-100-A-T")
        assert r is not None
        assert r.max_pop == "eas"  # sas has AN < 2000; founders/remaining are ineligible
        assert {p.id for p in r.populations} >= {"asj", "fin", "mid", "ami", "remaining"}

    @pytest.mark.unit
    def test_no_eligible_group_gives_none(self):
        r = _parse_response(_variant(genome=_block(1, 1000, _pops(nfe=(1, 1000)))), "x")
        assert r is not None
        assert r.max_pop_af is None and r.max_pop is None
        assert r.frequency_reliable is False
        assert any("low allele number" in w for w in r.warnings)

    @pytest.mark.unit
    def test_an_zero_is_unknown_not_zero(self):
        r = _parse_response(_variant(genome=_block(0, 0, _pops(nfe=(0, 0)))), "x")
        assert r is not None
        assert r.global_af is None
        assert r.populations[0].af is None
        assert r.frequency_reliable is False


class TestFrequencySources:
    @pytest.mark.unit
    def test_joint_preferred(self):
        joint = _block(14, 1613594, _pops(nfe=(14, 1200000)), faf=5e-06)
        exome = _block(14, 1461570, _pops(nfe=(14, 1100000)), faf=5.31e-06)
        r = _parse_response(_variant(exome=exome, joint=joint), "x")
        assert r is not None
        assert r.source == "joint"
        assert r.an == 1613594
        assert r.faf95_popmax == pytest.approx(5e-06)

    @pytest.mark.unit
    def test_exome_only_site_in_r4_still_counts_exome(self):
        """r4 genome is null for an exome-only site; joint carries the exome alleles."""
        joint = _block(1, 1613594, _pops(nfe=(1, 1200000)))
        exome = _block(1, 1461570, _pops(nfe=(1, 1100000)))
        r = _parse_response(_variant(genome=None, exome=exome, joint=joint), "x")
        assert r is not None and r.ac == 1

    @pytest.mark.unit
    def test_genome_and_exome_combined_without_joint(self):
        genome = _block(10, 30000, _pops(nfe=(10, 15000)), faf=1e-4)
        exome = _block(90, 250000, _pops(nfe=(90, 110000)), faf=5e-4)
        r = _parse_response(_variant(genome=genome, exome=exome), "x")
        assert r is not None
        assert r.source == "genome+exome"
        assert (r.ac, r.an) == (100, 280000)
        nfe = r.populations[0]
        assert (nfe.ac, nfe.an) == (100, 125000)
        assert r.faf95_popmax == pytest.approx(5e-4)  # the larger-AN block's faf

    @pytest.mark.unit
    def test_filtered_block_excluded_and_surfaced(self):
        genome = _block(10, 30000, _pops(nfe=(10, 15000)))
        exome = _block(900, 250000, _pops(nfe=(900, 110000)), filters=["RF"])
        r = _parse_response(_variant(genome=genome, exome=exome), "x")
        assert r is not None
        assert r.source == "genome"
        assert r.ac == 10
        assert r.filters == ["exome:RF"]
        assert any("exome data excluded" in w for w in r.warnings)

    @pytest.mark.unit
    def test_all_sources_filtered_is_not_reliable(self):
        genome = _block(1, 150000, _pops(nfe=(1, 60000)), filters=["AC0"])
        r = _parse_response(_variant(genome=genome), "x")
        assert r is not None
        assert r.filters == ["genome:AC0"]
        assert r.frequency_reliable is False


class TestBuildRouting:
    @pytest.mark.unit
    def test_resolve_dataset(self):
        assert resolve_dataset("auto", "GRCh37") == ("gnomad_r2_1", "GRCh37")
        assert resolve_dataset("auto", "hg19") == ("gnomad_r2_1", "GRCh37")
        assert resolve_dataset("auto", "hg38") == ("gnomad_r4", "GRCh38")
        assert resolve_dataset("gnomad_r3", "GRCh38") == ("gnomad_r3", "GRCh38")

    @pytest.mark.unit
    def test_unknown_build_refused(self):
        with pytest.raises(GnomadBuildError):
            resolve_dataset("auto", "unknown")
        with pytest.raises(GnomadBuildError):
            resolve_dataset("auto", None)

    @pytest.mark.unit
    def test_explicit_dataset_in_other_build_refused(self):
        with pytest.raises(GnomadBuildError, match="GRCh38"):
            resolve_dataset("gnomad_r4", "GRCh37")

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_grch37_lookup_queries_r2_1(self, httpx_mock):
        import json

        httpx_mock.add_response(url="https://gnomad.broadinstitute.org/api", json=GNOMAD_RESPONSE)
        client = GnomadClient()
        result = await client.lookup("17", 7579472, "G", "C", build="GRCh37")
        sent = json.loads(httpx_mock.get_requests()[0].content)
        assert sent["variables"]["dataset"] == "gnomad_r2_1"
        assert result is not None
        assert (result.dataset, result.build) == ("gnomad_r2_1", "GRCh37")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_unknown_build_lookup_raises_before_network(self):
        client = GnomadClient()
        with pytest.raises(GnomadBuildError):
            await client.lookup("17", 1, "G", "C", build="unknown")
        await client.close()


async def _gnomad_lookup_with_retry(client, *args, retries=3, backoff=2.0, **kwargs):
    """Retry gnomAD lookups with exponential backoff for CI resilience.

    Waits at least 2s between attempts to stay under 0.5 req/s.
    """
    import asyncio

    for attempt in range(retries):
        try:
            return await client.lookup(*args, **kwargs)
        except (httpx.ReadTimeout, httpx.ConnectTimeout):
            if attempt == retries - 1:
                raise
            await asyncio.sleep(backoff * (attempt + 1))


class TestGnomadIntegration:
    """Integration tests that hit the real gnomAD API.

    Uses retry with backoff and inter-test delays to stay under 0.5 req/s.
    """

    @pytest.mark.integration
    @pytest.mark.network
    @pytest.mark.asyncio
    async def test_real_api_known_variant(self):
        """Test lookup of a known common variant (rs1042522 in TP53).

        This is a well-characterized SNP (Pro72Arg) in TP53 with high
        population frequency, so it should always be present in gnomAD.
        GRCh38 coordinates: chr17:7676154 G>C
        """
        client = GnomadClient(dataset="gnomad_r4", timeout=60.0)
        result = await _gnomad_lookup_with_retry(client, "chr17", 7676154, "G", "C")

        # Variant should be found
        assert result is not None, "Known variant rs1042522 not found in gnomAD"
        assert isinstance(result, GnomadResult)

        # Should have frequency data
        assert result.global_af > 0, "Expected non-zero allele frequency"
        assert result.ac > 0, "Expected non-zero allele count"
        assert result.an > 0, "Expected non-zero allele number"

        # Should have population data
        assert len(result.populations) > 0, "Expected population frequency data"

        # Population af should be calculated correctly
        for pop in result.populations:
            if pop.an > 0:
                expected_af = pop.ac / pop.an
                assert pop.af == pytest.approx(expected_af), (
                    f"Population {pop.id} af mismatch: {pop.af} != {expected_af}"
                )

    @pytest.mark.integration
    @pytest.mark.network
    @pytest.mark.asyncio
    async def test_real_api_not_found(self):
        """Test lookup of a variant that shouldn't exist."""
        import asyncio

        await asyncio.sleep(2.0)  # Rate limit: stay under 0.5 req/s
        client = GnomadClient(dataset="gnomad_r4", timeout=60.0)
        # Made-up variant at a position unlikely to have this specific change
        result = await _gnomad_lookup_with_retry(client, "chr1", 1, "A", "T")

        # Should return None for non-existent variant
        assert result is None

    @pytest.mark.integration
    @pytest.mark.network
    @pytest.mark.asyncio
    async def test_real_api_grch37_routes_to_r2_1(self):
        """rs1042522 at its GRCh37 coordinate (17:7579472) is found via gnomad_r2_1."""
        import asyncio

        await asyncio.sleep(2.0)
        client = GnomadClient(timeout=60.0)
        result = await _gnomad_lookup_with_retry(client, "17", 7579472, "G", "C", build="GRCh37")
        assert result is not None and result.dataset == "gnomad_r2_1"
        assert result.max_pop in {"afr", "amr", "eas", "nfe", "sas"}

    @pytest.mark.integration
    @pytest.mark.network
    @pytest.mark.asyncio
    async def test_real_api_variant_not_found_is_none(self):
        """The live "Variant not found" error is a genuine not-found, not unavailable."""
        import asyncio

        await asyncio.sleep(2.0)
        client = GnomadClient(timeout=60.0)
        result = await _gnomad_lookup_with_retry(client, "17", 7676154, "G", "A", build="GRCh37")
        assert result is None

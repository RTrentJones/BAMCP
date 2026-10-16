"""Unit tests for bamcp.clinvar module."""

import re

import httpx
import pytest

from bamcp.clients.clinvar import (
    REVIEW_STARS,
    ClinVarClient,
    ClinVarResult,
    ClinVarUnavailableError,
    TokenBucketRateLimiter,
    _build_search_term,
    _parse_entry,
    _select_entry,
    _span_end,
    _spdi_matches,
)
from bamcp.clients.ncbi import NCBIUnavailableError

# -- Sample API response fixtures -------------------------------------------

ESEARCH_RESPONSE = {
    "esearchresult": {
        "count": "1",
        "retmax": "20",
        "idlist": ["12345"],
    }
}

ESEARCH_EMPTY = {
    "esearchresult": {
        "count": "0",
        "retmax": "20",
        "idlist": [],
    }
}


def _entry(
    uid,
    spdi,
    *,
    description="Pathogenic",
    review="criteria provided, multiple submitters, no conflicts",
    gene="TP53",
    traits=("Li-Fraumeni syndrome", "Hereditary cancer-predisposing syndrome"),
    title="NM_000546.6(TP53):c.743G>A (p.Arg248Gln)",
    loc38=None,
    loc37=None,
    oncogenicity="",
):
    """Build an esummary entry in ClinVar's current (germline_classification) schema.

    Mirrors the live esummary shape: the somatic blocks are always present (empty
    descriptions when unclassified) and ``variation_loc`` lists the record's span on both
    assemblies (``(chr, start, stop)``; 1-based, as ClinVar reports them).
    """
    locs = []
    for assembly, loc in (("GRCh38", loc38), ("GRCh37", loc37)):
        if loc is not None:
            chrom, start, stop = loc
            locs.append(
                {
                    "status": "current" if assembly == "GRCh38" else "previous",
                    "assembly_name": assembly,
                    "chr": chrom,
                    "start": str(start),
                    "stop": str(stop),
                    "ref": "",
                    "alt": "",
                }
            )
    return {
        "uid": uid,
        "title": title,
        "germline_classification": {
            "description": description,
            "review_status": review,
            "trait_set": [{"trait_name": t} for t in traits],
            "last_evaluated": "2023-01-15",
        },
        "clinical_impact_classification": {
            "description": "",
            "last_evaluated": "1/01/01 00:00",
            "review_status": "",
        },
        "oncogenicity_classification": {
            "description": oncogenicity,
            "last_evaluated": "1/01/01 00:00",
            "review_status": "criteria provided, single submitter" if oncogenicity else "",
        },
        "genes": [{"symbol": gene}] if gene else [],
        "variation_set": [{"canonical_spdi": spdi, "variant_type": "snv", "variation_loc": locs}],
    }


# TP53 c.743G>A at GRCh38 chr17:7674220 → SPDI 0-based 7674219.
ESUMMARY_RESPONSE = {
    "result": {"uids": ["12345"], "12345": _entry("12345", "NC_000017.11:7674219:G:A")}
}

# A position search returning an overlapping non-matching record plus the exact SNV.
ESUMMARY_MULTI = {
    "result": {
        "uids": ["111", "222"],
        "111": _entry("111", "NC_000017.11:7674218:GG:AA", description="Benign", title="delins"),
        "222": _entry("222", "NC_000017.11:7674219:G:A", description="Pathogenic"),
    }
}

ESUMMARY_EMPTY = {
    "result": {
        "uids": [],
    }
}


class TestBuildSearchTerm:
    """Tests for _build_search_term (position-only; allele resolved from SPDI)."""

    @pytest.mark.unit
    def test_basic_search_term(self):
        term = _build_search_term("chr17", 7674220)
        assert term == "17[Chromosome] AND 7674220[chrpos38]"

    @pytest.mark.unit
    def test_no_variant_name_clause(self):
        # ClinVar's [Variant name] is HGVS, not REF>ALT; a ref>alt clause matches nothing.
        assert "[Variant name]" not in _build_search_term("17", 100)

    @pytest.mark.unit
    def test_strips_chr_prefix(self):
        assert _build_search_term("chr1", 100).startswith("1[")

    @pytest.mark.unit
    def test_numeric_chrom(self):
        assert _build_search_term("17", 100).startswith("17[")


class TestSpdiMatching:
    """Tests for _spdi_matches (canonical SPDI seq:0based_pos:del:ins)."""

    @pytest.mark.unit
    def test_snv_matches_at_1based_pos(self):
        assert _spdi_matches("NC_000017.11:7674219:G:A", 7674220, "G", "A")

    @pytest.mark.unit
    def test_position_off_by_one_does_not_match(self):
        assert not _spdi_matches("NC_000017.11:7674219:G:A", 7674219, "G", "A")

    @pytest.mark.unit
    def test_allele_mismatch_does_not_match(self):
        assert not _spdi_matches("NC_000017.11:7674219:G:A", 7674220, "G", "T")

    @pytest.mark.unit
    def test_case_insensitive(self):
        assert _spdi_matches("NC_000017.11:7674219:g:a", 7674220, "G", "A")

    @pytest.mark.unit
    def test_vcf_insertion_matches_canonical_spdi(self):
        # VCF pos=100 ref=A alt=AT (leading anchor) ↔ canonical SPDI 100::T.
        assert _spdi_matches("NC_000001.11:100::T", 100, "A", "AT")

    @pytest.mark.unit
    def test_vcf_deletion_matches_canonical_spdi(self):
        # VCF pos=100 ref=AT alt=A ↔ canonical SPDI 100:T:.
        assert _spdi_matches("NC_000001.11:100:T:", 100, "AT", "A")

    @pytest.mark.unit
    def test_insertion_wrong_base_does_not_match(self):
        assert not _spdi_matches("NC_000001.11:100::G", 100, "A", "AT")

    @pytest.mark.unit
    def test_malformed_spdi(self):
        assert not _spdi_matches("not-an-spdi", 7674220, "G", "A")
        assert not _spdi_matches("NC:x:G:A", 7674220, "G", "A")


class TestSelectEntry:
    """Tests for _select_entry (disambiguate the exact allele among position hits)."""

    @pytest.mark.unit
    def test_selects_allele_matching_record(self):
        entry = _select_entry(ESUMMARY_MULTI["result"], 7674220, "G", "A")
        assert entry is not None
        assert entry["uid"] == "222"  # the exact SNV, not the overlapping delins

    @pytest.mark.unit
    def test_none_when_no_allele_matches(self):
        assert _select_entry(ESUMMARY_MULTI["result"], 7674220, "G", "C") is None

    @pytest.mark.unit
    def test_none_for_empty_result(self):
        assert _select_entry(ESUMMARY_EMPTY["result"], 7674220, "G", "A") is None

    @pytest.mark.unit
    def test_haplotype_record_not_matched(self):
        # A multi-variant (haplotype/genotype) record must not match a single-allele
        # lookup, even when one component SPDI matches.
        result = {
            "uids": ["h1"],
            "h1": {
                "uid": "h1",
                "variation_set": [
                    {"canonical_spdi": "NC_000017.11:7674219:G:A"},
                    {"canonical_spdi": "NC_000017.11:7674300:C:T"},
                ],
            },
        }
        assert _select_entry(result, 7674220, "G", "A") is None


class TestParseEntry:
    """Tests for _parse_entry (current + legacy esummary schema)."""

    @pytest.mark.unit
    def test_parse_germline_classification(self):
        result = _parse_entry(_entry("12345", "NC_000017.11:7674219:G:A"))
        assert result is not None
        assert result.variation_id == 12345
        assert result.clinical_significance == "Pathogenic"
        assert result.review_status == "criteria provided, multiple submitters, no conflicts"
        assert result.stars == 2
        assert "Li-Fraumeni syndrome" in result.conditions
        assert result.gene == "TP53"
        assert result.variant_name == "NM_000546.6(TP53):c.743G>A (p.Arg248Gln)"

    @pytest.mark.unit
    def test_parse_no_uid(self):
        assert _parse_entry({"title": "x"}) is None

    @pytest.mark.unit
    def test_parse_legacy_schema(self):
        # Pre-2024 top-level clinical_significance / review_status / trait_set.
        entry = {
            "uid": "999",
            "clinical_significance": {"description": "Pathogenic"},
            "review_status": "reviewed by expert panel",
            "genes": [{"symbol": "BRCA1"}],
            "trait_set": [{"trait_name": "Breast cancer"}],
        }
        result = _parse_entry(entry)
        assert result is not None
        assert result.clinical_significance == "Pathogenic"
        assert result.stars == 3
        assert result.conditions == ["Breast cancer"]

    @pytest.mark.unit
    def test_parse_significance_as_string(self):
        entry = {
            "uid": "999",
            "clinical_significance": "Benign",
            "review_status": "no assertion criteria provided",
            "genes": [],
            "trait_set": [],
        }
        result = _parse_entry(entry)
        assert result is not None
        assert result.clinical_significance == "Benign"
        assert result.stars == 0

    @pytest.mark.unit
    def test_parse_gene_as_string(self):
        entry = {
            "uid": "100",
            "germline_classification": {"description": "VUS", "review_status": "x"},
            "genes": ["BRCA1"],
        }
        result = _parse_entry(entry)
        assert result is not None
        assert result.gene == "BRCA1"


class TestReviewStars:
    """Tests for review status to stars mapping."""

    @pytest.mark.unit
    def test_practice_guideline(self):
        assert REVIEW_STARS["practice guideline"] == 4

    @pytest.mark.unit
    def test_expert_panel(self):
        assert REVIEW_STARS["reviewed by expert panel"] == 3

    @pytest.mark.unit
    def test_multiple_submitters(self):
        assert REVIEW_STARS["criteria provided, multiple submitters, no conflicts"] == 2

    @pytest.mark.unit
    def test_no_criteria(self):
        assert REVIEW_STARS["no assertion criteria provided"] == 0


class TestClinVarClient:
    """Tests for ClinVarClient."""

    @pytest.mark.unit
    def test_init_without_api_key(self):
        client = ClinVarClient()
        assert client.api_key is None
        assert client._semaphore._value == 3  # type: ignore[attr-defined]

    @pytest.mark.unit
    def test_init_with_api_key(self):
        client = ClinVarClient(api_key="test_key")
        assert client.api_key == "test_key"
        assert client._semaphore._value == 10  # type: ignore[attr-defined]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_found(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json=ESEARCH_RESPONSE,
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json=ESUMMARY_RESPONSE,
        )

        client = ClinVarClient()
        result = await client.lookup("chr17", 7674220, "G", "A")

        assert result is not None
        assert isinstance(result, ClinVarResult)
        assert result.clinical_significance == "Pathogenic"
        assert result.gene == "TP53"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_not_found(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json=ESEARCH_EMPTY,
        )

        client = ClinVarClient()
        result = await client.lookup("chr99", 1, "A", "T")
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_selects_matching_allele(self, httpx_mock):
        # Position search returns an overlapping delins + the exact SNV; pick the SNV.
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), json=ESEARCH_RESPONSE)
        httpx_mock.add_response(url=re.compile(r".*/esummary\.fcgi.*"), json=ESUMMARY_MULTI)

        client = ClinVarClient()
        result = await client.lookup("chr17", 7674220, "G", "A")
        assert result is not None
        assert result.variation_id == 222
        assert result.clinical_significance == "Pathogenic"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_ids_pages_past_limit(self, httpx_mock):
        # count=3 across two pages; _search_ids must collect every ID, not just page 1.
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"count": "3", "retstart": "0", "idlist": ["1", "2"]}},
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"count": "3", "retstart": "2", "idlist": ["3"]}},
        )
        client = ClinVarClient()
        ids, truncated = await client._search_ids("chr17", 7674220)
        assert ids == ["1", "2", "3"]
        assert truncated is False

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_allele_mismatch_returns_none(self, httpx_mock):
        # Records exist at the position but none match the requested allele.
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), json=ESEARCH_RESPONSE)
        httpx_mock.add_response(url=re.compile(r".*/esummary\.fcgi.*"), json=ESUMMARY_RESPONSE)

        client = ClinVarClient()
        result = await client.lookup("chr17", 7674220, "G", "C")  # G>C, but record is G>A
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_caches_result(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json=ESEARCH_RESPONSE,
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json=ESUMMARY_RESPONSE,
        )

        client = ClinVarClient()
        result1 = await client.lookup("chr17", 7674220, "G", "A")
        result2 = await client.lookup("chr17", 7674220, "G", "A")

        assert result1 is result2
        # Only 2 requests should have been made (esearch + esummary), not 4
        assert len(httpx_mock.get_requests()) == 2

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_lookup_timeout(self, httpx_mock):
        httpx_mock.add_exception(
            httpx.ReadTimeout("Connection timed out"),
            url=re.compile(r".*/esearch\.fcgi.*"),
        )

        client = ClinVarClient(timeout=1.0)
        with pytest.raises(httpx.ReadTimeout):
            await client.lookup("chr17", 7674220, "G", "A")


class TestBoundedTTLCache:
    """Tests for BoundedTTLCache behavior."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_cache_basic_set_get(self):
        from bamcp.clients.clinvar import BoundedTTLCache

        cache: BoundedTTLCache[str] = BoundedTTLCache(maxsize=100, ttl=3600)
        await cache.set(("key1",), "value1")
        result = await cache.get(("key1",))
        assert result == "value1"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_cache_miss_returns_none(self):
        from bamcp.clients.clinvar import BoundedTTLCache

        cache: BoundedTTLCache[str] = BoundedTTLCache(maxsize=100, ttl=3600)
        result = await cache.get(("nonexistent",))
        assert result is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_get_entry_distinguishes_miss_from_cached_none(self):
        """A cached None (negative result) reads back as a hit, not a miss."""
        from bamcp.clients.ttl_cache import MISS, BoundedTTLCache

        cache: BoundedTTLCache = BoundedTTLCache(maxsize=100, ttl=3600)
        # Miss → sentinel.
        assert await cache.get_entry(("absent",)) is MISS
        # Cache a negative result and read it back as a hit.
        await cache.set(("absent",), None)
        assert await cache.get_entry(("absent",)) is None
        assert await cache.get_entry(("absent",)) is not MISS

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_negative_result_is_not_refetched(self):
        """A cached 'not found' is served from cache instead of re-hitting the API."""
        from bamcp.clients.gnomad import GnomadClient

        client = GnomadClient()
        calls = {"n": 0}

        async def _fake_lookup(*_args):
            calls["n"] += 1
            return None  # variant not found

        client._do_lookup = _fake_lookup  # type: ignore[method-assign]
        assert await client.lookup("chr1", 1, "A", "T") is None
        assert await client.lookup("chr1", 1, "A", "T") is None
        assert calls["n"] == 1  # second lookup served from the negative cache

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_cache_evicts_oldest_when_full(self):
        from bamcp.clients.clinvar import BoundedTTLCache

        cache: BoundedTTLCache[str] = BoundedTTLCache(maxsize=3, ttl=3600)
        await cache.set(("key1",), "value1")
        await cache.set(("key2",), "value2")
        await cache.set(("key3",), "value3")
        # Cache is now full

        # Add a fourth item - should evict key1 (oldest)
        await cache.set(("key4",), "value4")

        assert await cache.get(("key1",)) is None  # Evicted
        assert await cache.get(("key2",)) == "value2"
        assert await cache.get(("key3",)) == "value3"
        assert await cache.get(("key4",)) == "value4"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_cache_lru_access_updates_order(self):
        from bamcp.clients.clinvar import BoundedTTLCache

        cache: BoundedTTLCache[str] = BoundedTTLCache(maxsize=3, ttl=3600)
        await cache.set(("key1",), "value1")
        await cache.set(("key2",), "value2")
        await cache.set(("key3",), "value3")

        # Access key1 to make it most recently used
        await cache.get(("key1",))

        # Add a fourth item - should evict key2 (now oldest)
        await cache.set(("key4",), "value4")

        assert await cache.get(("key1",)) == "value1"  # Still there (was accessed)
        assert await cache.get(("key2",)) is None  # Evicted
        assert await cache.get(("key3",)) == "value3"
        assert await cache.get(("key4",)) == "value4"


class TestTokenBucketRateLimiter:
    """Tests for TokenBucketRateLimiter behavior."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_rate_limiter_allows_initial_requests(self):
        from bamcp.clients.clinvar import TokenBucketRateLimiter

        limiter = TokenBucketRateLimiter(rate=5.0)  # 5 requests per second
        # Should be able to make 5 requests immediately
        for _ in range(5):
            await limiter.acquire()
        # Test passes if no exception

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_rate_limiter_initialized_with_rate(self):
        from bamcp.clients.clinvar import TokenBucketRateLimiter

        limiter = TokenBucketRateLimiter(rate=10.0)
        assert limiter._rate == 10.0


# -- Build-aware lookups, honest classifications, in-band errors --------------

# rs1042522 (TP53 c.215C>G): GRCh38 chr17:7676154, GRCh37 chr17:7579472 — shapes from the
# live esummary (SPDI is GRCh38-only; variation_loc carries both assemblies).
RS1042522 = _entry(
    "12351",
    "NC_000017.11:7676153:G:C",
    description="Benign",
    review="reviewed by expert panel",
    loc38=("17", 7676154, 7676154),
    loc37=("17", 7579472, 7579472),
    oncogenicity="Benign",
)


class TestBuildAwareSearch:
    @pytest.mark.unit
    def test_grch37_uses_chrpos37(self):
        assert _build_search_term("chr17", 7579472, "GRCh37") == (
            "17[Chromosome] AND 7579472[chrpos37]"
        )

    @pytest.mark.unit
    def test_indel_searches_reference_span(self):
        # VCF deletion pos=100 ref=ACG: records sit on the deleted bases after the pad base.
        assert _span_end(100, "ACG", "A") == 103
        assert _build_search_term("1", 100, "GRCh38", 103) == "1[Chromosome] AND 100:103[chrpos38]"
        # insertion: the event lies between pos and pos+1
        assert _span_end(100, "A", "AT") == 101
        # SNV / MNV: their own span
        assert _span_end(100, "A", "G") == 100
        assert _span_end(100, "AC", "GT") == 101

    @pytest.mark.unit
    def test_grch37_entry_matched_via_variation_loc(self):
        section = {"uids": ["12351"], "12351": RS1042522}
        assert _select_entry(section, 7579472, "G", "C", build="GRCh37", chrom="17") is RS1042522
        # the GRCh38 coordinate is NOT a GRCh37 match, and vice versa
        assert _select_entry(section, 7676154, "G", "C", build="GRCh37", chrom="17") is None
        assert _select_entry(section, 7579472, "G", "C", build="GRCh38") is None
        assert _select_entry(section, 7676154, "G", "C", build="GRCh38") is RS1042522

    @pytest.mark.unit
    def test_grch37_record_without_grch37_location_not_guessed(self):
        entry = _entry("1", "NC_000017.11:7676153:G:C", loc38=("17", 7676154, 7676154))
        section = {"uids": ["1"], "1": entry}
        assert _select_entry(section, 7579472, "G", "C", build="GRCh37", chrom="17") is None

    @pytest.mark.unit
    def test_grch37_wrong_chromosome_not_matched(self):
        section = {"uids": ["12351"], "12351": RS1042522}
        assert _select_entry(section, 7579472, "G", "C", build="GRCh37", chrom="chr3") is None

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_grch37_lookup_queries_chrpos37_and_finds_record(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*chrpos37.*"),
            json={"esearchresult": {"count": "1", "idlist": ["12351"]}},
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={"result": {"uids": ["12351"], "12351": RS1042522}},
        )
        client = ClinVarClient()
        result = await client.lookup("chr17", 7579472, "G", "C", build="hg19")
        assert result is not None
        assert result.variation_id == 12351
        assert result.build == "GRCh37"
        assert result.clinical_significance == "Benign"
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_default_build_comes_from_client(self, httpx_mock):
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*chrpos37.*"), json=ESEARCH_EMPTY)
        client = ClinVarClient(genome_build="GRCh37")
        assert await client.lookup("17", 7579472, "G", "C") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_unknown_build_is_unavailable_not_none(self):
        client = ClinVarClient()
        with pytest.raises(ClinVarUnavailableError, match="build"):
            await client.lookup("17", 7579472, "G", "C", build="unknown")
        await client.close()


class TestHonestClassification:
    @pytest.mark.unit
    def test_no_germline_classification_is_none_not_vus(self):
        entry = _entry("9", "NC_000017.11:1:G:A", description="")
        result = _parse_entry(entry)
        assert result is not None
        assert result.clinical_significance is None
        assert result.classification_type is None
        assert result.stars == 0

    @pytest.mark.unit
    def test_somatic_classification_does_not_feed_germline(self):
        entry = _entry("9", "NC_000017.11:1:G:A", description="", oncogenicity="Oncogenic")
        result = _parse_entry(entry)
        assert result is not None
        assert result.clinical_significance is None
        assert result.somatic_classifications == {"oncogenicity": "Oncogenic"}

    @pytest.mark.unit
    def test_germline_and_somatic_kept_apart(self):
        result = _parse_entry(RS1042522)
        assert result is not None
        assert result.clinical_significance == "Benign"
        assert result.classification_type == "germline"
        assert result.somatic_classifications == {"oncogenicity": "Benign"}

    @pytest.mark.unit
    def test_legacy_schema_tagged_legacy(self):
        entry = {
            "uid": "5",
            "clinical_significance": {"description": "Pathogenic"},
            "review_status": "criteria provided, single submitter",
        }
        result = _parse_entry(entry)
        assert result is not None
        assert result.classification_type == "legacy"
        assert result.stars == 1


class TestInBandErrors:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_esearch_error_is_unavailable_and_not_cached(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"ERROR": "Search Backend failed"}},
        )
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), json=ESEARCH_EMPTY)
        client = ClinVarClient()
        with pytest.raises(NCBIUnavailableError):
            await client.lookup("17", 7674220, "G", "A")
        # not cached: the retry reaches NCBI and gets a real answer
        assert await client.lookup("17", 7674220, "G", "A") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_unknown_search_field_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "count": "0",
                    "idlist": [],
                    "errorlist": {"phrasesnotfound": [], "fieldsnotfound": ["chrpos38"]},
                }
            },
        )
        client = ClinVarClient()
        with pytest.raises(NCBIUnavailableError):
            await client.lookup("17", 7674220, "G", "A")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_phrases_not_found_is_an_ordinary_empty_result(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={
                "esearchresult": {
                    "count": "0",
                    "idlist": [],
                    "errorlist": {"phrasesnotfound": ["099[chrpos38]"], "fieldsnotfound": []},
                }
            },
        )
        client = ClinVarClient()
        assert await client.lookup("17", 99, "G", "A") is None
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_esummary_uid_error_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), json=ESEARCH_RESPONSE)
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "uids": ["12345"],
                    "12345": {"uid": "12345", "error": "cannot get document summary"},
                }
            },
        )
        client = ClinVarClient()
        with pytest.raises(NCBIUnavailableError):
            await client.lookup("17", 7674220, "G", "A")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_html_error_page_is_unavailable(self, httpx_mock):
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"), text="<html>Service Unavailable</html>"
        )
        client = ClinVarClient()
        with pytest.raises(NCBIUnavailableError, match="unparseable"):
            await client.lookup("17", 7674220, "G", "A")
        await client.close()

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_truncated_search_without_match_is_unavailable(self, httpx_mock, monkeypatch):
        import bamcp.clients.clinvar as cv

        monkeypatch.setattr(cv, "_MAX_POSITION_RECORDS", 2)
        httpx_mock.add_response(
            url=re.compile(r".*/esearch\.fcgi.*"),
            json={"esearchresult": {"count": "5", "idlist": ["1", "2"]}},
        )
        httpx_mock.add_response(
            url=re.compile(r".*/esummary\.fcgi.*"),
            json={
                "result": {
                    "uids": ["1", "2"],
                    "1": _entry("1", "NC_000017.11:7674219:G:C"),
                    "2": _entry("2", "NC_000017.11:7674219:G:T"),
                }
            },
        )
        client = ClinVarClient()
        with pytest.raises(ClinVarUnavailableError, match="truncated"):
            await client.lookup("17", 7674220, "G", "A")
        await client.close()


class TestSharedRateLimiter:
    @pytest.mark.unit
    def test_injected_limiter_is_used(self):
        limiter = TokenBucketRateLimiter(3)
        client = ClinVarClient(rate_limiter=limiter)
        assert client._rate_limiter is limiter


# -- NCBI throttling: bounded retry with backoff ---------------------------------


class TestNcbiRetry:
    @pytest.fixture
    def sleeps(self, monkeypatch):
        import bamcp.clients.ncbi as ncbi

        recorded: list[float] = []

        async def fake_sleep(seconds):
            recorded.append(seconds)

        monkeypatch.setattr(ncbi.asyncio, "sleep", fake_sleep)
        return recorded

    @staticmethod
    def _client():
        class NoWait(TokenBucketRateLimiter):
            async def acquire(self):  # the bucket's own waits are not retries
                return None

        return ClinVarClient(rate_limiter=NoWait(3))

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_429_is_retried_honoring_retry_after(self, httpx_mock, sleeps):
        esearch = re.compile(r".*/esearch\.fcgi.*")
        httpx_mock.add_response(url=esearch, status_code=429, headers={"Retry-After": "2"})
        httpx_mock.add_response(url=esearch, json=ESEARCH_RESPONSE)
        httpx_mock.add_response(url=re.compile(r".*/esummary\.fcgi.*"), json=ESUMMARY_RESPONSE)

        result = await self._client().lookup("chr17", 7674220, "G", "A")
        assert result is not None and result.gene == "TP53"
        assert sleeps == [2.0]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_retry_is_bounded_and_backs_off(self, httpx_mock, sleeps):
        from bamcp.clients.ncbi import NCBI_MAX_RETRIES

        for _ in range(NCBI_MAX_RETRIES + 1):
            httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), status_code=429)
        with pytest.raises(httpx.HTTPStatusError):
            await self._client().lookup("chr17", 7674220, "G", "A")
        assert sleeps == [0.5, 1.0, 2.0]  # exponential, then give up (-> unavailable)

    @pytest.mark.unit
    def test_retry_after_is_capped_and_date_form_ignored(self):
        from bamcp.clients.ncbi import NCBI_RETRY_MAX_DELAY, _retry_delay

        big = httpx.Response(429, headers={"Retry-After": "3600"})
        assert _retry_delay(big, 0) == NCBI_RETRY_MAX_DELAY
        date = httpx.Response(429, headers={"Retry-After": "Wed, 21 Oct 2026 07:28:00 GMT"})
        assert _retry_delay(date, 1) == 1.0

"""Unit tests for bamcp.clinvar module."""

import re

import httpx
import pytest

from bamcp.clients.clinvar import (
    REVIEW_STARS,
    ClinVarClient,
    ClinVarResult,
    _build_search_term,
    _parse_entry,
    _select_entry,
    _spdi_matches,
)

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
):
    """Build an esummary entry in ClinVar's current (germline_classification) schema."""
    return {
        "uid": uid,
        "title": title,
        "germline_classification": {
            "description": description,
            "review_status": review,
            "trait_set": [{"trait_name": t} for t in traits],
            "last_evaluated": "2023-01-15",
        },
        "genes": [{"symbol": gene}] if gene else [],
        "variation_set": [{"canonical_spdi": spdi, "variant_type": "snv"}],
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
        assert term == "17[Chromosome] AND 7674220[Base Position]"

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
        ids = await client._search_ids("chr17", 7674220)
        assert ids == ["1", "2", "3"]

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

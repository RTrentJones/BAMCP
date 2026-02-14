"""Unit tests for bamcp.clinvar module."""

import re

import httpx
import pytest

from bamcp.clients.clinvar import (
    REVIEW_STARS,
    ClinVarClient,
    ClinVarResult,
    _build_search_term,
    _parse_summary,
)

# -- Sample API response fixtures -------------------------------------------

ESEARCH_RESPONSE = {
    "esearchresult": {
        "count": "1",
        "retmax": "5",
        "idlist": ["12345"],
    }
}

ESEARCH_EMPTY = {
    "esearchresult": {
        "count": "0",
        "retmax": "5",
        "idlist": [],
    }
}

ESUMMARY_RESPONSE = {
    "result": {
        "uids": ["12345"],
        "12345": {
            "uid": "12345",
            "title": "NM_000546.6(TP53):c.743G>A (p.Arg248Gln)",
            "clinical_significance": {
                "description": "Pathogenic",
            },
            "review_status": "criteria provided, multiple submitters, no conflicts",
            "genes": [{"symbol": "TP53"}],
            "trait_set": [
                {"trait_name": "Li-Fraumeni syndrome"},
                {"trait_name": "Hereditary cancer-predisposing syndrome"},
            ],
            "last_evaluated": "2023-01-15",
        },
    }
}

ESUMMARY_EMPTY = {
    "result": {
        "uids": [],
    }
}


class TestBuildSearchTerm:
    """Tests for _build_search_term."""

    @pytest.mark.unit
    def test_basic_search_term(self):
        term = _build_search_term("chr17", 7674220, "G", "A")
        assert "17" in term
        assert "7674220" in term
        assert "G>A" in term
        assert "[Chromosome]" in term
        assert "[Base Position]" in term

    @pytest.mark.unit
    def test_strips_chr_prefix(self):
        term = _build_search_term("chr1", 100, "A", "T")
        # Should not have "chr" in the chromosome part
        assert term.startswith("1[")

    @pytest.mark.unit
    def test_numeric_chrom(self):
        term = _build_search_term("17", 100, "A", "T")
        assert term.startswith("17[")


class TestParseSummary:
    """Tests for _parse_summary."""

    @pytest.mark.unit
    def test_parse_valid_response(self):
        result = _parse_summary(ESUMMARY_RESPONSE)
        assert result is not None
        assert result.variation_id == 12345
        assert result.clinical_significance == "Pathogenic"
        assert result.review_status == "criteria provided, multiple submitters, no conflicts"
        assert result.stars == 2
        assert "Li-Fraumeni syndrome" in result.conditions
        assert result.gene == "TP53"
        assert result.variant_name == "NM_000546.6(TP53):c.743G>A (p.Arg248Gln)"

    @pytest.mark.unit
    def test_parse_empty_response(self):
        result = _parse_summary(ESUMMARY_EMPTY)
        assert result is None

    @pytest.mark.unit
    def test_parse_no_result_key(self):
        result = _parse_summary({})
        assert result is None

    @pytest.mark.unit
    def test_parse_significance_as_string(self):
        data = {
            "result": {
                "uids": ["999"],
                "999": {
                    "uid": "999",
                    "clinical_significance": "Benign",
                    "review_status": "no assertion criteria provided",
                    "genes": [],
                    "trait_set": [],
                },
            }
        }
        result = _parse_summary(data)
        assert result is not None
        assert result.clinical_significance == "Benign"
        assert result.stars == 0

    @pytest.mark.unit
    def test_parse_gene_as_string(self):
        data = {
            "result": {
                "uids": ["100"],
                "100": {
                    "uid": "100",
                    "clinical_significance": {"description": "VUS"},
                    "review_status": "criteria provided, single submitter",
                    "genes": ["BRCA1"],
                    "trait_set": [],
                },
            }
        }
        result = _parse_summary(data)
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

"""Shared NCBI E-utilities plumbing: rate limiting and in-band error detection.

ClinVar and Gene lookups both hit ``eutils.ncbi.nlm.nih.gov``, and NCBI's rate limit
(3 req/s without an API key, 10 with one) is per caller, not per database — so the two
clients must draw from ONE token bucket. :class:`TokenBucketRateLimiter` lives here so a
server can build a single instance and hand it to both.

E-utilities also report many failures *in band* with HTTP 200 (``esearchresult.ERROR``,
``errorlist.fieldsnotfound``, a per-uid ``error`` in esummary). Parsing those as "no
records" would turn an outage or a malformed query into a false "not in ClinVar" — the
helpers below raise :class:`NCBIUnavailableError` instead, which callers report as
``unavailable`` and never cache.
"""

from __future__ import annotations

import asyncio
import time
from typing import Any

import httpx

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# NCBI rate limits (requests per second)
NCBI_RATE_LIMIT_NO_KEY = 3  # 3 req/sec without API key
NCBI_RATE_LIMIT_WITH_KEY = 10  # 10 req/sec with API key


class NCBIUnavailableError(Exception):
    """An NCBI lookup could not be answered — distinct from a genuine not-found.

    Raised for in-band E-utilities errors, unparseable (e.g. HTML) responses, truncated
    result sets and unsupported genome builds. Never cached.
    """


class TokenBucketRateLimiter:
    """Async token bucket rate limiter for API requests.

    Enforces a maximum number of requests per second.
    """

    def __init__(self, rate: float):
        """Initialize rate limiter.

        Args:
            rate: Maximum requests per second.
        """
        self._rate = rate
        self._tokens = rate
        self._last_refill = time.monotonic()
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        """Acquire a token, waiting if necessary."""
        async with self._lock:
            now = time.monotonic()
            elapsed = now - self._last_refill
            self._tokens = min(self._rate, self._tokens + elapsed * self._rate)
            self._last_refill = now

            if self._tokens < 1:
                wait_time = (1 - self._tokens) / self._rate
                await asyncio.sleep(wait_time)
                self._tokens = 0
            else:
                self._tokens -= 1


# Bounded retry for NCBI's shared-quota throttling (429) and transient gateway errors. Under
# shared demo load several sessions draw on one NCBI quota, so a 429 is expected, not an outage.
NCBI_RETRY_STATUSES = frozenset({429, 502, 503, 504})
NCBI_MAX_RETRIES = 3
NCBI_RETRY_BASE_DELAY = 0.5  # seconds; doubles per attempt
NCBI_RETRY_MAX_DELAY = 8.0  # cap, also on a server-sent Retry-After


def _retry_delay(resp: httpx.Response, attempt: int) -> float:
    """Seconds to wait before retry ``attempt`` (0-based): Retry-After if numeric, else backoff."""
    header = str(resp.headers.get("retry-after") or "").strip()
    delay = float(NCBI_RETRY_BASE_DELAY * (2**attempt))
    # Numeric seconds only; the HTTP-date form falls back to exponential backoff.
    if header.replace(".", "", 1).isdigit():
        delay = float(header)
    return max(0.0, min(delay, NCBI_RETRY_MAX_DELAY))


async def ncbi_get(
    client: httpx.AsyncClient,
    url: str,
    params: dict[str, Any],
    rate_limiter: TokenBucketRateLimiter,
    *,
    max_retries: int = NCBI_MAX_RETRIES,
) -> httpx.Response:
    """GET an E-utilities URL under the shared token bucket, retrying 429/5xx a bounded number
    of times (honoring a numeric ``Retry-After``, capped). Returns the last response; the caller
    still checks its status, so an exhausted retry surfaces as ``unavailable``."""
    attempt = 0
    while True:
        await rate_limiter.acquire()
        resp = await client.get(url, params=params)
        if resp.status_code not in NCBI_RETRY_STATUSES or attempt >= max_retries:
            return resp
        await asyncio.sleep(_retry_delay(resp, attempt))
        attempt += 1


def ncbi_rate_limiter(api_key: str | None) -> TokenBucketRateLimiter:
    """Build a token bucket sized to NCBI's limit for this key (or no key)."""
    return TokenBucketRateLimiter(NCBI_RATE_LIMIT_WITH_KEY if api_key else NCBI_RATE_LIMIT_NO_KEY)


def response_json(resp: httpx.Response, what: str) -> Any:
    """Decode a response body as JSON, mapping an HTML/garbled body to unavailable."""
    try:
        return resp.json()
    except ValueError as e:
        raise NCBIUnavailableError(f"{what}: unparseable response from NCBI") from e


def check_esearch(payload: Any, what: str) -> dict:
    """Return ``esearchresult``, raising on NCBI's in-band (HTTP 200) error forms.

    ``errorlist.phrasesnotfound`` is NOT an error — it is how esearch reports an ordinary
    empty result (e.g. a position with no records) — but ``ERROR`` and
    ``fieldsnotfound`` mean the query itself was not executed as written.
    """
    if not isinstance(payload, dict):
        raise NCBIUnavailableError(f"{what}: malformed esearch response")
    if payload.get("error") or payload.get("ERROR"):
        raise NCBIUnavailableError(
            f"{what}: NCBI error: {payload.get('error') or payload['ERROR']}"
        )
    result = payload.get("esearchresult")
    if not isinstance(result, dict):
        raise NCBIUnavailableError(f"{what}: esearch response has no esearchresult")
    if result.get("ERROR"):
        raise NCBIUnavailableError(f"{what}: NCBI esearch error: {result['ERROR']}")
    errorlist = result.get("errorlist") or {}
    if isinstance(errorlist, dict) and errorlist.get("fieldsnotfound"):
        raise NCBIUnavailableError(
            f"{what}: NCBI did not recognize search field(s) {errorlist['fieldsnotfound']}"
        )
    return result


def check_esummary(payload: Any, what: str) -> dict:
    """Return esummary's ``result`` map, raising on in-band errors.

    Covers a top-level ``error`` / ``esummaryresult`` message and a per-uid ``error``
    entry ("cannot get document summary") — a record NCBI failed to serve is unknown,
    not absent.
    """
    if not isinstance(payload, dict):
        raise NCBIUnavailableError(f"{what}: malformed esummary response")
    if payload.get("error") or payload.get("ERROR"):
        raise NCBIUnavailableError(
            f"{what}: NCBI error: {payload.get('error') or payload['ERROR']}"
        )
    result = payload.get("result")
    if not isinstance(result, dict):
        detail = payload.get("esummaryresult")
        raise NCBIUnavailableError(f"{what}: esummary response has no result ({detail!r})")
    for uid in result.get("uids", []):
        entry = result.get(uid)
        if isinstance(entry, dict) and entry.get("error"):
            raise NCBIUnavailableError(f"{what}: NCBI esummary error for {uid}: {entry['error']}")
    return result

"""Bounded LRU cache with TTL eviction.

Shared by ClinVar and gnomAD clients to avoid duplicate cache implementations.
"""

from __future__ import annotations

import asyncio
import time
from collections import OrderedDict
from typing import Generic, TypeVar

T = TypeVar("T")

# Default cache configuration for API clients
API_CACHE_MAX_SIZE = 1000  # Maximum number of cached entries
API_CACHE_TTL_SECONDS = 3600  # 1 hour TTL


class _Miss:
    """Sentinel distinguishing a cache miss from a cached ``None`` value."""

    __slots__ = ()

    def __repr__(self) -> str:  # pragma: no cover - debug aid
        return "MISS"


# Public sentinel: ``get_entry`` returns this on a miss so ``None`` can be a real cached value
# (e.g. a negative "variant not found" result), avoiding a re-fetch on every subsequent lookup.
MISS = _Miss()


class BoundedTTLCache(Generic[T]):
    """Thread-safe bounded cache with TTL eviction.

    Implements LRU eviction when maxsize is exceeded and TTL-based expiry.
    """

    def __init__(self, maxsize: int = API_CACHE_MAX_SIZE, ttl: float = API_CACHE_TTL_SECONDS):
        self._cache: OrderedDict[tuple, tuple[T, float]] = OrderedDict()
        self._maxsize = maxsize
        self._ttl = ttl
        self._lock = asyncio.Lock()

    async def get(self, key: tuple) -> T | None:
        """Get value from cache, returning None if missing or expired.

        Cannot distinguish a cached ``None`` from a miss — use :meth:`get_entry` when
        ``None`` is a meaningful cached value (negative results).
        """
        result = await self.get_entry(key)
        return None if result is MISS else result  # type: ignore[return-value]

    async def get_entry(self, key: tuple) -> T | _Miss:
        """Get value from cache, returning the ``MISS`` sentinel if missing or expired.

        Lets callers cache ``None`` (an absence) and read it back as a hit, so a
        not-found lookup is served from cache instead of re-hitting the upstream API.
        """
        async with self._lock:
            if key not in self._cache:
                return MISS

            value, timestamp = self._cache[key]
            if time.monotonic() - timestamp > self._ttl:
                del self._cache[key]
                return MISS

            # Move to end (most recently used)
            self._cache.move_to_end(key)
            return value

    async def set(self, key: tuple, value: T) -> None:
        """Set value in cache, evicting oldest if at capacity."""
        async with self._lock:
            if key in self._cache:
                del self._cache[key]
            elif len(self._cache) >= self._maxsize:
                # Evict oldest (first) item
                self._cache.popitem(last=False)

            self._cache[key] = (value, time.monotonic())

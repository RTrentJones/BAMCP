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
        """Get value from cache, returning None if missing or expired."""
        async with self._lock:
            if key not in self._cache:
                return None

            value, timestamp = self._cache[key]
            if time.monotonic() - timestamp > self._ttl:
                del self._cache[key]
                return None

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

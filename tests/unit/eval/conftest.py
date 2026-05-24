"""Shared fixtures for eval-package unit tests.

Provides a FakeRenderer so most runner tests don't have to launch Chromium.
"""

from __future__ import annotations

from typing import Any

import pytest

PNG_SENTINEL = b"\x89PNG\r\n\x1a\n" + b"FAKE_PNG_BYTES_FOR_TEST"


class FakeRenderer:
    """Stand-in for ViewerRenderer that returns a fixed PNG sentinel.

    Records every capture() invocation so tests can assert on arguments.
    """

    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []
        self.entered = False
        self.exited = False
        self._fail_next = False

    def make_fail_next(self) -> None:
        self._fail_next = True

    async def __aenter__(self) -> FakeRenderer:
        self.entered = True
        return self

    async def __aexit__(self, *exc: Any) -> None:
        self.exited = True

    async def capture(
        self,
        ui_init: dict,
        mode: str,
        *,
        active_variant_position: int | None = None,
        active_variant_alt: str | None = None,
        clip_selector: str = "#container",
    ) -> bytes:
        if self._fail_next:
            self._fail_next = False
            raise RuntimeError("fake capture failure")
        self.calls.append(
            {
                "ui_init": ui_init,
                "mode": mode,
                "active_variant_position": active_variant_position,
                "active_variant_alt": active_variant_alt,
            }
        )
        return PNG_SENTINEL


@pytest.fixture
def fake_renderer():
    return FakeRenderer()

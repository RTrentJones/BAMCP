"""Headless Chromium that renders the BAMCP viewer for vision-eval captures.

Hydrates ``viewer.html`` (the same HTML the production server serves) with a
tool's ``ui/init`` payload, sets the chosen display mode, and returns a PNG of
the reads + coverage + ruler region — exactly what a human user would see in
Claude Desktop.

Used in two places:

- :mod:`bamcp.eval.runner` — visual eval captures one screenshot per UI tool
  call and embeds it as an image content block in the next provider turn.
- :mod:`scripts.render_viewer` — developer CLI for the LLM feedback loop:
  edit ``renderer.ts`` → rebuild → render → Read PNG → iterate.

The renderer is async-context-managed: ``async with ViewerRenderer() as r:``
launches Chromium on entry and closes it on exit. Each capture loads a fresh
page so successive calls do not share state.
"""

from __future__ import annotations

from contextlib import AsyncExitStack
from typing import Any

from ..resources import get_viewer_html

_DEFAULT_VIEWPORT_WIDTH = 1024
_DEFAULT_VIEWPORT_HEIGHT = 800

# JS that hydrates the viewer with a payload and applies a display mode.
# Mirrors the e2e test pattern at tests/e2e/test_viewer_e2e.py:96-110 so the
# headless behavior matches what the existing Playwright tests already verify.
_HYDRATE_JS = """
(({payload, mode, activeVariantPosition, activeVariantAlt}) => {
    viewer.state.loadData(payload);
    viewer.state.settings.displayMode = mode;
    if (activeVariantPosition !== null) {
        viewer.state.settings.activeVariantPosition = activeVariantPosition;
        viewer.state.settings.activeVariantAlt = activeVariantAlt;
    }
    viewer.state.resortAndRepack();
    viewer.renderVariantTable();
    viewer.renderer.resize();
})
"""


class ViewerRenderer:
    """Async context manager that captures viewer screenshots via Playwright.

    Usage::

        async with ViewerRenderer() as renderer:
            png_bytes = await renderer.capture(ui_init, "dv-strips")

    Raises a clear error on enter if Playwright or Chromium is not installed.
    """

    def __init__(
        self,
        *,
        viewport_width: int = _DEFAULT_VIEWPORT_WIDTH,
        viewport_height: int = _DEFAULT_VIEWPORT_HEIGHT,
        headless: bool = True,
    ) -> None:
        self.viewport_width = viewport_width
        self.viewport_height = viewport_height
        self.headless = headless
        self._stack: AsyncExitStack | None = None
        self._browser: Any = None
        self._html: str | None = None

    async def __aenter__(self) -> ViewerRenderer:
        try:
            from playwright.async_api import async_playwright
        except ImportError as e:  # pragma: no cover — install hint path
            raise RuntimeError(
                "playwright is required for the visual eval renderer. "
                "Install with 'pip install \".[eval]\"' and run "
                "'python -m playwright install chromium'."
            ) from e

        self._stack = AsyncExitStack()
        playwright = await self._stack.enter_async_context(async_playwright())
        self._browser = await playwright.chromium.launch(
            headless=self.headless, args=["--no-sandbox"]
        )
        self._html = get_viewer_html()
        return self

    async def __aexit__(self, *exc: Any) -> None:
        if self._stack is not None:
            await self._stack.aclose()
            self._stack = None
        self._browser = None
        self._html = None

    async def capture(
        self,
        ui_init: dict,
        mode: str,
        *,
        active_variant_position: int | None = None,
        active_variant_alt: str | None = None,
        clip_selector: str = "#container",
    ) -> bytes:
        """Render ``ui_init`` at ``mode`` and return PNG bytes.

        Args:
            ui_init: The viewer payload (the same dict returned by
                ``visualize_region`` / ``jump_to`` in ``_meta.ui/init``).
            mode: One of ``squished``/``compact``/``expanded``/``dv-strips``/
                ``dv-composite``.
            active_variant_position: Optional genomic position to mark active.
                Drives the "supports-variant" channel in the DV modes.
            active_variant_alt: Alt base for the active variant.
            clip_selector: CSS selector for the region to screenshot. Defaults
                to the entire viewer container.

        Returns:
            PNG bytes ready to base64-encode or write to disk.
        """
        if self._browser is None or self._html is None:
            raise RuntimeError("ViewerRenderer used outside of an `async with` block")

        context = await self._browser.new_context(
            viewport={"width": self.viewport_width, "height": self.viewport_height}
        )
        try:
            page = await context.new_page()
            await page.set_content(self._html, wait_until="domcontentloaded")
            await page.wait_for_function("() => typeof window.viewer !== 'undefined'", timeout=5000)
            await page.evaluate(
                _HYDRATE_JS,
                {
                    "payload": ui_init,
                    "mode": mode,
                    "activeVariantPosition": active_variant_position,
                    "activeVariantAlt": active_variant_alt,
                },
            )
            # Give the renderer a moment to paint after resize().
            await page.wait_for_timeout(150)
            locator = page.locator(clip_selector)
            png: bytes = await locator.screenshot(type="png")
            return png
        finally:
            await context.close()

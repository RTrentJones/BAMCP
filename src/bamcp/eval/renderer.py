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

Blinded captures (``capture(..., blind=True)``) are the stimulus path for the
visual-vs-text eval: the payload is reduced to what the pileup draws (no variant
calls, no evidence, hashed read names, no file path), overlays marked
``data-blind-hide`` plus the toolbar and variant/evidence panel are hidden, the
display settings are pinned, the shot is clipped to ``#viewer`` and the viewport
height grows to fit the full read stack.
"""

from __future__ import annotations

import copy
import hashlib
import os
from contextlib import AsyncExitStack
from typing import Any

from ..resources import get_viewer_html

_DEFAULT_VIEWPORT_WIDTH = 1024
_DEFAULT_VIEWPORT_HEIGHT = 800
# Upper bound for an auto-fit (depth-sized) viewport; well under Chromium's canvas cap.
_MAX_AUTO_VIEWPORT_HEIGHT = 8000

CLIP_VIEWER = "#viewer"
CLIP_CONTAINER = "#container"

# Payload keys a blinded render keeps — everything the pileup draws and nothing that
# carries a verdict (variants, variant_evidence, applied_filters, summaries) or a label
# (file_path, reference/vcf paths).
_BLIND_KEEP_KEYS = frozenset(
    {
        "schema_version",
        "contig",
        "start",
        "end",
        "reads",
        "coverage",
        "reference_sequence",
        "compact",
        "coordinate_system",
    }
)

# Display settings pinned for a blinded render, so every stimulus is drawn the same way
# (and no sort-by-base at a site hints where to look).
BLIND_SETTINGS: dict[str, Any] = {
    "colorBy": "strand",
    "sortBy": "position",
    "showSoftClips": False,
    "showMismatches": True,
    "activeVariantPosition": None,
    "activeVariantAlt": None,
}

# Hidden in a blinded render. ``[data-blind-hide]`` is the viewer's contract for
# verdict-bearing overlays inside #viewer; the ids below are belt-and-braces for the
# panels/overlays that predate that attribute.
_BLIND_CSS = """
[data-blind-hide], #toolbar, #bottom-panel, #variant-panel, #evidence-panel,
#tooltip, #sequence-notice, #fetch-error, #debug-overlay {
    display: none !important;
}
"""


def _hash_read_name(name: Any) -> Any:
    """Stable, opaque stand-in for a read name (mates keep a shared hash)."""
    if name is None:
        return None
    return "r" + hashlib.sha256(str(name).encode("utf-8")).hexdigest()[:12]


def blind_payload(ui_init: dict) -> dict:
    """Return a copy of ``ui_init`` stripped of everything that could leak a verdict.

    Keeps only the drawable keys, forces ``variants: []`` and replaces every read name with
    a stable hash (so mate pairing still works but labels like ``d1_alt_fwd_0`` vanish).
    Handles both the columnar read layout and a plain list of read objects.
    """
    out = {k: copy.deepcopy(v) for k, v in ui_init.items() if k in _BLIND_KEEP_KEYS}
    out["variants"] = []
    reads = out.get("reads")
    if isinstance(reads, dict) and isinstance(reads.get("name"), list):
        reads["name"] = [_hash_read_name(n) for n in reads["name"]]
    elif isinstance(reads, list):
        for read in reads:
            if isinstance(read, dict) and "name" in read:
                read["name"] = _hash_read_name(read["name"])
    return out


# JS that hydrates the viewer with a payload and applies a display mode.
# Mirrors the e2e test pattern at tests/e2e/test_viewer_e2e.py:96-110 so the
# headless behavior matches what the existing Playwright tests already verify.
_HYDRATE_JS = """
(({payload, mode, activeVariantPosition, activeVariantAlt, pinned}) => {
    viewer.state.loadData(payload);
    if (pinned) {
        Object.assign(viewer.state.settings, pinned);
    }
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

# How many extra pixels #viewer needs to show its whole content without scrolling.
_OVERFLOW_JS = """
(() => {
    const el = document.querySelector('#viewer');
    return el ? Math.max(0, el.scrollHeight - el.clientHeight) : 0;
})
"""

_RESIZE_JS = "(() => { viewer.renderer.resize(); })"


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
        # A sandbox may ship a Chromium whose revision differs from the pinned playwright's;
        # point at it explicitly instead of requiring `playwright install` (unset: no-op).
        executable = os.environ.get("BAMCP_CHROMIUM_EXECUTABLE") or os.environ.get(
            "PLAYWRIGHT_CHROMIUM_EXECUTABLE"
        )
        launch_kwargs: dict[str, Any] = {"headless": self.headless, "args": ["--no-sandbox"]}
        if executable:
            launch_kwargs["executable_path"] = executable
        self._browser = await playwright.chromium.launch(**launch_kwargs)
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
        blind: bool = False,
        viewport: tuple[int, int | None] | None = None,
        clip: str | None = None,
        active_variant_position: int | None = None,
        active_variant_alt: str | None = None,
        clip_selector: str | None = None,
    ) -> bytes:
        """Render ``ui_init`` at ``mode`` and return PNG bytes.

        Args:
            ui_init: The viewer payload (the same dict returned by
                ``visualize_region`` / ``jump_to`` in ``_meta.ui/init``).
            mode: One of ``squished``/``compact``/``expanded``/``dv-strips``/
                ``dv-composite``.
            blind: Render a verdict-free stimulus: the payload goes through
                :func:`blind_payload`, ``[data-blind-hide]`` overlays and the
                toolbar/variant/evidence panels are hidden, display settings are
                pinned to :data:`BLIND_SETTINGS`, the clip defaults to ``#viewer``
                and the viewport height defaults to auto-fit.
            viewport: ``(width, height)`` for this capture. ``height=None`` grows
                the viewport until ``#viewer`` shows every read row (capped at
                ``_MAX_AUTO_VIEWPORT_HEIGHT``). Omitted: the instance viewport
                (auto-fit height when ``blind``).
            clip: CSS selector to screenshot — ``"#viewer"`` (tracks only) or
                ``"#container"`` (whole app, incl. the variant table). Defaults to
                ``#viewer`` when ``blind`` else ``#container``. A blinded capture
                refuses ``#container``: the variant table leaks verdicts.
            active_variant_position: Optional genomic position to mark active.
                Drives the "supports-variant" channel in the DV modes. Not
                allowed with ``blind`` (it points at the answer).
            active_variant_alt: Alt base for the active variant.
            clip_selector: Deprecated alias for ``clip``.

        Returns:
            PNG bytes ready to base64-encode or write to disk.
        """
        if self._browser is None or self._html is None:
            raise RuntimeError("ViewerRenderer used outside of an `async with` block")

        selector = clip or clip_selector or (CLIP_VIEWER if blind else CLIP_CONTAINER)
        if blind:
            if selector == CLIP_CONTAINER:
                raise ValueError("blind captures must not clip to #container (variant table)")
            if active_variant_position is not None or active_variant_alt is not None:
                raise ValueError("blind captures must not mark an active variant")
            payload = blind_payload(ui_init)
        else:
            payload = ui_init

        if viewport is None:
            width = self.viewport_width
            height: int | None = None if blind else self.viewport_height
        else:
            width, height = viewport
        auto_fit = height is None
        current_height = self.viewport_height if height is None else height

        context = await self._browser.new_context(
            viewport={"width": width, "height": current_height}
        )
        try:
            page = await context.new_page()
            await page.set_content(self._html, wait_until="domcontentloaded")
            await page.wait_for_function("() => typeof window.viewer !== 'undefined'", timeout=5000)
            if blind:
                await page.add_style_tag(content=_BLIND_CSS)
            await page.evaluate(
                _HYDRATE_JS,
                {
                    "payload": payload,
                    "mode": mode,
                    "activeVariantPosition": active_variant_position,
                    "activeVariantAlt": active_variant_alt,
                    "pinned": BLIND_SETTINGS if blind else None,
                },
            )
            if auto_fit:
                # Grow the viewport until #viewer no longer scrolls, so no read rows are cut.
                for _ in range(3):
                    overflow = int(await page.evaluate(_OVERFLOW_JS) or 0)
                    if overflow <= 0 or current_height >= _MAX_AUTO_VIEWPORT_HEIGHT:
                        break
                    current_height = min(current_height + overflow, _MAX_AUTO_VIEWPORT_HEIGHT)
                    await page.set_viewport_size({"width": width, "height": current_height})
                    await page.evaluate(_RESIZE_JS)
            # Give the renderer a moment to paint after resize().
            await page.wait_for_timeout(150)
            locator = page.locator(selector)
            png: bytes = await locator.screenshot(type="png")
            return png
        finally:
            await context.close()

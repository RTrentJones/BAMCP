"""Unit tests for ViewerRenderer with a fake Playwright (no Chromium launched)."""

from __future__ import annotations

import sys
import types
from typing import Any

import pytest

from bamcp.eval import renderer as renderer_mod
from bamcp.eval.renderer import ViewerRenderer


class _FakeLocator:
    def __init__(self, page: _FakePage, selector: str) -> None:
        self.page = page
        self.selector = selector

    async def screenshot(self, **kw: Any) -> bytes:
        self.page.screenshots.append({"selector": self.selector, **kw})
        return b"\x89PNG-fake"


class _FakePage:
    def __init__(self, ctx: _FakeContext) -> None:
        self.ctx = ctx
        self.evaluated: list[tuple[str, Any]] = []
        self.styles: list[str] = []
        self.screenshots: list[dict[str, Any]] = []

    async def set_content(self, html: str, **kw: Any) -> None:
        self.html = html

    async def wait_for_function(self, *a: Any, **kw: Any) -> None:
        return None

    async def evaluate(self, js: str, arg: Any = None) -> Any:
        self.evaluated.append((js, arg))
        return None

    async def add_style_tag(self, *, content: str) -> None:
        self.styles.append(content)

    async def wait_for_timeout(self, ms: int) -> None:
        return None

    def locator(self, selector: str) -> _FakeLocator:
        return _FakeLocator(self, selector)


class _FakeContext:
    def __init__(self, browser: _FakeBrowser, kw: dict[str, Any]) -> None:
        self.kw = kw
        self.pages: list[_FakePage] = []
        browser.contexts.append(self)

    async def new_page(self) -> _FakePage:
        page = _FakePage(self)
        self.pages.append(page)
        return page

    async def close(self) -> None:
        return None


class _FakeBrowser:
    def __init__(self) -> None:
        self.contexts: list[_FakeContext] = []

    async def new_context(self, **kw: Any) -> _FakeContext:
        return _FakeContext(self, kw)


class _FakeChromium:
    def __init__(self) -> None:
        self.launch_kwargs: dict[str, Any] | None = None
        self.browser = _FakeBrowser()

    async def launch(self, **kw: Any) -> _FakeBrowser:
        self.launch_kwargs = kw
        return self.browser


class _FakePlaywrightCM:
    def __init__(self, chromium: _FakeChromium) -> None:
        self.chromium = chromium

    async def __aenter__(self) -> _FakePlaywrightCM:
        return self

    async def __aexit__(self, *exc: Any) -> None:
        return None


@pytest.fixture
def fake_chromium(monkeypatch) -> _FakeChromium:
    """Install a fake ``playwright.async_api`` so ``__aenter__`` never launches Chromium."""
    chromium = _FakeChromium()
    pkg = types.ModuleType("playwright")
    api = types.ModuleType("playwright.async_api")
    api.async_playwright = lambda: _FakePlaywrightCM(chromium)  # type: ignore[attr-defined]
    pkg.async_api = api  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "playwright", pkg)
    monkeypatch.setitem(sys.modules, "playwright.async_api", api)
    monkeypatch.setattr(renderer_mod, "get_viewer_html", lambda: "<html></html>")
    monkeypatch.delenv("BAMCP_CHROMIUM_EXECUTABLE", raising=False)
    monkeypatch.delenv("PLAYWRIGHT_CHROMIUM_EXECUTABLE", raising=False)
    return chromium


class TestExecutablePath:
    @pytest.mark.unit
    async def test_no_env_omits_executable_path(self, fake_chromium):
        async with ViewerRenderer():
            pass
        assert fake_chromium.launch_kwargs is not None
        assert "executable_path" not in fake_chromium.launch_kwargs
        assert fake_chromium.launch_kwargs["headless"] is True

    @pytest.mark.unit
    async def test_playwright_env_is_used(self, fake_chromium, monkeypatch):
        monkeypatch.setenv("PLAYWRIGHT_CHROMIUM_EXECUTABLE", "/opt/pw/chrome")
        async with ViewerRenderer():
            pass
        assert fake_chromium.launch_kwargs["executable_path"] == "/opt/pw/chrome"

    @pytest.mark.unit
    async def test_bamcp_env_takes_precedence(self, fake_chromium, monkeypatch):
        monkeypatch.setenv("PLAYWRIGHT_CHROMIUM_EXECUTABLE", "/opt/pw/chrome")
        monkeypatch.setenv("BAMCP_CHROMIUM_EXECUTABLE", "/opt/bamcp/chrome")
        async with ViewerRenderer():
            pass
        assert fake_chromium.launch_kwargs["executable_path"] == "/opt/bamcp/chrome"

    @pytest.mark.unit
    async def test_empty_env_is_ignored(self, fake_chromium, monkeypatch):
        monkeypatch.setenv("BAMCP_CHROMIUM_EXECUTABLE", "")
        async with ViewerRenderer():
            pass
        assert "executable_path" not in fake_chromium.launch_kwargs


def _payload() -> dict:
    return {
        "schema_version": 2,
        "contig": "chr1",
        "start": 1000,
        "end": 1300,
        "file_path": "/data/het_snv_truth.bam",
        "reads": {
            "count": 3,
            "name": ["d1_alt_fwd_0", "d1_ref_rev_1", "d1_alt_fwd_0"],
            "position": [1000, 1010, 1100],
        },
        "coverage": [3] * 300,
        "variants": [{"position": 1050, "ref": "A", "alt": "T", "confidence": "high"}],
        "variant_evidence": {"1050:A>T": {"artifact_risk": "low"}},
        "applied_filters": {"min_vaf": 0.1},
        "reference_sequence": "ACGT" * 75,
        "compact": False,
    }


class TestBlindPayload:
    @pytest.mark.unit
    def test_strips_verdicts_and_labels(self):
        out = renderer_mod.blind_payload(_payload())
        assert out["variants"] == []
        for leaked in ("variant_evidence", "applied_filters", "file_path"):
            assert leaked not in out
        assert out["coverage"] == [3] * 300
        assert out["reference_sequence"] == "ACGT" * 75

    @pytest.mark.unit
    def test_read_names_hashed_stably_and_mates_kept(self):
        names = renderer_mod.blind_payload(_payload())["reads"]["name"]
        assert not any("alt" in n or "ref" in n for n in names)
        assert names[0] == names[2] != names[1]
        assert renderer_mod.blind_payload(_payload())["reads"]["name"] == names

    @pytest.mark.unit
    def test_list_form_reads_and_input_untouched(self):
        original = _payload()
        original["reads"] = [{"name": "d1_alt_fwd_0", "position": 1}, {"position": 2}]
        out = renderer_mod.blind_payload(original)
        assert out["reads"][0]["name"].startswith("r")
        assert "alt" not in out["reads"][0]["name"]
        assert original["reads"][0]["name"] == "d1_alt_fwd_0"
        assert original["variants"]


def _last_page(chromium: _FakeChromium) -> _FakePage:
    return chromium.browser.contexts[-1].pages[-1]


class TestCaptureOptions:
    @pytest.mark.unit
    async def test_default_capture_is_backward_compatible(self, fake_chromium):
        async with ViewerRenderer() as r:
            await r.capture(_payload(), "expanded")
        ctx = fake_chromium.browser.contexts[-1]
        page = _last_page(fake_chromium)
        assert ctx.kw["viewport"] == {"width": 1024, "height": 800}
        assert page.screenshots[-1]["selector"] == "#container"
        assert page.styles == []
        hydrate_arg = page.evaluated[0][1]
        assert hydrate_arg["payload"]["variants"]  # untouched
        assert hydrate_arg["pinned"] is None

    @pytest.mark.unit
    async def test_blind_clips_viewer_hides_overlays_strips_payload(self, fake_chromium):
        async with ViewerRenderer() as r:
            await r.capture(_payload(), mode="expanded", blind=True)
        page = _last_page(fake_chromium)
        assert page.screenshots[-1]["selector"] == "#viewer"
        css = "".join(page.styles)
        assert "[data-blind-hide]" in css
        assert "#bottom-panel" in css
        assert "display: none" in css
        hydrate_arg = page.evaluated[0][1]
        assert hydrate_arg["payload"]["variants"] == []
        assert "variant_evidence" not in hydrate_arg["payload"]
        assert "d1_alt_fwd_0" not in str(hydrate_arg["payload"])
        assert hydrate_arg["pinned"] == renderer_mod.BLIND_SETTINGS
        assert hydrate_arg["mode"] == "expanded"

    @pytest.mark.unit
    async def test_blind_refuses_container_clip_and_active_variant(self, fake_chromium):
        async with ViewerRenderer() as r:
            with pytest.raises(ValueError, match="container"):
                await r.capture(_payload(), "expanded", blind=True, clip="#container")
            with pytest.raises(ValueError, match="active variant"):
                await r.capture(_payload(), "dv-strips", blind=True, active_variant_position=1050)

    @pytest.mark.unit
    async def test_explicit_viewport_and_clip(self, fake_chromium):
        async with ViewerRenderer() as r:
            await r.capture(_payload(), "compact", viewport=(640, 480), clip="#viewer")
        ctx = fake_chromium.browser.contexts[-1]
        assert ctx.kw["viewport"] == {"width": 640, "height": 480}
        assert _last_page(fake_chromium).screenshots[-1]["selector"] == "#viewer"

    @pytest.mark.unit
    async def test_clip_selector_alias_still_works(self, fake_chromium):
        async with ViewerRenderer() as r:
            await r.capture(_payload(), "compact", clip_selector="#viewer")
        assert _last_page(fake_chromium).screenshots[-1]["selector"] == "#viewer"

    @pytest.mark.unit
    async def test_auto_height_grows_to_fit_reads(self, fake_chromium, monkeypatch):
        overflows = iter([700, 0])
        resized: list[dict] = []

        async def evaluate(self, js, arg=None):
            self.evaluated.append((js, arg))
            if js == renderer_mod._OVERFLOW_JS:
                return next(overflows)
            return None

        async def set_viewport_size(self, size):
            resized.append(size)

        monkeypatch.setattr(_FakePage, "evaluate", evaluate)
        monkeypatch.setattr(_FakePage, "set_viewport_size", set_viewport_size, raising=False)
        async with ViewerRenderer() as r:
            await r.capture(_payload(), "expanded", blind=True, viewport=(900, None))
        assert resized == [{"width": 900, "height": 1500}]

    @pytest.mark.unit
    async def test_auto_height_is_capped(self, fake_chromium, monkeypatch):
        resized: list[dict] = []

        async def evaluate(self, js, arg=None):
            return 10**6 if js == renderer_mod._OVERFLOW_JS else None

        async def set_viewport_size(self, size):
            resized.append(size)

        monkeypatch.setattr(_FakePage, "evaluate", evaluate)
        monkeypatch.setattr(_FakePage, "set_viewport_size", set_viewport_size, raising=False)
        async with ViewerRenderer() as r:
            await r.capture(_payload(), "expanded", blind=True)
        assert resized[-1]["height"] == renderer_mod._MAX_AUTO_VIEWPORT_HEIGHT
        assert len(resized) == 1

    @pytest.mark.unit
    async def test_capture_outside_context_raises(self):
        with pytest.raises(RuntimeError, match="async with"):
            await ViewerRenderer().capture({}, "expanded")

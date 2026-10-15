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

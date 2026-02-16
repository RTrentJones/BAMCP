"""E2E test configuration - Playwright tests use sync API with their own event loop."""

import contextlib
import os

import pytest
from playwright.sync_api import Page

FIXTURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "fixtures")


@pytest.fixture(scope="session")
def comprehensive_data():
    """Generate comprehensive BAM and return a loader function for any region.

    Usage:
        data = comprehensive_data("chr1:1000-1100")
        send_init_data(page, data)
    """
    from bamcp.core.parsers import fetch_region
    from bamcp.core.serialization import serialize_region_data
    from tests.create_fixtures import (
        create_comprehensive_bam,
        create_comprehensive_reference,
    )

    ref_path = create_comprehensive_reference()
    bam_path = create_comprehensive_bam(ref_path)

    def load(region: str, min_vaf: float = 0.02, min_depth: int = 2) -> dict:
        data = fetch_region(bam_path, region, ref_path, min_vaf=min_vaf, min_depth=min_depth)
        return serialize_region_data(data)

    return load


# ── Canvas pixel helpers ─────────────────────────────────────────────


def has_pixels_in_region(page: Page, canvas_id: str, x1: int, y1: int, w: int, h: int) -> bool:
    """Check if any opaque pixels exist in a canvas rectangle."""
    return page.evaluate(
        f"""() => {{
        const canvas = document.getElementById('{canvas_id}');
        const ctx = canvas.getContext('2d');
        const data = ctx.getImageData({x1}, {y1}, {w}, {h}).data;
        for (let i = 3; i < data.length; i += 4) {{
            if (data[i] > 0) return true;
        }}
        return false;
    }}"""
    )


def count_opaque_pixels(page: Page, canvas_id: str, x1: int, y1: int, w: int, h: int) -> int:
    """Count opaque pixels in a canvas rectangle."""
    return page.evaluate(
        f"""() => {{
        const canvas = document.getElementById('{canvas_id}');
        const ctx = canvas.getContext('2d');
        const data = ctx.getImageData({x1}, {y1}, {w}, {h}).data;
        let count = 0;
        for (let i = 3; i < data.length; i += 4) {{
            if (data[i] > 0) count++;
        }}
        return count;
    }}"""
    )


def sample_pixel(page: Page, canvas_id: str, x: int, y: int) -> dict:
    """Sample a single pixel at (x, y) on a canvas. Returns {r, g, b, a}."""
    return page.evaluate(
        f"""() => {{
        const canvas = document.getElementById('{canvas_id}');
        const ctx = canvas.getContext('2d');
        const p = ctx.getImageData({x}, {y}, 1, 1).data;
        return {{ r: p[0], g: p[1], b: p[2], a: p[3] }};
    }}"""
    )


def sample_pixel_at_genomic(page: Page, canvas_id: str, genomic_pos: int, y: int) -> dict:
    """Sample a pixel at a genomic position (auto-converts via viewport scale)."""
    return page.evaluate(
        f"""() => {{
        const canvas = document.getElementById('{canvas_id}');
        const ctx = canvas.getContext('2d');
        const scale = canvas.width / (viewer.state.viewport.end - viewer.state.viewport.start);
        const x = Math.floor(({genomic_pos} - viewer.state.viewport.start) * scale);
        const p = ctx.getImageData(x, {y}, 1, 1).data;
        return {{ r: p[0], g: p[1], b: p[2], a: p[3] }};
    }}"""
    )


def get_scale(page: Page) -> float:
    """Get the current pixels-per-base scale factor."""
    return page.evaluate(
        """() => {
        const canvas = document.getElementById('reads-canvas');
        return canvas.width / (viewer.state.viewport.end - viewer.state.viewport.start);
    }"""
    )


def genomic_to_pixel_x(page: Page, genomic_pos: int) -> int:
    """Convert genomic position to canvas pixel x-coordinate."""
    return page.evaluate(
        f"""() => {{
        const canvas = document.getElementById('reads-canvas');
        const scale = canvas.width / (viewer.state.viewport.end - viewer.state.viewport.start);
        return Math.floor(({genomic_pos} - viewer.state.viewport.start) * scale);
    }}"""
    )


# ── Screenshot-on-failure hook ───────────────────────────────────────

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "test-results")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Capture screenshot on test failure for debugging."""
    outcome = yield
    report = outcome.get_result()
    if report.when == "call" and report.failed:
        # Try to find a Page fixture
        page = item.funcargs.get("viewer_page") or item.funcargs.get("page")
        if page is not None:
            os.makedirs(RESULTS_DIR, exist_ok=True)
            path = os.path.join(RESULTS_DIR, f"{item.name}.png")
            with contextlib.suppress(Exception):
                page.screenshot(path=path)

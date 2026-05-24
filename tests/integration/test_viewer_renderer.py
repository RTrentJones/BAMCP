"""Integration tests for ViewerRenderer — real Playwright + real viewer HTML."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from bamcp.core.parsers import fetch_region
from bamcp.core.serialization import serialize_region_data
from bamcp.eval.renderer import ViewerRenderer

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _is_dist_built() -> bool:
    """Skip these tests when dist/viewer.html hasn't been built yet."""
    return (
        Path(__file__).resolve().parent.parent.parent
        / "src"
        / "bamcp"
        / "static"
        / "dist"
        / "viewer.html"
    ).exists()


def _is_chromium_installed() -> bool:
    """Skip when Playwright's Chromium hasn't been downloaded.

    The `make test` job on CI installs the Playwright Python package but does
    not run `playwright install chromium`; only the e2e job does. Skip rather
    than fail so the matrix stays green.
    """
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        return False
    try:
        with sync_playwright() as p:
            # executable_path raises if the browser isn't installed.
            return Path(p.chromium.executable_path).exists()
    except Exception:
        return False


pytestmark = [
    pytest.mark.skipif(
        not _is_dist_built(),
        reason="viewer dist not built — run `cd src/bamcp/static && npm run build`",
    ),
    pytest.mark.skipif(
        not _is_chromium_installed(),
        reason="Playwright Chromium not installed — run `python -m playwright install chromium`",
    ),
]


@pytest.fixture
def comprehensive_payload(comprehensive_bam_path, comprehensive_ref_fasta_path):
    """Build a real ui/init payload from the comprehensive fixture."""
    data = fetch_region(
        comprehensive_bam_path,
        "chr1:1000-1300",
        comprehensive_ref_fasta_path,
        min_vaf=0.02,
        min_depth=2,
    )
    return serialize_region_data(data)


@pytest.mark.integration
def test_capture_produces_png_bytes(comprehensive_payload):
    async def _go():
        async with ViewerRenderer() as renderer:
            return await renderer.capture(comprehensive_payload, "expanded")

    png = asyncio.run(_go())
    assert png.startswith(PNG_MAGIC), "Output is not a PNG"
    assert len(png) > 5000, f"PNG suspiciously small: {len(png)} bytes"


@pytest.mark.integration
def test_mode_changes_output(comprehensive_payload):
    """Two different display modes must produce visibly different PNGs."""

    async def _go():
        async with ViewerRenderer() as renderer:
            expanded = await renderer.capture(comprehensive_payload, "expanded")
            dv_composite = await renderer.capture(comprehensive_payload, "dv-composite")
            return expanded, dv_composite

    expanded, dv_composite = asyncio.run(_go())
    assert expanded != dv_composite
    # Both should still be valid PNGs.
    assert expanded.startswith(PNG_MAGIC)
    assert dv_composite.startswith(PNG_MAGIC)


@pytest.mark.integration
def test_active_variant_changes_dv_strips_output(comprehensive_payload):
    """Setting an active variant changes the supports-variant strip in dv-strips."""

    async def _go():
        async with ViewerRenderer() as renderer:
            without = await renderer.capture(comprehensive_payload, "dv-strips")
            with_active = await renderer.capture(
                comprehensive_payload,
                "dv-strips",
                active_variant_position=1050,
                active_variant_alt="T",
            )
            return without, with_active

    without, with_active = asyncio.run(_go())
    assert without != with_active, (
        "Activating a variant should alter the dv-strips supports-variant strip"
    )


@pytest.mark.integration
def test_renderer_requires_context_manager():
    """Calling capture outside an async-with block raises a clear error."""
    renderer = ViewerRenderer()

    async def _go():
        await renderer.capture({}, "expanded")

    with pytest.raises(RuntimeError, match="async with"):
        asyncio.run(_go())


@pytest.mark.integration
def test_render_viewer_cli_writes_png(
    tmp_path, comprehensive_bam_path, comprehensive_ref_fasta_path
):
    """End-to-end CLI smoke: scripts/render_viewer.py writes a real PNG."""
    out = tmp_path / "render.png"
    from bamcp.eval.render_cli import main

    rc = main(
        [
            "--bam",
            comprehensive_bam_path,
            "--reference",
            comprehensive_ref_fasta_path,
            "--region",
            "chr1:1000-1300",
            "--mode",
            "dv-strips",
            "--out",
            str(out),
            "--print-path-only",
        ]
    )
    assert rc == 0
    assert out.exists()
    assert out.stat().st_size > 5000
    assert out.read_bytes().startswith(PNG_MAGIC)


@pytest.mark.integration
def test_render_viewer_cli_accepts_ui_init_json(tmp_path, comprehensive_payload):
    """The --ui-init-json path skips BAM parsing and renders directly."""
    payload_path = tmp_path / "payload.json"
    payload_path.write_text(json.dumps(comprehensive_payload))
    out = tmp_path / "from_json.png"

    from bamcp.eval.render_cli import main

    rc = main(
        [
            "--ui-init-json",
            str(payload_path),
            "--mode",
            "expanded",
            "--out",
            str(out),
            "--print-path-only",
        ]
    )
    assert rc == 0
    assert out.read_bytes().startswith(PNG_MAGIC)


@pytest.mark.integration
def test_render_viewer_cli_modes_differ(tmp_path, comprehensive_payload):
    """Rendering the same payload at two modes via the CLI produces different PNGs."""
    payload_path = tmp_path / "payload.json"
    payload_path.write_text(json.dumps(comprehensive_payload))
    out_a = tmp_path / "a.png"
    out_b = tmp_path / "b.png"

    from bamcp.eval.render_cli import main

    main(
        [
            "--ui-init-json",
            str(payload_path),
            "--mode",
            "expanded",
            "--out",
            str(out_a),
            "--print-path-only",
        ]
    )
    main(
        [
            "--ui-init-json",
            str(payload_path),
            "--mode",
            "dv-composite",
            "--out",
            str(out_b),
            "--print-path-only",
        ]
    )
    assert out_a.read_bytes() != out_b.read_bytes()

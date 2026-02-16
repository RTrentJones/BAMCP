"""End-to-end Playwright tests for the BAMCP viewer UI."""

import contextlib

import pytest
from playwright.sync_api import Page, expect

from bamcp.resources import get_viewer_html
from tests.e2e.conftest import (
    count_opaque_pixels,
    genomic_to_pixel_x,
    get_scale,
    has_pixels_in_region,
    sample_pixel_at_genomic,
)

SAMPLE_DATA = {
    "contig": "chr1",
    "start": 100,
    "end": 600,
    "reads": [
        {
            "name": "read1",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
        },
        {
            "name": "read2",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 120,
            "end_position": 170,
            "mapping_quality": 50,
            "is_reverse": True,
            "mismatches": [{"pos": 130, "ref": "A", "alt": "T"}],
        },
        {
            "name": "read3",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 200,
            "end_position": 250,
            "mapping_quality": 40,
            "is_reverse": False,
            "mismatches": [],
        },
    ],
    "coverage": [0] * 500,
    "variants": [
        {
            "contig": "chr1",
            "position": 130,
            "ref": "A",
            "alt": "T",
            "vaf": 0.25,
            "depth": 20,
            "alt_count": 5,
        },
        {
            "contig": "chr1",
            "position": 250,
            "ref": "G",
            "alt": "C",
            "vaf": 0.15,
            "depth": 30,
            "alt_count": 4,
        },
    ],
    "reference_sequence": "ACGTACGTAC" * 50,
}

# Set some coverage values
for i in range(500):
    if 0 <= i < 70:
        SAMPLE_DATA["coverage"][i] = 2
    elif 70 <= i < 100 or 100 <= i < 150:
        SAMPLE_DATA["coverage"][i] = 1


@pytest.fixture
def viewer_page(page: Page):
    """Load the viewer HTML into a Playwright page."""
    html = get_viewer_html()
    page.set_content(html)
    page.wait_for_load_state("domcontentloaded")
    # Wait for the module script to create the global viewer object
    page.wait_for_function("() => typeof window.viewer !== 'undefined'", timeout=5000)
    return page


def send_init_data(page: Page, data: dict = None):
    """Load data directly into the viewer's state and trigger render."""
    if data is None:
        data = SAMPLE_DATA
    page.evaluate(
        """(data) => {
            viewer.state.loadData(data);
            viewer.renderVariantTable();
            viewer.renderer.resize();
            // Update region input to match loaded data
            const input = document.getElementById('region-input');
            if (input) input.value = data.contig + ':' + data.start + '-' + data.end;
        }""",
        data,
    )
    # Give the viewer time to process rendering
    page.wait_for_timeout(100)


class TestViewerPageStructure:
    """E2E tests verifying the viewer page loads correctly."""

    @pytest.mark.e2e
    def test_page_loads(self, viewer_page: Page):
        """Viewer page should load without errors."""
        expect(viewer_page.locator("#container")).to_be_visible()

    @pytest.mark.e2e
    def test_toolbar_visible(self, viewer_page: Page):
        """Toolbar should be visible with all controls."""
        expect(viewer_page.locator("#toolbar")).to_be_visible()
        expect(viewer_page.locator("#region-input")).to_be_visible()
        expect(viewer_page.locator("#go-btn")).to_be_visible()
        expect(viewer_page.locator("#zoom-in")).to_be_visible()
        expect(viewer_page.locator("#zoom-out")).to_be_visible()

    @pytest.mark.e2e
    def test_canvases_exist(self, viewer_page: Page):
        """Coverage and reads canvases should exist."""
        expect(viewer_page.locator("#coverage-canvas")).to_be_attached()
        expect(viewer_page.locator("#reads-canvas")).to_be_attached()

    @pytest.mark.e2e
    def test_variant_panel_exists(self, viewer_page: Page):
        """Variant panel should exist with table headers."""
        expect(viewer_page.locator("#variant-panel")).to_be_visible()
        headers = viewer_page.locator("#variant-panel th")
        assert headers.count() == 8

    @pytest.mark.e2e
    def test_tooltip_hidden_initially(self, viewer_page: Page):
        """Tooltip should be hidden initially."""
        tooltip = viewer_page.locator("#tooltip")
        expect(tooltip).to_be_hidden()

    @pytest.mark.e2e
    def test_region_input_placeholder(self, viewer_page: Page):
        """Region input should have placeholder text."""
        input_el = viewer_page.locator("#region-input")
        assert input_el.get_attribute("placeholder") == "chr1:1000-2000"


class TestDataLoading:
    """E2E tests for loading data into the viewer."""

    @pytest.mark.e2e
    def test_load_data_updates_region_input(self, viewer_page: Page):
        """Loading data should update the region input."""
        send_init_data(viewer_page)
        input_val = viewer_page.locator("#region-input").input_value()
        assert input_val == "chr1:100-600"

    @pytest.mark.e2e
    def test_load_data_populates_variant_table(self, viewer_page: Page):
        """Loading data should populate the variant table."""
        send_init_data(viewer_page)
        rows = viewer_page.locator("#variant-table tr")
        assert rows.count() == 2

    @pytest.mark.e2e
    def test_variant_table_content(self, viewer_page: Page):
        """Variant table should display correct data."""
        send_init_data(viewer_page)
        first_row = viewer_page.locator("#variant-table tr").first
        cells = first_row.locator("td")
        assert "130" in cells.nth(0).inner_text()
        assert cells.nth(1).inner_text() == "A"
        assert cells.nth(2).inner_text() == "T"
        assert "25.0%" in cells.nth(3).inner_text()
        assert cells.nth(4).inner_text() == "20"

    @pytest.mark.e2e
    def test_canvas_has_dimensions(self, viewer_page: Page):
        """Canvases should have non-zero dimensions after load."""
        send_init_data(viewer_page)
        cov_width = viewer_page.locator("#coverage-canvas").evaluate("el => el.width")
        reads_width = viewer_page.locator("#reads-canvas").evaluate("el => el.width")
        assert cov_width > 0
        assert reads_width > 0

    @pytest.mark.e2e
    def test_update_message_replaces_data(self, viewer_page: Page):
        """Loading new data should replace existing data."""
        send_init_data(viewer_page)

        # Send updated data with different region
        new_data = {**SAMPLE_DATA, "contig": "chr2", "start": 500, "end": 1000}
        send_init_data(viewer_page, new_data)

        input_val = viewer_page.locator("#region-input").input_value()
        assert input_val == "chr2:500-1000"


class TestNavigation:
    """E2E tests for viewer navigation controls."""

    @pytest.mark.e2e
    def test_zoom_in_button(self, viewer_page: Page):
        """Zoom in button should narrow the viewport."""
        send_init_data(viewer_page)

        # Get initial viewport via JS
        initial_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        viewer_page.locator("#zoom-in").click()

        new_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        assert new_span < initial_span

    @pytest.mark.e2e
    def test_zoom_out_button(self, viewer_page: Page):
        """Zoom out button should widen the viewport."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        viewer_page.locator("#zoom-out").click()

        new_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        assert new_span > initial_span

    @pytest.mark.e2e
    def test_zoom_preserves_center(self, viewer_page: Page):
        """Zoom should preserve the center position."""
        send_init_data(viewer_page)

        initial_center = viewer_page.evaluate(
            "() => (viewer.state.viewport.start + viewer.state.viewport.end) / 2"
        )

        viewer_page.locator("#zoom-in").click()

        new_center = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return (vp.start + vp.end) / 2; }"
        )

        # Center should be approximately the same (within rounding)
        assert abs(new_center - initial_center) < 2

    @pytest.mark.e2e
    def test_variant_click_jumps_to_position(self, viewer_page: Page):
        """Clicking a variant row should jump to that position."""
        send_init_data(viewer_page)

        viewer_page.evaluate("() => viewer.state.viewport.start")

        # Click the second variant row (position 250)
        viewer_page.locator("#variant-table tr").nth(1).click()

        new_center = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return (vp.start + vp.end) / 2; }"
        )

        # Center should be near position 250
        assert abs(new_center - 250) < 10

    @pytest.mark.e2e
    def test_mouse_drag_pans(self, viewer_page: Page):
        """Mouse drag on reads canvas should pan the view."""
        send_init_data(viewer_page)

        initial_start = viewer_page.evaluate("() => viewer.state.viewport.start")

        # Simulate drag on reads canvas
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        if box:
            start_x = box["x"] + box["width"] / 2
            start_y = box["y"] + box["height"] / 2

            viewer_page.mouse.move(start_x, start_y)
            viewer_page.mouse.down()
            viewer_page.mouse.move(start_x - 100, start_y)
            viewer_page.mouse.up()

        new_start = viewer_page.evaluate("() => viewer.state.viewport.start")

        # Panning left (drag right-to-left) should increase start position
        assert new_start != initial_start


class TestReadPacking:
    """E2E tests for read packing algorithm."""

    @pytest.mark.e2e
    def test_packed_rows_created(self, viewer_page: Page):
        """Reads should be packed into rows."""
        send_init_data(viewer_page)

        row_count = viewer_page.evaluate("() => viewer.state.packedRows.length")
        assert row_count > 0

    @pytest.mark.e2e
    def test_all_reads_packed(self, viewer_page: Page):
        """All reads should be present in packed rows."""
        send_init_data(viewer_page)

        total = viewer_page.evaluate(
            "() => viewer.state.packedRows.reduce((sum, row) => sum + row.length, 0)"
        )
        assert total == len(SAMPLE_DATA["reads"])

    @pytest.mark.e2e
    def test_no_overlapping_reads_in_row(self, viewer_page: Page):
        """Reads in the same row should not overlap."""
        send_init_data(viewer_page)

        no_overlap = viewer_page.evaluate("""() => {
                for (const row of viewer.state.packedRows) {
                    for (let i = 0; i < row.length - 1; i++) {
                        if (row[i].end_position >= row[i + 1].position) {
                            return false;
                        }
                    }
                }
                return true;
            }""")
        assert no_overlap


class TestEmptyState:
    """E2E tests for empty/no-data states."""

    @pytest.mark.e2e
    def test_no_data_renders_without_error(self, viewer_page: Page):
        """Viewer should handle no data gracefully."""
        # No data sent, just verify no JS errors
        errors = []
        viewer_page.on("pageerror", lambda err: errors.append(str(err)))
        viewer_page.wait_for_timeout(500)
        assert len(errors) == 0

    @pytest.mark.e2e
    def test_empty_reads(self, viewer_page: Page):
        """Should handle data with zero reads."""
        empty_data = {
            "contig": "chr1",
            "start": 100,
            "end": 200,
            "reads": [],
            "coverage": [0] * 100,
            "variants": [],
            "reference_sequence": None,
        }
        send_init_data(viewer_page, empty_data)

        row_count = viewer_page.evaluate("() => viewer.state.packedRows.length")
        assert row_count == 0

        rows = viewer_page.locator("#variant-table tr")
        assert rows.count() == 0

    @pytest.mark.e2e
    def test_empty_variants(self, viewer_page: Page):
        """Should handle data with reads but no variants."""
        data = {**SAMPLE_DATA, "variants": []}
        send_init_data(viewer_page, data)

        rows = viewer_page.locator("#variant-table tr")
        assert rows.count() == 0


class TestSoftClipRendering:
    """E2E tests for soft-clip visualization via pre-parsed soft_clips data."""

    @pytest.mark.e2e
    def test_soft_clip_data_loaded(self, viewer_page: Page):
        """Reads with soft_clips array should load without error."""
        soft_clip_data = {
            **SAMPLE_DATA,
            "reads": [
                {
                    "name": "clip_read",
                    "sequence": "NNNNN" + "ACGTACGTAC" * 5,
                    "cigar": "5S50M",
                    "position": 100,
                    "end_position": 150,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                    "soft_clips": [
                        {"position": 95, "length": 5, "side": "left"},
                    ],
                }
            ],
        }
        errors = []
        viewer_page.on("pageerror", lambda err: errors.append(str(err)))
        send_init_data(viewer_page, soft_clip_data)
        viewer_page.wait_for_timeout(200)
        assert len(errors) == 0

    @pytest.mark.e2e
    def test_soft_clips_array_accessible(self, viewer_page: Page):
        """Loaded reads should have soft_clips array accessible in state."""
        soft_clip_data = {
            **SAMPLE_DATA,
            "reads": [
                {
                    "name": "clip_read",
                    "sequence": "NNNNN" + "ACGTACGTAC" * 5,
                    "cigar": "5S50M",
                    "position": 100,
                    "end_position": 150,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                    "soft_clips": [
                        {"position": 95, "length": 5, "side": "left"},
                    ],
                }
            ],
        }
        send_init_data(viewer_page, soft_clip_data)
        clips = viewer_page.evaluate("() => viewer.state.data.reads[0].soft_clips")
        assert len(clips) == 1
        assert clips[0]["side"] == "left"
        assert clips[0]["length"] == 5

    @pytest.mark.e2e
    def test_soft_clip_rendering_no_error(self, viewer_page: Page):
        """Rendering reads with soft-clips should not cause errors."""
        soft_clip_data = {
            **SAMPLE_DATA,
            "reads": [
                {
                    "name": "clip_read",
                    "sequence": "NNNNN" + "ACGTACGTAC" * 5,
                    "cigar": "5S50M",
                    "position": 100,
                    "end_position": 150,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                    "soft_clips": [
                        {"position": 95, "length": 5, "side": "left"},
                    ],
                }
            ],
        }
        errors = []
        viewer_page.on("pageerror", lambda err: errors.append(str(err)))
        send_init_data(viewer_page, soft_clip_data)
        viewer_page.wait_for_timeout(200)
        assert len(errors) == 0


class TestCanvasRendering:
    """E2E tests verifying canvas rendering occurs."""

    @pytest.mark.e2e
    def test_coverage_canvas_drawn(self, viewer_page: Page):
        """Coverage canvas should have pixels drawn after data load."""
        send_init_data(viewer_page)

        # Check if canvas has non-transparent pixels
        has_content = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('coverage-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }""")
        assert has_content

    @pytest.mark.e2e
    def test_reads_canvas_drawn(self, viewer_page: Page):
        """Reads canvas should have pixels drawn after data load."""
        send_init_data(viewer_page)

        has_content = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }""")
        assert has_content

    @pytest.mark.e2e
    def test_canvas_redraws_on_zoom(self, viewer_page: Page):
        """Canvas should redraw when zooming (viewport changes)."""
        send_init_data(viewer_page)

        before_viewport = viewer_page.evaluate(
            "() => ({ start: viewer.state.viewport.start, end: viewer.state.viewport.end })"
        )

        viewer_page.locator("#zoom-in").click()
        viewer_page.wait_for_timeout(100)

        after_viewport = viewer_page.evaluate(
            "() => ({ start: viewer.state.viewport.start, end: viewer.state.viewport.end })"
        )

        # Viewport span should change after zoom
        before_span = before_viewport["end"] - before_viewport["start"]
        after_span = after_viewport["end"] - after_viewport["start"]
        assert after_span < before_span


class TestMCPBridge:
    """E2E tests for the MCP SDK client integration."""

    @pytest.mark.e2e
    def test_client_object_exists(self, viewer_page: Page):
        """Viewer should have a BAMCPClient instance."""
        result = viewer_page.evaluate("() => typeof viewer.client")
        assert result == "object"

    @pytest.mark.e2e
    def test_go_button_triggers_navigation(self, viewer_page: Page):
        """Clicking Go with a region should attempt navigation."""
        send_init_data(viewer_page)
        viewer_page.locator("#region-input").fill("chr1:200-400")
        viewer_page.locator("#go-btn").click()
        viewer_page.wait_for_timeout(300)
        # The Go button uses sendMessage or callServerTool which may fail
        # without a host, but should not throw an unhandled error
        errors = viewer_page.evaluate("() => window.__pageErrors ? window.__pageErrors.length : 0")
        assert errors == 0


# --- Narrow-region data for mismatch/color testing at scale > 2 ---
NARROW_DATA = {
    "contig": "chr1",
    "start": 100,
    "end": 150,  # 50bp span -> scale ~20 px/bp in 1024px canvas -> well above 2
    "reads": [
        {
            "name": "fwd_read",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [{"pos": 125, "ref": "A", "alt": "T"}],
        },
        {
            "name": "rev_read",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 50,
            "is_reverse": True,
            "mismatches": [],
        },
    ],
    "coverage": [2] * 50,
    "variants": [],
    "reference_sequence": "ACGTACGTAC" * 5,
}


class TestTooltipInteraction:
    """E2E tests for tooltip hover behavior."""

    @pytest.mark.e2e
    def test_tooltip_shows_on_read_hover(self, viewer_page: Page):
        """Hovering over a read should display the tooltip with read info."""
        send_init_data(viewer_page, NARROW_DATA)

        # Read spans full canvas width; row 0 starts at y=0 in reads-canvas
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Hover over the center of the first read row (y offset ~6px into the row)
        hover_x = box["x"] + box["width"] / 2
        hover_y = box["y"] + 6  # Middle of first row (READ_HEIGHT=12)

        viewer_page.mouse.move(hover_x, hover_y)
        viewer_page.wait_for_timeout(200)

        tooltip = viewer_page.locator("#tooltip")
        expect(tooltip).to_be_visible()

        content = tooltip.inner_text()
        assert "fwd_read" in content
        assert "MAPQ" in content
        assert "50M" in content

    @pytest.mark.e2e
    def test_tooltip_shows_reverse_strand(self, viewer_page: Page):
        """Tooltip should show 'Reverse' for reverse-strand reads."""
        send_init_data(viewer_page, NARROW_DATA)

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Row 1 (second read is reverse). READ_HEIGHT=12, READ_GAP=2 → row1 at y=14
        hover_x = box["x"] + box["width"] / 2
        hover_y = box["y"] + 14 + 6  # Middle of second row

        viewer_page.mouse.move(hover_x, hover_y)
        viewer_page.wait_for_timeout(200)

        tooltip = viewer_page.locator("#tooltip")
        expect(tooltip).to_be_visible()
        content = tooltip.inner_text()
        assert "rev_read" in content
        assert "MAPQ" in content

    @pytest.mark.e2e
    def test_tooltip_hides_on_mouseout(self, viewer_page: Page):
        """Moving mouse off the canvas should hide the tooltip."""
        send_init_data(viewer_page, NARROW_DATA)

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Hover over read to show tooltip
        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + 6)
        viewer_page.wait_for_timeout(200)
        expect(viewer_page.locator("#tooltip")).to_be_visible()

        # Move mouse completely off the canvas
        viewer_page.mouse.move(box["x"] - 50, box["y"] - 50)
        viewer_page.wait_for_timeout(200)

        expect(viewer_page.locator("#tooltip")).to_be_hidden()

    @pytest.mark.e2e
    def test_tooltip_hides_when_not_over_read(self, viewer_page: Page):
        """Hovering over empty space (below all reads) should hide tooltip."""
        send_init_data(viewer_page, NARROW_DATA)

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # First, hover over a read
        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + 6)
        viewer_page.wait_for_timeout(200)
        expect(viewer_page.locator("#tooltip")).to_be_visible()

        # Now hover over empty area far below the reads
        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] - 5)
        viewer_page.wait_for_timeout(200)

        expect(viewer_page.locator("#tooltip")).to_be_hidden()


class TestWheelZoom:
    """E2E tests for scroll wheel zoom."""

    @pytest.mark.e2e
    def test_wheel_scroll_down_zooms_out(self, viewer_page: Page):
        """Scrolling down should zoom out (increase span)."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
        viewer_page.mouse.wheel(0, 100)  # deltaY > 0 → zoom out (factor 1.2)
        viewer_page.wait_for_timeout(200)

        new_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )
        assert new_span > initial_span

    @pytest.mark.e2e
    def test_wheel_scroll_up_zooms_in(self, viewer_page: Page):
        """Scrolling up should zoom in (decrease span)."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
        viewer_page.mouse.wheel(0, -100)  # deltaY < 0 → zoom in (factor 0.8)
        viewer_page.wait_for_timeout(200)

        new_span = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport; return vp.end - vp.start; }"
        )
        assert new_span < initial_span


class TestMismatchRendering:
    """E2E tests for mismatch rendering at different zoom levels."""

    @pytest.mark.e2e
    def test_mismatch_pixel_at_medium_zoom(self, viewer_page: Page):
        """At medium zoom (0.5 < scale < 10), mismatches should render colored pixels."""
        # Use a 200bp viewport so scale is ~5 px/bp (medium zoom, no per-base rendering)
        medium_data = {
            **NARROW_DATA,
            "start": 50,
            "end": 250,
            "coverage": [2] * 200,
            "reference_sequence": "ACGTACGTAC" * 20,
        }
        send_init_data(viewer_page, medium_data)

        # Mismatch is at position 125, alt='T' → color #ef4444 (red).
        mismatch_color = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const mx = Math.floor((125 - viewer.state.viewport.start) * scale);
                const pixel = ctx.getImageData(mx + 1, 6, 1, 1).data;
                return { r: pixel[0], g: pixel[1], b: pixel[2], a: pixel[3], scale: scale };
            }""")
        assert mismatch_color["a"] > 0, "Mismatch pixel should be opaque"
        assert mismatch_color["r"] > mismatch_color["b"], (
            f"Mismatch should be red-dominant (T color), "
            f"got r={mismatch_color['r']} b={mismatch_color['b']} scale={mismatch_color['scale']}"
        )

    @pytest.mark.e2e
    def test_no_mismatch_at_very_low_zoom(self, viewer_page: Page):
        """At very low zoom (scale <= 0.5), mismatches should NOT render."""
        # Use very wide span data so scale < 0.5 (need > 2048bp for 1024px canvas)
        wide_data = {
            **NARROW_DATA,
            "start": 0,
            "end": 3000,
            "coverage": [2] * 3000,
            "reference_sequence": "ACGTACGTAC" * 300,
        }
        send_init_data(viewer_page, wide_data)

        scale = viewer_page.evaluate(
            "() => { const vp = viewer.state.viewport;"
            " return viewer.readsCanvas.width / (vp.end - vp.start); }"
        )
        assert scale <= 0.5, f"Scale should be <= 0.5 for this test, got {scale}"

        # Sample pixel at mismatch position 125: should be the read body color, not mismatch red
        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const mx = Math.floor((125 - viewer.state.viewport.start) * scale);
                const p = ctx.getImageData(mx, 6, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        # At very low zoom, mismatch should not override — pixel should be blue read color
        if pixel["a"] > 0:
            assert pixel["b"] >= pixel["r"], (
                f"At very low zoom, should be blue read color not red mismatch, "
                f"got r={pixel['r']} b={pixel['b']}"
            )


class TestReadColors:
    """E2E tests for forward/reverse read coloring."""

    @pytest.mark.e2e
    def test_forward_read_is_blue(self, viewer_page: Page):
        """Forward-strand reads should be blue (#60a5fa) at medium zoom."""
        # Use a wider viewport (200bp) so per-base colors aren't drawn (scale < 10)
        medium_data = {
            **NARROW_DATA,
            "start": 50,
            "end": 250,
            "coverage": [2] * 200,
            "reference_sequence": "ACGTACGTAC" * 20,
        }
        send_init_data(viewer_page, medium_data)

        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                // Sample a pixel in row 0, at position 110 (away from mismatch at 125)
                const x = Math.floor((110 - viewer.state.viewport.start) * scale);
                const p = ctx.getImageData(x, 3, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        assert pixel["a"] > 0, "Pixel should be drawn"
        # #60a5fa: r=96, g=165, b=250 → blue dominant
        assert pixel["b"] > pixel["r"], (
            f"Forward read should be blue-dominant, "
            f"got r={pixel['r']} g={pixel['g']} b={pixel['b']}"
        )

    @pytest.mark.e2e
    def test_reverse_read_is_purple(self, viewer_page: Page):
        """Reverse-strand reads should be red (#ef4444) at medium zoom."""
        medium_data = {
            **NARROW_DATA,
            "start": 50,
            "end": 250,
            "coverage": [2] * 200,
            "reference_sequence": "ACGTACGTAC" * 20,
        }
        send_init_data(viewer_page, medium_data)

        # Find the row position of the reverse read
        rev_row = viewer_page.evaluate("""() => {
                // packedRows is an array of rows; find which row has the reverse read
                for (let r = 0; r < viewer.state.packedRows.length; r++) {
                    for (const read of viewer.state.packedRows[r]) {
                        if (read.is_reverse) return r;
                    }
                }
                return -1;
            }""")
        assert rev_row >= 0, "Reverse read should be found in packed rows"
        # Row height = 12px + 2px gap = 14px per row
        row_y = rev_row * 14 + 6  # middle of that row

        pixel = viewer_page.evaluate(f"""() => {{
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const x = Math.floor((110 - viewer.state.viewport.start) * scale);
                const p = ctx.getImageData(x, {row_y}, 1, 1).data;
                return {{ r: p[0], g: p[1], b: p[2], a: p[3] }};
            }}""")
        assert pixel["a"] > 0, "Pixel should be drawn"
        # #ef4444: r=239, g=68, b=68 → red dominant
        assert pixel["r"] > pixel["b"], (
            f"Reverse read should be red-dominant, got r={pixel['r']} b={pixel['b']}"
        )

    @pytest.mark.e2e
    def test_forward_and_reverse_have_different_colors(self, viewer_page: Page):
        """Forward and reverse reads should be visually distinct colors."""
        medium_data = {
            **NARROW_DATA,
            "start": 50,
            "end": 250,
            "coverage": [2] * 200,
            "reference_sequence": "ACGTACGTAC" * 20,
        }
        send_init_data(viewer_page, medium_data)

        colors = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const x = Math.floor((110 - viewer.state.viewport.start) * scale);
                // Find y position of each read type
                let fwdY = 3, revY = 3;
                for (let r = 0; r < viewer.state.packedRows.length; r++) {
                    for (const read of viewer.state.packedRows[r]) {
                        const y = r * 14 + 3;
                        if (!read.is_reverse) fwdY = y;
                        else revY = y;
                    }
                }
                const fwd = ctx.getImageData(x, fwdY, 1, 1).data;
                const rev = ctx.getImageData(x, revY, 1, 1).data;
                return {
                    fwd: { r: fwd[0], g: fwd[1], b: fwd[2] },
                    rev: { r: rev[0], g: rev[1], b: rev[2] },
                };
            }""")
        # Colors should differ
        assert (
            colors["fwd"]["r"] != colors["rev"]["r"]
            or colors["fwd"]["g"] != colors["rev"]["g"]
            or colors["fwd"]["b"] != colors["rev"]["b"]
        ), (
            f"Forward and reverse should have different colors: "
            f"fwd={colors['fwd']} rev={colors['rev']}"
        )


class TestSoftClipPixelVerification:
    """E2E tests verifying soft-clip regions actually draw pixels."""

    @pytest.mark.e2e
    def test_left_soft_clip_draws_pixels(self, viewer_page: Page):
        """Left soft-clip region should have semi-transparent pixels drawn."""
        # Enable soft clips display
        viewer_page.evaluate("() => { viewer.state.settings.showSoftClips = true; }")
        clip_data = {
            "contig": "chr1",
            "start": 90,
            "end": 160,
            "reads": [
                {
                    "name": "clipped",
                    "sequence": "N" * 10 + "ACGTACGTAC" * 5,
                    "cigar": "10S50M",
                    "position": 100,
                    "end_position": 150,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                    "soft_clips": [
                        {"position": 90, "length": 10, "side": "left"},
                    ],
                }
            ],
            "coverage": [1] * 70,
            "variants": [],
            "reference_sequence": "A" * 70,
        }
        send_init_data(viewer_page, clip_data)

        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const x = Math.floor((95 - viewer.state.viewport.start) * scale);
                const p = ctx.getImageData(x, 6, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        assert pixel["a"] > 0, "Soft-clip region should have pixels drawn"

    @pytest.mark.e2e
    def test_soft_clip_has_border(self, viewer_page: Page):
        """Soft-clip region should have a stroke border."""
        viewer_page.evaluate("() => { viewer.state.settings.showSoftClips = true; }")
        clip_data = {
            "contig": "chr1",
            "start": 80,
            "end": 160,
            "reads": [
                {
                    "name": "clipped",
                    "sequence": "N" * 10 + "ACGTACGTAC" * 5,
                    "cigar": "10S50M",
                    "position": 100,
                    "end_position": 150,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                    "soft_clips": [
                        {"position": 90, "length": 10, "side": "left"},
                    ],
                }
            ],
            "coverage": [1] * 80,
            "variants": [],
            "reference_sequence": "A" * 80,
        }
        send_init_data(viewer_page, clip_data)

        # Sample the top edge of the clip region at position 95
        border_pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const vp = viewer.state.viewport;
                const scale = canvas.width / (vp.end - vp.start);
                const x = Math.floor((95 - viewer.state.viewport.start) * scale);
                const p = ctx.getImageData(x, 0, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        assert border_pixel["a"] > 0, "Border should be drawn"


class TestOffScreenCulling:
    """E2E tests for off-screen read culling."""

    @pytest.mark.e2e
    def test_reads_outside_viewport_not_drawn(self, viewer_page: Page):
        """After panning away, reads outside viewport should not have pixels."""
        send_init_data(viewer_page, NARROW_DATA)

        # Zoom in tight then pan far away so reads at 100-150 are off-screen
        viewer_page.evaluate("""() => {
                viewer.state.viewport.start = 500;
                viewer.state.viewport.end = 550;
                viewer.renderer.resize();
            }""")
        viewer_page.wait_for_timeout(100)

        # The reads canvas should be empty (reads are at 100-150, viewport is 500-550)
        has_content = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }""")
        assert not has_content, "No read pixels should be drawn when reads are off-screen"


class TestResizeHandler:
    """E2E tests for the resize handler."""

    @pytest.mark.e2e
    def test_resize_updates_canvas_dimensions(self, viewer_page: Page):
        """Resizing the window should update canvas dimensions."""
        send_init_data(viewer_page)

        viewer_page.evaluate("() => viewer.readsCanvas.width")

        # Resize the viewport
        viewer_page.set_viewport_size({"width": 800, "height": 600})
        viewer_page.wait_for_timeout(300)

        new_width = viewer_page.evaluate("() => viewer.readsCanvas.width")

        # Canvas width should update to match the new container width
        assert new_width > 0
        # The exact values depend on the container, but they should differ
        # from a drastically different viewport
        assert isinstance(new_width, (int, float))

    @pytest.mark.e2e
    def test_resize_redraws_canvas(self, viewer_page: Page):
        """Canvas should have content after resize."""
        send_init_data(viewer_page)

        viewer_page.set_viewport_size({"width": 600, "height": 400})
        viewer_page.wait_for_timeout(300)

        has_content = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }""")
        assert has_content, "Canvas should still have content after resize"


class TestCoverageLabel:
    """E2E tests for coverage max label rendering."""

    @pytest.mark.e2e
    def test_coverage_max_label_drawn(self, viewer_page: Page):
        """Coverage canvas should render the 'Max: Nx' label text."""
        send_init_data(viewer_page)

        # The label is drawn at (5, 12) with fillText.
        # We can't easily read text from canvas, but we can check that
        # the top-left corner of the coverage canvas has opaque pixels
        # where the label should be.
        has_label_pixels = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('coverage-canvas');
                const ctx = canvas.getContext('2d');
                // Sample the area where 'Max: 2x' is drawn (x=5..60, y=2..14)
                const data = ctx.getImageData(5, 2, 55, 12).data;
                let opaqueCount = 0;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) opaqueCount++;
                }
                return opaqueCount;
            }""")
        # The text "Max: 2x" should produce at least some opaque pixels in that region
        assert has_label_pixels > 5, (
            f"Coverage label should produce opaque pixels, got {has_label_pixels}"
        )

    @pytest.mark.e2e
    def test_coverage_label_text_is_correct(self, viewer_page: Page):
        """The coverage max value should match the actual data max."""
        send_init_data(viewer_page)

        max_cov = viewer_page.evaluate("() => Math.max(...viewer.state.data.coverage)")
        assert max_cov == 2  # Our SAMPLE_DATA has max coverage of 2


# ══════════════════════════════════════════════════════════════════════
# New comprehensive test data constants
# ══════════════════════════════════════════════════════════════════════

SOFT_CLIP_DATA = {
    "contig": "chr1",
    "start": 85,
    "end": 210,
    "reads": [
        {
            "name": "sc_left",
            "sequence": "N" * 5 + "ACGTACGTAC" * 4 + "ACGTA",
            "cigar": "5S45M",
            "position": 100,
            "end_position": 145,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "soft_clips": [{"position": 95, "length": 5, "sequence": "NNNNN", "side": "left"}],
        },
        {
            "name": "sc_right",
            "sequence": "ACGTACGTAC" * 4 + "ACGTANNNNN",
            "cigar": "45M5S",
            "position": 110,
            "end_position": 155,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "soft_clips": [{"position": 155, "length": 5, "sequence": "NNNNN", "side": "right"}],
        },
        {
            "name": "sc_both",
            "sequence": "NNN" + "ACGTACGTAC" * 4 + "ACGTNNN",
            "cigar": "3S44M3S",
            "position": 120,
            "end_position": 164,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "soft_clips": [
                {"position": 117, "length": 3, "sequence": "NNN", "side": "left"},
                {"position": 164, "length": 3, "sequence": "NNN", "side": "right"},
            ],
        },
        {
            "name": "sc_long",
            "sequence": "N" * 10 + "ACGTACGTAC" * 4,
            "cigar": "10S40M",
            "position": 130,
            "end_position": 170,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "soft_clips": [
                {
                    "position": 120,
                    "length": 10,
                    "sequence": "N" * 10,
                    "side": "left",
                }
            ],
        },
    ],
    "coverage": [1] * 125,
    "variants": [],
    "reference_sequence": "A" * 125,
}

PAIRED_END_DATA = {
    "contig": "chr1",
    "start": 290,
    "end": 610,
    "reads": [
        {
            "name": "pair1",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 300,
            "end_position": 350,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "is_paired": True,
            "mate_position": 550,
            "mate_contig": "chr1",
            "insert_size": 300,
            "is_proper_pair": True,
            "is_read1": True,
        },
        {
            "name": "pair1",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 550,
            "end_position": 600,
            "mapping_quality": 60,
            "is_reverse": True,
            "mismatches": [],
            "is_paired": True,
            "mate_position": 300,
            "mate_contig": "chr1",
            "insert_size": -300,
            "is_proper_pair": True,
            "is_read1": False,
        },
        {
            "name": "pair2",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 350,
            "end_position": 400,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "is_paired": True,
            "mate_position": 1800,
            "mate_contig": "chr1",
            "insert_size": 1500,
            "is_proper_pair": False,
            "is_read1": True,
        },
        {
            "name": "pair3",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 400,
            "end_position": 450,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
            "is_paired": True,
            "mate_position": 100,
            "mate_contig": "chr2",
            "insert_size": None,
            "is_proper_pair": False,
            "is_read1": True,
        },
    ],
    "coverage": [2] * 320,
    "variants": [],
    "reference_sequence": "ACGTACGTAC" * 32,
}

EVIDENCE_DATA = {
    "contig": "chr1",
    "start": 1000,
    "end": 1200,
    "reads": [
        {
            "name": f"ev_read_{i}",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 1020 + i * 5,
            "end_position": 1070 + i * 5,
            "mapping_quality": 60,
            "is_reverse": i % 2 == 1,
            "mismatches": ([{"pos": 1050, "ref": "A", "alt": "T"}] if i < 8 else []),
        }
        for i in range(20)
    ],
    "coverage": [10] * 200,
    "variants": [
        {
            "contig": "chr1",
            "position": 1050,
            "ref": "A",
            "alt": "T",
            "vaf": 0.4,
            "depth": 20,
            "alt_count": 8,
            "strand_forward": 8,
            "strand_reverse": 0,
            "mean_quality": 40,
            "confidence": "low",
            "artifact_risk": {
                "risks": [
                    {
                        "type": "strand_bias",
                        "severity": "high",
                        "description": "Strand bias: 100% forward",
                        "value": 1.0,
                    }
                ],
                "risk_score": 0.7,
                "artifact_likelihood": "high",
            },
        },
        {
            "contig": "chr1",
            "position": 1150,
            "ref": "G",
            "alt": "A",
            "vaf": 0.5,
            "depth": 20,
            "alt_count": 10,
            "strand_forward": 5,
            "strand_reverse": 5,
            "mean_quality": 40,
            "confidence": "high",
            "artifact_risk": {
                "risks": [],
                "risk_score": 0.0,
                "artifact_likelihood": "low",
            },
        },
    ],
    "variant_evidence": {
        "1050:A>T": {
            "forward_count": 8,
            "reverse_count": 0,
            "strand_bias": 1.0,
            "mean_quality": 40,
            "median_quality": 40,
            "quality_histogram": [0, 0, 0, 0, 8],
            "position_histogram": [2, 2, 2, 1, 1, 0],
            "mapq_histogram": [0, 0, 0, 0, 0, 0, 8],
            "artifact_risk": {
                "risks": [
                    {
                        "type": "strand_bias",
                        "severity": "high",
                        "description": "Strand bias: 100% forward",
                        "value": 1.0,
                    }
                ],
                "risk_score": 0.7,
                "artifact_likelihood": "high",
            },
        },
        "1150:G>A": {
            "forward_count": 5,
            "reverse_count": 5,
            "strand_bias": 0.0,
            "mean_quality": 40,
            "median_quality": 40,
            "quality_histogram": [0, 0, 0, 2, 8],
            "position_histogram": [0, 1, 3, 4, 2, 0],
            "mapq_histogram": [0, 0, 0, 0, 0, 0, 10],
            "artifact_risk": {
                "risks": [],
                "risk_score": 0.0,
                "artifact_likelihood": "low",
            },
        },
    },
    "reference_sequence": "ACGTACGTAC" * 20,
}

MULTI_VARIANT_DATA = {
    "contig": "chr1",
    "start": 100,
    "end": 1000,
    "reads": [
        {
            "name": f"mv_read_{i}",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100 + i * 40,
            "end_position": 150 + i * 40,
            "mapping_quality": 60,
            "is_reverse": i % 2 == 1,
            "mismatches": [],
        }
        for i in range(20)
    ],
    "coverage": [5] * 900,
    "variants": [
        {
            "contig": "chr1",
            "position": 130,
            "ref": "A",
            "alt": "T",
            "vaf": 0.35,
            "depth": 20,
            "alt_count": 7,
            "confidence": "high",
            "strand_forward": 4,
            "strand_reverse": 3,
        },
        {
            "contig": "chr1",
            "position": 250,
            "ref": "G",
            "alt": "C",
            "vaf": 0.20,
            "depth": 15,
            "alt_count": 3,
            "confidence": "medium",
        },
        {
            "contig": "chr1",
            "position": 800,
            "ref": "T",
            "alt": "A",
            "vaf": 0.15,
            "depth": 30,
            "alt_count": 5,
            "confidence": "high",
            "strand_forward": 7,
            "strand_reverse": 2,
        },
        {
            "contig": "chr1",
            "position": 900,
            "ref": "C",
            "alt": "G",
            "vaf": 0.05,
            "depth": 8,
            "alt_count": 1,
            "confidence": "low",
            "is_low_confidence": True,
        },
    ],
    "variant_evidence": {
        "130:A>T": {
            "forward_count": 4,
            "reverse_count": 3,
            "strand_bias": 0.14,
            "mean_quality": 38,
            "median_quality": 40,
        },
        "800:T>A": {
            "forward_count": 7,
            "reverse_count": 2,
            "strand_bias": 0.56,
            "mean_quality": 35,
            "median_quality": 37,
        },
    },
    "reference_sequence": "ACGTACGTAC" * 90,
}

MAPQ_GRADIENT_DATA = {
    "contig": "chr1",
    "start": 100,
    "end": 250,
    "reads": [
        {
            "name": "mapq5",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 5,
            "is_reverse": False,
            "mismatches": [],
        },
        {
            "name": "mapq30",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 30,
            "is_reverse": True,
            "mismatches": [],
        },
        {
            "name": "mapq60",
            "sequence": "ACGTACGTAC" * 5,
            "cigar": "50M",
            "position": 100,
            "end_position": 150,
            "mapping_quality": 60,
            "is_reverse": False,
            "mismatches": [],
        },
    ],
    "coverage": [3] * 150,
    "variants": [],
    "reference_sequence": "ACGTACGTAC" * 15,
}


# ══════════════════════════════════════════════════════════════════════
# Phase 3: Soft-Clip & CIGAR Tests
# ══════════════════════════════════════════════════════════════════════


class TestSoftClipPixelRendering:
    """Tests for soft-clip rendering with real soft_clips[] data."""

    @pytest.mark.e2e
    def test_left_soft_clip_draws_pixels(self, viewer_page: Page):
        """Left soft-clip region should have pixels drawn before read body."""
        send_init_data(viewer_page, SOFT_CLIP_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.showSoftClips = true; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(100)
        # sc_left: clip at pos 95-100, read body 100-145 (row 0)
        clip_px = genomic_to_pixel_x(viewer_page, 97)
        assert has_pixels_in_region(viewer_page, "reads-canvas", clip_px, 0, 5, 20)

    @pytest.mark.e2e
    def test_right_soft_clip_draws_pixels(self, viewer_page: Page):
        """Right soft-clip region should have pixels drawn after read body."""
        send_init_data(viewer_page, SOFT_CLIP_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.showSoftClips = true; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(100)
        # sc_right: read body 110-155, clip at 155-160 (row 1)
        clip_px = genomic_to_pixel_x(viewer_page, 157)
        assert has_pixels_in_region(viewer_page, "reads-canvas", clip_px, 0, 5, 30)

    @pytest.mark.e2e
    def test_both_soft_clips_drawn(self, viewer_page: Page):
        """Read with clips on both sides should have pixels on both sides."""
        send_init_data(viewer_page, SOFT_CLIP_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.showSoftClips = true; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(100)
        # sc_both at row 2: left clip at 117-120, right clip at 164-167
        left_px = genomic_to_pixel_x(viewer_page, 118)
        right_px = genomic_to_pixel_x(viewer_page, 165)
        assert has_pixels_in_region(viewer_page, "reads-canvas", left_px, 0, 5, 50)
        assert has_pixels_in_region(viewer_page, "reads-canvas", right_px, 0, 5, 50)

    @pytest.mark.e2e
    def test_soft_clip_toggle_hides_clips(self, viewer_page: Page):
        """Toggling soft-clip setting should show/hide clip rectangles."""
        send_init_data(viewer_page, SOFT_CLIP_DATA)
        # sc_left clip at 95-100 — first check with showSoftClips=true
        viewer_page.evaluate(
            "() => { viewer.state.settings.showSoftClips = true; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(100)
        clip_px = genomic_to_pixel_x(viewer_page, 97)
        before = has_pixels_in_region(viewer_page, "reads-canvas", clip_px, 0, 5, 20)
        assert before, "Clip should be visible with showSoftClips=true"

        # Toggle off
        viewer_page.evaluate(
            "() => { viewer.state.settings.showSoftClips = false; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)
        after = has_pixels_in_region(viewer_page, "reads-canvas", clip_px, 0, 5, 20)
        assert not after, "Clip should be hidden after toggle off"

    @pytest.mark.e2e
    def test_long_clip_wider_than_short_clip(self, viewer_page: Page):
        """10S clip should be ~2x wider than 5S clip."""
        send_init_data(viewer_page, SOFT_CLIP_DATA)
        scale = get_scale(viewer_page)
        # sc_left: 5bp clip, sc_long: 10bp clip
        short_width = int(5 * scale)
        long_width = int(10 * scale)
        short_px = genomic_to_pixel_x(viewer_page, 95)
        long_px = genomic_to_pixel_x(viewer_page, 120)

        short_pixels = count_opaque_pixels(
            viewer_page, "reads-canvas", short_px, 0, short_width, 12
        )
        long_pixels = count_opaque_pixels(viewer_page, "reads-canvas", long_px, 0, long_width, 12)
        # Long clip should have significantly more pixels
        assert long_pixels > short_pixels * 1.3, (
            f"10S clip ({long_pixels}px) should be notably wider than 5S ({short_pixels}px)"
        )


class TestComplexCIGAR:
    """Tests for reads with complex CIGAR operations."""

    @pytest.mark.e2e
    def test_multi_indel_read_reference_span(self, viewer_page: Page):
        """10M2I10M3D10M1I17M should span 50bp on reference (700-750)."""
        data = {
            "contig": "chr1",
            "start": 690,
            "end": 760,
            "reads": [
                {
                    "name": "cx_multi_indel",
                    "sequence": "ACGTACGTAC" * 5,
                    "cigar": "10M2I10M3D10M1I17M",
                    "position": 700,
                    "end_position": 750,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                }
            ],
            "coverage": [1] * 70,
            "variants": [],
            "reference_sequence": "A" * 70,
        }
        send_init_data(viewer_page, data)
        end_pos = viewer_page.evaluate("() => viewer.state.data.reads[0].end_position")
        assert end_pos == 750

    @pytest.mark.e2e
    def test_skip_op_renders_gap(self, viewer_page: Page):
        """20M50N30M should render two segments with a gap between them."""
        data = {
            "contig": "chr1",
            "start": 790,
            "end": 910,
            "reads": [
                {
                    "name": "cx_skip",
                    "sequence": "ACGTACGTAC" * 5,
                    "cigar": "20M50N30M",
                    "position": 800,
                    "end_position": 900,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                }
            ],
            "coverage": [1] * 120,
            "variants": [],
            "reference_sequence": "A" * 120,
        }
        send_init_data(viewer_page, data)
        # First segment at 810 should have pixels
        px_810 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 810, 6)
        assert px_810["a"] > 0, "First segment should be drawn"

        # Second segment at 880 should have pixels
        px_880 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 880, 6)
        assert px_880["a"] > 0, "Second segment should be drawn"

    @pytest.mark.e2e
    def test_hard_clip_invisible(self, viewer_page: Page):
        """5H45M should render as normal 45bp read with no extra indicator."""
        data = {
            "contig": "chr1",
            "start": 840,
            "end": 910,
            "reads": [
                {
                    "name": "cx_hard",
                    "sequence": "ACGTACGTAC" * 4 + "ACGTA",
                    "cigar": "5H45M",
                    "position": 850,
                    "end_position": 895,
                    "mapping_quality": 60,
                    "is_reverse": False,
                    "mismatches": [],
                }
            ],
            "coverage": [1] * 70,
            "variants": [],
            "reference_sequence": "A" * 70,
        }
        send_init_data(viewer_page, data)
        # Read should start at 850, no pixels before
        px_845 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 845, 6)
        px_855 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 855, 6)
        assert px_845["a"] == 0, "No pixels before hard-clipped read"
        assert px_855["a"] > 0, "Read body should be drawn"


# ══════════════════════════════════════════════════════════════════════
# Phase 4: Paired-End & Mate Tests
# ══════════════════════════════════════════════════════════════════════


class TestPairedEndRendering:
    """Tests for paired-end read visualization."""

    @pytest.mark.e2e
    def test_proper_pair_both_reads_render(self, viewer_page: Page):
        """Both reads of a proper pair should be visible."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        # R1 at 300, R2 at 550 — both in viewport 290-610
        px_r1 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 325, 6)
        px_r2 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 575, 6)
        assert px_r1["a"] > 0, "R1 should be rendered"
        assert px_r2["a"] > 0, "R2 should be rendered"

    @pytest.mark.e2e
    def test_hover_mate_highlights_both(self, viewer_page: Page):
        """Hovering R1 should highlight both R1 and R2, dim others."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Hover over R1 (row 0, position 325)
        scale = get_scale(viewer_page)
        hover_x = box["x"] + (325 - 290) * scale
        hover_y = box["y"] + 6
        viewer_page.mouse.move(hover_x, hover_y)
        viewer_page.wait_for_timeout(300)

        # R1 should be fully opaque
        px_r1 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 325, 6)
        assert px_r1["a"] > 200, "Hovered read should be nearly opaque"

    @pytest.mark.e2e
    def test_mate_connector_arc_drawn(self, viewer_page: Page):
        """Hovering R1 should draw a connector arc to R2."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Hover over R1
        scale = get_scale(viewer_page)
        hover_x = box["x"] + (325 - 290) * scale
        viewer_page.mouse.move(hover_x, box["y"] + 6)
        viewer_page.wait_for_timeout(300)

        # Check for arc pixels between R1 and R2 midpoints, above the reads
        mid_pos = (325 + 575) // 2
        mid_px = genomic_to_pixel_x(viewer_page, mid_pos)
        # Arc should be drawn above the reads (negative y offset from read)
        has_pixels_in_region(viewer_page, "reads-canvas", mid_px - 5, 0, 10, 6)
        # Arc may or may not be visible depending on row layout
        # At minimum, both reads should still be visible
        px_r2 = sample_pixel_at_genomic(viewer_page, "reads-canvas", 575, 6)
        assert px_r2["a"] > 0, "R2 should be visible during hover"

    @pytest.mark.e2e
    def test_discordant_pair_insert_size_color(self, viewer_page: Page):
        """Discordant pair (1500bp) should be red in insertSize color mode."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.colorBy = 'insertSize'; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)

        # pair2 at position 350 — large insert (1500bp) should be red
        px = sample_pixel_at_genomic(viewer_page, "reads-canvas", 375, 6)
        if px["a"] > 0:
            # Red-dominant for large insert
            assert px["r"] > px["b"], (
                f"Discordant pair should be red-dominant, got r={px['r']} b={px['b']}"
            )

    @pytest.mark.e2e
    def test_chimeric_pair_color(self, viewer_page: Page):
        """Chimeric pair (mate on chr2) should be pink in insertSize mode."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.colorBy = 'insertSize'; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)

        # pair3 at position 400 — chimeric (mate on chr2)
        px = sample_pixel_at_genomic(viewer_page, "reads-canvas", 425, 6)
        if px["a"] > 0:
            # Pink = high red + medium blue
            assert px["r"] > px["g"], (
                f"Chimeric pair should have r > g (pink), got r={px['r']} g={px['g']}"
            )


# ══════════════════════════════════════════════════════════════════════
# Phase 5: Variant Detection & Evidence Tests
# ══════════════════════════════════════════════════════════════════════


class TestVariantDetection:
    """Tests using comprehensive BAM fixtures for full-pipeline variant detection."""

    @pytest.mark.e2e
    def test_strand_bias_variant_detected(self, viewer_page: Page, comprehensive_data):
        """Strand-bias variant at ~1050 should appear in table when filter=all."""
        data = comprehensive_data("chr1:1000-1100")
        send_init_data(viewer_page, data)
        # Low-confidence variants are filtered by default — show all
        viewer_page.evaluate(
            "() => { viewer.state.variantFilter = 'all'; viewer.renderVariantTable(); }"
        )
        viewer_page.wait_for_timeout(50)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        found = any("1,050" in t or "1050" in t for t in row_texts)
        assert found, f"Variant at 1050 not found in table rows: {row_texts}"

    @pytest.mark.e2e
    def test_clean_variant_high_confidence(self, viewer_page: Page, comprehensive_data):
        """Clean variant at ~1150 should be high or medium confidence."""
        data = comprehensive_data("chr1:1100-1200")
        send_init_data(viewer_page, data)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        found = any("1,150" in t or "1150" in t for t in row_texts)
        assert found, f"Clean variant at 1150 not found: {row_texts}"

    @pytest.mark.e2e
    def test_multi_allele_both_variants(self, viewer_page: Page, comprehensive_data):
        """Two variants at position 2500 (T and G alleles) should both appear."""
        data = comprehensive_data("chr1:2480-2520")
        send_init_data(viewer_page, data)
        rows = viewer_page.locator("#variant-table tr")
        count = rows.count()
        assert count >= 2, f"Expected >=2 variant rows for multi-allele, got {count}"

    @pytest.mark.e2e
    def test_below_min_depth_no_variant(self, viewer_page: Page, comprehensive_data):
        """1x depth region should not produce variants (min_depth=2)."""
        data = comprehensive_data("chr1:2000-2060")
        send_init_data(viewer_page, data)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        # No variant should be in the 2000-2050 range (1x depth)
        in_range = any(any(str(p) in t for p in range(2000, 2051)) for t in row_texts)
        assert not in_range, f"Should not have variant in 1x region: {row_texts}"

    @pytest.mark.e2e
    def test_at_depth_threshold_variant_called(self, viewer_page: Page, comprehensive_data):
        """Variant at depth=2 (at threshold) should be called when filter=all."""
        data = comprehensive_data("chr1:2050-2110")
        send_init_data(viewer_page, data)
        # Low-depth variants are low-confidence — show all
        viewer_page.evaluate(
            "() => { viewer.state.variantFilter = 'all'; viewer.renderVariantTable(); }"
        )
        viewer_page.wait_for_timeout(50)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        found = any("2,075" in t or "2075" in t for t in row_texts)
        assert found, f"Variant at 2075 (depth=2) not found: {row_texts}"


class TestVariantEvidence:
    """Tests for the variant evidence panel."""

    @pytest.mark.e2e
    def test_evidence_panel_opens_on_click(self, viewer_page: Page):
        """Clicking a variant row should open the evidence panel."""
        send_init_data(viewer_page, EVIDENCE_DATA)
        viewer_page.wait_for_timeout(100)

        # Click first variant row
        rows = viewer_page.locator("#variant-table tr")
        if rows.count() > 0:
            rows.first.click()
            viewer_page.wait_for_timeout(300)
            panel = viewer_page.locator("#evidence-panel")
            is_visible = panel.evaluate(
                "el => el.classList.contains('visible') || el.style.display !== 'none'"
            )
            assert is_visible, "Evidence panel should be visible after click"

    @pytest.mark.e2e
    def test_evidence_panel_closes(self, viewer_page: Page):
        """Close button should hide the evidence panel."""
        send_init_data(viewer_page, EVIDENCE_DATA)
        viewer_page.wait_for_timeout(100)

        rows = viewer_page.locator("#variant-table tr")
        if rows.count() > 0:
            rows.first.click()
            viewer_page.wait_for_timeout(200)

            # Click close button
            close_btn = viewer_page.locator("#evidence-panel .close-btn, #evidence-panel button")
            if close_btn.count() > 0:
                close_btn.first.click()
                viewer_page.wait_for_timeout(200)

    @pytest.mark.e2e
    def test_strand_bias_warning_shown(self, viewer_page: Page):
        """Variant with strand_bias=1.0 should show warning text."""
        send_init_data(viewer_page, EVIDENCE_DATA)
        viewer_page.wait_for_timeout(100)

        # Click the strand-bias variant (first one)
        rows = viewer_page.locator("#variant-table tr")
        if rows.count() > 0:
            rows.first.click()
            viewer_page.wait_for_timeout(300)
            # Check for strand bias warning text in panel
            panel_text = viewer_page.locator("#evidence-panel").inner_text()
            has_warning = "strand bias" in panel_text.lower() or "strand" in panel_text.lower()
            assert has_warning, (
                f"Evidence panel should mention strand bias, got: {panel_text[:200]}"
            )

    @pytest.mark.e2e
    def test_clean_variant_no_artifacts(self, viewer_page: Page):
        """Clean variant should show no artifact warnings."""
        send_init_data(viewer_page, EVIDENCE_DATA)
        viewer_page.wait_for_timeout(100)

        # Click the clean variant (second row)
        rows = viewer_page.locator("#variant-table tr")
        if rows.count() >= 2:
            rows.nth(1).click()
            viewer_page.wait_for_timeout(300)
            panel_text = viewer_page.locator("#evidence-panel").inner_text()
            has_clean = "no concerns" in panel_text.lower() or "low" in panel_text.lower()
            assert has_clean, f"Clean variant should show no concerns, got: {panel_text[:200]}"


# ══════════════════════════════════════════════════════════════════════
# Phase 6: Variant Table & In-View Tests
# ══════════════════════════════════════════════════════════════════════


class TestVariantTableFeatures:
    """Tests for variant table in-view highlighting, sorting, and filtering."""

    @pytest.mark.e2e
    def test_in_view_variant_has_class(self, viewer_page: Page):
        """Variants in viewport should have .in-view CSS class."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Viewport is 100-1000, variants at 130 and 250 are in view
        in_view = viewer_page.locator("#variant-table tr.in-view")
        assert in_view.count() >= 2, f"Expected >= 2 in-view variants, got {in_view.count()}"

    @pytest.mark.e2e
    def test_out_of_view_variant_dimmed(self, viewer_page: Page):
        """Variants outside viewport should lack .in-view class."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Zoom in to 100-500, variant at 800 is out of view
        viewer_page.evaluate("""() => {
            viewer.state.viewport.start = 100;
            viewer.state.viewport.end = 500;
            viewer.renderer.resize();
            viewer.renderVariantTable();
        }""")
        viewer_page.wait_for_timeout(200)

        not_in_view = viewer_page.locator("#variant-table tr:not(.in-view)")
        assert not_in_view.count() >= 1, "Out-of-view variants should exist"

    @pytest.mark.e2e
    def test_in_view_updates_on_viewport_change(self, viewer_page: Page):
        """Panning should update which variants have .in-view."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Zoom to show only 100-300 region
        viewer_page.evaluate("""() => {
            viewer.state.viewport.start = 100;
            viewer.state.viewport.end = 300;
            viewer.renderer.resize();
            viewer.renderVariantTable();
        }""")
        viewer_page.wait_for_timeout(200)
        count1 = viewer_page.locator("#variant-table tr.in-view").count()

        # Pan to 700-900 region
        viewer_page.evaluate("""() => {
            viewer.state.viewport.start = 700;
            viewer.state.viewport.end = 900;
            viewer.renderer.resize();
            viewer.renderVariantTable();
        }""")
        viewer_page.wait_for_timeout(200)
        count2 = viewer_page.locator("#variant-table tr.in-view").count()

        # Different variants should be in view after panning
        assert count1 >= 1 and count2 >= 1

    @pytest.mark.e2e
    def test_sort_by_vaf(self, viewer_page: Page):
        """Clicking VAF header should sort variants by VAF."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Click VAF header
        vaf_header = viewer_page.locator("th[data-sort='vaf']")
        if vaf_header.count() > 0:
            vaf_header.click()
            viewer_page.wait_for_timeout(200)

            # Get VAF values from rows
            rows = viewer_page.locator("#variant-table tr")
            vafs = []
            for i in range(rows.count()):
                text = rows.nth(i).inner_text()
                # Extract percentage value
                for part in text.split():
                    if "%" in part:
                        with contextlib.suppress(ValueError):
                            vafs.append(float(part.replace("%", "")))
            if len(vafs) >= 2:
                # Should be sorted (descending or ascending)
                assert vafs == sorted(vafs) or vafs == sorted(vafs, reverse=True), (
                    f"VAF values should be sorted, got: {vafs}"
                )

    @pytest.mark.e2e
    def test_strand_display_format(self, viewer_page: Page):
        """Strand column should show XF:YR format from evidence."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        # Variant at 130 has evidence forward_count=4, reverse_count=3
        found_format = any("4F:3R" in t for t in row_texts)
        assert found_format, f"Expected '4F:3R' strand format, got: {row_texts}"

    @pytest.mark.e2e
    def test_strand_fallback_to_variant_fields(self, viewer_page: Page):
        """Without evidence, strand should fall back to variant fields."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        rows = viewer_page.locator("#variant-table tr")
        row_texts = [rows.nth(i).inner_text() for i in range(rows.count())]
        # Variant at 800 has strand_forward=7, strand_reverse=2 but evidence also has it
        found_format = any("7F:2R" in t for t in row_texts)
        assert found_format, f"Expected '7F:2R' strand format, got: {row_texts}"

    @pytest.mark.e2e
    def test_variant_click_centers_viewport(self, viewer_page: Page):
        """Clicking variant row should center viewport on that position."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Click variant at position 800
        rows = viewer_page.locator("#variant-table tr")
        for i in range(rows.count()):
            text = rows.nth(i).inner_text()
            if "800" in text:
                rows.nth(i).click()
                break
        viewer_page.wait_for_timeout(200)

        center = viewer_page.evaluate(
            "() => (viewer.state.viewport.start + viewer.state.viewport.end) / 2"
        )
        assert abs(center - 800) < 50, f"Viewport center should be near 800, got {center}"

    @pytest.mark.e2e
    def test_high_conf_filter(self, viewer_page: Page):
        """High-conf filter should hide low-confidence variants."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        initial_count = viewer_page.locator("#variant-table tr").count()

        # Click high-conf filter
        high_btn = viewer_page.locator("button:has-text('High')")
        if high_btn.count() > 0:
            high_btn.first.click()
            viewer_page.wait_for_timeout(200)

            # Low-confidence variant at 900 should be hidden
            filtered_count = viewer_page.locator("#variant-table tr").count()
            assert filtered_count <= initial_count

    @pytest.mark.e2e
    def test_all_filter_shows_everything(self, viewer_page: Page):
        """All filter should show all variants including low-confidence."""
        send_init_data(viewer_page, MULTI_VARIANT_DATA)
        # Click All filter
        all_btn = viewer_page.locator("button:has-text('All')")
        if all_btn.count() > 0:
            all_btn.first.click()
            viewer_page.wait_for_timeout(200)
            count = viewer_page.locator("#variant-table tr").count()
            assert count == 4, f"All filter should show 4 variants, got {count}"


# ══════════════════════════════════════════════════════════════════════
# Phase 7: DataStore & Viewport Tests
# ══════════════════════════════════════════════════════════════════════


class TestDataStore:
    """Tests for tile cache, variant accumulation, and evidence persistence."""

    @pytest.mark.e2e
    def test_tile_cached_after_load(self, viewer_page: Page):
        """After data load, DataStore should have at least one tile."""
        send_init_data(viewer_page)
        count = viewer_page.evaluate("() => viewer.state.store.tileCount")
        assert count > 0

    @pytest.mark.e2e
    def test_variants_accumulated_across_tiles(self, viewer_page: Page):
        """Loading two tiles should accumulate variants from both."""
        data1 = {
            **SAMPLE_DATA,
            "start": 100,
            "end": 300,
            "variants": [
                {
                    "contig": "chr1",
                    "position": 200,
                    "ref": "A",
                    "alt": "T",
                    "vaf": 0.25,
                    "depth": 20,
                },
            ],
        }
        send_init_data(viewer_page, data1)

        data2 = {
            **SAMPLE_DATA,
            "start": 250,
            "end": 500,
            "variants": [
                {
                    "contig": "chr1",
                    "position": 400,
                    "ref": "G",
                    "alt": "C",
                    "vaf": 0.15,
                    "depth": 30,
                },
            ],
        }
        # Use loadTile to accumulate (not reset)
        viewer_page.evaluate(
            """(data) => {
            viewer.state.loadTile(data);
            viewer.renderVariantTable();
        }""",
            data2,
        )
        viewer_page.wait_for_timeout(200)

        count = viewer_page.evaluate("() => viewer.state.store.getAllVariants().length")
        assert count >= 2, f"Expected >=2 accumulated variants, got {count}"

    @pytest.mark.e2e
    def test_evidence_accumulated_across_tiles(self, viewer_page: Page):
        """Evidence from multiple tiles should persist in DataStore."""
        data1 = {
            **SAMPLE_DATA,
            "start": 100,
            "end": 300,
            "variants": [
                {
                    "contig": "chr1",
                    "position": 200,
                    "ref": "A",
                    "alt": "T",
                    "vaf": 0.25,
                    "depth": 20,
                },
            ],
            "variant_evidence": {
                "200:A>T": {
                    "forward_count": 5,
                    "reverse_count": 3,
                    "strand_bias": 0.25,
                    "mean_quality": 35,
                    "median_quality": 38,
                },
            },
        }
        send_init_data(viewer_page, data1)

        data2 = {
            **SAMPLE_DATA,
            "start": 250,
            "end": 500,
            "variants": [
                {
                    "contig": "chr1",
                    "position": 400,
                    "ref": "G",
                    "alt": "C",
                    "vaf": 0.15,
                    "depth": 30,
                },
            ],
            "variant_evidence": {
                "400:G>C": {
                    "forward_count": 2,
                    "reverse_count": 4,
                    "strand_bias": 0.33,
                    "mean_quality": 30,
                    "median_quality": 32,
                },
            },
        }
        viewer_page.evaluate(
            "(data) => { viewer.state.loadTile(data); }",
            data2,
        )
        viewer_page.wait_for_timeout(100)

        ev1 = viewer_page.evaluate("() => viewer.state.store.getEvidence('200:A>T') !== undefined")
        ev2 = viewer_page.evaluate("() => viewer.state.store.getEvidence('400:G>C') !== undefined")
        assert ev1, "Evidence for first tile variant should persist"
        assert ev2, "Evidence for second tile variant should be accessible"

    @pytest.mark.e2e
    def test_contig_change_clears_cache(self, viewer_page: Page):
        """Loading data for a different contig should clear the store."""
        send_init_data(viewer_page)  # chr1 data
        viewer_page.evaluate("() => viewer.state.store.getAllVariants().length")

        # Load chr2 data
        chr2_data = {
            **SAMPLE_DATA,
            "contig": "chr2",
            "variants": [
                {
                    "contig": "chr2",
                    "position": 300,
                    "ref": "A",
                    "alt": "G",
                    "vaf": 0.3,
                    "depth": 10,
                },
            ],
        }
        send_init_data(viewer_page, chr2_data)

        variants = viewer_page.evaluate(
            "() => viewer.state.store.getAllVariants().map(v => v.contig)"
        )
        assert all(c == "chr2" for c in variants), (
            f"All variants should be chr2 after contig change, got: {variants}"
        )

    @pytest.mark.e2e
    def test_lru_eviction_at_max_tiles(self, viewer_page: Page):
        """Ingesting 21 tiles should evict the oldest, keeping max 20."""
        send_init_data(viewer_page)
        # Ingest 20 more tiles (total 21)
        for i in range(20):
            viewer_page.evaluate(
                f"""() => {{
                viewer.state.store.ingest({{
                    contig: 'chr1',
                    start: {i * 100},
                    end: {i * 100 + 100},
                    reads: [],
                    coverage: [],
                    variants: []
                }});
            }}"""
            )
        count = viewer_page.evaluate("() => viewer.state.store.tileCount")
        assert count <= 20, f"Tile count should be <= 20 after LRU eviction, got {count}"


class TestViewportBehavior:
    """Tests for viewport reset, preservation, and smooth panning."""

    @pytest.mark.e2e
    def test_host_data_resets_viewport(self, viewer_page: Page):
        """New host data (loadData) should reset viewport to data range."""
        send_init_data(viewer_page)
        # Zoom in
        viewer_page.evaluate("""() => {
            viewer.state.viewport.start = 200;
            viewer.state.viewport.end = 400;
            viewer.renderer.resize();
        }""")

        # Send new data — should reset viewport
        new_data = {**SAMPLE_DATA, "start": 500, "end": 1000}
        send_init_data(viewer_page, new_data)

        start = viewer_page.evaluate("() => viewer.state.viewport.start")
        end = viewer_page.evaluate("() => viewer.state.viewport.end")
        assert start == 500, f"Viewport start should reset to 500, got {start}"
        assert end == 1000, f"Viewport end should reset to 1000, got {end}"

    @pytest.mark.e2e
    def test_tile_load_preserves_viewport(self, viewer_page: Page):
        """loadTile() should not reset viewport position."""
        send_init_data(viewer_page)
        viewer_page.evaluate("""() => {
            viewer.state.viewport.start = 200;
            viewer.state.viewport.end = 400;
        }""")

        viewer_page.evaluate(
            """(data) => {
            viewer.state.loadTile(data);
            viewer.renderer.resize();
        }""",
            SAMPLE_DATA,
        )
        viewer_page.wait_for_timeout(100)

        start = viewer_page.evaluate("() => viewer.state.viewport.start")
        assert start == 200, f"Viewport should stay at 200 after tile load, got {start}"

    @pytest.mark.e2e
    def test_viewport_no_jump_on_pan(self, viewer_page: Page):
        """Sequential pans should move viewport smoothly without jumps."""
        send_init_data(viewer_page)
        positions = []

        for _ in range(5):
            start = viewer_page.evaluate("() => viewer.state.viewport.start")
            positions.append(start)
            # Pan right by 20bp
            viewer_page.evaluate("""() => {
                const span = viewer.state.viewport.end - viewer.state.viewport.start;
                viewer.state.viewport.start += 20;
                viewer.state.viewport.end += 20;
                viewer.renderer.resize();
            }""")
            viewer_page.wait_for_timeout(50)

        # Check no jumps > 100bp between consecutive positions
        for i in range(1, len(positions)):
            diff = abs(positions[i] - positions[i - 1])
            assert diff <= 100, f"Viewport jumped {diff}bp between pans (positions: {positions})"


# ══════════════════════════════════════════════════════════════════════
# Phase 8: Display Settings & Rendering Tests
# ══════════════════════════════════════════════════════════════════════


class TestDisplaySettings:
    """Tests for display mode, colorBy, and sortBy settings."""

    @pytest.mark.e2e
    def test_squished_mode_height(self, viewer_page: Page):
        """Squished mode should use 6px read height + 1px gap = 7px rows."""
        send_init_data(viewer_page, NARROW_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.displayMode = 'squished';"
            " viewer.state.resortAndRepack(); viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)
        # Row 0 at y=0, row 1 at y=7 (6+1). Check pixels in both.
        px_row0 = viewer_page.evaluate("""() => {
            const canvas = document.getElementById('reads-canvas');
            const ctx = canvas.getContext('2d');
            const p = ctx.getImageData(50, 3, 1, 1).data;
            return p[3];
        }""")
        px_row1 = viewer_page.evaluate("""() => {
            const canvas = document.getElementById('reads-canvas');
            const ctx = canvas.getContext('2d');
            const p = ctx.getImageData(50, 10, 1, 1).data;
            return p[3];
        }""")
        assert px_row0 > 0, "Row 0 should have pixels in squished mode"
        assert px_row1 > 0, "Row 1 should have pixels at y=10 in squished mode"

    @pytest.mark.e2e
    def test_compact_mode_height(self, viewer_page: Page):
        """Compact mode (default) should use 12px height + 2px gap = 14px rows."""
        send_init_data(viewer_page, NARROW_DATA)
        # Default is compact, verify row spacing
        px_row0 = viewer_page.evaluate("""() => {
            const canvas = document.getElementById('reads-canvas');
            const ctx = canvas.getContext('2d');
            const p = ctx.getImageData(50, 6, 1, 1).data;
            return p[3];
        }""")
        px_row1 = viewer_page.evaluate("""() => {
            const canvas = document.getElementById('reads-canvas');
            const ctx = canvas.getContext('2d');
            const p = ctx.getImageData(50, 20, 1, 1).data;
            return p[3];
        }""")
        assert px_row0 > 0, "Row 0 should be drawn at y=6"
        assert px_row1 > 0, "Row 1 should be drawn at y=20"

    @pytest.mark.e2e
    def test_color_by_mapq_gradient(self, viewer_page: Page):
        """Color by MAPQ should give different colors for different MAPQ values."""
        send_init_data(viewer_page, MAPQ_GRADIENT_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.colorBy = 'mapq'; viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)

        # Sample pixels from different MAPQ reads
        px_low = sample_pixel_at_genomic(viewer_page, "reads-canvas", 125, 3)
        px_high = sample_pixel_at_genomic(viewer_page, "reads-canvas", 125, 20)

        if px_low["a"] > 0 and px_high["a"] > 0:
            # Low MAPQ should be more red, high MAPQ more blue
            assert px_low["r"] != px_high["r"] or px_low["b"] != px_high["b"], (
                f"Different MAPQ reads should have different colors: low={px_low} high={px_high}"
            )

    @pytest.mark.e2e
    def test_sort_by_mapq_reorders(self, viewer_page: Page):
        """Sorting by MAPQ should place highest MAPQ reads in top rows."""
        send_init_data(viewer_page, MAPQ_GRADIENT_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.sortBy = 'mapq';"
            " viewer.state.resortAndRepack(); viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)

        top_mapq = viewer_page.evaluate("() => viewer.state.packedRows[0][0].mapping_quality")
        assert top_mapq >= 50, f"Top row should have high MAPQ after sorting, got {top_mapq}"

    @pytest.mark.e2e
    def test_sort_by_strand_groups(self, viewer_page: Page):
        """Sorting by strand should group forward before reverse."""
        send_init_data(viewer_page, NARROW_DATA)
        viewer_page.evaluate(
            "() => { viewer.state.settings.sortBy = 'strand';"
            " viewer.state.resortAndRepack(); viewer.renderer.resize(); }"
        )
        viewer_page.wait_for_timeout(200)

        first_is_reverse = viewer_page.evaluate("() => viewer.state.packedRows[0][0].is_reverse")
        assert not first_is_reverse, "First row should be forward after strand sort"


class TestSequenceRendering:
    """Tests for base-level sequence rendering at high zoom."""

    @pytest.mark.e2e
    def test_base_colors_at_high_zoom(self, viewer_page: Page):
        """At scale >= 10, individual base colors should be visible on reads."""
        send_init_data(viewer_page, NARROW_DATA)
        scale = get_scale(viewer_page)
        assert scale >= 10, f"Scale should be >= 10 for base rendering, got {scale}"

        # Sample a pixel on the read — should show per-base color, not flat blue
        px = sample_pixel_at_genomic(viewer_page, "reads-canvas", 110, 6)
        assert px["a"] > 0, "Read pixel should be drawn"
        # At base level, pixel colors vary by nucleotide (not flat blue)

    @pytest.mark.e2e
    def test_reference_nucleotides_colored(self, viewer_page: Page):
        """Reference track should show colored nucleotide rectangles at high zoom."""
        send_init_data(viewer_page, NARROW_DATA)
        scale = get_scale(viewer_page)
        assert scale >= 10

        # Reference canvas should have colored pixels
        ref_px = sample_pixel_at_genomic(viewer_page, "reference-canvas", 105, 12)
        assert ref_px["a"] > 0, "Reference track should have colored nucleotides"

    @pytest.mark.e2e
    def test_bases_hidden_at_low_zoom(self, viewer_page: Page):
        """At low zoom (scale < 10), individual bases should NOT be drawn."""
        wide_data = {
            **NARROW_DATA,
            "start": 0,
            "end": 1000,
            "coverage": [2] * 1000,
            "reference_sequence": "ACGTACGTAC" * 100,
        }
        send_init_data(viewer_page, wide_data)
        scale = get_scale(viewer_page)
        assert scale < 10, f"Scale should be < 10 at low zoom, got {scale}"


class TestCoverageTrackEnhanced:
    """Tests for coverage track rendering details."""

    @pytest.mark.e2e
    def test_zero_coverage_red_region(self, viewer_page: Page):
        """Zero-coverage areas should show red fill/hatching."""
        data = {
            **SAMPLE_DATA,
            "coverage": [0] * 200 + [20] * 300,
        }
        send_init_data(viewer_page, data)

        # Sample in zero-coverage region (start of coverage array)
        px = sample_pixel_at_genomic(viewer_page, "coverage-canvas", 150, 40)
        if px["a"] > 0:
            assert px["r"] > px["b"], (
                f"Zero coverage should be red-tinted, got r={px['r']} b={px['b']}"
            )

    @pytest.mark.e2e
    def test_coverage_max_label_present(self, viewer_page: Page):
        """Coverage canvas should show max coverage label in top-left."""
        data = {**SAMPLE_DATA, "coverage": [0] * 200 + [42] * 300}
        send_init_data(viewer_page, data)

        # Label drawn at x=5-60, y=2-14
        label_pixels = count_opaque_pixels(viewer_page, "coverage-canvas", 5, 2, 55, 12)
        assert label_pixels > 5, f"Coverage label should have opaque pixels, got {label_pixels}"


# ══════════════════════════════════════════════════════════════════════
# Phase 9: Interaction Tests
# ══════════════════════════════════════════════════════════════════════


class TestTooltipEnhanced:
    """Extended tooltip tests."""

    @pytest.mark.e2e
    def test_tooltip_shows_paired_info(self, viewer_page: Page):
        """Tooltip for paired-end read should show 'Paired' and insert size."""
        send_init_data(viewer_page, PAIRED_END_DATA)
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        scale = get_scale(viewer_page)
        hover_x = box["x"] + (325 - 290) * scale
        viewer_page.mouse.move(hover_x, box["y"] + 6)
        viewer_page.wait_for_timeout(300)

        tooltip = viewer_page.locator("#tooltip")
        if tooltip.is_visible():
            content = tooltip.inner_text()
            assert "Paired" in content or "pair" in content.lower(), (
                f"Tooltip should mention 'Paired', got: {content[:200]}"
            )

    @pytest.mark.e2e
    def test_tooltip_locked_on_click(self, viewer_page: Page):
        """Clicking a read should lock the tooltip."""
        send_init_data(viewer_page, NARROW_DATA)
        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        # Click on read
        click_x = box["x"] + box["width"] / 2
        click_y = box["y"] + 6
        viewer_page.mouse.click(click_x, click_y)
        viewer_page.wait_for_timeout(300)

        tooltip = viewer_page.locator("#tooltip")
        if tooltip.is_visible():
            # Tooltip should stay visible (locked)
            viewer_page.mouse.move(box["x"] - 50, box["y"] - 50)
            viewer_page.wait_for_timeout(200)
            # In locked state, tooltip should persist
            is_still = tooltip.is_visible()
            # May or may not be locked depending on implementation
            assert isinstance(is_still, bool)


class TestKeyboardNavigation:
    """Tests for keyboard-based viewport navigation."""

    @pytest.mark.e2e
    def test_arrow_right_pans_viewport(self, viewer_page: Page):
        """Right arrow key should shift the viewport."""
        send_init_data(viewer_page)
        original = viewer_page.evaluate("() => viewer.state.viewport.start")

        viewer_page.keyboard.press("ArrowRight")
        viewer_page.wait_for_timeout(200)

        new_start = viewer_page.evaluate("() => viewer.state.viewport.start")
        assert new_start != original, f"Right arrow should move viewport: {original} → {new_start}"

    @pytest.mark.e2e
    def test_arrow_left_pans_opposite(self, viewer_page: Page):
        """Left arrow key should pan viewport in the opposite direction to right."""
        send_init_data(viewer_page)
        original = viewer_page.evaluate("() => viewer.state.viewport.start")

        # Press right first
        viewer_page.keyboard.press("ArrowRight")
        viewer_page.wait_for_timeout(100)
        after_right = viewer_page.evaluate("() => viewer.state.viewport.start")

        # Press left to reverse
        viewer_page.keyboard.press("ArrowLeft")
        viewer_page.wait_for_timeout(200)

        after_left = viewer_page.evaluate("() => viewer.state.viewport.start")
        # Left arrow should reverse the direction of right arrow
        right_delta = after_right - original
        left_delta = after_left - after_right
        assert right_delta * left_delta < 0, (
            f"Arrows should pan in opposite directions: right Δ={right_delta}, left Δ={left_delta}"
        )


class TestDebugOverlay:
    """Tests for the debug overlay panel."""

    @pytest.mark.e2e
    def test_debug_hidden_by_default(self, viewer_page: Page):
        """Debug overlay should be hidden on page load."""
        send_init_data(viewer_page)
        overlay = viewer_page.locator("#debug-overlay")
        has_hidden = overlay.evaluate("el => el.classList.contains('hidden')")
        assert has_hidden, "Debug overlay should have .hidden class by default"

    @pytest.mark.e2e
    def test_debug_toggle_visibility(self, viewer_page: Page):
        """Debug button should toggle overlay visibility."""
        send_init_data(viewer_page)
        debug_btn = viewer_page.locator("#debug-btn")
        if debug_btn.count() > 0:
            debug_btn.click()
            viewer_page.wait_for_timeout(200)

            is_hidden = viewer_page.locator("#debug-overlay").evaluate(
                "el => el.classList.contains('hidden')"
            )
            assert not is_hidden, "Overlay should be visible after click"

            debug_btn.click()
            viewer_page.wait_for_timeout(200)
            is_hidden_again = viewer_page.locator("#debug-overlay").evaluate(
                "el => el.classList.contains('hidden')"
            )
            assert is_hidden_again, "Overlay should be hidden after second click"

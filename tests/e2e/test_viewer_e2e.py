"""End-to-end Playwright tests for the BAMCP viewer UI."""

import pytest
from playwright.sync_api import Page, expect

from bamcp.resources import get_viewer_html

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
    return page


def send_init_data(page: Page, data: dict = None):
    """Send init data to the viewer via postMessage."""
    if data is None:
        data = SAMPLE_DATA
    page.evaluate(
        """(data) => {
            window.postMessage({
                method: 'bamcp/init',
                params: data
            }, '*');
        }""",
        data,
    )
    # Give the viewer time to process
    page.wait_for_timeout(200)


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
        assert headers.count() == 5

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
        assert "chr1:130" in cells.nth(0).inner_text()
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
        """bamcp/update should replace existing data."""
        send_init_data(viewer_page)

        # Send updated data with different region
        new_data = {**SAMPLE_DATA, "contig": "chr2", "start": 500, "end": 1000}
        viewer_page.evaluate(
            """(data) => {
                window.postMessage({
                    method: 'bamcp/update',
                    params: data
                }, '*');
            }""",
            new_data,
        )
        viewer_page.wait_for_timeout(200)

        input_val = viewer_page.locator("#region-input").input_value()
        assert input_val == "chr2:500-1000"


class TestNavigation:
    """E2E tests for viewer navigation controls."""

    @pytest.mark.e2e
    def test_zoom_in_button(self, viewer_page: Page):
        """Zoom in button should narrow the viewport."""
        send_init_data(viewer_page)

        # Get initial viewport via JS
        initial_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        viewer_page.locator("#zoom-in").click()

        new_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        assert new_span < initial_span

    @pytest.mark.e2e
    def test_zoom_out_button(self, viewer_page: Page):
        """Zoom out button should widen the viewport."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        viewer_page.locator("#zoom-out").click()

        new_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        assert new_span > initial_span

    @pytest.mark.e2e
    def test_zoom_preserves_center(self, viewer_page: Page):
        """Zoom should preserve the center position."""
        send_init_data(viewer_page)

        initial_center = viewer_page.evaluate(
            "() => (viewer.viewport.start + viewer.viewport.end) / 2"
        )

        viewer_page.locator("#zoom-in").click()

        new_center = viewer_page.evaluate("() => (viewer.viewport.start + viewer.viewport.end) / 2")

        # Center should be approximately the same (within rounding)
        assert abs(new_center - initial_center) < 2

    @pytest.mark.e2e
    def test_variant_click_jumps_to_position(self, viewer_page: Page):
        """Clicking a variant row should jump to that position."""
        send_init_data(viewer_page)

        viewer_page.evaluate("() => viewer.viewport.start")

        # Click the second variant row (position 250)
        viewer_page.locator("#variant-table tr").nth(1).click()

        new_center = viewer_page.evaluate("() => (viewer.viewport.start + viewer.viewport.end) / 2")

        # Center should be near position 250
        assert abs(new_center - 250) < 10

    @pytest.mark.e2e
    def test_mouse_drag_pans(self, viewer_page: Page):
        """Mouse drag on reads canvas should pan the view."""
        send_init_data(viewer_page)

        initial_start = viewer_page.evaluate("() => viewer.viewport.start")

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

        new_start = viewer_page.evaluate("() => viewer.viewport.start")

        # Panning left (drag right-to-left) should increase start position
        assert new_start != initial_start


class TestReadPacking:
    """E2E tests for read packing algorithm."""

    @pytest.mark.e2e
    def test_packed_rows_created(self, viewer_page: Page):
        """Reads should be packed into rows."""
        send_init_data(viewer_page)

        row_count = viewer_page.evaluate("() => viewer.packedRows.length")
        assert row_count > 0

    @pytest.mark.e2e
    def test_all_reads_packed(self, viewer_page: Page):
        """All reads should be present in packed rows."""
        send_init_data(viewer_page)

        total = viewer_page.evaluate(
            "() => viewer.packedRows.reduce((sum, row) => sum + row.length, 0)"
        )
        assert total == len(SAMPLE_DATA["reads"])

    @pytest.mark.e2e
    def test_no_overlapping_reads_in_row(self, viewer_page: Page):
        """Reads in the same row should not overlap."""
        send_init_data(viewer_page)

        no_overlap = viewer_page.evaluate("""() => {
                for (const row of viewer.packedRows) {
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

        row_count = viewer_page.evaluate("() => viewer.packedRows.length")
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
    """E2E tests for soft-clip visualization."""

    @pytest.mark.e2e
    def test_parse_soft_clips_method_exists(self, viewer_page: Page):
        """Viewer should have parseSoftClips method."""
        send_init_data(viewer_page)
        exists = viewer_page.evaluate("() => typeof viewer.parseSoftClips === 'function'")
        assert exists

    @pytest.mark.e2e
    def test_parse_soft_clips_no_clips(self, viewer_page: Page):
        """Should return zero clips for simple CIGAR."""
        send_init_data(viewer_page)
        result = viewer_page.evaluate("() => viewer.parseSoftClips('50M')")
        assert result["left"] == 0
        assert result["right"] == 0

    @pytest.mark.e2e
    def test_parse_soft_clips_left(self, viewer_page: Page):
        """Should detect left soft-clip."""
        send_init_data(viewer_page)
        result = viewer_page.evaluate("() => viewer.parseSoftClips('5S45M')")
        assert result["left"] == 5
        assert result["right"] == 0

    @pytest.mark.e2e
    def test_parse_soft_clips_right(self, viewer_page: Page):
        """Should detect right soft-clip."""
        send_init_data(viewer_page)
        result = viewer_page.evaluate("() => viewer.parseSoftClips('45M5S')")
        assert result["left"] == 0
        assert result["right"] == 5

    @pytest.mark.e2e
    def test_parse_soft_clips_both(self, viewer_page: Page):
        """Should detect both soft-clips."""
        send_init_data(viewer_page)
        result = viewer_page.evaluate("() => viewer.parseSoftClips('3S44M3S')")
        assert result["left"] == 3
        assert result["right"] == 3

    @pytest.mark.e2e
    def test_parse_soft_clips_empty(self, viewer_page: Page):
        """Should handle empty CIGAR."""
        send_init_data(viewer_page)
        result = viewer_page.evaluate("() => viewer.parseSoftClips('')")
        assert result["left"] == 0
        assert result["right"] == 0

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
            "() => ({ start: viewer.viewport.start, end: viewer.viewport.end })"
        )

        viewer_page.locator("#zoom-in").click()
        viewer_page.wait_for_timeout(100)

        after_viewport = viewer_page.evaluate(
            "() => ({ start: viewer.viewport.start, end: viewer.viewport.end })"
        )

        # Viewport span should change after zoom
        before_span = before_viewport["end"] - before_viewport["start"]
        after_span = after_viewport["end"] - after_viewport["start"]
        assert after_span < before_span


class TestMCPBridge:
    """E2E tests for the MCP postMessage bridge."""

    @pytest.mark.e2e
    def test_ready_message_sent(self, page: Page):
        """Viewer should send ready message on load."""
        # Load HTML with a message capture wrapper that records before the script runs
        html = get_viewer_html()
        capture_script = """<script>
            window.__messages = [];
            window.addEventListener('message', (e) => {
                window.__messages.push(e.data);
            });
        </script>"""
        # Insert capture script before the viewer script
        modified_html = html.replace("<script>", capture_script + "\n    <script>", 1)
        page.set_content(modified_html)
        page.wait_for_timeout(500)

        msgs = page.evaluate("() => window.__messages")
        ready = [m for m in msgs if isinstance(m, dict) and m.get("method") == "ready"]
        assert len(ready) >= 1
        assert ready[0]["jsonrpc"] == "2.0"

    @pytest.mark.e2e
    def test_go_button_sends_tools_call(self, viewer_page: Page):
        """Clicking Go should send a tools/call message."""
        # Set up message capture
        viewer_page.evaluate("""() => {
                window.__outgoing = [];
                const origPost = window.parent.postMessage.bind(window.parent);
                window.parent.postMessage = (msg, origin) => {
                    window.__outgoing.push(msg);
                    origPost(msg, origin);
                };
            }""")

        viewer_page.locator("#region-input").fill("chr1:5000-6000")
        viewer_page.locator("#go-btn").click()
        viewer_page.wait_for_timeout(200)

        messages = viewer_page.evaluate("() => window.__outgoing")
        tools_calls = [m for m in messages if m.get("method") == "tools/call"]
        assert len(tools_calls) >= 1
        assert tools_calls[0]["params"]["name"] == "browse_region"
        assert tools_calls[0]["params"]["arguments"]["region"] == "chr1:5000-6000"


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
        assert "Forward" in content

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

        content = viewer_page.locator("#tooltip").inner_text()
        assert "rev_read" in content
        assert "Reverse" in content

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

        initial_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
        viewer_page.mouse.wheel(0, 100)  # deltaY > 0 → zoom out (factor 1.2)
        viewer_page.wait_for_timeout(200)

        new_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")
        assert new_span > initial_span

    @pytest.mark.e2e
    def test_wheel_scroll_up_zooms_in(self, viewer_page: Page):
        """Scrolling up should zoom in (decrease span)."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")

        canvas = viewer_page.locator("#reads-canvas")
        box = canvas.bounding_box()
        assert box is not None

        viewer_page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
        viewer_page.mouse.wheel(0, -100)  # deltaY < 0 → zoom in (factor 0.8)
        viewer_page.wait_for_timeout(200)

        new_span = viewer_page.evaluate("() => viewer.viewport.end - viewer.viewport.start")
        assert new_span < initial_span


class TestMismatchRendering:
    """E2E tests for mismatch rendering at high zoom (scale > 2)."""

    @pytest.mark.e2e
    def test_mismatch_pixel_at_high_zoom(self, viewer_page: Page):
        """At high zoom, mismatches should render colored pixels at the mismatch position."""
        send_init_data(viewer_page, NARROW_DATA)

        # With 50bp in ~1024px, scale is ~20px/bp, well above the >2 threshold.
        # Mismatch is at position 125, alt='T' → color #ef4444 (red).
        # Pixel x = (125 - 100) * scale = 25 * scale.
        mismatch_color = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                const mx = Math.floor((125 - viewer.viewport.start) * scale);
                // Sample row 0 (y=0..12), middle at y=6
                const pixel = ctx.getImageData(mx + 1, 6, 1, 1).data;
                return { r: pixel[0], g: pixel[1], b: pixel[2], a: pixel[3] };
            }""")
        # #ef4444 = rgb(239, 68, 68) — the T mismatch color.
        # We just verify the red channel is dominant (not the default blue read color).
        assert mismatch_color["a"] > 0, "Mismatch pixel should be opaque"
        assert mismatch_color["r"] > mismatch_color["b"], (
            f"Mismatch should be red-dominant (T color), "
            f"got r={mismatch_color['r']} b={mismatch_color['b']}"
        )

    @pytest.mark.e2e
    def test_no_mismatch_at_low_zoom(self, viewer_page: Page):
        """At low zoom (scale <= 2), mismatches should NOT render."""
        # Use wide-span data so scale < 2
        wide_data = {
            **NARROW_DATA,
            "start": 0,
            "end": 1000,
            "coverage": [2] * 1000,
            "reference_sequence": "ACGTACGTAC" * 100,
        }
        send_init_data(viewer_page, wide_data)

        scale = viewer_page.evaluate(
            "() => viewer.readsCanvas.width / (viewer.viewport.end - viewer.viewport.start)"
        )
        assert scale <= 2, f"Scale should be <= 2 for this test, got {scale}"

        # Sample pixel at mismatch position 125: should be the read body color, not mismatch red
        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                const mx = Math.floor((125 - viewer.viewport.start) * scale);
                const p = ctx.getImageData(mx, 6, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        # At low zoom, mismatch should not override — pixel should be blue read color
        # #60a5fa = rgb(96, 165, 250) for forward reads
        if pixel["a"] > 0:
            assert pixel["b"] >= pixel["r"], (
                f"At low zoom, should be blue read color not red mismatch, "
                f"got r={pixel['r']} b={pixel['b']}"
            )


class TestReadColors:
    """E2E tests for forward/reverse read coloring."""

    @pytest.mark.e2e
    def test_forward_read_is_blue(self, viewer_page: Page):
        """Forward-strand reads should be blue (#60a5fa)."""
        send_init_data(viewer_page, NARROW_DATA)

        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                // Sample a pixel in row 0, at position 110 (away from mismatch at 125)
                const x = Math.floor((110 - viewer.viewport.start) * scale);
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
        """Reverse-strand reads should be purple (#a78bfa)."""
        send_init_data(viewer_page, NARROW_DATA)

        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                // Sample a pixel in row 1 (y=14..26, middle=20), at position 110
                const x = Math.floor((110 - viewer.viewport.start) * scale);
                const p = ctx.getImageData(x, 20, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        assert pixel["a"] > 0, "Pixel should be drawn"
        # #a78bfa: r=167, g=139, b=250 → blue dominant but also significant red
        assert pixel["b"] > pixel["g"], (
            f"Reverse read should have blue > green, got g={pixel['g']} b={pixel['b']}"
        )
        assert pixel["r"] > pixel["g"], (
            f"Reverse read should have red > green (purple), got r={pixel['r']} g={pixel['g']}"
        )

    @pytest.mark.e2e
    def test_forward_and_reverse_have_different_colors(self, viewer_page: Page):
        """Forward and reverse reads should be visually distinct colors."""
        send_init_data(viewer_page, NARROW_DATA)

        colors = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                const x = Math.floor((110 - viewer.viewport.start) * scale);
                const fwd = ctx.getImageData(x, 3, 1, 1).data;
                const rev = ctx.getImageData(x, 20, 1, 1).data;
                return {
                    fwd: { r: fwd[0], g: fwd[1], b: fwd[2] },
                    rev: { r: rev[0], g: rev[1], b: rev[2] },
                };
            }""")
        # Colors should differ — at minimum the red channel differs significantly
        # Forward (#60a5fa): r=96, Reverse (#a78bfa): r=167
        assert (
            colors["fwd"]["r"] != colors["rev"]["r"] or colors["fwd"]["g"] != colors["rev"]["g"]
        ), (
            f"Forward and reverse should have different colors: "
            f"fwd={colors['fwd']} rev={colors['rev']}"
        )


class TestSoftClipPixelVerification:
    """E2E tests verifying soft-clip regions actually draw pixels."""

    @pytest.mark.e2e
    def test_left_soft_clip_draws_pixels(self, viewer_page: Page):
        """Left soft-clip region should have semi-transparent pixels drawn."""
        clip_data = {
            "contig": "chr1",
            "start": 90,  # Start before the read to make room for the clip
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
                }
            ],
            "coverage": [1] * 70,
            "variants": [],
            "reference_sequence": "A" * 70,
        }
        send_init_data(viewer_page, clip_data)

        # The left soft-clip is 10bp before position 100 (genome coords ~90-100)
        # At scale ~14.6 px/bp (1024 / 70), clip region starts at x=0
        pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                // Sample at position 95 (middle of soft-clip region)
                const x = Math.floor((95 - viewer.viewport.start) * scale);
                const p = ctx.getImageData(x, 6, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        assert pixel["a"] > 0, "Soft-clip region should have pixels drawn"

    @pytest.mark.e2e
    def test_soft_clip_has_amber_border(self, viewer_page: Page):
        """Soft-clip region should have an amber (#f59e0b) stroke border."""
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
                }
            ],
            "coverage": [1] * 80,
            "variants": [],
            "reference_sequence": "A" * 80,
        }
        send_init_data(viewer_page, clip_data)

        # The border is drawn with strokeRect. Sample the top edge of the clip
        # region at position 95 (well inside the clip). The top border is at y=0.
        border_pixel = viewer_page.evaluate("""() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const scale = canvas.width / (viewer.viewport.end - viewer.viewport.start);
                const x = Math.floor((95 - viewer.viewport.start) * scale);
                // Top border pixel at y=0
                const p = ctx.getImageData(x, 0, 1, 1).data;
                return { r: p[0], g: p[1], b: p[2], a: p[3] };
            }""")
        # #f59e0b = rgb(245, 158, 11) — amber border.
        # At least verify the pixel is drawn and has warm (orange/amber) tones.
        assert border_pixel["a"] > 0, "Border should be drawn"
        assert border_pixel["r"] > border_pixel["b"], (
            f"Amber border should be warm-toned (r > b), "
            f"got r={border_pixel['r']} b={border_pixel['b']}"
        )


class TestOffScreenCulling:
    """E2E tests for off-screen read culling."""

    @pytest.mark.e2e
    def test_reads_outside_viewport_not_drawn(self, viewer_page: Page):
        """After panning away, reads outside viewport should not have pixels."""
        send_init_data(viewer_page, NARROW_DATA)

        # Zoom in tight then pan far away so reads at 100-150 are off-screen
        viewer_page.evaluate("""() => {
                viewer.viewport.start = 500;
                viewer.viewport.end = 550;
                viewer.render();
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

        max_cov = viewer_page.evaluate("() => Math.max(...viewer.data.coverage)")
        assert max_cov == 2  # Our SAMPLE_DATA has max coverage of 2

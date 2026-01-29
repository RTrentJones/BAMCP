"""End-to-end Playwright tests for the BAMCP viewer UI."""

import json
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
    elif 70 <= i < 100:
        SAMPLE_DATA["coverage"][i] = 1
    elif 100 <= i < 150:
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
        cov_width = viewer_page.locator("#coverage-canvas").evaluate(
            "el => el.width"
        )
        reads_width = viewer_page.locator("#reads-canvas").evaluate(
            "el => el.width"
        )
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
        initial_span = viewer_page.evaluate(
            "() => viewer.viewport.end - viewer.viewport.start"
        )

        viewer_page.locator("#zoom-in").click()

        new_span = viewer_page.evaluate(
            "() => viewer.viewport.end - viewer.viewport.start"
        )

        assert new_span < initial_span

    @pytest.mark.e2e
    def test_zoom_out_button(self, viewer_page: Page):
        """Zoom out button should widen the viewport."""
        send_init_data(viewer_page)

        initial_span = viewer_page.evaluate(
            "() => viewer.viewport.end - viewer.viewport.start"
        )

        viewer_page.locator("#zoom-out").click()

        new_span = viewer_page.evaluate(
            "() => viewer.viewport.end - viewer.viewport.start"
        )

        assert new_span > initial_span

    @pytest.mark.e2e
    def test_zoom_preserves_center(self, viewer_page: Page):
        """Zoom should preserve the center position."""
        send_init_data(viewer_page)

        initial_center = viewer_page.evaluate(
            "() => (viewer.viewport.start + viewer.viewport.end) / 2"
        )

        viewer_page.locator("#zoom-in").click()

        new_center = viewer_page.evaluate(
            "() => (viewer.viewport.start + viewer.viewport.end) / 2"
        )

        # Center should be approximately the same (within rounding)
        assert abs(new_center - initial_center) < 2

    @pytest.mark.e2e
    def test_variant_click_jumps_to_position(self, viewer_page: Page):
        """Clicking a variant row should jump to that position."""
        send_init_data(viewer_page)

        initial_start = viewer_page.evaluate("() => viewer.viewport.start")

        # Click the second variant row (position 250)
        viewer_page.locator("#variant-table tr").nth(1).click()

        new_center = viewer_page.evaluate(
            "() => (viewer.viewport.start + viewer.viewport.end) / 2"
        )

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

        no_overlap = viewer_page.evaluate(
            """() => {
                for (const row of viewer.packedRows) {
                    for (let i = 0; i < row.length - 1; i++) {
                        if (row[i].end_position >= row[i + 1].position) {
                            return false;
                        }
                    }
                }
                return true;
            }"""
        )
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
        has_content = viewer_page.evaluate(
            """() => {
                const canvas = document.getElementById('coverage-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }"""
        )
        assert has_content

    @pytest.mark.e2e
    def test_reads_canvas_drawn(self, viewer_page: Page):
        """Reads canvas should have pixels drawn after data load."""
        send_init_data(viewer_page)

        has_content = viewer_page.evaluate(
            """() => {
                const canvas = document.getElementById('reads-canvas');
                const ctx = canvas.getContext('2d');
                const data = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
                for (let i = 3; i < data.length; i += 4) {
                    if (data[i] > 0) return true;
                }
                return false;
            }"""
        )
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
        viewer_page.evaluate(
            """() => {
                window.__outgoing = [];
                const origPost = window.parent.postMessage.bind(window.parent);
                window.parent.postMessage = (msg, origin) => {
                    window.__outgoing.push(msg);
                    origPost(msg, origin);
                };
            }"""
        )

        viewer_page.locator("#region-input").fill("chr1:5000-6000")
        viewer_page.locator("#go-btn").click()
        viewer_page.wait_for_timeout(200)

        messages = viewer_page.evaluate("() => window.__outgoing")
        tools_calls = [m for m in messages if m.get("method") == "tools/call"]
        assert len(tools_calls) >= 1
        assert tools_calls[0]["params"]["name"] == "browse_region"
        assert tools_calls[0]["params"]["arguments"]["region"] == "chr1:5000-6000"

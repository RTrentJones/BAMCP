"""Unit tests for bamcp.resources module."""

import pytest

from bamcp.resources import get_viewer_html


class TestGetViewerHtml:
    """Tests for the viewer HTML resource."""

    @pytest.mark.unit
    def test_returns_html_string(self):
        """Should return an HTML string."""
        html = get_viewer_html()
        assert isinstance(html, str)
        assert len(html) > 0

    @pytest.mark.unit
    def test_valid_html_structure(self):
        """HTML should have basic structure."""
        html = get_viewer_html()
        assert "<!DOCTYPE html>" in html
        assert "<html>" in html
        assert "</html>" in html
        assert "<head>" in html
        assert "</head>" in html
        assert "<body>" in html
        assert "</body>" in html

    @pytest.mark.unit
    def test_contains_toolbar(self):
        """HTML should contain the toolbar elements."""
        html = get_viewer_html()
        assert 'id="toolbar"' in html
        assert 'id="region-input"' in html
        assert 'id="go-btn"' in html
        assert 'id="zoom-in"' in html
        assert 'id="zoom-out"' in html

    @pytest.mark.unit
    def test_contains_canvases(self):
        """HTML should contain the canvas elements."""
        html = get_viewer_html()
        assert 'id="coverage-canvas"' in html
        assert 'id="reads-canvas"' in html

    @pytest.mark.unit
    def test_contains_variant_panel(self):
        """HTML should contain the variant panel."""
        html = get_viewer_html()
        assert 'id="variant-panel"' in html
        assert 'id="variant-table"' in html
        assert "<th>Position</th>" in html
        assert "<th>Ref</th>" in html
        assert "<th>Alt</th>" in html
        assert "<th>VAF</th>" in html
        assert "<th>Depth</th>" in html

    @pytest.mark.unit
    def test_contains_base_colors(self):
        """HTML should define base colors for A, T, G, C, N."""
        html = get_viewer_html()
        assert "BASE_COLORS" in html
        assert "'A'" in html
        assert "'T'" in html
        assert "'G'" in html
        assert "'C'" in html
        assert "'N'" in html

    @pytest.mark.unit
    def test_contains_viewer_class(self):
        """HTML should define the BAMCPViewer class."""
        html = get_viewer_html()
        assert "class BAMCPViewer" in html

    @pytest.mark.unit
    def test_contains_mcp_bridge(self):
        """HTML should have the MCP postMessage bridge."""
        html = get_viewer_html()
        assert "setupMCPBridge" in html
        assert "postMessage" in html
        assert "bamcp/init" in html

    @pytest.mark.unit
    def test_contains_rendering_methods(self):
        """HTML should have rendering methods."""
        html = get_viewer_html()
        assert "renderCoverage" in html
        assert "renderReads" in html
        assert "renderRead" in html
        assert "renderVariantTable" in html

    @pytest.mark.unit
    def test_contains_navigation_methods(self):
        """HTML should have navigation methods."""
        html = get_viewer_html()
        assert "zoom" in html
        assert "pan" in html
        assert "jumpTo" in html
        assert "packReads" in html

    @pytest.mark.unit
    def test_contains_tooltip(self):
        """HTML should contain tooltip elements."""
        html = get_viewer_html()
        assert 'id="tooltip"' in html
        assert "showTooltip" in html

    @pytest.mark.unit
    def test_contains_soft_clip_rendering(self):
        """HTML should contain soft-clip rendering logic."""
        html = get_viewer_html()
        assert "parseSoftClips" in html
        assert "softClips" in html

    @pytest.mark.unit
    def test_css_styles_present(self):
        """HTML should include CSS styles."""
        html = get_viewer_html()
        assert "<style>" in html
        assert "#container" in html
        assert "#toolbar" in html
        assert "#viewer" in html

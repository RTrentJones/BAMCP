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
        # Column headers (may have class/data attributes)
        assert "Position" in html
        assert ">Ref<" in html
        assert ">Alt<" in html
        assert "VAF" in html
        assert "Depth" in html

    @pytest.mark.unit
    def test_contains_tooltip(self):
        """HTML should contain tooltip elements."""
        html = get_viewer_html()
        assert 'id="tooltip"' in html

    @pytest.mark.unit
    def test_css_styles_present(self):
        """HTML should include CSS styles."""
        html = get_viewer_html()
        assert "<style>" in html
        assert "#container" in html
        assert "#toolbar" in html
        assert "#viewer" in html

    @pytest.mark.unit
    def test_has_inline_script(self):
        """Bundled HTML should have inline JavaScript (not external CDN)."""
        html = get_viewer_html()
        assert '<script type="module"' in html
        # Should NOT have external CDN references
        assert "esm.sh" not in html
        assert "cdn.jsdelivr.net" not in html
        assert "unpkg.com" not in html

    @pytest.mark.unit
    def test_no_external_script_src(self):
        """Bundled HTML should not load scripts from external URLs."""
        html = get_viewer_html()
        # Check that we don't have src= pointing to http/https URLs
        import re

        external_scripts = re.findall(r'<script[^>]+src=["\']https?://', html)
        assert len(external_scripts) == 0, f"Found external script sources: {external_scripts}"

    @pytest.mark.unit
    def test_contains_mcp_app_functionality(self):
        """Bundled HTML should contain MCP Apps SDK functionality (even if minified)."""
        html = get_viewer_html()
        # These strings should survive minification since they're API method names
        assert "ontoolresult" in html
        assert "structuredContent" in html
        # connect() should be present (may be minified but method name preserved)
        assert ".connect(" in html

    @pytest.mark.unit
    def test_bundled_size_reasonable(self):
        """Bundled HTML should be a reasonable size (SDK is ~100KB bundled)."""
        html = get_viewer_html()
        # Should be > 50KB (SDK + viewer code) but < 1MB
        assert len(html) > 50_000, "HTML seems too small, SDK may not be bundled"
        assert len(html) < 1_000_000, "HTML seems too large"

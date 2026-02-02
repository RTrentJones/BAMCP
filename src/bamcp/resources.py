"""UI resource provider for the MCP Apps extension."""

import importlib.resources


def get_viewer_html() -> str:
    """Return the HTML content for the alignment viewer."""
    return (
        importlib.resources.files("bamcp")
        .joinpath("static/viewer.html")
        .read_text(encoding="utf-8")
    )

"""UI resource provider for the MCP Apps extension."""

import importlib.resources
from pathlib import Path

# Path to the static directory when running from source
_STATIC_DIR = Path(__file__).parent / "static"


def get_viewer_html() -> str:
    """Return the HTML content for the alignment viewer.

    Prefers the bundled version from dist/ if it exists (built with Vite),
    otherwise falls back to the source viewer.html.
    """
    # First try bundled version (has SDK inlined)
    dist_path = _STATIC_DIR / "dist" / "viewer.html"
    if dist_path.exists():
        return dist_path.read_text(encoding="utf-8")

    # Fall back to importlib.resources for installed package
    try:
        bundled = importlib.resources.files("bamcp").joinpath("static/dist/viewer.html")
        if bundled.is_file():
            return bundled.read_text(encoding="utf-8")
    except (TypeError, FileNotFoundError):
        pass

    # Final fallback: source viewer.html (won't work in sandboxed iframe)
    return (
        importlib.resources.files("bamcp")
        .joinpath("static/viewer.html")
        .read_text(encoding="utf-8")
    )

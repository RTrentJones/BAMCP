"""Viewer build stamp: content hash of the viewer sources recorded next to dist/.

Stdlib-only on purpose — ``make build-viewer`` runs this right after ``npm run build``, before
the Python dependencies are installed (CI's ``make install`` builds the viewer first), so it
must not import the rest of the package. ``python -m bamcp.viewer_build`` writes the stamp.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

_STATIC_DIR = Path(__file__).resolve().parent / "static"


# Written next to dist/viewer.html after a build (``--record-build-hash``): the content hash
# of the sources that produced it, plus the dist hash so a later un-stamped rebuild is detected.
BUILD_STAMP_NAME = ".viewer-sources.json"


def _viewer_source_files(static_dir: Path) -> list[Path]:
    """The files whose content determines dist/viewer.html (tests excluded)."""
    files = [f for f in static_dir.glob("*.ts") if not f.name.endswith(".test.ts")]
    files += [static_dir / "viewer.html", static_dir / "package-lock.json"]
    return sorted((f for f in files if f.is_file()), key=lambda f: f.name)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def viewer_source_hash(static_dir: Path | None = None) -> str:
    """Content hash of the viewer sources (names + bytes) — independent of mtimes."""
    static_dir = static_dir or _STATIC_DIR
    h = hashlib.sha256()
    for f in _viewer_source_files(static_dir):
        h.update(f.name.encode("utf-8") + b"\0")
        h.update(f.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def record_build_hash(static_dir: Path | None = None) -> Path:
    """Stamp dist/ with the source + dist hashes. Run right after ``npm run build``."""
    static_dir = static_dir or _STATIC_DIR
    dist = static_dir / "dist" / "viewer.html"
    if not dist.exists():
        raise SystemExit(f"{dist} does not exist; build the viewer before recording its hash.")
    stamp = static_dir / "dist" / BUILD_STAMP_NAME
    stamp.write_text(
        json.dumps({"sources": viewer_source_hash(static_dir), "dist": _sha256_file(dist)}) + "\n",
        encoding="utf-8",
    )
    return stamp


def _read_stamp(static_dir: Path, dist: Path) -> str | None:
    """The recorded source hash, or None when there is no stamp matching this dist build."""
    stamp = static_dir / "dist" / BUILD_STAMP_NAME
    try:
        data = json.loads(stamp.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or data.get("dist") != _sha256_file(dist):
        return None  # dist was rebuilt after the stamp was written — the stamp says nothing.
    sources = data.get("sources")
    return sources if isinstance(sources, str) else None


if __name__ == "__main__":  # pragma: no cover — exercised by `make build-viewer`
    print(f"Recorded viewer build hash -> {record_build_hash()}")

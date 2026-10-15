"""Developer CLI implementation for the visual feedback loop.

The thin shell at ``scripts/render_viewer.py`` delegates to ``main()`` here.
Pulled into the package so tests can ``from bamcp.eval.render_cli import main``
without sys.path gymnastics.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import re
import sys
import time
from pathlib import Path

from ..config import BAMCPConfig
from ..core.serialization import serialize_region_data
from ..core.tools import _fetch_region_with_timeout
from .renderer import ViewerRenderer

_VARIANT_PATTERN = re.compile(r"^([^:]+):(\d+):([ACGTN]+)>([ACGTN]+)$", re.IGNORECASE)

_DEFAULT_HEIGHT = 800

_STATIC_DIR = Path(__file__).resolve().parent.parent / "static"


# Re-exported: the stamp logic lives in the stdlib-only bamcp.viewer_build (it runs before
# the package's dependencies are installed).
from ..viewer_build import (  # noqa: E402
    BUILD_STAMP_NAME as BUILD_STAMP_NAME,
)
from ..viewer_build import (  # noqa: E402
    _read_stamp,
    _viewer_source_files,
    record_build_hash,
    viewer_source_hash,
)


def _check_rebuild_needed() -> str | None:
    """Return a warning string when dist/viewer.html was not built from the current sources.

    Compares a content hash of the sources against the hash recorded at build time
    (``--record-build-hash``), so a fresh checkout or a branch switch that only touches
    mtimes does not warn. Without a stamp it falls back to the old mtime comparison.
    """
    dist = _STATIC_DIR / "dist" / "viewer.html"
    rebuild = "Run 'cd src/bamcp/static && npm run build'"
    if not dist.exists():
        return f"WARNING: {dist} does not exist. {rebuild}."
    recorded = _read_stamp(_STATIC_DIR, dist)
    if recorded is not None:
        if recorded == viewer_source_hash(_STATIC_DIR):
            return None
        return (
            "WARNING: viewer sources changed since dist/viewer.html was built "
            f"(content hash mismatch). {rebuild} to pick up your changes."
        )
    dist_mtime = dist.stat().st_mtime
    newer = [f for f in _viewer_source_files(_STATIC_DIR) if f.stat().st_mtime > dist_mtime]
    if newer:
        names = ", ".join(f.name for f in newer)
        return (
            f"WARNING: viewer sources newer than dist/viewer.html ({names}). "
            f"{rebuild} to pick up your changes (then '--record-build-hash' "
            "to switch this check to content hashes)."
        )
    return None


async def _build_payload_from_bam(
    bam_path: str, reference: str | None, region: str, config: BAMCPConfig
) -> dict:
    """Run the same parse + serialize the MCP server uses."""
    data = await _fetch_region_with_timeout(bam_path, region, reference, config)
    payload = serialize_region_data(data)
    payload["file_path"] = bam_path
    return payload


def _load_ui_init_json(path: str) -> dict:
    raw = Path(path).read_text(encoding="utf-8")
    payload = json.loads(raw)
    if not isinstance(payload, dict):
        raise ValueError(
            f"--ui-init-json file must contain a JSON object, got {type(payload).__name__}"
        )
    return payload


def _parse_active_variant(spec: str | None) -> tuple[int | None, str | None]:
    if not spec:
        return None, None
    match = _VARIANT_PATTERN.match(spec)
    if not match:
        raise SystemExit(f"--active-variant must look like 'chrom:pos:ref>alt', got {spec!r}")
    _chrom, pos, _ref, alt = match.groups()
    return int(pos), alt.upper()


def _resolve_region(bam_path: str, region: str | None) -> str:
    if region:
        return region
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contig = bam.references[0]
        length = bam.lengths[0]
    end = min(length, 1000)
    return f"{contig}:1-{end}"


def _parse_height(value: str) -> int | str:
    if value.lower() == "auto":
        return "auto"
    try:
        height = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--height must be an integer or 'auto', got {value!r}"
        ) from None
    if height <= 0:
        raise argparse.ArgumentTypeError("--height must be positive")
    return height


def _resolve_height(args: argparse.Namespace) -> int | None:
    """Explicit --height wins; otherwise 800, or auto-fit (None) for a blinded render."""
    if args.height == "auto":
        return None
    if args.height is not None:
        return int(args.height)
    return None if args.blind else _DEFAULT_HEIGHT


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="render_viewer",
        description="Render the BAMCP viewer to a PNG for visual verification.",
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--bam", help="Path to a BAM/CRAM file. The region is parsed live.")
    src.add_argument(
        "--ui-init-json",
        help="Path to a JSON file containing a pre-serialized ui/init payload.",
    )
    p.add_argument("--reference", help="Reference FASTA (only used with --bam).")
    p.add_argument(
        "--region",
        help=(
            "Genomic region 'chrom:start-end' (only used with --bam). "
            "Defaults to the first 1000bp of the first contig."
        ),
    )
    p.add_argument(
        "--mode",
        default="expanded",
        choices=["squished", "compact", "expanded", "dv-strips", "dv-composite"],
        help="Display mode to render at.",
    )
    p.add_argument(
        "--active-variant",
        help=(
            "Mark a variant as active for the supports-variant channel, format 'chrom:pos:ref>alt'."
        ),
    )
    p.add_argument("--width", type=int, default=1024, help="Viewport width in pixels.")
    p.add_argument(
        "--height",
        type=_parse_height,
        default=None,
        help=(
            "Viewport height in pixels, or 'auto' to grow until every read row fits. "
            "Default: 800, or 'auto' with --blind."
        ),
    )
    p.add_argument(
        "--blind",
        action="store_true",
        help=(
            "Blinded eval stimulus: no variant table/overlays/verdicts, hashed read names, "
            "pinned display settings, clipped to #viewer, depth-sized height."
        ),
    )
    p.add_argument(
        "--clip",
        choices=["#viewer", "#container"],
        default=None,
        help="Element to screenshot. Default: #container, or #viewer with --blind.",
    )
    p.add_argument("--out", default=".render-dev/last.png", help="Output PNG path.")
    p.epilog = (
        "Run 'render_viewer --record-build-hash' right after 'npm run build' to stamp dist/ "
        "with a content hash of the sources (enables the hash-based staleness check)."
    )
    p.add_argument(
        "--print-path-only",
        action="store_true",
        help="Print only the output path on stdout. Useful for piping into Read.",
    )
    return p.parse_args(argv)


async def _async_main(args: argparse.Namespace) -> int:
    config = BAMCPConfig(reference=args.reference)

    if args.bam:
        region = _resolve_region(args.bam, args.region)
        payload = await _build_payload_from_bam(args.bam, args.reference, region, config)
    else:
        payload = _load_ui_init_json(args.ui_init_json)

    pos, alt = _parse_active_variant(args.active_variant)
    if args.blind and (pos is not None or args.clip == "#container"):
        raise SystemExit("--blind cannot be combined with --active-variant or --clip '#container'")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.monotonic()
    height = _resolve_height(args)
    async with ViewerRenderer(
        viewport_width=args.width, viewport_height=height or _DEFAULT_HEIGHT
    ) as renderer:
        png_bytes = await renderer.capture(
            payload,
            args.mode,
            blind=args.blind,
            viewport=(args.width, height),
            clip=args.clip,
            active_variant_position=pos,
            active_variant_alt=alt,
        )
    out_path.write_bytes(png_bytes)
    elapsed_ms = (time.monotonic() - start) * 1000

    if args.print_path_only:
        print(str(out_path))
    else:
        print(
            f"Rendered {len(png_bytes):,} bytes to {out_path} "
            f"(mode={args.mode}{', blind' if args.blind else ''}, {elapsed_ms:.0f}ms)"
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    raw = sys.argv[1:] if argv is None else argv
    if raw == ["--record-build-hash"]:
        print(record_build_hash())
        return 0
    args = _parse_args(argv)
    if not args.print_path_only:
        warning = _check_rebuild_needed()
        if warning:
            print(warning, file=sys.stderr)
    return asyncio.run(_async_main(args))

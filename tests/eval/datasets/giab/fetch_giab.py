"""Fetch a GIAB truth slice and emit a BAMCP truth-set manifest.

Scaffold for the real-data eval path. Network-gated and intentionally NOT run
in CI — the synthetic_v1 gate covers regressions on every PR; GIAB is a
manual / nightly job. See ./README.md for the full rationale.

The functions below define the contract (URLs, checksum verification, manifest
shape) so the milestone work is "fill in the slice + VCF parse", not "design
the pipeline". Each step that requires the actual download raises
NotImplementedError with a precise description of what to implement.
"""

from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass
from pathlib import Path

# GIAB v4.2.1 release for HG001 / NA12878 on GRCh38.
GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/"
)
TRUTH_VCF = GIAB_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HIGHCONF_BED = GIAB_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"


@dataclass(frozen=True)
class GiabSlice:
    """A small region to score, keeping the run fast and cheap."""

    region: str = "chr20:1000000-1300000"
    output_dir: Path = Path("tests/eval/datasets/giab")


def sha256(path: Path) -> str:
    """Stream a file through SHA-256 (used to verify downloaded artifacts)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> Path:  # pragma: no cover — network
    """Download ``url`` to ``dest`` if missing. Implement with urllib/requests."""
    raise NotImplementedError(
        f"Download {url} -> {dest}. Stream to disk, then verify with sha256(); "
        "GIAB publishes checksums alongside each release artifact."
    )


def build_manifest(slice_: GiabSlice) -> Path:  # pragma: no cover — network
    """Emit a truth-set manifest for the slice, matching synthetic_v1's schema.

    Implementation outline:
      1. download(TRUTH_VCF) and download(HIGHCONF_BED); verify checksums.
      2. Restrict to ``slice_.region`` (e.g. via pysam.VariantFile fetch).
      3. variant_sites <- truth VCF records inside the high-confidence BED
         (expected_artifact / is_clean_control left null).
      4. negative_regions <- high-confidence stretches with no truth variant.
      5. Write ``slice_.output_dir / "manifest.yaml"``.
    """
    raise NotImplementedError(
        "Wire the slice + VCF parse per the steps in this function's docstring "
        "and ./README.md. The scorer (bamcp.eval.truthset) then runs unchanged."
    )


def main(argv: list[str] | None = None) -> int:  # pragma: no cover — scaffold
    p = argparse.ArgumentParser(description="Fetch a GIAB truth slice (scaffold).")
    p.add_argument("--region", default=GiabSlice.region)
    args = p.parse_args(argv)
    build_manifest(GiabSlice(region=args.region))
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())

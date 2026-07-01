"""Generate a real Genome-in-a-Bottle benchmark slice for the BAMCP eval.

This turns the eval's real-data path from a scaffold into something that runs.
It is **semi-synthetic over real biology**:

- Real GRCh38 reference (UCSC hg38 chr20).
- Real NIST GIAB v4.2.1 truth variants — real coordinates, alleles, and
  zygosity — restricted to the high-confidence BED.
- Simulated short reads at a controlled depth with an Illumina-like base-error
  rate, honoring each truth site's genotype (het ~50% alt, hom ~100% alt).

What it measures: whether BAMCP's caller recovers real variant structures on
real sequence context (recall) and rejects sequencing noise (precision on the
high-confidence intervals). What it does NOT measure: robustness to real
alignment/mapping artifacts, mismapping in segmental duplications, or real
error profiles — those require real reads (a ~300 GB alignment), which is out of
scope for a portable, reproducible benchmark. See GIAB_RESULTS.md for the
honest framing of what the numbers do and don't say.

The emitted manifest uses the same schema as ``synthetic_v1/manifest.yaml``, so
the existing scorer runs it unchanged::

    python tests/eval/datasets/giab/fetch_giab.py --region chr20:1000000-1060000
    python -m bamcp.eval.truthset \\
      --manifest tests/eval/datasets/giab/manifest.yaml \\
      --min-recall 0.90 --min-precision 0.90

Network-gated (streams the truth VCF, downloads ~15 MB of reference + BED).
Not run in CI — the synthetic_v1 gate covers per-PR regressions.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import os
import random
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

# GIAB v4.2.1 benchmark for HG001 / NA12878 on GRCh38 (chr-prefixed).
_GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/"
)
TRUTH_VCF = _GIAB_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HIGHCONF_BED = _GIAB_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"

# Per-chromosome GRCh38 reference (much smaller than the whole genome).
_UCSC_CHROM = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"

_HERE = Path(__file__).resolve().parent
_BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class TruthSNV:
    """A single real GIAB truth SNV, 0-based (matching BAMCP's coordinates)."""

    pos0: int  # 0-based reference position
    ref: str
    alt: str
    zygosity: str  # "het" or "hom"


# ── Downloads (curl: it already works through the sandbox proxy) ───────────


def _curl(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url} -> {dest.name}")
    curl = shutil.which("curl") or "curl"
    # Fixed argument vector, no shell; URLs are the module's own constants.
    subprocess.run(  # noqa: S603
        [curl, "-sSfL", "--retry", "3", "-o", str(dest), url],
        check=True,
        timeout=1800,
    )
    return dest


# ── High-confidence BED ────────────────────────────────────────────────────


def load_bed_intervals(bed_path: Path, chrom: str, start: int, end: int) -> list[tuple[int, int]]:
    """Return sorted [start,end) intervals on ``chrom`` overlapping [start,end)."""
    intervals: list[tuple[int, int]] = []
    with bed_path.open() as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            parts = line.split("\t")
            b_start, b_end = int(parts[1]), int(parts[2])
            if b_end <= start or b_start >= end:
                continue
            intervals.append((max(b_start, start), min(b_end, end)))
    intervals.sort()
    return intervals


def _in_intervals(pos: int, starts: list[int], intervals: list[tuple[int, int]]) -> bool:
    """True if ``pos`` falls in any [start,end) interval (starts is the sorted keys)."""
    i = bisect.bisect_right(starts, pos) - 1
    return 0 <= i < len(intervals) and intervals[i][0] <= pos < intervals[i][1]


# ── Real truth variants (streamed) ─────────────────────────────────────────


def fetch_truth_snvs(
    chrom: str, start: int, end: int, intervals: list[tuple[int, int]]
) -> list[TruthSNV]:
    """Stream the GIAB truth VCF and return PASS biallelic SNVs inside the BED."""
    import pysam  # local import: pysam is heavy and only needed here

    starts = [iv[0] for iv in intervals]
    out: list[TruthSNV] = []
    vf = pysam.VariantFile(TRUTH_VCF)
    try:
        for rec in vf.fetch(chrom, start, end):
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref, alt = rec.ref, rec.alts[0]
            if len(ref) != 1 or len(alt) != 1:
                continue  # SNVs only — indel read-simulation is future work
            if ref.upper() not in _BASES or alt.upper() not in _BASES:
                continue
            pos0 = rec.pos - 1  # VCF is 1-based; BAMCP reports 0-based
            if not _in_intervals(pos0, starts, intervals):
                continue
            gt = rec.samples[0].get("GT") if rec.samples else None
            zygosity = "hom" if gt == (1, 1) else "het"
            out.append(TruthSNV(pos0=pos0, ref=ref.upper(), alt=alt.upper(), zygosity=zygosity))
    finally:
        vf.close()
    return out


# ── Reference ──────────────────────────────────────────────────────────────


def ensure_reference(chrom: str, data_dir: Path) -> Path:
    """Download + index the GRCh38 chromosome FASTA. Returns the .fa path."""
    import pysam

    fa = data_dir / f"{chrom}.fa"
    if not (fa.with_suffix(".fa.fai")).exists():
        gz = _curl(_UCSC_CHROM.format(chrom=chrom), data_dir / f"{chrom}.fa.gz")
        print(f"  decompressing {gz.name}")
        with gzip.open(gz, "rb") as src, fa.open("wb") as dst:
            dst.write(src.read())
        pysam.faidx(str(fa))
    return fa


# ── Read simulation ─────────────────────────────────────────────────────────


def simulate_reads(
    ref_fa: Path,
    chrom: str,
    start: int,
    end: int,
    truth: list[TruthSNV],
    out_bam: Path,
    *,
    depth: int,
    read_len: int,
    error_rate: float,
    seed: int,
) -> Path:
    """Write an indexed BAM of simulated reads honoring the truth genotypes.

    Deterministic given ``seed``. Het sites carry the alt on ~50% of covering
    reads; hom sites on ~100%. A uniform per-base substitution error is applied
    at ``error_rate`` so the caller's noise-rejection is exercised.
    """
    import pysam

    rng = random.Random(seed)  # noqa: S311 — read simulation, not cryptographic
    fa = pysam.FastaFile(str(ref_fa))
    chrom_len = fa.get_reference_length(chrom)
    truth_by_pos = {t.pos0: t for t in truth}

    qual = _phred_for_error_rate(error_rate)
    n_reads = depth * (end - start) // read_len

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": chrom_len}]}
    segments = []
    for i in range(n_reads):
        r_start = rng.randint(start, max(start, end - read_len))
        seq = list(fa.fetch(chrom, r_start, r_start + read_len).upper())
        # Apply real truth alleles.
        for j in range(len(seq)):
            ref_pos = r_start + j
            t = truth_by_pos.get(ref_pos)
            if t is None:
                continue
            carries_alt = t.zygosity == "hom" or rng.random() < 0.5
            if carries_alt:
                seq[j] = t.alt
        # Apply uniform sequencing error.
        if error_rate > 0:
            for j in range(len(seq)):
                if rng.random() < error_rate:
                    seq[j] = rng.choice([b for b in _BASES if b != seq[j]])
        seg = pysam.AlignedSegment()
        seg.query_name = f"giab_sim_{i}"
        seg.query_sequence = "".join(seq)
        seg.flag = 0 if i % 2 == 0 else 16
        seg.reference_id = 0
        seg.reference_start = r_start
        seg.mapping_quality = 60
        seg.cigartuples = [(0, read_len)]
        seg.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * read_len)
        segments.append(seg)
    fa.close()

    segments.sort(key=lambda s: s.reference_start)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(out_bam), "wb", header=header) as outf:
        for seg in segments:
            outf.write(seg)
    pysam.index(str(out_bam))
    return out_bam


def _phred_for_error_rate(error_rate: float) -> int:
    """Phred quality consistent with the injected error rate (capped 2..40)."""
    if error_rate <= 0:
        return 40
    import math

    return max(2, min(40, round(-10 * math.log10(error_rate))))


# ── Manifest ────────────────────────────────────────────────────────────────


def _pick_negative_regions(
    intervals: list[tuple[int, int]], truth: list[TruthSNV], want: int = 3, span: int = 2000
) -> list[tuple[int, int]]:
    """Find high-confidence sub-regions containing no truth SNV (for FP control)."""
    truth_positions = sorted(t.pos0 for t in truth)
    negatives: list[tuple[int, int]] = []
    for iv_start, iv_end in intervals:
        cursor = iv_start
        while cursor + span <= iv_end and len(negatives) < want:
            lo = bisect.bisect_left(truth_positions, cursor)
            hi = bisect.bisect_left(truth_positions, cursor + span)
            if lo == hi:  # no truth variant in [cursor, cursor+span)
                negatives.append((cursor, cursor + span))
                cursor += span
            else:
                cursor = truth_positions[hi - 1] + 1
        if len(negatives) >= want:
            break
    return negatives


def write_manifest(
    out_path: Path,
    chrom: str,
    start: int,
    end: int,
    truth: list[TruthSNV],
    negatives: list[tuple[int, int]],
    ref_path: Path,
    bam_path: Path,
    *,
    depth: int,
    error_rate: float,
    seed: int,
) -> Path:
    def _rel(p: Path) -> str:
        # Repo-relative when possible so the committed manifest is portable.
        try:
            return os.path.relpath(p, Path.cwd())
        except ValueError:  # pragma: no cover — different drive on Windows
            return str(p)

    lines = [
        "# GIAB benchmark slice — GENERATED, do not hand-edit.",
        "# Real GRCh38 reference + real NIST v4.2.1 truth SNVs (real coords,",
        "# alleles, zygosity) + simulated reads. See fetch_giab.py / GIAB_RESULTS.md.",
        "",
        "dataset: giab_hg001_grch38_v4.2.1",
        "version: 1",
        "description: >-",
        f"  {len(truth)} real NIST truth SNVs on {chrom}:{start}-{end} inside the",
        f"  high-confidence BED, simulated at {depth}x with base-error rate {error_rate}.",
        "",
        "provenance:",
        f"  generator: tests/eval/datasets/giab/fetch_giab.py (seed={seed})",
        f"  bam: {_rel(bam_path)}",
        f"  reference: {_rel(ref_path)}",
        f"  truth_vcf: {TRUTH_VCF}",
        f"  highconf_bed: {HIGHCONF_BED}",
        "  note: >-",
        "    Semi-synthetic. Measures recovery of real variant structures and",
        "    noise rejection; NOT robustness to real alignment/mapping artifacts.",
        "",
        f"detection_region: {chrom}:{start}-{end}",
        "",
        "variant_sites:",
    ]
    for t in truth:
        lines.append(
            f"  - {{chrom: {chrom}, pos: {t.pos0}, ref: {t.ref}, alt: {t.alt}, "
            f"category: {t.zygosity}}}"
        )
    lines.append("")
    lines.append("negative_regions:")
    if negatives:
        for n_start, n_end in negatives:
            lines.append(
                f"  - {{region: '{chrom}:{n_start}-{n_end}', note: high-confidence, no truth SNV}}"
            )
    else:
        lines.append("  []")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


# ── Orchestration ───────────────────────────────────────────────────────────


def _parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":")
    start, end = span.split("-")
    return chrom, int(start), int(end)


def generate(
    region: str,
    *,
    depth: int = 30,
    read_len: int = 100,
    error_rate: float = 0.005,
    seed: int = 1,
    data_dir: Path | None = None,
) -> dict[str, object]:
    """Build the reference, reads, and manifest for ``region``. Returns a summary."""
    chrom, start, end = _parse_region(region)
    data_dir = data_dir or (_HERE / "data")
    data_dir.mkdir(parents=True, exist_ok=True)

    print("1/5 high-confidence BED")
    bed = _curl(HIGHCONF_BED, data_dir / "highconf.bed")
    intervals = load_bed_intervals(bed, chrom, start, end)
    if not intervals:
        raise SystemExit(f"No high-confidence intervals in {region}")

    print("2/5 truth SNVs (streaming VCF)")
    truth = fetch_truth_snvs(chrom, start, end, intervals)
    if not truth:
        raise SystemExit(f"No truth SNVs in {region}")

    print("3/5 reference")
    ref = ensure_reference(chrom, data_dir)

    print(f"4/5 simulating reads ({depth}x, err={error_rate})")
    bam = simulate_reads(
        ref,
        chrom,
        start,
        end,
        truth,
        data_dir / "reads.bam",
        depth=depth,
        read_len=read_len,
        error_rate=error_rate,
        seed=seed,
    )

    print("5/5 manifest")
    negatives = _pick_negative_regions(intervals, truth)
    manifest = write_manifest(
        _HERE / "manifest.yaml",
        chrom,
        start,
        end,
        truth,
        negatives,
        ref,
        bam,
        depth=depth,
        error_rate=error_rate,
        seed=seed,
    )

    summary = {
        "region": region,
        "truth_snvs": len(truth),
        "het": sum(1 for t in truth if t.zygosity == "het"),
        "hom": sum(1 for t in truth if t.zygosity == "hom"),
        "negative_regions": len(negatives),
        "manifest": str(manifest),
        "bam": str(bam),
        "reference": str(ref),
    }
    print("done:", summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Generate a real GIAB benchmark slice.")
    p.add_argument("--region", default="chr20:1000000-1060000")
    p.add_argument("--depth", type=int, default=30)
    p.add_argument("--read-len", type=int, default=100)
    p.add_argument("--error-rate", type=float, default=0.005)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args(argv)
    generate(
        args.region,
        depth=args.depth,
        read_len=args.read_len,
        error_rate=args.error_rate,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover — network-gated CLI
    raise SystemExit(main())

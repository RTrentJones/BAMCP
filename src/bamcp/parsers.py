"""BAM/CRAM file parsing using pysam."""

from dataclasses import dataclass, field

import pysam


@dataclass
class AlignedRead:
    """A single aligned read extracted from a BAM/CRAM file."""

    name: str
    sequence: str
    qualities: list[int]
    cigar: str
    position: int
    end_position: int
    mapping_quality: int
    is_reverse: bool
    mismatches: list[dict] = field(default_factory=list)


@dataclass
class RegionData:
    """Data for a genomic region including reads, coverage, and variants."""

    contig: str
    start: int
    end: int
    reads: list[AlignedRead]
    coverage: list[int]
    variants: list[dict]
    reference_sequence: str | None = None


def parse_region(region: str) -> tuple[str, int, int]:
    """
    Parse a genomic region string into contig, start, end.

    Supports formats:
        - chr1:1000-2000
        - chr1:1,000-2,000
        - 1:1000-2000

    Returns:
        Tuple of (contig, start, end)

    Raises:
        ValueError: If region format is invalid.
    """
    region = region.replace(",", "")
    try:
        contig, coords = region.split(":")
        start_str, end_str = coords.split("-")
        start = int(start_str)
        end = int(end_str)
    except (ValueError, AttributeError) as e:
        raise ValueError(
            f"Invalid region format: '{region}'. Expected format: 'chr1:1000-2000'"
        ) from e

    if start < 0:
        raise ValueError(f"Start position must be non-negative, got {start}")
    if end <= start:
        raise ValueError(f"End position ({end}) must be greater than start ({start})")

    return contig, start, end


def fetch_region(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    max_reads: int = 10000,
    min_mapq: int = 0,
    index_filename: str | None = None,
) -> RegionData:
    """
    Fetch reads from a BAM/CRAM file for a given region.

    Args:
        bam_path: Path to BAM/CRAM file (local or remote).
        region: Genomic region (e.g., "chr1:1000-2000").
        reference_path: Path to reference FASTA (required for CRAM).
        max_reads: Maximum reads to return (stop after this many).
        min_mapq: Minimum mapping quality filter.
        index_filename: Path to index file (.bai/.crai). Used to redirect
            remote BAM index downloads to a cache directory.

    Returns:
        RegionData with reads, coverage, and detected variants.
    """
    contig, start, end = parse_region(region)

    mode = "rc" if bam_path.endswith(".cram") else "rb"
    samfile = pysam.AlignmentFile(
        bam_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference_path,
        index_filename=index_filename,
    )

    reads: list[AlignedRead] = []
    coverage = [0] * (end - start)
    base_counts: list[dict[str, int]] = [{} for _ in range(end - start)]

    ref_seq: str | None = None
    if reference_path:
        with pysam.FastaFile(reference_path) as fasta:
            ref_seq = fasta.fetch(contig, start, end)

    read_count = 0
    for read in samfile.fetch(contig, start, end):
        if read.mapping_quality < min_mapq:
            continue
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        read_count += 1
        if read_count > max_reads:
            break

        mismatches: list[dict] = []
        if ref_seq and read.query_sequence:
            try:
                aligned_pairs = read.get_aligned_pairs(with_seq=True)
            except ValueError:
                # MD tag not present; fall back to coordinate-based comparison
                aligned_pairs = []
                for qpos, rpos in read.get_aligned_pairs():
                    if rpos is not None and start <= rpos < end:
                        idx = rpos - start
                        ref_base_chr = ref_seq[idx] if 0 <= idx < len(ref_seq) else None
                        aligned_pairs.append((qpos, rpos, ref_base_chr))
                    else:
                        aligned_pairs.append((qpos, rpos, None))
            for qpos, rpos, ref_base in aligned_pairs:
                if rpos is None or qpos is None:
                    continue
                if start <= rpos < end:
                    query_base = read.query_sequence[qpos]
                    if ref_base and query_base != ref_base.upper():
                        mismatches.append({"pos": rpos, "ref": ref_base.upper(), "alt": query_base})
                    idx = rpos - start
                    if 0 <= idx < len(base_counts):
                        base_counts[idx][query_base] = base_counts[idx].get(query_base, 0) + 1

        ref_start = max(read.reference_start, start)
        ref_end = min(read.reference_end or read.reference_start, end)
        for pos in range(ref_start, ref_end):
            idx = pos - start
            if 0 <= idx < len(coverage):
                coverage[idx] += 1

        reads.append(
            AlignedRead(
                name=read.query_name or "",
                sequence=read.query_sequence or "",
                qualities=list(read.query_qualities or []),
                cigar=read.cigarstring or "",
                position=read.reference_start,
                end_position=read.reference_end or read.reference_start,
                mapping_quality=read.mapping_quality,
                is_reverse=read.is_reverse,
                mismatches=mismatches,
            )
        )

    samfile.close()

    variants = detect_variants(base_counts, ref_seq, contig, start)

    return RegionData(
        contig=contig,
        start=start,
        end=end,
        reads=reads,
        coverage=coverage,
        variants=variants,
        reference_sequence=ref_seq,
    )


def detect_variants(
    base_counts: list[dict[str, int]],
    ref_seq: str | None,
    contig: str,
    start: int,
    min_vaf: float = 0.1,
    min_depth: int = 10,
) -> list[dict]:
    """
    Detect variants from base counts.

    Args:
        base_counts: List of dicts mapping base -> count for each position.
        ref_seq: Reference sequence for the region.
        contig: Contig name.
        start: Start position of the region.
        min_vaf: Minimum variant allele frequency threshold.
        min_depth: Minimum read depth to consider.

    Returns:
        List of variant dicts with contig, position, ref, alt, vaf, depth, alt_count.
    """
    variants: list[dict] = []

    if not ref_seq:
        return variants

    for i, counts in enumerate(base_counts):
        if not counts:
            continue

        total = sum(counts.values())
        if total < min_depth:
            continue

        ref_base = ref_seq[i].upper()

        for base, count in counts.items():
            if base == ref_base:
                continue
            vaf = count / total
            if vaf >= min_vaf:
                variants.append(
                    {
                        "contig": contig,
                        "position": start + i,
                        "ref": ref_base,
                        "alt": base,
                        "vaf": round(vaf, 3),
                        "depth": total,
                        "alt_count": count,
                    }
                )

    return variants

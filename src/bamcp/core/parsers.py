"""BAM/CRAM file parsing using pysam."""

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pysam


@dataclass
class SoftClip:
    """A soft-clipped region from a read."""

    position: int  # Reference position where clip starts
    length: int  # Number of bases clipped
    sequence: str | None  # Clipped sequence if available
    side: str  # 'left' or 'right'


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
    soft_clips: list[SoftClip] = field(default_factory=list)

    # Paired-end fields
    mate_position: int | None = None
    mate_contig: str | None = None
    insert_size: int | None = None
    is_proper_pair: bool = False
    is_read1: bool = False
    is_paired: bool = False


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


# Maximum region size to prevent DoS via unbounded memory allocation
MAX_REGION_SIZE = 1_000_000  # 1 Mbp

# CIGAR operation codes
CIGAR_SOFT_CLIP = 4  # S operation


def extract_soft_clips(read: pysam.AlignedSegment) -> list[SoftClip]:
    """
    Extract soft clip positions and sequences from CIGAR.

    Args:
        read: A pysam aligned segment.

    Returns:
        List of SoftClip objects for left and right soft clips.
    """
    clips: list[SoftClip] = []
    cigar = read.cigartuples

    if not cigar:
        return clips

    query_seq = read.query_sequence
    ref_start = read.reference_start or 0
    query_pos = 0

    for op, length in cigar:
        if op == CIGAR_SOFT_CLIP:
            # Get clipped sequence
            clip_seq = None
            if query_seq and query_pos < len(query_seq):
                clip_seq = query_seq[query_pos : query_pos + length]

            # Determine if left or right clip
            side = "left" if query_pos == 0 else "right"

            clips.append(
                SoftClip(
                    position=ref_start,
                    length=length,
                    sequence=clip_seq,
                    side=side,
                )
            )

        # Update positions based on CIGAR operation
        # Operations that consume query: M, I, S, =, X (0, 1, 4, 7, 8)
        if op in (0, 1, 4, 7, 8):
            query_pos += length
        # Operations that consume reference: M, D, N, =, X (0, 2, 3, 7, 8)
        if op in (0, 2, 3, 7, 8):
            ref_start += length

    return clips


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
        ValueError: If region format is invalid or exceeds MAX_REGION_SIZE.
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

    region_size = end - start
    if region_size > MAX_REGION_SIZE:
        raise ValueError(
            f"Region size {region_size:,}bp exceeds maximum allowed {MAX_REGION_SIZE:,}bp. "
            f"Please request a smaller region."
        )

    return contig, start, end


def fetch_region(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    max_reads: int = 10000,
    min_mapq: int = 0,
    index_filename: str | None = None,
    min_vaf: float = 0.1,
    min_depth: int = 3,
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

    # Use context manager to ensure file handles are properly closed on exception
    with pysam.AlignmentFile(
        bam_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference_path,
        index_filename=index_filename,
    ) as samfile:
        # 1. Calculate coverage and base counts using pysam's optimized C engine
        # Returns tuple of arrays (A, C, G, T)
        # Use read_callback to respect min_mapq filter for consistency with read filtering
        cov_A, cov_C, cov_G, cov_T = samfile.count_coverage(
            contig,
            start,
            end,
            quality_threshold=0,
            read_callback=lambda r: (
                r.mapping_quality >= min_mapq
                and not r.is_unmapped
                and not r.is_secondary
                and not r.is_supplementary
            ),
        )

        # Vectorized sum using numpy for O(n) instead of Python loop
        # count_coverage returns numpy arrays, so we can add directly
        coverage = (np.array(cov_A) + np.array(cov_C) + np.array(cov_G) + np.array(cov_T)).tolist()

        reads: list[AlignedRead] = []

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
                            mismatches.append(
                                {
                                    "pos": rpos,
                                    "ref": ref_base.upper(),
                                    "alt": query_base,
                                }
                            )

            # Extract paired-end information
            mate_position = None
            mate_contig = None
            insert_size = None
            is_proper_pair = False
            is_read1 = False

            if read.is_paired:
                is_proper_pair = read.is_proper_pair
                is_read1 = read.is_read1
                insert_size = read.template_length if read.template_length != 0 else None

                if read.next_reference_id >= 0:
                    mate_position = read.next_reference_start
                    mate_contig = samfile.get_reference_name(read.next_reference_id)

            # Extract soft clips from CIGAR
            soft_clips = extract_soft_clips(read)

            reads.append(
                AlignedRead(
                    name=read.query_name or "",
                    sequence=read.query_sequence or "",
                    qualities=list(read.query_qualities or []),
                    cigar=read.cigarstring or "",
                    position=read.reference_start if read.reference_start is not None else 0,
                    end_position=(
                        read.reference_end
                        if read.reference_end is not None
                        else (read.reference_start or 0)
                    ),
                    mapping_quality=read.mapping_quality,
                    is_reverse=read.is_reverse,
                    mismatches=mismatches,
                    soft_clips=soft_clips,
                    mate_position=mate_position,
                    mate_contig=mate_contig,
                    insert_size=insert_size,
                    is_proper_pair=is_proper_pair,
                    is_read1=is_read1,
                    is_paired=read.is_paired,
                )
            )

    # Context manager ensures samfile is closed before we continue

    variants = detect_variants(
        (cov_A, cov_C, cov_G, cov_T), ref_seq, contig, start, min_vaf, min_depth
    )

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
    coverage_counts: tuple[Sequence[int], Sequence[int], Sequence[int], Sequence[int]],
    ref_seq: str | None,
    contig: str,
    start: int,
    min_vaf: float = 0.1,
    min_depth: int = 10,
) -> list[dict]:
    """
    Detect variants from coverage counts (A, C, G, T).

    Args:
        coverage_counts: Tuple of 4 arrays (A, C, G, T) from count_coverage.
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

    cov_A, cov_C, cov_G, cov_T = coverage_counts

    # Iterate over positions
    for i in range(len(cov_A)):
        # Total depth at this position (excluding N/Del)
        counts = {"A": cov_A[i], "C": cov_C[i], "G": cov_G[i], "T": cov_T[i]}
        total = sum(counts.values())

        if total < min_depth:
            continue

        if i >= len(ref_seq):
            break

        ref_base = ref_seq[i].upper()

        # Check for alternatives
        for base, count in counts.items():
            if base == ref_base:
                continue

            if count == 0:
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

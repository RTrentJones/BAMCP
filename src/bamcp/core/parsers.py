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


def read_passes_filters(read: pysam.AlignedSegment, min_mapq: int = 0) -> bool:
    """Return whether a read should contribute to BAMCP coverage/read evidence."""
    return (
        read.mapping_quality >= min_mapq
        and not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and not read.is_duplicate
        and not read.is_qcfail
    )


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


def _open_alignment(
    bam_path: str,
    reference_path: str | None = None,
    index_filename: str | None = None,
) -> pysam.AlignmentFile:
    """Open a BAM/CRAM with the correct mode and shared options."""
    mode = "rc" if bam_path.endswith(".cram") else "rb"
    return pysam.AlignmentFile(
        bam_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference_path,
        index_filename=index_filename,
    )


def fetch_reference_sequence(reference_path: str, region: str) -> str:
    """Fetch the reference sequence for a region (0-based half-open)."""
    contig, start, end = parse_region(region)
    with pysam.FastaFile(reference_path) as fasta:
        return fasta.fetch(contig, start, end)


def count_coverage_arrays(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
) -> tuple[str, int, int, tuple[Sequence[int], Sequence[int], Sequence[int], Sequence[int]]]:
    """Return base-count coverage arrays for a region without extracting reads."""
    contig, start, end = parse_region(region)
    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        counts = samfile.count_coverage(
            contig,
            start,
            end,
            quality_threshold=min_baseq,
            read_callback=lambda r: read_passes_filters(r, min_mapq),
        )
    return contig, start, end, counts


def fetch_coverage_only(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
) -> RegionData:
    """Fetch only coverage for a region, skipping read serialization."""
    contig, start, end, counts = count_coverage_arrays(
        bam_path, region, reference_path, min_mapq, min_baseq, index_filename
    )
    cov_A, cov_C, cov_G, cov_T = counts
    coverage = (np.array(cov_A) + np.array(cov_C) + np.array(cov_G) + np.array(cov_T)).tolist()
    return RegionData(contig=contig, start=start, end=end, reads=[], coverage=coverage, variants=[])


def fetch_candidate_variants_only(
    bam_path: str,
    region: str,
    reference_path: str,
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
    min_vaf: float = 0.1,
    min_depth: int = 3,
) -> RegionData:
    """Fetch coverage-derived candidate variants without serializing reads."""
    contig, start, end, counts = count_coverage_arrays(
        bam_path, region, reference_path, min_mapq, min_baseq, index_filename
    )
    cov_A, cov_C, cov_G, cov_T = counts
    coverage = (np.array(cov_A) + np.array(cov_C) + np.array(cov_G) + np.array(cov_T)).tolist()
    with pysam.FastaFile(reference_path) as fasta:
        ref_seq = fasta.fetch(contig, start, end)
    indels, gap_depth = _detect_indels_if_small(
        bam_path,
        region,
        reference_path,
        ref_seq,
        min_mapq,
        min_baseq,
        index_filename,
        min_vaf,
        min_depth,
    )
    variants = detect_variants(counts, ref_seq, contig, start, min_vaf, min_depth, gap_depth)
    variants.extend(indels)
    return RegionData(
        contig=contig,
        start=start,
        end=end,
        reads=[],
        coverage=coverage,
        variants=variants,
        reference_sequence=ref_seq,
    )


def load_vcf_variants(vcf_path: str, region: str) -> list[dict]:
    """Load VCF/BCF records for a region as internal 0-based candidate variants."""
    contig, start, end = parse_region(region)
    variants: list[dict] = []
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.fetch(contig, start, end):
            pos0 = int(record.pos) - 1
            if pos0 < start or pos0 >= end:
                continue
            depth = int(record.info.get("DP", 0) or 0)
            af_value = record.info.get("AF", 0.0)
            af_values = af_value if isinstance(af_value, tuple) else (af_value,)
            sv_type = record.info.get("SVTYPE")
            sv_end = record.info.get("END")
            sv_len = record.info.get("SVLEN")
            is_structural = bool(sv_type) or any(
                str(alt).startswith("<") for alt in record.alts or []
            )
            samples = {
                sample_name: {
                    key: list(value) if isinstance(value, tuple) else value
                    for key, value in sample_data.items()
                }
                for sample_name, sample_data in record.samples.items()
            }
            for idx, alt in enumerate(record.alts or []):
                af = float(af_values[idx]) if idx < len(af_values) else 0.0
                alt_len = len(str(alt))
                ref_len = len(str(record.ref))
                variant_kind = (
                    "sv" if is_structural else ("snv" if ref_len == 1 and alt_len == 1 else "indel")
                )
                variants.append(
                    {
                        "contig": record.contig,
                        "position": pos0,
                        "ref": record.ref,
                        "alt": alt,
                        "variant_kind": variant_kind,
                        "vaf": round(af, 3),
                        "depth": depth,
                        "alt_count": round(af * depth) if depth else 0,
                        "source": "vcf",
                        "samples": samples,
                        "sample_names": list(samples.keys()),
                        "sv_type": sv_type,
                        "sv_end": sv_end,
                        "sv_len": sv_len,
                    }
                )
    return variants


def annotate_vcf_snv_support(
    bam_path: str,
    region: str,
    reference_path: str | None,
    variants: list[dict],
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
) -> list[dict]:
    """Attach truncation-free read-level support + per-read metrics to VCF SNVs (in place).

    Pileups each VCF SNV site — htslib sees every read there, so this is immune to
    the ``max_reads`` cap that limits the serialized read list — and, for reads
    carrying the alt allele, collects strand, base quality, position-in-read, and
    MAPQ. Adds ``strand_forward``/``strand_reverse``/``read_depth``/
    ``read_support_vaf`` plus raw ``_alt_qualities``/``_alt_positions``/``_alt_mapqs``
    arrays that :func:`bamcp.analysis.evidence._vcf_evidence` turns into the
    quality/position/MAPQ histograms and mean quality. Indels/SVs are untouched
    (single-base support can't corroborate them).
    """
    snvs = [v for v in variants if v.get("source") == "vcf" and len(v["ref"]) == 1 == len(v["alt"])]
    if not snvs:
        return variants

    contig, start, end = parse_region(region)
    # A position may host several alt alleles, so index the targets by position.
    by_pos: dict[int, list[dict]] = {}
    for v in snvs:
        by_pos.setdefault(v["position"], []).append(v)

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        for column in samfile.pileup(
            contig, start, end, truncate=True, min_base_quality=min_baseq, stepper="nofilter"
        ):
            targets = by_pos.get(column.reference_pos)
            if not targets:
                continue

            depth = 0
            # Per-alt accumulators for alt-supporting reads at this column.
            acc: dict[str, dict] = {
                v["alt"].upper(): {"fwd": 0, "rev": 0, "quals": [], "positions": [], "mapqs": []}
                for v in targets
            }
            for pread in column.pileups:
                aln = pread.alignment
                if not read_passes_filters(aln, min_mapq) or pread.is_del or pread.is_refskip:
                    continue
                qpos = pread.query_position
                seq = aln.query_sequence
                if qpos is None or seq is None:
                    continue
                depth += 1
                support = acc.get(seq[qpos].upper())
                if support is None:
                    continue
                if aln.is_reverse:
                    support["rev"] += 1
                else:
                    support["fwd"] += 1
                quals = aln.query_qualities
                if quals is not None and qpos < len(quals):
                    support["quals"].append(int(quals[qpos]))
                    support["positions"].append(min(qpos, len(quals) - qpos - 1))
                support["mapqs"].append(aln.mapping_quality)

            for v in targets:
                support = acc[v["alt"].upper()]
                alt = support["fwd"] + support["rev"]
                v["strand_forward"] = support["fwd"]
                v["strand_reverse"] = support["rev"]
                v["read_depth"] = depth
                v["read_support_vaf"] = round(alt / depth, 3) if depth else 0.0
                v["_alt_qualities"] = support["quals"]
                v["_alt_positions"] = support["positions"]
                v["_alt_mapqs"] = support["mapqs"]
                # Fill any DP/AF/alt_count the VCF omitted from the counted support,
                # keeping the displayed depth/VAF/alt-reads mutually consistent.
                if not v.get("depth"):
                    v["depth"] = depth
                if not v.get("vaf"):
                    v["vaf"] = v["read_support_vaf"]
                if not v.get("alt_count"):
                    v["alt_count"] = alt

    return variants


def fetch_region(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    max_reads: int = 10000,
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
    min_vaf: float = 0.1,
    min_depth: int = 3,
    detect: bool = True,
) -> RegionData:
    """
    Fetch reads from a BAM/CRAM file for a given region.

    Args:
        bam_path: Path to BAM/CRAM file (local or remote).
        region: Genomic region (e.g., "chr1:1000-2000").
        reference_path: Path to reference FASTA (required for CRAM).
        max_reads: Maximum reads to return (stop after this many).
        min_mapq: Minimum mapping quality filter.
        min_baseq: Minimum base quality for coverage/candidate variant counts.
        index_filename: Path to index file (.bai/.crai). Used to redirect
            remote BAM index downloads to a cache directory.

    Returns:
        RegionData with reads, coverage, and detected variants.
    """
    contig, start, end = parse_region(region)

    # Use context manager to ensure file handles are properly closed on exception
    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        # 1. Calculate coverage and base counts using pysam's optimized C engine
        # Returns tuple of arrays (A, C, G, T)
        # Use read_callback to respect min_mapq filter for consistency with read filtering
        cov_A, cov_C, cov_G, cov_T = samfile.count_coverage(
            contig,
            start,
            end,
            quality_threshold=min_baseq,
            read_callback=lambda r: read_passes_filters(r, min_mapq),
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
            if not read_passes_filters(read, min_mapq):
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

    # Candidate detection can be skipped (e.g. VCF-primary, where local candidates
    # would be discarded) to avoid scanning a broad region for nothing.
    variants: list[dict] = []
    if detect:
        indels, gap_depth = _detect_indels_if_small(
            bam_path,
            region,
            reference_path,
            ref_seq,
            min_mapq,
            min_baseq,
            index_filename,
            min_vaf,
            min_depth,
        )
        variants = detect_variants(
            (cov_A, cov_C, cov_G, cov_T), ref_seq, contig, start, min_vaf, min_depth, gap_depth
        )
        variants.extend(indels)

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
    gap_depth: Sequence[int] | None = None,
) -> list[dict]:
    """
    Detect candidate SNVs from coverage counts (A, C, G, T).

    Args:
        coverage_counts: Tuple of 4 arrays (A, C, G, T) from count_coverage.
        ref_seq: Reference sequence for the region.
        contig: Contig name.
        start: Start position of the region.
        min_vaf: Minimum variant allele frequency threshold.
        min_depth: Minimum read depth to consider.
        gap_depth: Optional per-position count of reads carrying a deletion or
            reference-skip at each column (from :func:`detect_indels`). count_coverage
            only tallies A/C/G/T, so at deletion-heavy loci the naive ``sum(A,C,G,T)``
            denominator understates depth and inflates VAF. Adding gap_depth restores
            the true spanning depth.

    Returns:
        List of candidate variant dicts (contig, position, ref, alt, vaf, depth,
        alt_count, variant_kind="snv", source="bamcp").
    """
    variants: list[dict] = []

    if not ref_seq:
        return variants

    cov_A, cov_C, cov_G, cov_T = coverage_counts

    # Iterate over positions
    for i in range(len(cov_A)):
        # True spanning depth: matched bases (A/C/G/T) plus reads deleted/skipped here.
        counts = {"A": cov_A[i], "C": cov_C[i], "G": cov_G[i], "T": cov_T[i]}
        gap = gap_depth[i] if gap_depth is not None and i < len(gap_depth) else 0
        total = sum(counts.values()) + gap

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
                        "variant_kind": "snv",
                        "source": "bamcp",
                    }
                )

    return variants


def detect_indels(
    bam_path: str,
    region: str,
    reference_path: str | None = None,
    ref_seq: str | None = None,
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
    min_vaf: float = 0.1,
    min_depth: int = 3,
) -> tuple[list[dict], list[int]]:
    """Detect candidate insertions/deletions via a pysam pileup pass.

    count_coverage (used for SNVs/coverage) only tallies A/C/G/T, so it can neither
    surface indels nor account for reads spanning a deletion. This pileup pass fills
    both gaps in one traversal.

    Returns a ``(indels, gap_depth)`` tuple where:
      - ``indels`` are candidate records (``variant_kind="indel"``, ``indel_type``
        of ``"ins"``/``"del"``) in VCF-style anchored representation (ref/alt share a
        leading anchor base), 0-based ``position`` like :func:`detect_variants`.
      - ``gap_depth`` is a per-position count of reads carrying a deletion or
        reference-skip at each column, fed back into :func:`detect_variants` to
        correct the SNV depth denominator.

    ``gap_depth`` always has length ``end - start``. Without ``ref_seq`` no alleles can
    be built, so ``indels`` is empty (SNV detection also no-ops without a reference).
    """
    contig, start, end = parse_region(region)
    length = max(end - start, 0)
    gap_depth = [0] * length
    indels: list[dict] = []

    if not ref_seq or length == 0:
        return indels, gap_depth

    ref_len = len(ref_seq)

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        # nofilter + explicit read_passes_filters mirrors count_coverage's read set
        # (no BAQ, no mate-overlap dedup) so indel depth stays consistent with SNV depth.
        for column in samfile.pileup(
            contig,
            start,
            end,
            truncate=True,
            min_base_quality=min_baseq,
            stepper="nofilter",
        ):
            idx = column.reference_pos - start
            if not (0 <= idx < length):
                continue

            matched = 0
            del_here = 0
            del_events: dict[int, int] = {}
            ins_events: dict[str, int] = {}

            for pread in column.pileups:
                aln = pread.alignment
                if not read_passes_filters(aln, min_mapq):
                    continue
                if pread.is_refskip:
                    del_here += 1
                    continue
                if pread.is_del:
                    del_here += 1
                    continue
                matched += 1
                # pread.indel is the indel BETWEEN this column and the next.
                if pread.indel > 0:
                    qpos = pread.query_position
                    seq = aln.query_sequence
                    if qpos is not None and seq is not None:
                        ins_seq = seq[qpos + 1 : qpos + 1 + pread.indel].upper()
                        if ins_seq:
                            ins_events[ins_seq] = ins_events.get(ins_seq, 0) + 1
                elif pread.indel < 0:
                    del_len = -pread.indel
                    del_events[del_len] = del_events.get(del_len, 0) + 1

            gap_depth[idx] = del_here
            # Guard against a reference shorter than the BAM contig (mismatched
            # reference): pileup can yield a column past ref_seq, and detect_variants
            # makes the same check with `i >= len(ref_seq)`.
            if idx >= ref_len:
                continue
            depth = matched + del_here  # reads spanning this column
            if depth < min_depth:
                continue

            anchor = ref_seq[idx].upper()

            for del_len, support in del_events.items():
                vaf = support / depth
                if vaf < min_vaf:
                    continue
                # VCF-style: ref = anchor + deleted bases, alt = anchor. If the
                # deletion runs past the fetched reference slice (anchor near the
                # region edge) we can't build the real deleted allele, so skip it
                # rather than emit a non-variant with ref == alt == anchor.
                del_end = idx + 1 + del_len
                if del_end > ref_len:
                    continue
                ref_allele = ref_seq[idx:del_end].upper()
                indels.append(
                    {
                        "contig": contig,
                        "position": column.reference_pos,
                        "ref": ref_allele,
                        "alt": anchor,
                        "variant_kind": "indel",
                        "indel_type": "del",
                        "vaf": round(vaf, 3),
                        "depth": depth,
                        "alt_count": support,
                        "source": "bamcp",
                    }
                )

            for ins_seq, support in ins_events.items():
                vaf = support / depth
                if vaf < min_vaf:
                    continue
                # VCF-style: ref = anchor, alt = anchor + inserted bases.
                indels.append(
                    {
                        "contig": contig,
                        "position": column.reference_pos,
                        "ref": anchor,
                        "alt": anchor + ins_seq,
                        "variant_kind": "indel",
                        "indel_type": "ins",
                        "vaf": round(vaf, 3),
                        "depth": depth,
                        "alt_count": support,
                        "source": "bamcp",
                    }
                )

    return indels, gap_depth


def _detect_indels_if_small(
    bam_path: str,
    region: str,
    reference_path: str | None,
    ref_seq: str | None,
    min_mapq: int,
    min_baseq: int,
    index_filename: str | None,
    min_vaf: float,
    min_depth: int,
) -> tuple[list[dict], list[int] | None]:
    """Run :func:`detect_indels` only when the region is small enough for a pileup.

    Returns ``(indels, gap_depth)``. For oversized regions returns ``([], None)`` so
    the SNV path keeps its fast count_coverage-only behavior and stays responsive.
    """
    from ..constants import INDEL_DETECTION_MAX_REGION

    _, start, end = parse_region(region)
    if (end - start) > INDEL_DETECTION_MAX_REGION:
        return [], None
    return detect_indels(
        bam_path,
        region,
        reference_path,
        ref_seq,
        min_mapq,
        min_baseq,
        index_filename,
        min_vaf,
        min_depth,
    )


def scan_variants_chunked(
    bam_path: str,
    contig: str,
    reference_path: str,
    chunk_size: int = 50_000,
    min_vaf: float = 0.1,
    min_depth: int = 3,
    min_mapq: int = 0,
    min_baseq: int = 0,
    max_region: int = 250_000_000,
    index_filename: str | None = None,
) -> list[dict]:
    """Scan an entire contig for variants using fast coverage-based detection.

    Uses pysam's count_coverage (C implementation) in chunks — no read
    extraction — so scanning a full human chromosome is feasible.

    Args:
        bam_path: Path to BAM/CRAM file.
        contig: Contig to scan (e.g. "chr13").
        reference_path: Path to reference FASTA (required).
        chunk_size: Window size per iteration.
        min_vaf: Minimum variant allele frequency.
        min_depth: Minimum read depth.
        min_mapq: Minimum mapping quality filter.
        min_baseq: Minimum base quality for coverage/candidate variant counts.
        max_region: Maximum bases to scan (safety limit).
        index_filename: Path to cached index file for remote BAMs.

    Returns:
        List of variant dicts sorted by VAF descending, capped at 500.
    """
    from ..constants import SCAN_VARIANTS_MAX_RESULTS

    all_variants: list[dict] = []

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        # Get contig length from header
        contig_lengths = dict(zip(samfile.references, samfile.lengths, strict=False))
        contig_length = contig_lengths.get(contig)
        if contig_length is None:
            return []

        scan_end = min(contig_length, max_region)

        with pysam.FastaFile(reference_path) as fasta:
            for chunk_start in range(0, scan_end, chunk_size):
                chunk_end = min(chunk_start + chunk_size, scan_end)

                ref_seq = fasta.fetch(contig, chunk_start, chunk_end)

                cov_A, cov_C, cov_G, cov_T = samfile.count_coverage(
                    contig,
                    chunk_start,
                    chunk_end,
                    quality_threshold=min_baseq,
                    read_callback=lambda r: read_passes_filters(r, min_mapq),
                )

                chunk_variants = detect_variants(
                    (cov_A, cov_C, cov_G, cov_T),
                    ref_seq,
                    contig,
                    chunk_start,
                    min_vaf,
                    min_depth,
                )
                all_variants.extend(chunk_variants)

    # Sort by VAF descending and cap results
    all_variants.sort(key=lambda v: v["vaf"], reverse=True)
    return all_variants[:SCAN_VARIANTS_MAX_RESULTS]

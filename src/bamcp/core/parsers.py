"""BAM/CRAM file parsing using pysam."""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pysam

from .validation import is_cram


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
    # True when the max_reads cap stopped read extraction early: ``reads`` is then a
    # prefix of the region's reads, so read-derived evidence is incomplete (coverage
    # and candidate VAF/depth come from uncapped passes and are unaffected).
    reads_truncated: bool = False


# Maximum region size to prevent DoS via unbounded memory allocation
MAX_REGION_SIZE = 1_000_000  # 1 Mbp

# pysam.pileup defaults to max_depth=8000, which silently truncates high-depth
# (amplicon/targeted) sites. The count_coverage path this replaces was uncapped, so
# lift the ceiling to effectively unbounded for the VCF-support pileup — a genomic
# region is already capped at MAX_REGION_SIZE, and htslib allocates per-column lazily.
PILEUP_MAX_DEPTH = 1_000_000_000

# The VCF-support pileup runs a Python per-read loop, so it must not walk the whole
# request window (up to MAX_REGION_SIZE) when only a few SNVs are present. SNV positions
# within this many bases are grouped into one pileup interval — small enough to skip the
# empty coverage between sparse sites, wide enough that clustered sites share a pass.
SNV_PILEUP_MERGE_GAP = 1000

# CIGAR operation codes
CIGAR_SOFT_CLIP = 4  # S operation

# Reads BAMCP never counts — for coverage, candidates, VCF support, or display. Echoed
# in every variant response's applied_filters so the exclusions are never silent.
EXCLUDED_READ_FLAGS = ("unmapped", "secondary", "supplementary", "duplicate", "qcfail")

# Overlapping mates policy: when both mates of a fragment cover a base, BAMCP counts
# BOTH (like ``samtools depth`` without ``-s``). count_coverage cannot de-duplicate
# overlaps, so every pileup pass sets ``ignore_overlaps=False`` explicitly to agree with
# it — pysam's default (True) silently zeroes one mate's base quality, which made depth
# and VAF differ between the coverage, indel and VCF-support passes.
OVERLAPPING_MATES_POLICY = "both_counted"

_ACGT = frozenset("ACGT")


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


def _pileup(samfile: pysam.AlignmentFile, contig: str, start: int, stop: int) -> Any:
    """The one pileup configuration every BAMCP pass uses.

    ``stepper="nofilter"`` + explicit :func:`read_passes_filters` mirrors count_coverage's
    read set; ``min_base_quality=0`` because pileup's own base-quality filter also drops
    deletion/ref-skip entries (which have no base) — callers threshold real bases
    themselves via :func:`_counted_base`; ``ignore_overlaps=False`` per
    :data:`OVERLAPPING_MATES_POLICY`; no BAQ; ``max_depth`` lifted off pysam's 8000.
    """
    return samfile.pileup(
        contig,
        start,
        stop,
        truncate=True,
        min_base_quality=0,
        stepper="nofilter",
        ignore_overlaps=False,
        compute_baq=False,
        max_depth=PILEUP_MAX_DEPTH,
    )


def _counted_base(pread: Any, min_baseq: int) -> str | None:
    """The upper-case base a pileup read contributes, or None if count_coverage wouldn't.

    count_coverage tallies only A/C/G/T whose stored quality meets the threshold, so
    under a positive ``min_baseq`` a no-QUAL ("*") read — or one whose quality array
    is too short — does NOT count either. Deletion/ref-skip entries carry no base.
    """
    if pread.is_del or pread.is_refskip:
        return None
    qpos = pread.query_position
    seq = pread.alignment.query_sequence
    if qpos is None or seq is None or qpos >= len(seq):
        return None
    base = seq[qpos].upper()
    if base not in _ACGT:
        return None
    if min_baseq:
        quals = pread.alignment.query_qualities
        if quals is None or qpos >= len(quals) or quals[qpos] < min_baseq:
            return None
    return str(base)


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
    Parse a 1-based, inclusive (samtools-style) region string into internal coordinates.

    The tool-boundary contract: ``chr1:1051-1051`` is the single base at 1-based
    position 1051, and ``chr1:1000-2000`` spans 1000..2000 inclusive (1001 bases) —
    exactly what ``samtools view``/IGV mean by the same string. Internally BAMCP works
    in pysam's 0-based half-open frame, so this is the ONE place the conversion
    happens: the returned ``start`` is ``int(start) - 1`` and ``end`` is unchanged
    (1-based inclusive end == 0-based exclusive end).

    Supports formats:
        - chr1:1000-2000
        - chr1:1,000-2,000
        - 1:1000-2000
        - HLA-A*01:01:01:01:100-200 (contig names containing ':' split on the last ':')

    The end is NOT clamped here (no header is available); the fetchers clamp it to
    the contig length once the file is open (see :func:`_clamp_end`).

    Returns:
        Tuple of (contig, start0, end) — 0-based half-open.

    Raises:
        ValueError: If region format is invalid or exceeds MAX_REGION_SIZE.
    """
    region = region.replace(",", "")
    try:
        contig, coords = region.rsplit(":", 1)
        start_str, end_str = coords.split("-")
        start1 = int(start_str)
        end = int(end_str)
    except (ValueError, AttributeError) as e:
        raise ValueError(
            f"Invalid region format: '{region}'. Expected 1-based inclusive "
            f"'chr1:1000-2000' (samtools-style)"
        ) from e
    if not contig:
        raise ValueError(f"Invalid region format: '{region}'. Missing contig name")

    if start1 < 1:
        raise ValueError(
            f"Start position must be >= 1 (regions are 1-based, inclusive), got {start1}"
        )
    if end < start1:
        raise ValueError(f"End position ({end}) must be >= start ({start1})")

    start = start1 - 1
    region_size = end - start
    if region_size > MAX_REGION_SIZE:
        raise ValueError(
            f"Region size {region_size:,}bp exceeds maximum allowed {MAX_REGION_SIZE:,}bp. "
            f"Please request a smaller region."
        )

    return contig, start, end


def format_region(contig: str, start0: int, end: int) -> str:
    """Render internal 0-based half-open ``[start0, end)`` as a 1-based inclusive region.

    The inverse of :func:`parse_region`: every internal caller that builds a region
    string from 0-based numbers must go through here so the string it emits means
    the same bases at the tool boundary.
    """
    return f"{contig}:{start0 + 1}-{end}"


def contig_not_found_message(contig: str, names: Iterable[str]) -> str:
    """A caller-safe "contig not found" message that names the other-style spelling.

    A BAM using ``13`` rejects ``chr13`` (and vice versa) — the commonest wrong-contig error,
    usually from coordinates copied out of a differently-styled source. Built only from the
    caller's contig and the file's own header names (never operator config).
    """
    from .reference import map_contig_to_style

    msg = f"Contig '{contig}' not found in the alignment header"
    alt = map_contig_to_style(contig, frozenset(names))
    if alt is not None and alt != contig:
        msg += f"; this file names it '{alt}' — retry with '{alt}'"
    else:
        msg += "; check contig names with list_contigs"
    return msg


def _clamp_end(lengths: dict[str, int], contig: str, start: int, end: int) -> int:
    """Clamp ``end`` to the contig length; raise a clear error for unknown contigs/starts."""
    length = lengths.get(contig)
    if length is None:
        raise ValueError(contig_not_found_message(contig, lengths))
    if start >= length:
        raise ValueError(
            f"Region start {start + 1} is past the end of {contig} (length {length:,})"
        )
    return min(end, length)


def _header_lengths(samfile: pysam.AlignmentFile) -> dict[str, int]:
    return dict(zip(samfile.references, samfile.lengths, strict=False))


def _open_alignment(
    bam_path: str,
    reference_path: str | None = None,
    index_filename: str | None = None,
) -> pysam.AlignmentFile:
    """Open a BAM/CRAM with the correct mode and shared options."""
    mode = "rc" if is_cram(bam_path) else "rb"
    return pysam.AlignmentFile(
        bam_path,
        mode,  # type: ignore[arg-type]
        reference_filename=reference_path,
        index_filename=index_filename,
    )


def fetch_reference_sequence(reference_path: str, region: str) -> str:
    """Fetch the reference sequence for a 1-based inclusive region string.

    Returns the bases at internal 0-based half-open ``[start, end)``; pysam truncates
    a slice that runs past the contig end, matching the alignment-side clamp.
    """
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
        end = _clamp_end(_header_lengths(samfile), contig, start, end)
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


def _per_alt_value(value: Any, idx: int) -> Any:
    """Select the ``idx``-th element of a per-ALT INFO value.

    Number=A INFO fields (e.g. ``SVLEN``) carry one value per ALT as a tuple; scalar
    (Number=1 / broadcast) values apply to every ALT. Returns ``None`` when a per-ALT
    tuple is shorter than the ALT it's asked for, so a mismatched record can't report
    another ALT's value.
    """
    if isinstance(value, (list, tuple)):
        return value[idx] if idx < len(value) else None
    return value


def _info_value(record: Any, key: str, default: Any = None) -> Any:
    """Read an INFO field without tripping pysam's undeclared-key error.

    ``record.info.get(key)`` raises ``ValueError("Invalid header")`` when the VCF
    header does not declare ``key`` — the norm for bcftools/GATK/DeepVariant small-
    variant VCFs, which never declare SVTYPE/END/SVLEN (and often not AF/DP).
    """
    if key not in record.header.info:
        return default
    return record.info.get(key, default)


_NUCLEOTIDES = frozenset("ACGTNacgtn")


def _trim_allele(pos0: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Minimally represent one ALT against REF (multi-allelic records are padded).

    A multi-allelic record shares one REF across ALTs (``REF=ACG ALT=ATG,A``), so
    each allele is trimmed on its own: shared suffix first, then shared prefix,
    always leaving at least one base on each side (VCF anchor convention). The
    result matches BAMCP's own candidate representation, so the ``auto`` overlay
    de-duplicates a VCF call against the same local candidate.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1].upper() == alt[-1].upper():
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0].upper() == alt[0].upper():
        ref, alt, pos0 = ref[1:], alt[1:], pos0 + 1
    return pos0, ref, alt


def _classify_allele(ref: str, alt: str, structural: bool) -> str:
    """Variant kind of a (trimmed) ref/alt pair: snv, mnp, indel, or sv."""
    if structural:
        return "sv"
    if len(ref) == len(alt):
        return "snv" if len(ref) == 1 else "mnp"
    return "indel"


def _sample_support(record: Any, idx: int) -> tuple[int, int] | None:
    """``(depth, alt_count)`` for ALT ``idx`` from a single-sample FORMAT AD/DP.

    The sample's own allele depths are what VAF means for the BAM being viewed;
    INFO AF/DP are cohort/site-level and only a fallback. Used only when the VCF has
    exactly one sample (with several, which one matches the BAM is unknowable) and
    AD carries a value for this ALT.
    """
    if len(record.samples) != 1 or "AD" not in record.format:
        return None
    sample = record.samples[0]
    ad = sample.get("AD")
    if not isinstance(ad, (list, tuple)) or idx + 1 >= len(ad) or ad[idx + 1] is None:
        return None
    alt_count = int(ad[idx + 1])
    dp = sample.get("DP") if "DP" in record.format else None
    depth = int(dp) if dp is not None else sum(int(x) for x in ad if x is not None)
    return depth, alt_count


def load_vcf_variants(vcf_path: str, region: str) -> list[dict]:
    """Load VCF/BCF records for a region as internal 0-based candidate variants.

    One variant per ALT allele, each trimmed to its minimal representation and
    classified as snv / mnp / indel / sv. Spanning-deletion ``*`` alleles are not
    variants at this site (the upstream deletion record carries them) and are
    skipped. Depth/VAF come from the single sample's FORMAT AD/DP when present,
    else from INFO AF/DP; a missing (``.``) value leaves VAF 0 for the read-support
    annotation to fill.

    Raises:
        ValueError: For a VCF that cannot be region-queried (plain/unindexed .vcf).
    """
    contig, start, end = parse_region(region)
    variants: list[dict] = []
    with pysam.VariantFile(vcf_path) as vcf:
        try:
            records = vcf.fetch(contig, start, end)
        except ValueError as e:
            raise ValueError(
                "vcf_path must be bgzip-compressed and indexed (.tbi/.csi) or an indexed "
                "BCF to be queried by region; run `bgzip calls.vcf && tabix -p vcf "
                f"calls.vcf.gz` ({e})"
            ) from e
        for record in records:
            pos0 = int(record.pos) - 1
            if pos0 < start or pos0 >= end:
                continue
            info_dp = _info_value(record, "DP", 0) or 0
            af_value = _info_value(record, "AF", 0.0)
            af_values = af_value if isinstance(af_value, tuple) else (af_value,)
            sv_type = _info_value(record, "SVTYPE")
            # pysam reserves INFO/END: ``record.info.get("END")`` is always None for a
            # real file, and the value surfaces as ``record.stop`` instead. Only report
            # it when the record actually carries an END (span differs from REF's).
            sv_end = (
                int(record.stop)
                if "END" in record.header.info and record.rlen != len(record.ref or "")
                else None
            )
            sv_len = _info_value(record, "SVLEN")
            samples = {
                sample_name: {
                    key: list(value) if isinstance(value, tuple) else value
                    for key, value in sample_data.items()
                }
                for sample_name, sample_data in record.samples.items()
            }
            for idx, alt in enumerate(record.alts or []):
                alt = str(alt)
                if alt == "*":
                    continue
                ref = str(record.ref)
                is_structural = bool(sv_type) or not (
                    set(alt) <= _NUCLEOTIDES and set(ref) <= _NUCLEOTIDES
                )
                var_pos0, var_ref, var_alt = (
                    (pos0, ref, alt) if is_structural else _trim_allele(pos0, ref, alt)
                )

                support = _sample_support(record, idx)
                if support is not None:
                    depth, alt_count = support
                    vaf = alt_count / depth if depth else 0.0
                else:
                    depth = int(info_dp)
                    raw_af = af_values[idx] if idx < len(af_values) else None
                    # AF=.,0.3 arrives as (None, 0.3): a missing value is "unknown",
                    # not a crash.
                    vaf = float(raw_af) if raw_af is not None else 0.0
                    alt_count = round(vaf * depth) if depth else 0
                variants.append(
                    {
                        "contig": record.contig,
                        "position": var_pos0,
                        "ref": var_ref,
                        "alt": var_alt,
                        "variant_kind": _classify_allele(var_ref, var_alt, is_structural),
                        "vaf": round(vaf, 3),
                        "depth": depth,
                        "alt_count": alt_count,
                        "source": "vcf",
                        "samples": samples,
                        "sample_names": list(samples.keys()),
                        "sv_type": sv_type,
                        "sv_end": sv_end,
                        # SVLEN is Number=A (one value per ALT); pick this ALT's.
                        "sv_len": _per_alt_value(sv_len, idx),
                    }
                )
    return variants


def _merge_snv_intervals(positions: list[int], gap: int) -> list[tuple[int, int]]:
    """Merge SNV positions into half-open ``[start, stop)`` pileup intervals.

    A position within ``gap`` bases of the current interval's end is folded into it, so
    clustered SNVs share one pileup pass while sparse SNVs stay narrow — keeping the
    per-read pileup loop off the unrelated coverage between distant sites.
    """
    intervals: list[tuple[int, int]] = []
    for p in sorted(positions):
        if intervals and p - intervals[-1][1] <= gap:
            intervals[-1] = (intervals[-1][0], p + 1)
        else:
            intervals.append((p, p + 1))
    return intervals


def annotate_vcf_snv_support(
    bam_path: str,
    region: str,
    reference_path: str | None,
    variants: list[dict],
    min_mapq: int = 0,
    min_baseq: int = 0,
    index_filename: str | None = None,
    sources: tuple[str, ...] = ("vcf",),
) -> list[dict]:
    """Attach truncation-free read-level support + per-read metrics to VCF SNVs (in place).

    ``sources`` selects which variants are annotated; :func:`fetch_region` also passes
    ``("bamcp",)`` to score its own SNV candidates from the pileup when the
    ``max_reads`` cap truncated the read list their evidence would otherwise use.

    Pileups each VCF SNV site — htslib sees every read there, so this is immune to
    the ``max_reads`` cap that limits the serialized read list — and, for reads
    carrying the alt allele, collects strand, base quality, position-in-read, and
    MAPQ. Adds ``strand_forward``/``strand_reverse``/``read_depth``/
    ``read_support_vaf`` plus raw ``_alt_qualities``/``_alt_positions``/``_alt_mapqs``
    arrays that :func:`bamcp.analysis.evidence._vcf_evidence` turns into the
    quality/position/MAPQ histograms and mean quality. Indels/SVs are untouched
    (single-base support can't corroborate them).
    """
    snvs = [
        v for v in variants if v.get("source") in sources and len(v["ref"]) == 1 == len(v["alt"])
    ]
    if not snvs:
        return variants

    contig = parse_region(region)[0]
    # A position may host several alt alleles, so index the targets by position.
    by_pos: dict[int, list[dict]] = {}
    for v in snvs:
        by_pos.setdefault(v["position"], []).append(v)
        # Zero-init every in-region SNV up front. pileup() emits no column for a
        # position with no overlapping reads, so covered columns overwrite these while
        # uncovered SNVs keep explicit zeros — letting callers distinguish true zero
        # BAM support from an unannotated variant (as the old count_coverage loop did).
        v["strand_forward"] = 0
        v["strand_reverse"] = 0
        v["read_depth"] = 0
        v["read_support_vaf"] = 0.0
        v["_alt_qualities"] = []
        v["_alt_positions"] = []
        v["_alt_mapqs"] = []

    # Pileup only short intervals around the SNVs — not the whole request window, which
    # (up to MAX_REGION_SIZE) would make this Python per-read loop walk unrelated coverage.
    intervals = _merge_snv_intervals(list(by_pos), SNV_PILEUP_MERGE_GAP)

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        # See _pileup: no pileup-side base filtering, both overlapping mates counted,
        # and max_depth lifted so high-depth sites aren't silently truncated.
        columns = (
            column
            for iv_start, iv_stop in intervals
            for column in _pileup(samfile, contig, iv_start, iv_stop)
        )
        for column in columns:
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
                if not read_passes_filters(aln, min_mapq):
                    continue
                # Honor min_baseq (and A/C/G/T-only) for both depth and support so
                # VCF-primary VAF/confidence agrees with count_coverage used elsewhere.
                base = _counted_base(pread, min_baseq)
                if base is None:
                    continue
                qpos = pread.query_position
                quals = aln.query_qualities
                depth += 1
                support = acc.get(base)
                if support is None:
                    continue
                if aln.is_reverse:
                    support["rev"] += 1
                else:
                    support["fwd"] += 1
                if quals is not None and qpos < len(quals):
                    support["quals"].append(int(quals[qpos]))
                    support["positions"].append(min(qpos, len(quals) - qpos - 1))
                support["mapqs"].append(aln.mapping_quality)

            for v in targets:
                support = acc[v["alt"].upper()]
                v["strand_forward"] = support["fwd"]
                v["strand_reverse"] = support["rev"]
                v["read_depth"] = depth
                v["read_support_vaf"] = (
                    round((support["fwd"] + support["rev"]) / depth, 3) if depth else 0.0
                )
                v["_alt_qualities"] = support["quals"]
                v["_alt_positions"] = support["positions"]
                v["_alt_mapqs"] = support["mapqs"]

    # Fill any DP/AF/alt_count the VCF omitted from the counted support, for every SNV
    # (covered or not), keeping the displayed depth/VAF/alt-reads mutually consistent.
    for v in snvs:
        if not v.get("depth"):
            v["depth"] = v["read_depth"]
        if not v.get("vaf"):
            v["vaf"] = v["read_support_vaf"]
        if not v.get("alt_count"):
            v["alt_count"] = v["strand_forward"] + v["strand_reverse"]

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
        end = _clamp_end(_header_lengths(samfile), contig, start, end)
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
        reads_truncated = False
        for read in samfile.fetch(contig, start, end):
            if not read_passes_filters(read, min_mapq):
                continue

            read_count += 1
            if read_count > max_reads:
                reads_truncated = True
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
        if reads_truncated:
            # The candidates' VAF/depth come from uncapped count_coverage, but their
            # strand/quality evidence would come from the capped read list — silently
            # zero for any site past the cap. Pileup the SNV sites instead (htslib sees
            # every read); serialization scores them from this support.
            annotate_vcf_snv_support(
                bam_path,
                region,
                reference_path,
                variants,
                min_mapq,
                min_baseq,
                index_filename,
                sources=("bamcp",),
            )

    return RegionData(
        contig=contig,
        start=start,
        end=end,
        reads=reads,
        coverage=coverage,
        variants=variants,
        reference_sequence=ref_seq,
        reads_truncated=reads_truncated,
    )


def detect_variants(
    coverage_counts: tuple[Sequence[int], Sequence[int], Sequence[int], Sequence[int]],
    ref_seq: str | None,
    contig: str,
    start: int,
    min_vaf: float = 0.1,
    min_depth: int = 10,
    gap_depth: Sequence[int] | None = None,
    min_alt_count: int = 1,
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
        gap_depth: Optional per-position count of reads carrying a deletion at each
            column (from :func:`detect_indels`; spliced ref-skips are excluded).
            count_coverage only tallies A/C/G/T, so at deletion-heavy loci the naive
            ``sum(A,C,G,T)`` denominator understates depth and inflates VAF. Adding
            gap_depth restores the true spanning depth.
        min_alt_count: Minimum alt-supporting reads (1 = no floor). The tool layer applies
            the configured floor itself; scan_contig passes it here so its result cap
            counts only candidates that survive it.

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
        # True spanning depth: matched bases (A/C/G/T) plus reads deleted here.
        counts = {"A": cov_A[i], "C": cov_C[i], "G": cov_G[i], "T": cov_T[i]}
        gap = gap_depth[i] if gap_depth is not None and i < len(gap_depth) else 0
        total = sum(counts.values()) + gap

        if total < min_depth:
            continue

        if i >= len(ref_seq):
            break

        ref_base = ref_seq[i].upper()
        # An N/IUPAC reference base has no "reference allele" to differ from: every
        # covered read would otherwise read as a 100%-VAF "variant".
        if ref_base not in _ACGT:
            continue

        # Check for alternatives
        for base, count in counts.items():
            if base == ref_base:
                continue

            if count == 0 or count < min_alt_count:
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
      - ``gap_depth`` is a per-position count of reads carrying a deletion at each
        column (a spliced ``N`` skip does not span the base and is not counted),
        fed back into :func:`detect_variants` to correct the SNV depth denominator.

    ``gap_depth`` always has length ``end - start`` (``end`` clamped to the contig
    length, as for the coverage arrays). Without ``ref_seq`` no alleles can
    be built, so ``indels`` is empty (SNV detection also no-ops without a reference).
    """
    contig, start, end = parse_region(region)
    indels: list[dict] = []

    if not ref_seq:
        return indels, [0] * max(end - start, 0)

    ref_len = len(ref_seq)

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        # Same clamp as count_coverage_arrays/fetch_region so gap_depth stays aligned
        # with the coverage arrays it corrects.
        end = _clamp_end(_header_lengths(samfile), contig, start, end)
        length = max(end - start, 0)
        gap_depth = [0] * length
        # _pileup mirrors count_coverage's read set (no BAQ, both overlapping mates
        # counted) so indel depth stays consistent with SNV depth.
        for column in _pileup(samfile, contig, start, end):
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
                    # A spliced (N) read skips this base entirely — it neither covers it
                    # nor carries a deletion, so it must not enter any depth denominator.
                    continue
                if pread.is_del:
                    # A deletion has no base quality; count it regardless of min_baseq
                    # (pileup's own min_base_quality used to drop it here).
                    del_here += 1
                    continue
                # Matched bases follow count_coverage's rule (A/C/G/T, qual >= min_baseq).
                if _counted_base(pread, min_baseq) is None:
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
            # As for SNVs: an N/IUPAC anchor can't form a real VCF-style allele.
            if anchor not in _ACGT:
                continue

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
                if not set(ref_allele) <= _ACGT:
                    continue  # deleted bases include N/IUPAC — no real ref allele
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


@dataclass
class ScanResult:
    """A contig scan's candidates plus what was (and wasn't) scanned.

    ``scanned_start``/``scanned_end`` are internal 0-based half-open; the tool layer
    renders them 1-based. ``total_candidates`` counts candidates before the
    ``SCAN_VARIANTS_MAX_RESULTS`` cap, so ``results_capped`` is self-describing.
    """

    contig: str
    contig_length: int
    scanned_start: int
    scanned_end: int
    variants: list[dict]
    total_candidates: int
    results_capped: bool
    truncated_by_max_region: bool


def _deletion_depth_at(
    samfile: pysam.AlignmentFile, contig: str, positions: list[int], min_mapq: int
) -> dict[int, int]:
    """Count reads carrying a deletion at each of ``positions`` (spliced skips excluded).

    The scan's SNV candidates come from count_coverage, whose A/C/G/T-only denominator
    understates depth inside deletions. Rather than pileup a whole chromosome, pileup
    only short intervals around the candidates and add the deletion depth back — the
    same gap_depth correction get_variants applies via detect_indels.
    """
    wanted = set(positions)
    counts: dict[int, int] = {}
    for iv_start, iv_stop in _merge_snv_intervals(positions, SNV_PILEUP_MERGE_GAP):
        for column in _pileup(samfile, contig, iv_start, iv_stop):
            pos = column.reference_pos
            if pos not in wanted:
                continue
            counts[pos] = sum(
                1
                for pread in column.pileups
                if pread.is_del and read_passes_filters(pread.alignment, min_mapq)
            )
    return counts


def scan_contig(
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
    min_alt_count: int = 1,
) -> ScanResult:
    """Scan a contig for SNV candidates using fast coverage-based detection.

    Uses pysam's count_coverage (C implementation) in chunks — no read
    extraction — so scanning a full human chromosome is feasible. SNV-only: indels
    need a per-read pileup that is too slow at chromosome scale. Candidate VAF/depth
    include reads deleted at the site (as get_variants does), and N/IUPAC reference
    bases are skipped.

    Raises:
        ValueError: If ``contig`` is not in the alignment header.
    """
    from ..constants import SCAN_VARIANTS_MAX_RESULTS

    all_variants: list[dict] = []

    with _open_alignment(bam_path, reference_path, index_filename) as samfile:
        header_lengths = _header_lengths(samfile)
        contig_length = header_lengths.get(contig)
        if contig_length is None:
            raise ValueError(contig_not_found_message(contig, header_lengths))

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
                    min_alt_count=min_alt_count,
                )
                if chunk_variants:
                    gaps = _deletion_depth_at(
                        samfile, contig, sorted({v["position"] for v in chunk_variants}), min_mapq
                    )
                    for v in chunk_variants:
                        gap = gaps.get(v["position"], 0)
                        if gap:
                            v["depth"] += gap
                            v["vaf"] = round(v["alt_count"] / v["depth"], 3)
                    # Adding deletion depth only lowers VAF: re-apply the threshold.
                    chunk_variants = [v for v in chunk_variants if v["vaf"] >= min_vaf]
                all_variants.extend(chunk_variants)

    # Sort by VAF descending and cap results
    all_variants.sort(key=lambda v: v["vaf"], reverse=True)
    return ScanResult(
        contig=contig,
        contig_length=contig_length,
        scanned_start=0,
        scanned_end=scan_end,
        variants=all_variants[:SCAN_VARIANTS_MAX_RESULTS],
        total_candidates=len(all_variants),
        results_capped=len(all_variants) > SCAN_VARIANTS_MAX_RESULTS,
        truncated_by_max_region=scan_end < contig_length,
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
    """Candidates from :func:`scan_contig`, sorted by VAF descending, capped at 500."""
    return scan_contig(
        bam_path,
        contig,
        reference_path,
        chunk_size,
        min_vaf,
        min_depth,
        min_mapq,
        min_baseq,
        max_region,
        index_filename,
    ).variants

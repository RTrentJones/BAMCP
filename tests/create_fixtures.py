"""Script to create test BAM and FASTA fixtures for testing."""

import os

import pysam

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


def create_reference_fasta():
    """Create a small reference FASTA file."""
    os.makedirs(FIXTURES_DIR, exist_ok=True)
    ref_path = os.path.join(FIXTURES_DIR, "ref.fa")

    # 1000bp reference sequence for chr1
    # Use a known sequence so we can create reads with mismatches
    seq = "ACGTACGTAC" * 100  # 1000bp

    with open(ref_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")

    # Index the FASTA
    pysam.faidx(ref_path)
    return ref_path


def create_small_bam(ref_path: str):
    """Create a small BAM file with test reads."""
    bam_path = os.path.join(FIXTURES_DIR, "small.bam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 500}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Read 1: Forward strand, perfect match at position 100
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "ACGTACGTAC" * 5  # 50bp
        a.flag = 0  # forward strand
        a.reference_id = 0  # chr1
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [(0, 50)]  # 50M
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # Read 2: Reverse strand at position 120
        a = pysam.AlignedSegment()
        a.query_name = "read2"
        a.query_sequence = "ACGTACGTAC" * 5
        a.flag = 16  # reverse strand
        a.reference_id = 0
        a.reference_start = 120
        a.mapping_quality = 50
        a.cigartuples = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # Read 3: With a mismatch (T->G at position 133 relative to ref)
        seq3 = list("ACGTACGTAC" * 5)
        seq3[3] = "G"  # Introduce mismatch at query pos 3 (ref pos 133)
        a = pysam.AlignedSegment()
        a.query_name = "read3"
        a.query_sequence = "".join(seq3)
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 130
        a.mapping_quality = 40
        a.cigartuples = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # Read 4: With insertion (2I at position 10 in query)
        a = pysam.AlignedSegment()
        a.query_name = "read4"
        a.query_sequence = "ACGTACGTAC" + "TT" + "ACGTACGTAC" * 4  # 52bp
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 200
        a.mapping_quality = 55
        a.cigartuples = [(0, 10), (1, 2), (0, 40)]  # 10M2I40M
        a.query_qualities = pysam.qualitystring_to_array("I" * 52)
        outf.write(a)

        # Read 5: With deletion (3D)
        a = pysam.AlignedSegment()
        a.query_name = "read5"
        a.query_sequence = "ACGTACGTAC" * 5  # 50bp
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 300
        a.mapping_quality = 45
        a.cigartuples = [(0, 25), (2, 3), (0, 25)]  # 25M3D25M
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # Read 6: Low quality mapping
        a = pysam.AlignedSegment()
        a.query_name = "read6_lowq"
        a.query_sequence = "ACGTACGTAC" * 3  # 30bp
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 5
        a.cigartuples = [(0, 30)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 30)
        outf.write(a)

        # Read 7: Secondary alignment (should be filtered)
        a = pysam.AlignedSegment()
        a.query_name = "read7_secondary"
        a.query_sequence = "ACGTACGTAC" * 3
        a.flag = 256  # secondary
        a.reference_id = 0
        a.reference_start = 400
        a.mapping_quality = 30
        a.cigartuples = [(0, 30)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 30)
        outf.write(a)

        # Read 8: On chr2
        a = pysam.AlignedSegment()
        a.query_name = "read8_chr2"
        a.query_sequence = "ACGTACGTAC" * 3
        a.flag = 0
        a.reference_id = 1  # chr2
        a.reference_start = 50
        a.mapping_quality = 60
        a.cigartuples = [(0, 30)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 30)
        outf.write(a)

        # Reads 9-18: Multiple overlapping reads for coverage and variant testing
        # at position 500-600 region
        for i in range(10):
            seq = list("ACGTACGTAC" * 5)
            # Reads 9-13: introduce a T->C mismatch at position 525 (query pos 25+offset)
            if i < 5:
                offset = 500 + i * 5
                mismatch_query_pos = 525 - offset
                if 0 <= mismatch_query_pos < 50:
                    seq[mismatch_query_pos] = "C"

            a = pysam.AlignedSegment()
            a.query_name = f"read_cov_{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0 if i % 2 == 0 else 16
            a.reference_id = 0
            a.reference_start = 500 + i * 5
            a.mapping_quality = 60
            a.cigartuples = [(0, 50)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 50)
            outf.write(a)

    # Sort and index
    sorted_bam = os.path.join(FIXTURES_DIR, "small_sorted.bam")
    pysam.sort("-o", sorted_bam, bam_path)

    # Replace unsorted with sorted
    os.replace(sorted_bam, bam_path)

    pysam.index(bam_path)
    return bam_path


def create_comprehensive_reference():
    """Create a 10,000bp chr1 + 5,000bp chr2 reference with targeted modifications.

    Modifications for testing:
    - Position 800-803: AAAA (4bp homopolymer, at artifact threshold)
    - Position 1800-1805: TTTTTT (6bp homopolymer, above threshold)
    - Position 2500: A (known base for multi-allele tests)
    """
    os.makedirs(FIXTURES_DIR, exist_ok=True)
    ref_path = os.path.join(FIXTURES_DIR, "comprehensive_ref.fa")

    # Build chr1 (10,000bp)
    chr1 = list("ACGTACGTAC" * 1000)  # 10,000bp
    # Homopolymer at 800-803
    for i in range(800, 804):
        chr1[i] = "A"
    # Homopolymer at 1800-1805
    for i in range(1800, 1806):
        chr1[i] = "T"
    # Known base at 2500
    chr1[2500] = "A"
    chr1_seq = "".join(chr1)

    # Build chr2 (5,000bp)
    chr2_seq = "ACGTACGTAC" * 500

    with open(ref_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(chr1_seq), 80):
            f.write(chr1_seq[i : i + 80] + "\n")
        f.write(">chr2\n")
        for i in range(0, len(chr2_seq), 80):
            f.write(chr2_seq[i : i + 80] + "\n")

    pysam.faidx(ref_path)
    return ref_path


def _make_read(
    name,
    seq,
    ref_id,
    pos,
    cigar,
    mapq=60,
    flag=0,
    quals=None,
    next_ref_id=-1,
    next_ref_start=0,
    template_length=0,
):
    """Helper to construct a pysam AlignedSegment."""
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigartuples = cigar
    if quals is None:
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    else:
        a.query_qualities = quals
    if next_ref_id >= 0:
        a.next_reference_id = next_ref_id
        a.next_reference_start = next_ref_start
        a.template_length = template_length
    return a


def _ref_seq_at(ref_str, pos, length):
    """Get reference sequence at a given 0-based position."""
    return ref_str[pos : pos + length]


def create_comprehensive_bam(ref_path: str):
    """Create a comprehensive BAM with reads covering all CIGAR ops, paired-end,
    variant evidence patterns, depth extremes, homopolymers, and multi-allele sites.

    See plan phases 1A-1G for the complete read group specification.
    """
    bam_path = os.path.join(FIXTURES_DIR, "comprehensive.bam")

    # Read the reference so we can construct matching/mismatching sequences
    fa = pysam.FastaFile(ref_path)
    chr1_ref = fa.fetch("chr1")
    fa.close()

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 5000}],
    }

    reads = []

    # ── Group A: Soft-clipped reads (positions 100-200) ──────────────

    # sc_left: 5S45M at position 100
    reads.append(
        _make_read(
            "sc_left",
            "NNNNN" + _ref_seq_at(chr1_ref, 100, 45),  # 5 clipped + 45 aligned
            0,
            100,
            [(4, 5), (0, 45)],
        )
    )

    # sc_right: 45M5S at position 110
    reads.append(
        _make_read(
            "sc_right",
            _ref_seq_at(chr1_ref, 110, 45) + "NNNNN",  # 45 aligned + 5 clipped
            0,
            110,
            [(0, 45), (4, 5)],
        )
    )

    # sc_both: 3S44M3S at position 120
    reads.append(
        _make_read(
            "sc_both",
            "NNN" + _ref_seq_at(chr1_ref, 120, 44) + "NNN",  # 3+44+3=50bp
            0,
            120,
            [(4, 3), (0, 44), (4, 3)],
        )
    )

    # sc_long: 10S40M at position 130
    reads.append(
        _make_read(
            "sc_long",
            "N" * 10 + _ref_seq_at(chr1_ref, 130, 40),  # 10+40=50bp
            0,
            130,
            [(4, 10), (0, 40)],
        )
    )

    # ── Group B: Paired-end reads (positions 300-600) ────────────────

    # pair1: Proper pair, insert 300bp
    # R1: forward, position 300, mate at 550
    reads.append(
        _make_read(
            "pair1",
            _ref_seq_at(chr1_ref, 300, 50),
            0,
            300,
            [(0, 50)],
            flag=0x63,  # paired+proper+mate_reverse+read1
            next_ref_id=0,
            next_ref_start=550,
            template_length=300,
        )
    )
    # R2: reverse, position 550, mate at 300
    reads.append(
        _make_read(
            "pair1",
            _ref_seq_at(chr1_ref, 550, 50),
            0,
            550,
            [(0, 50)],
            flag=0xA3,  # paired+proper+reverse+mate_forward+read2
            next_ref_id=0,
            next_ref_start=300,
            template_length=-300,
        )
    )

    # pair2: Discordant pair, insert 1500bp
    reads.append(
        _make_read(
            "pair2",
            _ref_seq_at(chr1_ref, 350, 50),
            0,
            350,
            [(0, 50)],
            flag=0x63,
            next_ref_id=0,
            next_ref_start=1800,
            template_length=1500,
        )
    )
    reads.append(
        _make_read(
            "pair2",
            _ref_seq_at(chr1_ref, 1800, 50),
            0,
            1800,
            [(0, 50)],
            flag=0xA3,
            next_ref_id=0,
            next_ref_start=350,
            template_length=-1500,
        )
    )

    # pair3: Chimeric (mate on chr2)
    reads.append(
        _make_read(
            "pair3",
            _ref_seq_at(chr1_ref, 400, 50),
            0,
            400,
            [(0, 50)],
            flag=0x41,  # paired+read1 (no proper)
            next_ref_id=1,  # chr2
            next_ref_start=100,
            template_length=0,
        )
    )
    reads.append(
        _make_read(
            "pair3",
            "ACGTACGTAC" * 5,
            1,  # chr2
            100,
            [(0, 50)],
            flag=0xA1,  # paired+reverse+read2
            next_ref_id=0,  # chr1
            next_ref_start=400,
            template_length=0,
        )
    )

    # ── Group C: Complex CIGAR reads (positions 700-900) ─────────────

    # cx_multi_indel: 10M2I10M3D10M1I17M at position 700
    # Ref span: 10+10+3+10+17 = 50bp (pos 700-750)
    # Query: 10+2+10+10+1+17 = 50bp
    cx_seq = (
        _ref_seq_at(chr1_ref, 700, 10)
        + "TT"  # 2bp insertion
        + _ref_seq_at(chr1_ref, 710, 10)
        # 3bp deletion (no query bases)
        + _ref_seq_at(chr1_ref, 723, 10)
        + "G"  # 1bp insertion
        + _ref_seq_at(chr1_ref, 733, 17)
    )
    reads.append(
        _make_read(
            "cx_multi_indel",
            cx_seq,
            0,
            700,
            [(0, 10), (1, 2), (0, 10), (2, 3), (0, 10), (1, 1), (0, 17)],
        )
    )

    # cx_clip_indel: 5S20M2I23M at position 750
    # Ref span: 20+23 = 43bp (pos 750-793)
    # Query: 5+20+2+23 = 50bp
    cx2_seq = (
        "NNNNN"
        + _ref_seq_at(chr1_ref, 750, 20)
        + "AA"  # 2bp insertion
        + _ref_seq_at(chr1_ref, 770, 23)
    )
    reads.append(
        _make_read(
            "cx_clip_indel",
            cx2_seq,
            0,
            750,
            [(4, 5), (0, 20), (1, 2), (0, 23)],
        )
    )

    # cx_skip: 20M50N30M at position 800 (intron skip)
    # Ref span: 20+50+30 = 100bp (pos 800-900)
    # Query: 20+30 = 50bp
    reads.append(
        _make_read(
            "cx_skip",
            _ref_seq_at(chr1_ref, 800, 20) + _ref_seq_at(chr1_ref, 870, 30),
            0,
            800,
            [(0, 20), (3, 50), (0, 30)],
        )
    )

    # cx_hard: 5H45M at position 850
    # Ref span: 45bp (pos 850-895)
    # Query: 45bp (hard clip doesn't appear in query)
    reads.append(
        _make_read(
            "cx_hard",
            _ref_seq_at(chr1_ref, 850, 45),
            0,
            850,
            [(5, 5), (0, 45)],
        )
    )

    # ── Group D: Variant evidence reads (positions 1000-1300) ────────
    # 20x coverage across 1000-1300 region

    # D1: Strand-bias variant at position 1050 (ref base from chr1_ref)
    ref_base_1050 = chr1_ref[1050]
    alt_base_1050 = "T" if ref_base_1050 != "T" else "G"
    # 8 forward reads with alt at 1050
    for i in range(8):
        seq = list(_ref_seq_at(chr1_ref, 1020 + i * 2, 50))
        qpos = 1050 - (1020 + i * 2)
        if 0 <= qpos < 50:
            seq[qpos] = alt_base_1050
        reads.append(
            _make_read(
                f"d1_alt_fwd_{i}",
                "".join(seq),
                0,
                1020 + i * 2,
                [(0, 50)],
                flag=0,  # forward
            )
        )
    # 6 forward + 6 reverse reference reads
    for i in range(12):
        reads.append(
            _make_read(
                f"d1_ref_{i}",
                _ref_seq_at(chr1_ref, 1020 + i * 2, 50),
                0,
                1020 + i * 2,
                [(0, 50)],
                flag=0 if i < 6 else 16,
            )
        )

    # D2: Position-bias variant at position 1100
    ref_base_1100 = chr1_ref[1100]
    alt_base_1100 = "C" if ref_base_1100 != "C" else "G"
    # 6 reads with mismatch at pos 1100, all at query positions 0-5 (near read start)
    for i in range(6):
        start = 1100 - i  # query pos = i (0 to 5, near start)
        seq = list(_ref_seq_at(chr1_ref, start, 50))
        qpos = 1100 - start
        if 0 <= qpos < 50:
            seq[qpos] = alt_base_1100
        reads.append(
            _make_read(
                f"d2_alt_{i}",
                "".join(seq),
                0,
                start,
                [(0, 50)],
                flag=0,
            )
        )
    # 14 reference reads spanning position 1100
    for i in range(14):
        reads.append(
            _make_read(
                f"d2_ref_{i}",
                _ref_seq_at(chr1_ref, 1080 + i * 2, 50),
                0,
                1080 + i * 2,
                [(0, 50)],
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # D3: Clean variant at position 1150
    ref_base_1150 = chr1_ref[1150]
    alt_base_1150 = "A" if ref_base_1150 != "A" else "G"
    # 5 forward + 5 reverse with mismatch at mid-read positions
    for i in range(10):
        start = 1150 - 25 + i  # query pos ~25 (mid-read)
        seq = list(_ref_seq_at(chr1_ref, start, 50))
        qpos = 1150 - start
        if 0 <= qpos < 50:
            seq[qpos] = alt_base_1150
        reads.append(
            _make_read(
                f"d3_alt_{i}",
                "".join(seq),
                0,
                start,
                [(0, 50)],
                flag=0 if i < 5 else 16,
            )
        )
    # 10 reference reads
    for i in range(10):
        reads.append(
            _make_read(
                f"d3_ref_{i}",
                _ref_seq_at(chr1_ref, 1130 + i * 2, 50),
                0,
                1130 + i * 2,
                [(0, 50)],
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # D4: Low-MAPQ variant at position 1200
    ref_base_1200 = chr1_ref[1200]
    alt_base_1200 = "T" if ref_base_1200 != "T" else "C"
    # 4 reads MAPQ 5-15 with mismatch
    for i in range(4):
        start = 1175 + i * 3
        seq = list(_ref_seq_at(chr1_ref, start, 50))
        qpos = 1200 - start
        if 0 <= qpos < 50:
            seq[qpos] = alt_base_1200
        reads.append(
            _make_read(
                f"d4_alt_lowq_{i}",
                "".join(seq),
                0,
                start,
                [(0, 50)],
                mapq=5 + i * 3,  # MAPQ 5, 8, 11, 14
            )
        )
    # 1 read MAPQ 60 with mismatch
    seq = list(_ref_seq_at(chr1_ref, 1175, 50))
    qpos = 1200 - 1175
    if 0 <= qpos < 50:
        seq[qpos] = alt_base_1200
    reads.append(
        _make_read(
            "d4_alt_highq",
            "".join(seq),
            0,
            1175,
            [(0, 50)],
            mapq=60,
        )
    )
    # 15 reference reads MAPQ 60
    for i in range(15):
        reads.append(
            _make_read(
                f"d4_ref_{i}",
                _ref_seq_at(chr1_ref, 1175 + i * 2, 50),
                0,
                1175 + i * 2,
                [(0, 50)],
                mapq=60,
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # ── Group E: Depth extremes (positions 2000-2200) ────────────────

    # 2000-2050: 1x depth (single read, no variant possible at min_depth=2)
    reads.append(
        _make_read(
            "depth_1x",
            _ref_seq_at(chr1_ref, 2000, 50),
            0,
            2000,
            [(0, 50)],
        )
    )

    # 2050-2100: 2x depth, both with A→G mismatch at 2075
    ref_base_2075 = chr1_ref[2075]
    alt_base_2075 = "G" if ref_base_2075 != "G" else "C"
    for i in range(2):
        seq = list(_ref_seq_at(chr1_ref, 2050, 50))
        qpos = 2075 - 2050
        seq[qpos] = alt_base_2075
        reads.append(
            _make_read(
                f"depth_2x_{i}",
                "".join(seq),
                0,
                2050,
                [(0, 50)],
                flag=0 if i == 0 else 16,
            )
        )

    # 2100-2150: 50x depth
    for i in range(50):
        reads.append(
            _make_read(
                f"depth_50x_{i}",
                _ref_seq_at(chr1_ref, 2100, 50),
                0,
                2100,
                [(0, 50)],
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # 2150-2200: 0x depth (no reads — gap)

    # ── Group F: Homopolymer context (positions 1800-1900) ───────────
    # Reference has TTTTTT at 1800-1805. 3 of 10 reads have T→A at 1802
    for i in range(10):
        seq = list(_ref_seq_at(chr1_ref, 1780 + i * 2, 50))
        if i < 3:
            qpos = 1802 - (1780 + i * 2)
            if 0 <= qpos < 50:
                seq[qpos] = "A"
        reads.append(
            _make_read(
                f"homo_{i}",
                "".join(seq),
                0,
                1780 + i * 2,
                [(0, 50)],
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # ── Group G: Multi-allele at position 2500 ──────────────────────
    # ref=A at 2500, 3 reads T, 3 reads G, 14 reads A
    for i in range(20):
        seq = list(_ref_seq_at(chr1_ref, 2475, 50))
        qpos = 2500 - 2475  # = 25
        if i < 3:
            seq[qpos] = "T"
        elif i < 6:
            seq[qpos] = "G"
        # else: keep reference A
        reads.append(
            _make_read(
                f"multi_{i}",
                "".join(seq),
                0,
                2475,
                [(0, 50)],
                flag=0 if i % 2 == 0 else 16,
            )
        )

    # ── Write BAM ────────────────────────────────────────────────────
    # Sort reads by (reference_id, reference_start) before writing
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)

    # Index
    pysam.index(bam_path)
    return bam_path


def create_empty_bam():
    """Create an empty BAM file (no reads, just header)."""
    bam_path = os.path.join(FIXTURES_DIR, "empty.bam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header):
        pass  # No reads

    pysam.index(bam_path)
    return bam_path


if __name__ == "__main__":
    ref = create_reference_fasta()
    print(f"Created reference: {ref}")

    bam = create_small_bam(ref)
    print(f"Created BAM: {bam}")

    empty = create_empty_bam()
    print(f"Created empty BAM: {empty}")

    comp_ref = create_comprehensive_reference()
    print(f"Created comprehensive reference: {comp_ref}")

    comp_bam = create_comprehensive_bam(comp_ref)
    print(f"Created comprehensive BAM: {comp_bam}")

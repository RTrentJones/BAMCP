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

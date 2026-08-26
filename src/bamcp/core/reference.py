"""Reference genome registry and detection."""

from __future__ import annotations

# Contig naming styles. A reference is only usable for a BAM that names its
# contigs the same way: fetching "1" from a "chr1"-style FASTA raises KeyError.
CONTIG_STYLE_CHR = "chr"  # UCSC-style: chr1, chr2, chrM
CONTIG_STYLE_NOCHR = "nochr"  # Ensembl/1000G-style: 1, 2, MT

# Known genome builds with detection signatures and public references.
#
# Every URL here MUST be openable by pysam.FastaFile for RANDOM ACCESS, which
# requires both:
#   1. an uncompressed FASTA, or one compressed with BGZF ("bgzip") — plain gzip
#      is not seekable and htslib rejects it; and
#   2. a faidx index served next to it (`.fai`, plus `.gzi` for bgzip).
#
# The UCSC bigZips files this registry used to advertise
# (hg19.fa.gz / hg38.fa.gz) satisfy NEITHER — they are plain gzip (gzip FLG byte
# 0x00, no FEXTRA/BC subfield) and have no published `.fai`. Handing one to
# pysam fails with a bare "error when opening file <url>", so `list_contigs`
# was steering callers straight into an error that no argument could fix.
# Re-verify (HTTP 206 on all three) before changing any URL below.
GENOME_BUILDS: dict[str, dict] = {
    "GRCh38": {
        "aliases": ["hg38", "grch38", "grch38.p14", "grch38.p13"],
        "chr1_length": 248956422,
        "description": "Human genome build 38 (Dec 2013)",
        "fasta_urls": {
            # bgzip + .fai + .gzi
            CONTIG_STYLE_NOCHR: (
                "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna_index/"
                "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
            ),
            # uncompressed FASTA + .fai
            CONTIG_STYLE_CHR: (
                "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
                "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
            ),
        },
    },
    "GRCh37": {
        "aliases": ["hg19", "grch37", "b37", "hs37d5"],
        "chr1_length": 249250621,
        "description": "Human genome build 37 (Feb 2009)",
        "fasta_urls": {
            # bgzip + .fai + .gzi — the reference the 1000 Genomes phase-3 BAMs
            # were actually aligned against, so contigs line up exactly.
            CONTIG_STYLE_NOCHR: (
                "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
                "phase2_reference_assembly_sequence/hs37d5.fa.gz"
            ),
            # No public UCSC-style GRCh37 FASTA publishes a faidx index beside it
            # (hg19.fa.gz and the analysisSet .fa.gz are both plain gzip, no .fai),
            # so there is nothing here that would actually open. None is honest;
            # a URL that always fails is what caused this bug.
            CONTIG_STYLE_CHR: None,
        },
    },
}


def normalize_build_name(name: str) -> str | None:
    """Normalize build aliases to canonical names.

    Args:
        name: Build name or alias (e.g., "hg38", "GRCh38", "hg19")

    Returns:
        Canonical build name (e.g., "GRCh38") or None if not recognized.

    Examples:
        >>> normalize_build_name("hg38")
        "GRCh38"
        >>> normalize_build_name("hg19")
        "GRCh37"
        >>> normalize_build_name("unknown")
        None
    """
    name_lower = name.lower()

    for canonical, info in GENOME_BUILDS.items():
        if name_lower == canonical.lower():
            return canonical
        if name_lower in info["aliases"]:
            return canonical

    return None


def contig_style(contigs: list[dict]) -> str:
    """Return the naming style a BAM's contigs use.

    A suggested reference is only usable if its contig names match the BAM's:
    ``fasta.fetch("1", ...)`` against a UCSC-style FASTA raises KeyError. Judge by
    majority so a stray unprefixed decoy (hs37d5, NC_007605) in an otherwise
    chr-prefixed header does not flip the answer.
    """
    prefixed = sum(1 for c in contigs if str(c.get("name", "")).startswith("chr"))
    return CONTIG_STYLE_CHR if prefixed * 2 > len(contigs) else CONTIG_STYLE_NOCHR


def get_public_reference_url(build: str, style: str | None = None) -> str | None:
    """Get a public, index-backed FASTA URL for a genome build.

    Args:
        build: Build name or alias (e.g., "GRCh38", "hg38").
        style: Contig naming style the BAM uses — ``CONTIG_STYLE_CHR`` or
            ``CONTIG_STYLE_NOCHR`` (see :func:`contig_style`). Defaults to
            ``CONTIG_STYLE_NOCHR``.

    Returns:
        URL to a public FASTA that pysam can actually open (see GENOME_BUILDS),
        or None if the build is unrecognized or no verified reference exists for
        that naming style. None means "nothing safe to suggest" — callers should
        say so rather than falling back to a URL that will fail to open.
    """
    canonical = normalize_build_name(build)
    if not canonical or canonical not in GENOME_BUILDS:
        return None
    url = GENOME_BUILDS[canonical]["fasta_urls"].get(style or CONTIG_STYLE_NOCHR)
    return str(url) if url else None


def detect_genome_build(contigs: list[dict]) -> dict:
    """Detect genome build from contig names and lengths.

    Uses chromosome 1 length as the primary detection method, which differs
    significantly between GRCh37 (249,250,621 bp) and GRCh38 (248,956,422 bp).

    Args:
        contigs: List of {"name": str, "length": int} dicts from BAM header.

    Returns:
        Dictionary with:
        - build: "GRCh38", "GRCh37", or "unknown"
        - confidence: "high", "medium", or "low"
        - evidence: List of evidence strings explaining the detection
    """
    if not contigs:
        return {
            "build": "unknown",
            "confidence": "low",
            "evidence": ["No contigs in BAM header"],
        }

    # Build a lookup of contig lengths
    contig_lengths: dict[str, int] = {}
    for contig in contigs:
        name = contig.get("name", "")
        length = contig.get("length", 0)
        contig_lengths[name] = length
        # Also store normalized name (without chr prefix)
        if name.startswith("chr"):
            contig_lengths[name[3:]] = length

    # Try to find chr1/1 length
    chr1_length = contig_lengths.get("chr1") or contig_lengths.get("1")

    if chr1_length:
        # Check against known builds
        for build_name, info in GENOME_BUILDS.items():
            if chr1_length == info["chr1_length"]:
                return {
                    "build": build_name,
                    "confidence": "high",
                    "evidence": [f"chr1 length ({chr1_length:,}) matches {build_name}"],
                }

        # Chr1 exists but doesn't match known builds
        return {
            "build": "unknown",
            "confidence": "medium",
            "evidence": [
                f"chr1 length ({chr1_length:,}) does not match known human builds",
                "GRCh38 chr1 = 248,956,422; GRCh37 chr1 = 249,250,621",
            ],
        }

    # No chr1 found - might be non-human or unusual naming
    evidence = ["chr1/1 not found in contigs"]

    # Check for common contig names
    has_chr_prefix = any(c.get("name", "").startswith("chr") for c in contigs)
    if has_chr_prefix:
        evidence.append("Contigs use 'chr' prefix (common in UCSC-style references)")
    else:
        evidence.append("Contigs do not use 'chr' prefix (may be NCBI-style or non-human)")

    return {
        "build": "unknown",
        "confidence": "low",
        "evidence": evidence,
    }

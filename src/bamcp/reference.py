"""Reference genome registry and detection."""

from __future__ import annotations

# Known genome builds with detection signatures and public references
GENOME_BUILDS: dict[str, dict] = {
    "GRCh38": {
        "aliases": ["hg38", "grch38", "grch38.p14", "grch38.p13"],
        "chr1_length": 248956422,
        "description": "Human genome build 38 (Dec 2013)",
        "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    },
    "GRCh37": {
        "aliases": ["hg19", "grch37", "b37", "hs37d5"],
        "chr1_length": 249250621,
        "description": "Human genome build 37 (Feb 2009)",
        "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
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


def get_public_reference_url(build: str) -> str | None:
    """Get public FASTA URL for a genome build.

    Args:
        build: Build name or alias (e.g., "GRCh38", "hg38")

    Returns:
        URL to a public FASTA file, or None if build not recognized.
    """
    canonical = normalize_build_name(build)
    if canonical and canonical in GENOME_BUILDS:
        return str(GENOME_BUILDS[canonical]["fasta_url"])
    return None


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

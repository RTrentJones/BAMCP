"""Input validation for BAMCP tool handlers.

Provides validation for genomic regions, file paths, remote URLs (SSRF prevention),
and variant parameters (chromosome, position, alleles).
"""

from __future__ import annotations

import ipaddress
import re
import socket
from pathlib import Path
from urllib.parse import urlparse

from ..config import BAMCPConfig
from ..constants import REMOTE_FILE_SCHEMES

# Pattern for valid chromosome names (chr prefix optional, allows chr1-chr99 or 1-99, X, Y, M, MT)
CHROM_PATTERN = re.compile(r"^(chr)?(\d{1,2}|[XYM]|MT)$", re.IGNORECASE)

# Pattern for valid allele strings (only ACGTN)
ALLELE_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)

# Pattern for valid genomic region strings (e.g., chr1:1000-2000 or chr1:1,000-2,000)
REGION_PATTERN = re.compile(
    r"^(chr)?(\d{1,2}|[XYM]|MT):(\d{1,3}(?:,\d{3})*|\d{1,11})-"
    r"(\d{1,3}(?:,\d{3})*|\d{1,11})$",
    re.IGNORECASE,
)

# Input length limits
MAX_FILE_PATH_LENGTH = 2048
MAX_REGION_LENGTH = 100

# Allowed file extensions
ALLOWED_FILE_EXTENSIONS = (".bam", ".cram")
ALLOWED_VARIANT_FILE_EXTENSIONS = (".vcf", ".vcf.gz", ".bcf")
ALLOWED_REFERENCE_EXTENSIONS = (".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz")


def _ext_target(path_or_url: str) -> str:
    """Lower-cased string to run the file-extension check against.

    For remote URLs, use only the URL *path* so a query string (e.g. a signed S3/GCS
    ``?X-Amz-Signature=…`` URL) does not defeat the ``endswith`` extension check.
    """
    if "://" in path_or_url:
        return urlparse(path_or_url).path.lower()
    return path_or_url.lower()


def _is_private_ip(addr: str) -> bool:
    """Check if an IP address is private, loopback, or link-local.

    Blocks access to internal networks and cloud metadata endpoints
    (e.g., 169.254.169.254) to prevent SSRF attacks.
    """
    try:
        ip = ipaddress.ip_address(addr)
    except ValueError:
        return True  # If we can't parse it, block it

    return ip.is_private or ip.is_loopback or ip.is_link_local or ip.is_reserved


def validate_remote_url(url: str, config: BAMCPConfig) -> None:
    """Validate a remote URL is safe to access (SSRF prevention).

    Resolves the hostname to IP addresses and blocks private/internal ranges.

    Args:
        url: The remote URL to validate.
        config: Server configuration.

    Raises:
        ValueError: If the URL targets a private/internal IP or is not allowed.

    Note:
        This resolves DNS and the real fetch happens separately, so a hostname
        could be rebound to an internal address between the two calls (a
        TOCTOU / DNS-rebinding gap). Callers re-validate derived URLs (e.g. index
        downloads) to narrow the window, but the only robust defense for
        untrusted input is the ``allowed_remote_hosts`` allow-list. Host matching
        is exact per-host (no implicit subdomain wildcards).
    """
    parsed = urlparse(url)
    hostname = parsed.hostname

    if not hostname:
        raise ValueError("Remote URL has no hostname")

    # Check domain allowlist if configured. Hostnames are case-insensitive
    # (urlparse already lower-cases parsed.hostname), so normalize the allow-list
    # too — otherwise a config entry like "Example.com" would never match.
    if config.allowed_remote_hosts:
        allowed_hosts = {h.lower() for h in config.allowed_remote_hosts}
        if hostname.lower() not in allowed_hosts:
            raise ValueError(f"Host '{hostname}' is not in the allowed remote hosts list")

    # Resolve hostname to IP addresses and check each one
    try:
        addr_infos = socket.getaddrinfo(hostname, parsed.port or 443, proto=socket.IPPROTO_TCP)
    except socket.gaierror as e:
        raise ValueError(f"Cannot resolve hostname '{hostname}': {e}") from e

    if not addr_infos:
        raise ValueError(f"No addresses found for hostname '{hostname}'")

    for addr_info in addr_infos:
        ip_str = str(addr_info[4][0])
        if _is_private_ip(ip_str):
            raise ValueError("Remote URL resolves to private/internal address (blocked)")


def validate_region(region: str) -> None:
    """Validate a genomic region string format.

    Args:
        region: Region string like 'chr1:1000-2000'.

    Raises:
        ValueError: If the region string is malformed.
    """
    if len(region) > MAX_REGION_LENGTH:
        raise ValueError(f"Region string too long (max {MAX_REGION_LENGTH} characters)")

    if not REGION_PATTERN.match(region):
        raise ValueError(f"Invalid region format: '{region}'. Expected format: chr1:1000-2000")


def validate_path(file_path: str, config: BAMCPConfig) -> None:
    """Validate that the file path is allowed by configuration.

    Args:
        file_path: Path or URL to validate.
        config: Server configuration.

    Raises:
        ValueError: If the path is not allowed.
    """
    if len(file_path) > MAX_FILE_PATH_LENGTH:
        raise ValueError(f"File path too long (max {MAX_FILE_PATH_LENGTH} characters)")

    # Check for remote URLs
    if "://" in file_path:
        if not config.allow_remote_files:
            raise ValueError("Remote files are disabled")

        if not file_path.startswith(REMOTE_FILE_SCHEMES):
            raise ValueError(f"Scheme not supported for remote file: {file_path}")

        # SSRF prevention: validate the URL target
        validate_remote_url(file_path, config)
        return

    # Validate file extension for local files
    lower_path = file_path.lower()
    if not any(lower_path.endswith(ext) for ext in ALLOWED_FILE_EXTENSIONS):
        raise ValueError(f"Unsupported file type. Allowed extensions: {ALLOWED_FILE_EXTENSIONS}")

    # Check local files if restrictions are configured
    if config.allowed_directories:
        try:
            abs_path = Path(file_path).resolve()
        except OSError as e:
            raise ValueError(f"Invalid path: {file_path}") from e

        allowed = False
        for d in config.allowed_directories:
            try:
                allowed_dir = Path(d).resolve()
                if abs_path.is_relative_to(allowed_dir):
                    allowed = True
                    break
            except OSError:
                continue

        if not allowed:
            raise ValueError("Path is not in allowed directories")

    # Surface a clean error for a missing local file instead of a cryptic pysam
    # failure deeper in the call stack. Checked last so allow-list violations take
    # precedence and never leak the existence of files outside allowed directories.
    if not Path(file_path).is_file():
        raise FileNotFoundError(f"File not found: {file_path}")


def validate_variant_file_path(file_path: str, config: BAMCPConfig) -> None:
    """Validate an optional VCF/BCF path used to overlay called variants."""
    if len(file_path) > MAX_FILE_PATH_LENGTH:
        raise ValueError(f"File path too long (max {MAX_FILE_PATH_LENGTH} characters)")

    ext_target = _ext_target(file_path)
    if "://" in file_path:
        if not config.allow_remote_files:
            raise ValueError("Remote files are disabled")
        if not file_path.startswith(REMOTE_FILE_SCHEMES):
            raise ValueError(f"Scheme not supported for remote file: {file_path}")
        if not any(ext_target.endswith(ext) for ext in ALLOWED_VARIANT_FILE_EXTENSIONS):
            raise ValueError(
                "Unsupported variant file type. "
                f"Allowed extensions: {ALLOWED_VARIANT_FILE_EXTENSIONS}"
            )
        validate_remote_url(file_path, config)
        return

    if not any(ext_target.endswith(ext) for ext in ALLOWED_VARIANT_FILE_EXTENSIONS):
        raise ValueError(
            f"Unsupported variant file type. Allowed extensions: {ALLOWED_VARIANT_FILE_EXTENSIONS}"
        )

    if config.allowed_directories:
        try:
            abs_path = Path(file_path).resolve()
        except OSError as e:
            raise ValueError(f"Invalid path: {file_path}") from e

        if not any(abs_path.is_relative_to(Path(d).resolve()) for d in config.allowed_directories):
            raise ValueError("Path is not in allowed directories")


def validate_reference_path(reference_path: str, config: BAMCPConfig) -> None:
    """Validate a reference FASTA path/URL with the same policy as BAM/VCF inputs.

    Closes an SSRF gap: the ``reference`` argument is handed straight to
    ``pysam.AlignmentFile(reference_filename=…)`` and htslib will fetch it, so an
    unvalidated remote reference could target a private/internal host or a cloud
    metadata endpoint. Remote references are a first-class use case (public UCSC
    FASTA for CRAM decode), so this is **allow-list aware** rather than blocking —
    it shares the SSRF core (:func:`validate_remote_url`) with the BAM/CRAM/VCF and
    index paths, so one policy governs every remote resource BAMCP hands to pysam.
    """
    if len(reference_path) > MAX_FILE_PATH_LENGTH:
        raise ValueError(f"Reference path too long (max {MAX_FILE_PATH_LENGTH} characters)")

    ext_target = _ext_target(reference_path)
    if "://" in reference_path:
        if not config.allow_remote_files:
            raise ValueError("Remote files are disabled")
        if not reference_path.startswith(REMOTE_FILE_SCHEMES):
            raise ValueError(f"Scheme not supported for remote reference: {reference_path}")
        if not any(ext_target.endswith(ext) for ext in ALLOWED_REFERENCE_EXTENSIONS):
            raise ValueError(
                "Unsupported reference file type. "
                f"Allowed extensions: {ALLOWED_REFERENCE_EXTENSIONS}"
            )
        validate_remote_url(reference_path, config)
        return

    if not any(ext_target.endswith(ext) for ext in ALLOWED_REFERENCE_EXTENSIONS):
        raise ValueError(
            f"Unsupported reference file type. Allowed extensions: {ALLOWED_REFERENCE_EXTENSIONS}"
        )

    if config.allowed_directories:
        try:
            abs_path = Path(reference_path).resolve()
        except OSError as e:
            raise ValueError(f"Invalid reference path: {reference_path}") from e
        if not any(abs_path.is_relative_to(Path(d).resolve()) for d in config.allowed_directories):
            raise ValueError("Reference path is not in allowed directories")

    if not Path(reference_path).is_file():
        raise FileNotFoundError(f"Reference file not found: {reference_path}")


def validate_data_sources(
    file_path: str,
    config: BAMCPConfig,
    *,
    reference: str | None = None,
    vcf_path: str | None = None,
) -> None:
    """Validate every data source a tool call will hand to pysam, in one place.

    The single choke point for the remote-resource policy: alignment file, optional
    reference FASTA, and optional variant (VCF/BCF) overlay all pass through the same
    SSRF-aware validators so no input reaches htslib unchecked.
    """
    validate_path(file_path, config)
    if reference:
        validate_reference_path(reference, config)
    if vcf_path:
        validate_variant_file_path(vcf_path, config)


def validate_variant_input(chrom: str, pos: int, ref: str, alt: str) -> str | None:
    """Validate variant lookup input parameters.

    Returns:
        Error message if validation fails, None if valid.
    """
    if not CHROM_PATTERN.match(chrom):
        return f"Invalid chromosome: {chrom}"
    if pos < 1:
        return f"Position must be positive, got {pos}"
    if not ALLELE_PATTERN.match(ref):
        return f"Invalid reference allele: {ref}"
    if not ALLELE_PATTERN.match(alt):
        return f"Invalid alternate allele: {alt}"
    if len(ref) > 1000 or len(alt) > 1000:
        return "Allele length exceeds maximum (1000bp)"
    return None

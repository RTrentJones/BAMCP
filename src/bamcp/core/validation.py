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

# Pattern for valid genomic region strings (e.g., chr1:1000-2000)
REGION_PATTERN = re.compile(r"^(chr)?(\d{1,2}|[XYM]|MT):(\d{1,11})-(\d{1,11})$", re.IGNORECASE)

# Input length limits
MAX_FILE_PATH_LENGTH = 2048
MAX_REGION_LENGTH = 100

# Allowed BAM/CRAM file extensions
ALLOWED_FILE_EXTENSIONS = (".bam", ".cram")


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
    """
    parsed = urlparse(url)
    hostname = parsed.hostname

    if not hostname:
        raise ValueError("Remote URL has no hostname")

    # Check domain allowlist if configured
    if config.allowed_remote_hosts and hostname not in config.allowed_remote_hosts:
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

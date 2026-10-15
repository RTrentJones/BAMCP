"""Input validation for BAMCP tool handlers.

Provides validation for genomic regions, file paths, remote URLs (SSRF prevention),
and variant parameters (chromosome, position, alleles).
"""

from __future__ import annotations

import ipaddress
import logging
import re
import socket
from pathlib import Path
from urllib.parse import urlparse

from ..config import BAMCPConfig
from ..constants import REMOTE_FILE_SCHEMES

logger = logging.getLogger(__name__)

# Pattern for valid chromosome names (chr prefix optional, allows chr1-chr99 or 1-99, X, Y, M, MT)
CHROM_PATTERN = re.compile(r"^(chr)?(\d{1,2}|[XYM]|MT)$", re.IGNORECASE)

# Pattern for valid allele strings (only ACGTN)
ALLELE_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)

# Pattern for valid genomic region strings (e.g., chr1:1000-2000 or chr1:1,000-2,000).
# The contig follows the SAM spec's reference-name grammar ([0-9A-Za-z!#$%&+./:;?@^_|~-]
# then [0-9A-Za-z!#$%&*+./:;=?@^_|~-]*), so non-primary and non-human contigs (chrUn_*,
# GL000*, chrEBV, HLA-A*01:01:01:01, ...) are accepted; whitespace, quotes, brackets,
# backslashes and control characters are not. The contig may itself contain ':' (HLA), so the
# coordinates are anchored to the END of the string — split on the last ':' when parsing.
REGION_PATTERN = re.compile(
    r"^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*:"
    r"(\d{1,3}(?:,\d{3})*|\d{1,11})(?:-(\d{1,3}(?:,\d{3})*|\d{1,11}))?\Z"
)

# Input length limits
MAX_FILE_PATH_LENGTH = 2048
MAX_REGION_LENGTH = 100

# Allowed file extensions
ALLOWED_FILE_EXTENSIONS = (".bam", ".cram")
ALLOWED_VARIANT_FILE_EXTENSIONS = (".vcf", ".vcf.gz", ".bcf")
ALLOWED_REFERENCE_EXTENSIONS = (".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz")

# Remote access is restricted to the default web ports (an explicit :80/:443 is fine).
ALLOWED_REMOTE_PORTS = (80, 443)
_REMOTE_SCHEMES = tuple(s.split("://", 1)[0] for s in REMOTE_FILE_SCHEMES)


def _ext_target(path_or_url: str) -> str:
    """Lower-cased string to run the file-extension check against.

    For remote URLs, use only the URL *path* so a query string (e.g. a signed S3/GCS
    ``?X-Amz-Signature=…`` URL) does not defeat the ``endswith`` extension check.
    """
    if "://" in path_or_url:
        return urlparse(path_or_url).path.lower()
    return path_or_url.lower()


def is_cram(path_or_url: str) -> bool:
    """True when a local path or remote URL names a CRAM (vs BAM) file.

    Judged on the URL *path* (not the query string of a signed URL) and
    case-insensitively, so ``https://h/S1.CRAM?sig=…`` is recognized as CRAM.
    """
    return _ext_target(path_or_url).endswith(".cram")


def _is_private_ip(addr: str) -> bool:
    """Check if an IP address is anything other than a public unicast address.

    Blocks private, loopback, link-local (incl. the 169.254.169.254 cloud metadata
    endpoint), CGNAT (100.64.0.0/10), multicast, unspecified, and reserved ranges —
    i.e. anything that is not ``is_global`` — to prevent SSRF attacks. IPv4-mapped
    IPv6 addresses (``::ffff:10.0.0.1``) are judged by their embedded IPv4 address.
    """
    try:
        ip = ipaddress.ip_address(addr)
    except ValueError:
        return True  # If we can't parse it, block it

    if isinstance(ip, ipaddress.IPv6Address) and ip.ipv4_mapped is not None:
        ip = ip.ipv4_mapped

    return (
        not ip.is_global
        or ip.is_private
        or ip.is_loopback
        or ip.is_link_local
        or ip.is_reserved
        or ip.is_multicast
        or ip.is_unspecified
    )


def validate_remote_url(url: str, config: BAMCPConfig) -> None:
    """Validate a remote URL is safe to access (SSRF prevention).

    Rejects fragments (htslib ``##idx##`` index smuggling), userinfo, parser-
    differential characters, non-http(s) schemes and non-web ports, then resolves
    the hostname and blocks every address that is not public unicast.

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

        Redirects: BAMCP's own index downloads follow redirects by hand and
        re-validate every hop, but the BAM/CRAM, reference and VCF URLs are fetched
        by htslib/libcurl, whose redirect handling cannot be hooked from Python. An
        allowed host that redirects (or is rebound) to an internal address is
        therefore only stopped by the allow-list — see
        :func:`warn_if_remote_unrestricted`.
    """
    # htslib reads ``<url>##idx##<index-url>`` as "load the index from <index-url>".
    # urlparse sees everything after the first '#' as a fragment, so only the first
    # host would be validated and htslib would happily fetch an arbitrary second URL
    # (e.g. the cloud metadata endpoint). A fragment is meaningless for a data file,
    # so reject '#' outright rather than trying to parse htslib's syntax.
    if "#" in url:
        raise ValueError("Remote URL must not contain '#' (fragments are not allowed)")
    # Characters that URL parsers disagree on (urlparse vs libcurl) are a classic
    # SSRF parser-differential vector; none belong in a legitimate data URL.
    if "\\" in url or any(c.isspace() or ord(c) < 0x20 or ord(c) == 0x7F for c in url):
        raise ValueError("Remote URL contains disallowed characters")

    parsed = urlparse(url)
    if parsed.scheme.lower() not in _REMOTE_SCHEMES:
        raise ValueError(f"Scheme not supported for remote URL: {parsed.scheme or '(none)'}")

    hostname = parsed.hostname

    if not hostname:
        raise ValueError("Remote URL has no hostname")

    # Userinfo ("user:pass@host") is another parser-differential vector (which '@'
    # delimits the host?) and is never needed for public genomic data.
    if "@" in parsed.netloc:
        raise ValueError("Remote URL must not contain credentials (userinfo)")

    # Only the standard web ports: an open demo must not be usable to probe
    # arbitrary services on an allowed host.
    try:
        port = parsed.port
    except ValueError as e:
        raise ValueError("Remote URL has an invalid port") from e
    if port is not None and port not in ALLOWED_REMOTE_PORTS:
        raise ValueError(f"Remote URL port {port} is not allowed (only 80/443)")

    # Check domain allowlist if configured. Hostnames are case-insensitive
    # (urlparse already lower-cases parsed.hostname), so normalize the allow-list
    # too — otherwise a config entry like "Example.com" would never match.
    if config.allowed_remote_hosts:
        allowed_hosts = {h.lower() for h in config.allowed_remote_hosts}
        if hostname.lower() not in allowed_hosts:
            # The allow-list is public policy on an open demo (it names public genomics
            # hosts, never credentials), so naming it tells the caller what WILL work.
            raise ValueError(
                f"Host '{hostname}' is not in the allowed remote hosts list. Allowed hosts: "
                + ", ".join(sorted(allowed_hosts))
            )

    # Resolve hostname to IP addresses and check each one
    try:
        addr_infos = socket.getaddrinfo(hostname, port or 443, proto=socket.IPPROTO_TCP)
    except socket.gaierror as e:
        raise ValueError(f"Cannot resolve hostname '{hostname}': {e}") from e

    if not addr_infos:
        raise ValueError(f"No addresses found for hostname '{hostname}'")

    for addr_info in addr_infos:
        ip_str = str(addr_info[4][0])
        if _is_private_ip(ip_str):
            raise ValueError("Remote URL resolves to private/internal address (blocked)")


def warn_if_remote_unrestricted(config: BAMCPConfig) -> bool:
    """Log a startup warning when remote files are on without a host allow-list.

    Per-URL validation cannot see redirects followed inside htslib/libcurl or DNS
    rebinding between check and fetch, so on an internet-facing server the
    ``BAMCP_ALLOWED_REMOTE_HOSTS`` allow-list is the real SSRF control. This only
    warns (startup must not fail); returns True when the warning was emitted.
    """
    if config.allow_remote_files and not config.allowed_remote_hosts:
        logger.warning(
            "BAMCP_ALLOW_REMOTE_FILES is enabled without BAMCP_ALLOWED_REMOTE_HOSTS: "
            "remote URLs are only checked for private/internal addresses, which does "
            "not cover redirects followed by htslib or DNS rebinding. Set an exact-host "
            "allow-list for any internet-facing deployment."
        )
        return True
    return False


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
        raise ValueError(
            f"Invalid region format: '{region}'. Expected 1-based inclusive 'chr1:1000-2000' "
            "(or a single base, 'chr1:1000')"
        )


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

        # Same extension policy as local files, judged on the URL *path* so a signed
        # URL's query string neither defeats nor breaks the check.
        if not any(_ext_target(file_path).endswith(ext) for ext in ALLOWED_FILE_EXTENSIONS):
            raise ValueError(
                f"Unsupported file type. Allowed extensions: {ALLOWED_FILE_EXTENSIONS}"
            )

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

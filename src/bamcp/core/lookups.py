"""ClinVar / gnomAD lookup handlers and their shared answer vocabulary.

Every outcome is an answer with an explicit ``status`` (found / not_found / unavailable /
disabled) and the assembly queried; only invalid input is a tool error.
"""

from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import asdict
from typing import Any

import httpx

from ..clients.gnomad import GnomadBuildError, GnomadQueryError
from ..clients.ncbi import NCBIUnavailableError
from ..config import BAMCPConfig
from ..middleware.telemetry import telemetry_wrapper
from .errors import ToolInputError
from .services import get_clinvar_client, get_gnomad_client
from .validation import check_variant_input

logger = logging.getLogger(__name__)

_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)

# Lookup outcome vocabulary for the external-database tools. The load-bearing distinction
# is not_found (genuine absence — the variant is truly not in the database) vs. unavailable
# (the lookup errored/timed out — absence is UNKNOWN and must NOT be read as "rare/benign").
# Collapsing the two is the circular trap that turns an empty gnomAD result into false
# rarity (PM2) evidence; the explicit status keeps callers honest.
_LOOKUP_FOUND = "found"
_LOOKUP_NOT_FOUND = "not_found"
_LOOKUP_UNAVAILABLE = "unavailable"
_LOOKUP_DISABLED = "disabled"


def _lookup_result(
    payload: dict[str, Any], *, status: str, disclaimer: str, build_info: dict[str, Any]
) -> dict:
    """Wrap an external-database lookup payload with an explicit ``status`` field so a
    caller can tell a genuine not-found from a failed lookup (see the vocabulary above).

    Every outcome — found, not_found, unavailable, disabled — carries ``build_info`` (the
    assembly the lookup would be / was made in), so a caller never has to guess it.
    """
    payload = {**payload, **build_info, "status": status}
    payload.setdefault("disclaimer", disclaimer)
    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


# Fixed, caller-safe reasons for an unavailable lookup. Upstream (NCBI / gnomAD) error text is
# logged server-side only — it is unbounded third-party content and never echoed to callers.
_LOOKUP_REASONS: dict[str, str] = {
    "timeout": "the upstream service did not answer in time; absence is not established",
    "network_error": "the upstream service could not be reached or returned an HTTP error; "
    "absence is not established",
    "upstream_error": "the upstream service returned an error or an incomplete/unparseable "
    "result; absence is not established",
    "unsupported_build": "the upstream service has no dataset for this genome build",
}


def _lookup_failure(
    source: str,
    reason_code: str,
    *,
    error: str,
    disclaimer: str,
    build_info: dict[str, Any],
) -> dict:
    """An ``unavailable`` lookup outcome with a fixed reason code + short message."""
    return _lookup_result(
        {
            "error": error,
            "reason_code": reason_code,
            "reason": f"{source}: {_LOOKUP_REASONS[reason_code]}",
        },
        status=_LOOKUP_UNAVAILABLE,
        disclaimer=disclaimer,
        build_info=build_info,
    )


def _lookup_disabled(source: str, *, disclaimer: str, build_info: dict[str, Any]) -> dict:
    return _lookup_result(
        {"error": f"{source} lookup is disabled"},
        status=_LOOKUP_DISABLED,
        disclaimer=disclaimer,
        build_info=build_info,
    )


def _lookup_not_found(
    source: str, args: dict[str, Any], *, disclaimer: str, build_info: dict[str, Any]
) -> dict:
    variant = f"{args['chrom']}:{args['pos']} {args['ref']}>{args['alt']}"
    return _lookup_result(
        {
            "found": False,
            "message": f"No {source} entry found for {variant} ({build_info['genome_build']})",
        },
        status=_LOOKUP_NOT_FOUND,
        disclaimer=disclaimer,
        build_info=build_info,
    )


def _lookup_build(args: dict[str, Any], config: BAMCPConfig) -> tuple[str | None, str]:
    """Resolve a lookup's genome build: explicit ``build`` arg, else the server default.

    Returns ``(canonical_build, source)``; ``canonical_build`` is None for an unrecognized
    explicit build (the caller reports invalid input rather than querying the wrong assembly).
    """
    from ..constants import normalize_genome_build

    if args.get("build"):
        return normalize_genome_build(args["build"]), "explicit"
    return normalize_genome_build(config.genome_build), "server-default"


def _validate_lookup_input(args: dict[str, Any], config: BAMCPConfig) -> tuple[str, str]:
    """Validate a variant lookup's chrom/pos/ref/alt/build; raise ToolInputError if invalid.

    Returns the resolved ``(build, build_source)``.
    """
    check_variant_input(args["chrom"], args["pos"], args["ref"], args["alt"])
    build, build_source = _lookup_build(args, config)
    if build is None:
        raise ToolInputError(
            f"Unrecognized build {args.get('build')!r}; use GRCh37/hg19 or GRCh38/hg38"
        )
    return build, build_source


@telemetry_wrapper("lookup_clinvar")
async def handle_lookup_clinvar(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in ClinVar via NCBI E-utilities.

    Returns clinical significance, review status, and associated conditions.
    """
    chrom, pos, ref, alt = args["chrom"], args["pos"], args["ref"], args["alt"]

    # Invalid input is a tool ERROR (isError), not an answer; not_found/unavailable/disabled
    # stay answers. Validate before the enabled check so bad input fails the same either way.
    build, build_source = _validate_lookup_input(args, config)
    build_info = {"genome_build": build, "build_source": build_source}
    common: dict[str, Any] = {"disclaimer": _CLINVAR_DISCLAIMER, "build_info": build_info}

    if not config.clinvar_enabled:
        return _lookup_disabled("ClinVar", **common)

    client = get_clinvar_client(config)
    variant = f"{chrom}:{pos} {ref}>{alt}"
    try:
        result = await client.lookup(chrom, pos, ref, alt, build=build)
    # TimeoutError first: on 3.11+ asyncio.TimeoutError IS the builtin TimeoutError, an
    # OSError subclass, so the network clause below would otherwise swallow it.
    except (asyncio.TimeoutError, TimeoutError, httpx.TimeoutException):
        logger.warning("ClinVar lookup timed out for %s", variant)
        return _lookup_failure("ClinVar", "timeout", error="ClinVar lookup timed out", **common)
    except NCBIUnavailableError as e:
        # In-band NCBI error, truncated search, unsupported build: absence NOT established.
        logger.warning("ClinVar lookup unavailable for %s: %s", variant, e)
        return _lookup_failure(
            "ClinVar", "upstream_error", error="ClinVar lookup unavailable", **common
        )
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError) as e:
        # Network/HTTP failure — the lookup is UNAVAILABLE, not an authoritative not-found.
        logger.warning("ClinVar lookup failed for %s: %s", variant, e)
        return _lookup_failure("ClinVar", "network_error", error="ClinVar lookup failed", **common)
    except ValueError as e:
        # An unparseable body.
        logger.warning("ClinVar lookup failed for %s: %s", variant, e)
        return _lookup_failure("ClinVar", "upstream_error", error="ClinVar lookup failed", **common)

    if result is None:
        return _lookup_not_found("ClinVar", args, **common)
    return _lookup_result(asdict(result), status=_LOOKUP_FOUND, **common)


@telemetry_wrapper("lookup_gnomad")
async def handle_lookup_gnomad(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom, pos, ref, alt = args["chrom"], args["pos"], args["ref"], args["alt"]

    # Invalid input is a tool ERROR (isError), not an answer; not_found/unavailable/disabled
    # stay answers. Validate before the enabled check so bad input fails the same either way.
    build, build_source = _validate_lookup_input(args, config)
    build_info = {"genome_build": build, "build_source": build_source}
    common: dict[str, Any] = {"disclaimer": _GNOMAD_DISCLAIMER, "build_info": build_info}

    if not config.gnomad_enabled:
        return _lookup_disabled("gnomAD", **common)

    client = get_gnomad_client(config)
    variant = f"{chrom}:{pos} {ref}>{alt}"
    try:
        result = await client.lookup(chrom, pos, ref, alt, build=build)
    # TimeoutError first: on 3.11+ asyncio.TimeoutError IS the builtin TimeoutError, an
    # OSError subclass, so the network clause below would otherwise swallow it.
    except (asyncio.TimeoutError, TimeoutError, httpx.TimeoutException):
        logger.warning("gnomAD lookup timed out for %s", variant)
        return _lookup_failure("gnomAD", "timeout", error="gnomAD lookup timed out", **common)
    except GnomadBuildError as e:
        # No gnomAD dataset for this build: unavailable, never a not-found.
        logger.warning("gnomAD build refused for %s: %s", variant, e)
        return _lookup_failure(
            "gnomAD", "unsupported_build", error="gnomAD lookup failed", **common
        )
    except GnomadQueryError as e:
        # Transient GraphQL error — the client propagates it (rather than returning None) so we
        # report "unavailable" here instead of a misleading "not found in gnomAD".
        logger.warning("gnomAD query error for %s: %s", variant, e)
        return _lookup_failure("gnomAD", "upstream_error", error="gnomAD lookup failed", **common)
    except (httpx.HTTPStatusError, httpx.RequestError, ConnectionError, OSError) as e:
        # Network/HTTP failure — the lookup is UNAVAILABLE, not an authoritative not-found.
        logger.warning("gnomAD lookup failed for %s: %s", variant, e)
        return _lookup_failure("gnomAD", "network_error", error="gnomAD lookup failed", **common)
    except ValueError as e:
        logger.warning("gnomAD lookup failed for %s: %s", variant, e)
        return _lookup_failure("gnomAD", "upstream_error", error="gnomAD lookup failed", **common)

    if result is None:
        return _lookup_not_found("gnomAD", args, **common)
    return _lookup_result(asdict(result), status=_LOOKUP_FOUND, **common)

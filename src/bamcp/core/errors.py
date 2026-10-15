"""One error envelope for every BAMCP tool.

A tool *failure* (bad input, unreadable file, unknown contig, timeout) must reach the
client as an MCP error result (``isError=true``) with a clean, caller-actionable message.
Two things go wrong without a single choke point:

- **Inconsistent envelopes.** Some handlers returned ``{"error": ...}`` as a *successful*
  result while others raised, so a client could not tell a failure from an answer.
- **Leaked operator config.** pysam/OS exceptions carry server-side paths
  (``BAMCP_REFERENCE``, the cache dir, ...) and raw ``KeyError`` reprs; FastMCP renders
  ``str(exc)`` straight to the client.

Handlers raise :class:`ToolInputError` for caller-facing failures; the server wraps every
tool in :func:`sanitize_tool_exception`, which passes those through and maps everything
else to a sanitized message. Note the vocabulary boundary: a lookup that legitimately finds
nothing (``status: not_found``) is a successful *answer*, not an error — this module is for
failures only and does not change lookup-tool semantics.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

from ..config import BAMCPConfig

logger = logging.getLogger(__name__)


class ToolInputError(Exception):
    """A caller-facing tool failure whose message is safe to show the client verbatim.

    Raise it with a message built only from caller-supplied values and static text —
    never from operator config (paths, hosts, tokens).
    """


def _config_secrets(config: BAMCPConfig | None) -> list[str]:
    """Operator-sourced strings that must never be echoed to a client."""
    if config is None:
        return []
    values = [config.reference, str(config.cache_dir) if config.cache_dir else None]
    return [v for v in values if v]


def _scrub(message: str, config: BAMCPConfig | None) -> str:
    for secret in _config_secrets(config):
        message = message.replace(secret, "<server-configured path>")
    return message


def sanitize_tool_exception(
    exc: BaseException, tool: str, config: BAMCPConfig | None = None
) -> ToolInputError:
    """Map any exception raised by a tool handler to a client-safe :class:`ToolInputError`.

    - ``ToolInputError`` passes through (already safe).
    - Timeouts keep their message (built by BAMCP from the caller's region), or a generic one.
    - ``ValueError`` is BAMCP's input-validation vocabulary (caller values only); it is kept,
      with any operator-configured path scrubbed defensively.
    - ``KeyError`` (pysam: unknown contig/reference) becomes a plain "not found" message
      instead of the raw ``KeyError 'unknown reference chrZZ'`` repr.
    - ``OSError`` (missing/unreadable file — often the operator's reference or cache) becomes
      a generic message: its text names server paths.
    - Anything else is logged server-side and reported generically.
    """
    if isinstance(exc, ToolInputError):
        return exc
    # TimeoutError is an OSError subclass (3.10+), so test it first; asyncio.TimeoutError is
    # a distinct class on 3.10.
    if isinstance(exc, (TimeoutError, asyncio.TimeoutError)):
        detail = _scrub(str(exc), config)
        return ToolInputError(detail or f"{tool} timed out; retry, or narrow the request.")
    if isinstance(exc, KeyError):
        detail = _scrub(str(exc.args[0]), config) if exc.args else ""
        found = f" ({detail})" if detail else ""
        return ToolInputError(
            f"Contig not found in the alignment or reference{found}. "
            "Check contig names with list_contigs."
        )
    if isinstance(exc, ValueError):
        return ToolInputError(_scrub(str(exc), config) or f"Invalid input for {tool}.")
    if isinstance(exc, OSError):
        logger.warning("%s: I/O error: %s", tool, exc)
        hint = ""
        text = str(exc)
        if config is not None and config.reference and config.reference in text:
            hint = " The server's default reference could not be opened."
        return ToolInputError(
            "Could not open or read an input file (missing, unreadable, unindexed, or "
            f"wrong format).{hint}"
        )
    logger.exception("%s: unexpected error", tool, exc_info=exc)
    return ToolInputError(f"Internal error while running {tool}.")


def format_argument_error(
    tool: str, exc: Exception, schema: dict[str, Any], arguments: dict[str, Any] | None
) -> str:
    """A short, caller-actionable message for arguments that failed schema validation.

    FastMCP renders a pydantic ``ValidationError`` verbatim — a multi-line dump ending in an
    ``errors.pydantic.dev`` URL that never names the parameters the tool actually takes. List
    the missing / unexpected / mistyped arguments and the tool's expected parameters instead.
    """
    props: dict[str, Any] = schema.get("properties") or {}
    required = set(schema.get("required") or [])
    problems: list[str] = []
    errors = exc.errors() if hasattr(exc, "errors") else []
    for err in errors:
        loc = ".".join(str(p) for p in err.get("loc", ()) if p is not None) or "arguments"
        if err.get("type") == "missing":
            problems.append(f"missing required '{loc}'")
        else:
            problems.append(f"'{loc}': {err.get('msg', 'invalid value')}")
    unexpected = sorted(set(arguments or {}) - set(props))
    if unexpected:
        problems.append("unexpected " + ", ".join(f"'{u}'" for u in unexpected))
    expected = ", ".join(f"{name}{' (required)' if name in required else ''}" for name in props)
    detail = "; ".join(problems) or "invalid arguments"
    return f"Invalid arguments for {tool}: {detail}. Expected parameters: {expected or 'none'}."

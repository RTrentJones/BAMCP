"""Tool-use telemetry: structured event emission for every tool invocation.

Captures one event per tool call with timing, sanitized args, result summary,
and (optionally) an OpenTelemetry span. The JSONL schema matches the
``tool_invocations`` table in archived/BAMCP_EVAL_HARNESS.md so the eval
harness can grade runs deterministically from telemetry output.

Usage:

    from bamcp.middleware.telemetry import telemetry_wrapper, configure_telemetry

    configure_telemetry(enabled=True, jsonl_path="/tmp/bamcp.jsonl")

    @telemetry_wrapper("get_variants")
    async def handle_get_variants(args: dict, config) -> dict:
        ...

When telemetry is disabled the decorator is a fast no-op pass-through.
"""

from __future__ import annotations

import functools
import hashlib
import json
import logging
import time
import uuid
from collections.abc import Awaitable, Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# Process-global telemetry config — set by configure_telemetry() at startup.
_enabled: bool = False
_jsonl_path: str = ""
_otel_tracer: Any = None
_session_id: str = ""

# Argument keys that always contain sensitive paths and must be hashed.
_HASH_ARG_KEYS = frozenset({"file_path", "reference", "index_filename", "bam_path"})

# Argument keys that should be dropped entirely (auth/credentials).
_DROP_ARG_KEYS = frozenset({"token", "access_token", "auth", "authorization", "api_key"})

# Cap on the serialized size of a single result_summary blob.
_RESULT_SUMMARY_MAX_BYTES = 1024

# Keys to surface from a JSON tool result, if present.
_RESULT_SUMMARY_KEYS = ("count", "variants", "reads", "contigs", "error", "found", "message")


def configure_telemetry(
    *,
    enabled: bool,
    jsonl_path: str = "",
    otel_enabled: bool = False,
    session_id: str | None = None,
) -> None:
    """Configure telemetry sinks for this process.

    Called once at server startup. When ``enabled`` is False the decorator
    becomes a fast pass-through.
    """
    global _enabled, _jsonl_path, _otel_tracer, _session_id
    _enabled = enabled
    _jsonl_path = jsonl_path or ""
    _session_id = session_id or uuid.uuid4().hex

    if not enabled:
        _otel_tracer = None
        return

    if jsonl_path:
        Path(jsonl_path).parent.mkdir(parents=True, exist_ok=True)

    _otel_tracer = _init_otel_tracer() if otel_enabled else None


def reset_telemetry() -> None:
    """Reset all telemetry state. Intended for tests."""
    global _enabled, _jsonl_path, _otel_tracer, _session_id
    _enabled = False
    _jsonl_path = ""
    _otel_tracer = None
    _session_id = ""


def _init_otel_tracer() -> Any:
    """Initialize an OpenTelemetry tracer. Returns None if the SDK is unavailable."""
    try:
        from opentelemetry import trace
        from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter
        from opentelemetry.sdk.resources import Resource
        from opentelemetry.sdk.trace import TracerProvider
        from opentelemetry.sdk.trace.export import BatchSpanProcessor
    except ImportError:
        logger.warning(
            "Telemetry OTel requested but opentelemetry SDK is not installed. "
            "Install with: pip install '.[telemetry]'"
        )
        return None

    current = trace.get_tracer_provider()
    if not isinstance(current, TracerProvider):
        provider = TracerProvider(resource=Resource.create({"service.name": "bamcp"}))
        provider.add_span_processor(BatchSpanProcessor(OTLPSpanExporter()))
        trace.set_tracer_provider(provider)

    return trace.get_tracer("bamcp")


def is_enabled() -> bool:
    """Whether telemetry capture is currently enabled."""
    return _enabled


def get_session_id() -> str:
    """Return the session ID for this server process."""
    return _session_id


def sanitize_args(args: dict[str, Any]) -> dict[str, Any]:
    """Sanitize tool arguments for telemetry capture.

    File-path-like values are hashed; auth-like keys are dropped; all other
    values pass through unchanged.
    """
    if not isinstance(args, dict):
        return {}
    out: dict[str, Any] = {}
    for k, v in args.items():
        key_lower = k.lower()
        if key_lower in _DROP_ARG_KEYS:
            continue
        if key_lower in _HASH_ARG_KEYS and isinstance(v, str) and v:
            out[k] = _hash_path(v)
        else:
            out[k] = v
    return out


def _hash_path(path: str) -> str:
    """Hash a file path / URL to ``sha256[:12]:filename``."""
    digest = hashlib.sha256(path.encode()).hexdigest()[:12]
    name = Path(path).name or path
    return f"{digest}:{name}"


def summarize_result(result: Any) -> dict[str, Any]:
    """Build a small result summary for the telemetry event.

    Extracts structured signals (counts, errors, UI presence) from the standard
    tool content envelope. Caps the total serialized size at
    ``_RESULT_SUMMARY_MAX_BYTES``.
    """
    summary: dict[str, Any] = {}
    if isinstance(result, dict):
        content = result.get("content")
        if isinstance(content, list) and content:
            first = content[0]
            text = first.get("text") if isinstance(first, dict) else None
            if isinstance(text, str):
                _enrich_summary_from_text(summary, text)
        if "_meta" in result:
            summary["has_ui"] = True
    elif isinstance(result, str):
        summary["text_preview"] = result[:200]

    serialized = json.dumps(summary, default=str)
    if len(serialized) > _RESULT_SUMMARY_MAX_BYTES:
        return {"truncated": True, "size": len(serialized)}
    return summary


def _enrich_summary_from_text(summary: dict[str, Any], text: str) -> None:
    """Pull count/error fields out of a JSON-shaped text payload."""
    try:
        parsed = json.loads(text)
    except (ValueError, TypeError):
        summary["text_preview"] = text[:200]
        return
    if not isinstance(parsed, dict):
        summary["text_preview"] = text[:200]
        return
    for k in _RESULT_SUMMARY_KEYS:
        if k not in parsed:
            continue
        v = parsed[k]
        if isinstance(v, list):
            summary[f"{k}_len"] = len(v)
        elif isinstance(v, (str, int, float, bool)) or v is None:
            summary[k] = v


def telemetry_wrapper(
    tool_name: str,
) -> Callable[[Callable[..., Awaitable[Any]]], Callable[..., Awaitable[Any]]]:
    """Decorate an async tool handler to capture telemetry on every call.

    Captures ``tool_name``, ``session_id``, sanitized ``args``, ``result_summary``,
    ``duration_ms``, ``ok``, ``error``, and (when OTel is configured) ``trace_id``
    / ``span_id``. Writes a JSONL row when a path is configured.
    """

    def decorator(handler: Callable[..., Awaitable[Any]]) -> Callable[..., Awaitable[Any]]:
        @functools.wraps(handler)
        async def wrapper(*args: Any, **kwargs: Any) -> Any:
            if not _enabled:
                return await handler(*args, **kwargs)

            tool_args = args[0] if args and isinstance(args[0], dict) else {}
            start = time.monotonic()
            ts = datetime.now(timezone.utc).isoformat()

            span = None
            trace_id = None
            span_id = None
            if _otel_tracer is not None:
                span = _otel_tracer.start_span(f"bamcp.{tool_name}")
                ctx = span.get_span_context()
                trace_id = f"{ctx.trace_id:032x}"
                span_id = f"{ctx.span_id:016x}"

            err: str | None = None
            ok = True
            result: Any = None
            try:
                result = await handler(*args, **kwargs)
                return result
            except Exception as e:
                ok = False
                err = f"{type(e).__name__}: {e}"
                raise
            finally:
                duration_ms = round((time.monotonic() - start) * 1000, 3)
                event = {
                    "ts": ts,
                    "session_id": _session_id,
                    "tool_name": tool_name,
                    "args": sanitize_args(tool_args),
                    "result_summary": summarize_result(result) if ok else {},
                    "duration_ms": duration_ms,
                    "ok": ok,
                    "error": err,
                    "trace_id": trace_id,
                    "span_id": span_id,
                }
                _write_event(event)
                if span is not None:
                    span.set_attribute("bamcp.tool", tool_name)
                    span.set_attribute("bamcp.duration_ms", duration_ms)
                    if not ok:
                        span.set_attribute("error", True)
                        span.set_attribute("error.message", err or "")
                    span.end()

        return wrapper

    return decorator


def _write_event(event: dict[str, Any]) -> None:
    """Append a telemetry event to the configured JSONL file, if any."""
    if not _jsonl_path:
        return
    try:
        line = json.dumps(event, default=str) + "\n"
        # O_APPEND on POSIX guarantees atomic appends for short writes.
        with open(_jsonl_path, "a", encoding="utf-8") as fh:
            fh.write(line)
    except OSError as e:
        logger.warning("Telemetry JSONL write failed: %s", e)

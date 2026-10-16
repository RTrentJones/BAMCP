"""Best-effort startup prewarm of the demo datasets (HTTP transports only).

The first call to a remote BAM is a cold start: the index (.bai) is downloaded, the remote
header is read, and the reference FASTA's .fai/.gzi are fetched. On the hosted demo that made
a newcomer's very first ``get_variants`` time out. When the server starts on an HTTP transport,
this module walks ``bamcp.demos.DEMO_DATASETS`` in the background and does that work up front:
cache each demo's index and its reference's .fai/.gzi (reused by every later open of that
FASTA), then open the BAM header and the reference in the parse pool (which also marks them
warm, so the cold-start budget isn't needed later).

It never blocks readiness, never raises, and is a no-op when the demos module is absent, when
``config.prewarm`` is false, or on the stdio transport (a single local client gains nothing
from warming someone else's demo files).
"""

from __future__ import annotations

import asyncio
import importlib
import logging
from typing import Any

from ..config import BAMCPConfig

logger = logging.getLogger(__name__)


def should_prewarm(config: BAMCPConfig) -> bool:
    """Prewarm only when enabled and serving over HTTP (never stdio)."""
    return bool(config.prewarm) and config.transport != "stdio"


def _demo_datasets() -> dict[str, Any]:
    """The demo registry, or ``{}`` when the demos module isn't available."""
    try:
        # Imported by name so this module neither requires nor type-depends on bamcp.demos.
        demos = importlib.import_module("bamcp.demos")
    except ImportError:
        return {}
    datasets = getattr(demos, "DEMO_DATASETS", None)
    return dict(datasets) if isinstance(datasets, dict) else {}


async def _prewarm_one(name: str, dataset: Any, config: BAMCPConfig) -> bool:
    """Cache one demo's index and open its header + reference .fai. True on success."""
    from .tools import _ensure_cached_index, _warm_remote_inputs
    from .validation import validate_path, validate_reference_path

    file_path = getattr(dataset, "file_path", None)
    reference = getattr(dataset, "reference_url", None)
    if not file_path:
        return False
    # Same gates a tool call would apply: a demo on a host this server doesn't allow is skipped.
    await asyncio.to_thread(validate_path, file_path, config)
    if reference:
        await asyncio.to_thread(validate_reference_path, reference, config)
    index_path = await _ensure_cached_index(file_path, config)
    await _warm_remote_inputs(file_path, reference, index_path, config)
    logger.info("Prewarmed demo %s", name)
    return True


async def prewarm_demos(config: BAMCPConfig) -> int:
    """Prewarm every demo dataset; return how many succeeded. Never raises."""
    datasets = _demo_datasets()
    if not datasets:
        logger.debug("Prewarm: no demo datasets available")
        return 0
    warmed = 0
    for name, dataset in datasets.items():
        try:
            if await _prewarm_one(name, dataset, config):
                warmed += 1
        except asyncio.CancelledError:
            raise
        except Exception as exc:  # noqa: BLE001 — best effort: log and move on
            logger.warning("Prewarm of demo %s failed (ignored): %s", name, exc)
    return warmed


def start_prewarm(config: BAMCPConfig) -> asyncio.Task[int] | None:
    """Schedule :func:`prewarm_demos` in the running loop without awaiting it.

    Returns the task (callers should keep a reference), or ``None`` when prewarm is disabled,
    on stdio, or when there is no running loop.
    """
    if not should_prewarm(config):
        return None
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        return None
    task = loop.create_task(prewarm_demos(config), name="bamcp-prewarm")
    task.add_done_callback(_log_task_failure)
    return task


def _log_task_failure(task: asyncio.Task[int]) -> None:
    if task.cancelled():
        return
    exc = task.exception()
    if exc is not None:  # pragma: no cover — prewarm_demos swallows its own errors
        logger.warning("Prewarm task failed (ignored): %s", exc)

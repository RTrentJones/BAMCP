"""Shared fixtures: one small (1 instance per case type) run set per test session."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import pytest

from bamcp.eval.visual import stimuli
from bamcp.eval.visual.stimuli import build_run_set


def tiny_png(text: str | None = None) -> bytes:
    """A real 4x4 PNG (readable by matplotlib), optionally with a tEXt chunk."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    plt.imsave(
        buf, [[0.0, 1.0], [1.0, 0.0]], format="png", metadata={"Comment": text} if text else None
    )
    return buf.getvalue()


class FakeViewerRenderer:
    calls: list[dict[str, Any]] = []

    def __init__(self, *a: Any, **k: Any) -> None:
        pass

    async def __aenter__(self) -> FakeViewerRenderer:
        return self

    async def __aexit__(self, *exc: Any) -> None:
        return None

    async def capture(self, payload: dict, mode: str, **kw: Any) -> bytes:
        FakeViewerRenderer.calls.append({"payload": payload, "mode": mode, **kw})
        return tiny_png()


@pytest.fixture(scope="session")
def fixture_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("vt_fixtures")


@pytest.fixture(scope="session")
def run_set(tmp_path_factory: pytest.TempPathFactory, fixture_dir: Path) -> Path:
    """Text + visual conditions (PNG step faked), 18 items, 2 models x 2 repeats."""
    import bamcp.eval.renderer as renderer_mod

    out = tmp_path_factory.mktemp("vt_run")
    mp = pytest.MonkeyPatch()
    mp.setattr(renderer_mod, "ViewerRenderer", FakeViewerRenderer)

    async def _available() -> None:
        return None

    mp.setattr(stimuli, "chromium_available", _available)
    mp.setattr(stimuli, "viewer_bundle_problem", lambda: None)
    try:
        build_run_set(
            out,
            seed=7,
            instances=1,
            conditions=stimuli.CONDITIONS,
            models=["m1", "m2"],
            repeats=2,
            fixture_dir=fixture_dir,
            log=lambda *_: None,
        )
    finally:
        mp.undo()
    return out

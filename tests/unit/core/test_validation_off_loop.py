"""Path/URL validation resolves DNS (socket.getaddrinfo); it must not block the event loop."""

import threading

import pytest

import bamcp.core.tools as tools
from bamcp.config import BAMCPConfig


@pytest.mark.unit
async def test_validation_runs_off_the_event_loop(monkeypatch):
    loop_thread = threading.get_ident()
    seen: list[int] = []

    def fake_validate_path(file_path, config):
        seen.append(threading.get_ident())
        raise ValueError("stop here")

    monkeypatch.setattr(tools, "validate_path", fake_validate_path)
    with pytest.raises(ValueError, match="stop here"):
        await tools.handle_get_coverage(
            {"file_path": "https://example.org/x.bam", "region": "chr1:1-10"}, BAMCPConfig()
        )
    assert seen and seen[0] != loop_thread


@pytest.mark.unit
async def test_resolve_reference_validates_off_loop(monkeypatch):
    loop_thread = threading.get_ident()
    seen: list[int] = []

    def fake_validate_reference(path, config):
        seen.append(threading.get_ident())

    monkeypatch.setattr(tools, "validate_reference_path", fake_validate_reference)
    explicit, effective = await tools._resolve_reference(
        {"reference": "https://example.org/ref.fa"}, BAMCPConfig()
    )
    assert explicit == effective == "https://example.org/ref.fa"
    assert seen and seen[0] != loop_thread
    # A config default is trusted and not re-validated.
    seen.clear()
    _, effective = await tools._resolve_reference({}, BAMCPConfig(reference="/ref.fa"))
    assert effective == "/ref.fa" and not seen

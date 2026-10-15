"""Cold-start budget for remote inputs, and the best-effort startup prewarm."""

import asyncio
import logging
import sys
import time
import types
from unittest.mock import MagicMock

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core import prewarm, tools
from bamcp.core.parsers import RegionData

pytestmark = pytest.mark.unit

BAM = "https://example.org/demo.bam"
REF = "https://example.org/ref.fa"


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    """Fresh warm set; no network for index downloads."""
    monkeypatch.setattr(tools, "_warm_inputs", type(tools._warm_inputs)())

    async def no_index(file_path, config):
        return None

    monkeypatch.setattr(tools, "_ensure_cached_index", no_index)


@pytest.fixture
def opened(monkeypatch):
    calls = []

    def fake_open(file_path, reference, index_path):
        calls.append((file_path, reference))

    monkeypatch.setattr(tools, "_open_remote_inputs_sync", fake_open)
    return calls


@pytest.fixture
def fake_fetch(monkeypatch):
    async def fetch(file_path, region, reference, config, **kwargs):
        return RegionData(
            contig="chr1",
            start=0,
            end=10,
            reads=[],
            coverage=[],
            variants=[],
            reference_sequence=None,
        )

    monkeypatch.setattr(tools, "_fetch_region_data", fetch)


# --- cold-start budget -------------------------------------------------------------------


async def test_cold_open_once_then_normal_budget(opened, fake_fetch):
    config = BAMCPConfig()
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, config)
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, config)
    assert opened == [(BAM, REF)]  # the cold open happened exactly once


async def test_only_the_cold_input_is_reopened(opened, fake_fetch):
    config = BAMCPConfig()
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, config)
    other = "https://example.org/other.fa"
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", other, config)
    assert opened == [(BAM, REF), (None, other)]


async def test_cold_budget_is_separate_from_parse_budget(monkeypatch, fake_fetch):
    """A cold open slower than the parse budget succeeds under the larger cold budget."""
    monkeypatch.setattr(tools, "BAM_PARSE_TIMEOUT_SECONDS", 0.01)
    monkeypatch.setattr(tools, "COLD_START_TIMEOUT_SECONDS", 2.0)
    monkeypatch.setattr(tools, "_open_remote_inputs_sync", lambda *a: time.sleep(0.1))
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, BAMCPConfig())
    assert BAM in tools._warm_inputs and REF in tools._warm_inputs


async def test_cold_timeout_is_diagnosable_and_stays_cold(monkeypatch, fake_fetch):
    monkeypatch.setattr(tools, "COLD_START_TIMEOUT_SECONDS", 0.05)
    monkeypatch.setattr(tools, "_open_remote_inputs_sync", lambda *a: time.sleep(0.3))
    with pytest.raises(TimeoutError, match="first use") as exc:
        await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, BAMCPConfig())
    assert "reference FASTA" in str(exc.value) and "retrying" in str(exc.value)
    assert BAM not in tools._warm_inputs
    await asyncio.sleep(0.3)


async def test_cold_open_error_is_left_to_the_real_parse(monkeypatch, fake_fetch):
    def boom(*a):
        raise OSError("error when opening file")

    monkeypatch.setattr(tools, "_open_remote_inputs_sync", boom)
    await tools._fetch_region_with_timeout(BAM, "chr1:1-10", REF, BAMCPConfig())
    assert BAM not in tools._warm_inputs


async def test_local_inputs_are_never_cold_opened(opened, fake_fetch, small_bam_path):
    await tools._fetch_region_with_timeout(small_bam_path, "chr1:1-10", None, BAMCPConfig())
    assert opened == []


def test_warm_set_is_bounded(monkeypatch):
    monkeypatch.setattr(tools, "WARM_INPUTS_MAX_ENTRIES", 2)
    for i in range(5):
        tools._mark_warm(f"https://example.org/{i}.bam")
    assert list(tools._warm_inputs) == ["https://example.org/3.bam", "https://example.org/4.bam"]


def test_open_remote_inputs_sync_opens_local_files(small_bam_path, ref_fasta_path):
    """The real opener works against local inputs (the same pysam calls as remote)."""
    tools._open_remote_inputs_sync(small_bam_path, ref_fasta_path, None)


# --- prewarm -----------------------------------------------------------------------------


def _demo(file_path, reference_url):
    return types.SimpleNamespace(file_path=file_path, reference_url=reference_url)


@pytest.fixture
def demos(monkeypatch):
    def install(datasets):
        module = types.ModuleType("bamcp.demos")
        module.DEMO_DATASETS = datasets
        monkeypatch.setitem(sys.modules, "bamcp.demos", module)

    monkeypatch.setattr("bamcp.core.validation.validate_path", lambda *a: None)
    monkeypatch.setattr("bamcp.core.validation.validate_reference_path", lambda *a: None)
    return install


def test_should_prewarm():
    assert prewarm.should_prewarm(BAMCPConfig(transport="streamable-http"))
    assert prewarm.should_prewarm(BAMCPConfig(transport="sse"))
    assert not prewarm.should_prewarm(BAMCPConfig(transport="stdio"))
    assert not prewarm.should_prewarm(BAMCPConfig(transport="streamable-http", prewarm=False))


async def test_prewarm_warms_each_demo(demos, opened):
    demos({"HG002": _demo(BAM, REF), "B": _demo("https://example.org/b.bam", None)})
    assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 2
    assert opened == [(BAM, REF), ("https://example.org/b.bam", None)]
    assert BAM in tools._warm_inputs and REF in tools._warm_inputs


async def test_prewarm_failures_are_logged_and_swallowed(demos, opened, monkeypatch, caplog):
    async def flaky_index(file_path, config):
        if "bad" in file_path:
            raise RuntimeError("index server down")
        return None

    monkeypatch.setattr(tools, "_ensure_cached_index", flaky_index)
    demos({"bad": _demo("https://example.org/bad.bam", REF), "good": _demo(BAM, REF)})
    with caplog.at_level(logging.WARNING, logger="bamcp.core.prewarm"):
        assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 1
    assert "Prewarm of demo bad failed" in caplog.text
    assert opened == [(BAM, REF)]


async def test_prewarm_skips_disallowed_demo(demos, opened, monkeypatch, caplog):
    def reject(*a):
        raise ValueError("Remote host not allowed")

    demos({"x": _demo(BAM, REF)})
    monkeypatch.setattr("bamcp.core.validation.validate_path", reject)
    with caplog.at_level(logging.WARNING, logger="bamcp.core.prewarm"):
        assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 0
    assert opened == []


async def test_prewarm_noop_when_demos_module_missing(monkeypatch, opened):
    monkeypatch.setitem(sys.modules, "bamcp.demos", None)  # forces ImportError
    assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 0
    assert opened == []


async def test_prewarm_noop_when_registry_malformed(demos, opened):
    demos(["not", "a", "dict"])
    assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 0
    demos({"no-path": _demo(None, None)})
    assert await prewarm.prewarm_demos(BAMCPConfig(transport="sse")) == 0


async def test_start_prewarm_schedules_in_background(demos, opened):
    demos({"HG002": _demo(BAM, REF)})
    task = prewarm.start_prewarm(BAMCPConfig(transport="streamable-http"))
    assert task is not None
    assert await task == 1


async def test_start_prewarm_never_on_stdio_or_disabled(demos, opened):
    demos({"HG002": _demo(BAM, REF)})
    assert prewarm.start_prewarm(BAMCPConfig(transport="stdio")) is None
    assert prewarm.start_prewarm(BAMCPConfig(transport="sse", prewarm=False)) is None
    assert opened == []


def test_start_prewarm_without_running_loop():
    assert prewarm.start_prewarm(BAMCPConfig(transport="sse")) is None


# --- __main__ wiring + config -------------------------------------------------------------


def _run_main(monkeypatch, transport):
    import uvicorn

    import bamcp.__main__ as entry

    for key in list(__import__("os").environ):
        if key.startswith("BAMCP_"):
            monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("BAMCP_TRANSPORT", transport)
    server = MagicMock()
    monkeypatch.setattr(entry, "create_server", lambda config: server)
    monkeypatch.setattr(entry.atexit, "register", lambda fn: None)
    started = []
    monkeypatch.setattr(prewarm, "start_prewarm", lambda config: started.append(config))

    async def no_serve(self):
        return None

    monkeypatch.setattr(uvicorn.Server, "serve", no_serve)
    entry.main()
    return server, started


def test_main_stdio_never_prewarms(monkeypatch):
    server, started = _run_main(monkeypatch, "stdio")
    server.run.assert_called_once_with(transport="stdio")
    assert started == []


def test_main_http_starts_prewarm(monkeypatch):
    _, started = _run_main(monkeypatch, "streamable-http")
    assert len(started) == 1 and started[0].transport == "streamable-http"


def test_config_prewarm_env(monkeypatch):
    monkeypatch.delenv("BAMCP_PREWARM", raising=False)
    assert BAMCPConfig.from_env().prewarm is True
    monkeypatch.setenv("BAMCP_PREWARM", "false")
    assert BAMCPConfig.from_env().prewarm is False
    monkeypatch.setenv("BAMCP_PREWARM", "nope")
    with pytest.raises(ValueError, match="BAMCP_PREWARM"):
        BAMCPConfig.from_env()

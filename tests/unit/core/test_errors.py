"""Tests for the shared tool error envelope (core/errors.py)."""

from __future__ import annotations

import asyncio

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core.errors import ToolInputError, sanitize_tool_exception

SECRET_REF = "/opt/secret/x.fa"


@pytest.fixture
def cfg(tmp_path):
    return BAMCPConfig(reference=SECRET_REF, cache_dir=str(tmp_path / "secret-cache"))


@pytest.mark.unit
def test_tool_input_error_passes_through(cfg):
    err = ToolInputError("bad input")
    assert sanitize_tool_exception(err, "t", cfg) is err


@pytest.mark.unit
def test_os_error_never_echoes_operator_reference(cfg):
    out = sanitize_tool_exception(OSError(f"file `{SECRET_REF}` not found"), "get_variants", cfg)
    assert SECRET_REF not in str(out)
    assert "/opt/secret" not in str(out)
    assert "default reference" in str(out)


@pytest.mark.unit
def test_os_error_hides_cache_paths(cfg):
    path = f"{cfg.cache_dir}/abcd1234/x.bai"
    out = sanitize_tool_exception(FileNotFoundError(2, "No such file", path), "t", cfg)
    assert cfg.cache_dir not in str(out)
    assert "Could not open or read an input file" in str(out)


@pytest.mark.unit
def test_key_error_is_clean_not_a_repr(cfg):
    out = str(sanitize_tool_exception(KeyError("unknown reference chrZZ"), "jump_to", cfg))
    assert "KeyError" not in out
    assert "'unknown reference" not in out  # no raw repr quoting
    assert "chrZZ" in out and "list_contigs" in out


@pytest.mark.unit
def test_key_error_without_args(cfg):
    assert "Contig not found" in str(sanitize_tool_exception(KeyError(), "t", cfg))


@pytest.mark.unit
def test_value_error_kept_but_scrubbed(cfg):
    out = sanitize_tool_exception(ValueError(f"bad thing near {SECRET_REF}"), "t", cfg)
    assert "bad thing" in str(out) and SECRET_REF not in str(out)
    assert "Invalid input for t" in str(sanitize_tool_exception(ValueError(), "t", cfg))


@pytest.mark.unit
@pytest.mark.parametrize("exc", [TimeoutError(), asyncio.TimeoutError()])
def test_bare_timeout_gets_a_message(cfg, exc):
    # str(asyncio.TimeoutError()) is "" — the client must still see something actionable.
    out = str(sanitize_tool_exception(exc, "visualize_region", cfg))
    assert "visualize_region timed out" in out


@pytest.mark.unit
def test_timeout_message_kept(cfg):
    out = sanitize_tool_exception(TimeoutError("Timed out after 30s reading chr1:1-2"), "t", cfg)
    assert str(out) == "Timed out after 30s reading chr1:1-2"


@pytest.mark.unit
def test_unexpected_error_is_generic(cfg):
    out = sanitize_tool_exception(RuntimeError(f"boom {SECRET_REF}"), "t", cfg)
    assert str(out) == "Internal error while running t."


@pytest.mark.unit
def test_no_config(cfg):
    assert "Could not open" in str(sanitize_tool_exception(OSError("x"), "t", None))

"""The shared PyYAML helper: lazy import, install hint, JSON fallback for dumps."""

from __future__ import annotations

import json
import sys

import pytest

from bamcp.eval.yaml_io import dump_yaml, load_yaml_file

pytest.importorskip("yaml")


@pytest.mark.unit
def test_round_trip(tmp_path):
    p = tmp_path / "x.yaml"
    p.write_text(dump_yaml({"b": 1, "a": [1, 2]}), encoding="utf-8")
    assert load_yaml_file(p) == {"b": 1, "a": [1, 2]}
    assert p.read_text().index("b:") < p.read_text().index("a:")  # key order kept


@pytest.mark.unit
def test_without_pyyaml(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "yaml", None)
    assert json.loads(dump_yaml({"a": 1})) == {"a": 1}
    with pytest.raises(RuntimeError, match="pyyaml is required for the truth-set scorer"):
        load_yaml_file(tmp_path / "x.yaml", purpose="the truth-set scorer")

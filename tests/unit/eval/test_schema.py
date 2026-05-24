"""Unit tests for eval schema loading."""

from __future__ import annotations

from pathlib import Path

import pytest

from bamcp.eval.schema import EvalCase, apply_subset, load_cases


@pytest.mark.unit
def test_loads_marrvel_style_yaml(tmp_path: Path):
    yaml_text = """\
- case:
    name: "Gene for variant"
    category: "Gene Utilities"
    input: "What gene is associated with this variant?"
    expected: "Genes include NUTM2G."
"""
    p = tmp_path / "cases.yaml"
    p.write_text(yaml_text)
    cases = load_cases(p)
    assert len(cases) == 1
    assert cases[0].name == "Gene for variant"
    assert cases[0].category == "Gene Utilities"
    assert cases[0].tools_expected == []


@pytest.mark.unit
def test_loads_flat_yaml_without_case_wrapper(tmp_path: Path):
    yaml_text = """\
- name: "Flat case"
  input: "Q?"
  expected: "A"
  tools_expected:
    - get_variants
"""
    p = tmp_path / "cases.yaml"
    p.write_text(yaml_text)
    cases = load_cases(p)
    assert len(cases) == 1
    assert cases[0].tools_expected == ["get_variants"]


@pytest.mark.unit
def test_loads_bamcp_extra_fields(tmp_path: Path):
    yaml_text = """\
- case:
    name: "Render comparison"
    input: "Render this region"
    expected: "Reads visible"
    bam_fixture: "tests/fixtures/x.bam"
    rendering_modes: ["expanded", "dv-strips"]
"""
    p = tmp_path / "cases.yaml"
    p.write_text(yaml_text)
    cases = load_cases(p)
    assert cases[0].bam_fixture == "tests/fixtures/x.bam"
    assert cases[0].rendering_modes == ["expanded", "dv-strips"]


@pytest.mark.unit
def test_top_level_must_be_list(tmp_path: Path):
    p = tmp_path / "cases.yaml"
    p.write_text("name: not_a_list\ninput: x\nexpected: y\n")
    with pytest.raises(ValueError, match="must be a list"):
        load_cases(p)


@pytest.mark.unit
def test_apply_subset_inclusive():
    cases = [EvalCase(name=f"c{i}", input="q", expected="a") for i in range(5)]
    s = apply_subset(cases, (2, 4))
    assert [c.name for c in s] == ["c1", "c2", "c3"]


@pytest.mark.unit
def test_apply_subset_none_returns_all():
    cases = [EvalCase(name="x", input="q", expected="a")]
    assert apply_subset(cases, None) == cases


@pytest.mark.unit
def test_seed_test_cases_file_is_valid():
    """The shipped test_cases.yaml must parse cleanly so `make eval` works."""
    path = Path(__file__).parent.parent.parent / "eval" / "test_cases.yaml"
    cases = load_cases(path)
    assert cases, "Seed cases file must not be empty"
    for c in cases:
        assert c.name
        assert c.input
        assert c.expected

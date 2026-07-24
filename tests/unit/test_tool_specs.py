"""Tests for the canonical tool-spec registry."""

from __future__ import annotations

import pytest

from bamcp.eval.router import InProcessRouter
from bamcp.tool_specs import TOOL_SPECS, get_tool_spec


@pytest.mark.unit
def test_every_routed_tool_has_a_spec():
    """The eval router's tools must all be described in the registry.

    This is the fidelity guard: if a tool is added to the router without a spec,
    the eval would silently fall back to a placeholder schema.
    """
    router_tools = set(InProcessRouter().list_tools())
    missing = router_tools - set(TOOL_SPECS)
    assert not missing, f"tools missing a spec: {sorted(missing)}"


@pytest.mark.unit
def test_specs_are_wellformed_json_schema():
    for name, spec in TOOL_SPECS.items():
        assert spec["description"], f"{name} has an empty description"
        schema = spec["input_schema"]
        assert schema["type"] == "object"
        assert "properties" in schema
        # Every required field must be a declared property.
        props = set(schema["properties"])
        for req in schema.get("required", []):
            assert req in props, f"{name}: required field {req!r} not in properties"


@pytest.mark.unit
def test_get_tool_spec_returns_none_for_unknown():
    assert get_tool_spec("does_not_exist") is None
    assert get_tool_spec("get_variants") is not None

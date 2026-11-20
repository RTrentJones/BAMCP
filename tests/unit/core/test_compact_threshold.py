"""The compact-mode cut-off is one constant, sent to the viewer in the payload."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from bamcp.constants import COMPACT_REGION_THRESHOLD_BP
from bamcp.core.parsers import RegionData
from bamcp.core.serialization import serialize_region_data


def _data(span: int) -> RegionData:
    return RegionData(contig="chr1", start=0, end=span, reads=[], coverage=[0] * span, variants=[])


@pytest.mark.unit
def test_payload_carries_threshold():
    payload = serialize_region_data(_data(10))
    assert payload["compact_threshold_bp"] == COMPACT_REGION_THRESHOLD_BP


@pytest.mark.unit
def test_auto_compact_uses_the_constant():
    assert serialize_region_data(_data(COMPACT_REGION_THRESHOLD_BP))["compact"] is False
    assert serialize_region_data(_data(COMPACT_REGION_THRESHOLD_BP + 1))["compact"] is True


@pytest.mark.unit
def test_viewer_fallback_matches_server():
    """The viewer's fallback (for payloads without the field) tracks the server constant."""
    constants_ts = Path(__file__).parents[3] / "src/bamcp/static/constants.ts"
    m = re.search(r"DEFAULT_COMPACT_THRESHOLD_BP\s*=\s*(\d+)", constants_ts.read_text())
    assert m, "DEFAULT_COMPACT_THRESHOLD_BP not found in constants.ts"
    assert int(m.group(1)) == COMPACT_REGION_THRESHOLD_BP

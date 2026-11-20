"""The one PyYAML entry point for the eval harness.

PyYAML ships in the ``[eval]`` extra only, so it is imported lazily: the server and the
CI lint job (``[dev]`` only) never need it. Its missing stubs are covered by one
``[[tool.mypy.overrides]]`` entry in ``pyproject.toml``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _require_yaml(purpose: str) -> Any:
    try:
        import yaml
    except ImportError as e:  # pragma: no cover — install hint
        raise RuntimeError(
            f"pyyaml is required for {purpose}. Install with 'pip install \".[eval]\"'."
        ) from e
    return yaml


def load_yaml_file(path: str | Path, *, purpose: str = "the eval harness") -> Any:
    """``yaml.safe_load`` a UTF-8 file; a missing PyYAML raises RuntimeError naming ``purpose``."""
    return _require_yaml(purpose).safe_load(Path(path).read_text(encoding="utf-8"))


def dump_yaml(obj: object) -> str:
    """Block-style YAML for ``obj`` (key order kept); JSON if PyYAML is not installed."""
    try:
        import yaml
    except ImportError:
        # Hard fallback so a report is still written without pyyaml.
        return json.dumps(obj, indent=2)
    return str(yaml.safe_dump(obj, sort_keys=False, default_flow_style=False))

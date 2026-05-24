#!/usr/bin/env python3
"""Entry point for ``make eval``. Delegates to bamcp.eval.cli."""

from __future__ import annotations

import sys

from bamcp.eval.cli import main

if __name__ == "__main__":
    sys.exit(main())

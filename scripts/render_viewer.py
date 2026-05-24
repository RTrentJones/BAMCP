#!/usr/bin/env python3
"""Entry point for ``make render-viewer``. Delegates to bamcp.eval.render_cli."""

from __future__ import annotations

import sys

from bamcp.eval.render_cli import main

if __name__ == "__main__":
    sys.exit(main())

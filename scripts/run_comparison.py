#!/usr/bin/env python3
"""Entry point for ``make eval-matrix``. Delegates to bamcp.eval.compare.

Model comparison + tool-use ablation. Needs ANTHROPIC_API_KEY / OPENAI_API_KEY
for real providers; ``--provider mock`` runs the whole pipeline offline.
"""

from __future__ import annotations

import sys

from bamcp.eval.compare import main

if __name__ == "__main__":
    sys.exit(main())

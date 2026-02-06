#!/usr/bin/env python3
"""Dev launcher for BAMCP.

Usage:
    # MCP Inspector (web UI on http://localhost:6274):
    ./dev_server.py

    # Or directly:
    python dev_server.py
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(ROOT, "tests", "fixtures")
BAM = os.path.join(FIXTURES, "small.bam")
REF = os.path.join(FIXTURES, "ref.fa")

# Ensure fixtures exist
if not os.path.exists(BAM):
    print("Creating test fixtures...")
    subprocess.check_call([sys.executable, os.path.join(ROOT, "tests", "create_fixtures.py")])

# Set defaults for manual testing
os.environ.setdefault("BAMCP_REFERENCE", REF)

if __name__ == "__main__":
    print(f"BAM fixture:  {BAM}")
    print(f"Reference:    {REF}")
    print(f"Test region:  chr1:500-600  (10 reads, variant T>C at pos 525)")
    print(f"All contigs:  chr1 (1000bp), chr2 (500bp)")
    print()
    print("Starting MCP Inspector...")
    print()
    subprocess.check_call(
        ["mcp", "dev", "dev_server.py:mcp", "-e", ROOT, "--with", "typer"],
        cwd=ROOT,
    )
else:
    # Imported by `mcp dev dev_server.py:mcp`
    from bamcp.server import create_server

    mcp = create_server()

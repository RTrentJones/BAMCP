"""Create the visual-vs-text eval fixtures: tests/fixtures/vt_eval.bam (+ .bai) and
vt_eval_ref.fa (+ .fai). Generated, not committed (see .gitignore).

    python tests/create_visual_eval_fixtures.py [--seed 1] [--instances 3] [--out-dir DIR]

The generator lives in :mod:`bamcp.eval.visual.fixtures` so the stimulus builder can
re-derive the (label-bearing) item plan from the seed; this script is a thin shell.
Same seed → byte-identical files.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from bamcp.eval.visual.fixtures import create_visual_eval_fixtures  # noqa: E402

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")

__all__ = ["create_visual_eval_fixtures"]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--instances", type=int, default=3)
    p.add_argument("--out-dir", default=FIXTURES_DIR)
    args = p.parse_args(argv)
    for name, path in create_visual_eval_fixtures(args.seed, args.out_dir, args.instances).items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

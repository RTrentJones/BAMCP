"""``python -m bamcp.eval.visual`` — build, lint, self-test and grade the visual-vs-text eval.

Subcommands::

    fixtures  write tests/fixtures/vt_eval.{bam,bam.bai} + vt_eval_ref.fa{,.fai}
    build     fixtures + stimuli + labels + runs.jsonl for a run set
    lint      leak-lint a run set's stim/ tree
    selftest  oracle = 100%, random ≈ chance, invalid = wrong
    grade     grade a directory of answers → results.csv, summary.json, charts
    dry       build (all conditions) + lint + selftest — the hermetic CI gate
    select-real  scan GIAB HG002 (network) and freeze the real-tier site manifest
    build-real   blinded fixtures + stimuli + labels + runs.jsonl from a site manifest
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .stimuli import CONDITIONS, PILOT_CONDITIONS


def _add_build_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--instances", type=int, default=3, help="seeded instances per case type")
    p.add_argument("--fixture-dir", default=None, help="default: tests/fixtures")
    png = p.add_mutually_exclusive_group()
    png.add_argument(
        "--png",
        dest="png",
        action="store_true",
        default=None,
        help="require the PNG step (fail without Chromium)",
    )
    png.add_argument(
        "--no-png",
        dest="png",
        action="store_false",
        help="skip the PNG step and the visual conditions",
    )


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m bamcp.eval.visual",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fixtures", help="write the vt_eval fixtures")
    f.add_argument("--seed", type=int, default=1)
    f.add_argument("--instances", type=int, default=3)
    f.add_argument("--fixture-dir", default=None)

    b = sub.add_parser("build", help="build a run set")
    b.add_argument("--out", required=True)
    b.add_argument("--conditions", nargs="+", default=list(PILOT_CONDITIONS), choices=CONDITIONS)
    b.add_argument("--models", nargs="+", default=["model"])
    b.add_argument("--repeats", type=int, default=1)
    b.add_argument(
        "--n-items",
        type=int,
        default=None,
        help="stratified subset (all case types' instance 0 first)",
    )
    _add_build_args(b)

    lint = sub.add_parser("lint", help="leak-lint a run set")
    lint.add_argument("--run-dir", required=True)
    lint.add_argument("--bam", default=None, help="default: the BAM recorded in build.json")

    st = sub.add_parser("selftest", help="oracle / random / invalid self-tests")
    st.add_argument("--run-dir", required=True)
    st.add_argument("--no-charts", action="store_true")

    g = sub.add_parser("grade", help="grade answers")
    g.add_argument("--run-dir", required=True)
    g.add_argument(
        "--answers-dir",
        default=None,
        help="default: answer_path from runs.jsonl, relative to --run-dir",
    )
    g.add_argument("--out", default=None, help="default: <run-dir>/results")
    g.add_argument(
        "--exclude",
        default=None,
        help="file with one run_id per line to drop (pre-registered exclusions)",
    )
    g.add_argument("--no-charts", action="store_true")

    sr = sub.add_parser("select-real", help="freeze the real-tier site manifest (network)")
    sr.add_argument("--out", required=True, help="sites .jsonl to write")
    sr.add_argument("--seed", type=int, required=True)
    sr.add_argument("--confident-bed", required=True, help="local copy of the GIAB benchmark BED")
    sr.add_argument("--repeat-bed", required=True, help="local copy of the repeat stratification")
    sr.add_argument("--chrom", default="chr20")

    br = sub.add_parser("build-real", help="build the real-tier run set from a site manifest")
    br.add_argument("--out", required=True)
    br.add_argument("--sites", required=True)
    br.add_argument("--seed", type=int, required=True)
    br.add_argument("--conditions", nargs="+", default=list(PILOT_CONDITIONS), choices=CONDITIONS)
    br.add_argument("--models", nargs="+", default=["model"])
    br.add_argument("--repeats", type=int, default=1)

    d = sub.add_parser("dry", help="hermetic gate: build + lint + selftest")
    d.add_argument("--out", default=".eval-results/visual-dry")
    _add_build_args(d)
    return p


def _lint(run_dir: str, bam: str | None = None) -> bool:
    from .leaklint import lint_run_set

    rep = lint_run_set(run_dir, bam)
    for prob in rep.problems:
        print(f"LEAK: {prob}")
    print(f"leak-lint: {rep.files_checked} files checked, {len(rep.problems)} problem(s)")
    return rep.ok


def _selftest(run_dir: str, charts: bool) -> bool:
    from .selftest import run_selftests

    ok = True
    for r in run_selftests(run_dir, charts=charts):
        print(f"{'PASS' if r.ok else 'FAIL'} {r.name}: {r.detail}")
        ok &= r.ok
    return ok


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)

    if args.cmd == "fixtures":
        from .fixtures import create_visual_eval_fixtures

        for k, v in create_visual_eval_fixtures(
            args.seed, args.fixture_dir, args.instances
        ).items():
            print(f"{k}: {v}")
        return 0

    if args.cmd in ("build", "dry"):
        from .stimuli import build_run_set

        conditions = CONDITIONS if args.cmd == "dry" else args.conditions
        models = ["selftest"] if args.cmd == "dry" else args.models
        res = build_run_set(
            args.out,
            seed=args.seed,
            conditions=conditions,
            models=models,
            repeats=1 if args.cmd == "dry" else args.repeats,
            n_items=None if args.cmd == "dry" else args.n_items,
            instances=args.instances,
            png=args.png,
            fixture_dir=args.fixture_dir,
        )
        print(
            f"built {res.n_stimuli} stimuli for {res.n_items} items, {res.n_runs} runs → "
            f"{res.out} (PNGs: {'yes' if res.pngs else 'no'})"
        )
        if args.cmd == "build":
            return 0 if _lint(str(res.out)) else 1
        ok = _lint(str(res.out))
        ok &= _selftest(str(res.out), charts=True)
        print("eval-visual-dry: OK" if ok else "eval-visual-dry: FAILED")
        return 0 if ok else 1

    if args.cmd == "select-real":
        from .real import select_sites, sites_to_jsonl

        sites = select_sites(
            args.seed,
            chrom=args.chrom,
            confident_bed=args.confident_bed,
            repeat_bed=args.repeat_bed,
        )
        digest = sites_to_jsonl(sites, args.out)
        print(f"{len(sites)} sites → {args.out} (sha256 {digest})")
        return 0

    if args.cmd == "build-real":
        from .real import build_real_run_set

        res = build_real_run_set(
            args.out,
            args.sites,
            seed=args.seed,
            conditions=args.conditions,
            models=args.models,
            repeats=args.repeats,
        )
        print(
            f"built {res.n_stimuli} stimuli for {res.n_items} items, {res.n_runs} runs → {res.out}"
        )
        return 0 if _lint(str(res.out)) else 1

    if args.cmd == "lint":
        return 0 if _lint(args.run_dir, args.bam) else 1

    if args.cmd == "selftest":
        return 0 if _selftest(args.run_dir, charts=not args.no_charts) else 1

    if args.cmd == "grade":
        from .grade import grade_run_set

        exclude = set()
        if args.exclude:
            exclude = {x.strip() for x in Path(args.exclude).read_text().splitlines() if x.strip()}
        s = grade_run_set(
            args.run_dir,
            answers_dir=args.answers_dir,
            out_dir=args.out,
            exclude=exclude,
            charts=not args.no_charts,
        )
        for key, g in s["groups"].items():
            lo, hi = g["primary_joint_ci95"]
            acc = g["primary_joint_accuracy"]
            print(
                f"{key:30s} joint={_fmt(acc)} [{_fmt(lo)}, {_fmt(hi)}] "
                f"valid={g['n_valid']}/{g['n_runs']}"
            )
        for key, c in s["contrasts"].items():
            if key.startswith("primary|"):
                lo, hi = c["ci95"]
                print(
                    f"{key:40s} delta={_fmt(c['delta_joint'])} [{_fmt(lo)}, {_fmt(hi)}] "
                    f"McNemar p={_fmt(c['mcnemar']['p'])}"
                )
        return 0
    return 2  # pragma: no cover — argparse enforces a subcommand


def _fmt(x: float | None) -> str:
    return "n/a" if x is None or x != x else f"{x:.3f}"


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

"""Round-4 analysis subcommands for ``python -m bamcp.eval.visual``.

Wire into the main parser with :func:`add_r4_analysis_commands` and dispatch with
:func:`run_r4_analysis` (returns ``None`` for a command it does not own)::

    grade-r4     grade a round-4 run set: balanced primary, C1/C2 with Holm, guardrail, TOST
    calibrate    B-only pilot → shrunken cells, knob rule, G1-G8 → calibration.json
                 (--simulate: the rule's operating characteristic under the section-7 prior)
    sufficiency  G8: cross-validated logistic AUROC per pair / mosaic band on a candidate pool
    selftest-r4  oracle = 1.0, always-real / always-artifact = 0.5 per pair, random ≈ chance

Standalone: ``python -m bamcp.eval.visual.cli_r4 <command> ...``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

R4_ANALYSIS_COMMANDS = ("grade-r4", "calibrate", "sufficiency", "selftest-r4")


def add_r4_analysis_commands(sub: Any) -> None:
    """Register the round-4 analysis subcommands on an ``add_subparsers()`` object."""
    g = sub.add_parser("grade-r4", help="grade a round-4 run set (design sections 5, 8)")
    g.add_argument("--run-dir", required=True)
    g.add_argument("--answers-dir", default=None)
    g.add_argument("--out", default=None, help="default: <run-dir>/results")
    g.add_argument("--exclude", default=None, help="file: one run_id per line to drop")
    g.add_argument(
        "--exclude-items",
        default=None,
        help="file: one item id per line to drop (adjudication sensitivity analysis)",
    )
    g.add_argument("--reference", default=None, help="default: build.json's reference")
    g.add_argument("--n-flips", type=int, default=None)
    g.add_argument("--seed", type=int, default=0)

    c = sub.add_parser("calibrate", help="B-only pilot calibration → calibration.json")
    c.add_argument("--run-dir", action="append", default=[], help="pilot run set (repeatable)")
    c.add_argument("--answers-dir", action="append", default=None)
    c.add_argument("--exclude", default=None, help="file: one run_id per line to drop")
    c.add_argument(
        "--gates", default=None, help='JSON {"G4": {"pass": ...}, "G5": ..., "G6": ..., "G8": ...}'
    )
    c.add_argument("--g8", default=None, help="a `sufficiency` output file (fills G8)")
    c.add_argument("--previous", default=None, help="the first calibration.json (re-pilot)")
    c.add_argument("--repilot", action="store_true", help="this is the B-only re-pilot")
    c.add_argument("--opus-model", default=None)
    c.add_argument("--out", default=None, help="default: <first run-dir>/calibration.json")
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("--simulate", action="store_true", help="print the operating characteristic")
    c.add_argument("--n-sims", type=int, default=500)

    s = sub.add_parser("sufficiency", help="G8 recoverability on a candidate pool")
    s.add_argument("--pool", required=True, help="candidate-pool JSONL (manifest rows)")
    s.add_argument("--out", required=True, help="G8 JSON to write")
    s.add_argument("--partition", default="pilot", help='"pilot" (default), "main" or "all"')
    s.add_argument("--items", default=None, help="items JSONL for the adjudication proxy")
    s.add_argument("--proxy-out", default=None, help="proxy posteriors JSONL to write")
    s.add_argument("--seed", type=int, default=0)

    t = sub.add_parser("selftest-r4", help="round-4 grader self-tests")
    t.add_argument("--run-dir", required=True)
    t.add_argument("--seed", type=int, default=0)


def _ids(path: str | None) -> set[str]:
    if not path:
        return set()
    return {x.strip() for x in Path(path).read_text().splitlines() if x.strip()}


def _fmt(x: Any) -> str:
    return "n/a" if x is None or x != x else f"{x:.3f}"


def _grade(args: argparse.Namespace) -> int:
    from .grade import fasta_ref_fetch, grade_run_set_r4
    from .stats import N_FLIPS

    s = grade_run_set_r4(
        args.run_dir,
        answers_dir=args.answers_dir,
        out_dir=args.out,
        exclude=_ids(args.exclude),
        exclude_items=_ids(args.exclude_items),
        n_flips=args.n_flips or N_FLIPS,
        seed=args.seed,
        ref_fetch=fasta_ref_fetch(args.reference) if args.reference else None,
    )
    for key, g in s["groups"].items():
        print(
            f"{key:34s} primary={_fmt(g['primary'])} sens={_fmt(g['sensitivity'])} "
            f"spec={_fmt(g['specificity'])} utility={_fmt(g['utility'])} "
            f"valid={g['n_valid']}/{g['n_runs']}"
        )
    for name, c in s["contrasts"].items():
        lo, hi = c["ci95"]
        print(
            f"{name} {c['a']} vs {c['b']}: delta={_fmt(c['delta'])} [{_fmt(lo)}, {_fmt(hi)}] "
            f"p={_fmt(c['p'])} holm={_fmt(c.get('holm_p'))} guardrail="
            f"{'ok' if c['guardrail']['holds'] else 'FAIL'} → {c.get('verdict', 'n/a')}"
        )
    st = s["directional_sign_test"]
    print(
        f"directional sign test: {st['n_matched']}/{st['n_counted']} matched, "
        f"p={_fmt(st['p_one_sided'])}"
    )
    return 0


def _calibrate(args: argparse.Namespace) -> int:
    from .calibrate import (
        calibrate,
        format_operating_characteristic,
        load_pilot_runs,
        simulate_operating_characteristic,
        write_calibration,
    )

    if args.simulate:
        print(
            format_operating_characteristic(
                simulate_operating_characteristic(args.n_sims, seed=args.seed)
            )
        )
        return 0
    if not args.run_dir:
        print("calibrate: --run-dir is required (or --simulate)", file=sys.stderr)
        return 2
    gates: dict[str, Any] = json.loads(Path(args.gates).read_text()) if args.gates else {}
    if args.g8:
        gates["G8"] = json.loads(Path(args.g8).read_text())
    previous = json.loads(Path(args.previous).read_text()) if args.previous else None
    runs, caption = load_pilot_runs(
        args.run_dir, answers_dirs=args.answers_dir, exclude=_ids(args.exclude)
    )
    rec = calibrate(
        runs,
        gates=gates,
        caption_runs=caption,
        previous=previous,
        repilot=args.repilot,
        opus_model=args.opus_model,
        seed=args.seed,
    )
    out = Path(args.out) if args.out else Path(args.run_dir[0]) / "calibration.json"
    write_calibration(rec, out)
    p = rec["prediction"]
    print(
        f"predicted main primary B = {p['mean']:.3f} (90% PI {p['pi90'][0]:.3f}-"
        f"{p['pi90'][1]:.3f}); per model: "
        + ", ".join(f"{m}={v:.3f}" for m, v in p["per_model"].items())
    )
    for mv in rec["gates"]["G2"]["moves"]:
        print(f"knob move: {mv['stratum']} {mv['from']} → {mv['to']} ({mv['reason']})")
    for g in ("G1", "G3", "G4", "G5", "G6", "G8"):
        print(f"{g}: {rec['gates'][g]['pass']}")
    print(
        f"decision: {rec['decision']}"
        + (f" — {'; '.join(rec['reasons'])}" if rec["reasons"] else "")
    )
    print(f"→ {out}")
    return 0


def _sufficiency(args: argparse.Namespace) -> int:
    from .grade import _clean
    from .sufficiency import adjudication_proxy, load_jsonl, sufficiency

    pool = load_jsonl(args.pool)
    partition = None if args.partition == "all" else args.partition
    res = sufficiency(pool, partition=partition, seed=args.seed)
    Path(args.out).write_text(json.dumps(_clean(res), indent=2) + "\n")
    for pair, r in res["pairs"].items():
        print(f"{pair:10s} AUROC={_fmt(r['auroc'])} n={r['n_pos']}+/{r['n_neg']}- pass={r['pass']}")
    for band, r in res["mosaic_bands"].items():
        print(f"mosaic band {band:4s} AUROC={_fmt(r['auroc'])} drop={r['drop']}")
    print(f"G8 pass={res['g8']['pass']} → {args.out}")
    if args.items:
        rows = adjudication_proxy(pool, load_jsonl(args.items), partition=partition)
        dest = (
            Path(args.proxy_out) if args.proxy_out else Path(args.out).with_suffix(".proxy.jsonl")
        )
        with dest.open("w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        print(
            f"adjudication proxy: {sum(r['undecidable'] for r in rows)}/{len(rows)} "
            f"undecidable → {dest}"
        )
    return 0


def _selftest(args: argparse.Namespace) -> int:
    from .selftest import run_selftests_r4

    ok = True
    for r in run_selftests_r4(args.run_dir, seed=args.seed):
        print(f"{'PASS' if r.ok else 'FAIL'} {r.name}: {r.detail}")
        ok &= r.ok
    return 0 if ok else 1


_DISPATCH = {
    "grade-r4": _grade,
    "calibrate": _calibrate,
    "sufficiency": _sufficiency,
    "selftest-r4": _selftest,
}


def run_r4_analysis(args: argparse.Namespace) -> int | None:
    """Run ``args.cmd`` if it is a round-4 analysis command, else ``None``."""
    fn = _DISPATCH.get(getattr(args, "cmd", ""))
    return fn(args) if fn else None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m bamcp.eval.visual.cli_r4",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_r4_analysis_commands(p.add_subparsers(dest="cmd", required=True))
    rc = run_r4_analysis(p.parse_args(argv))
    return 2 if rc is None else rc  # pragma: no branch — argparse enforces a command


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

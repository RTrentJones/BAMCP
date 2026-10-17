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
    pools-r4     round 4: scan candidate pools (network) and report them per stratum
    select-r4    round 4: draw the pilot or main site manifest from the pools (network)
    build-r4     round 4: blinded paired-end fixtures + r4 stimuli from a manifest
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .cli_r4 import add_r4_analysis_commands, run_r4_analysis
from .stimuli import CONDITIONS, PILOT_CONDITIONS


def _add_build_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--instances", type=int, default=3, help="seeded instances per case type")
    p.add_argument("--fixture-dir", default=None, help="default: tests/fixtures")
    _add_png_args(p)


def _add_png_args(p: argparse.ArgumentParser) -> None:
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
    add_r4_analysis_commands(sub)

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
    _add_png_args(br)

    pr4 = sub.add_parser("pools-r4", help="round 4: candidate pools per stratum (network)")
    pr4.add_argument("--out", required=True, help="pools report .json")
    _add_r4_source_args(pr4)

    sr4 = sub.add_parser("select-r4", help="round 4: draw a split's site manifest (network)")
    sr4.add_argument("--split", required=True, choices=("pilot", "main"))
    sr4.add_argument("--seed", type=int, required=True, help="draw seed (pools use --pool-seed)")
    sr4.add_argument("--out", required=True, help="manifest .jsonl to write")
    sr4.add_argument("--pools-report", default=None, help="also write the pools report here")
    sr4.add_argument("--draw-report", default=None, help="default: <out>.report.json")
    sr4.add_argument(
        "--knob",
        action="append",
        default=[],
        metavar="STRATUM=alt",
        help="main split: use the piloted alternative knob for a stratum (repeatable)",
    )
    sr4.add_argument(
        "--contamination",
        type=float,
        default=None,
        help="measured contamination; above 0.01 moves stratum 2's VAF floor to 0.08",
    )
    sr4.add_argument("--no-hifi", action="store_true", help="skip the HiFi label check")
    sr4.add_argument("--no-trio", action="store_true", help="skip trio recurrence")
    sr4.add_argument("--no-gnomad", action="store_true", help="skip the gnomAD AF filter")
    sr4.add_argument("--strict", action="store_true", help="exit 1 if any stratum is short")
    _add_r4_source_args(sr4)

    br4 = sub.add_parser("build-r4", help="round 4: run set from an r4 site manifest")
    br4.add_argument("--out", required=True)
    br4.add_argument("--manifest", required=True)
    br4.add_argument("--seed", type=int, required=True)
    br4.add_argument("--conditions", nargs="+", default=None, help="default: the r4 conditions")
    br4.add_argument("--models", nargs="+", default=["model"])
    br4.add_argument("--repeats", type=int, default=1)
    br4.add_argument("--bam-index", default=None, help="local .bai for the remote HG002 BAM")
    _add_png_args(br4)

    d = sub.add_parser("dry", help="hermetic gate: build + lint + selftest")
    d.add_argument("--out", default=".eval-results/visual-dry")
    _add_build_args(d)
    return p


def _add_r4_source_args(p: argparse.ArgumentParser) -> None:
    """Where round-4 selection reads from; every source defaults to its public URL."""
    p.add_argument("--chroms", nargs="+", default=["chr20"], help="scan-driven strata")
    p.add_argument(
        "--truth-chroms",
        nargs="+",
        default=None,
        help="truth-list strata (mosaic, CMRG/PSV); default: all autosomes",
    )
    p.add_argument(
        "--pairs",
        nargs="+",
        default=["mosaic", "tr", "segdup", "indel_adj", "anchor"],
        choices=("mosaic", "tr", "segdup", "indel_adj", "anchor"),
    )
    p.add_argument("--pool-seed", type=int, default=4, help="frozen seed of the pool scans")
    p.add_argument("--candidates", default=None, help="candidate-rows .jsonl: reused if present")
    p.add_argument("--cache", default=None, help="scan cache .jsonl (default: user cache dir)")
    p.add_argument("--bam-index", default=None, help="local .bai for the HG002 60x BAM")
    p.add_argument("--hifi-bai", default=None, help="local .bai for the HiFi BAM")
    for name in ("v421", "v5q", "cmrg", "mosaic"):
        p.add_argument(f"--{name}-vcf", default=None)
        p.add_argument(f"--{name}-bed", default=None)
    p.add_argument("--strat-dir", default=None, help="local v3.3 stratification BEDs")
    p.add_argument("--superdups", default=None, help="local genomicSuperDups.txt.gz")
    p.add_argument(
        "--superdup-min-identity",
        type=float,
        default=0.0,
        help="cluster only segdup pairings with fracMatch >= this (default: all, as designed)",
    )
    p.add_argument("--mosaic-database", default=None, help="mosaic database .tsv (optional)")
    p.add_argument("--limit-tr", type=int, default=1500, help="TR loci read per polarity")
    p.add_argument("--limit-adj", type=int, default=1500, help="truth indels read per polarity")
    p.add_argument("--limit-anchor", type=int, default=60)
    p.add_argument("--limit-mosaic-fp", type=int, default=3000, help="caller FPs read")
    p.add_argument("--psv-windows", type=int, default=800, help="2 kb PSV scan windows")


def _r4_rows(args: argparse.Namespace, *, checks: bool) -> tuple[list, object]:
    """Open round-4 resources and return (candidate rows, resources)."""
    from . import real, truth

    def pair(name: str) -> tuple[str, str] | None:
        v, b = getattr(args, f"{name}_vcf"), getattr(args, f"{name}_bed")
        if v is None and b is None:
            return None
        defaults = {
            "v421": (truth.V421_VCF, truth.V421_BED),
            "v5q": (truth.V5Q_VCF, truth.V5Q_BED),
            "cmrg": (truth.CMRG_VCF, truth.CMRG_BED),
            "mosaic": (truth.MOSAIC_VCF, truth.MOSAIC_BED),
        }[name]
        return (v or defaults[0], b or defaults[1])

    strat = None
    if args.strat_dir:
        strat = {
            k: str(Path(args.strat_dir) / u.rsplit("/", 1)[-1])
            for k, u in truth.STRAT_URLS.items()
            if (Path(args.strat_dir) / u.rsplit("/", 1)[-1]).exists()
        }
    truth_chroms = args.truth_chroms or list(truth.AUTOSOMES)
    use_hifi = checks and not getattr(args, "no_hifi", False)
    use_trio = checks and not getattr(args, "no_trio", False)
    res = real.open_r4_resources(
        chroms=args.chroms,
        truth_chroms=truth_chroms,
        bam_index=args.bam_index,
        v421=pair("v421"),
        v5q=pair("v5q"),
        cmrg=pair("cmrg"),
        mosaic=pair("mosaic"),
        strat=strat,
        superdups=args.superdups,
        superdup_min_identity=args.superdup_min_identity,
        hifi=(truth.HIFI_BAM, args.hifi_bai) if use_hifi else None,
        parents={n: (u, None) for n, u in truth.PARENT_BAMS.items()} if use_trio else None,
        parent_truth=dict(truth.PARENT_TRUTH) if use_trio else None,
        cache_path=args.cache or truth.cache_dir() / "r4_scan_cache.jsonl",
    )
    if args.candidates and Path(args.candidates).exists():
        rows = real.manifest_from_jsonl(args.candidates)
        print(f"{len(rows)} candidate rows from {args.candidates}")
        return rows, res
    rows = real.collect_pools_r4(
        res,
        chroms=args.chroms,
        truth_chroms=truth_chroms,
        seed=args.pool_seed,
        pairs=args.pairs,
        limits={
            "tr": args.limit_tr,
            "indel_adj": args.limit_adj,
            "anchor": args.limit_anchor,
            "mosaic_fp": args.limit_mosaic_fp,
            "psv_windows": args.psv_windows,
        },
        mosaic_database=args.mosaic_database,
    )
    if args.candidates:
        real.manifest_to_jsonl(rows, args.candidates)
    return rows, res


def _write_json(path: str | Path, obj: object) -> None:
    import json

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")


def _r4_command(args: argparse.Namespace) -> int:
    from . import real

    if args.cmd == "build-r4":
        res = real.build_r4_run_set(
            args.out,
            args.manifest,
            seed=args.seed,
            conditions=args.conditions or list(real.R4_CONDITIONS),
            models=args.models,
            repeats=args.repeats,
            png=args.png,
            bam_index=args.bam_index,
        )
        print(
            f"built {res.n_stimuli} stimuli for {res.n_items} items, {res.n_runs} runs → {res.out}"
        )
        return 0 if _lint(str(res.out)) else 1

    rows, resources = _r4_rows(args, checks=args.cmd == "select-r4")
    if args.cmd == "pools-r4":
        rep = real.pools_report_r4(rows)
        _write_json(args.out, rep)
        for s, r in rep.items():
            print(f"{s:30s} candidates={r['candidates']:5d} eligible={r['eligible']}")
        return 0

    from .strata_r4 import CONTAMINATION_VAF_FLOOR

    knobs = {}
    for kv in args.knob:
        k, _, v = kv.partition("=")
        knobs[k] = v or "alt"
    params = {}
    if args.contamination is not None and args.contamination > 0.01:
        params["mosaic_caller_fp_neg"] = {"vaf_floor_override": CONTAMINATION_VAF_FLOOR}
    if args.pools_report:
        _write_json(args.pools_report, real.pools_report_r4(rows, overrides=params))
    checks = real.LabelChecks(resources, gnomad=not args.no_gnomad)  # type: ignore[arg-type]
    chosen, report = real.select_r4(
        rows,
        split=args.split,
        seed=args.seed,
        checks=checks,
        knob_overrides=knobs,
        param_overrides=params,
        pairs=args.pairs,
    )
    report["check_failures"] = dict(checks.counts)
    digest = real.manifest_to_jsonl(chosen, args.out)
    report["manifest_sha256"] = digest
    _write_json(args.draw_report or str(args.out) + ".report.json", report)
    print(f"{len(chosen)} sites → {args.out} (sha256 {digest})")
    return 1 if args.strict and report["short"] else 0


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
    rc = run_r4_analysis(args)
    if rc is not None:
        return rc

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
            png=args.png,
        )
        print(
            f"built {res.n_stimuli} stimuli for {res.n_items} items, {res.n_runs} runs → {res.out}"
        )
        return 0 if _lint(str(res.out)) else 1

    if args.cmd in ("pools-r4", "select-r4", "build-r4"):
        return _r4_command(args)

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

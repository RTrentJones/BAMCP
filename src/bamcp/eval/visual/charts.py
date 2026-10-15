"""Static PNG charts for the visual-vs-text eval (matplotlib, Agg backend).

- ``accuracy_dots.png``  primary-endpoint accuracy per condition x model, 95% cluster CIs
- ``forest_<model>.png`` paired deltas by case type for each pre-registered contrast
- ``reliability.png``    reliability diagrams per condition (pooled over models)
- ``confusion_<condition>.png`` artifact-type confusion matrices
- ``gallery.png``        stimulus gallery (wide V screenshot per item, labelled by case type)

Palette: the dataviz reference categorical order (first three slots, validated all-pairs)
plus one sequential blue ramp; text in ink tokens, recessive grid.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")


def _plt() -> Any:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK_2,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "text.color": INK,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9,
        }
    )
    return plt


def _nan(x: Any) -> float:
    return float("nan") if x is None else float(x)


def accuracy_dots(summary: dict[str, Any], out: Path) -> Path:
    plt = _plt()
    conds = summary["conditions"]
    models = summary["models"]
    fig, ax = plt.subplots(figsize=(max(5.0, 1.3 * len(conds) + 2), 3.6))
    width = 0.6 / max(1, len(models))
    for mi, m in enumerate(models):
        xs, ys, los, his = [], [], [], []
        for ci, c in enumerate(conds):
            g = summary["groups"].get(f"{m}|{c}")
            if not g:
                continue
            y = _nan(g["primary_joint_accuracy"])
            lo, hi = (_nan(v) for v in g["primary_joint_ci95"])
            xs.append(ci - 0.3 + width * (mi + 0.5))
            ys.append(y)
            los.append(max(0.0, y - lo) if not math.isnan(lo) else 0.0)
            his.append(max(0.0, hi - y) if not math.isnan(hi) else 0.0)
        color = SERIES[mi % len(SERIES)]
        ax.errorbar(xs, ys, yerr=[los, his], fmt="none", ecolor=color, elinewidth=2, capsize=0)
        ax.plot(
            xs,
            ys,
            MARKERS[mi % len(MARKERS)],
            color=color,
            markersize=8,
            label=m,
            markeredgecolor=SURFACE,
            markeredgewidth=1.5,
            linestyle="none",
        )
    ax.set_xticks(range(len(conds)), conds)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("primary endpoint accuracy")
    ax.set_title("Joint (variant_present, artifact_verdict) accuracy, 95% cluster CI", loc="left")
    if len(models) > 1:
        ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    path = out / "accuracy_dots.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def forest_plots(summary: dict[str, Any], out: Path) -> list[Path]:
    plt = _plt()
    paths: list[Path] = []
    by_model: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for key, res in summary["contrasts"].items():
        kind, model, name = key.split("|", 2)
        if kind == "primary":
            by_model.setdefault(model, []).append((name, res))
    for model, contrasts in by_model.items():
        cases = sorted({c for _, r in contrasts for c in r["per_case_type"]})
        fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(cases) + 1.6))
        n = len(contrasts)
        for k, (name, res) in enumerate(contrasts):
            color = SERIES[k % len(SERIES)]
            ys, xs, lo, hi = [], [], [], []
            for i, c in enumerate(cases):
                pc = res["per_case_type"].get(c)
                if not pc or pc["delta"] is None or math.isnan(_nan(pc["delta"])):
                    continue
                d = _nan(pc["delta"])
                ys.append(i + (k - (n - 1) / 2) * 0.25)
                xs.append(d)
                lo.append(max(0.0, d - _nan(pc["ci95"][0])))
                hi.append(max(0.0, _nan(pc["ci95"][1]) - d))
            ax.errorbar(
                xs,
                ys,
                xerr=[lo, hi],
                fmt=MARKERS[k],
                color=color,
                markersize=6,
                elinewidth=1.5,
                capsize=0,
                label=f"{name}  (overall {_nan(res['delta_joint']):+.2f})",
            )
        ax.axvline(0, color=INK_2, linewidth=1)
        ax.set_yticks(range(len(cases)), cases)
        ax.invert_yaxis()
        ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel("paired delta in primary endpoint accuracy (A − B)")
        ax.set_title(f"Paired deltas by case type — {model}", loc="left")
        ax.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.12 - 0.6 / len(cases)),
            ncol=1,
        )
        fig.tight_layout()
        path = out / f"forest_{_slug(model)}.png"
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths


def reliability(graded: list[Any], out: Path) -> Path:
    from .stats import reliability_bins

    plt = _plt()
    conds = sorted({g.run["condition"] for g in graded})
    cols = min(3, max(1, len(conds)))
    rows = math.ceil(len(conds) / cols) or 1
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 3.0 * rows), squeeze=False)
    for i, c in enumerate(conds):
        ax = axes[i // cols][i % cols]
        valid = [g for g in graded if g.run["condition"] == c and g.answer]
        bins = reliability_bins(
            [float(g.answer["confidence"]) for g in valid], [g.joint for g in valid]
        )
        ax.plot([0, 1], [0, 1], color=INK_2, linewidth=1, linestyle="--")
        if bins:
            xs, ys, ns = zip(*bins, strict=True)
            ax.plot(
                xs, ys, "-o", color=SERIES[0], linewidth=2, markersize=6, markeredgecolor=SURFACE
            )
            for x, y, n in bins:
                ax.annotate(
                    str(n),
                    (x, y),
                    textcoords="offset points",
                    xytext=(4, -10),
                    fontsize=7,
                    color=INK_2,
                )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(c, loc="left")
        ax.set_xlabel("stated confidence")
        ax.set_ylabel("observed accuracy")
    for j in range(len(conds), rows * cols):
        axes[j // cols][j % cols].axis("off")
    fig.suptitle("Reliability (valid answers, pooled over models; n per bin)", x=0.02, ha="left")
    fig.tight_layout()
    path = out / "reliability.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def confusion(summary: dict[str, Any], out: Path) -> list[Path]:
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("seq_blue", ["#eef4fc", "#2a78d6", "#0d3a75"])
    paths: list[Path] = []
    pooled: dict[str, dict[str, dict[str, int]]] = {}
    for key, g in summary["groups"].items():
        _, cond = key.split("|", 1)
        tgt = pooled.setdefault(cond, {})
        for t, preds in g["artifact_type_confusion"].items():
            for p, n in preds.items():
                tgt.setdefault(t, {})
                tgt[t][p] = tgt[t].get(p, 0) + n
    for cond, mat in sorted(pooled.items()):
        truths = sorted(mat)
        preds = sorted({p for row in mat.values() for p in row} | set(truths))
        grid = [[mat[t].get(p, 0) for p in preds] for t in truths]
        fig, ax = plt.subplots(figsize=(0.6 * len(preds) + 2.6, 0.45 * len(truths) + 1.8))
        ax.grid(False)
        ax.imshow(grid, cmap=cmap, aspect="auto")
        vmax = max((max(r) for r in grid), default=1) or 1
        for i, row in enumerate(grid):
            for j, v in enumerate(row):
                if v:
                    ax.text(
                        j,
                        i,
                        str(v),
                        ha="center",
                        va="center",
                        fontsize=8,
                        color=SURFACE if v > vmax * 0.55 else INK,
                    )
        ax.set_xticks(range(len(preds)), preds, rotation=40, ha="right")
        ax.set_yticks(range(len(truths)), truths)
        ax.set_xlabel("answered artifact_type")
        ax.set_ylabel("true artifact_type")
        ax.set_title(f"artifact_type confusion — {cond} (pooled over models)", loc="left")
        fig.tight_layout()
        path = out / f"confusion_{_slug(cond)}.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        paths.append(path)
    return paths


def gallery(run_dir: Path, out: Path) -> Path | None:
    """Wide V screenshot per item labelled with its case type (for humans, post-hoc)."""
    import json

    runs = [json.loads(x) for x in (run_dir / "runs.jsonl").read_text().splitlines() if x]
    labels = {
        json.loads(x)["item"]: json.loads(x)
        for x in (run_dir / "labels" / "labels.jsonl").read_text().splitlines()
        if x
    }
    pngs: dict[str, Path] = {}
    for r in runs:
        p = run_dir / r["stim_dir"] / "view_wide.png"
        if r["item"] not in pngs and p.exists():
            pngs[r["item"]] = p
    if not pngs:
        return None
    plt = _plt()
    import matplotlib.image as mpimg

    items = sorted(pngs, key=lambda i: (labels[i]["case_type"], labels[i]["instance"]))
    cols = 6
    rows = math.ceil(len(items) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.9 * rows), squeeze=False)
    for k, ax in enumerate(ax for row in axes for ax in row):
        ax.axis("off")
        if k < len(items):
            img = mpimg.imread(pngs[items[k]])
            ax.imshow(img[: min(len(img), 900)])
            lab = labels[items[k]]
            ax.set_title(f"{lab['case_type']} #{lab['instance']}", fontsize=8, loc="left")
    fig.suptitle(
        "Stimulus gallery — wide V screenshot per item (labels shown post hoc)", x=0.01, ha="left"
    )
    fig.tight_layout()
    path = out / "gallery.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def _slug(s: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_") or "x"


def write_charts(
    summary: dict[str, Any], graded: list[Any], run_dir: Path, out: Path
) -> list[Path]:
    paths = [accuracy_dots(summary, out), reliability(graded, out)]
    paths += forest_plots(summary, out)
    paths += confusion(summary, out)
    g = gallery(run_dir, out)
    if g:
        paths.append(g)
    return paths

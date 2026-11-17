"""v2 figures (PNG + PDF) from saved result tables only. Palette: reference categorical slots 1-3 (blue, orange,
aqua; validated all-pairs) with neutral text ink and recessive grid.

  python model_dev/src/v2_figures.py phase1|phase2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
FIG = C.RESULTS / "figures"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 9, "axes.titleweight": "bold", "text.color": INK})


def style(ax, grid_axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=170 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("saved", FIG / f"{name}.png")


def ci_points(ax, xs, rows, color, label, val="Q_caf_expr", lo="Q_ci_low", hi="Q_ci_high", off=0.0):
    v = np.array([r[val] for r in rows]); l_ = np.array([r[lo] for r in rows]); h = np.array([r[hi] for r in rows])
    ax.errorbar(np.array(xs) + off, v, yerr=[v - l_, h - v], fmt="o", ms=5, color=color, lw=1.6, capsize=0, label=label)


def phase1():
    d = pd.read_csv(C.RESULTS / "v2_phase1" / "phase1_summary.csv")
    dec = json.loads((C.RESULTS / "v2_phase1" / "phase1_decision.json").read_text())
    b = pd.read_csv(C.RESULTS / "v2_phase1" / "b_contamination_by_adjacency.csv")
    f = pd.read_csv(C.RESULTS / "v2_phase1" / "f_target_genes_caf_vs_macrophage.csv")
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    # A: Q by distance bin
    ax = axs[0, 0]
    bins = ["touching", "near", "separated"]
    for t, color, lab, off in (("myeloid:C1Q_TAM", BLUE, "C1Q_TAM (primary)", -0.15),
                               ("S_mean", AQUA, "spillover-free macrophage genes S", 0.0),
                               ("M_mean", ORANGE, "CAF marker genes M (pure spillover)", 0.15)):
        rows = [d[(d.group == "ref") & (d.stratum == "bin") & (d.level == lv) & (d.target == t) & (d.model == "xgb")].iloc[0] for lv in bins]
        ci_points(ax, range(3), rows, color, lab, off=off)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(3), ["touching\n(<= 12.1 um)", "near\n(12.1-24.3 um)", "separated\n(> 24.3 um)"])
    ax.set_ylabel("CAF-expression contribution Q (R2)")
    ax.set_title("A  CAF-expression contribution by distance to nearest CAF")
    ax.legend(frameon=False, fontsize=7); style(ax)
    # B: contamination by bin
    ax = axs[0, 1]
    bb = b[b.stratifier == "bin"]
    x = np.arange(3)
    ax.bar(x, bb.marker_cp10k, color=BLUE, width=0.55)
    ax.errorbar(x, bb.marker_cp10k, yerr=[bb.marker_cp10k - bb.ci_low, bb.ci_high - bb.marker_cp10k], fmt="none", ecolor=INK, lw=1)
    caf = float(b[b.stratifier == "CAF cells"].marker_cp10k.iloc[0])
    ax.axhline(caf, color=ORANGE, lw=1.6, ls="--")
    ax.text(2.35, caf, f"CAFs {caf:,.0f}", color=INK, va="bottom", ha="right", fontsize=7)
    for xi, v in zip(x, bb.marker_cp10k):
        ax.text(xi, v + 60, f"{v:,.0f}", ha="center", fontsize=7, color=INK)
    rate = pd.read_csv(C.RESULTS / "v2_phase1" / "b_spillover_rate.csv")
    rt = rate[rate.bin == "touching"].iloc[0]
    ax.set_xticks(x, ["touching", "near", "separated"])
    ax.set_ylabel("CP10k, sum over the 20 marker genes M")
    ax.set_title("B  CAF-marker transcripts in macrophages\n"
                 f"(touching: {rt.spillover_rate:.0%} [{rt.ci_low:.0%}, {rt.ci_high:.0%}] of the CAF excess)")
    style(ax)
    # C: checks
    ax = axs[1, 0]
    items = [("reference\nC1Q_TAM", ("ref", "myeloid:C1Q_TAM")), ("adjusted\nC1Q_adj+s_c", ("sc", "C1Q_TAM_adj")),
             ("S genes\nspill.-free", ("ref", "S_mean")), ("M genes\nCAF markers", ("ref", "M_mean")),
             ("CD45+ q75\nC1Q_TAM", ("cd45_gt_caf75", "myeloid:C1Q_TAM")), ("CD45+ q90\nC1Q_TAM", ("cd45_gt_caf90", "myeloid:C1Q_TAM"))]
    rows = [d[(d.group == g) & (d.stratum == "all") & (d.target == t) & (d.model == "xgb")].iloc[0] for _, (g, t) in items]
    colors = [BLUE, BLUE, AQUA, ORANGE, INK2, INK2]
    for i, (r, c) in enumerate(zip(rows, colors)):
        ci_points(ax, [i], [r], c, None)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(items)), [i[0] for i in items], fontsize=7)
    ax.set_ylabel("Q (R2), all receivers")
    ax.set_title(f"C  Checks (XGBoost). Verdict: {dec['verdict']}")
    style(ax)
    # D: target genes CAF vs isolated macrophages
    ax = axs[1, 1]
    f = f.sort_values("ratio_caf_over_iso")
    y = np.arange(len(f))
    ax.scatter(f.caf_mean_cp10k_median, y, color=ORANGE, s=30, label="CAFs", zorder=3)
    ax.scatter(f.iso_mac_mean_cp10k_median, y, color=BLUE, s=30, label="isolated macrophages", zorder=3)
    ax.scatter(f.touch_mac_mean_cp10k_median, y, color=AQUA, s=30, marker="D", label="touching macrophages", zorder=3)
    for yi, r in zip(y, f.itertuples()):
        ax.plot([r.caf_mean_cp10k_median, r.iso_mac_mean_cp10k_median], [yi, yi], color=GRID, lw=2, zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(y, [f"{g}{' *' if h else ''}  ({r:.2f})" for g, h, r in zip(f.gene, f.hypothesis_gene, f.ratio_caf_over_iso)], fontsize=7)
    ax.set_xlabel("mean CP10k (median over 13 tumours; log scale)")
    ax.set_title("D  Target genes in CAFs vs macrophages\n(label: ratio CAF / isolated macrophage; * hypothesis gene)")
    ax.legend(frameon=False, fontsize=7, loc="upper right"); style(ax, "x")
    fig.suptitle("v2 Phase 1: segmentation-spillover control for the CAF -> macrophage signal (development set, 13 tumours; tumour-bootstrap 95% CIs)", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig1_spillover")


def phase2():
    q = pd.read_csv(C.RESULTS / "v2_phase2" / "phase2_quantities.csv")
    show = [("myeloid_orig", "myeloid:C1Q_TAM"), ("myeloid_orig", "myeloid:SPP1_TAM"), ("tumour_orig", "tumour:classical"),
            ("tumour_orig", "tumour:basal_like"), ("myeloid_adj_sc", "C1Q_TAM_adj"), ("myeloid_adj_sc", "SPP1_TAM_adj")]
    show = [s for s in show if ((q.group == s[0]) & (q.target == s[1])).any()]
    qs = [("neighbour gain (all - own)", "neighbour gain", INK2), ("cross-type contribution (all - same)", "cross-type", BLUE),
          ("who is nearby: other-type composition (same_plus_other_comp - same)", "other-type composition", AQUA),
          ("what they express: other-type expression (all - same_plus_other_comp)", "other-type expression", ORANGE)]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axs[0]
    for j, (qq, lab, color) in enumerate(qs):
        rows = []
        for g, t in show:
            r = q[(q.group == g) & (q.target == t) & (q.model == "xgb (genes set)") & (q.quantity == qq)]
            rows.append(r.iloc[0] if len(r) else None)
        xs = [i for i, r in enumerate(rows) if r is not None]
        ci_points(ax, xs, [r for r in rows if r is not None], color, lab, "value", "ci_low", "ci_high", off=-0.24 + 0.16 * j)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(show)), [f"{t.split(':')[-1]}\n({g.split('_')[0]}{', adjusted' if 'adj' in g else ''})" for g, t in show], fontsize=7)
    ax.set_ylabel("held-out R2 difference")
    ax.set_title("A  XGBoost (genes set): neighbour gain and its other-type parts")
    ax.legend(frameon=False, fontsize=7); style(ax)
    ax = axs[1]
    rows = []
    for g, t in show:
        r = q[(q.group == g) & (q.target == t) & (q.model == "xgb (genes set)") & q.quantity.str.startswith("unique contribution")]
        rows.append(r)
    types = ["CAF", "tumour", "myeloid", "T", "other"]
    cols = {"CAF": ORANGE, "tumour": BLUE, "myeloid": BLUE, "T": AQUA, "other": INK2}
    for i, r in enumerate(rows):
        for k, t in enumerate(types):
            rr = r[r.quantity.str.contains(f" of {t} ")]
            if rr.empty:
                continue
            rr = rr.iloc[0]
            ax.errorbar(rr.value, i + (k - 2) * 0.15, xerr=[[rr.value - rr.ci_low], [rr.ci_high - rr.value]], fmt="o", ms=4,
                        color=cols[t], lw=1.3, label=t if i == 0 or (t == "myeloid" and i == 2) else None)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(range(len(show)), [f"{t.split(':')[-1]} ({g.split('_')[0]})" for g, t in show], fontsize=7)
    ax.set_xlabel("unique contribution (all - minus type), R2")
    ax.set_title("B  Unique contribution of each other cell type")
    h, l_ = ax.get_legend_handles_labels()
    ax.legend(h, l_, frameon=False, fontsize=7, title="neighbour type (other-type only)"); style(ax, "x")
    fig.suptitle("v2 Phase 2: cross-type share on the aggregated baselines (development set; paired tumour-bootstrap 95% CIs)", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig2_crosstype")


if __name__ == "__main__":
    {"phase1": phase1, "phase2": phase2}[sys.argv[1]]()

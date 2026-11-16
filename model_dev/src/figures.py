"""Phase 7 figures for the report and the dashboard (PNG + PDF in model_dev/results/figures/).

Every figure is drawn from a generated results file (paths in each function); nothing is typed in by hand.
Also writes a checkpoint manifest (paths, sizes, SHA-256) for the git-ignored model_cache/checkpoints.
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

R = C.RESULTS
OUT = R / "figures"
TARGETS = {"myeloid": ["myeloid:C1Q_TAM", "myeloid:SPP1_TAM", "genes_pooled"],
           "tumour": ["tumour:classical", "tumour:basal_like", "genes_pooled"]}
PRETTY = {"genes_pooled": "gene vector (pooled)"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def model_comparison() -> None:
    b = pd.read_csv(R / "phase2" / "baseline_summary.csv")
    g = pd.read_csv(R / "phase3" / "graph_summary.csv")
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5))
    order = ["ridge/own", "ridge/programs", "xgb/programs", "ridge/genes", "xgb/genes", "mlp/genes", "graph/sage", "graph/attention"]
    colors = {"ridge/own": "#9e9e9e", "ridge/programs": "#6baed6", "xgb/programs": "#2171b5", "ridge/genes": "#9ecae1",
              "xgb/genes": "#08519c", "mlp/genes": "#fd8d3c", "graph/sage": "#d62728", "graph/attention": "#e377c2"}
    for i, (rec, tlist) in enumerate(TARGETS.items()):
        for jx, t in enumerate(tlist):
            ax = axes[i, jx]
            rows = []
            for r in b[(b.receiver == rec) & (b.target == t) & (b.model != "mean")].itertuples():
                if r.model == "xgb" and t == "genes_pooled":
                    continue  # XGBoost was fitted to program targets only (gene columns = mean)
                rows.append((f"{r.model}/{r.feature_set}", r.r2, r.ci_low, r.ci_high))
            for r in g[(g.receiver == rec) & (g.target == t)].itertuples():
                rows.append((f"graph/{r.arch}", r.r2, r.ci_low, r.ci_high))
            d = pd.DataFrame(rows, columns=["model", "r2", "lo", "hi"]).set_index("model").reindex([o for o in order if o in {x[0] for x in rows}])
            y = np.arange(len(d))
            ax.errorbar(d.r2, y, xerr=[d.r2 - d.lo, d.hi - d.r2], fmt="none", ecolor="#555", capsize=2, lw=1)
            ax.scatter(d.r2, y, c=[colors[m] for m in d.index], s=40, zorder=3)
            own = d.loc["ridge/own", "r2"]
            ax.axvline(own, color="#9e9e9e", ls=":", lw=1)
            ax.set_yticks(y, d.index)
            ax.set_title(f"{rec}: {PRETTY.get(t, t)}", fontsize=10, fontweight="bold")
            ax.set_xlabel("held-out R² (grouped 5-fold; 95% tumour CI)")
    fig.suptitle("Neighbourhood models vs baselines (development set, 13 tumours; dotted line = own covariates only)", fontsize=11)
    fig.tight_layout()
    save(fig, "fig1_model_comparison")


def ablations() -> None:
    a = pd.read_csv(R / "phase3" / "ablation_summary.csv")
    labels = {"same_type": "same-type only\n(= cross-type contribution)", "other_type": "other-type only",
              "composition": "composition only", "shuffled": "tile-shuffled\nneighbour expression", "no_graph": "no graph"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=False)
    for ax, rec in zip(axes, ["myeloid", "tumour"]):
        d = a[(a.receiver == rec) & a.target.isin(TARGETS[rec])]
        modes = [m for m in labels if m in set(d["mode"])]
        w = 0.8 / len(TARGETS[rec])
        for k, t in enumerate(TARGETS[rec]):
            e = d[d.target == t].set_index("mode").reindex(modes)
            x = np.arange(len(modes)) + (k - 1) * w
            ax.bar(x, e.full_minus_mode, w, yerr=[e.full_minus_mode - e.diff_ci_low, e.diff_ci_high - e.full_minus_mode],
                   capsize=2, label=PRETTY.get(t, t), color=["#d62728", "#2171b5", "#7f7f7f"][k], alpha=.85)
        ax.axhline(0, color="k", lw=.7)
        ax.set_xticks(np.arange(len(modes)), [labels[m] for m in modes], fontsize=8)
        ax.set_ylabel("R²(full graph model) − R²(ablated)\n(95% paired tumour CI)")
        ax.set_title(f"{rec} receivers", fontweight="bold")
        ax.legend(fontsize=7)
    fig.suptitle("Ablations: what the graph model loses when information is removed (SAGE, lr 3e-4, 3 seeds)", fontsize=11)
    fig.tight_layout()
    save(fig, "fig2_ablations_cross_type")


def external() -> None:
    e = pd.read_csv(R / "phase4" / "external_results.csv").drop_duplicates(["receiver", "model", "target"])
    e = e[e.target.isin(TARGETS["myeloid"] + TARGETS["tumour"]) & (e.model != "mean")]
    e = e[~((e.model.str.startswith("xgb")) & (e.target == "genes_pooled"))]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    mk = {"myeloid": "o", "tumour": "s"}
    for r in e.itertuples():
        col = "#d62728" if r.model.startswith("graph") else ("#9e9e9e" if r.model == "ridge/own" else "#2171b5")
        ax.errorbar(r.internal_r2, r.external_r2, yerr=[[r.external_r2 - r.unit_ci_low], [r.unit_ci_high - r.external_r2]],
                    fmt=mk[r.receiver], color=col, ms=5, capsize=1.5, alpha=.8)
    lim = [min(-0.1, e.external_r2.min() - .05), max(e.internal_r2.max(), e.external_r2.max()) + .05]
    ax.plot(lim, lim, "k:", lw=.8)
    ax.set_xlabel("internal R² (development, grouped CV)")
    ax.set_ylabel("external R² (GSE310352; 95% unit CI)")
    ax.set_title("Internal vs external (circles myeloid, squares tumour;\nred graph, blue aggregated baselines, grey own covariates)", fontsize=9)
    ax = axes[1]
    d = e[(e.receiver == "myeloid") & (e.target == "myeloid:C1Q_TAM")]
    for r in d.itertuples():
        ps = json.loads(r.per_slide_r2)
        col = "#d62728" if r.model.startswith("graph") else ("#9e9e9e" if r.model == "ridge/own" else None)
        ax.plot(list(ps), list(ps.values()), marker="o", label=r.model, color=col, lw=1.5 if col else 1)
    ax.set_ylabel("external R² per slide")
    ax.set_title("myeloid C1Q_TAM: per-slide external R²\n(slides 5/5b2/5b3 likely one specimen; independence UNVERIFIED)", fontsize=9)
    ax.legend(fontsize=7)
    plt.setp(ax.get_xticklabels(), rotation=30)
    fig.tight_layout()
    save(fig, "fig3_external_test")


def training_curves() -> None:
    sel = json.loads((R / "phase3" / "architecture_selection.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, rec in zip(axes, ["myeloid", "tumour"]):
        for f in sorted((R / "phase3" / "curves").glob(f"graph_{rec}_full_{sel['selected_arch']}_lr{sel['selected_lr']}_seed*_fold*.csv")):
            c = pd.read_csv(f)
            fold = int(f.stem.split("fold")[-1])
            ax.plot(c.epoch, c.val_mse, color=plt.cm.tab10(fold), alpha=.6, lw=1)
            ax.plot(c.epoch, c.train_mse, color=plt.cm.tab10(fold), alpha=.25, lw=1, ls="--")
        ax.set_title(f"{rec}: {sel['selected_arch']} lr {sel['selected_lr']} (5 folds x 3 seeds)", fontsize=10)
        ax.set_xlabel("epoch")
        ax.set_ylabel("MSE (solid = inner validation, dashed = train)")
    fig.suptitle("Graph-model training curves with early stopping (patience 8)", fontsize=11)
    fig.tight_layout()
    save(fig, "fig4_training_curves")


def importance_and_perturbation() -> None:
    p = R / "phase5" / "permutation_importance.csv"
    if not p.exists():
        return
    imp = pd.read_csv(p).sort_values("r2_drop")
    pert = pd.read_csv(R / "phase5" / "in_silico_perturbation.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.3]})
    ax = axes[0]
    d = imp[imp.r2_drop.abs() > 0.0001]
    col = ["#d62728" if g.startswith("class:myeloid") or "myeloid:" in g else ("#2171b5" if g.startswith("class:") else "#7f7f7f") for g in d.group]
    ax.barh(range(len(d)), d.r2_drop, xerr=[d.r2_drop - d.ci_low, d.ci_high - d.r2_drop], color=col, capsize=2)
    ax.set_yticks(range(len(d)), d.group.str.replace("class:", "neighbour class: ").str.replace("genes:", "neighbour genes: "), fontsize=7.5)
    ax.set_xlabel("R² lost when the group is permuted/zeroed (95% tumour CI)")
    ax.set_title("Permutation importance, myeloid gene vector (SAGE)", fontsize=10)
    ax = axes[1]
    genes = [g for g in pert.target.unique() if g != "genes_pooled" and ":" not in g]
    progs = list(pert.neighbour_program.unique())
    M = np.full((len(genes), len(progs)), np.nan)
    for i, gname in enumerate(genes):
        for j, pn in enumerate(progs):
            hi = pert[(pert.target == gname) & (pert.neighbour_program == pn) & (pert["quantile"] == 0.9)].predicted_change
            lo = pert[(pert.target == gname) & (pert.neighbour_program == pn) & (pert["quantile"] == 0.1)].predicted_change
            if len(hi) and len(lo):
                M[i, j] = float(hi.iloc[0]) - float(lo.iloc[0])
    v = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(progs)), [f"neighbour {p_}" for p_ in progs], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)), genes, fontsize=8)
    plt.colorbar(im, ax=ax, label="predicted change (log1p CP10k): 90th − 10th percentile setting")
    ax.set_title("In-silico perturbation of other-type neighbours (model-based hypothesis)", fontsize=10)
    fig.tight_layout()
    save(fig, "fig5_importance_perturbation")


def learning_and_pretraining() -> None:
    lc = pd.read_csv(R / "phase6" / "learning_curve_runs.csv")
    pr = pd.read_csv(R / "phase6" / "probe_summary.csv")
    g = pd.read_csv(R / "phase3" / "graph_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    for rec, c in (("myeloid", "#d62728"), ("tumour", "#2171b5")):
        d = lc[lc.receiver == rec]
        ax.scatter(d.n_train_cells, d.test_r2_fold, color=c, alpha=.5, s=18, label=f"{rec} ({d.target.iloc[0]}) per fold")
        agg = d.groupby("n_train_tumours").agg(x=("n_train_cells", "mean"), y=("test_r2_fold", "mean"))
        b_, a_ = np.polyfit(np.log10(agg.x), agg.y, 1)
        xs = np.logspace(np.log10(agg.x.min()), np.log10(agg.x.max() * 10), 50)
        ax.plot(xs, a_ + b_ * np.log10(xs), color=c, ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("training receiver cells (log)")
    ax.set_ylabel("held-out R² (fold)")
    ax.set_title("Learning curve (graph model, 1 seed); dashed = log-linear\nextrapolation to 10x data (rough)", fontsize=9)
    ax.legend(fontsize=7)
    ax = axes[1]
    rows = []
    for rec in ("myeloid", "tumour"):
        for t in TARGETS[rec][:2]:
            sel = g[(g.receiver == rec) & (g.target == t) & (g.arch == "sage")].iloc[0]
            rows.append((f"{rec}\n{t.split(':')[1]}", "supervised graph", sel.r2, sel.ci_low, sel.ci_high))
            for pb in ("neighbourhood_only", "target_masked"):
                r = pr[(pr.receiver == rec) & (pr.target == t) & (pr.probe == pb)].iloc[0]
                rows.append((f"{rec}\n{t.split(':')[1]}", {"neighbourhood_only": "probe: neighbourhood only (like-for-like)",
                                                         "target_masked": "probe: + own non-target genes (not like-for-like)"}[pb],
                             r.r2, r.ci_low, r.ci_high))
    d = pd.DataFrame(rows, columns=["t", "m", "r2", "lo", "hi"])
    ts, ms = list(dict.fromkeys(d.t)), list(dict.fromkeys(d.m))
    for k, m in enumerate(ms):
        e = d[d.m == m].set_index("t").reindex(ts)
        x = np.arange(len(ts)) + (k - 1) * 0.27
        ax.bar(x, e.r2, 0.27, yerr=[e.r2 - e.lo, e.hi - e.r2], capsize=2, label=m, color=["#d62728", "#7f7f7f", "#c7c7c7"][k])
    ax.set_xticks(range(len(ts)), ts, fontsize=8)
    ax.set_ylabel("held-out R² (95% tumour CI)")
    ax.set_title("Masked-expression pretraining: linear probe vs supervised", fontsize=10)
    ax.legend(fontsize=7)
    fig.tight_layout()
    save(fig, "fig6_learning_curve_pretraining")


def spatial_example() -> None:
    """Illustration only: one development section, myeloid cells coloured by C1Q_TAM score with CAF context."""
    import assemble as A
    path = C.CACHE / "datasets" / "mendeley_pdac_smi" / "U6-a_1.npz"
    d = A.load_section(path)
    xy, lab = d["xy"], d["label"].astype(str)
    s = d["scores"][:, list(d["score_names"]).index("myeloid:C1Q_TAM")]
    my = lab == "myeloid"
    cx, cy = np.median(xy[my], axis=0)
    win = (np.abs(xy[:, 0] - cx) < 600) & (np.abs(xy[:, 1] - cy) < 600)
    fig, ax = plt.subplots(figsize=(7.5, 7))
    for cls, col, sz in (("epithelial_ductal_tumor", "#fdd0a2", 3), ("CAF_fibroblast", "#9ecae1", 3)):
        m = win & (lab == cls)
        ax.scatter(xy[m, 0], xy[m, 1], s=sz, c=col, label=cls.replace("epithelial_ductal_tumor", "tumour").replace("CAF_fibroblast", "CAF"), lw=0)
    m = win & ~np.isin(lab, ["epithelial_ductal_tumor", "CAF_fibroblast", "myeloid"])
    ax.scatter(xy[m, 0], xy[m, 1], s=2, c="#e0e0e0", lw=0, label="other")
    m = win & my
    sc = ax.scatter(xy[m, 0], xy[m, 1], s=9, c=s[m], cmap="magma", lw=0, label="myeloid (colour = C1Q_TAM score)")
    plt.colorbar(sc, ax=ax, shrink=.7, label="myeloid C1Q_TAM program score")
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ext_x, ext_y = np.ptp(xy[win, 0]), np.ptp(xy[win, 1])
    ax.set_title(f"Section U6-a_1 ({ext_x:.0f} x {ext_y:.0f} µm window): macrophage state in its tissue neighbourhood", fontsize=10)
    ax.legend(fontsize=7, markerscale=3, loc="lower left")
    save(fig, "fig0_spatial_example")


def checkpoint_manifest() -> None:
    rows = []
    for p in sorted((C.CACHE / "checkpoints").rglob("*.pt")):
        rows.append({"path": C.rel(p), "bytes": p.stat().st_size, "sha256": C.sha256(p)})
    pd.DataFrame(rows).to_csv(R / "phase7" / "checkpoint_manifest.csv", index=False)


def main() -> None:
    log = C.start_run("phase7", "figures")
    only = sys.argv[1:]
    for f in (spatial_example, model_comparison, ablations, external, training_curves, importance_and_perturbation,
              learning_and_pretraining, checkpoint_manifest):
        if only and f.__name__ not in only:
            continue
        f()
        log.info("done: %s", f.__name__)
    log.info("done")


if __name__ == "__main__":
    main()

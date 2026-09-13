"""Runner for spatial edge-score null permutation testing and comparison.

THEORETICAL AND SCIENTIFIC CONTEXT:
-----------------------------------
This runner addresses a core limitation of cluster-marginal CCC scoring (CellPhoneDB /
Squidpy ligrec). Standard cluster-mean scoring evaluates (mean_L + mean_R) / 2 across
all cells of a cluster without ever consulting spatial coordinates or edges, making it
mathematically identical to dissociated scRNA-seq.

Here, we evaluate an edge-based spatial LR score over the actual 6-NN spatial connectivity
graph (adata.obsp['spatial_connectivities']):
  edge_score = mean_{ (i->j) in graph : type[i]==sender, type[j]==receiver } ( expr_L[i] * expr_R[j] )

We reuse the exact same permutation generator (seed=42, N=500 label permutations) from
spatial_null.py, holding the spatial graph topology and single-cell expression profiles
fixed while permuting cell-type identities.

This produces:
  1. results/spatial_edge_null.csv
  2. results/null_model_comparison.csv (with direct verdict_flip tracking vs compositional null)
  3. figures/spatial_edge_null_comparison.png / .pdf
  4. figures/null_model_flip_comparison.png / .pdf
  5. logs/spatial_edge_null.log
"""
from __future__ import annotations

import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURES_DIR, LOGS_DIR, PROCESSED_DIR, RESULTS_DIR, SEED
from spatial_edge_score import (
    compute_edge_scores_for_labels,
    get_spatial_graph_edges,
    precompute_edge_products,
)
from spatial_null import generate_label_permutations


def setup_logger(log_path: Path) -> logging.Logger:
    """Configure logger with standard format (asctime levelname message)."""
    logger = logging.getLogger("spatial_edge_null")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def save_figure(fig: plt.Figure, path: Path) -> None:
    """Save figure in both PNG and PDF formats."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_edge_null_figure(
    edge_df: pd.DataFrame,
    null_matrix: np.ndarray,
    baseline_df: pd.DataFrame,
    out_path: Path,
    top_n: int = 10,
) -> None:
    """Plot null distributions against observed edge scores for top N interactions."""
    top_df = edge_df.head(top_n).copy()
    if top_df.empty:
        return

    interaction_labels = [
        f"{row['ligand']}→{row['receptor']} ({row['sender']}→{row['receiver']})"
        for _, row in top_df.iterrows()
    ]

    baseline_lookup = {
        (r["sender"], r["receiver"], r["ligand"], r["receptor"]): idx
        for idx, r in baseline_df.iterrows()
    }

    plot_records = []
    observed_points = []

    for rank, (_, row) in enumerate(top_df.iterrows()):
        key = (row["sender"], row["receiver"], row["ligand"], row["receptor"])
        base_idx = baseline_lookup[key]
        null_vals = null_matrix[:, base_idx]
        label = interaction_labels[rank]

        for val in null_vals:
            plot_records.append({"interaction": label, "score": val})
        observed_points.append({
            "interaction": label,
            "score": row["edge_observed_score"],
            "significant": row["edge_significant"],
            "rank": rank,
        })

    dist_df = pd.DataFrame(plot_records)
    obs_df = pd.DataFrame(observed_points)

    sns.set_theme(style="whitegrid", font="sans-serif")
    fig, ax = plt.subplots(figsize=(11, 7))

    sns.boxplot(
        data=dist_df,
        y="interaction",
        x="score",
        ax=ax,
        color="#c6dbef",
        showmeans=True,
        meanprops={"marker": "o", "markerfacecolor": "#3182bd", "markeredgecolor": "#3182bd", "markersize": 5},
        flierprops={"marker": ".", "markersize": 3, "alpha": 0.3, "color": "#9ecae1"},
        boxprops={"alpha": 0.8},
        width=0.55,
    )

    sig_labeled = False
    nonsig_labeled = False
    for _, pt in obs_df.iterrows():
        marker_color = "#de2d26" if pt["significant"] else "#636363"
        edge_color = "#a50f15" if pt["significant"] else "#252525"
        marker_symbol = "*" if pt["significant"] else "X"
        lbl = None
        if pt["significant"] and not sig_labeled:
            lbl = "Observed Score (Sig p < 0.05)"
            sig_labeled = True
        elif not pt["significant"] and not nonsig_labeled:
            lbl = "Observed Score (Non-Sig)"
            nonsig_labeled = True

        ax.scatter(
            pt["score"],
            pt["interaction"],
            color=marker_color,
            edgecolors=edge_color,
            s=160,
            marker=marker_symbol,
            zorder=10,
            label=lbl,
        )

    ax.set_title(
        f"Top {top_n} Spatial Edge-Based Scores vs. Spatial-Null Distribution (N={null_matrix.shape[0]} Permutations)",
        fontsize=13,
        weight="bold",
        pad=14,
    )
    ax.set_xlabel("Edge-Based Spatial Score: Mean(expr_L(i) * expr_R(j)) over 6-NN Edges", fontsize=11)
    ax.set_ylabel("", fontsize=11)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc")

    plt.tight_layout()
    save_figure(fig, out_path)


def make_flip_comparison_figure(
    merged_df: pd.DataFrame,
    out_path: Path,
) -> None:
    """Generate a clear comparison figure highlighting verdict flips between null models."""
    df = merged_df.copy()
    df["interaction_label"] = (
        df["ligand"] + "→" + df["receptor"] + " (" + df["sender"] + "→" + df["receiver"] + ")"
    )

    # Sort by compositional observed score descending
    df = df.sort_values("compositional_observed_score", ascending=True).reset_index(drop=True)

    sns.set_theme(style="whitegrid", font="sans-serif")
    fig, axes = plt.subplots(1, 2, figsize=(15, 9), sharey=True)

    y_pos = np.arange(len(df))

    # Panel 1: Compositional Z-score vs Edge Z-score
    ax1 = axes[0]
    # Highlight flipped vs non-flipped
    colors = ["#e41a1c" if flip else "#377eb8" for flip in df["verdict_flip"]]

    ax1.scatter(
        df["compositional_z_score"],
        y_pos,
        color="#999999",
        s=80,
        label="Compositional Z-Score",
        marker="o",
        alpha=0.8,
    )
    ax1.scatter(
        df["edge_z_score"],
        y_pos,
        c=colors,
        s=100,
        label="Edge Z-Score (Red = Flipped)",
        marker="s",
        edgecolors="black",
        linewidths=0.8,
    )

    for i in range(len(df)):
        comp_z = df.loc[i, "compositional_z_score"]
        edge_z = df.loc[i, "edge_z_score"]
        flip = df.loc[i, "verdict_flip"]
        line_color = "#e41a1c" if flip else "#bbbbbb"
        line_style = "-" if flip else ":"
        line_width = 1.6 if flip else 1.0
        ax1.plot([comp_z, edge_z], [y_pos[i], y_pos[i]], color=line_color, linestyle=line_style, linewidth=line_width)

    ax1.axvline(1.96, color="#4daf4a", linestyle="--", linewidth=1.2, label="Sig Threshold (z ~ 1.96)")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(df["interaction_label"], fontsize=9)
    ax1.set_xlabel("Standardized Z-Score", fontsize=11)
    ax1.set_title("Z-Score Shift: Compositional vs. Spatial Edge", fontsize=12, weight="bold")
    ax1.legend(loc="lower right", frameon=True, facecolor="white", fontsize=9)

    # Panel 2: Verdict comparison matrix / summary
    ax2 = axes[1]
    comp_sig = df["compositional_significant"].to_numpy()
    edge_sig = df["edge_significant"].to_numpy()
    flips = df["verdict_flip"].to_numpy()

    ax2.scatter(
        np.zeros(len(df)),
        y_pos,
        c=["#4daf4a" if s else "#999999" for s in comp_sig],
        s=120,
        marker="o",
        edgecolors="black",
        label="Compositional Verdict",
    )
    ax2.scatter(
        np.ones(len(df)),
        y_pos,
        c=["#4daf4a" if s else "#999999" for s in edge_sig],
        s=120,
        marker="s",
        edgecolors="black",
        label="Edge Verdict",
    )

    for i in range(len(df)):
        if flips[i]:
            ax2.plot([0, 1], [y_pos[i], y_pos[i]], color="#e41a1c", linewidth=2.5, linestyle="-")
            ax2.text(
                1.15,
                y_pos[i],
                "FLIPPED (Sig → Non-Sig)",
                va="center",
                color="#e41a1c",
                fontweight="bold",
                fontsize=8.5,
            )
        else:
            ax2.plot([0, 1], [y_pos[i], y_pos[i]], color="#e0e0e0", linewidth=1.0, linestyle=":")

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["Compositional Null", "Spatial Edge Null"], fontsize=11, weight="bold")
    ax2.set_xlim(-0.3, 2.2)
    ax2.set_title(
        f"Verdict Flips ({flips.sum()} Flipped out of {len(df)} Interactions)",
        fontsize=12,
        weight="bold",
    )

    # Legend for verdicts
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#4daf4a", markeredgecolor="k", markersize=9, label="Significant (p < 0.05)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999", markeredgecolor="k", markersize=9, label="Non-Significant"),
        Line2D([0], [0], color="#e41a1c", lw=2.5, label="Verdict Flip"),
    ]
    ax2.legend(handles=legend_elements, loc="lower right", frameon=True, facecolor="white", fontsize=9)

    plt.tight_layout()
    save_figure(fig, out_path)


def main() -> None:
    for directory in (LOGS_DIR, RESULTS_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "spatial_edge_null.log"
    logger = setup_logger(log_path)

    start_iso = datetime.now(timezone.utc).isoformat()
    t_start = time.perf_counter()
    n_perms = 500

    logger.info("spatial_edge_null start=%s seed=%d N=%d", start_iso, SEED, n_perms)

    # 1. Load data
    h5ad_path = PROCESSED_DIR / "xenium_breast_baseline.h5ad"
    baseline_csv = RESULTS_DIR / "baseline_lr_ranked.csv"
    comp_null_csv = RESULTS_DIR / "spatial_null_corrected.csv"

    if not h5ad_path.exists():
        raise FileNotFoundError(f"Processed AnnData missing: {h5ad_path}")
    if not baseline_csv.exists():
        raise FileNotFoundError(f"Baseline CSV missing: {baseline_csv}")
    if not comp_null_csv.exists():
        raise FileNotFoundError(f"Prior spatial_null_corrected.csv missing: {comp_null_csv}")

    adata = ad.read_h5ad(h5ad_path)
    baseline_df = pd.read_csv(baseline_csv)
    comp_df = pd.read_csv(comp_null_csv)

    n_interactions = len(baseline_df)
    logger.info("Loaded AnnData (%d cells) and baseline (%d interactions).", adata.n_obs, n_interactions)

    # 2. Extract spatial edges & precompute edge products
    rows, cols = get_spatial_graph_edges(adata)
    n_total_edges = len(rows)
    logger.info("Extracted %d directed spatial edges from adata.obsp['spatial_connectivities'].", n_total_edges)

    interactions: list[tuple[str, str, str, str]] = [
        (row["sender"], row["receiver"], row["ligand"], row["receptor"])
        for _, row in baseline_df.iterrows()
    ]
    lr_pairs = list(set([(lig, rec) for _, _, lig, rec in interactions]))
    edge_products = precompute_edge_products(adata, lr_pairs, rows, cols)

    # 3. Compute observed edge scores
    orig_labels = np.asarray(adata.obs["cell_type"].astype(str).to_numpy())
    obs_scores, obs_counts = compute_edge_scores_for_labels(
        labels=orig_labels,
        row_indices=rows,
        col_indices=cols,
        edge_products=edge_products,
        interactions=interactions,
    )

    # 4. Generate permutations and evaluate null distribution
    logger.info("Generating %d label permutations (seed=%d) and evaluating edge scores...", n_perms, SEED)
    t_perm_start = time.perf_counter()
    perms = generate_label_permutations(orig_labels, n_perms=n_perms, seed=SEED)

    null_matrix = np.zeros((n_perms, n_interactions), dtype=np.float64)
    null_edge_counts = np.zeros((n_perms, n_interactions), dtype=np.int64)

    for p, p_labels in enumerate(perms):
        p_scores, p_counts = compute_edge_scores_for_labels(
            labels=p_labels,
            row_indices=rows,
            col_indices=cols,
            edge_products=edge_products,
            interactions=interactions,
        )
        # In permutations with zero qualifying edges, edge signaling potential is 0.0
        p_scores_clean = np.where(np.isnan(p_scores), 0.0, p_scores)
        null_matrix[p, :] = p_scores_clean
        null_edge_counts[p, :] = p_counts

    perm_runtime = time.perf_counter() - t_perm_start
    logger.info("Permutations completed in %.4f s (%.5f s/permutation)", perm_runtime, perm_runtime / n_perms)

    # 5. Compute statistics per interaction
    edge_records = []
    for i, (sender, receiver, ligand, receptor) in enumerate(interactions):
        obs_score = obs_scores[i]
        n_qualifying = int(obs_counts[i])

        if n_qualifying == 0:
            null_mean = np.nan
            null_std = np.nan
            empirical_p = np.nan
            z_score = np.nan
            sig = False
            status = "insufficient_edges"
        else:
            null_dist = null_matrix[:, i]
            null_mean = float(np.mean(null_dist))
            null_std = float(np.std(null_dist, ddof=1)) if n_perms > 1 else 0.0
            empirical_p = float(np.mean(null_dist >= obs_score))
            z_score = float((obs_score - null_mean) / null_std) if null_std > 1e-12 else 0.0
            sig = bool(empirical_p < 0.05)
            status = "evaluated"

        edge_records.append({
            "sender": sender,
            "receiver": receiver,
            "ligand": ligand,
            "receptor": receptor,
            "n_qualifying_edges": n_qualifying,
            "edge_observed_score": obs_score,
            "edge_null_mean": null_mean,
            "edge_null_std": null_std,
            "edge_empirical_pvalue": empirical_p,
            "edge_z_score": z_score,
            "edge_significant": sig,
        })

    edge_df = pd.DataFrame(edge_records)
    # Sort by edge_observed_score descending
    edge_df = edge_df.sort_values("edge_observed_score", ascending=False).reset_index(drop=True)

    edge_csv = RESULTS_DIR / "spatial_edge_null.csv"
    edge_df.to_csv(edge_csv, index=False)
    logger.info("Saved spatial edge null results to %s", edge_csv)

    # 6. Merge with compositional null for direct comparison
    comp_renamed = comp_df.rename(columns={
        "observed_score": "compositional_observed_score",
        "null_mean": "compositional_null_mean",
        "null_std": "compositional_null_std",
        "empirical_pvalue": "compositional_empirical_pvalue",
        "z_score": "compositional_z_score",
        "significant": "compositional_significant",
    })

    merged_df = pd.merge(
        comp_renamed,
        edge_df,
        on=["sender", "receiver", "ligand", "receptor"],
        how="inner",
    )

    merged_df["verdict_flip"] = merged_df["compositional_significant"] != merged_df["edge_significant"]
    # Sort merged table by compositional_observed_score descending (matching baseline order)
    merged_df = merged_df.sort_values("compositional_observed_score", ascending=False).reset_index(drop=True)

    comparison_csv = RESULTS_DIR / "null_model_comparison.csv"
    merged_df.to_csv(comparison_csv, index=False)
    logger.info("Saved null model comparison to %s", comparison_csv)

    n_edge_sig = int(edge_df["edge_significant"].sum())
    n_flips = int(merged_df["verdict_flip"].sum())
    logger.info(
        "Summary: %d/%d interactions significant under edge null. Verdict flips vs compositional null: %d.",
        n_edge_sig,
        n_interactions,
        n_flips,
    )

    # 7. Generate comparison figures
    edge_fig_path = FIGURES_DIR / "spatial_edge_null_comparison"
    make_edge_null_figure(
        edge_df=edge_df,
        null_matrix=null_matrix,
        baseline_df=baseline_df,
        out_path=edge_fig_path,
        top_n=10,
    )
    logger.info("Saved spatial edge null figure to %s.png / .pdf", edge_fig_path)

    flip_fig_path = FIGURES_DIR / "null_model_flip_comparison"
    make_flip_comparison_figure(
        merged_df=merged_df,
        out_path=flip_fig_path,
    )
    logger.info("Saved flip comparison figure to %s.png / .pdf", flip_fig_path)

    total_runtime = time.perf_counter() - t_start
    logger.info(
        "spatial_edge_null finished: runtime=%.2f s, seed=%d, N=%d, tested=%d, edge_significant=%d, flips=%d",
        total_runtime,
        SEED,
        n_perms,
        n_interactions,
        n_edge_sig,
        n_flips,
    )

    print("\n=== Spatial Edge-Score Null Test Complete ===")
    print(f"Total runtime: {total_runtime:.2f} s")
    print(f"Permutations: {n_perms}")
    print(f"Interactions tested: {n_interactions}")
    print(f"Edge-significant (p < 0.05): {n_edge_sig}")
    print(f"Verdict flips vs compositional null: {n_flips}")
    print(f"Edge Results CSV: {edge_csv}")
    print(f"Comparison CSV: {comparison_csv}")
    print(f"Figures: {edge_fig_path}.png, {flip_fig_path}.png (+ .pdf)")
    print(f"Log: {log_path}\n")


if __name__ == "__main__":
    main()

"""Run the spatial-null permutation test and generate corrected CCC outputs.

CITED LITERATURE AND METHODOLOGICAL CONTEXT:
--------------------------------------------
This script implements a spatial-null test to determine whether the top-ranked
ligand-receptor interactions identified in results/baseline_lr_ranked.csv represent
true spatial cell-cell communication or are artifacts of spatial autocorrelation and
global expression abundance.

The method follows the spatial-null evaluation principles established in:
  - CONCISE (Zhao et al., 2026, bioRxiv / PMC13320749): "Spatial co-expression and
    cell-cell communication inference from spatially resolved transcriptomics with CONCISE"
  - SOAAR (Khatri et al., 2026, bioRxiv): "Spatial Autocorrelation Aware Resampling
    Improves Cell-Cell Interaction Inference in Spatial Transcriptomics Data"

Variant implemented:
  Spatial Label-Permutation Null Model (preserving spatial graph topology and single-cell
  expression profiles while randomly shuffling cell-type labels across cells).
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
from spatial_null import calculate_null_statistics, run_spatial_null_permutations


def setup_logger(log_path: Path) -> logging.Logger:
    """Configure logger with prototype.log format (asctime levelname message)."""
    logger = logging.getLogger("spatial_null")
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


def make_comparison_figure(
    corrected_df: pd.DataFrame,
    null_matrix: np.ndarray,
    baseline_df: pd.DataFrame,
    out_path: Path,
    top_n: int = 10,
) -> None:
    """Plot null distributions against observed scores for top N interactions."""
    top_df = corrected_df.head(top_n).copy()
    if top_df.empty:
        return

    # Build interaction display labels
    interaction_labels = [
        f"{row['ligand']}→{row['receptor']} ({row['sender']}→{row['receiver']})"
        for _, row in top_df.iterrows()
    ]

    # Map baseline index for each row in top_df to extract matching null distribution
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
            plot_records.append({"interaction": label, "score": val, "type": "Null Distribution"})
        observed_points.append({
            "interaction": label,
            "score": row["observed_score"],
            "significant": row["significant"],
            "rank": rank,
        })

    dist_df = pd.DataFrame(plot_records)
    obs_df = pd.DataFrame(observed_points)

    sns.set_theme(style="whitegrid", font="sans-serif")
    fig, ax = plt.subplots(figsize=(11, 7))

    # Boxplot of null distribution for each interaction
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
        f"Top {top_n} Baseline Interactions vs. Spatial-Null Distribution (N={null_matrix.shape[0]} Permutations)",
        fontsize=13,
        weight="bold",
        pad=14,
    )
    ax.set_xlabel("Cell-Cell Interaction Score: (Mean Ligand + Mean Receptor) / 2", fontsize=11)
    ax.set_ylabel("", fontsize=11)

    # Custom legend
    handles, labels = ax.get_legend_handles_labels()
    # Deduplicate legend items
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc")

    plt.tight_layout()
    save_figure(fig, out_path)


def main() -> None:
    for directory in (LOGS_DIR, RESULTS_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "spatial_null.log"
    logger = setup_logger(log_path)

    start_iso = datetime.now(timezone.utc).isoformat()
    t_start = time.perf_counter()
    logger.info("spatial_null start=%s seed=%d", start_iso, SEED)

    # Load baseline artifacts
    h5ad_path = PROCESSED_DIR / "xenium_breast_baseline.h5ad"
    baseline_csv = RESULTS_DIR / "baseline_lr_ranked.csv"

    if not h5ad_path.exists():
        raise FileNotFoundError(f"Processed AnnData missing: {h5ad_path}")
    if not baseline_csv.exists():
        raise FileNotFoundError(f"Baseline CSV missing: {baseline_csv}")

    adata = ad.read_h5ad(h5ad_path)
    baseline_df = pd.read_csv(baseline_csv)
    n_interactions = len(baseline_df)
    logger.info("Loaded AnnData: %d cells, %d genes. Testing %d baseline interactions.", adata.n_obs, adata.n_vars, n_interactions)

    # Step 1: Benchmark timing on first 10 permutations
    logger.info("Timing initial benchmark of 10 permutations...")
    t0 = time.perf_counter()
    _ = run_spatial_null_permutations(adata, baseline_df, n_perms=10, seed=SEED)
    t1 = time.perf_counter()
    time_10 = t1 - t0
    rate_per_perm = time_10 / 10.0
    logger.info("10 permutations completed in %.4f s (%.5f s/permutation)", time_10, rate_per_perm)

    # Extrapolate runtime
    est_200 = rate_per_perm * 200
    est_500 = rate_per_perm * 500
    logger.info("Extrapolated runtime: N=200 -> %.2f s, N=500 -> %.2f s", est_200, est_500)

    # Decide N: choose N=500 since runtime is negligible (< 1-2 seconds)
    final_n_perms = 500
    logger.info("Timing decision: Selecting N=%d permutations (est. runtime %.2f s, within recommended 200-500 range)", final_n_perms, est_500)

    # Step 2: Run full permutation set with fixed seed
    logger.info("Executing full spatial-null permutation test (N=%d, seed=%d)...", final_n_perms, SEED)
    t_perm_start = time.perf_counter()
    null_matrix = run_spatial_null_permutations(adata, baseline_df, n_perms=final_n_perms, seed=SEED)
    t_perm_end = time.perf_counter()
    perm_runtime = t_perm_end - t_perm_start
    logger.info("Permutations completed in %.4f s", perm_runtime)

    # Save aggregated null distributions to .npz
    npz_path = RESULTS_DIR / "spatial_null_distributions.npz"
    np.savez_compressed(
        npz_path,
        null_scores=null_matrix,
        interactions=baseline_df[["sender", "receiver", "ligand", "receptor"]].to_numpy(),
        seed=SEED,
        n_perms=final_n_perms,
    )
    logger.info("Saved aggregated null distribution matrix to %s", npz_path)

    # Step 3: Compute statistics and generate corrected CSV
    corrected_df = calculate_null_statistics(baseline_df, null_matrix)
    corrected_csv = RESULTS_DIR / "spatial_null_corrected.csv"
    corrected_df.to_csv(corrected_csv, index=False)
    logger.info("Saved spatial-null corrected results to %s", corrected_csv)

    n_sig = int(corrected_df["significant"].sum())
    n_nonsig = n_interactions - n_sig
    logger.info("Results summary: %d/%d interactions significant (empirical p < 0.05), %d non-significant", n_sig, n_interactions, n_nonsig)

    # Step 4: Generate comparison figure
    fig_path = FIGURES_DIR / "spatial_null_comparison"
    make_comparison_figure(
        corrected_df=corrected_df,
        null_matrix=null_matrix,
        baseline_df=baseline_df,
        out_path=fig_path,
        top_n=10,
    )
    logger.info("Saved comparison figure to %s.png and .pdf", fig_path)

    total_runtime = time.perf_counter() - t_start
    logger.info("spatial_null finished: runtime=%.2f s, seed=%d, N=%d, tested=%d, significant=%d", total_runtime, SEED, final_n_perms, n_interactions, n_sig)

    print("\n=== Spatial-Null Permutation Test Complete ===")
    print(f"Total runtime: {total_runtime:.2f} s")
    print(f"Permutations: {final_n_perms}")
    print(f"Interactions tested: {n_interactions}")
    print(f"Significant interactions (p < 0.05): {n_sig}")
    print(f"Non-significant interactions: {n_nonsig}")
    print(f"Results CSV: {corrected_csv}")
    print(f"Figure: {fig_path}.png (+ .pdf)")
    print(f"Log: {log_path}\n")


if __name__ == "__main__":
    main()

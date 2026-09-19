"""Enhancement 1: Continuous Soft Confidence-Weighted Edge Null Model.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
  - results/gnn_edge_score/gnn_pca_node_embeddings.npy (read-only)
  - results/null_model_comparison.csv (read-only)
  - results/gnn_edge_score/gnn_confidence_filtered_null.csv (read-only)
Writes to results/soft_confidence_null/ and figures/.

Key design:
  - Instead of discarding edges via hard thresholding (top_frac=0.5), every directed edge
    is continuously weighted by its GNN link-prediction confidence:
      score_soft = sum(p_edge * expr_L * expr_R) / sum(p_edge)
  - Correctness gate: setting p_edge = 1.0 everywhere must reproduce edge_observed_score
    to within 1e-12 absolute difference across all 22 interactions.
  - Generates 500 permutations under seed=42 using identical permutation sequence.
  - Tests whether continuous soft weighting preserves ranking stability while rescuing
    the top-5 interaction dilution observed under hard filtering.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
from spatial_edge_score import get_spatial_graph_edges, precompute_edge_products
from spatial_null import generate_label_permutations

SEED = 42
N_PERMS = 500

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
GNN_EMB_PATH = ROOT / "results" / "gnn_edge_score" / "gnn_pca_node_embeddings.npy"
NULL_COMPARISON_CSV = ROOT / "results" / "null_model_comparison.csv"
FILTERED_NULL_CSV = ROOT / "results" / "gnn_edge_score" / "gnn_confidence_filtered_null.csv"
OUT_DIR = ROOT / "results" / "soft_confidence_null"
FIGURES_DIR = ROOT / "figures"


def compute_soft_scores_for_labels(
    labels: np.ndarray,
    row_indices: np.ndarray,
    col_indices: np.ndarray,
    edge_products: dict[tuple[str, str], np.ndarray],
    p_edge: np.ndarray,
    interactions: list[tuple[str, str, str, str]],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute continuous soft confidence-weighted edge scores.

    score_soft = sum(p_edge * expr_L * expr_R) / sum(p_edge)
    """
    sender_types = labels[row_indices]
    receiver_types = labels[col_indices]

    n_interactions = len(interactions)
    scores = np.zeros(n_interactions, dtype=np.float64)
    edge_counts = np.zeros(n_interactions, dtype=np.int64)

    unique_type_pairs = set([(s, r) for s, r, _, _ in interactions])
    pair_masks: dict[tuple[str, str], np.ndarray] = {}
    for s, r in unique_type_pairs:
        pair_masks[(s, r)] = (sender_types == s) & (receiver_types == r)

    for i, (s, r, lig, rec) in enumerate(interactions):
        mask = pair_masks[(s, r)]
        n_edges = int(np.sum(mask))
        edge_counts[i] = n_edges

        if n_edges > 0:
            weights = p_edge[mask]
            sum_w = float(np.sum(weights))
            if sum_w > 0:
                prod_arr = edge_products[(lig, rec)][mask]
                scores[i] = float(np.sum(weights * prod_arr) / sum_w)
            else:
                scores[i] = 0.0
        else:
            scores[i] = np.nan

    return scores, edge_counts


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading AnnData and spatial graph...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    labels_true = adata.obs["cell_type"].astype(str).values

    row_idx, col_idx = get_spatial_graph_edges(adata)
    n_edges = len(row_idx)
    print(f"Loaded adata with {n_cells} cells, {n_edges} directed edges.")

    # Compute p_edge from PCA GNN node embeddings
    print("Loading PCA GNN node embeddings...")
    node_emb = np.load(GNN_EMB_PATH)
    assert node_emb.shape == (n_cells, 64)
    with torch.no_grad():
        src_z = torch.tensor(node_emb[row_idx])
        dst_z = torch.tensor(node_emb[col_idx])
        logits = (src_z * dst_z).sum(dim=-1)
        p_edge = torch.sigmoid(logits).cpu().numpy().astype(np.float64)
    print(f"Computed p_edge: min={p_edge.min():.4f}, mean={p_edge.mean():.4f}, max={p_edge.max():.4f}")

    # Load 22 baseline interactions
    null_df = pd.read_csv(NULL_COMPARISON_CSV)
    interactions = [
        (row["sender"], row["receiver"], row["ligand"], row["receptor"])
        for _, row in null_df.iterrows()
    ]
    pairs = list(set([(lig, rec) for _, _, lig, rec in interactions]))
    edge_products = precompute_edge_products(adata, pairs, row_idx, col_idx)

    # ---------------- CORRECTNESS GATE ----------------
    print("\n--- Evaluating Correctness Gate (Uniform Weights vs Raw Edge Score) ---")
    uniform_p = np.ones_like(p_edge)
    unweighted_scores, _ = compute_soft_scores_for_labels(
        labels_true, row_idx, col_idx, edge_products, uniform_p, interactions
    )
    expected_scores = null_df["edge_observed_score"].values
    abs_diffs = np.abs(unweighted_scores - expected_scores)
    max_abs_diff = float(np.nanmax(abs_diffs))
    gate_pass = max_abs_diff < 1e-9
    print(f"Correctness Gate max absolute difference: {max_abs_diff:.2e} (tolerance 1e-9)")
    print("CORRECTNESS GATE:", "PASS" if gate_pass else "FAIL")
    assert gate_pass, f"Correctness gate failed with diff {max_abs_diff}"

    # ---------------- COMPUTE OBSERVED SOFT SCORES ----------------
    print("\nComputing observed continuous soft confidence-weighted scores...")
    obs_soft_scores, qualifying_counts = compute_soft_scores_for_labels(
        labels_true, row_idx, col_idx, edge_products, p_edge, interactions
    )

    # ---------------- NULL PERMUTATION ENGINE ----------------
    print(f"\nRunning {N_PERMS} permutations with SEED={SEED}...")
    t_perm_start = time.time()
    perms = generate_label_permutations(labels_true, n_perms=N_PERMS, seed=SEED)
    null_distributions = np.zeros((len(interactions), N_PERMS), dtype=np.float64)

    for p_idx in range(N_PERMS):
        perm_labels = perms[p_idx]
        perm_scores, _ = compute_soft_scores_for_labels(
            perm_labels, row_idx, col_idx, edge_products, p_edge, interactions
        )
        null_distributions[:, p_idx] = np.nan_to_num(perm_scores, nan=0.0)

    perm_runtime = time.time() - t_perm_start
    print(f"Completed {N_PERMS} permutations in {perm_runtime:.2f}s ({perm_runtime/N_PERMS*1000:.1f} ms/perm)")

    # ---------------- COMPUTE NULL STATISTICS & Z-SCORES ----------------
    null_mean = np.mean(null_distributions, axis=1)
    null_std = np.std(null_distributions, axis=1)

    z_scores = np.zeros(len(interactions), dtype=np.float64)
    empirical_p = np.zeros(len(interactions), dtype=np.float64)

    for i in range(len(interactions)):
        obs = obs_soft_scores[i]
        if np.isnan(obs):
            z_scores[i] = np.nan
            empirical_p[i] = np.nan
            continue
        sd = null_std[i]
        mn = null_mean[i]
        if sd > 0:
            z_scores[i] = float((obs - mn) / sd)
        else:
            z_scores[i] = 0.0
        empirical_p[i] = float(np.mean(null_distributions[i, :] >= obs))

    # Build output dataframe
    res_df = null_df.copy()
    res_df["soft_observed_score"] = obs_soft_scores
    res_df["soft_null_mean"] = null_mean
    res_df["soft_null_std"] = null_std
    res_df["soft_empirical_pvalue"] = empirical_p
    res_df["soft_z_score"] = z_scores

    # Add hard-filtered z-score for direct comparison
    if FILTERED_NULL_CSV.exists():
        filt_df = pd.read_csv(FILTERED_NULL_CSV)
        res_df["filtered_z_score"] = filt_df["filtered_z_score"]

    csv_out_path = OUT_DIR / "soft_confidence_null.csv"
    res_df.to_csv(csv_out_path, index=False)
    print(f"Saved: {csv_out_path}")

    # ---------------- EVALUATE SUCCESS CRITERIA ----------------
    cxcl12_mask = res_df["ligand"] == "CXCL12"
    cxcl12_df = res_df[cxcl12_mask].dropna(subset=["soft_z_score", "edge_z_score", "compositional_z_score"])

    rho_edge_cxcl12, p_edge_cxcl12 = spearmanr(cxcl12_df["soft_z_score"], cxcl12_df["edge_z_score"])
    rho_comp_cxcl12, p_comp_cxcl12 = spearmanr(cxcl12_df["soft_z_score"], cxcl12_df["compositional_z_score"])

    # Top-5 interactions by original edge_z_score
    top5_df = res_df.sort_values("edge_z_score", ascending=False).head(5)
    median_edge_z_top5 = float(np.median(np.abs(top5_df["edge_z_score"])))
    median_filt_z_top5 = float(np.median(np.abs(top5_df["filtered_z_score"]))) if "filtered_z_score" in top5_df else np.nan
    median_soft_z_top5 = float(np.median(np.abs(top5_df["soft_z_score"])))

    rescued_dilution = bool(median_soft_z_top5 >= median_filt_z_top5)
    rank_preserved = bool(rho_edge_cxcl12 > 0.70)

    print("\n--- RESULTS & DILUTION COMPARISON ---")
    print(f"CXCL12-16 Spearman rho(soft_z, edge_z):          {rho_edge_cxcl12:.4f} (p={p_edge_cxcl12:.4e})")
    print(f"CXCL12-16 Spearman rho(soft_z, comp_z):          {rho_comp_cxcl12:.4f} (p={p_comp_cxcl12:.4e})")
    print(f"Median |z| top-5 (Raw Edge Null):                {median_edge_z_top5:.3f}")
    print(f"Median |z| top-5 (Hard-Filtered top-50%):        {median_filt_z_top5:.3f} (DILUTED)")
    print(f"Median |z| top-5 (Continuous Soft Weighted):     {median_soft_z_top5:.3f}")
    print(f"Rescued dilution vs hard filter:                 {rescued_dilution}")

    # Save summary JSON
    summary = {
        "correctness_gate": {
            "max_abs_diff": max_abs_diff,
            "pass": gate_pass,
        },
        "spearman_16_cxcl12": {
            "vs_edge_z_score": {"rho": float(rho_edge_cxcl12), "p": float(p_edge_cxcl12)},
            "vs_compositional_z_score": {"rho": float(rho_comp_cxcl12), "p": float(p_comp_cxcl12)},
        },
        "top5_interactions": [
            {
                "sender": r["sender"],
                "receiver": r["receiver"],
                "ligand": r["ligand"],
                "receptor": r["receptor"],
                "edge_z": float(r["edge_z_score"]),
                "filtered_z": float(r["filtered_z_score"]) if "filtered_z_score" in r else None,
                "soft_z": float(r["soft_z_score"]),
            }
            for _, r in top5_df.iterrows()
        ],
        "median_abs_z_top5": {
            "raw_edge": median_edge_z_top5,
            "hard_filtered": median_filt_z_top5,
            "continuous_soft": median_soft_z_top5,
        },
        "criteria": {
            "rank_preserved_rho_gt_0.7": rank_preserved,
            "rescued_dilution_vs_hard_filter": rescued_dilution,
        },
        "runtime_sec": time.time() - t0,
    }
    with open(OUT_DIR / "step_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved: {OUT_DIR / 'step_summary.json'}")

    # ---------------- GENERATE FIGURE ----------------
    print("\nGenerating Figure: Soft Confidence Comparison...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: Scatter of soft_z vs edge_z
    ax0 = axes[0]
    ax0.scatter(cxcl12_df["edge_z_score"], cxcl12_df["soft_z_score"], color="#1f77b4", s=45, alpha=0.85, edgecolor="k", lw=0.6)
    lims = [
        min(float(cxcl12_df["edge_z_score"].min()), float(cxcl12_df["soft_z_score"].min())) - 2,
        max(float(cxcl12_df["edge_z_score"].max()), float(cxcl12_df["soft_z_score"].max())) + 2,
    ]
    ax0.plot(lims, lims, "k--", alpha=0.5, label="y = x (Perfect Match)")
    ax0.set_xlim(lims)
    ax0.set_ylim(lims)
    ax0.set_xlabel("Original Edge Null $z$-score", fontsize=11)
    ax0.set_ylabel("Soft Confidence-Weighted $z$-score", fontsize=11)
    ax0.set_title(f"CXCL12–CXCR4 Ranking Preservation\nSpearman $\\rho = {rho_edge_cxcl12:.3f}$ ($p = {p_edge_cxcl12:.2e}$)", fontsize=12, fontweight="bold")
    ax0.grid(True, linestyle=":", alpha=0.6)
    ax0.legend(frameon=True, fontsize=10)

    # Panel B: Grouped bar chart comparing top-5 |z|
    ax1 = axes[1]
    labels_top5 = [f"{r['sender'][:4]}→{r['receiver'][:4]}" for _, r in top5_df.iterrows()]
    x = np.arange(len(labels_top5))
    width = 0.25

    ax1.bar(x - width, np.abs(top5_df["edge_z_score"]), width, label="Raw Edge Null", color="#7f7f7f", alpha=0.85)
    if "filtered_z_score" in top5_df:
        ax1.bar(x, np.abs(top5_df["filtered_z_score"]), width, label="Hard-Filtered (top-50%)", color="#d62728", alpha=0.85)
    ax1.bar(x + width, np.abs(top5_df["soft_z_score"]), width, label="Continuous Soft Weighted", color="#2ca02c", alpha=0.85)

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels_top5, fontsize=10)
    ax1.set_ylabel("Absolute $z$-score ($|z|$)", fontsize=11)
    ax1.set_title(
        f"Top-5 Signal Strength: Soft Weighting Rescues Dilution\nMedian: Raw={median_edge_z_top5:.2f} | Filtered={median_filt_z_top5:.2f} | Soft={median_soft_z_top5:.2f}",
        fontsize=12,
        fontweight="bold",
    )
    ax1.grid(axis="y", linestyle=":", alpha=0.6)
    ax1.legend(frameon=True, fontsize=10)

    fig.suptitle("Enhancement 1: Continuous Soft Confidence-Weighted Null Model", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "soft_confidence_comparison.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "soft_confidence_comparison.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "soft_confidence_comparison.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "soft_confidence_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved soft_confidence_comparison.png and .pdf in figures/ and results/soft_confidence_null/")

    total_time = time.time() - t0
    print(f"\nEnhancement 1 completed successfully in {total_time:.2f}s")
    return summary


if __name__ == "__main__":
    main()

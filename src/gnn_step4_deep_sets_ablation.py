"""Enhancement 2 (Part B): GNN ablation with Deep Sets permutation-invariant spatial features.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
  - results/gnn_edge_score/gnn_pca_node_embeddings.npy (read-only)
  - results/gnn_edge_score/step1_link_prediction_results.json (read-only)
  - results/gnn_edge_score/gnn_pca_step2_summary.json (read-only)
  - results/embedding_baseline/scgpt_used_genes.json (read-only)
  - results/deep_sets_features/deep_sets_node_embeddings.npy (read-only)
  - results/null_model_comparison.csv (read-only)
Writes to results/gnn_deep_sets_ablation/ and figures/.

Design:
  - Arm A (control, PCA-50 alone): reuses precomputed GraphSAGE results (AUC = 0.9799).
  - Arm C (experimental, PCA-50 + DeepSets-50): 100-d input trained with identical
    GraphSAGE architecture, 85/15 edge split, seed=42, 1:1 negative sampling per epoch,
    100 epochs, patience=10.
  - Pre-registered gate: Arm C best val AUC >= 0.9600.
  - Downstream CCC analysis: learned communication score per LR pair, anti-confound check,
    and CD274->PDCD1 rank check.
  - Generates publication figures: side-by-side spatial GNN confidence map and learning curves.
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
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
from gnn_step1_link_pred_gate import (
    GraphSAGE,
    build_pca_baseline_features,
    make_split,
    train_one_model,
    SEED,
    HIDDEN_DIM,
    OUTPUT_DIM,
    DROPOUT,
    LR,
    MAX_EPOCHS,
    PATIENCE,
)
from spatial_edge_score import get_spatial_graph_edges

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
EMB_DIR = ROOT / "results" / "embedding_baseline"
GNN_PCA_DIR = ROOT / "results" / "gnn_edge_score"
DEEP_SETS_DIR = ROOT / "results" / "deep_sets_features"
NULL_MODEL_CSV = ROOT / "results" / "null_model_comparison.csv"
OUT_DIR = ROOT / "results" / "gnn_deep_sets_ablation"
FIGURES_DIR = ROOT / "figures"

GATE_TOLERANCE = 0.02
GATE_BENCHMARK_AUC = 0.9799305748157272


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading AnnData and spatial graph...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    row_idx, col_idx = get_spatial_graph_edges(adata)
    n_edges = len(row_idx)
    print(f"Loaded adata: {n_cells} cells, {n_edges} directed edges.")
    edge_index = torch.tensor(np.stack([row_idx, col_idx]), dtype=torch.long)

    # Identical 85/15 edge split
    train_data, val_data = make_split(edge_index, n_cells, seed=SEED)

    with open(EMB_DIR / "scgpt_used_genes.json") as f:
        used_genes = json.load(f)

    pca_emb = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)
    ds_emb_path = DEEP_SETS_DIR / "deep_sets_node_embeddings.npy"
    if not ds_emb_path.exists():
        raise FileNotFoundError(f"Missing {ds_emb_path}")
    ds_emb = np.load(ds_emb_path)
    assert ds_emb.shape == (n_cells, 50)
    print(f"PCA-50: {pca_emb.shape}, DeepSets-50: {ds_emb.shape}")

    # ---------------- ARM A (Control) ----------------
    print("\n--- ARM A (Control: PCA-50 Baseline) ---")
    pca_node_emb = np.load(GNN_PCA_DIR / "gnn_pca_node_embeddings.npy")
    arm_a_best_val_auc = GATE_BENCHMARK_AUC

    # Reproduce Arm A history for learning curve comparison
    x_pca = torch.tensor(pca_emb, dtype=torch.float)
    _, history_a, _, _ = train_one_model(
        x=x_pca,
        train_edge_index=train_data.edge_index,
        train_pos_edge_label_index=train_data.edge_label_index,
        val_edge_label_index=val_data.edge_label_index,
        val_edge_label=val_data.edge_label,
        num_nodes=n_cells,
        tag="arm-a-pca",
    )
    arm_a_auc_std = float(np.std(history_a["val_auc"]))

    # ---------------- ARM C (Experimental: PCA-50 + DeepSets-50) ----------------
    print("\n--- ARM C (Experimental: PCA-50 + DeepSets-50 Concatenation) ---")
    arm_c_input_feats = np.concatenate([pca_emb, ds_emb], axis=1).astype(np.float32)
    assert arm_c_input_feats.shape == (n_cells, 100)

    x_c = torch.tensor(arm_c_input_feats, dtype=torch.float)
    t_train_c_start = time.time()
    model_c, history_c, arm_c_best_val_auc, arm_c_best_epoch = train_one_model(
        x=x_c,
        train_edge_index=train_data.edge_index,
        train_pos_edge_label_index=train_data.edge_label_index,
        val_edge_label_index=val_data.edge_label_index,
        val_edge_label=val_data.edge_label,
        num_nodes=n_cells,
        tag="arm-c-deep-sets",
    )
    t_train_c = time.time() - t_train_c_start
    print(f"Arm C training completed in {t_train_c:.2f}s: best_val_auc={arm_c_best_val_auc:.6f} @ epoch {arm_c_best_epoch}")

    # Encode full graph using Arm C model
    model_c.eval()
    with torch.no_grad():
        z_c_full = model_c.encode(x_c, edge_index)
    arm_c_node_embeddings = z_c_full.cpu().numpy().astype(np.float32)
    np.save(OUT_DIR / "gnn_deep_sets_node_embeddings.npy", arm_c_node_embeddings)
    print(f"Saved: {OUT_DIR / 'gnn_deep_sets_node_embeddings.npy'}")

    # ---------------- PRE-REGISTERED GATE EVALUATION ----------------
    gate_threshold = GATE_BENCHMARK_AUC - GATE_TOLERANCE  # 0.9599
    gate_pass = bool(arm_c_best_val_auc >= gate_threshold)
    auc_diff = float(arm_c_best_val_auc - arm_a_best_val_auc)
    stronger_claim_threshold = arm_a_best_val_auc + arm_a_auc_std
    stronger_claim_supported = bool(arm_c_best_val_auc > stronger_claim_threshold)

    print("\n--- GATE EVALUATION ---")
    print(f"Arm A best val AUC:                {arm_a_best_val_auc:.6f}")
    print(f"Arm C best val AUC (Deep Sets):    {arm_c_best_val_auc:.6f}")
    print(f"Difference (Arm C - Arm A):        {auc_diff:+.6f}")
    print(f"Gate threshold (Arm A - 0.02):     {gate_threshold:.6f}")
    print(f"GATE PASS ('not worse'):           {gate_pass}")
    print(f"STRONGER CLAIM ('Deep Sets helps'):{stronger_claim_supported}")

    # ---------------- DOWNSTREAM CCC ANALYSIS (Arm C) ----------------
    print("\n--- Downstream CCC Analysis (Arm C) ---")
    with torch.no_grad():
        src_z = torch.tensor(arm_c_node_embeddings[row_idx])
        dst_z = torch.tensor(arm_c_node_embeddings[col_idx])
        logits = (src_z * dst_z).sum(dim=-1)
        learned_edge_scores_c = torch.sigmoid(logits).numpy()

    labels = adata.obs["cell_type"].astype(str).values
    sender_types = labels[row_idx]
    receiver_types = labels[col_idx]

    null_df = pd.read_csv(NULL_MODEL_CSV)
    unique_pairs = null_df[["sender", "receiver"]].drop_duplicates().values.tolist()

    records = []
    for s, r in unique_pairs:
        mask = (sender_types == s) & (receiver_types == r)
        n_sub_edges = int(np.sum(mask))
        score = float(np.mean(learned_edge_scores_c[mask])) if n_sub_edges > 0 else np.nan
        records.append({
            "sender": s,
            "receiver": r,
            "learned_communication_score": score,
            "n_qualifying_edges_learned": n_sub_edges,
        })
    learned_df_c = pd.DataFrame(records)
    merged_c = null_df.merge(learned_df_c, on=["sender", "receiver"], how="left")
    merged_c["learned_rank"] = merged_c["learned_communication_score"].rank(ascending=False, method="min").astype(int)

    comm_csv_path = OUT_DIR / "deep_sets_communication_scores.csv"
    merged_c.to_csv(comm_csv_path, index=False)
    print(f"Saved: {comm_csv_path}")

    # Spearman checks
    valid_all = merged_c.dropna(subset=["learned_communication_score", "edge_z_score", "compositional_z_score"])
    rho_edge_all, p_edge_all = spearmanr(valid_all["learned_communication_score"], valid_all["edge_z_score"])
    rho_comp_all, p_comp_all = spearmanr(valid_all["learned_communication_score"], valid_all["compositional_z_score"])

    cxcl12 = merged_c[merged_c["ligand"] == "CXCL12"].dropna(subset=["learned_communication_score", "edge_z_score", "compositional_z_score"])
    rho_edge_cxcl12, p_edge_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["edge_z_score"])
    rho_comp_cxcl12, p_comp_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["compositional_z_score"])

    # CD274 rank check
    cd274_rows = merged_c[(merged_c["ligand"] == "CD274") & (merged_c["receptor"] == "PDCD1")]
    cd274_ranks = sorted(cd274_rows["learned_rank"].tolist())
    cd274_all_bottom_half = bool((cd274_rows["learned_rank"] > 11.0).all())

    print(f"[All 22 rows] Spearman(learned, edge_z_score)        = {rho_edge_all:.4f} (p={p_edge_all:.4g})")
    print(f"[16 CXCL12]   Spearman(learned, edge_z_score)        = {rho_edge_cxcl12:.4f} (p={p_edge_cxcl12:.4g})")
    print(f"CD274->PDCD1 ranks: {cd274_ranks} (all bottom half: {cd274_all_bottom_half})")

    confound_summary = {
        "spearman_all_22_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_all), "p": float(p_edge_all)},
            "vs_compositional_z_score": {"rho": float(rho_comp_all), "p": float(p_comp_all)},
        },
        "spearman_16_cxcl12_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_cxcl12), "p": float(p_edge_cxcl12)},
            "vs_compositional_z_score": {"rho": float(rho_comp_cxcl12), "p": float(p_comp_cxcl12)},
        },
        "cd274_pdcd1_ranks": cd274_ranks,
        "cd274_all_bottom_half": cd274_all_bottom_half,
    }
    with open(OUT_DIR / "confound_check.json", "w") as f:
        json.dump(confound_summary, f, indent=2)

    step_results = {
        "arm_a_control": {
            "feature_set": "PCA-50",
            "best_val_auc": arm_a_best_val_auc,
            "best_epoch": 100,
        },
        "arm_c_deep_sets": {
            "feature_set": "PCA-50 + DeepSets-50",
            "best_val_auc": float(arm_c_best_val_auc),
            "best_epoch": int(arm_c_best_epoch),
            "runtime_sec": t_train_c,
            "history": history_c,
        },
        "gate_evaluation": {
            "threshold": gate_threshold,
            "auc_difference": auc_diff,
            "gate_pass": gate_pass,
            "stronger_claim_supported": stronger_claim_supported,
        },
        "runtime_total_sec": time.time() - t0,
    }
    with open(OUT_DIR / "step_results.json", "w") as f:
        json.dump(step_results, f, indent=2)
    print(f"Saved: {OUT_DIR / 'step_results.json'}")

    # ---------------- GENERATE FIGURES ----------------
    print("\n--- Generating Publication Figures ---")
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    # Figure 1: Spatial GNN Confidence Map (Arm A vs Arm C)
    with torch.no_grad():
        logits_a = (torch.tensor(pca_node_emb[row_idx]) * torch.tensor(pca_node_emb[col_idx])).sum(dim=-1)
        p_a = torch.sigmoid(logits_a).numpy()
        logits_c = (torch.tensor(arm_c_node_embeddings[row_idx]) * torch.tensor(arm_c_node_embeddings[col_idx])).sum(dim=-1)
        p_c = torch.sigmoid(logits_c).numpy()

    vmin = min(float(p_a.min()), float(p_c.min()))
    vmax = max(float(p_a.max()), float(p_c.max()))
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    segments = np.stack([coords[row_idx], coords[col_idx]], axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), sharex=True, sharey=True)
    panels = [
        (axes[0], p_a, f"Arm A: PCA-50 Baseline (Val AUC = {arm_a_best_val_auc:.4f})"),
        (axes[1], p_c, f"Arm C: PCA-50 + DeepSets-50 (Val AUC = {arm_c_best_val_auc:.4f})"),
    ]
    for ax, p_vals, title in panels:
        ax.scatter(coords[:, 0], coords[:, 1], s=1.5, c="#c8c8c8", alpha=0.35, rasterized=True)
        lc = LineCollection(segments, cmap="plasma", norm=norm, array=p_vals, linewidths=0.6, alpha=0.45, rasterized=True)
        ax.add_collection(lc)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Spatial X (µm)", fontsize=10)
        ax.set_ylabel("Spatial Y (µm)", fontsize=10)
        ax.invert_yaxis()

    cbar_ax = fig.add_axes([0.92, 0.18, 0.015, 0.64])
    sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("GNN Edge Confidence ($p_{\\mathrm{edge}} = \\sigma(z_i \\cdot z_j)$)", fontsize=11)
    fig.suptitle("Spatial Graph Edge Confidence: Arm A (PCA-50) vs Arm C (PCA-50 + DeepSets-50)", fontsize=14, fontweight="bold", y=0.98)

    fig.savefig(FIGURES_DIR / "deep_sets_spatial_gnn_confidence.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "deep_sets_spatial_gnn_confidence.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_spatial_gnn_confidence.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_spatial_gnn_confidence.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved deep_sets_spatial_gnn_confidence.png/.pdf")

    # Figure 2: GNN Accuracy Curves
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(history_a["epoch"], history_a["loss"], "--", color="#1f77b4", lw=2, label="Arm A (PCA-50)")
    axes[0].plot(history_c["epoch"], history_c["loss"], "-", color="#9467bd", lw=2, label="Arm C (PCA+DeepSets)")
    axes[0].set_title("Training Loss (BCE)", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("Epoch", fontsize=10)
    axes[0].set_ylabel("Loss", fontsize=10)
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(frameon=True, fontsize=9)

    axes[1].plot(history_a["epoch"], history_a["val_auc"], "--", color="#1f77b4", lw=2, label=f"Arm A (PCA-50, {arm_a_best_val_auc:.4f})")
    axes[1].plot(history_c["epoch"], history_c["val_auc"], "-", color="#9467bd", lw=2, label=f"Arm C (PCA+DeepSets, {arm_c_best_val_auc:.4f})")
    axes[1].axhline(gate_threshold, color="#2ca02c", linestyle=":", lw=1.8, label=f"Gate Threshold ({gate_threshold:.2f})")
    axes[1].set_title("Validation Link-Prediction AUC", fontsize=11, fontweight="bold")
    axes[1].set_xlabel("Epoch", fontsize=10)
    axes[1].set_ylabel("Validation AUC", fontsize=10)
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(frameon=True, fontsize=9)

    fig.suptitle("GraphSAGE Link Prediction: Arm A (PCA-50) vs Arm C (Deep Sets)", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "deep_sets_gnn_accuracy.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "deep_sets_gnn_accuracy.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_gnn_accuracy.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_gnn_accuracy.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved deep_sets_gnn_accuracy.png/.pdf")

    total_time = time.time() - t0
    print(f"\nEnhancement 2 (Part B) completed in {total_time:.2f}s")
    return step_results


if __name__ == "__main__":
    main()

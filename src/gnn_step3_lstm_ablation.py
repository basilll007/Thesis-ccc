"""Task 2: GNN ablation — with vs. without LSTM features.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
  - results/gnn_edge_score/step1_link_prediction_results.json (read-only)
  - results/gnn_edge_score/gnn_pca_node_embeddings.npy (read-only)
  - results/gnn_edge_score/gnn_pca_step2_summary.json (read-only)
  - results/embedding_baseline/scgpt_used_genes.json (read-only)
  - results/lstm_features/lstm_node_embeddings.npy (read-only)
  - results/null_model_comparison.csv (read-only)
Writes only to results/gnn_lstm_ablation/.

Design:
  - Arm A (control, 'without'): Re-load precomputed PCA-50 GraphSAGE results
    from results/gnn_edge_score/step1_link_prediction_results.json and
    gnn_pca_node_embeddings.npy. Re-executes the deterministic 100-epoch history
    to compute per-epoch AUC fluctuation std without modifying saved embeddings.
  - Arm B ('with LSTM'): Concatenate PCA-50 (50-d) + LSTM-50 (50-d) -> 100-d node features.
    Train identical GraphSAGE architecture/hyperparameters/split/seed as
    gnn_step1_link_pred_gate.py (RandomLinkSplit seed=42, 85/15 split, 1:1 negative sampling,
    100 epochs, patience=10).
  - GATE (pre-registered): Arm B best val AUC >= 0.9799 - 0.02 (>= 0.9599).
    Stronger claim: Arm B AUC > Arm A AUC + std(Arm A per-epoch AUC fluctuation).
  - Downstream steps:
    (a) learned communication score per LR pair via decoder sigmoid(dot product)
    (b) anti-confound Spearman check vs edge_z_score and compositional_z_score
    (c) CD274->PDCD1 rank check (bottom half check)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
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
LSTM_DIR = ROOT / "results" / "lstm_features"
NULL_MODEL_CSV = ROOT / "results" / "null_model_comparison.csv"
OUT_DIR = ROOT / "results" / "gnn_lstm_ablation"

GATE_TOLERANCE = 0.02
GATE_BENCHMARK_AUC = 0.9799305748157272  # from Step 2 reproduction / 0.9799


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading AnnData and spatial graph...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs

    row_idx, col_idx = get_spatial_graph_edges(adata)
    n_edges = len(row_idx)
    print(f"Loaded adata with {n_cells} cells, {n_edges} directed spatial edges.")
    edge_index = torch.tensor(np.stack([row_idx, col_idx]), dtype=torch.long)

    # Identical 85/15 edge split
    train_data, val_data = make_split(edge_index, n_cells, seed=SEED)
    print(f"Train message passing edges: {train_data.edge_index.size(1)}, "
          f"Val supervision edges: {val_data.edge_label_index.size(1)}")

    # Load gene list
    with open(EMB_DIR / "scgpt_used_genes.json") as f:
        used_genes = json.load(f)

    # Reconstruct PCA-50 features
    pca_emb = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)
    print(f"PCA-50 baseline shape: {pca_emb.shape}")

    # Load LSTM features
    lstm_emb_path = LSTM_DIR / "lstm_node_embeddings.npy"
    if not lstm_emb_path.exists():
        raise FileNotFoundError(f"Missing {lstm_emb_path}. Run Task 1 first.")
    lstm_emb = np.load(lstm_emb_path)
    print(f"LSTM-50 features shape: {lstm_emb.shape}")
    assert lstm_emb.shape == (n_cells, 50)

    # ---------------- ARM A: Re-load PCA-50 control ----------------
    print("\n--- ARM A (Control: PCA-50 features without LSTM) ---")
    pca_emb_saved_path = GNN_PCA_DIR / "gnn_pca_node_embeddings.npy"
    if not pca_emb_saved_path.exists():
        raise FileNotFoundError(f"Missing {pca_emb_saved_path}")
    arm_a_node_embeddings = np.load(pca_emb_saved_path)
    print(f"Loaded existing Arm A node embeddings: {arm_a_node_embeddings.shape}")

    # Load existing Step 1 / Step 2 summary
    step1_json_path = GNN_PCA_DIR / "step1_link_prediction_results.json"
    step2_json_path = GNN_PCA_DIR / "gnn_pca_step2_summary.json"
    with open(step1_json_path) as f:
        step1_res = json.load(f)
    with open(step2_json_path) as f:
        step2_res = json.load(f)

    arm_a_best_val_auc = float(step2_res["step1_reproduction"]["this_run_best_val_auc"])
    arm_a_best_epoch = int(step2_res["step1_reproduction"]["best_epoch"])
    print(f"Arm A best val AUC (from Step 2): {arm_a_best_val_auc:.6f} @ epoch {arm_a_best_epoch}")

    # Reproduce Arm A's deterministic training history to compute per-epoch AUC fluctuation std
    print("Computing Arm A per-epoch AUC trajectory for std fluctuation benchmark...")
    x_pca = torch.tensor(pca_emb, dtype=torch.float)
    _, history_a, _, _ = train_one_model(
        x=x_pca,
        train_edge_index=train_data.edge_index,
        train_pos_edge_label_index=train_data.edge_label_index,
        val_edge_label_index=val_data.edge_label_index,
        val_edge_label=val_data.edge_label,
        num_nodes=n_cells,
        tag="arm-a-history",
    )
    arm_a_auc_curve = history_a["val_auc"]
    arm_a_auc_std = float(np.std(arm_a_auc_curve))
    print(f"Arm A per-epoch AUC fluctuation std: {arm_a_auc_std:.6f} (across {len(arm_a_auc_curve)} epochs)")

    # ---------------- ARM B: Train with LSTM features ----------------
    print("\n--- ARM B (Experimental: PCA-50 + LSTM-50 concatenated features) ---")
    arm_b_input_feats = np.concatenate([pca_emb, lstm_emb], axis=1).astype(np.float32)
    print(f"Arm B input features shape: {arm_b_input_feats.shape}")
    assert arm_b_input_feats.shape == (n_cells, 100)

    x_b = torch.tensor(arm_b_input_feats, dtype=torch.float)
    t_train_b_start = time.time()
    model_b, history_b, arm_b_best_val_auc, arm_b_best_epoch = train_one_model(
        x=x_b,
        train_edge_index=train_data.edge_index,
        train_pos_edge_label_index=train_data.edge_label_index,
        val_edge_label_index=val_data.edge_label_index,
        val_edge_label=val_data.edge_label,
        num_nodes=n_cells,
        tag="arm-b-lstm",
    )
    t_train_b = time.time() - t_train_b_start
    print(f"Arm B training finished in {t_train_b:.2f}s: best_val_auc={arm_b_best_val_auc:.6f} @ epoch {arm_b_best_epoch}")

    # Encode full graph using Arm B trained model
    model_b.eval()
    with torch.no_grad():
        z_b_full = model_b.encode(x_b, edge_index)
    arm_b_node_embeddings = z_b_full.cpu().numpy().astype(np.float32)
    print(f"Arm B full-graph node embeddings shape: {arm_b_node_embeddings.shape}")
    np.save(OUT_DIR / "gnn_lstm_node_embeddings.npy", arm_b_node_embeddings)
    print(f"Saved: {OUT_DIR / 'gnn_lstm_node_embeddings.npy'}")

    # ---------------- PRE-REGISTERED GATE EVALUATION ----------------
    gate_threshold = GATE_BENCHMARK_AUC - GATE_TOLERANCE  # 0.97993 - 0.02 = 0.95993
    gate_pass = bool(arm_b_best_val_auc >= gate_threshold)
    auc_diff = float(arm_b_best_val_auc - arm_a_best_val_auc)
    stronger_claim_threshold = arm_a_best_val_auc + arm_a_auc_std
    stronger_claim_supported = bool(arm_b_best_val_auc > stronger_claim_threshold)

    print("\n--- GATE EVALUATION ---")
    print(f"Arm A best val AUC:                {arm_a_best_val_auc:.6f}")
    print(f"Arm B best val AUC:                {arm_b_best_val_auc:.6f}")
    print(f"Difference (Arm B - Arm A):        {auc_diff:+.6f}")
    print(f"Gate threshold (Arm A - 0.02):     {gate_threshold:.6f}")
    print(f"GATE PASS ('not worse'):           {gate_pass}")
    print(f"Arm A per-epoch AUC std:           {arm_a_auc_std:.6f}")
    print(f"Stronger claim threshold:          {stronger_claim_threshold:.6f}")
    print(f"STRONGER CLAIM ('LSTM helps'):     {stronger_claim_supported}")

    # ---------------- DOWNSTREAM CCC ANALYSIS FOR ARM B ----------------
    print("\n--- DOWNSTREAM CCC ANALYSIS (Arm B) ---")
    with torch.no_grad():
        src_z = torch.tensor(arm_b_node_embeddings[row_idx])
        dst_z = torch.tensor(arm_b_node_embeddings[col_idx])
        logits = (src_z * dst_z).sum(dim=-1)
        learned_edge_scores_b = torch.sigmoid(logits).numpy()

    labels = adata.obs["cell_type"].astype(str).values
    sender_types = labels[row_idx]
    receiver_types = labels[col_idx]

    null_df = pd.read_csv(NULL_MODEL_CSV)
    unique_pairs = null_df[["sender", "receiver"]].drop_duplicates().values.tolist()

    records = []
    for s, r in unique_pairs:
        mask = (sender_types == s) & (receiver_types == r)
        n_sub_edges = int(np.sum(mask))
        score = float(np.mean(learned_edge_scores_b[mask])) if n_sub_edges > 0 else np.nan
        records.append({
            "sender": s,
            "receiver": r,
            "learned_communication_score": score,
            "n_qualifying_edges_learned": n_sub_edges,
        })
    learned_df_b = pd.DataFrame(records)

    merged_b = null_df.merge(learned_df_b, on=["sender", "receiver"], how="left")
    merged_b["learned_rank"] = merged_b["learned_communication_score"].rank(
        ascending=False, method="min"
    ).astype(int)

    comm_csv_path = OUT_DIR / "lstm_communication_scores.csv"
    merged_b.to_csv(comm_csv_path, index=False)
    print(f"Saved: {comm_csv_path}")

    # --- (b) Anti-confound Spearman check ---
    valid_all = merged_b.dropna(subset=["learned_communication_score", "edge_z_score", "compositional_z_score"])
    rho_edge_all, p_edge_all = spearmanr(valid_all["learned_communication_score"], valid_all["edge_z_score"])
    rho_comp_all, p_comp_all = spearmanr(valid_all["learned_communication_score"], valid_all["compositional_z_score"])

    cxcl12 = merged_b[merged_b["ligand"] == "CXCL12"].dropna(
        subset=["learned_communication_score", "edge_z_score", "compositional_z_score"]
    )
    rho_edge_cxcl12, p_edge_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["edge_z_score"])
    rho_comp_cxcl12, p_comp_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["compositional_z_score"])

    def meaningfully_weaker(rho_edge, rho_comp, tol=0.05):
        return bool(rho_edge < (rho_comp - tol))

    confound_all_weaker = meaningfully_weaker(rho_edge_all, rho_comp_all)
    confound_cxcl12_weaker = meaningfully_weaker(rho_edge_cxcl12, rho_comp_cxcl12)

    print(f"[All 22 rows] Spearman(learned, edge_z_score)        = {rho_edge_all:.4f} (p={p_edge_all:.4g})")
    print(f"[All 22 rows] Spearman(learned, compositional_z_score) = {rho_comp_all:.4f} (p={p_comp_all:.4g})")
    print(f"[16 CXCL12]   Spearman(learned, edge_z_score)        = {rho_edge_cxcl12:.4f} (p={p_edge_cxcl12:.4g})")
    print(f"[16 CXCL12]   Spearman(learned, compositional_z_score) = {rho_comp_cxcl12:.4f} (p={p_comp_cxcl12:.4g})")

    # --- (c) CD274->PDCD1 rank check ---
    cd274_rows = merged_b[(merged_b["ligand"] == "CD274") & (merged_b["receptor"] == "PDCD1")]
    cd274_ranks = sorted(cd274_rows["learned_rank"].tolist())
    n_ranked = len(merged_b)
    bottom_half_threshold = n_ranked / 2.0
    cd274_all_bottom_half = bool((cd274_rows["learned_rank"] > bottom_half_threshold).all())

    print(f"CD274->PDCD1 ranks among {n_ranked} rows: {cd274_ranks}")
    print(f"All CD274->PDCD1 in bottom half (> {bottom_half_threshold}): {cd274_all_bottom_half}")

    # Save confound check JSON
    confound_summary = {
        "spearman_all_22_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_all), "p": float(p_edge_all)},
            "vs_compositional_z_score": {"rho": float(rho_comp_all), "p": float(p_comp_all)},
            "edge_meaningfully_weaker": confound_all_weaker,
        },
        "spearman_16_cxcl12_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_cxcl12), "p": float(p_edge_cxcl12)},
            "vs_compositional_z_score": {"rho": float(rho_comp_cxcl12), "p": float(p_comp_cxcl12)},
            "edge_meaningfully_weaker": confound_cxcl12_weaker,
        },
        "cd274_pdcd1_rank_check": {
            "ranks": cd274_ranks,
            "total_ranked": n_ranked,
            "bottom_half_threshold": bottom_half_threshold,
            "all_bottom_half": cd274_all_bottom_half,
            "scores": {
                f"{r.sender}->{r.receiver}": float(r.learned_communication_score)
                for r in cd274_rows.itertuples()
            },
        },
    }
    with open(OUT_DIR / "confound_check.json", "w") as f:
        json.dump(confound_summary, f, indent=2)
    print(f"Saved: {OUT_DIR / 'confound_check.json'}")

    # Save overall step_results.json
    step_results = {
        "arm_a_control": {
            "feature_set": "PCA-50",
            "feature_dim": 50,
            "best_val_auc": arm_a_best_val_auc,
            "best_epoch": arm_a_best_epoch,
            "auc_fluctuation_std": arm_a_auc_std,
            "history": {
                "epoch": history_a["epoch"],
                "loss": history_a["loss"],
                "val_auc": history_a["val_auc"],
            },
        },
        "arm_b_with_lstm": {
            "feature_set": "PCA-50 + LSTM-50",
            "feature_dim": 100,
            "best_val_auc": float(arm_b_best_val_auc),
            "best_epoch": int(arm_b_best_epoch),
            "runtime_sec": t_train_b,
            "history": {
                "epoch": history_b["epoch"],
                "loss": history_b["loss"],
                "val_auc": history_b["val_auc"],
            },
        },
        "gate_evaluation": {
            "tolerance": GATE_TOLERANCE,
            "threshold": gate_threshold,
            "auc_difference": auc_diff,
            "gate_pass": gate_pass,
            "stronger_claim_threshold": stronger_claim_threshold,
            "stronger_claim_supported": stronger_claim_supported,
        },
        "runtime_total_sec": time.time() - t0,
    }
    with open(OUT_DIR / "step_results.json", "w") as f:
        json.dump(step_results, f, indent=2)
    print(f"Saved: {OUT_DIR / 'step_results.json'}")

    total_time = time.time() - t0
    print(f"\nTask 2 completed successfully in {total_time:.2f}s")
    return step_results


if __name__ == "__main__":
    main()

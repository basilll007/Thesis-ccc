"""GNN Step 2 continuation: PCA-features branch only (scGPT branch failed its
Step 1 gate last round and is closed).

Reproduces (does not retune) the exact GraphSAGE / split / training procedure
from src/gnn_step1_link_pred_gate.py using the PCA-50 node features that
already passed the val-AUC gate (0.9786), then uses the trained decoder's
sigmoid(dot product) as a learned edge score over the full 42,978-edge
spatial graph, aggregates it per sender-receiver cell-type pair, and compares
it against results/null_model_comparison.csv (edge_z_score and
compositional_z_score) including the mandatory anti-confound check.

Standalone script for the `scgpt-zeroshot` conda env. Reads
data/processed/xenium_breast_baseline.h5ad and
results/embedding_baseline/scgpt_used_genes.json (both read-only). Imports
(does not modify) src/gnn_step1_link_pred_gate.py and src/spatial_edge_score.py.
Writes only to results/gnn_edge_score/.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import anndata as ad
import torch
from scipy.stats import spearmanr

sys.path.insert(0, r"F:\Thesis\src")
from gnn_step1_link_pred_gate import (  # noqa: E402
    GraphSAGE, build_pca_baseline_features, make_split, train_one_model,
    SEED, HIDDEN_DIM, OUTPUT_DIM, DROPOUT, LR, MAX_EPOCHS, PATIENCE,
)
from spatial_edge_score import get_spatial_graph_edges  # noqa: E402

H5AD_PATH = r"F:\Thesis\data\processed\xenium_breast_baseline.h5ad"
EMB_DIR = r"F:\Thesis\results\embedding_baseline"
OUT_DIR = r"F:\Thesis\results\gnn_edge_score"
NULL_MODEL_CSV = r"F:\Thesis\results\null_model_comparison.csv"
PRIOR_RESULTS_JSON = os.path.join(OUT_DIR, "step1_link_prediction_results.json")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---------------- STEP 1: reproduce the PCA-features model ----------------
    t0 = time.time()

    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    print("Loaded adata:", adata.shape)

    with open(os.path.join(EMB_DIR, "scgpt_used_genes.json")) as f:
        used_genes = json.load(f)
    print("used genes:", len(used_genes))

    pca_emb = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)
    print("PCA baseline features:", pca_emb.shape)

    row_idx, col_idx = get_spatial_graph_edges(adata)
    print("n_directed_edges:", len(row_idx))
    edge_index = torch.tensor(np.stack([row_idx, col_idx]), dtype=torch.long)

    train_data, val_data = make_split(edge_index, n_cells, seed=SEED)
    print("train edges (message passing):", train_data.edge_index.size(1))
    print("val supervision edges (pos+neg):", val_data.edge_label_index.size(1),
          "| positives:", int(val_data.edge_label.sum().item()))

    x = torch.tensor(pca_emb, dtype=torch.float)
    model, history, best_val_auc, best_epoch = train_one_model(
        x=x,
        train_edge_index=train_data.edge_index,
        train_pos_edge_label_index=train_data.edge_label_index,
        val_edge_label_index=val_data.edge_label_index,
        val_edge_label=val_data.edge_label,
        num_nodes=n_cells,
        tag="pca-reproduce",
    )
    print(f"Reproduced PCA-features model: best_val_auc={best_val_auc:.6f} @ epoch {best_epoch}, "
          f"ran {len(history['epoch'])} epochs")

    with open(PRIOR_RESULTS_JSON) as f:
        prior = json.load(f)
    prior_auc = prior["pca_features"]["best_val_auc"]
    auc_diff = abs(best_val_auc - prior_auc)
    print(f"Prior run's PCA best_val_auc: {prior_auc:.6f} | this run: {best_val_auc:.6f} | "
          f"abs diff: {auc_diff:.6f}")
    reproduction_close = auc_diff < 0.01
    print("Reproduction close to prior run (abs diff < 0.01):", reproduction_close)

    # Re-encode ALL 7,163 nodes on the FULL real graph (not just the train-split
    # subgraph) -- Step 2 needs the learned score over the complete existing graph.
    model.eval()
    with torch.no_grad():
        z_full = model.encode(x, edge_index)
    node_embeddings = z_full.numpy()
    print("Full-graph node embeddings shape:", node_embeddings.shape)

    emb_path = os.path.join(OUT_DIR, "gnn_pca_node_embeddings.npy")
    np.save(emb_path, node_embeddings)
    print("Saved:", emb_path)

    t1 = time.time()
    step1_runtime = t1 - t0
    print(f"Step 1 (reproduce + save embeddings) runtime: {step1_runtime:.2f}s")

    # ---------------- STEP 2: learned communication score vs ground truth ----------------
    t2 = time.time()

    with torch.no_grad():
        src_z = torch.tensor(node_embeddings[row_idx])
        dst_z = torch.tensor(node_embeddings[col_idx])
        logits = (src_z * dst_z).sum(dim=-1)
        learned_edge_scores = torch.sigmoid(logits).numpy()

    labels = adata.obs["cell_type"].astype(str).values
    sender_types = labels[row_idx]
    receiver_types = labels[col_idx]

    null_df = pd.read_csv(NULL_MODEL_CSV)
    unique_pairs = null_df[["sender", "receiver"]].drop_duplicates().values.tolist()

    records = []
    for s, r in unique_pairs:
        mask = (sender_types == s) & (receiver_types == r)
        n_edges = int(np.sum(mask))
        score = float(np.mean(learned_edge_scores[mask])) if n_edges > 0 else np.nan
        records.append({
            "sender": s,
            "receiver": r,
            "learned_communication_score": score,
            "n_qualifying_edges_learned": n_edges,
        })
    learned_df = pd.DataFrame(records)

    merged = null_df.merge(learned_df, on=["sender", "receiver"], how="left")
    out_csv = os.path.join(OUT_DIR, "gnn_pca_communication_scores.csv")
    merged.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

    # --- Spearman correlations: all 22 rows ---
    valid_all = merged.dropna(subset=["learned_communication_score", "edge_z_score", "compositional_z_score"])
    rho_edge_all, p_edge_all = spearmanr(valid_all["learned_communication_score"], valid_all["edge_z_score"])
    rho_comp_all, p_comp_all = spearmanr(valid_all["learned_communication_score"], valid_all["compositional_z_score"])
    print(f"\n[All {len(valid_all)} rows] Spearman(learned, edge_z_score)        = {rho_edge_all:.4f} (p={p_edge_all:.4g})")
    print(f"[All {len(valid_all)} rows] Spearman(learned, compositional_z_score) = {rho_comp_all:.4f} (p={p_comp_all:.4g})")

    # --- Spearman correlations: 16 CXCL12 rows only ---
    cxcl12 = merged[merged["ligand"] == "CXCL12"].dropna(
        subset=["learned_communication_score", "edge_z_score", "compositional_z_score"])
    rho_edge_cxcl12, p_edge_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["edge_z_score"])
    rho_comp_cxcl12, p_comp_cxcl12 = spearmanr(cxcl12["learned_communication_score"], cxcl12["compositional_z_score"])
    print(f"\n[16 CXCL12 rows] Spearman(learned, edge_z_score)        = {rho_edge_cxcl12:.4f} (p={p_edge_cxcl12:.4g})")
    print(f"[16 CXCL12 rows] Spearman(learned, compositional_z_score) = {rho_comp_cxcl12:.4f} (p={p_comp_cxcl12:.4g})")

    # --- CD274->PDCD1 rows: values + rank within full 22-row ranking (descending by learned score) ---
    merged_ranked = merged.dropna(subset=["learned_communication_score"]).copy()
    merged_ranked["learned_rank"] = merged_ranked["learned_communication_score"].rank(
        ascending=False, method="min").astype(int)
    n_ranked = len(merged_ranked)

    cd274 = merged_ranked[(merged_ranked["ligand"] == "CD274") & (merged_ranked["receptor"] == "PDCD1")]
    print(f"\nCD274->PDCD1 rows (edge_observed_score=0.000 in ground truth), ranked among {n_ranked} rows "
          "(rank 1 = highest learned score):")
    print(cd274[["sender", "receiver", "learned_communication_score", "learned_rank",
                 "edge_z_score", "compositional_z_score"]].to_string(index=False))

    bottom_half_threshold = n_ranked / 2.0  # ranks > this are "bottom half"
    cd274_bottom_half = (cd274["learned_rank"] > bottom_half_threshold).all()
    print(f"\nAll 6 CD274->PDCD1 rows in bottom half (rank > {bottom_half_threshold}): {cd274_bottom_half}")
    print("CD274->PDCD1 ranks:", sorted(cd274["learned_rank"].tolist()))

    # --- Full ranking printout for transparency ---
    print("\nFull 22-row ranking by learned_communication_score (descending):")
    print(merged_ranked.sort_values("learned_rank")[
        ["sender", "receiver", "ligand", "receptor", "learned_communication_score",
         "learned_rank", "edge_z_score", "compositional_z_score"]
    ].to_string(index=False))

    t3 = time.time()
    step2_runtime = t3 - t2
    print(f"\nStep 2 (scoring/aggregation/comparison) runtime: {step2_runtime:.2f}s")

    # ---------------- STEP 3: anti-confound check ----------------
    def meaningfully_weaker(rho_edge, rho_comp, tol=0.05):
        # "meaningfully weaker" = edge-rho is more than `tol` below comp-rho.
        return rho_edge < (rho_comp - tol)

    confound_all_weaker = meaningfully_weaker(rho_edge_all, rho_comp_all)
    confound_cxcl12_weaker = meaningfully_weaker(rho_edge_cxcl12, rho_comp_cxcl12)
    edge_at_least_as_strong_cxcl12 = rho_edge_cxcl12 >= (rho_comp_cxcl12 - 0.05)

    print("\n--- STEP 3: anti-confound check ---")
    print(f"All rows:    edge_z_score rho={rho_edge_all:.4f}  vs  compositional_z_score rho={rho_comp_all:.4f}  "
          f"-> edge meaningfully weaker than compositional: {confound_all_weaker}")
    print(f"16 CXCL12:   edge_z_score rho={rho_edge_cxcl12:.4f}  vs  compositional_z_score rho={rho_comp_cxcl12:.4f}  "
          f"-> edge meaningfully weaker than compositional: {confound_cxcl12_weaker}")
    if confound_all_weaker or confound_cxcl12_weaker:
        print("CONFOUND WARNING: the learned score correlates with compositional_z_score at least as "
              "strongly as with edge_z_score in at least one grouping -- this suggests the GNN may have "
              "re-derived plain expression/cell-type composition from the PCA features rather than "
              "genuine spatially-restricted contact structure.")
    else:
        print("No confound detected by this check: the learned score is not meaningfully more correlated "
              "with compositional_z_score than with edge_z_score in either grouping.")

    # ---------------- STEP 4: pre-registered success criterion ----------------
    criterion_a = bool(cd274_bottom_half)
    criterion_b = bool(rho_edge_cxcl12 > 0.6)
    criterion_c = bool(edge_at_least_as_strong_cxcl12)  # "not meaningfully weaker" at the 16-CXCL12 level

    n_met = sum([criterion_a, criterion_b, criterion_c])
    if n_met == 3:
        verdict = "MET"
    elif n_met == 0:
        verdict = "NOT MET"
    else:
        verdict = "PARTIAL"

    print("\n--- STEP 4: pre-registered success criterion ---")
    print(f"(a) CD274->PDCD1 all in bottom half of 22-row ranking: {criterion_a}")
    print(f"(b) Spearman rho(learned, edge_z_score) > 0.6 across 16 CXCL12 rows: {criterion_b} "
          f"(actual rho={rho_edge_cxcl12:.4f})")
    print(f"(c) edge_z_score correlation not meaningfully weaker than compositional_z_score "
          f"(16-CXCL12 level, tol=0.05): {criterion_c}")
    print(f"VERDICT: {verdict}")

    # ---------------- save summary artifacts ----------------
    summary = {
        "step1_reproduction": {
            "prior_run_best_val_auc": prior_auc,
            "this_run_best_val_auc": float(best_val_auc),
            "abs_diff": float(auc_diff),
            "reproduction_close": bool(reproduction_close),
            "best_epoch": int(best_epoch),
            "n_epochs_run": len(history["epoch"]),
        },
        "spearman_all_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_all), "p": float(p_edge_all), "n": int(len(valid_all))},
            "vs_compositional_z_score": {"rho": float(rho_comp_all), "p": float(p_comp_all), "n": int(len(valid_all))},
        },
        "spearman_16_cxcl12_rows": {
            "vs_edge_z_score": {"rho": float(rho_edge_cxcl12), "p": float(p_edge_cxcl12), "n": int(len(cxcl12))},
            "vs_compositional_z_score": {"rho": float(rho_comp_cxcl12), "p": float(p_comp_cxcl12), "n": int(len(cxcl12))},
        },
        "cd274_pdcd1": {
            "ranks": [int(r) for r in sorted(cd274["learned_rank"].tolist())],
            "n_total_ranked": int(n_ranked),
            "bottom_half_threshold_rank": float(bottom_half_threshold),
            "all_bottom_half": bool(cd274_bottom_half),
            "scores": {
                f"{row.sender}->{row.receiver}": row.learned_communication_score
                for row in cd274.itertuples()
            },
        },
        "step3_confound_check": {
            "all_rows_edge_meaningfully_weaker_than_compositional": bool(confound_all_weaker),
            "cxcl12_rows_edge_meaningfully_weaker_than_compositional": bool(confound_cxcl12_weaker),
        },
        "step4_criteria": {
            "a_cd274_bottom_half": criterion_a,
            "b_rho_gt_0.6": criterion_b,
            "c_edge_not_meaningfully_weaker": criterion_c,
            "verdict": verdict,
        },
        "runtime_sec_step1": step1_runtime,
        "runtime_sec_step2": step2_runtime,
    }
    with open(os.path.join(OUT_DIR, "gnn_pca_step2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved:", os.path.join(OUT_DIR, "gnn_pca_step2_summary.json"))

    paragraph = (
        f"The PCA-features GraphSAGE link-prediction model was reproduced with best validation AUC "
        f"{best_val_auc:.4f} (prior run: {prior_auc:.4f}, abs diff {auc_diff:.4f}), confirming the earlier "
        f"gate-passing result. Its trained decoder's learned edge score correlates with edge_z_score at "
        f"rho={rho_edge_cxcl12:.3f} and with compositional_z_score at rho={rho_comp_cxcl12:.3f} across the 16 "
        f"CXCL12-CXCR4 pairs (all-22-row level: edge rho={rho_edge_all:.3f}, compositional rho={rho_comp_all:.3f}). "
        f"{'The score is NOT meaningfully more correlated with compositional structure than with the spatial-edge ground truth.' if not (confound_all_weaker or confound_cxcl12_weaker) else 'The score correlates with compositional_z_score at least as strongly as with edge_z_score, indicating the model likely re-derived expression/cell-type composition rather than genuine spatial-contact structure.'} "
        f"All 6 CD274->PDCD1 rows rank {sorted(cd274['learned_rank'].tolist())} of {n_ranked} by learned score "
        f"({'all in the bottom half' if cd274_bottom_half else 'NOT all in the bottom half'}). "
        f"Pre-registered success criterion: (a) CD274->PDCD1 bottom-half = {criterion_a}; "
        f"(b) rho>0.6 on 16 CXCL12 rows = {criterion_b} (actual {rho_edge_cxcl12:.3f}); "
        f"(c) edge-vs-compositional not meaningfully weaker = {criterion_c}. "
        f"Overall verdict: {verdict}."
    )
    with open(os.path.join(OUT_DIR, "gnn_pca_comparison_summary.txt"), "w") as f:
        f.write(paragraph + "\n")
    print("\nSaved comparison paragraph:\n", paragraph)

    return summary


if __name__ == "__main__":
    main()

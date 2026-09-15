"""Step 2 (only run if Step 1 gate passed): embedding-based communication score.

Reuses get_spatial_graph_edges from src/spatial_edge_score.py (imported, not modified).
Reads data/processed/xenium_breast_baseline.h5ad (read-only) and the scGPT embeddings
saved by scgpt_step1_embedding_gate.py. Writes only to results/embedding_baseline/.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr

sys.path.insert(0, r"F:\Thesis\src")
from spatial_edge_score import get_spatial_graph_edges  # noqa: E402

H5AD_PATH = r"F:\Thesis\data\processed\xenium_breast_baseline.h5ad"
OUT_DIR = r"F:\Thesis\results\embedding_baseline"
NULL_MODEL_CSV = r"F:\Thesis\results\null_model_comparison.csv"


def main():
    t0 = time.time()

    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    print("Loaded adata:", adata.shape)

    scgpt_emb = np.load(os.path.join(OUT_DIR, "scgpt_cell_embeddings.npy"))
    assert scgpt_emb.shape[0] == n_cells, "embedding cell count mismatch with adata"
    print("Loaded scGPT embeddings:", scgpt_emb.shape)

    row_idx, col_idx = get_spatial_graph_edges(adata)
    print("n_directed_edges:", len(row_idx))

    # Cosine similarity per directed edge i -> j
    norms = np.linalg.norm(scgpt_emb, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    emb_unit = scgpt_emb / norms
    edge_cos_sim = np.sum(emb_unit[row_idx] * emb_unit[col_idx], axis=1)

    labels = adata.obs["cell_type"].astype(str).values
    sender_types = labels[row_idx]
    receiver_types = labels[col_idx]

    null_df = pd.read_csv(NULL_MODEL_CSV)
    unique_pairs = null_df[["sender", "receiver"]].drop_duplicates().values.tolist()

    records = []
    for s, r in unique_pairs:
        mask = (sender_types == s) & (receiver_types == r)
        n_edges = int(np.sum(mask))
        emb_score = float(np.mean(edge_cos_sim[mask])) if n_edges > 0 else np.nan
        records.append({
            "sender": s,
            "receiver": r,
            "embedding_similarity_score": emb_score,
            "n_qualifying_edges_embedding": n_edges,
        })

    emb_df = pd.DataFrame(records)

    merged = null_df.merge(emb_df, on=["sender", "receiver"], how="left")
    out_csv = os.path.join(OUT_DIR, "embedding_communication_scores.csv")
    merged.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

    # --- Spearman correlation: embedding similarity ranking vs edge_z_score ranking ---
    # Use the per-(sender,receiver) unique-pair level (embedding score is identical
    # across ligand-receptor rows sharing the same sender-receiver pair; for ranking
    # against edge_z_score we correlate row-by-row across all 22 interactions, since
    # edge_z_score varies per LR pair even when sender-receiver repeats).
    valid = merged.dropna(subset=["embedding_similarity_score", "edge_z_score"])
    rho, pval = spearmanr(valid["embedding_similarity_score"], valid["edge_z_score"])
    print(f"Spearman rho (embedding_similarity_score vs edge_z_score), n={len(valid)}: rho={rho:.4f}, p={pval:.4g}")

    # Sender-receiver-pair-level ranking only (unique pairs, one row each)
    rho_pairs, pval_pairs = spearmanr(
        emb_df.merge(null_df.groupby(["sender", "receiver"])["edge_z_score"].mean().reset_index(),
                      on=["sender", "receiver"])["embedding_similarity_score"],
        emb_df.merge(null_df.groupby(["sender", "receiver"])["edge_z_score"].mean().reset_index(),
                      on=["sender", "receiver"])["edge_z_score"],
    )
    print(f"Spearman rho (unique sender-receiver pairs, mean edge_z_score): rho={rho_pairs:.4f}, p={pval_pairs:.4g}")

    # --- CD274->PDCD1 check ---
    cd274 = merged[(merged["ligand"] == "CD274") & (merged["receptor"] == "PDCD1")]
    print("\nCD274->PDCD1 rows (edge_observed_score is 0.000 for all in ground truth):")
    print(cd274[["sender", "receiver", "edge_observed_score", "edge_z_score", "embedding_similarity_score"]].to_string(index=False))

    all_emb_scores = emb_df["embedding_similarity_score"].dropna()
    cd274_emb_scores = cd274["embedding_similarity_score"].dropna()
    print(f"\nAll sender-receiver pairs embedding score: mean={all_emb_scores.mean():.4f}, "
          f"min={all_emb_scores.min():.4f}, max={all_emb_scores.max():.4f}")
    print(f"CD274->PDCD1 pairs embedding score: mean={cd274_emb_scores.mean():.4f}, "
          f"min={cd274_emb_scores.min():.4f}, max={cd274_emb_scores.max():.4f}")
    cd274_rank_pct = [
        (emb_df["embedding_similarity_score"] <= v).mean() for v in cd274_emb_scores
    ]
    print("CD274->PDCD1 embedding-score percentile rank among all sender-receiver pairs:",
          [f"{p:.2f}" for p in cd274_rank_pct])

    # --- fibroblast->epithelial / fibroblast->immune (CXCL12-CXCR4) top-rank check ---
    cxcl12 = merged[(merged["ligand"] == "CXCL12") & (merged["receptor"] == "CXCR4")].copy()
    cxcl12_sorted_by_emb = cxcl12.sort_values("embedding_similarity_score", ascending=False).reset_index(drop=True)
    print("\nCXCL12->CXCR4 rows ranked by embedding_similarity_score (descending):")
    print(cxcl12_sorted_by_emb[["sender", "receiver", "embedding_similarity_score", "edge_z_score"]].to_string(index=False))

    fib_immune_rank = cxcl12_sorted_by_emb.index[
        (cxcl12_sorted_by_emb["sender"] == "fibroblast") & (cxcl12_sorted_by_emb["receiver"] == "immune")
    ].tolist()
    fib_epi_rank = cxcl12_sorted_by_emb.index[
        (cxcl12_sorted_by_emb["sender"] == "fibroblast") & (cxcl12_sorted_by_emb["receiver"] == "epithelial")
    ].tolist()
    print(f"fibroblast->immune rank by embedding similarity (0=top): {fib_immune_rank}")
    print(f"fibroblast->epithelial rank by embedding similarity (0=top): {fib_epi_rank}")

    t1 = time.time()
    step2_runtime = t1 - t0
    print(f"\nStep 2 runtime: {step2_runtime:.2f} s")

    summary = {
        "n_directed_edges": int(len(row_idx)),
        "n_unique_sender_receiver_pairs": int(len(unique_pairs)),
        "spearman_rho_all_rows": float(rho),
        "spearman_pvalue_all_rows": float(pval),
        "spearman_rho_unique_pairs": float(rho_pairs),
        "spearman_pvalue_unique_pairs": float(pval_pairs),
        "cd274_pdcd1_embedding_scores": {
            f"{row.sender}->{row.receiver}": row.embedding_similarity_score
            for row in cd274.itertuples()
        },
        "all_pairs_embedding_score_mean": float(all_emb_scores.mean()),
        "all_pairs_embedding_score_min": float(all_emb_scores.min()),
        "all_pairs_embedding_score_max": float(all_emb_scores.max()),
        "fibroblast_immune_rank_by_embedding": fib_immune_rank,
        "fibroblast_epithelial_rank_by_embedding": fib_epi_rank,
        "runtime_sec": step2_runtime,
    }
    with open(os.path.join(OUT_DIR, "step2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    return summary


if __name__ == "__main__":
    main()

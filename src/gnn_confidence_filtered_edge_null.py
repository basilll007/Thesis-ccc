"""GNN-confidence-filtered edge null model: making the GNN communication score
ligand-receptor-aware.

Prior round's `learned_communication_score` (results/gnn_edge_score/gnn_pca_communication_scores.csv)
was identical across LR pairs whenever sender/receiver cell types matched (e.g.
fibroblast->fibroblast scored 0.936226 for both CXCL12-CXCR4 and CD274-PDCD1),
because it only used the trained GNN's per-edge link-prediction confidence
(a function of node identity/structure only) with no reference to which genes
were being tested.

This module fixes that by using the GNN's per-edge confidence purely as a
FILTER on which physical edges count toward the existing exact edge-product
score (expr_L[i] * expr_R[j], from spatial_edge_score.py), rather than as the
score itself. Restricting the edge-product average to only the top-confidence
fraction of each (sender,receiver) pair's qualifying edges makes the result
ligand-receptor-aware again (since edge_products differs per LR pair), while
weighting toward edges the trained GNN considers more "real" spatial contacts.

Reuses, unmodified:
  - get_spatial_graph_edges, precompute_edge_products (spatial_edge_score.py)
  - generate_label_permutations (spatial_null.py)
  - the exact z-score / empirical-p-value convention from run_spatial_edge_null.py

Does NOT retrain the GNN or recompute PCA features -- loads the already-saved
results/gnn_edge_score/gnn_pca_node_embeddings.npy only.
"""
from __future__ import annotations

import json
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import PROCESSED_DIR, RESULTS_DIR, SEED  # noqa: E402
from spatial_edge_score import get_spatial_graph_edges, precompute_edge_products  # noqa: E402
from spatial_null import generate_label_permutations  # noqa: E402

H5AD_PATH = PROCESSED_DIR / "xenium_breast_baseline.h5ad"
NULL_MODEL_CSV = RESULTS_DIR / "null_model_comparison.csv"
EMB_PATH = RESULTS_DIR / "gnn_edge_score" / "gnn_pca_node_embeddings.npy"
OUT_DIR = RESULTS_DIR / "gnn_edge_score"

N_PERMS = 500
TOP_FRAC = 0.5
CORRECTNESS_TOL = 1e-9


def filtered_score_for_labels(
    labels: np.ndarray,
    row_indices: np.ndarray,
    col_indices: np.ndarray,
    edge_products: dict[tuple[str, str], np.ndarray],
    interactions: list[tuple[str, str, str, str]],
    p_edge: np.ndarray,
    top_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Mirrors compute_edge_scores_for_labels's pair-masking logic exactly,
    with one addition: within each (sender,receiver) qualifying-edge mask,
    keep only the top `top_frac` fraction of edges by p_edge (computed
    within that qualifying subset only) before averaging expr_L*expr_R.

    At top_frac=1.0 this reproduces compute_edge_scores_for_labels exactly.
    """
    sender_types = labels[row_indices]
    receiver_types = labels[col_indices]

    n_interactions = len(interactions)
    scores = np.full(n_interactions, np.nan, dtype=np.float64)
    edge_counts = np.zeros(n_interactions, dtype=np.int64)

    unique_type_pairs = set([(s, r) for s, r, _, _ in interactions])
    pair_masks: dict[tuple[str, str], np.ndarray] = {}
    for s, r in unique_type_pairs:
        pair_masks[(s, r)] = (sender_types == s) & (receiver_types == r)

    for i, (s, r, lig, rec) in enumerate(interactions):
        mask = pair_masks[(s, r)]
        qualifying_idx = np.where(mask)[0]
        n_qual = len(qualifying_idx)
        edge_counts[i] = n_qual
        if n_qual == 0:
            continue

        p_qual = p_edge[qualifying_idx]
        k = max(1, int(round(top_frac * n_qual)))
        if k >= n_qual:
            selected_idx = qualifying_idx
        else:
            top_local = np.argpartition(p_qual, -k)[-k:]
            selected_idx = qualifying_idx[top_local]

        prod_arr = edge_products[(lig, rec)]
        scores[i] = float(np.mean(prod_arr[selected_idx]))

    return scores, edge_counts


def edge_null_stats(observed_score: float, null_dist: np.ndarray) -> tuple[float, float, float]:
    """Exact convention from run_spatial_edge_null.py."""
    null_dist_clean = np.where(np.isnan(null_dist), 0.0, null_dist)
    null_mean = float(np.mean(null_dist_clean))
    null_std = float(np.std(null_dist_clean, ddof=1)) if len(null_dist_clean) > 1 else 0.0
    empirical_p = float(np.mean(null_dist_clean >= observed_score))
    z_score = float((observed_score - null_mean) / null_std) if null_std > 1e-12 else 0.0
    return null_mean, null_std, empirical_p, z_score


def main() -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    print("Loaded adata:", adata.shape)

    embeddings = np.load(EMB_PATH)
    assert embeddings.shape[0] == n_cells, "embedding cell count mismatch with adata"
    print("Loaded GNN PCA-features node embeddings:", embeddings.shape)

    row_idx, col_idx = get_spatial_graph_edges(adata)
    print("n_directed_edges:", len(row_idx))

    # p_edge: sigmoid(dot product) between endpoint embeddings, vectorized, no loop.
    dot = np.sum(embeddings[row_idx] * embeddings[col_idx], axis=1)
    p_edge = 1.0 / (1.0 + np.exp(-dot))
    print(f"p_edge: min={p_edge.min():.4f} max={p_edge.max():.4f} mean={p_edge.mean():.4f}")

    null_df = pd.read_csv(NULL_MODEL_CSV)
    interactions: list[tuple[str, str, str, str]] = [
        (row["sender"], row["receiver"], row["ligand"], row["receptor"])
        for _, row in null_df.iterrows()
    ]
    lr_pairs = list(set([(lig, rec) for _, _, lig, rec in interactions]))
    edge_products = precompute_edge_products(adata, lr_pairs, row_idx, col_idx)

    orig_labels = np.asarray(adata.obs["cell_type"].astype(str).to_numpy())

    t_step01 = time.time()
    step0_runtime = t_step01 - t0

    # ---------------- STEP 1: CORRECTNESS GATE (top_frac=1.0) ----------------
    reproduced_scores, _ = filtered_score_for_labels(
        labels=orig_labels,
        row_indices=row_idx,
        col_indices=col_idx,
        edge_products=edge_products,
        interactions=interactions,
        p_edge=p_edge,
        top_frac=1.0,
    )
    existing_scores = null_df["edge_observed_score"].to_numpy(dtype=np.float64)
    diffs = np.abs(reproduced_scores - existing_scores)

    print("\n=== STEP 1 CORRECTNESS GATE (top_frac=1.0 vs existing edge_observed_score) ===")
    gate_table = pd.DataFrame({
        "sender": null_df["sender"], "receiver": null_df["receiver"],
        "ligand": null_df["ligand"], "receptor": null_df["receptor"],
        "existing_edge_observed_score": existing_scores,
        "reproduced_score_top_frac_1.0": reproduced_scores,
        "abs_diff": diffs,
    })
    print(gate_table.to_string(index=False))

    max_diff = float(np.max(diffs))
    gate_pass = bool(max_diff < CORRECTNESS_TOL)
    print(f"\nMax abs diff across all 22 rows: {max_diff:.3e}")
    print(f"CORRECTNESS GATE (all rows < {CORRECTNESS_TOL:.0e}): {'PASS' if gate_pass else 'FAIL'}")

    t_step1 = time.time()
    step1_runtime = t_step1 - t0

    result: dict = {
        "step0_runtime_sec": step0_runtime,
        "step1_runtime_sec_total": step1_runtime,
        "correctness_gate": {
            "max_abs_diff": max_diff,
            "tolerance": CORRECTNESS_TOL,
            "gate_pass": gate_pass,
            "per_row_diffs": gate_table.to_dict(orient="records"),
        },
    }

    if not gate_pass:
        mismatched = gate_table[gate_table["abs_diff"] >= CORRECTNESS_TOL]
        print("\nMISMATCHED ROWS:")
        print(mismatched.to_string(index=False))
        result["gate_pass"] = False
        with open(OUT_DIR / "gnn_confidence_filtered_CORRECTNESS_GATE_FAILED.json", "w") as f:
            json.dump(result, f, indent=2, default=str)
        print("\nSTOPPED: correctness gate failed. No Step 2-4 artifacts produced.")
        return result

    print("\nCORRECTNESS GATE PASSED. Proceeding to Step 2.")

    # ---------------- STEP 2: filtered observed scores (top_frac=0.5) + null ----------------
    t2 = time.time()

    filtered_observed, filtered_edge_counts = filtered_score_for_labels(
        labels=orig_labels,
        row_indices=row_idx,
        col_indices=col_idx,
        edge_products=edge_products,
        interactions=interactions,
        p_edge=p_edge,
        top_frac=TOP_FRAC,
    )
    print(f"\nComputed {len(interactions)} filtered observed scores at top_frac={TOP_FRAC}.")

    perms = generate_label_permutations(orig_labels, n_perms=N_PERMS, seed=SEED)
    null_matrix = np.zeros((N_PERMS, len(interactions)), dtype=np.float64)
    for p, p_labels in enumerate(perms):
        p_scores, _ = filtered_score_for_labels(
            labels=p_labels,
            row_indices=row_idx,
            col_indices=col_idx,
            edge_products=edge_products,
            interactions=interactions,
            p_edge=p_edge,
            top_frac=TOP_FRAC,
        )
        null_matrix[p, :] = p_scores

    filtered_null_mean = np.zeros(len(interactions))
    filtered_null_std = np.zeros(len(interactions))
    filtered_p = np.zeros(len(interactions))
    filtered_z = np.zeros(len(interactions))
    for i in range(len(interactions)):
        nm, ns, ep, z = edge_null_stats(filtered_observed[i], null_matrix[:, i])
        filtered_null_mean[i] = nm
        filtered_null_std[i] = ns
        filtered_p[i] = ep
        filtered_z[i] = z

    t3 = time.time()
    step2_runtime = t3 - t2
    print(f"Step 2 (500-permutation filtered null) runtime: {step2_runtime:.2f}s")

    # ---------------- STEP 3: comparison ----------------
    merged = null_df.copy()
    merged["filtered_observed_score"] = filtered_observed
    merged["filtered_null_mean"] = filtered_null_mean
    merged["filtered_null_std"] = filtered_null_std
    merged["filtered_empirical_pvalue"] = filtered_p
    merged["filtered_z_score"] = filtered_z
    merged["n_qualifying_edges_filtered"] = filtered_edge_counts

    out_csv = OUT_DIR / "gnn_confidence_filtered_null.csv"
    merged.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

    cxcl12 = merged[merged["ligand"] == "CXCL12"].copy()
    rho, pval = spearmanr(cxcl12["filtered_z_score"], cxcl12["edge_z_score"])
    print(f"\n[16 CXCL12 rows] Spearman(filtered_z_score, edge_z_score) = {rho:.4f} (p={pval:.4g})")

    top5 = merged.reindex(merged["edge_z_score"].abs().sort_values(ascending=False).index).head(5)
    # Task says "5 rows with the highest original edge_z_score" -> use signed edge_z_score descending
    top5 = merged.sort_values("edge_z_score", ascending=False).head(5).copy()
    med_filtered_abs = float(top5["filtered_z_score"].abs().median())
    med_edge_abs = float(top5["edge_z_score"].abs().median())
    print("\nTop-5 rows by highest edge_z_score:")
    print(top5[["sender", "receiver", "ligand", "receptor", "edge_z_score", "filtered_z_score"]].to_string(index=False))
    print(f"\nMedian |edge_z_score| (top-5) = {med_edge_abs:.4f}")
    print(f"Median |filtered_z_score| (top-5) = {med_filtered_abs:.4f}")
    sharpened = med_filtered_abs > med_edge_abs
    print(f"Filtering {'SHARPENED' if sharpened else 'DID NOT sharpen (diluted or unchanged)'} the top signal.")

    cd274 = merged[(merged["ligand"] == "CD274") & (merged["receptor"] == "PDCD1")].copy()
    print("\n=== CD274->PDCD1 rows -- SANITY CHECK ONLY, NOT SCORED ===")
    print(cd274[["sender", "receiver", "filtered_observed_score", "filtered_z_score",
                 "filtered_empirical_pvalue", "edge_observed_score", "edge_z_score"]].to_string(index=False))
    cd274_stayed_near_zero = bool((cd274["filtered_observed_score"].abs() < 1e-6).all())
    cd274_nonsignificant = bool((cd274["filtered_empirical_pvalue"] >= 0.05).all() or
                                 (cd274["filtered_z_score"].abs() < 1.96).all())
    print(f"CD274->PDCD1 filtered_observed_score all ~0 (< 1e-6): {cd274_stayed_near_zero}")
    print(f"CD274->PDCD1 stayed non-significant (p>=0.05 or |z|<1.96): {cd274_nonsignificant}")

    # ---------------- STEP 4: pre-registered success criterion ----------------
    criterion_a = bool(rho > 0.7)
    criterion_b = bool(sharpened)
    n_met = sum([criterion_a, criterion_b])
    verdict = "MET" if n_met == 2 else ("NOT MET" if n_met == 0 else "PARTIAL")

    print("\n=== STEP 4: pre-registered success criterion (CXCL12-CXCR4 only) ===")
    print(f"(a) Spearman rho > 0.7 (16 CXCL12 rows): {criterion_a} (actual rho={rho:.4f})")
    print(f"(b) median|filtered_z| > median|edge_z| (top-5 rows): {criterion_b} "
          f"({med_filtered_abs:.4f} vs {med_edge_abs:.4f})")
    print(f"VERDICT: {verdict}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    ax0 = axes[0]
    ax0.scatter(cxcl12["edge_z_score"], cxcl12["filtered_z_score"], s=70, color="#3182bd",
                edgecolors="black", linewidths=0.6, zorder=3)
    for _, r in cxcl12.iterrows():
        ax0.annotate(f"{r['sender'][:3]}→{r['receiver'][:3]}",
                     (r["edge_z_score"], r["filtered_z_score"]), fontsize=7,
                     xytext=(3, 3), textcoords="offset points", alpha=0.8)
    lims = [min(cxcl12["edge_z_score"].min(), cxcl12["filtered_z_score"].min()) - 2,
            max(cxcl12["edge_z_score"].max(), cxcl12["filtered_z_score"].max()) + 2]
    ax0.plot(lims, lims, linestyle="--", color="#999999", linewidth=1.0, label="y = x", zorder=1)
    ax0.set_xlim(lims); ax0.set_ylim(lims)
    ax0.set_xlabel("edge_z_score (original, all qualifying edges)")
    ax0.set_ylabel("filtered_z_score (top-50% GNN-confidence edges)")
    ax0.set_title(f"CXCL12→CXCR4 (n=16): Spearman ρ={rho:.3f}, p={pval:.3g}")
    ax0.legend(loc="upper left", fontsize=9)
    ax0.grid(alpha=0.3)

    ax1 = axes[1]
    labels_bar = [f"{r['sender']}→{r['receiver']}\n({r['ligand']}-{r['receptor']})" for _, r in top5.iterrows()]
    x = np.arange(len(top5))
    width = 0.35
    ax1.bar(x - width / 2, top5["edge_z_score"].abs(), width, label="|edge_z_score| (original)",
            color="#999999", edgecolor="black")
    ax1.bar(x + width / 2, top5["filtered_z_score"].abs(), width, label="|filtered_z_score| (GNN-filtered)",
            color="#3182bd", edgecolor="black")
    ax1.axhline(med_edge_abs, color="#636363", linestyle=":", linewidth=1.0)
    ax1.axhline(med_filtered_abs, color="#08519c", linestyle=":", linewidth=1.0)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels_bar, fontsize=8)
    ax1.set_ylabel("|z-score|")
    ax1.set_title(f"Top-5 by edge_z_score: median |edge_z|={med_edge_abs:.2f} vs median |filtered_z|={med_filtered_abs:.2f}")
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3, axis="y")

    fig.suptitle("GNN-Confidence-Filtered Edge Null vs. Original Edge Null (top_frac=0.5)", fontsize=13)
    fig.tight_layout()
    fig_path_png = OUT_DIR / "gnn_confidence_filtered_comparison.png"
    fig_path_pdf = OUT_DIR / "gnn_confidence_filtered_comparison.pdf"
    fig.savefig(fig_path_png, dpi=150, bbox_inches="tight")
    fig.savefig(fig_path_pdf, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved figure:", fig_path_png, fig_path_pdf)

    # ---------------- save summary ----------------
    summary = {
        "step0_runtime_sec": step0_runtime,
        "step1_runtime_sec_total": step1_runtime,
        "step2_runtime_sec": step2_runtime,
        "correctness_gate": {"max_abs_diff": max_diff, "gate_pass": gate_pass},
        "top_frac": TOP_FRAC,
        "n_perms": N_PERMS,
        "spearman_16_cxcl12": {"rho": float(rho), "p": float(pval)},
        "top5_by_edge_z_score": top5[["sender", "receiver", "ligand", "receptor",
                                       "edge_z_score", "filtered_z_score"]].to_dict(orient="records"),
        "median_abs_edge_z_top5": med_edge_abs,
        "median_abs_filtered_z_top5": med_filtered_abs,
        "filtering_sharpened_top_signal": sharpened,
        "cd274_pdcd1_sanity_check_not_scored": cd274[["sender", "receiver", "filtered_observed_score",
                                                        "filtered_z_score", "filtered_empirical_pvalue"]].to_dict(orient="records"),
        "cd274_stayed_near_zero": cd274_stayed_near_zero,
        "cd274_stayed_nonsignificant": cd274_nonsignificant,
        "step4_criteria": {
            "a_rho_gt_0.7": criterion_a,
            "b_filtering_sharpens_top5": criterion_b,
            "verdict": verdict,
        },
    }
    with open(OUT_DIR / "gnn_confidence_filtered_step_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print("\nSaved:", OUT_DIR / "gnn_confidence_filtered_step_summary.json")

    paragraph = (
        f"Step 1 correctness gate: at top_frac=1.0, the new filtered_score_for_labels function "
        f"reproduces the existing edge_observed_score for all 22 interactions with a maximum absolute "
        f"difference of {max_diff:.2e} (tolerance 1e-9) -- PASSED. At top_frac=0.5 (fixed, GNN-confidence-filtered "
        f"edges only), filtered_z_score correlates with the original edge_z_score at Spearman rho={rho:.3f} "
        f"(p={pval:.3g}) across the 16 CXCL12-CXCR4 rows. Among the 5 rows with the highest original edge_z_score, "
        f"median |edge_z_score|={med_edge_abs:.3f} vs median |filtered_z_score|={med_filtered_abs:.3f} -- filtering "
        f"{'sharpened' if sharpened else 'did not sharpen'} the already-established top signal. "
        f"CD274->PDCD1 (SANITY CHECK -- NOT SCORED): filtered_observed_score stayed at essentially 0 for all 6 rows "
        f"({'confirmed' if cd274_stayed_near_zero else 'NOT confirmed'}, all < 1e-6) and remained non-significant "
        f"({'confirmed' if cd274_nonsignificant else 'NOT confirmed'}), as expected since this is a panel-detection "
        f"limit (near-zero PDCD1 expression) rather than a spatial-contact limit that GNN-confidence filtering could "
        f"ever change. Pre-registered success criterion (CXCL12-CXCR4 only): (a) rho>0.7 = {criterion_a} "
        f"(actual {rho:.3f}); (b) filtering sharpens the top-5 signal = {criterion_b}. Overall verdict: {verdict}."
    )
    with open(OUT_DIR / "gnn_confidence_filtered_comparison_summary.txt", "w") as f:
        f.write(paragraph + "\n")
    print("\nSaved comparison paragraph:\n", paragraph)

    result.update(summary)
    return result


if __name__ == "__main__":
    main()

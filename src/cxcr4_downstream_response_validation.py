"""Enhancement 3: Downstream Receiver Response Gene Activation Validation.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
Writes to results/cxcr4_response_validation/ and figures/.

Key design:
  - Provides mechanistic biological validation for the primary spatial communication axis:
    fibroblast -> immune (CXCL12 -> CXCR4).
  - Evaluates whether receiver cells in direct physical contact with CXCL12+ sender cells
    exhibit upregulation of downstream signaling response genes (MAPK feedback, motility,
    cytoskeleton, and cell cycle: DUSP2, DUSP5, MAP3K8, PIM1, CDC42EP1, CTTN, S100A4,
    MMP2, SNAI1, CCND1, MKI67).
  - Partitions receiver cells into Contact-Engaged (incoming CXCL12 exposure > 0) vs
    Unengaged (incoming CXCL12 exposure == 0).
  - Computes log2 fold changes, Welch's t-test, Mann-Whitney U test, Benjamini-Hochberg
    FDR q-values, Cohen's d effect sizes, and composite pathway activation scores.
  - Runs a 500-permutation spatial null test to confirm spatial specificity.
  - Generates publication-grade figures (spatial activation map and volcano/effect size plot).
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
from scipy.stats import mannwhitneyu, ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
from spatial_edge_score import get_spatial_graph_edges

SEED = 42
N_PERMS = 500

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
OUT_DIR = ROOT / "results" / "cxcr4_response_validation"
FIGURES_DIR = ROOT / "figures"

# Candidate downstream response genes present in Xenium 280-gene panel
RESPONSE_GENES = [
    "DUSP2",     # Dual-specificity phosphatase 2 (MAPK/ERK negative feedback)
    "DUSP5",     # Dual-specificity phosphatase 5 (nuclear MAPK phosphatase)
    "MAP3K8",    # MAP kinase kinase kinase 8 / Tpl2 (direct MAPK activator)
    "PIM1",      # Serine/threonine kinase (chemokine survival signaling)
    "CDC42EP1",  # CDC42 effector protein 1 (actin cytoskeletal remodeling)
    "CTTN",      # Cortactin (chemokine-induced migration/motility)
    "S100A4",    # Metastasis/motility protein downstream of chemokine motility
    "MMP2",      # Matrix metalloproteinase 2 (ECM invasion)
    "SNAI1",     # Snail family transcriptional repressor (motility/EMT)
    "CCND1",     # Cyclin D1 (cell cycle entry downstream of CXCR4-ERK)
    "MKI67",     # Ki-67 (proliferation marker)
]


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Compute Benjamini-Hochberg False Discovery Rate q-values."""
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    q_values = np.zeros(n, dtype=np.float64)

    current_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        q = (sorted_p[i] * n) / rank
        current_min = min(current_min, q)
        q_values[i] = current_min

    out_q = np.zeros(n, dtype=np.float64)
    out_q[sorted_indices] = q_values
    return out_q


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Compute Cohen's d effect size between two groups."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / max(n1 + n2 - 2, 1))
    if pooled_sd <= 1e-12:
        return 0.0
    return float((np.mean(group1) - np.mean(group2)) / pooled_sd)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    row_idx, col_idx = get_spatial_graph_edges(adata)
    print(f"Loaded spatial graph: {n_cells} cells, {len(row_idx)} directed edges.")

    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)

    # Extract CXCL12 and CXCR4 expression
    cxcl12_col = raw_source.X[:, var_names.index("CXCL12")]
    if hasattr(cxcl12_col, "toarray"):
        cxcl12_col = cxcl12_col.toarray()
    cxcl12_expr = np.asarray(cxcl12_col).ravel().astype(np.float64)

    cxcr4_col = raw_source.X[:, var_names.index("CXCR4")]
    if hasattr(cxcr4_col, "toarray"):
        cxcr4_col = cxcr4_col.toarray()
    cxcr4_expr = np.asarray(cxcr4_col).ravel().astype(np.float64)

    # Filter available response genes in panel
    valid_response_genes = [g for g in RESPONSE_GENES if g in var_names]
    print(f"Validated response genes in panel ({len(valid_response_genes)}/{len(RESPONSE_GENES)}): {valid_response_genes}")

    # Extract response gene expression matrix
    resp_matrix = np.zeros((n_cells, len(valid_response_genes)), dtype=np.float64)
    for g_idx, g_name in enumerate(valid_response_genes):
        c = raw_source.X[:, var_names.index(g_name)]
        if hasattr(c, "toarray"):
            c = c.toarray()
        resp_matrix[:, g_idx] = np.asarray(c).ravel().astype(np.float64)

    # Calculate incoming CXCL12 exposure for each receiver cell j
    # row_idx = sender (i), col_idx = receiver (j)
    incoming_ligand_exposure = np.zeros(n_cells, dtype=np.float64)
    np.add.at(incoming_ligand_exposure, col_idx, cxcl12_expr[row_idx])

    # Receiver cells of interest: immune cells or CXCR4-expressing cells
    cell_types = adata.obs["cell_type"].astype(str).values
    is_immune = cell_types == "immune"
    cxcr4_pos = cxcr4_expr > 0

    # Primary analysis subset: CXCR4-expressing immune and stromal receiver cells
    # We analyze receiver cells with CXCR4 expression to test receptor-dependent downstream response
    receiver_mask = is_immune | cxcr4_pos
    receiver_indices = np.where(receiver_mask)[0]
    n_receivers = len(receiver_indices)
    print(f"Total receiver cells evaluated (Immune or CXCR4+): {n_receivers} / {n_cells}")

    rec_exposure = incoming_ligand_exposure[receiver_indices]
    engaged_mask = rec_exposure > 0
    n_engaged = int(np.sum(engaged_mask))
    n_unengaged = n_receivers - n_engaged
    print(f"Contact-Engaged receivers (incoming CXCL12 > 0): {n_engaged} ({n_engaged/n_receivers*100:.1f}%)")
    print(f"Unengaged receivers (incoming CXCL12 == 0):     {n_unengaged} ({n_unengaged/n_receivers*100:.1f}%)")

    engaged_indices = receiver_indices[engaged_mask]
    unengaged_indices = receiver_indices[~engaged_mask]

    # Compute composite z-scored activation score
    z_resp_matrix = np.zeros_like(resp_matrix)
    for g_idx in range(len(valid_response_genes)):
        col_vals = resp_matrix[:, g_idx]
        sd = np.std(col_vals)
        mn = np.mean(col_vals)
        z_resp_matrix[:, g_idx] = (col_vals - mn) / max(sd, 1e-12)
    activation_signature = np.mean(z_resp_matrix, axis=1)

    # ---------------- DIFFERENTIAL EXPRESSION TESTING ----------------
    records = []
    p_values_ttest = []
    p_values_mwu = []

    for g_idx, g_name in enumerate(valid_response_genes):
        vals_eng = resp_matrix[engaged_indices, g_idx]
        vals_uneng = resp_matrix[unengaged_indices, g_idx]

        mean_eng = float(np.mean(vals_eng))
        mean_uneng = float(np.mean(vals_uneng))
        pct_eng = float(100.0 * np.mean(vals_eng > 0))
        pct_uneng = float(100.0 * np.mean(vals_uneng > 0))

        eps = 1e-3
        log2_fc = float(np.log2((mean_eng + eps) / (mean_uneng + eps)))
        d_val = cohens_d(vals_eng, vals_uneng)

        # Welch's t-test
        t_stat, t_pval = ttest_ind(vals_eng, vals_uneng, equal_var=False)
        # Mann-Whitney U test
        u_stat, u_pval = mannwhitneyu(vals_eng, vals_uneng, alternative="two-sided")

        p_values_ttest.append(float(t_pval))
        p_values_mwu.append(float(u_pval))

        records.append({
            "gene": g_name,
            "mean_engaged": mean_eng,
            "mean_unengaged": mean_uneng,
            "pct_detected_engaged": pct_eng,
            "pct_detected_unengaged": pct_uneng,
            "log2_fold_change": log2_fc,
            "cohens_d": d_val,
            "welch_t_stat": float(t_stat),
            "welch_pvalue": float(t_pval),
            "mwu_stat": float(u_stat),
            "mwu_pvalue": float(u_pval),
        })

    diff_df = pd.DataFrame(records)
    diff_df["fdr_qvalue_welch"] = benjamini_hochberg(np.asarray(p_values_ttest))
    diff_df["fdr_qvalue_mwu"] = benjamini_hochberg(np.asarray(p_values_mwu))
    diff_df["significant_upregulated"] = (diff_df["log2_fold_change"] > 0) & (diff_df["fdr_qvalue_welch"] < 0.05)

    diff_csv_path = OUT_DIR / "response_genes_differential.csv"
    diff_df.to_csv(diff_csv_path, index=False)
    print(f"\nSaved differential expression table: {diff_csv_path}")
    print(diff_df[["gene", "mean_engaged", "mean_unengaged", "log2_fold_change", "cohens_d", "fdr_qvalue_welch", "significant_upregulated"]].to_string(index=False))

    # Composite activation score comparison
    eng_act = activation_signature[engaged_indices]
    uneng_act = activation_signature[unengaged_indices]
    act_diff = float(np.mean(eng_act) - np.mean(uneng_act))
    act_d = cohens_d(eng_act, uneng_act)
    act_t, act_p = ttest_ind(eng_act, uneng_act, equal_var=False)
    print(f"\nComposite Pathway Activation: Engaged={np.mean(eng_act):.3f} vs Unengaged={np.mean(uneng_act):.3f} | Diff={act_diff:+.3f}, Cohen's d={act_d:.3f}, p={act_p:.4e}")

    # ---------------- SPATIAL PERMUTATION NULL TEST ----------------
    print(f"\nRunning {N_PERMS} spatial coordinate shuffles to test spatial specificity of activation...")
    rng = np.random.default_rng(SEED)
    null_act_diffs = np.zeros(N_PERMS, dtype=np.float64)

    for p_idx in range(N_PERMS):
        shuffled_cxcl12 = rng.permutation(cxcl12_expr)
        perm_incoming = np.zeros(n_cells, dtype=np.float64)
        np.add.at(perm_incoming, col_idx, shuffled_cxcl12[row_idx])

        perm_rec_exp = perm_incoming[receiver_indices]
        perm_eng = perm_rec_exp > 0
        if np.sum(perm_eng) > 0 and np.sum(~perm_eng) > 0:
            null_act_diffs[p_idx] = float(np.mean(activation_signature[receiver_indices[perm_eng]]) - np.mean(activation_signature[receiver_indices[~perm_eng]]))
        else:
            null_act_diffs[p_idx] = 0.0

    spatial_p = float(np.mean(null_act_diffs >= act_diff))
    null_mean_diff = float(np.mean(null_act_diffs))
    null_std_diff = float(np.std(null_act_diffs))
    spatial_z = float((act_diff - null_mean_diff) / max(null_std_diff, 1e-12))
    print(f"Spatial Permutation Test: Null Mean={null_mean_diff:.4f}, Null Std={null_std_diff:.4f}, Spatial z={spatial_z:.2f}, p_spatial={spatial_p:.4f}")

    # Save summary JSON
    n_sig_up = int(diff_df["significant_upregulated"].sum())
    summary = {
        "n_receiver_cells": n_receivers,
        "n_contact_engaged": n_engaged,
        "n_unengaged": n_unengaged,
        "composite_activation": {
            "mean_engaged": float(np.mean(eng_act)),
            "mean_unengaged": float(np.mean(uneng_act)),
            "difference": act_diff,
            "cohens_d": act_d,
            "welch_pvalue": float(act_p),
        },
        "spatial_permutation_null": {
            "n_perms": N_PERMS,
            "spatial_z_score": spatial_z,
            "spatial_pvalue": spatial_p,
            "spatially_significant": bool(spatial_p < 0.05),
        },
        "response_genes_tested": len(valid_response_genes),
        "n_significantly_upregulated": n_sig_up,
        "top_upregulated_genes": diff_df.sort_values("cohens_d", ascending=False).head(5)["gene"].tolist(),
        "runtime_sec": time.time() - t0,
    }
    with open(OUT_DIR / "response_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved: {OUT_DIR / 'response_summary.json'}")

    # ---------------- GENERATE PUBLICATION FIGURES ----------------
    print("\n--- Generating Publication Figures ---")

    # Figure 1: Spatial Map of Contact Engagement & Pathway Activation
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), sharex=True, sharey=True)

    # Panel A: Spatial Engagement
    ax0 = axes[0]
    ax0.scatter(coords[:, 0], coords[:, 1], s=1.5, c="#e0e0e0", alpha=0.3, rasterized=True)
    ax0.scatter(coords[unengaged_indices, 0], coords[unengaged_indices, 1], s=8, c="#337ab7", alpha=0.75, label=f"Unengaged (E=0, n={n_unengaged})", rasterized=True)
    ax0.scatter(coords[engaged_indices, 0], coords[engaged_indices, 1], s=14, c="#d9534f", alpha=0.9, label=f"Contact-Engaged (E>0, n={n_engaged})", rasterized=True)
    ax0.set_aspect("equal")
    ax0.invert_yaxis()
    ax0.set_title(f"Spatial Receptor Engagement (CXCL12 → CXCR4)\n{n_engaged} Engaged vs. {n_unengaged} Unengaged Receivers", fontsize=11, fontweight="bold")
    ax0.set_xlabel("Spatial X (µm)", fontsize=10)
    ax0.set_ylabel("Spatial Y (µm)", fontsize=10)
    ax0.legend(frameon=True, fontsize=9, loc="upper right")

    # Panel B: Continuous Pathway Activation Score
    ax1 = axes[1]
    ax1.scatter(coords[:, 0], coords[:, 1], s=1.5, c="#e0e0e0", alpha=0.3, rasterized=True)
    sc = ax1.scatter(
        coords[receiver_indices, 0],
        coords[receiver_indices, 1],
        s=10,
        c=activation_signature[receiver_indices],
        cmap="plasma",
        alpha=0.85,
        rasterized=True,
    )
    cbar = fig.colorbar(sc, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("CXCR4 Pathway Activation Score (z)", fontsize=10)
    ax1.set_aspect("equal")
    ax1.invert_yaxis()
    ax1.set_title(f"Receiver Cell Pathway Activation Signature\nCohen's d = {act_d:+.3f} (p = {act_p:.2e}, Spatial z = {spatial_z:.1f})", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Spatial X (µm)", fontsize=10)

    fig.suptitle("Enhancement 3: Downstream CXCR4 Signaling Activation Map", fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "cxcr4_response_activation_spatial.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "cxcr4_response_activation_spatial.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "cxcr4_response_activation_spatial.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "cxcr4_response_activation_spatial.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved cxcr4_response_activation_spatial.png/.pdf")

    # Figure 2: Differential Expression & Effect Sizes
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))

    # Volcano Plot
    ax0 = axes[0]
    p_log = -np.log10(np.clip(diff_df["welch_pvalue"].values, 1e-30, 1.0))
    fc = diff_df["log2_fold_change"].values
    colors = ["#d9534f" if sig else "#7f7f7f" for sig in diff_df["significant_upregulated"]]

    ax0.scatter(fc, p_log, c=colors, s=50, alpha=0.85, edgecolor="k", lw=0.5)
    ax0.axhline(-np.log10(0.05), color="blue", linestyle=":", lw=1.2, label="p = 0.05")
    ax0.axvline(0, color="k", linestyle="--", lw=0.8, alpha=0.6)

    for _, row in diff_df.iterrows():
        if row["significant_upregulated"] or abs(row["log2_fold_change"]) > 0.2:
            ax0.annotate(
                row["gene"],
                (row["log2_fold_change"], -np.log10(max(row["welch_pvalue"], 1e-30))),
                xytext=(5, 3),
                textcoords="offset points",
                fontsize=9,
                fontweight="bold" if row["significant_upregulated"] else "normal",
            )
    ax0.set_xlabel("Log2 Fold Change (Engaged / Unengaged)", fontsize=10)
    ax0.set_ylabel("-Log10 (Welch p-value)", fontsize=10)
    ax0.set_title("Downstream Response Gene Volcano Plot", fontsize=11, fontweight="bold")
    ax0.grid(True, linestyle=":", alpha=0.6)
    ax0.legend(frameon=True, fontsize=9)

    # Bar Plot of Effect Sizes (Cohen's d)
    ax1 = axes[1]
    sorted_diff = diff_df.sort_values("cohens_d", ascending=True)
    bar_colors = ["#2ca02c" if d > 0 and q < 0.05 else ("#1f77b4" if d > 0 else "#7f7f7f") for d, q in zip(sorted_diff["cohens_d"], sorted_diff["fdr_qvalue_welch"])]
    bars = ax1.barh(sorted_diff["gene"], sorted_diff["cohens_d"], color=bar_colors, alpha=0.85, height=0.65)
    ax1.axvline(0, color="k", linestyle="-", lw=0.8)
    ax1.set_xlabel("Cohen's d Effect Size", fontsize=10)
    ax1.set_title("Effect Size by Gene (Engaged vs. Unengaged Receivers)", fontsize=11, fontweight="bold")
    ax1.grid(axis="x", linestyle=":", alpha=0.6)

    for bar, (_, row) in zip(bars, sorted_diff.iterrows()):
        q = row["fdr_qvalue_welch"]
        star = "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else ""))
        d = row["cohens_d"]
        ax1.text(d + (0.01 if d >= 0 else -0.01), bar.get_y() + bar.get_height() / 2.0, f"{d:+.2f} {star}", va="center", ha="left" if d >= 0 else "right", fontsize=8)

    fig.suptitle("Enhancement 3: Mechanistic Verification of Downstream CXCR4 Signaling", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "cxcr4_response_differential_expression.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "cxcr4_response_differential_expression.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "cxcr4_response_differential_expression.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "cxcr4_response_differential_expression.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved cxcr4_response_differential_expression.png/.pdf")

    total_time = time.time() - t0
    print(f"\nEnhancement 3 completed successfully in {total_time:.2f}s")
    return summary


if __name__ == "__main__":
    main()

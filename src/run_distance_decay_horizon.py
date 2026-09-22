"""Option 1: Multi-Scale Spatial Distance-Decay Horizon & Receiver Dose-Response Kinetics.

Analyzes the effective paracrine diffusion range of CXCL12-CXCR4 in 10x Xenium breast carcinoma:
1. Distance-Decay Horizon Sweep: Radii r in [15, 30, 50, 75, 100, 150, 200, 300] um.
   - Evaluates observed communication score, qualifying edge count, and 500-permutation spatial null z-score.
   - Evaluates continuous exponential diffusion kernel w(d) = exp(-d / lambda).
2. Single-Cell Receiver Dose-Response Kinetics:
   - Measures downstream response gene expression (S100A4, MMP2, MAP3K8, PIM1, CCND1, MKI67)
     in single receiver cells as a continuous function of Euclidean distance to the nearest CXCL12-secreting fibroblast.
   - Fits exponential decay kinetics to quantify empirical biological action half-distance (d_1/2).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

SEED = 42
N_PERMS = 500
RADII = [15.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0]  # in microns (um)
DIFFUSION_LAMBDAS = [15.0, 30.0, 50.0, 100.0, 200.0]

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
OUT_DIR = ROOT / "results" / "distance_decay_horizon"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def exp_decay_model(d, y_inf, a, lam):
    """Exponential decay: y(d) = y_inf + a * exp(-d / lam)."""
    return y_inf + a * np.exp(-d / np.maximum(lam, 1e-3))


def run_distance_decay_horizon():
    t0 = time.time()
    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values
    n_cells = adata.n_obs

    # Gene expression extraction
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)

    def get_gene_expr(gene_name: str) -> np.ndarray:
        if gene_name not in var_names:
            return np.zeros(n_cells, dtype=np.float32)
        idx = var_names.index(gene_name)
        col = raw_source.X[:, idx]
        if hasattr(col, "toarray"):
            col = col.toarray()
        return np.asarray(col, dtype=np.float32).flatten()

    cxcl12_expr = get_gene_expr("CXCL12")
    cxcr4_expr = get_gene_expr("CXCR4")
    cd274_expr = get_gene_expr("CD274")
    pdcd1_expr = get_gene_expr("PDCD1")

    # Downstream response genes
    response_gene_names = ["S100A4", "MMP2", "MAP3K8", "PIM1", "CCND1", "MKI67"]
    response_expr = {g: get_gene_expr(g) for g in response_gene_names}

    print(f"Dataset loaded: {n_cells} cells, coords range X: [{coords[:,0].min():.1f}, {coords[:,0].max():.1f}], Y: [{coords[:,1].min():.1f}, {coords[:,1].max():.1f}] um.")

    # -------------------------------------------------------------
    # 1. Distance-Decay Horizon Sweep across Radii r
    # -------------------------------------------------------------
    print("\n--- Part 1: Sweeping Distance Horizon Radii (r in 15 to 300 um) ---")
    horizon_records = []

    # Pre-generate label permutations
    rng = np.random.RandomState(SEED)
    perm_labels_list = [rng.permutation(cell_types) for _ in range(N_PERMS)]

    nn_engine = NearestNeighbors(metric="euclidean", n_jobs=-1)
    nn_engine.fit(coords)

    for r in RADII:
        t_r0 = time.time()
        # Find neighbors within radius r
        adj_radius = nn_engine.radius_neighbors(coords, radius=r, return_distance=True)
        
        row_list = []
        col_list = []
        dist_list = []
        for i, (dists, neighbors) in enumerate(zip(adj_radius[0], adj_radius[1])):
            mask = neighbors != i  # exclude self
            if np.any(mask):
                row_list.extend([i] * np.sum(mask))
                col_list.extend(neighbors[mask])
                dist_list.extend(dists[mask])

        rows = np.asarray(row_list, dtype=np.int64)
        cols = np.asarray(col_list, dtype=np.int64)
        dists = np.asarray(dist_list, dtype=np.float32)
        n_edges = len(rows)

        # Precompute edge products
        cxcl12_cxcr4_prod = cxcl12_expr[rows] * cxcr4_expr[cols]
        cd274_pdcd1_prod = cd274_expr[rows] * pdcd1_expr[cols]

        # Evaluate target axis: fibroblast -> immune (CXCL12 -> CXCR4)
        fib_imm_mask_obs = (cell_types[rows] == "fibroblast") & (cell_types[cols] == "immune")
        n_fib_imm = int(np.sum(fib_imm_mask_obs))
        obs_score_cxcl12 = float(np.mean(cxcl12_cxcr4_prod[fib_imm_mask_obs])) if n_fib_imm > 0 else 0.0

        # Evaluate tumor -> immune (CD274 -> PDCD1)
        epi_imm_mask_obs = (cell_types[rows] == "epithelial") & (cell_types[cols] == "immune")
        n_epi_imm = int(np.sum(epi_imm_mask_obs))
        obs_score_cd274 = float(np.mean(cd274_pdcd1_prod[epi_imm_mask_obs])) if n_epi_imm > 0 else 0.0

        # 500 permutations for CXCL12-CXCR4
        null_cxcl12 = np.zeros(N_PERMS, dtype=np.float64)
        null_cd274 = np.zeros(N_PERMS, dtype=np.float64)

        for p_idx, perm_labels in enumerate(perm_labels_list):
            p_fib_imm = (perm_labels[rows] == "fibroblast") & (perm_labels[cols] == "immune")
            if np.any(p_fib_imm):
                null_cxcl12[p_idx] = np.mean(cxcl12_cxcr4_prod[p_fib_imm])

            p_epi_imm = (perm_labels[rows] == "epithelial") & (perm_labels[cols] == "immune")
            if np.any(p_epi_imm):
                null_cd274[p_idx] = np.mean(cd274_pdcd1_prod[p_epi_imm])

        # CXCL12 statistics
        null_mean_12 = float(np.mean(null_cxcl12))
        null_std_12 = float(np.std(null_cxcl12, ddof=1))
        z_score_12 = float((obs_score_cxcl12 - null_mean_12) / null_std_12) if null_std_12 > 1e-9 else 0.0
        p_val_12 = float(np.mean(null_cxcl12 >= obs_score_cxcl12))

        # CD274 statistics
        null_mean_274 = float(np.mean(null_cd274))
        null_std_274 = float(np.std(null_cd274, ddof=1))
        z_score_274 = float((obs_score_cd274 - null_mean_274) / null_std_274) if null_std_274 > 1e-9 else 0.0
        p_val_274 = float(np.mean(null_cd274 >= obs_score_cd274))

        # Continuous Exponential Diffusion Kernel over this graph: lambda = r / 2
        lam = r / 2.0
        w_d = np.exp(-dists / lam)
        w_fib_imm = w_d[fib_imm_mask_obs]
        diff_score_cxcl12 = float(np.sum(w_fib_imm * cxcl12_cxcr4_prod[fib_imm_mask_obs]) / np.maximum(np.sum(w_fib_imm), 1e-9)) if n_fib_imm > 0 else 0.0

        rec = {
            "radius_um": r,
            "total_edges": n_edges,
            "mean_edges_per_cell": float(n_edges / n_cells),
            "fib_imm_edges": n_fib_imm,
            "epi_imm_edges": n_epi_imm,
            "cxcl12_obs_score": obs_score_cxcl12,
            "cxcl12_null_mean": null_mean_12,
            "cxcl12_null_std": null_std_12,
            "cxcl12_z_score": z_score_12,
            "cxcl12_pvalue": p_val_12,
            "cxcl12_diff_kernel_score": diff_score_cxcl12,
            "cd274_obs_score": obs_score_cd274,
            "cd274_null_mean": null_mean_274,
            "cd274_null_std": null_std_274,
            "cd274_z_score": z_score_274,
            "cd274_pvalue": p_val_274,
            "runtime_sec": round(time.time() - t_r0, 2),
        }
        horizon_records.append(rec)
        print(f"Radius {r:5.1f} um: edges={n_edges:6d} | fib->imm edges={n_fib_imm:5d} | CXCL12 z={z_score_12:6.2f} (p={p_val_12:.4f}) | CD274 z={z_score_274:5.2f}")

    df_horizon = pd.DataFrame(horizon_records)
    df_horizon.to_csv(OUT_DIR / "distance_decay_summary.csv", index=False)
    print(f"Saved horizon summary: {OUT_DIR / 'distance_decay_summary.csv'}")

    # -------------------------------------------------------------
    # 2. Single-Cell Receiver Dose-Response Kinetics
    # -------------------------------------------------------------
    print("\n--- Part 2: Single-Cell Receiver Distance Dose-Response Kinetics ---")
    
    # Identify CXCL12-secreting fibroblasts
    fib_mask = (cell_types == "fibroblast") & (cxcl12_expr > 0)
    if not np.any(fib_mask):
        fib_mask = cell_types == "fibroblast"
    fib_coords = coords[fib_mask]
    print(f"Number of active CXCL12-producing fibroblasts: {len(fib_coords)}")

    # Candidate receiver cells (immune cells or CXCR4-expressing cells)
    receiver_mask = (cell_types == "immune") | (cxcr4_expr > 0)
    receiver_indices = np.where(receiver_mask)[0]
    receiver_coords = coords[receiver_indices]
    n_receivers = len(receiver_indices)
    print(f"Number of receiver cells evaluated: {n_receivers}")

    # Compute Euclidean distance from each receiver to the nearest CXCL12+ fibroblast
    nn_fib = NearestNeighbors(n_neighbors=1, metric="euclidean")
    nn_fib.fit(fib_coords)
    dists_to_nearest_fib, _ = nn_fib.kneighbors(receiver_coords)
    dists_to_nearest_fib = dists_to_nearest_fib.flatten()  # in microns

    # Distance bins
    bin_edges = [0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 1000.0]
    bin_labels = [
        "0-25 um", "25-50 um", "50-75 um", "75-100 um",
        "100-150 um", "150-200 um", "200-300 um", ">300 um"
    ]
    bin_idx = np.digitize(dists_to_nearest_fib, bin_edges) - 1
    bin_idx = np.clip(bin_idx, 0, len(bin_labels) - 1)

    # Build per-cell dataframe
    df_cell_dose = pd.DataFrame({
        "cell_id": adata.obs_names[receiver_indices],
        "cell_type": cell_types[receiver_indices],
        "cxcr4_expr": cxcr4_expr[receiver_indices],
        "dist_to_cxcl12_fib_um": dists_to_nearest_fib,
        "dist_bin": [bin_labels[b] for b in bin_idx],
        "dist_bin_idx": bin_idx,
    })
    for g in response_gene_names:
        df_cell_dose[f"{g}_expr"] = response_expr[g][receiver_indices]

    # Aggregate by distance bin
    bin_records = []
    for b_i, b_lbl in enumerate(bin_labels):
        m_b = bin_idx == b_i
        n_b = int(np.sum(m_b))
        if n_b == 0:
            continue
        d_mean = float(np.mean(dists_to_nearest_fib[m_b]))
        d_med = float(np.median(dists_to_nearest_fib[m_b]))
        rec = {
            "bin_index": b_i,
            "bin_label": b_lbl,
            "n_cells": n_b,
            "mean_dist_um": round(d_mean, 2),
            "median_dist_um": round(d_med, 2),
        }
        for g in response_gene_names:
            vals = response_expr[g][receiver_indices][m_b]
            rec[f"{g}_mean"] = float(np.mean(vals))
            rec[f"{g}_sem"] = float(np.std(vals, ddof=1) / np.sqrt(n_b)) if n_b > 1 else 0.0
            rec[f"{g}_pct_detected"] = float(np.mean(vals > 0) * 100.0)
        bin_records.append(rec)

    df_bin_dose = pd.DataFrame(bin_records)
    df_bin_dose.to_csv(OUT_DIR / "receiver_dose_response.csv", index=False)
    print(f"Saved binned dose-response: {OUT_DIR / 'receiver_dose_response.csv'}")

    # Fit exponential decay for S100A4 and MMP2
    kinetics_results = {}
    for g in ["S100A4", "MMP2", "MAP3K8"]:
        y_vals = df_bin_dose[f"{g}_mean"].values
        d_vals = df_bin_dose["mean_dist_um"].values
        # Spearman correlation with continuous distance
        rho, p_rho = spearmanr(dists_to_nearest_fib, response_expr[g][receiver_indices])

        try:
            # Fit: y(d) = y_inf + a * exp(-d / lam)
            p0 = [float(y_vals[-1]), float(y_vals[0] - y_vals[-1]), 50.0]
            bounds = ([0.0, 0.0, 1.0], [np.inf, np.inf, 500.0])
            popt, _ = curve_fit(exp_decay_model, d_vals, y_vals, p0=p0, bounds=bounds, maxfev=5000)
            y_inf, a, lam = popt
            half_dist = float(lam * np.log(2.0))
            kinetics_results[g] = {
                "baseline_y_inf": round(float(y_inf), 4),
                "amplitude_a": round(float(a), 4),
                "decay_length_lambda_um": round(float(lam), 2),
                "half_distance_d12_um": round(half_dist, 2),
                "spearman_rho_vs_dist": round(float(rho), 4),
                "spearman_pvalue": float(p_rho),
                "fit_success": True,
            }
            print(f"Kinetic Fit for {g:6s}: Half-Distance d_1/2 = {half_dist:5.1f} um | Decay lambda = {lam:5.1f} um | Spearman rho = {rho:+.4f} (p={p_rho:.2e})")
        except Exception as e:
            kinetics_results[g] = {
                "fit_success": False,
                "error": str(e),
                "spearman_rho_vs_dist": round(float(rho), 4),
                "spearman_pvalue": float(p_rho),
            }

    # Find Critical Signaling Radius r* where CXCL12 z-score peaks
    best_row = df_horizon.loc[df_horizon["cxcl12_z_score"].idxmax()]
    r_star = float(best_row["radius_um"])
    max_z = float(best_row["cxcl12_z_score"])

    summary_out = {
        "critical_signaling_radius_r_star_um": r_star,
        "max_spatial_z_score": max_z,
        "radii_evaluated_um": RADII,
        "kinetics": kinetics_results,
        "cd274_isolated_radii_um": [float(r) for r in df_horizon.loc[df_horizon["cd274_obs_score"] == 0.0, "radius_um"].values],
        "n_receivers_evaluated": n_receivers,
        "n_producers_evaluated": len(fib_coords),
        "runtime_total_sec": round(time.time() - t0, 2),
    }

    with open(OUT_DIR / "horizon_summary.json", "w") as f:
        json.dump(summary_out, f, indent=2)
    print(f"Saved complete horizon summary: {OUT_DIR / 'horizon_summary.json'}")
    print(f"Critical Signaling Radius r* = {r_star} um with peak z-score = {max_z:.2f}")


if __name__ == "__main__":
    run_distance_decay_horizon()

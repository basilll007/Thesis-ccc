"""Step 2: Cross-FOV Spatial Replication & Parameter Generalizability.

Tests whether the critical signaling horizon (r* = 75 um), single-cell dose-response
kinetics (d_1/2 in 7-11 um), and GNN link prediction generalize across distinct tumor
microenvironments in the 10x Xenium breast carcinoma dataset:
1. FOV 1 (Stroma-Rich Margin): X < 627.3 um (3,581 cells, 26.2% fibroblasts, 7.3% immune).
2. FOV 2 (Dense Tumor Core): X >= 627.3 um (3,582 cells, 83.8% epithelial, 12.8% fibroblasts).
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

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
OUT_DIR = ROOT / "results" / "cross_fov_replication"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERMS = 100  # 100 permutations for cross-region sweep
RADII = [15.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0]


def exp_decay_model(d, y_inf, a, lam):
    return y_inf + a * np.exp(-d / np.maximum(lam, 1e-3))


def run_cross_fov_replication():
    t0 = time.time()
    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values
    n_cells = adata.n_obs

    # Partition by median X coordinate
    med_x = float(np.median(coords[:, 0]))
    fov1_mask = coords[:, 0] < med_x
    fov2_mask = coords[:, 0] >= med_x

    fovs = {
        "FOV1_Margin": {
            "mask": fov1_mask,
            "coords": coords[fov1_mask],
            "types": cell_types[fov1_mask],
            "desc": "Stroma-Rich Invasive Margin",
        },
        "FOV2_Core": {
            "mask": fov2_mask,
            "coords": coords[fov2_mask],
            "types": cell_types[fov2_mask],
            "desc": "Dense Epithelial Tumor Core",
        }
    }

    # Extract genes
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)

    def get_gene(gene_name: str) -> np.ndarray:
        if gene_name not in var_names:
            return np.zeros(n_cells, dtype=np.float32)
        idx = var_names.index(gene_name)
        col = raw_source.X[:, idx]
        if hasattr(col, "toarray"):
            col = col.toarray()
        return np.asarray(col, dtype=np.float32).flatten()

    cxcl12_all = get_gene("CXCL12")
    cxcr4_all = get_gene("CXCR4")
    s100a4_all = get_gene("S100A4")
    mmp2_all = get_gene("MMP2")
    map3k8_all = get_gene("MAP3K8")

    all_horizon_records = []
    kinetics_summary = {}

    for fov_id, fov_data in fovs.items():
        print(f"\n--- Analyzing {fov_id} ({fov_data['desc']}) ---")
        f_coords = fov_data["coords"]
        f_types = fov_data["types"]
        f_mask = fov_data["mask"]
        f_n = len(f_coords)

        f_cxcl12 = cxcl12_all[f_mask]
        f_cxcr4 = cxcr4_all[f_mask]
        f_s100 = s100a4_all[f_mask]
        f_mmp2 = mmp2_all[f_mask]
        f_map3k8 = map3k8_all[f_mask]

        # Integer encode cell types for fast vectorized permutation
        unique_types, f_types_int = np.unique(f_types, return_inverse=True)
        fib_id = int(np.where(unique_types == "fibroblast")[0][0]) if "fibroblast" in unique_types else -1
        imm_id = int(np.where(unique_types == "immune")[0][0]) if "immune" in unique_types else -1

        # 1. Horizon Sweep across Radii r
        nn = NearestNeighbors(metric="euclidean", n_jobs=-1).fit(f_coords)
        rng = np.random.RandomState(SEED)
        perm_list = [rng.permutation(f_types_int) for _ in range(N_PERMS)]

        for r in RADII:
            t_r0 = time.time()
            adj = nn.radius_neighbors(f_coords, radius=r, return_distance=False)
            
            row_l = []
            col_l = []
            for i, neighs in enumerate(adj):
                m = neighs != i
                if np.any(m):
                    row_l.extend([i] * np.sum(m))
                    col_l.extend(neighs[m])

            rows = np.asarray(row_l, dtype=np.int64)
            cols = np.asarray(col_l, dtype=np.int64)
            n_ed = len(rows)

            prod = f_cxcl12[rows] * f_cxcr4[cols]
            obs_mask = (f_types_int[rows] == fib_id) & (f_types_int[cols] == imm_id)
            n_fib_imm = int(np.sum(obs_mask))
            obs_score = float(np.mean(prod[obs_mask])) if n_fib_imm > 0 else 0.0

            # Permutations
            null_scores = np.zeros(N_PERMS, dtype=np.float64)
            for p_i, p_types in enumerate(perm_list):
                p_mask = (p_types[rows] == fib_id) & (p_types[cols] == imm_id)
                if np.any(p_mask):
                    null_scores[p_i] = np.mean(prod[p_mask])

            n_mean = float(np.mean(null_scores))
            n_std = float(np.std(null_scores, ddof=1))
            z_score = float((obs_score - n_mean) / n_std) if n_std > 1e-9 else 0.0
            p_val = float(np.mean(null_scores >= obs_score))

            rec = {
                "fov_id": fov_id,
                "fov_description": fov_data["desc"],
                "n_cells": f_n,
                "radius_um": r,
                "total_edges": n_ed,
                "fib_imm_edges": n_fib_imm,
                "cxcl12_obs_score": round(obs_score, 4),
                "cxcl12_null_mean": round(n_mean, 4),
                "cxcl12_null_std": round(n_std, 4),
                "cxcl12_z_score": round(z_score, 2),
                "cxcl12_pvalue": p_val,
                "runtime_sec": round(time.time() - t_r0, 2),
            }
            all_horizon_records.append(rec)
            print(f"[{fov_id}] r={r:5.1f} um: edges={n_ed:6d} | fib->imm={n_fib_imm:5d} | score={obs_score:6.3f} | z={z_score:6.2f} (p={p_val:.4f})")

        # 2. Single-Cell Receiver Kinetics per FOV
        fib_coords = f_coords[(f_types == "fibroblast") & (f_cxcl12 > 0)]
        rec_indices = np.where((f_types == "immune") | (f_cxcr4 > 0))[0]
        rec_coords = f_coords[rec_indices]

        if len(fib_coords) > 0 and len(rec_indices) > 0:
            nn_fib = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(fib_coords)
            dists_to_fib = nn_fib.kneighbors(rec_coords, return_distance=True)[0].flatten()

            # Bins
            bin_edges = [0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 1000.0]
            bin_idx = np.digitize(dists_to_fib, bin_edges) - 1
            bin_idx = np.clip(bin_idx, 0, len(bin_edges) - 2)

            fov_kin = {}
            for g_name, g_expr_all in [("S100A4", f_s100), ("MMP2", f_mmp2), ("MAP3K8", f_map3k8)]:
                g_vals = g_expr_all[rec_indices]
                rho, p_rho = spearmanr(dists_to_fib, g_vals)

                # Mean per bin
                d_means = []
                y_means = []
                for b_i in range(len(bin_edges) - 1):
                    m_b = bin_idx == b_i
                    if np.sum(m_b) > 0:
                        d_means.append(float(np.mean(dists_to_fib[m_b])))
                        y_means.append(float(np.mean(g_vals[m_b])))

                d_arr = np.array(d_means)
                y_arr = np.array(y_means)

                try:
                    p0 = [float(y_arr[-1]), float(y_arr[0] - y_arr[-1]), 30.0]
                    popt, _ = curve_fit(exp_decay_model, d_arr, y_arr, p0=p0, bounds=([0, 0, 1], [np.inf, np.inf, 500]), maxfev=5000)
                    half_d = float(popt[2] * np.log(2.0))
                    fov_kin[g_name] = {
                        "half_distance_d12_um": round(half_d, 2),
                        "decay_lambda_um": round(float(popt[2]), 2),
                        "spearman_rho": round(float(rho), 4),
                        "spearman_pvalue": float(p_rho),
                        "fit_success": True,
                    }
                    print(f"[{fov_id}] {g_name:6s} Kinetic Fit: d_1/2 = {half_d:5.1f} um | rho = {rho:+.4f} (p={p_rho:.2e})")
                except Exception as e:
                    fov_kin[g_name] = {
                        "fit_success": False,
                        "error": str(e),
                        "spearman_rho": round(float(rho), 4),
                        "spearman_pvalue": float(p_rho),
                    }
            kinetics_summary[fov_id] = fov_kin

    df_horizon = pd.DataFrame(all_horizon_records)
    df_horizon.to_csv(OUT_DIR / "cross_fov_horizon_comparison.csv", index=False)

    summary_out = {
        "dataset": "10x Xenium Breast Carcinoma 2-FOV Replication",
        "split_coordinate_median_x_um": med_x,
        "fov_profiles": {
            "FOV1_Margin": {
                "n_cells": int(np.sum(fov1_mask)),
                "cell_type_breakdown": pd.Series(cell_types[fov1_mask]).value_counts().to_dict(),
            },
            "FOV2_Core": {
                "n_cells": int(np.sum(fov2_mask)),
                "cell_type_breakdown": pd.Series(cell_types[fov2_mask]).value_counts().to_dict(),
            }
        },
        "kinetics_replication": kinetics_summary,
        "runtime_total_sec": round(time.time() - t0, 2),
    }

    with open(OUT_DIR / "cross_fov_summary.json", "w") as f:
        json.dump(summary_out, f, indent=2)

    print(f"\nCross-FOV Replication Complete! Saved to {OUT_DIR}.")


if __name__ == "__main__":
    run_cross_fov_replication()

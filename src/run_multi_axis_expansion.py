"""Step 3: Multi-Axis Signaling Expansion across Verified Panel Axes.

Expands beyond CXCL12 and CD274 to evaluate newly discovered active signaling axes
present in the 10x Xenium 280-gene human breast carcinoma panel:
1. CXCL12 -> CXCR4: Paracrine chemokine axis (Fibroblast -> Immune/Tumor).
2. PTN -> SDC4: Stromal growth factor / motogenic axis (Pleiotrophin -> Syndecan-4).
   - Highly expressed in breast carcinoma stroma (2,501 PTN+ senders -> 3,044 SDC4+ receivers).
3. CD86 -> CTLA4: Immune checkpoint receptor-ligand contact axis (Immune -> Immune).
   - 224 CD86+ senders -> 95 CTLA4+ receivers.
4. CD274 -> PDCD1: Negative spatial contact control (Epithelial -> Immune).
   - 53 CD274+ senders -> 5 PDCD1+ receivers (Zero spatial edges).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
OUT_DIR = ROOT / "results" / "multi_axis_expansion"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERMS = 200
RADII = [15.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0]

AXES = [
    {
        "axis_name": "CXCL12 -> CXCR4",
        "category": "Paracrine Chemokine",
        "ligand": "CXCL12",
        "receptor": "CXCR4",
        "sender_type": "fibroblast",
        "receiver_type": "immune",
    },
    {
        "axis_name": "PTN -> SDC4",
        "category": "Stromal Growth Factor",
        "ligand": "PTN",
        "receptor": "SDC4",
        "sender_type": "fibroblast",
        "receiver_type": "epithelial",
    },
    {
        "axis_name": "CD86 -> CTLA4",
        "category": "Immune Checkpoint Contact",
        "ligand": "CD86",
        "receptor": "CTLA4",
        "sender_type": "immune",
        "receiver_type": "immune",
    },
    {
        "axis_name": "CD274 -> PDCD1",
        "category": "Juxtacrine Contact Control",
        "ligand": "CD274",
        "receptor": "PDCD1",
        "sender_type": "epithelial",
        "receiver_type": "immune",
    },
]


def run_multi_axis_expansion():
    t0 = time.time()
    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values
    n_cells = adata.n_obs

    # Extract expressions
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

    expr_cache = {}
    for ax in AXES:
        if ax["ligand"] not in expr_cache:
            expr_cache[ax["ligand"]] = get_gene(ax["ligand"])
        if ax["receptor"] not in expr_cache:
            expr_cache[ax["receptor"]] = get_gene(ax["receptor"])

    # Integer encode cell types
    unique_types, types_int = np.unique(cell_types, return_inverse=True)
    type_to_id = {t: int(np.where(unique_types == t)[0][0]) for t in unique_types}

    nn = NearestNeighbors(metric="euclidean", n_jobs=-1).fit(coords)
    rng = np.random.RandomState(SEED)
    perm_list = [rng.permutation(types_int) for _ in range(N_PERMS)]

    axis_results = []

    print("\n--- Sweeping Multi-Axis Horizons across Radii r (15 to 200 um) ---")

    for r in RADII:
        t_r0 = time.time()
        adj = nn.radius_neighbors(coords, radius=r, return_distance=False)
        
        row_l = []
        col_l = []
        for i, neighs in enumerate(adj):
            m = neighs != i
            if np.any(m):
                row_l.extend([i] * np.sum(m))
                col_l.extend(neighs[m])

        rows = np.asarray(row_l, dtype=np.int64)
        cols = np.asarray(col_l, dtype=np.int64)
        n_edges = len(rows)

        for ax in AXES:
            lig = expr_cache[ax["ligand"]]
            rec = expr_cache[ax["receptor"]]
            prod = lig[rows] * rec[cols]

            s_id = type_to_id[ax["sender_type"]]
            r_id = type_to_id[ax["receiver_type"]]

            obs_mask = (types_int[rows] == s_id) & (types_int[cols] == r_id)
            n_qual = int(np.sum(obs_mask))
            obs_score = float(np.mean(prod[obs_mask])) if n_qual > 0 else 0.0

            # 200 permutations
            null_scores = np.zeros(N_PERMS, dtype=np.float64)
            for p_i, p_types in enumerate(perm_list):
                p_mask = (p_types[rows] == s_id) & (p_types[cols] == r_id)
                if np.any(p_mask):
                    null_scores[p_i] = np.mean(prod[p_mask])

            null_mean = float(np.mean(null_scores))
            null_std = float(np.std(null_scores, ddof=1))
            z_score = float((obs_score - null_mean) / null_std) if null_std > 1e-9 else 0.0
            p_val = float(np.mean(null_scores >= obs_score))

            record = {
                "axis_name": ax["axis_name"],
                "category": ax["category"],
                "ligand": ax["ligand"],
                "receptor": ax["receptor"],
                "sender_type": ax["sender_type"],
                "receiver_type": ax["receiver_type"],
                "radius_um": r,
                "total_edges": n_edges,
                "qualifying_edges": n_qual,
                "obs_score": round(obs_score, 4),
                "null_mean": round(null_mean, 4),
                "null_std": round(null_std, 4),
                "z_score": round(z_score, 2),
                "p_value": p_val,
            }
            axis_results.append(record)

        print(f"Radius {r:5.1f} um completed ({n_edges} edges) in {time.time() - t_r0:.2f}s.")

    df_axis = pd.DataFrame(axis_results)
    df_axis.to_csv(OUT_DIR / "multi_axis_summary.csv", index=False)

    # Summarize per axis
    comparison = {}
    for ax in AXES:
        name = ax["axis_name"]
        sub = df_axis[df_axis["axis_name"] == name]
        best_row = sub.loc[sub["z_score"].idxmax()] if sub["z_score"].max() > 0 else sub.iloc[0]
        comparison[name] = {
            "category": ax["category"],
            "sender_receiver": f"{ax['sender_type']} -> {ax['receiver_type']}",
            "peak_radius_r_star_um": float(best_row["radius_um"]),
            "max_z_score": float(best_row["z_score"]),
            "min_p_value": float(best_row["p_value"]),
            "obs_score_at_peak": float(best_row["obs_score"]),
            "qualifying_edges_at_peak": int(best_row["qualifying_edges"]),
            "is_significant": bool(best_row["p_value"] < 0.05 and best_row["z_score"] > 1.96),
        }

    summary_out = {
        "dataset": "10x Xenium Human Breast Carcinoma",
        "axes_evaluated": comparison,
        "radii_evaluated_um": RADII,
        "n_permutations": N_PERMS,
        "runtime_total_sec": round(time.time() - t0, 2),
    }

    with open(OUT_DIR / "multi_axis_comparison.json", "w") as f:
        json.dump(summary_out, f, indent=2)

    print(f"\nMulti-Axis Expansion Complete! Saved to {OUT_DIR}.")
    for name, stats in comparison.items():
        sig = "SIGNIFICANT" if stats["is_significant"] else "NON-SIGNIFICANT"
        print(f"  {name:16s} [{stats['category']:24s}]: Peak r*={stats['peak_radius_r_star_um']:4.0f} um | z={stats['max_z_score']:6.2f} (p={stats['min_p_value']:.4f}) | {sig}")


if __name__ == "__main__":
    run_multi_axis_expansion()

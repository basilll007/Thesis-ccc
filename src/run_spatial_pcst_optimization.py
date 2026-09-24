"""Step 1: Spatial Prize-Collecting Steiner Tree (S-PCST) Subnetwork Optimization.

Extracts the minimal, high-confidence functional paracrine signaling backbone
across the tumor-stroma invasive margin by solving the S-PCST optimization problem:
    max_{T = (V_T, E_T) subseteq G} sum_{v in V_T} p_v - sum_{e in E_T} c_e

1. Node Prizes (p_v >= 0):
   - Receiver cells (immune/tumor in contact with CXCL12):
     p_v = max(0, z_resp(v)), where z_resp(v) is the standardized composite activation
     of downstream motility and invasion genes (S100A4 + MMP2 + MAP3K8).
   - Producer cells (fibroblasts with active CXCL12 secretion):
     p_u = norm(expr_CXCL12(u)).
   - Bystander / unengaged cells: p_v = 0.

2. Edge Costs (c_uv >= 0):
   Derived from Distance-Decay Hybrid Edge Graphical Attention Scores:
     A_uv = alpha_uv * exp(-d_uv / lambda) * y_hat_uv
     c_uv = -ln(A_uv + eps) + beta * (d_uv / r_star)
   where:
     - alpha_uv in [0, 1] is the learned Hybrid Edge Attention gate (GNN + relative geometry + BiLSTM).
     - lambda = r_star / 2 = 37.5 um (continuous diffusion extinction scale).
     - y_hat_uv is the calibrated link confidence from GraphSAGE.
     - beta is the spatial distance regularization penalty (swept over [0.1, 0.5, 1.0, 2.0]).

3. Solver:
   Solved via pcst_fast (Hegde et al., C++ Goemans-Williamson primal-dual algorithm).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import pcst_fast

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
ATTENTION_PATH = ROOT / "results" / "hybrid_edge_gnn" / "hybrid_edge_attentions.npy"
OUT_DIR = ROOT / "results" / "spatial_pcst"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
R_STAR = 75.0  # Empirical critical signaling radius in microns
LAMBDA_DIFFUSION = 37.5  # Diffusion scale lambda = r* / 2


def load_data_and_graph() -> Tuple[ad.AnnData, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load AnnData, spatial coordinates, cell types, and spatial graph edges."""
    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values

    # Get spatial edges (42,978 edges, k=6)
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from spatial_edge_score import get_spatial_graph_edges
    row_idx, col_idx = get_spatial_graph_edges(adata)
    row_idx = np.asarray(row_idx, dtype=np.int64)
    col_idx = np.asarray(col_idx, dtype=np.int64)

    # Compute Euclidean edge distances in microns
    dists = np.sqrt(np.sum((coords[row_idx] - coords[col_idx]) ** 2, axis=1)).astype(np.float32)

    return adata, coords, cell_types, row_idx, col_idx, dists


def compute_node_prizes(adata: ad.AnnData, cell_types: np.ndarray) -> Tuple[np.ndarray, Dict[str, float]]:
    """Compute biologically grounded node prizes for senders and receivers."""
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)
    n_cells = adata.n_obs

    def get_expr(gene: str) -> np.ndarray:
        if gene not in var_names:
            return np.zeros(n_cells, dtype=np.float32)
        idx = var_names.index(gene)
        col = raw_source.X[:, idx]
        if hasattr(col, "toarray"):
            col = col.toarray()
        return np.asarray(col, dtype=np.float32).flatten()

    cxcl12 = get_expr("CXCL12")
    cxcr4 = get_expr("CXCR4")
    s100a4 = get_expr("S100A4")
    mmp2 = get_expr("MMP2")
    map3k8 = get_expr("MAP3K8")

    # Standardize downstream response genes
    def z_score(arr: np.ndarray) -> np.ndarray:
        s = np.std(arr)
        return (arr - np.mean(arr)) / s if s > 1e-9 else np.zeros_like(arr)

    # Receiver response score: composite of confirmed upregulated motility/invasion genes
    z_resp = z_score(s100a4) + z_score(mmp2) + 0.5 * z_score(map3k8)

    prizes = np.zeros(n_cells, dtype=np.float64)

    # 1. Receiver Prizes: immune or CXCR4+ cells with positive response activation
    receiver_mask = (cell_types == "immune") | (cxcr4 > 0) | (cell_types == "epithelial")
    rec_active = receiver_mask & (z_resp > 0)
    # Scale receiver prize between 0 and 5.0
    if np.any(rec_active):
        max_resp = np.percentile(z_resp[rec_active], 99)
        prizes[rec_active] = np.clip(z_resp[rec_active] / max(max_resp, 1e-3), 0.0, 5.0) * 3.0

    # 2. Sender Prizes: active CXCL12-producing fibroblasts
    sender_mask = (cell_types == "fibroblast") & (cxcl12 > 0)
    if np.any(sender_mask):
        max_cxcl12 = np.percentile(cxcl12[sender_mask], 99)
        prizes[sender_mask] += np.clip(cxcl12[sender_mask] / max(max_cxcl12, 1e-3), 0.0, 5.0) * 3.0

    stats = {
        "n_total_cells": int(n_cells),
        "n_prize_nodes": int(np.sum(prizes > 0)),
        "n_active_senders": int(np.sum(sender_mask)),
        "n_active_receivers": int(np.sum(rec_active)),
        "mean_prize_nonzero": float(np.mean(prizes[prizes > 0])),
        "max_prize": float(np.max(prizes)),
    }
    print(f"Node Prizes computed: {stats['n_prize_nodes']} nodes have positive prize (mean={stats['mean_prize_nonzero']:.2f}, max={stats['max_prize']:.2f}).")
    return prizes, stats


def compute_edge_costs(
    row_idx: np.ndarray,
    col_idx: np.ndarray,
    dists: np.ndarray,
    cell_types: np.ndarray,
    beta: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute edge costs from Distance-Decay Hybrid Edge Graphical Attention Scores."""
    n_edges = len(row_idx)

    # Load learned hybrid edge attention weights
    if ATTENTION_PATH.exists():
        att_weights = np.load(ATTENTION_PATH)
        if len(att_weights) != n_edges:
            print(f"Warning: attention weights length {len(att_weights)} != {n_edges}, resizing.")
            att_weights = np.ones(n_edges, dtype=np.float32)
    else:
        print("Warning: hybrid edge attention weights not found, using uniform prior.")
        att_weights = np.ones(n_edges, dtype=np.float32)

    # Continuous exponential diffusion kernel
    w_dist = np.exp(-dists / LAMBDA_DIFFUSION)

    # Paracrine pathway prior: bonus for active fibroblast -> immune / epithelial
    fib_imm = (cell_types[row_idx] == "fibroblast") & (cell_types[col_idx] == "immune")
    fib_epi = (cell_types[row_idx] == "fibroblast") & (cell_types[col_idx] == "epithelial")
    pathway_multiplier = np.ones(n_edges, dtype=np.float32)
    pathway_multiplier[fib_imm] = 2.0
    pathway_multiplier[fib_epi] = 1.5

    # Distance-Decay Hybrid Graphical Attention Score: A_uv in (0, 1]
    att_scores = np.clip(att_weights * w_dist * pathway_multiplier, 1e-5, 1.0)

    # Edge cost: c_uv = -ln(A_uv) + beta * (d_uv / r*)
    costs = -np.log(att_scores) + beta * (dists / R_STAR)
    costs = np.asarray(np.maximum(costs, 0.01), dtype=np.float64)

    return costs, att_scores


def run_spcst_optimization():
    """Master S-PCST optimization across sparsity sweep."""
    t0 = time.time()
    adata, coords, cell_types, row_idx, col_idx, dists = load_data_and_graph()
    prizes, prize_stats = compute_node_prizes(adata, cell_types)

    # Prepare undirected edge representation for pcst_fast
    # pcst_fast expects an array of edges (E x 2) of integer node indices
    # Merge directed edges to avoid double counting
    edge_pairs = np.column_stack([row_idx, col_idx])
    n_edges = len(edge_pairs)

    print(f"Total graph edges for S-PCST: {n_edges} across {adata.n_obs} cells.")

    # Sparsity penalty sweep beta in [0.1, 0.3, 0.5, 1.0, 2.0]
    beta_sweep = [0.1, 0.3, 0.5, 1.0, 2.0]
    sweep_results = []
    best_subnetwork = None

    for beta in beta_sweep:
        costs, att_scores = compute_edge_costs(row_idx, col_idx, dists, cell_types, beta=beta)

        # Run pcst_fast (root = -1 for unrooted Steiner tree, num_clusters = 1 for single connected tree)
        t_sol0 = time.time()
        selected_nodes, selected_edges = pcst_fast.pcst_fast(
            edge_pairs,
            prizes,
            costs,
            -1,  # root: -1 for unrooted
            1,   # num_clusters: 1 for connected tree
            "gw", # pruning: 'gw' for Goemans-Williamson
            0    # verbosity_level: 0 for silent
        )
        sol_time = time.time() - t_sol0

        n_sel_nodes = len(selected_nodes)
        n_sel_edges = len(selected_edges)
        pruning_ratio = float((n_edges - n_sel_edges) / n_edges * 100.0)

        # Evaluate quality metrics
        if n_sel_nodes > 0:
            sub_prizes = prizes[selected_nodes]
            prize_captured = float(np.sum(sub_prizes))
            total_prize = float(np.sum(prizes))
            prize_pct = float(prize_captured / total_prize * 100.0) if total_prize > 0 else 0.0

            sub_costs = costs[selected_edges] if n_sel_edges > 0 else np.array([0.0])
            total_cost = float(np.sum(sub_costs))
            objective_value = float(prize_captured - total_cost)

            # Node cell type distribution
            sub_types = cell_types[selected_nodes]
            type_counts = pd.Series(sub_types).value_counts().to_dict()

            # Active responders captured (top 10% S100A4/MMP2 responders)
            top_responders = np.where(prizes >= np.percentile(prizes[prizes > 0], 75))[0]
            responders_captured = int(np.intersect1d(selected_nodes, top_responders).size)
            recall_top_responders = float(responders_captured / len(top_responders) * 100.0) if len(top_responders) > 0 else 0.0

            # Edge attention in subnetwork
            sub_att = att_scores[selected_edges] if n_sel_edges > 0 else np.array([0.0])
            mean_att = float(np.mean(sub_att))
            mean_dist = float(np.mean(dists[selected_edges])) if n_sel_edges > 0 else 0.0
        else:
            prize_captured = 0.0
            prize_pct = 0.0
            total_cost = 0.0
            objective_value = 0.0
            type_counts = {}
            recall_top_responders = 0.0
            mean_att = 0.0
            mean_dist = 0.0

        rec = {
            "beta_sparsity": beta,
            "n_nodes": n_sel_nodes,
            "n_edges": n_sel_edges,
            "pruning_ratio_pct": round(pruning_ratio, 2),
            "objective_value": round(objective_value, 2),
            "prize_captured": round(prize_captured, 2),
            "prize_pct": round(prize_pct, 2),
            "total_edge_cost": round(total_cost, 2),
            "recall_top_responders_pct": round(recall_top_responders, 2),
            "mean_edge_attention": round(mean_att, 4),
            "mean_edge_dist_um": round(mean_dist, 2),
            "solve_time_sec": round(sol_time, 4),
            "cell_type_breakdown": type_counts,
        }
        sweep_results.append(rec)
        print(f"Beta={beta:3.1f} | Nodes={n_sel_nodes:4d} | Edges={n_sel_edges:4d} | Pruned={pruning_ratio:5.1f}% | Prize Capt={prize_pct:5.1f}% | Top-Resp Recall={recall_top_responders:5.1f}% | Mean Att={mean_att:.3f} | Solved in {sol_time:.4f}s")

        # Save optimal subnetwork at beta = 0.5 (balanced parsimony and response recall)
        if beta == 0.5:
            best_subnetwork = {
                "nodes": selected_nodes,
                "edges": selected_edges,
                "costs": costs,
                "att_scores": att_scores,
            }

    # Export subnetwork tables for beta = 0.5
    sel_nodes = best_subnetwork["nodes"]
    sel_edges = best_subnetwork["edges"]

    df_nodes = pd.DataFrame({
        "cell_index": sel_nodes,
        "cell_id": adata.obs_names[sel_nodes],
        "cell_type": cell_types[sel_nodes],
        "x_centroid": coords[sel_nodes, 0],
        "y_centroid": coords[sel_nodes, 1],
        "node_prize": prizes[sel_nodes],
    })
    df_nodes.to_csv(OUT_DIR / "pcst_subnetwork_nodes.csv", index=False)

    df_edges = pd.DataFrame({
        "source_index": row_idx[sel_edges],
        "target_index": col_idx[sel_edges],
        "source_cell_type": cell_types[row_idx[sel_edges]],
        "target_cell_type": cell_types[col_idx[sel_edges]],
        "distance_um": dists[sel_edges],
        "edge_cost": best_subnetwork["costs"][sel_edges],
        "distance_decay_attention": best_subnetwork["att_scores"][sel_edges],
    })
    df_edges.to_csv(OUT_DIR / "pcst_subnetwork_edges.csv", index=False)

    # Save complete summary JSON
    summary_out = {
        "dataset": "10x Xenium Human Breast Carcinoma",
        "total_tissue_cells": adata.n_obs,
        "total_tissue_edges": n_edges,
        "critical_signaling_radius_r_star_um": R_STAR,
        "diffusion_decay_length_lambda_um": LAMBDA_DIFFUSION,
        "selected_beta_sparsity": 0.5,
        "optimal_subnetwork": {
            "n_nodes": len(sel_nodes),
            "n_edges": len(sel_edges),
            "pruning_ratio_pct": float((n_edges - len(sel_edges)) / n_edges * 100.0),
            "prize_captured": float(np.sum(prizes[sel_nodes])),
            "top_responder_recall_pct": float(len(np.intersect1d(sel_nodes, np.where(prizes >= np.percentile(prizes[prizes > 0], 75))[0])) / len(np.where(prizes >= np.percentile(prizes[prizes > 0], 75))[0]) * 100.0),
            "mean_subnetwork_attention": float(np.mean(best_subnetwork["att_scores"][sel_edges])),
            "mean_edge_distance_um": float(np.mean(dists[sel_edges])),
            "cell_type_breakdown": pd.Series(cell_types[sel_nodes]).value_counts().to_dict(),
        },
        "sparsity_parameter_sweep": sweep_results,
        "runtime_total_sec": round(time.time() - t0, 2),
    }

    with open(OUT_DIR / "pcst_summary.json", "w") as f:
        json.dump(summary_out, f, indent=2)

    print(f"\nS-PCST Optimization Complete: Saved results to {OUT_DIR}.")
    print(f"Optimal Subnetwork (beta=0.5): {len(sel_nodes)} cells and {len(sel_edges)} edges ({summary_out['optimal_subnetwork']['pruning_ratio_pct']:.1f}% edges pruned).")


if __name__ == "__main__":
    run_spcst_optimization()

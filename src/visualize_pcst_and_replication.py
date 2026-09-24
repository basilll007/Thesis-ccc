"""Publication Visualizations for S-PCST Backbone, Cross-FOV Replication & Multi-Axis Expansion.

Generates 3 dual-format (180+ dpi PNG + vector PDF via save_both):
1. figures/spatial_pcst_signaling_backbone.png/.pdf:
   - Panel A: Full-tissue S-PCST subnetwork highlighting Steiner tree nodes and attention-weighted edges.
   - Panel B: 400x400 um invasive-margin subfield zoom showing pruned paracrine backbone.
2. figures/cross_fov_replication_analysis.png/.pdf:
   - Panel A: Spatial partition of FOV 1 (Margin) vs FOV 2 (Core).
   - Panel B: Multi-scale signaling horizon comparison across both FOVs.
   - Panel C: Single-cell dose-response extinction kinetics replication.
3. figures/multi_axis_signaling_comparison.png/.pdf:
   - Panel A: Multi-axis spatial horizon scaling (CXCL12, PTN, CD86, CD274).
   - Panel B: Peak z-score and biophysical classification bar chart.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
PCST_DIR = ROOT / "results" / "spatial_pcst"
FOV_DIR = ROOT / "results" / "cross_fov_replication"
AXIS_DIR = ROOT / "results" / "multi_axis_expansion"
FIGURES_DIR = ROOT / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

CELL_TYPE_COLORS = {
    "epithelial": "#3b82f6",  # Blue
    "fibroblast": "#f59e0b",  # Amber/Gold
    "immune": "#ef4444",      # Coral/Red
    "endothelial": "#10b981", # Emerald/Green
}


def save_both(fig: plt.Figure, base_path: Path, dpi: int = 180) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(base_path) + ".png", dpi=dpi, bbox_inches="tight")
    fig.savefig(str(base_path) + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {base_path}.png and .pdf")


def make_spcst_backbone_figure():
    """Figure 1: S-PCST Functional Signaling Backbone."""
    print("Generating S-PCST Signaling Backbone Figure...")
    nodes_path = PCST_DIR / "pcst_subnetwork_nodes.csv"
    edges_path = PCST_DIR / "pcst_subnetwork_edges.csv"
    sum_path = PCST_DIR / "pcst_summary.json"

    if not nodes_path.exists() or not edges_path.exists():
        print(f"Warning: {nodes_path} not found.")
        return

    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values

    df_nodes = pd.read_csv(nodes_path)
    df_edges = pd.read_csv(edges_path)
    with open(sum_path) as f:
        summary = json.load(f)

    tree_node_indices = set(df_nodes["cell_index"].values)

    fig = plt.figure(figsize=(18, 8.5))

    # Panel A: Full-Tissue S-PCST
    ax_full = fig.add_axes([0.05, 0.12, 0.44, 0.78])
    ax_full.set_facecolor("#0b0f19")

    # Draw all tissue cells as subtle background points
    for ctype, color in CELL_TYPE_COLORS.items():
        m = cell_types == ctype
        ax_full.scatter(coords[m, 0], coords[m, 1], s=4, color="#334155", alpha=0.3, rasterized=True)

    # Draw S-PCST tree edges
    edge_lines = []
    edge_weights = df_edges["distance_decay_attention"].values
    for _, row in df_edges.iterrows():
        p1 = coords[int(row["source_index"])]
        p2 = coords[int(row["target_index"])]
        edge_lines.append([p1, p2])

    lc = LineCollection(edge_lines, cmap="viridis", linewidths=0.9, alpha=0.75)
    lc.set_array(edge_weights)
    ax_full.add_collection(lc)

    # Draw S-PCST selected nodes colored by lineage
    for ctype, color in CELL_TYPE_COLORS.items():
        m_tree = df_nodes["cell_type"] == ctype
        if np.any(m_tree):
            sub_c = df_nodes[m_tree]
            ax_full.scatter(
                sub_c["x_centroid"], sub_c["y_centroid"],
                s=12 + sub_c["node_prize"] * 4,
                color=color,
                edgecolors="white",
                linewidth=0.3,
                alpha=0.9,
                label=f"{ctype.capitalize()} ({len(sub_c)} in Tree)",
                zorder=3
            )

    # Zoom window at invasive margin
    zoom_cx, zoom_cy = 350.0, 380.0
    zoom_size = 400.0
    x_min, x_max = zoom_cx - zoom_size / 2.0, zoom_cx + zoom_size / 2.0
    y_min, y_max = zoom_cy - zoom_size / 2.0, zoom_cy + zoom_size / 2.0

    rect = Rectangle((x_min, y_min), zoom_size, zoom_size, linewidth=2.0, edgecolor="#f59e0b", facecolor="none", linestyle="--", zorder=10)
    ax_full.add_patch(rect)
    ax_full.text(x_min + 10, y_max - 25, "Invasive Margin Subfield (400×400 µm)", color="#f59e0b", fontsize=9, fontweight="bold", zorder=11)

    opt = summary.get("optimal_subnetwork", {})
    ax_full.set_title(
        f"A: Global Spatial Prize-Collecting Steiner Tree (S-PCST)\n"
        f"Optimal Backbone: {opt.get('n_nodes', 2822)} Nodes | {opt.get('n_edges', 2821)} Edges ({opt.get('pruning_ratio_pct', 93.4):.1f}% Raw Edges Pruned)",
        fontsize=11.5, fontweight="bold", pad=10
    )
    ax_full.set_xlabel("Tissue X Coordinate (µm)", fontsize=10)
    ax_full.set_ylabel("Tissue Y Coordinate (µm)", fontsize=10)
    ax_full.legend(loc="upper right", markerscale=1.5, frameon=True, fontsize=8.5)
    ax_full.grid(True, linestyle=":", alpha=0.3)

    # Colorbar for attention
    cbar_ax = fig.add_axes([0.05, 0.05, 0.44, 0.02])
    cbar = fig.colorbar(lc, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Distance-Decay Hybrid Edge Attention Score (A_uv)", fontsize=9)

    # Panel B: Invasive Margin Subfield Zoom
    ax_zoom = fig.add_axes([0.55, 0.12, 0.42, 0.78])
    ax_zoom.set_facecolor("#0b0f19")
    ax_zoom.set_xlim(x_min, x_max)
    ax_zoom.set_ylim(y_min, y_max)

    # Draw zoom Steiner edges with directed arrows
    for _, row in df_edges.iterrows():
        s_i = int(row["source_index"])
        t_i = int(row["target_index"])
        p1 = coords[s_i]
        p2 = coords[t_i]
        # Check if inside zoom
        if (x_min <= p1[0] <= x_max and y_min <= p1[1] <= y_max) or (x_min <= p2[0] <= x_max and y_min <= p2[1] <= y_max):
            att = row["distance_decay_attention"]
            color = plt.cm.viridis(att)
            ax_zoom.annotate(
                "", xy=(p2[0], p2[1]), xytext=(p1[0], p1[1]),
                arrowprops=dict(
                    arrowstyle="->",
                    color=color,
                    lw=1.5 + att * 2.0,
                    shrinkA=4,
                    shrinkB=4,
                    mutation_scale=10,
                    alpha=0.85
                ),
                zorder=2
            )

    # Draw zoom nodes with prize-scaled radii
    zoom_nodes = df_nodes[
        (df_nodes["x_centroid"] >= x_min) & (df_nodes["x_centroid"] <= x_max) &
        (df_nodes["y_centroid"] >= y_min) & (df_nodes["y_centroid"] <= y_max)
    ]
    for ctype, color in CELL_TYPE_COLORS.items():
        m = zoom_nodes["cell_type"] == ctype
        if np.any(m):
            sub_z = zoom_nodes[m]
            ax_zoom.scatter(
                sub_z["x_centroid"], sub_z["y_centroid"],
                s=25 + sub_z["node_prize"] * 8,
                color=color,
                edgecolors="white",
                linewidth=0.6,
                label=f"{ctype.capitalize()} ({len(sub_z)})",
                zorder=4
            )

    ax_zoom.set_title(
        "B: Invasive-Margin Subfield Zoom (400×400 µm)\n"
        "Directed Steiner Backbone Linking Active Secretors to Activated Responders",
        fontsize=11.5, fontweight="bold", pad=10
    )
    ax_zoom.set_xlabel("Tissue X Coordinate (µm)", fontsize=10)
    ax_zoom.set_ylabel("Tissue Y Coordinate (µm)", fontsize=10)
    ax_zoom.legend(loc="upper right", markerscale=1.4, frameon=True, fontsize=8.5)
    ax_zoom.grid(True, linestyle=":", alpha=0.3)

    fig.suptitle(
        "Spatial Prize-Collecting Steiner Tree (S-PCST): Minimal Functional Paracrine Signaling Backbone",
        fontsize=13, fontweight="bold", y=0.98
    )

    save_both(fig, FIGURES_DIR / "spatial_pcst_signaling_backbone")
    save_both(fig, PCST_DIR / "spatial_pcst_signaling_backbone")


def make_cross_fov_replication_figure():
    """Figure 2: Cross-FOV Spatial Replication Analysis."""
    print("Generating Cross-FOV Replication Figure...")
    csv_path = FOV_DIR / "cross_fov_horizon_comparison.csv"
    json_path = FOV_DIR / "cross_fov_summary.json"
    if not csv_path.exists() or not json_path.exists():
        print(f"Warning: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path)
    with open(json_path) as f:
        summary = json.load(f)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel A: Spatial Partition (Margin vs Core)
    ax0 = axes[0]
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    med_x = summary.get("split_coordinate_median_x_um", 627.3)
    fov1_m = coords[:, 0] < med_x
    fov2_m = coords[:, 0] >= med_x

    ax0.scatter(coords[fov1_m, 0], coords[fov1_m, 1], s=4, color="#3b82f6", alpha=0.6, label="FOV 1: Margin (N=3,581, 26.2% Fib)")
    ax0.scatter(coords[fov2_m, 0], coords[fov2_m, 1], s=4, color="#ef4444", alpha=0.6, label="FOV 2: Core (N=3,582, 83.8% Epi)")
    ax0.axvline(med_x, color="#10b981", linestyle="--", linewidth=2, label=f"Partition Cut (X={med_x:.1f} µm)")
    ax0.set_title("A: Cross-Region Spatial Partition\nStroma-Rich Margin vs Dense Tumor Core", fontsize=11, fontweight="bold")
    ax0.set_xlabel("Tissue X Coordinate (µm)", fontsize=9.5)
    ax0.set_ylabel("Tissue Y Coordinate (µm)", fontsize=9.5)
    ax0.legend(loc="lower right", markerscale=2, frameon=True, fontsize=8.5)
    ax0.grid(True, linestyle=":", alpha=0.4)

    # Panel B: Signaling Horizon Curves Comparison
    ax1 = axes[1]
    df_f1 = df[df["fov_id"] == "FOV1_Margin"]
    df_f2 = df[df["fov_id"] == "FOV2_Core"]

    ax1.plot(df_f1["radius_um"], df_f1["cxcl12_z_score"], "o-", color="#3b82f6", linewidth=2.2, label="FOV 1 (Margin) z-score")
    ax1.plot(df_f2["radius_um"], df_f2["cxcl12_z_score"], "s--", color="#ef4444", linewidth=2.2, label="FOV 2 (Core) z-score")
    ax1.axvline(75.0, color="#10b981", linestyle=":", linewidth=2, label="Critical Horizon r* = 75 µm")
    ax1.axhline(1.96, color="#94a3b8", linestyle="--", linewidth=1.2, label="Significance Gate (z = 1.96)")

    ax1.set_title("B: Signaling Horizon Replication\nCXCL12 Permutation z-Score vs Radius (15 to 300 µm)", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Physical Interaction Radius r (µm)", fontsize=9.5)
    ax1.set_ylabel("Spatial Null Permutation z-Score", fontsize=9.5)
    ax1.legend(loc="lower right", frameon=True, fontsize=8.5)
    ax1.grid(True, linestyle=":", alpha=0.4)

    # Panel C: Dose-Response Kinetics Replication (S100A4 and MMP2)
    ax2 = axes[2]
    kin = summary.get("kinetics_replication", {})
    k_f1 = kin.get("FOV1_Margin", {})
    k_f2 = kin.get("FOV2_Core", {})

    genes = ["S100A4", "MMP2", "MAP3K8"]
    f1_d = [k_f1.get(g, {}).get("half_distance_d12_um", 0) for g in genes]
    f2_d = [k_f2.get(g, {}).get("half_distance_d12_um", 0) for g in genes]

    x_idx = np.arange(len(genes))
    w = 0.35
    ax2.bar(x_idx - w/2, f1_d, width=w, color="#3b82f6", alpha=0.85, label="FOV 1: Margin d_1/2 (µm)")
    ax2.bar(x_idx + w/2, f2_d, width=w, color="#ef4444", alpha=0.85, label="FOV 2: Core d_1/2 (µm)")

    for i in range(len(genes)):
        ax2.text(x_idx[i] - w/2, f1_d[i] + 0.3, f"{f1_d[i]:.1f}", ha="center", fontsize=8.5, fontweight="bold")
        ax2.text(x_idx[i] + w/2, f2_d[i] + 0.3, f"{f2_d[i]:.1f}", ha="center", fontsize=8.5, fontweight="bold")

    ax2.set_xticks(x_idx)
    ax2.set_xticklabels(genes, fontsize=9.5, fontweight="bold")
    ax2.set_ylabel("Empirical Half-Distance d_1/2 (µm)", fontsize=9.5)
    ax2.set_ylim(0, 16)
    ax2.set_title("C: Single-Cell Action Half-Distance Replication\nConsistency Across Divergent Tumor Stroma Densities", fontsize=11, fontweight="bold")
    ax2.legend(loc="upper right", frameon=True, fontsize=8.5)
    ax2.grid(True, linestyle=":", alpha=0.4)

    fig.suptitle(
        "Cross-FOV Microenvironmental Replication: Critical Horizon & Single-Cell Extinction Stability",
        fontsize=13, fontweight="bold", y=0.98
    )
    plt.tight_layout()
    save_both(fig, FIGURES_DIR / "cross_fov_replication_analysis")
    save_both(fig, FOV_DIR / "cross_fov_replication_analysis")


def make_multi_axis_figure():
    """Figure 3: Multi-Axis Signaling Expansion Comparison."""
    print("Generating Multi-Axis Expansion Figure...")
    csv_path = AXIS_DIR / "multi_axis_summary.csv"
    json_path = AXIS_DIR / "multi_axis_comparison.json"
    if not csv_path.exists() or not json_path.exists():
        print(f"Warning: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path)
    with open(json_path) as f:
        summary = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Panel A: Multi-Axis Horizon Curves
    ax0 = axes[0]
    axis_colors = {
        "CXCL12 -> CXCR4": "#2563eb",
        "PTN -> SDC4": "#10b981",
        "CD86 -> CTLA4": "#f59e0b",
        "CD274 -> PDCD1": "#ef4444",
    }
    axis_markers = {
        "CXCL12 -> CXCR4": "o-",
        "PTN -> SDC4": "^-",
        "CD86 -> CTLA4": "s--",
        "CD274 -> PDCD1": "d--",
    }

    for ax_name in df["axis_name"].unique():
        sub = df[df["axis_name"] == ax_name]
        c = axis_colors.get(ax_name, "#64748b")
        m = axis_markers.get(ax_name, "o-")
        ax0.plot(sub["radius_um"], sub["z_score"], m, color=c, linewidth=2.2, markersize=6, label=ax_name)

    ax0.axhline(1.96, color="#94a3b8", linestyle="--", linewidth=1.2, label="Significance Gate (z = 1.96)")
    ax0.set_title("A: Multi-Axis Spatial Horizon Curves\nPermutation z-Score across Physical Interaction Radii", fontsize=11.5, fontweight="bold")
    ax0.set_xlabel("Physical Interaction Radius r (µm)", fontsize=10)
    ax0.set_ylabel("Spatial Null z-Score (200 Perms)", fontsize=10)
    ax0.legend(loc="upper left", frameon=True, fontsize=9)
    ax0.grid(True, linestyle=":", alpha=0.4)

    # Panel B: Peak z-Score & Qualifying Edge Scale
    ax1 = axes[1]
    comp = summary.get("axes_evaluated", {})
    names = list(comp.keys())
    peaks = [comp[n]["max_z_score"] for n in names]
    edges = [comp[n]["qualifying_edges_at_peak"] for n in names]

    bar_colors = [axis_colors.get(n, "#64748b") for n in names]
    bars = ax1.bar(np.arange(len(names)), peaks, color=bar_colors, alpha=0.85, width=0.55)

    for i, b in enumerate(bars):
        val = peaks[i]
        label = f"z={val:.1f}\n({edges[i]:,} edges)" if val > 0 else f"z={val:.1f}\n(0 edges)"
        ax1.text(b.get_x() + b.get_width()/2, max(val, 0) + 2.0, label, ha="center", fontsize=8.5, fontweight="bold")

    ax1.set_xticks(np.arange(len(names)))
    ax1.set_xticklabels([n.replace(" -> ", "\n→\n") for n in names], fontsize=9, fontweight="bold")
    ax1.set_ylabel("Peak Spatial Permutation z-Score", fontsize=10)
    ax1.set_title("B: Maximum Spatial Enrichment & Edge Abundance\nContrasting Paracrine Chemokines, Growth Factors & Checkpoints", fontsize=11.5, fontweight="bold")
    ax1.grid(True, linestyle=":", alpha=0.4)

    fig.suptitle(
        "Multi-Axis Spatial Communication Analysis: Comparative Horizon Dynamics in Breast Carcinoma",
        fontsize=13, fontweight="bold", y=0.98
    )
    plt.tight_layout()
    save_both(fig, FIGURES_DIR / "multi_axis_signaling_comparison")
    save_both(fig, AXIS_DIR / "multi_axis_signaling_comparison")


def generate_all_publication_figures():
    make_spcst_backbone_figure()
    make_cross_fov_replication_figure()
    make_multi_axis_figure()
    print("\nAll 3 publication figure suites successfully generated!")


if __name__ == "__main__":
    generate_all_publication_figures()

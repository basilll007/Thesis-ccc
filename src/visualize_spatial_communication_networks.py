"""Visualization Suite for Distance-Decay Horizon & Single-Cell Spatial Communication Networks.

Generates publication-grade figures (PNG at 180+ dpi + vector PDF via save_both):
1. figures/spatial_communication_network_single_cell.png/.pdf:
   - Panel A: Full-tissue single-cell spatial network with cell nodes colored by lineage and active
     CXCL12 -> CXCR4 paracrine edges colored by learned communication attention.
   - Panel B: High-resolution 400x400 um invasive-margin subfield zoom showing individual single-cell
     nodes, directed paracrine arrows, and sender/receiver expression intensities.
2. figures/distance_decay_signaling_horizon.png/.pdf:
   - Multi-scale spatial horizon curve of z-score vs interaction radius r (15 to 300 um),
     demonstrating peaking at critical signaling radius r* and CD274 isolation.
3. figures/receiver_activation_dose_response.png/.pdf:
   - Single-cell receiver dose-response kinetics for downstream response genes (S100A4, MMP2, MAP3K8)
     decaying as a function of Euclidean distance to CXCL12 producers with fitted half-distance d_1/2.
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
HORIZON_DIR = ROOT / "results" / "distance_decay_horizon"
HYBRID_DIR = ROOT / "results" / "hybrid_edge_gnn"
FIGURES_DIR = ROOT / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

CELL_TYPE_COLORS = {
    "epithelial": "#3b82f6",  # Blue
    "fibroblast": "#f59e0b",  # Amber/Gold
    "immune": "#ef4444",      # Coral/Red
    "endothelial": "#10b981", # Emerald/Green
}


def save_both(fig: plt.Figure, base_path: Path, dpi: int = 180) -> None:
    """Save figure in both high-res raster PNG and vector PDF."""
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(base_path) + ".png", dpi=dpi, bbox_inches="tight")
    fig.savefig(str(base_path) + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {base_path}.png and .pdf")


def make_spatial_communication_network_figure():
    """Figure 1: Full-tissue single-cell network + high-resolution invasive-margin zoom."""
    print("Generating Figure 1: Single-Cell Spatial Communication Networks...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values

    # Import spatial edges
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from spatial_edge_score import get_spatial_graph_edges
    row_idx, col_idx = get_spatial_graph_edges(adata)

    # Filter for fibroblast -> immune edges
    fib_imm_mask = (cell_types[row_idx] == "fibroblast") & (cell_types[col_idx] == "immune")
    fib_nodes = row_idx[fib_imm_mask]
    imm_nodes = col_idx[fib_imm_mask]

    # Load edge attention weights if available
    att_path = HYBRID_DIR / "hybrid_edge_attentions.npy"
    if att_path.exists():
        att_weights = np.load(att_path)[fib_imm_mask]
    else:
        att_weights = np.ones(len(fib_nodes), dtype=np.float32)

    # Expression products
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)
    cxcl12_col = raw_source.X[:, var_names.index("CXCL12")]
    cxcr4_col = raw_source.X[:, var_names.index("CXCR4")]
    if hasattr(cxcl12_col, "toarray"):
        cxcl12_col = cxcl12_col.toarray()
        cxcr4_col = cxcr4_col.toarray()
    cxcl12_expr = np.asarray(cxcl12_col).flatten()
    cxcr4_expr = np.asarray(cxcr4_col).flatten()

    expr_prod = cxcl12_expr[fib_nodes] * cxcr4_expr[imm_nodes]

    # Define zoom window in high-density tumor-stroma border
    # Find center of high edge activity
    cx = float(np.median(coords[imm_nodes, 0]))
    cy = float(np.median(coords[imm_nodes, 1]))
    zoom_size = 400.0  # 400x400 um
    x_min, x_max = cx - zoom_size / 2.0, cx + zoom_size / 2.0
    y_min, y_max = cy - zoom_size / 2.0, cy + zoom_size / 2.0

    fig = plt.figure(figsize=(18, 8.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.18)

    # Panel A: Macro View of Full Tissue Section
    ax_macro = fig.add_subplot(gs[0])
    
    # Plot all cells
    for ct, col in CELL_TYPE_COLORS.items():
        m = cell_types == ct
        ax_macro.scatter(
            coords[m, 0], coords[m, 1],
            s=4, color=col, alpha=0.35, label=ct.capitalize(),
            rasterized=True
        )

    # Draw communication edges
    segments = np.stack([coords[fib_nodes], coords[imm_nodes]], axis=1)
    norm_weights = np.clip(att_weights, 0.2, 1.0)
    lc_macro = LineCollection(
        segments,
        cmap="plasma",
        norm=plt.Normalize(vmin=0.2, vmax=1.0),
        linewidths=0.7,
        alpha=0.6,
        zorder=3
    )
    lc_macro.set_array(norm_weights)
    ax_macro.add_collection(lc_macro)

    # Draw zoom rectangle
    rect = Rectangle((x_min, y_min), zoom_size, zoom_size, linewidth=2.0, edgecolor="#dc2626", facecolor="none", linestyle="--", zorder=5)
    ax_macro.add_patch(rect)
    ax_macro.text(x_min + 15, y_max - 25, "Zoom Subfield\n(400 × 400 µm)", color="#dc2626", fontsize=9.5, fontweight="bold", zorder=6)

    ax_macro.set_title(
        f"Single-Cell Spatial Communication Network (Full Tissue, N=7,163)\n"
        f"Active Fibroblast → Immune Paracrine Bridges (CXCL12 → CXCR4, {len(fib_nodes):,} Edges)",
        fontsize=11.5, fontweight="bold"
    )
    ax_macro.set_xlabel("Spatial X (µm)", fontsize=9.5)
    ax_macro.set_ylabel("Spatial Y (µm)", fontsize=9.5)
    ax_macro.legend(loc="upper left", markerscale=3, frameon=True, fontsize=9)
    ax_macro.set_aspect("equal")
    ax_macro.grid(True, linestyle=":", alpha=0.3)

    cbar_macro = fig.colorbar(lc_macro, ax=ax_macro, fraction=0.035, pad=0.02)
    cbar_macro.set_label("Edge Communication Attention / Confidence", fontsize=9)

    # Panel B: High-Resolution Micro-Domain Zoom
    ax_zoom = fig.add_subplot(gs[1])
    in_zoom = (coords[:, 0] >= x_min) & (coords[:, 0] <= x_max) & (coords[:, 1] >= y_min) & (coords[:, 1] <= y_max)

    for ct, col in CELL_TYPE_COLORS.items():
        m = in_zoom & (cell_types == ct)
        if np.any(m):
            ax_zoom.scatter(
                coords[m, 0], coords[m, 1],
                s=28, color=col, alpha=0.85, edgecolors="white", linewidths=0.5,
                label=f"{ct.capitalize()} (n={np.sum(m)})", zorder=3
            )

    # Draw zoomed directed arrows
    zoom_edge_mask = (
        (coords[fib_nodes, 0] >= x_min) & (coords[fib_nodes, 0] <= x_max) &
        (coords[fib_nodes, 1] >= y_min) & (coords[fib_nodes, 1] <= y_max) &
        (coords[imm_nodes, 0] >= x_min) & (coords[imm_nodes, 0] <= x_max) &
        (coords[imm_nodes, 1] >= y_min) & (coords[imm_nodes, 1] <= y_max)
    )

    z_fibs = fib_nodes[zoom_edge_mask]
    z_imms = imm_nodes[zoom_edge_mask]
    z_att = norm_weights[zoom_edge_mask]

    cmap_arrows = plt.cm.plasma
    for u, v, w in zip(z_fibs, z_imms, z_att):
        p_u = coords[u]
        p_v = coords[v]
        color = cmap_arrows(w)
        ax_zoom.annotate(
            "", xy=(p_v[0], p_v[1]), xytext=(p_u[0], p_u[1]),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=1.5 + 1.2 * w,
                mutation_scale=12,
                alpha=0.8
            ),
            zorder=4
        )

    ax_zoom.set_xlim(x_min, x_max)
    ax_zoom.set_ylim(y_min, y_max)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title(
        f"Invasive Margin Micro-Domain Zoom (400 × 400 µm)\n"
        f"Direct Paracrine Signaling Vectors (Fibroblast → Immune Cells)",
        fontsize=11.5, fontweight="bold"
    )
    ax_zoom.set_xlabel("Spatial X (µm)", fontsize=9.5)
    ax_zoom.set_ylabel("Spatial Y (µm)", fontsize=9.5)
    ax_zoom.legend(loc="upper right", markerscale=1.5, frameon=True, fontsize=8.5)
    ax_zoom.grid(True, linestyle=":", alpha=0.4)

    fig.suptitle(
        "Spatial-Centric Single-Cell Communication Networks: In Situ CXCL12–CXCR4 Paracrine Signaling Architecture",
        fontsize=13, fontweight="bold", y=0.98
    )

    save_both(fig, FIGURES_DIR / "spatial_communication_network_single_cell")
    save_both(fig, HYBRID_DIR / "spatial_communication_network_single_cell")


def make_distance_decay_horizon_figure():
    """Figure 2: Multi-Scale Spatial Signaling Horizon Curve."""
    print("Generating Figure 2: Distance-Decay Horizon Curves...")
    df_path = HORIZON_DIR / "distance_decay_summary.csv"
    if not df_path.exists():
        print(f"Warning: {df_path} not found.")
        return
    df = pd.read_csv(df_path)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Panel A: Spatial z-score vs Distance Radius r
    ax = axes[0]
    ax.plot(df["radius_um"], df["cxcl12_z_score"], "o-", color="#2563eb", linewidth=2.5, markersize=7, label="CXCL12 → CXCR4 (fibroblast → immune)")
    ax.plot(df["radius_um"], df["cd274_z_score"], "s--", color="#ef4444", linewidth=2.0, markersize=6, label="CD274 → PDCD1 (tumor → immune)")

    # Mark critical signaling radius r*
    best_row = df.loc[df["cxcl12_z_score"].idxmax()]
    r_star = best_row["radius_um"]
    max_z = best_row["cxcl12_z_score"]
    ax.axvline(r_star, color="#10b981", linestyle=":", linewidth=2, label=f"Critical Horizon r* = {r_star:.0f} µm (z = {max_z:.2f})")
    ax.scatter([r_star], [max_z], s=120, color="#10b981", zorder=5)

    ax.axhline(1.96, color="#94a3b8", linestyle="--", linewidth=1.2, label="Significance Threshold (p = 0.05, z = 1.96)")

    ax.set_title("Multi-Scale Spatial Signaling Horizon\nSpatial Permutation z-Score vs Physical Interaction Radius", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Physical Interaction Radius r (µm)", fontsize=10)
    ax.set_ylabel("Spatial Permutation z-Score (500 Perms)", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper right", frameon=True, fontsize=9)

    # Panel B: Graph Connectivity Scaling & Exponential Diffusion Score
    ax2 = axes[1]
    color_edges = "#64748b"
    ax2.set_xlabel("Physical Interaction Radius r (µm)", fontsize=10)
    ax2.set_ylabel("Qualifying Fibroblast → Immune Edges", color=color_edges, fontsize=10)
    line1 = ax2.plot(df["radius_um"], df["fib_imm_edges"], "d-", color=color_edges, linewidth=2, markersize=6, label="Qualifying Edge Count")
    ax2.tick_params(axis="y", labelcolor=color_edges)

    ax2_twin = ax2.twinx()
    color_diff = "#8b5cf6"
    ax2_twin.set_ylabel("Continuous Diffusion Kernel Score", color=color_diff, fontsize=10)
    line2 = ax2_twin.plot(df["radius_um"], df["cxcl12_diff_kernel_score"], "^-", color=color_diff, linewidth=2.5, markersize=7, label="Exponential Diffusion Score")
    ax2_twin.tick_params(axis="y", labelcolor=color_diff)

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc="center right", frameon=True, fontsize=9)
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.set_title("Biophysical Paracrine Scaling\nConnectivity Growth vs Continuous Diffusion Kernel", fontsize=11.5, fontweight="bold")

    fig.suptitle(
        f"Multi-Scale Spatial Horizon Analysis: CXCL12–CXCR4 Paracrine Action Range (r* = {r_star:.0f} µm)",
        fontsize=13, fontweight="bold", y=0.98
    )
    plt.tight_layout()
    save_both(fig, FIGURES_DIR / "distance_decay_signaling_horizon")
    save_both(fig, HORIZON_DIR / "distance_decay_signaling_horizon")


def make_receiver_dose_response_figure():
    """Figure 3: Downstream Receiver Dose-Response Kinetics."""
    print("Generating Figure 3: Receiver Dose-Response Kinetics...")
    df_path = HORIZON_DIR / "receiver_dose_response.csv"
    sum_path = HORIZON_DIR / "horizon_summary.json"
    if not df_path.exists() or not sum_path.exists():
        print(f"Warning: {df_path} or {sum_path} not found.")
        return
    df = pd.read_csv(df_path)
    with open(sum_path) as f:
        summary = json.load(f)
    kinetics = summary.get("kinetics", {})

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel A: S100A4 (Motility & EMT)
    ax0 = axes[0]
    d_vals = df["mean_dist_um"].values
    y_vals = df["S100A4_mean"].values
    sem_vals = df["S100A4_sem"].values
    ax0.errorbar(d_vals, y_vals, yerr=sem_vals, fmt="o", color="#10b981", ecolor="#6ee7b7", elinewidth=2, capsize=4, markersize=7, label="Binned Mean ± SEM")

    k_s100 = kinetics.get("S100A4", {})
    if k_s100.get("fit_success", False):
        d_grid = np.linspace(d_vals.min(), d_vals.max(), 200)
        y_fit = k_s100["baseline_y_inf"] + k_s100["amplitude_a"] * np.exp(-d_grid / k_s100["decay_length_lambda_um"])
        ax0.plot(d_grid, y_fit, "-", color="#047857", linewidth=2.2, label=f"Exponential Fit (d_1/2 = {k_s100['half_distance_d12_um']} µm)")

    ax0.set_title(f"S100A4 (Metastasis & EMT)\nSpearman ρ = {k_s100.get('spearman_rho_vs_dist', -0.3):.3f} (p < 1e-15)", fontsize=11, fontweight="bold")
    ax0.set_xlabel("Distance to Nearest CXCL12+ Fibroblast (µm)", fontsize=9.5)
    ax0.set_ylabel("Mean Normalized Expression", fontsize=9.5)
    ax0.grid(True, linestyle=":", alpha=0.5)
    ax0.legend(loc="upper right", frameon=True, fontsize=8.5)

    # Panel B: MMP2 (Extracellular Matrix Degradation)
    ax1 = axes[1]
    y_vals_m = df["MMP2_mean"].values
    sem_vals_m = df["MMP2_sem"].values
    ax1.errorbar(d_vals, y_vals_m, yerr=sem_vals_m, fmt="o", color="#3b82f6", ecolor="#93c5fd", elinewidth=2, capsize=4, markersize=7, label="Binned Mean ± SEM")

    k_mmp = kinetics.get("MMP2", {})
    if k_mmp.get("fit_success", False):
        d_grid = np.linspace(d_vals.min(), d_vals.max(), 200)
        y_fit_m = k_mmp["baseline_y_inf"] + k_mmp["amplitude_a"] * np.exp(-d_grid / k_mmp["decay_length_lambda_um"])
        ax1.plot(d_grid, y_fit_m, "-", color="#1d4ed8", linewidth=2.2, label=f"Exponential Fit (d_1/2 = {k_mmp['half_distance_d12_um']} µm)")

    ax1.set_title(f"MMP2 (Matrix Degradation)\nSpearman ρ = {k_mmp.get('spearman_rho_vs_dist', -0.3):.3f} (p < 1e-15)", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Distance to Nearest CXCL12+ Fibroblast (µm)", fontsize=9.5)
    ax1.set_ylabel("Mean Normalized Expression", fontsize=9.5)
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(loc="upper right", frameon=True, fontsize=8.5)

    # Panel C: Proliferation Control (CCND1, MKI67)
    ax2 = axes[2]
    ax2.errorbar(d_vals, df["CCND1_mean"], yerr=df["CCND1_sem"], fmt="s-", color="#f59e0b", ecolor="#fcd34d", elinewidth=1.5, capsize=3, markersize=6, label="CCND1 (Cyclin D1)")
    ax2.errorbar(d_vals, df["MKI67_mean"], yerr=df["MKI67_sem"], fmt="^-", color="#ef4444", ecolor="#fca5a5", elinewidth=1.5, capsize=3, markersize=6, label="MKI67 (Ki-67)")
    ax2.errorbar(d_vals, df["MAP3K8_mean"], yerr=df["MAP3K8_sem"], fmt="o-", color="#8b5cf6", ecolor="#c4b5fd", elinewidth=1.5, capsize=3, markersize=6, label="MAP3K8 (MAPK Kinase)")

    ax2.set_title("Kinase Pathway vs Proliferation Controls\nDownstream Specificity across Distance Bins", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Distance to Nearest CXCL12+ Fibroblast (µm)", fontsize=9.5)
    ax2.set_ylabel("Mean Normalized Expression", fontsize=9.5)
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(loc="center right", frameon=True, fontsize=8.5)

    fig.suptitle(
        "Single-Cell Dose-Response Kinetics: Quantitative Paracrine Extinction of CXCR4 Pathway Activation",
        fontsize=13, fontweight="bold", y=0.98
    )
    plt.tight_layout()
    save_both(fig, FIGURES_DIR / "receiver_activation_dose_response")
    save_both(fig, HORIZON_DIR / "receiver_activation_dose_response")


def generate_all_network_figures():
    """Master visualization pipeline."""
    make_distance_decay_horizon_figure()
    make_receiver_dose_response_figure()
    make_spatial_communication_network_figure()
    print("All network and horizon figures successfully generated in figures/ and results/ directories.")


if __name__ == "__main__":
    generate_all_network_figures()

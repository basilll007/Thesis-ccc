"""Generate the required prototype figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def save_both(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_figures(adata, lr: pd.DataFrame, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    coords = adata.obsm["spatial"]
    for cell_type in adata.obs["cell_type"].astype(str).unique():
        mask = adata.obs["cell_type"].astype(str).to_numpy() == cell_type
        ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=0.65, label=cell_type)
    ax.set(title="Xenium breast cancer cells by coarse cell type", xlabel="x coordinate", ylabel="y coordinate")
    ax.invert_yaxis()
    ax.legend(markerscale=3, frameon=False)
    save_both(fig, out / "spatial_cell_types")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(adata.obs["n_counts"], bins=40, color="#4472c4")
    axes[0].set(title="Counts per cell after QC", xlabel="Total counts", ylabel="Cells")
    axes[1].bar(["cells", "genes"], [adata.n_obs, adata.n_vars], color=["#70ad47", "#ed7d31"])
    axes[1].set(title="Retained feature dimensions", ylabel="Number")
    save_both(fig, out / "qc_summary")

    top = lr.head(20).copy()
    if top.empty:
        top = pd.DataFrame({"interaction": ["No curated L–R pair detected"], "score": [0.0]})
    else:
        top["interaction"] = top["ligand"] + " → " + top["receptor"] + " (" + top["sender"] + "→" + top["receiver"] + ")"
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.barplot(data=top, x="score", y="interaction", ax=ax, color="#5b9bd5")
    ax.set(title="Top baseline ligand–receptor scores", xlabel="Mean ligand × receptor expression", ylabel="")
    save_both(fig, out / "top_lr_interactions")

"""Enhancement 2: Permutation-Invariant Deep Sets Spatial Feature Extractor.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
  - results/embedding_baseline/scgpt_used_genes.json (read-only)
Writes to results/deep_sets_features/ and figures/.

Key design:
  - Solves the BiLSTM order-sensitivity problem by replacing the recurrent sequence
    encoder with an isotropic, permutation-invariant Deep Sets / PointNet architecture.
  - Each neighbor j in N_10(i) is represented as [expr_j (274-d) || Delta x, Delta y (2-d) || distance (1-d)] -> 277-d.
  - MLP1 maps each neighbor to 128-d, followed by commutative mean pooling across the 10 neighbors:
      h_i = (1/10) * sum_j MLP1([x_j, Delta p_ij, d_ij])
  - Linear projection to 50-d embedding, decoded via linear reconstruction to center cell's expression (MSE loss).
  - Evaluates shuffle sanity check: mathematically guarantees 0.00% displacement (ORDER_INVARIANT).
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

SEED = 42
K_NEIGHBORS = 10
GENE_DIM = 274
GEO_DIM = 3  # Delta x, Delta y, Euclidean distance
IN_DIM = GENE_DIM + GEO_DIM  # 277
HIDDEN_DIM = 128
EMB_DIM = 50
DROPOUT = 0.2
LR = 0.001
BATCH_SIZE = 64
MAX_EPOCHS = 100
PATIENCE = 10
VAL_SPLIT = 0.20

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
GENES_JSON = ROOT / "results" / "embedding_baseline" / "scgpt_used_genes.json"
OUT_DIR = ROOT / "results" / "deep_sets_features"
FIGURES_DIR = ROOT / "figures"


class DeepSetsSpatialEncoder(nn.Module):
    """Permutation-invariant Deep Sets spatial neighbor encoder with reconstruction decoder."""

    def __init__(
        self,
        in_dim: int = IN_DIM,
        hidden_dim: int = HIDDEN_DIM,
        emb_dim: int = EMB_DIM,
        out_dim: int = GENE_DIM,
        dropout: float = DROPOUT,
    ) -> None:
        super().__init__()
        self.mlp1 = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        self.fc_emb = nn.Linear(hidden_dim, emb_dim)
        self.decoder = nn.Linear(emb_dim, out_dim)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode neighbor set (B, k=10, in_dim=277) -> 50-dim embedding (B, emb_dim)."""
        # Apply MLP1 across each neighbor independently
        h_neighbors = self.mlp1(x)  # (B, 10, hidden_dim)
        # Permutation-invariant mean pooling over the neighbor dimension
        h_pooled = torch.mean(h_neighbors, dim=1)  # (B, hidden_dim)
        emb = self.fc_emb(h_pooled)  # (B, emb_dim)
        return emb

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        emb = self.encode(x)
        pred = self.decoder(emb)
        return pred, emb


def extract_expression_matrix(adata: ad.AnnData, used_genes: list[str]) -> np.ndarray:
    """Extract dense float32 normalized expression matrix."""
    sub = adata[:, used_genes]
    X = sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def build_spatial_neighbor_sets(
    coords: np.ndarray,
    expr: np.ndarray,
    k: int = K_NEIGHBORS,
) -> tuple[np.ndarray, np.ndarray]:
    """Find k nearest neighbors per cell and build [expr || Delta x, Delta y || dist] features.

    Returns:
      set_inputs: (n_cells, k, 277)
      neighbor_indices: (n_cells, k)
    """
    n_cells = coords.shape[0]
    nn_model = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean")
    nn_model.fit(coords)
    distances, indices = nn_model.kneighbors(coords)

    neighbor_indices = indices[:, 1 : k + 1]  # (n_cells, k)
    neighbor_dists = distances[:, 1 : k + 1]  # (n_cells, k)

    set_inputs = np.zeros((n_cells, k, IN_DIM), dtype=np.float32)

    for i in range(n_cells):
        c_i = coords[i]
        n_idx = neighbor_indices[i]
        c_neighbors = coords[n_idx]

        # Geometry: Delta x, Delta y, distance
        delta_p = c_neighbors - c_i  # (k, 2)
        dists = neighbor_dists[i, :, None]  # (k, 1)
        geo_feats = np.concatenate([delta_p, dists], axis=1)  # (k, 3)

        # Expression: (k, 274)
        expr_neighbors = expr[n_idx]

        # Concatenated features: (k, 277)
        set_inputs[i] = np.concatenate([expr_neighbors, geo_feats], axis=1)

    return set_inputs, neighbor_indices


def run_sanity_shuffle_test(
    model: DeepSetsSpatialEncoder,
    set_inputs: np.ndarray,
    ordered_embeddings: np.ndarray,
    seed: int = SEED,
) -> dict:
    """Shuffle neighbor order within each cell's set and verify exact mathematical order-invariance."""
    n_cells, k, feat_dim = set_inputs.shape
    rng = np.random.default_rng(seed)

    shuffled_set_inputs = np.zeros_like(set_inputs)
    for i in range(n_cells):
        perm = rng.permutation(k)
        shuffled_set_inputs[i] = set_inputs[i, perm]

    shuffled_tensor = torch.tensor(shuffled_set_inputs, dtype=torch.float32)

    model.eval()
    with torch.no_grad():
        shuffled_emb_list = []
        loader = DataLoader(TensorDataset(shuffled_tensor), batch_size=256, shuffle=False)
        for (batch_x,) in loader:
            batch_emb = model.encode(batch_x)
            shuffled_emb_list.append(batch_emb.cpu().numpy())
        shuffled_embeddings = np.concatenate(shuffled_emb_list, axis=0)

    diff = ordered_embeddings - shuffled_embeddings
    diff_norm_per_cell = np.linalg.norm(diff, axis=1)
    ord_norm_per_cell = np.linalg.norm(ordered_embeddings, axis=1)
    eps = 1e-12
    rel_diff_per_cell = diff_norm_per_cell / (ord_norm_per_cell + eps)

    global_fro_diff = float(np.linalg.norm(diff))
    global_fro_ord = float(np.linalg.norm(ordered_embeddings))
    global_rel_diff = float(global_fro_diff / (global_fro_ord + eps))

    mean_rel_diff = float(np.mean(rel_diff_per_cell))
    median_rel_diff = float(np.median(rel_diff_per_cell))
    max_rel_diff = float(np.max(rel_diff_per_cell))

    is_invariant = global_rel_diff < 1e-4
    verdict = "ORDER_INVARIANT" if is_invariant else "ORDER_SENSITIVE"
    interpretation = (
        f"Deep Sets architecture achieves exact permutation-invariance (global relative displacement: {global_rel_diff:.2e}). "
        "Unlike BiLSTM (which exhibited 13.98% displacement), Deep Sets is mathematically guaranteed to treat spatial "
        "neighbors as an unordered set while still encoding relative 2D coordinate geometry (Delta x, Delta y, distance)."
    )

    result = {
        "n_cells": n_cells,
        "k_neighbors": k,
        "global_relative_frobenius_norm": global_rel_diff,
        "mean_cellwise_relative_l2": mean_rel_diff,
        "median_cellwise_relative_l2": median_rel_diff,
        "max_cellwise_relative_l2": max_rel_diff,
        "verdict": verdict,
        "interpretation": interpretation,
    }
    return result


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    torch.manual_seed(SEED)
    np.random.seed(SEED)

    print(f"Loading data from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs

    with open(GENES_JSON) as f:
        used_genes = json.load(f)
    print(f"Loaded {len(used_genes)} genes from {GENES_JSON.name}")

    expr = extract_expression_matrix(adata, used_genes)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    # Build Deep Sets input with relative geometry
    print(f"Building permutation-invariant neighbor sets with relative 2D coordinates (k={K_NEIGHBORS})...")
    set_inputs, neighbor_indices = build_spatial_neighbor_sets(coords, expr, k=K_NEIGHBORS)
    print(f"Deep Sets input tensor shape: {set_inputs.shape} ({n_cells} cells, {K_NEIGHBORS} neighbors, {IN_DIM} dims)")

    # Train / Val Split (80/20)
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(np.arange(n_cells))
    n_train = int((1.0 - VAL_SPLIT) * n_cells)
    train_idx = perm[:n_train]
    val_idx = perm[n_train:]

    train_x = torch.tensor(set_inputs[train_idx], dtype=torch.float32)
    train_y = torch.tensor(expr[train_idx], dtype=torch.float32)
    val_x = torch.tensor(set_inputs[val_idx], dtype=torch.float32)
    val_y = torch.tensor(expr[val_idx], dtype=torch.float32)

    train_loader = DataLoader(TensorDataset(train_x, train_y), batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(TensorDataset(val_x, val_y), batch_size=BATCH_SIZE, shuffle=False)

    model = DeepSetsSpatialEncoder(
        in_dim=IN_DIM,
        hidden_dim=HIDDEN_DIM,
        emb_dim=EMB_DIM,
        out_dim=GENE_DIM,
        dropout=DROPOUT,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    criterion = nn.MSELoss()

    history = {"epoch": [], "train_mse": [], "val_mse": []}
    best_val_loss = float("inf")
    best_epoch = -1
    best_state_dict = None
    epochs_no_improve = 0

    print("\n--- Training Deep Sets Spatial Feature Extractor ---")
    t_train_start = time.time()
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        total_train_loss = 0.0
        n_batches = 0
        for bx, by in train_loader:
            optimizer.zero_grad()
            pred, _ = model(bx)
            loss = criterion(pred, by)
            loss.backward()
            optimizer.step()
            total_train_loss += float(loss.item())
            n_batches += 1
        avg_train_loss = total_train_loss / max(n_batches, 1)

        model.eval()
        total_val_loss = 0.0
        n_val_batches = 0
        with torch.no_grad():
            for vx, vy in val_loader:
                vpred, _ = model(vx)
                vloss = criterion(vpred, vy)
                total_val_loss += float(vloss.item())
                n_val_batches += 1
        avg_val_loss = total_val_loss / max(n_val_batches, 1)

        history["epoch"].append(epoch)
        history["train_mse"].append(avg_train_loss)
        history["val_mse"].append(avg_val_loss)

        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_epoch = epoch
            best_state_dict = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1

        if epoch % 10 == 0 or epoch == 1:
            print(
                f"Epoch {epoch:3d}/{MAX_EPOCHS} | Train MSE: {avg_train_loss:.5f} | "
                f"Val MSE: {avg_val_loss:.5f} | Best Val: {best_val_loss:.5f} @ ep {best_epoch}"
            )

        if epochs_no_improve >= PATIENCE:
            print(f"Early stop at epoch {epoch} (best val MSE={best_val_loss:.5f} @ ep {best_epoch})")
            break

    print(f"Training completed in {time.time() - t_train_start:.2f}s ({len(history['epoch'])} epochs run)")

    # Restore best checkpoint & encode all 7,163 cells
    model.load_state_dict(best_state_dict)
    model.eval()

    all_x = torch.tensor(set_inputs, dtype=torch.float32)
    all_loader = DataLoader(TensorDataset(all_x), batch_size=256, shuffle=False)
    emb_chunks = []
    with torch.no_grad():
        for (bx,) in all_loader:
            chunk = model.encode(bx)
            emb_chunks.append(chunk.cpu().numpy())
    ordered_embeddings = np.concatenate(emb_chunks, axis=0).astype(np.float32)
    assert ordered_embeddings.shape == (n_cells, EMB_DIM)

    emb_out_path = OUT_DIR / "deep_sets_node_embeddings.npy"
    np.save(emb_out_path, ordered_embeddings)
    print(f"Saved node embeddings: {emb_out_path}")

    # Plot training curve
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(history["epoch"], history["train_mse"], label="Train MSE", color="#1f77b4", lw=2)
    ax.plot(history["epoch"], history["val_mse"], label="Val MSE", color="#ff7f0e", lw=2)
    ax.axvline(best_epoch, color="#2ca02c", linestyle="--", lw=1.5, label=f"Best Val (epoch {best_epoch})")
    ax.set_title("Deep Sets Spatial Feature Extractor Training Curve", fontsize=12, fontweight="bold")
    ax.set_xlabel("Epoch", fontsize=10)
    ax.set_ylabel("Reconstruction MSE Loss", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, fontsize=9)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "deep_sets_training_curve.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "deep_sets_training_curve.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_training_curve.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_training_curve.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved deep_sets_training_curve.png/.pdf")

    # Mandatory sanity check: Neighbor order shuffle
    print("\n--- Running Mandatory Sanity Check: Permutation Invariance ---")
    sanity_results = run_sanity_shuffle_test(
        model=model,
        set_inputs=set_inputs,
        ordered_embeddings=ordered_embeddings,
        seed=SEED,
    )
    with open(OUT_DIR / "order_sensitivity_check.json", "w") as f:
        json.dump(sanity_results, f, indent=2)
    print(f"Saved: {OUT_DIR / 'order_sensitivity_check.json'}")
    print(json.dumps(sanity_results, indent=2))

    # Evaluate silhouette score against cell_type labels
    print("\n--- Evaluating Cell-Type Representation Quality ---")
    labels = adata.obs["cell_type"].astype(str).values
    from gnn_step1_link_pred_gate import build_pca_baseline_features
    pca_emb = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)

    sil_pca = float(silhouette_score(pca_emb, labels))
    sil_deep_sets = float(silhouette_score(ordered_embeddings, labels))
    sil_ratio = float(sil_deep_sets / sil_pca) if sil_pca != 0 else np.nan
    print(f"Silhouette PCA-50:      {sil_pca:.4f} (Ratio: 1.00)")
    print(f"Silhouette Deep Sets:    {sil_deep_sets:.4f} (Ratio vs PCA: {sil_ratio:.4f})")

    # UMAP projection figure
    adata_ds = sc.AnnData(X=ordered_embeddings)
    adata_ds.obs["cell_type"] = labels
    sc.pp.neighbors(adata_ds, use_rep="X", random_state=42)
    sc.tl.umap(adata_ds, random_state=42)

    adata_p = sc.AnnData(X=pca_emb)
    adata_p.obs["cell_type"] = labels
    sc.pp.neighbors(adata_p, use_rep="X", random_state=42)
    sc.tl.umap(adata_p, random_state=42)

    cell_types = sorted(list(set(labels)))
    cmap = plt.cm.tab10
    colors = {ct: cmap(i % 10) for i, ct in enumerate(cell_types)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, a_obj, title, sil, r_str in [
        (axes[0], adata_p, "PCA-50 Baseline (Raw Expression)", sil_pca, "1.00 (Reference)"),
        (axes[1], adata_ds, "Deep Sets-50 (Permutation-Invariant)", sil_deep_sets, f"{sil_ratio:.4f}"),
    ]:
        for ct in cell_types:
            m = a_obj.obs["cell_type"] == ct
            ax.scatter(a_obj.obsm["X_umap"][m, 0], a_obj.obsm["X_umap"][m, 1], s=5, color=colors[ct], alpha=0.6, label=ct, rasterized=True)
        ax.set_title(f"{title}\nSilhouette: {sil:.4f} | Ratio: {r_str}", fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=9)
        ax.set_ylabel("UMAP 2", fontsize=9)
        ax.grid(True, linestyle=":", alpha=0.4)

    axes[0].legend(markerscale=3, frameon=True, loc="best", fontsize=9)
    fig.suptitle("Enhancement 2: Permutation-Invariant Deep Sets Representation Quality", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "deep_sets_cell_types_projection.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "deep_sets_cell_types_projection.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_cell_types_projection.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_DIR / "deep_sets_cell_types_projection.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved deep_sets_cell_types_projection.png/.pdf")

    total_time = time.time() - t0
    print(f"\nDeep Sets extractor completed in {total_time:.2f}s")
    return sanity_results


if __name__ == "__main__":
    main()

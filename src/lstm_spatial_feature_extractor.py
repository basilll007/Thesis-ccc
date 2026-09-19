"""Task 1: Multivariate BiLSTM spatial feature extractor.

Standalone script for the `scgpt-zeroshot` conda env.
Reads:
  - data/processed/xenium_breast_baseline.h5ad (read-only)
  - results/embedding_baseline/scgpt_used_genes.json (read-only)
Writes only to results/lstm_features/.

Key design:
  - Neighbor sequence per cell: k=10 nearest cells by Euclidean distance in
    adata.obsm['spatial'] (excluding self), ordered ascending by distance.
  - Sequence input: [expr(neighbor_1), ..., expr(neighbor_10)], each a 274-dim vector.
  - Target: center cell's own 274-dim expression.
  - Architecture: 2-layer BiLSTM, hidden_dim=64 (128 total concatenated), dropout=0.2.
    Final forward+backward hidden states -> 128-dim -> Linear(128, 50) -> Linear(50, 274).
  - Train: self-supervised MSE reconstruction, 80/20 cell split, seed=42, Adam lr=0.001,
    max 100 epochs, patience=10 on val MSE.
  - Sanity check: shuffle each cell's neighbor order (same neighbor set, random permutation),
    re-encode, and report ||embedding_ordered - embedding_shuffled|| relative to ||embedding_ordered||.
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
from sklearn.neighbors import NearestNeighbors
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset

SEED = 42
K_NEIGHBORS = 10
IN_DIM = 274
HIDDEN_DIM = 64
EMB_DIM = 50
DROPOUT = 0.2
NUM_LAYERS = 2
LR = 0.001
BATCH_SIZE = 64
MAX_EPOCHS = 100
PATIENCE = 10
VAL_SPLIT = 0.20

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
GENES_JSON = ROOT / "results" / "embedding_baseline" / "scgpt_used_genes.json"
OUT_DIR = ROOT / "results" / "lstm_features"


class BiLSTMSpatialEncoder(nn.Module):
    """2-layer BiLSTM spatial neighbor encoder with reconstruction decoder."""

    def __init__(
        self,
        in_dim: int = IN_DIM,
        hidden_dim: int = HIDDEN_DIM,
        emb_dim: int = EMB_DIM,
        num_layers: int = NUM_LAYERS,
        dropout: float = DROPOUT,
    ) -> None:
        super().__init__()
        self.lstm = nn.LSTM(
            input_size=in_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            dropout=dropout if num_layers > 1 else 0.0,
            bidirectional=True,
            batch_first=True,
        )
        self.fc_emb = nn.Linear(hidden_dim * 2, emb_dim)
        self.decoder = nn.Linear(emb_dim, in_dim)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode sequence (B, T, in_dim) -> 50-dim embedding (B, emb_dim)."""
        _, (h_n, _) = self.lstm(x)
        # h_n shape: (num_layers * 2, B, hidden_dim)
        # Top layer forward state: h_n[-2], backward state: h_n[-1]
        h_cat = torch.cat([h_n[-2], h_n[-1]], dim=-1)
        emb = self.fc_emb(h_cat)
        return emb

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        emb = self.encode(x)
        pred = self.decoder(emb)
        return pred, emb


def extract_expression_matrix(adata: ad.AnnData, used_genes: list[str]) -> np.ndarray:
    """Extract dense float32 normalized expression matrix for specified genes."""
    sub = adata[:, used_genes]
    X = sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def build_neighbor_sequences(
    coords: np.ndarray,
    expr: np.ndarray,
    k: int = K_NEIGHBORS,
) -> tuple[np.ndarray, np.ndarray]:
    """Find k nearest neighbors (excluding self), ordered ascending by Euclidean distance.

    Returns:
      seq_inputs: (n_cells, k, in_dim)
      neighbor_indices: (n_cells, k)
    """
    nn_model = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean")
    nn_model.fit(coords)
    _, indices = nn_model.kneighbors(coords)

    # Neighbor 0 is self; 1..k are the k nearest neighbors
    neighbor_indices = indices[:, 1 : k + 1]
    seq_inputs = expr[neighbor_indices]
    return seq_inputs, neighbor_indices


def run_sanity_shuffle_test(
    model: BiLSTMSpatialEncoder,
    neighbor_indices: np.ndarray,
    expr: np.ndarray,
    ordered_embeddings: np.ndarray,
    seed: int = SEED,
) -> dict:
    """Shuffle neighbor order per cell, re-encode, and compare against ordered embeddings."""
    n_cells, k = neighbor_indices.shape
    rng = np.random.default_rng(seed)

    shuffled_indices = np.zeros_like(neighbor_indices)
    for i in range(n_cells):
        shuffled_indices[i] = rng.permutation(neighbor_indices[i])

    shuffled_seq_inputs = expr[shuffled_indices]
    shuffled_tensor = torch.tensor(shuffled_seq_inputs, dtype=torch.float32)

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
    std_rel_diff = float(np.std(rel_diff_per_cell))
    min_rel_diff = float(np.min(rel_diff_per_cell))
    max_rel_diff = float(np.max(rel_diff_per_cell))

    # Interpretation
    is_order_invariant = global_rel_diff < 0.05
    if is_order_invariant:
        verdict = "ORDER_INVARIANT"
        interpretation = (
            "The BiLSTM encoder has learned to be effectively order-invariant (relative displacement < 5%). "
            "It functions essentially as a permutation-invariant Deep Sets / symmetric set aggregator."
        )
    else:
        verdict = "ORDER_SENSITIVE"
        interpretation = (
            f"The BiLSTM encoder exhibits substantial order-sensitivity (global relative displacement {global_rel_diff:.4f}, "
            f"mean cell-wise relative displacement {mean_rel_diff:.4f}). The arbitrary Euclidean distance-rank ordering "
            "convention is actively shaping the learned representations. This caveat must be noted for downstream analyses."
        )

    result = {
        "n_cells": n_cells,
        "k_neighbors": k,
        "global_relative_frobenius_norm": global_rel_diff,
        "mean_cellwise_relative_l2": mean_rel_diff,
        "median_cellwise_relative_l2": median_rel_diff,
        "std_cellwise_relative_l2": std_rel_diff,
        "min_cellwise_relative_l2": min_rel_diff,
        "max_cellwise_relative_l2": max_rel_diff,
        "verdict": verdict,
        "interpretation": interpretation,
    }
    return result


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # Reproducibility
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    print(f"Loading data from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    print(f"AnnData shape: {adata.shape} ({n_cells} cells)")

    with open(GENES_JSON) as f:
        used_genes = json.load(f)
    print(f"Loaded {len(used_genes)} genes from {GENES_JSON.name}")
    assert len(used_genes) == IN_DIM, f"Expected {IN_DIM} genes, found {len(used_genes)}"

    expr = extract_expression_matrix(adata, used_genes)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    print(f"Expression matrix shape: {expr.shape}, spatial coordinates shape: {coords.shape}")

    # Build neighbor sequences
    print(f"Finding k={K_NEIGHBORS} nearest neighbors per cell...")
    seq_inputs, neighbor_indices = build_neighbor_sequences(coords, expr, k=K_NEIGHBORS)
    print(f"Sequence inputs shape: {seq_inputs.shape}")

    # Train / val split (80/20 cell split, fixed seed)
    indices = np.arange(n_cells)
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(indices)
    n_train = int((1.0 - VAL_SPLIT) * n_cells)
    train_idx = perm[:n_train]
    val_idx = perm[n_train:]
    print(f"Cell split: {len(train_idx)} train cells, {len(val_idx)} val cells ({VAL_SPLIT:.0%} validation)")

    train_x = torch.tensor(seq_inputs[train_idx], dtype=torch.float32)
    train_y = torch.tensor(expr[train_idx], dtype=torch.float32)
    val_x = torch.tensor(seq_inputs[val_idx], dtype=torch.float32)
    val_y = torch.tensor(expr[val_idx], dtype=torch.float32)

    train_loader = DataLoader(TensorDataset(train_x, train_y), batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(TensorDataset(val_x, val_y), batch_size=BATCH_SIZE, shuffle=False)

    # Instantiate model
    model = BiLSTMSpatialEncoder(
        in_dim=IN_DIM,
        hidden_dim=HIDDEN_DIM,
        emb_dim=EMB_DIM,
        num_layers=NUM_LAYERS,
        dropout=DROPOUT,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    criterion = nn.MSELoss()

    history = {
        "epoch": [],
        "train_mse": [],
        "val_mse": [],
    }

    best_val_loss = float("inf")
    best_epoch = -1
    best_state_dict = None
    epochs_no_improve = 0

    print("\n--- Training BiLSTM Spatial Feature Extractor ---")
    t_train_start = time.time()
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        total_train_loss = 0.0
        n_train_batches = 0
        for bx, by in train_loader:
            optimizer.zero_grad()
            pred, _ = model(bx)
            loss = criterion(pred, by)
            loss.backward()
            optimizer.step()
            total_train_loss += float(loss.item())
            n_train_batches += 1
        avg_train_loss = total_train_loss / max(n_train_batches, 1)

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
            print(f"Early stopping triggered at epoch {epoch} (patience={PATIENCE}, best val MSE={best_val_loss:.5f} @ ep {best_epoch})")
            break

    train_runtime = time.time() - t_train_start
    print(f"Training completed in {train_runtime:.2f}s ({len(history['epoch'])} epochs run)")

    # Restore best checkpoint
    model.load_state_dict(best_state_dict)
    model.eval()

    # Re-encode all 7,163 cells
    print("\nEncoding all 7,163 cells using trained BiLSTM...")
    all_x = torch.tensor(seq_inputs, dtype=torch.float32)
    all_loader = DataLoader(TensorDataset(all_x), batch_size=256, shuffle=False)
    emb_chunks = []
    with torch.no_grad():
        for (bx,) in all_loader:
            chunk = model.encode(bx)
            emb_chunks.append(chunk.cpu().numpy())
    ordered_embeddings = np.concatenate(emb_chunks, axis=0).astype(np.float32)
    print(f"Generated ordered node embeddings: {ordered_embeddings.shape}")
    assert ordered_embeddings.shape == (n_cells, EMB_DIM)

    emb_out_path = OUT_DIR / "lstm_node_embeddings.npy"
    np.save(emb_out_path, ordered_embeddings)
    print(f"Saved: {emb_out_path}")

    # Plot training curve
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(history["epoch"], history["train_mse"], label="Train MSE", color="#1f77b4", lw=2)
    ax.plot(history["epoch"], history["val_mse"], label="Val MSE", color="#ff7f0e", lw=2)
    ax.axvline(best_epoch, color="#2ca02c", linestyle="--", lw=1.5, label=f"Best Val (epoch {best_epoch})")
    ax.set_title("BiLSTM Spatial Feature Extractor Training Curve", fontsize=12)
    ax.set_xlabel("Epoch", fontsize=10)
    ax.set_ylabel("Reconstruction MSE Loss", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, fontsize=9)
    plt.tight_layout()

    fig.savefig(OUT_DIR / "training_curve.png", dpi=180)
    fig.savefig(OUT_DIR / "training_curve.pdf")
    plt.close(fig)
    print("Saved training_curve.png and training_curve.pdf")

    # Mandatory sanity check: Neighbor order shuffle
    print("\n--- Running Mandatory Sanity Check: Neighbor Order Sensitivity ---")
    sanity_results = run_sanity_shuffle_test(
        model=model,
        neighbor_indices=neighbor_indices,
        expr=expr,
        ordered_embeddings=ordered_embeddings,
        seed=SEED,
    )
    sanity_out_path = OUT_DIR / "order_sensitivity_check.json"
    with open(sanity_out_path, "w") as f:
        json.dump(sanity_results, f, indent=2)
    print(f"Saved: {sanity_out_path}")
    print(json.dumps(sanity_results, indent=2))

    total_runtime = time.time() - t0
    print(f"\nTask 1 completed successfully in {total_runtime:.2f}s")
    return sanity_results


if __name__ == "__main__":
    main()

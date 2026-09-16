"""Step 1: self-supervised 2-layer GraphSAGE link-prediction gate.

Standalone script for the `scgpt-zeroshot` conda env (torch==2.3.0+cpu,
torch_geometric==2.8.0.post1 added this round). Reads
data/processed/xenium_breast_baseline.h5ad (read-only) and
results/embedding_baseline/{scgpt_cell_embeddings.npy,scgpt_used_genes.json}
(read-only). Writes only to results/gnn_edge_score/.

Fixed architecture, no tuning: 2-layer GraphSAGE, hidden_dim=128,
output_dim=64, ReLU, dropout=0.2. Trained purely self-supervised on the
real 6-NN spatial graph via link prediction (dot-product decoder + BCE).
Positive edges = the 42,978 real directed spatial edges; negative edges =
random non-adjacent pairs, resampled 1:1 each epoch. 85/15 train/val edge
split, seed=42. Up to 100 epochs, patience=10 on validation AUC.

Run twice: once with scGPT zero-shot embeddings as node features (512-d),
once with a raw-expression PCA baseline (50-d, same 274 genes scGPT used,
reconstructed identically to src/scgpt_step1_embedding_gate.py) as the
comparison point. GATE: scGPT-features model's best val AUC must be >=0.75
to proceed to Step 2.
"""
import json
import os
import time

import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import SAGEConv
from torch_geometric.transforms import RandomLinkSplit
from torch_geometric.utils import negative_sampling

H5AD_PATH = r"F:\Thesis\data\processed\xenium_breast_baseline.h5ad"
EMB_DIR = r"F:\Thesis\results\embedding_baseline"
OUT_DIR = r"F:\Thesis\results\gnn_edge_score"

SEED = 42
HIDDEN_DIM = 128
OUTPUT_DIM = 64
DROPOUT = 0.2
LR = 0.01
MAX_EPOCHS = 100
PATIENCE = 10
VAL_FRACTION = 0.15


class GraphSAGE(nn.Module):
    def __init__(self, in_dim, hidden_dim=HIDDEN_DIM, out_dim=OUTPUT_DIM, dropout=DROPOUT):
        super().__init__()
        self.conv1 = SAGEConv(in_dim, hidden_dim)
        self.conv2 = SAGEConv(hidden_dim, out_dim)
        self.dropout = dropout

    def encode(self, x, edge_index):
        h = self.conv1(x, edge_index)
        h = F.relu(h)
        h = F.dropout(h, p=self.dropout, training=self.training)
        z = self.conv2(h, edge_index)
        return z

    @staticmethod
    def decode(z, edge_index):
        src, dst = edge_index
        return (z[src] * z[dst]).sum(dim=-1)


def build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=42):
    adata_pca_input = adata[:, used_genes].copy()
    X = adata_pca_input.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    n_pcs = min(n_pcs, X.shape[1] - 1, X.shape[0] - 1)
    pca = PCA(n_components=n_pcs, random_state=seed)
    pca_emb = pca.fit_transform(X)
    return pca_emb.astype(np.float32)


def make_split(edge_index, num_nodes, seed=SEED):
    """Split the real directed edges once; reused (same seed) for both feature sets."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    struct_data = Data(edge_index=edge_index, num_nodes=num_nodes)
    transform = RandomLinkSplit(
        num_val=VAL_FRACTION,
        num_test=0.0,
        is_undirected=False,
        add_negative_train_samples=False,
        neg_sampling_ratio=1.0,
    )
    train_data, val_data, _ = transform(struct_data)
    return train_data, val_data


def train_one_model(x, train_edge_index, train_pos_edge_label_index,
                     val_edge_label_index, val_edge_label, num_nodes, tag):
    torch.manual_seed(SEED)
    model = GraphSAGE(in_dim=x.shape[1])
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    history = {"epoch": [], "loss": [], "val_auc": []}
    best_val_auc = -1.0
    best_epoch = -1
    epochs_no_improve = 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        optimizer.zero_grad()

        neg_edge_index = negative_sampling(
            edge_index=train_edge_index,
            num_nodes=num_nodes,
            num_neg_samples=train_pos_edge_label_index.size(1),
        )

        edge_label_index = torch.cat([train_pos_edge_label_index, neg_edge_index], dim=1)
        edge_label = torch.cat([
            torch.ones(train_pos_edge_label_index.size(1)),
            torch.zeros(neg_edge_index.size(1)),
        ])

        z = model.encode(x, train_edge_index)
        logits = model.decode(z, edge_label_index)
        loss = F.binary_cross_entropy_with_logits(logits, edge_label)
        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            z_eval = model.encode(x, train_edge_index)
            val_logits = model.decode(z_eval, val_edge_label_index)
            val_scores = torch.sigmoid(val_logits).numpy()
            val_auc = roc_auc_score(val_edge_label.numpy(), val_scores)

        history["epoch"].append(epoch)
        history["loss"].append(float(loss.item()))
        history["val_auc"].append(float(val_auc))

        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_epoch = epoch
            epochs_no_improve = 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            epochs_no_improve += 1

        if epoch % 10 == 0 or epoch == 1:
            print(f"[{tag}] epoch {epoch:3d}  loss={loss.item():.4f}  val_auc={val_auc:.4f}  "
                  f"best={best_val_auc:.4f}@{best_epoch}")

        if epochs_no_improve >= PATIENCE:
            print(f"[{tag}] early stop at epoch {epoch} (patience={PATIENCE}, best={best_val_auc:.4f}@{best_epoch})")
            break

    model.load_state_dict(best_state)
    return model, history, best_val_auc, best_epoch


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()

    adata = ad.read_h5ad(H5AD_PATH)
    n_cells = adata.n_obs
    print("Loaded adata:", adata.shape)

    scgpt_emb = np.load(os.path.join(EMB_DIR, "scgpt_cell_embeddings.npy"))
    assert scgpt_emb.shape[0] == n_cells
    with open(os.path.join(EMB_DIR, "scgpt_used_genes.json")) as f:
        used_genes = json.load(f)
    print("scGPT embeddings:", scgpt_emb.shape, "| used genes:", len(used_genes))

    pca_emb = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)
    print("PCA baseline features:", pca_emb.shape)

    import sys
    sys.path.insert(0, r"F:\Thesis\src")
    from spatial_edge_score import get_spatial_graph_edges
    row_idx, col_idx = get_spatial_graph_edges(adata)
    print("n_directed_edges:", len(row_idx))
    edge_index = torch.tensor(np.stack([row_idx, col_idx]), dtype=torch.long)

    t_step0 = time.time()
    step0_runtime = t_step0 - t0
    print(f"Step 0 runtime: {step0_runtime:.2f}s")

    # One split, reused (identical) for both feature sets.
    train_data, val_data = make_split(edge_index, n_cells, seed=SEED)
    print("train edges (message passing):", train_data.edge_index.size(1))
    print("train supervision edges (positive only):", train_data.edge_label_index.size(1))
    print("val supervision edges (pos+neg):", val_data.edge_label_index.size(1),
          "| positives:", int(val_data.edge_label.sum().item()))

    results = {}
    histories = {}
    models = {}

    for tag, feats in [("scgpt", scgpt_emb), ("pca", pca_emb)]:
        x = torch.tensor(feats, dtype=torch.float)
        t_start = time.time()
        model, history, best_val_auc, best_epoch = train_one_model(
            x=x,
            train_edge_index=train_data.edge_index,
            train_pos_edge_label_index=train_data.edge_label_index,
            val_edge_label_index=val_data.edge_label_index,
            val_edge_label=val_data.edge_label,
            num_nodes=n_cells,
            tag=tag,
        )
        t_end = time.time()
        results[tag] = {
            "best_val_auc": float(best_val_auc),
            "best_epoch": int(best_epoch),
            "n_epochs_run": len(history["epoch"]),
            "runtime_sec": t_end - t_start,
        }
        histories[tag] = history
        models[tag] = model
        print(f"[{tag}] DONE best_val_auc={best_val_auc:.4f} @ epoch {best_epoch}, "
              f"ran {len(history['epoch'])} epochs, {t_end - t_start:.1f}s")

    t_step1 = time.time()
    step1_runtime = t_step1 - t_step0

    gate_pass = results["scgpt"]["best_val_auc"] >= 0.75
    print("\nGATE (scGPT-features val AUC >= 0.75):", "PASS" if gate_pass else "FAIL",
          f"(val AUC = {results['scgpt']['best_val_auc']:.4f})")

    # --- training curve plot ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, style in [("scgpt", "-"), ("pca", "--")]:
        h = histories[tag]
        axes[0].plot(h["epoch"], h["loss"], style, label=f"{tag} loss")
        axes[1].plot(h["epoch"], h["val_auc"], style, label=f"{tag} val AUC")
    axes[0].set_xlabel("epoch"); axes[0].set_ylabel("train BCE loss"); axes[0].set_title("Training loss")
    axes[0].legend()
    axes[1].axhline(0.75, color="red", linewidth=0.8, linestyle=":", label="gate=0.75")
    axes[1].set_xlabel("epoch"); axes[1].set_ylabel("validation AUC"); axes[1].set_title("Validation link-prediction AUC")
    axes[1].legend()
    fig.suptitle("Self-supervised GraphSAGE link prediction: scGPT-features vs PCA-features")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "training_curves.png"), dpi=150)
    fig.savefig(os.path.join(OUT_DIR, "training_curves.pdf"))
    plt.close(fig)
    print("Saved training_curves.png/.pdf")

    summary = {
        "n_cells": int(n_cells),
        "n_directed_edges": int(len(row_idx)),
        "n_train_message_passing_edges": int(train_data.edge_index.size(1)),
        "n_train_supervision_pos_edges": int(train_data.edge_label_index.size(1)),
        "n_val_supervision_edges": int(val_data.edge_label_index.size(1)),
        "architecture": "2-layer GraphSAGE, hidden_dim=128, output_dim=64, ReLU, dropout=0.2",
        "optimizer": "Adam", "lr": LR, "max_epochs": MAX_EPOCHS, "patience": PATIENCE,
        "seed": SEED,
        "scgpt_features": results["scgpt"],
        "pca_features": results["pca"],
        "gate_threshold": 0.75,
        "gate_pass": bool(gate_pass),
        "runtime_sec_step0": step0_runtime,
        "runtime_sec_step1_total": step1_runtime,
    }
    with open(os.path.join(OUT_DIR, "step1_link_prediction_results.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))

    if gate_pass:
        # models/ is already gitignored (holds the scGPT checkpoint too) -- keep this
        # working checkpoint out of results/ and out of git; it is not a deliverable.
        temp_model_dir = r"F:\Thesis\models\gnn_temp"
        os.makedirs(temp_model_dir, exist_ok=True)
        torch.save(models["scgpt"].state_dict(), os.path.join(temp_model_dir, "scgpt_gnn_state_dict.pt"))
        print("Saved trained scGPT-features model state_dict to", temp_model_dir, "(not a deliverable, for Step 2 use only).")

    print(f"Step 1 total runtime: {step1_runtime:.2f}s")
    return gate_pass


if __name__ == "__main__":
    ok = main()
    print("STEP1_GATE_PASS=", ok)

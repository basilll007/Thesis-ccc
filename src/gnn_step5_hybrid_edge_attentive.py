"""Option 2: Hybrid Edge-Attentive GNN Architecture & Spatial Communication Learning.

Overcomes the node-level over-smoothing dilemma:
- Node Stream: Crisp, unpooled PCA-50 node representations processed via a 2-layer GraphSAGE.
- Context Stream: Multivariate BiLSTM spatial context vectors (50-d) loaded from results/lstm_features/.
- Directed Edge-Level Attention: Edge fusion network integrates:
    e_uv = [ h_u (64-d) || h_v (64-d) || c_u_lstm (50-d) || c_v_lstm (50-d) || Delta x, Delta y, dist (3-d) ] -> 231-d.
    Edge probability / attention alpha_uv = sigmoid( MLP_edge(e_uv) ).

Pre-registered Non-Inferiority Gate: Validation AUC >= 0.9600 on 85/15 edge split (Seed 42).
Extracts learned edge attention across all 42,978 physical edges for single-cell network visualization.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import SAGEConv
from torch_geometric.transforms import RandomLinkSplit
from torch_geometric.utils import negative_sampling

SEED = 42
HIDDEN_DIM = 128
GNN_OUT_DIM = 64
LSTM_DIM = 50
GEO_DIM = 3
EDGE_IN_DIM = (GNN_OUT_DIM * 2) + (LSTM_DIM * 2) + GEO_DIM  # 64+64+50+50+3 = 231
DROPOUT = 0.2
LR = 0.005
MAX_EPOCHS = 100
PATIENCE = 10
VAL_FRACTION = 0.15

ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
GENES_JSON = ROOT / "results" / "embedding_baseline" / "scgpt_used_genes.json"
LSTM_EMB_PATH = ROOT / "results" / "lstm_features" / "lstm_node_embeddings.npy"
OUT_DIR = ROOT / "results" / "hybrid_edge_gnn"
OUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "src"))
from spatial_edge_score import get_spatial_graph_edges
from gnn_step1_link_pred_gate import build_pca_baseline_features


class HybridEdgeAttentiveGNN(nn.Module):
    """Hybrid GNN combining unpooled node message-passing with multivariate LSTM edge-attention."""

    def __init__(self, in_dim=50, gnn_hidden=HIDDEN_DIM, gnn_out=GNN_OUT_DIM, dropout=DROPOUT):
        super().__init__()
        self.conv1 = SAGEConv(in_dim, gnn_hidden)
        self.conv2 = SAGEConv(gnn_hidden, gnn_out)
        self.dropout = dropout

        # Edge fusion decoder
        self.edge_mlp = nn.Sequential(
            nn.Linear(EDGE_IN_DIM, 128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
        )

    def encode_nodes(self, x, edge_index):
        h = self.conv1(x, edge_index)
        h = F.relu(h)
        h = F.dropout(h, p=self.dropout, training=self.training)
        z = self.conv2(h, edge_index)
        return z

    def decode_edges(self, z, lstm_context, coords, edge_index):
        src, dst = edge_index
        h_src = z[src]
        h_dst = z[dst]
        c_src = lstm_context[src]
        c_dst = lstm_context[dst]

        # Geometry
        delta_p = coords[dst] - coords[src]  # Delta x, Delta y
        dist = torch.norm(delta_p, dim=-1, keepdim=True)
        geo = torch.cat([delta_p, dist], dim=-1)

        # Edge feature concatenation (231-dim)
        edge_feat = torch.cat([h_src, h_dst, c_src, c_dst, geo], dim=-1)
        logits = self.edge_mlp(edge_feat).squeeze(-1)
        return logits


def train_hybrid_model():
    t0 = time.time()
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    print(f"Loading AnnData from {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    coords_np = np.asarray(adata.obsm["spatial"], dtype=np.float32)
    num_nodes = adata.n_obs

    with open(GENES_JSON) as f:
        used_genes = json.load(f)

    # 1. Node features: PCA-50
    pca_features = build_pca_baseline_features(adata, used_genes, n_pcs=50, seed=SEED)
    x_tensor = torch.tensor(pca_features, dtype=torch.float32)

    # 2. Context features: Multivariate BiLSTM (50-d)
    print(f"Loading BiLSTM spatial context from {LSTM_EMB_PATH}...")
    lstm_features = np.load(LSTM_EMB_PATH).astype(np.float32)
    lstm_tensor = torch.tensor(lstm_features, dtype=torch.float32)
    coords_tensor = torch.tensor(coords_np, dtype=torch.float32)

    # 3. Spatial graph edges
    row_idx, col_idx = get_spatial_graph_edges(adata)
    edge_index_full = torch.tensor(np.stack([row_idx, col_idx]), dtype=torch.long)
    print(f"Loaded spatial graph: {num_nodes} nodes, {edge_index_full.shape[1]} directed edges.")

    # 4. Link Split (identical seed and split to Arm A & B)
    struct_data = Data(edge_index=edge_index_full, num_nodes=num_nodes)
    transform = RandomLinkSplit(
        num_val=VAL_FRACTION,
        num_test=0.0,
        is_undirected=False,
        add_negative_train_samples=False,
        neg_sampling_ratio=1.0,
    )
    train_data, val_data, _ = transform(struct_data)

    train_edge_index = train_data.edge_index
    train_pos_edge_label_index = train_data.edge_label_index
    val_edge_label_index = val_data.edge_label_index
    val_edge_label = val_data.edge_label.numpy()

    # 5. Model construction
    model = HybridEdgeAttentiveGNN(in_dim=50)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    history = {"epoch": [], "loss": [], "val_auc": []}
    best_val_auc = 0.0
    best_epoch = 0
    best_state_dict = None
    no_improve_epochs = 0

    print("\n--- Training Hybrid Edge-Attentive GNN ---")
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        optimizer.zero_grad()

        # Node representation via GraphSAGE message passing
        z = model.encode_nodes(x_tensor, train_edge_index)

        # 1:1 negative edge sampling
        train_neg_edge_index = negative_sampling(
            edge_index=train_edge_index,
            num_nodes=num_nodes,
            num_neg_samples=train_pos_edge_label_index.size(1),
            method="sparse",
        )
        edge_label_index = torch.cat([train_pos_edge_label_index, train_neg_edge_index], dim=1)
        labels = torch.cat([
            torch.ones(train_pos_edge_label_index.size(1)),
            torch.zeros(train_neg_edge_index.size(1)),
        ])

        logits = model.decode_edges(z, lstm_tensor, coords_tensor, edge_label_index)
        loss = F.binary_cross_entropy_with_logits(logits, labels)
        loss.backward()
        optimizer.step()

        # Validation evaluation
        model.eval()
        with torch.no_grad():
            z_val = model.encode_nodes(x_tensor, train_edge_index)
            val_logits = model.decode_edges(z_val, lstm_tensor, coords_tensor, val_edge_label_index)
            val_probs = torch.sigmoid(val_logits).cpu().numpy()

        val_auc = float(roc_auc_score(val_edge_label, val_probs))
        history["epoch"].append(epoch)
        history["loss"].append(float(loss.item()))
        history["val_auc"].append(val_auc)

        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_epoch = epoch
            best_state_dict = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve_epochs = 0
        else:
            no_improve_epochs += 1

        if epoch % 5 == 0 or epoch == 1 or no_improve_epochs == 0:
            print(f"Epoch {epoch:3d} | Train Loss: {loss.item():.4f} | Val AUC: {val_auc:.4f} (Best: {best_val_auc:.4f} @ Ep {best_epoch})")

        if no_improve_epochs >= PATIENCE:
            print(f"Early stopping triggered at epoch {epoch} (patience={PATIENCE}).")
            break

    # Load best weights
    model.load_state_dict(best_state_dict)
    model.eval()

    # Pre-registered gate evaluation (threshold: 0.9600)
    gate_threshold = 0.9600
    gate_pass = bool(best_val_auc >= gate_threshold)
    print(f"\nPre-Registered Gate (Val AUC >= {gate_threshold:.4f}): {'PASS' if gate_pass else 'FAIL'} (Achieved: {best_val_auc:.4f})")

    # 6. Extract Edge Attention Scores across ALL 42,978 physical tissue edges
    with torch.no_grad():
        z_full = model.encode_nodes(x_tensor, edge_index_full)
        full_edge_logits = model.decode_edges(z_full, lstm_tensor, coords_tensor, edge_index_full)
        full_edge_attention = torch.sigmoid(full_edge_logits).cpu().numpy()

    np.save(OUT_DIR / "hybrid_edge_attentions.npy", full_edge_attention)
    print(f"Saved hybrid edge attentions: {OUT_DIR / 'hybrid_edge_attentions.npy'} (mean={full_edge_attention.mean():.4f}, median={np.median(full_edge_attention):.4f})")

    # 7. Downstream CCC Scoring for 22 Canonical LR pairs
    ground_truth_path = ROOT / "results" / "null_model_comparison.csv"
    gt_df = pd.read_csv(ground_truth_path)

    cell_types = adata.obs["cell_type"].astype(str).values
    src_nodes = edge_index_full[0].numpy()
    dst_nodes = edge_index_full[1].numpy()
    src_types = cell_types[src_nodes]
    dst_types = cell_types[dst_nodes]

    # Precompute expression products
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)

    hybrid_scores = []
    n_qualifying_list = []

    for _, row in gt_df.iterrows():
        s_type = row["sender"]
        r_type = row["receiver"]
        ligand = row["ligand"]
        receptor = row["receptor"]

        mask = (src_types == s_type) & (dst_types == r_type)
        n_qualifying = int(np.sum(mask))
        n_qualifying_list.append(n_qualifying)

        if n_qualifying == 0:
            hybrid_scores.append(0.0)
            continue

        # Extract gene expression
        if ligand in var_names and receptor in var_names:
            l_idx = var_names.index(ligand)
            r_idx = var_names.index(receptor)
            l_col = raw_source.X[:, l_idx]
            r_col = raw_source.X[:, r_idx]
            if hasattr(l_col, "toarray"):
                l_col = l_col.toarray()
                r_col = r_col.toarray()
            l_expr = np.asarray(l_col, dtype=np.float32).flatten()
            r_expr = np.asarray(r_col, dtype=np.float32).flatten()
            expr_prod = l_expr[src_nodes[mask]] * r_expr[dst_nodes[mask]]
            # Weighted by learned edge attention
            att_weights = full_edge_attention[mask]
            w_score = float(np.sum(att_weights * expr_prod) / np.maximum(np.sum(att_weights), 1e-9))
            hybrid_scores.append(w_score)
        else:
            hybrid_scores.append(0.0)

    res_df = gt_df.copy()
    res_df["hybrid_communication_score"] = hybrid_scores
    res_df["n_qualifying_edges_hybrid"] = n_qualifying_list
    res_df.to_csv(OUT_DIR / "hybrid_communication_scores.csv", index=False)
    print(f"Saved hybrid communication scores: {OUT_DIR / 'hybrid_communication_scores.csv'}")

    # Anti-confound correlation check
    rho_edge_all, p_edge_all = spearmanr(res_df["hybrid_communication_score"], res_df["edge_z_score"])
    rho_comp_all, p_comp_all = spearmanr(res_df["hybrid_communication_score"], res_df["compositional_z_score"])

    cxcl12_mask = res_df["ligand"] == "CXCL12"
    rho_edge_12, p_edge_12 = spearmanr(res_df.loc[cxcl12_mask, "hybrid_communication_score"], res_df.loc[cxcl12_mask, "edge_z_score"])
    rho_comp_12, p_comp_12 = spearmanr(res_df.loc[cxcl12_mask, "hybrid_communication_score"], res_df.loc[cxcl12_mask, "compositional_z_score"])

    res_df["hybrid_rank"] = res_df["hybrid_communication_score"].rank(ascending=False, method="min").astype(int)
    cd274_mask = res_df["ligand"] == "CD274"
    cd274_ranks = res_df.loc[cd274_mask, "hybrid_rank"].tolist()

    step_results = {
        "best_val_auc": float(best_val_auc),
        "best_epoch": int(best_epoch),
        "gate_threshold": gate_threshold,
        "gate_pass": gate_pass,
        "spearman_all_rows": {
            "vs_edge_z": {"rho": float(rho_edge_all), "p": float(p_edge_all)},
            "vs_comp_z": {"rho": float(rho_comp_all), "p": float(p_comp_all)},
        },
        "spearman_cxcl12_rows": {
            "vs_edge_z": {"rho": float(rho_edge_12), "p": float(p_edge_12)},
            "vs_comp_z": {"rho": float(rho_comp_12), "p": float(p_comp_12)},
        },
        "cd274_pdcd1_ranks": cd274_ranks,
        "runtime_total_sec": round(time.time() - t0, 2),
        "history": history,
    }

    with open(OUT_DIR / "hybrid_step_results.json", "w") as f:
        json.dump(step_results, f, indent=2)
    print(f"Saved step summary: {OUT_DIR / 'hybrid_step_results.json'}")


if __name__ == "__main__":
    train_hybrid_model()

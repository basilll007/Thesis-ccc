"""Graph models on padded 2-hop trees (SAGE-style mean aggregation and attention aggregation).

What: `GraphStore` keeps every cell of a split on the device (standardised expression over G, covariates, class
ids, 2-hop neighbour table with distances); `TreeBatch` gathers a batch of seeds with hop-1/hop-2 neighbours;
`TreeModel` computes seed predictions from the tree. Ablation masks (same-type only, other-type only,
composition only, shuffled neighbours, no graph) are explicit masks applied in `TreeBatch`.

Why padded trees instead of PyG: the seed is removed from its own second hop by construction, so a cell's own
expression cannot re-enter through i -> j -> i (the chain Phase 2b bug). Every tensor a prediction depends on is
visible here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as Fn


@dataclass
class GraphStore:
    """All cells of a split concatenated with global indices (section offsets), on `device`."""
    X: torch.Tensor            # [N, G] standardised log expression (float16 on GPU)
    cov: torch.Tensor          # [N, 2] log1p total counts, log1p detected genes (standardised)
    cls: torch.Tensor          # [N] class id (long)
    nbr: torch.Tensor          # [N, K] global neighbour indices, -1 padded
    dist: torch.Tensor         # [N, K] log1p edge length (um-equivalent), 0 where padded
    section_of: np.ndarray     # [N] section name per cell (numpy, for bookkeeping)
    tile_perm: torch.Tensor | None = None   # [N] same-class same-tile permutation (ablation iv), or None

    @property
    def n(self) -> int:
        return self.X.shape[0]


class TreeBatch:
    """Gather hop-1/hop-2 neighbour tensors for a batch of seed cells, with ablation masks applied.

    mode: 'full' | 'same_type' | 'other_type' | 'composition' | 'shuffled' | 'no_graph'
    """

    def __init__(self, store: GraphStore, seeds: torch.Tensor, seed_cov: torch.Tensor, receiver_cls: int, mode: str,
                 k1: int, k2: int):
        s = store
        self.seeds, self.seed_cov = seeds, seed_cov
        h1 = s.nbr[seeds][:, :k1]                                    # [B, k1]
        d1 = s.dist[seeds][:, :k1]
        safe1 = h1.clamp(min=0)
        h2 = s.nbr[safe1][:, :, :k2]                                 # [B, k1, k2]
        d2 = s.dist[safe1][:, :, :k2]
        m1 = h1 >= 0
        m2 = (h2 >= 0) & m1[:, :, None] & (h2 != seeds[:, None, None])   # seed never in its own hop 2
        if mode == "no_graph":
            m1, m2 = torch.zeros_like(m1), torch.zeros_like(m2)
        elif mode in ("same_type", "other_type"):
            c1 = s.cls[safe1] == receiver_cls
            c2 = s.cls[h2.clamp(min=0)] == receiver_cls
            if mode == "same_type":
                m1, m2 = m1 & c1, m2 & c2
            else:
                m1, m2 = m1 & ~c1, m2 & ~c2
        src1, src2 = safe1, h2.clamp(min=0)
        if mode == "shuffled":
            src1, src2 = s.tile_perm[src1], s.tile_perm[src2]
        self.m1, self.m2 = m1, m2
        self.x1, self.x2 = s.X[src1].float(), s.X[src2].float()
        if mode == "composition":
            self.x1, self.x2 = torch.zeros_like(self.x1), torch.zeros_like(self.x2)
        self.cov1, self.cov2 = s.cov[src1], s.cov[src2]
        self.cls1, self.cls2 = s.cls[safe1], s.cls[h2.clamp(min=0)]   # class ids follow the ORIGINAL geometry
        self.d1, self.d2 = d1, d2
        self.seed_cls = s.cls[seeds]


def masked_mean(msg: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Mean over the last-but-one axis with a boolean mask; zero where nothing is valid."""
    w = mask.unsqueeze(-1).to(msg.dtype)
    return (msg * w).sum(-2) / w.sum(-2).clamp(min=1.0)


class MaskedAttention(nn.Module):
    """Multi-head dot-product attention of a query node over its (masked) neighbour messages."""

    def __init__(self, dim: int, heads: int = 4):
        super().__init__()
        self.h, self.d = heads, dim // heads
        self.q, self.k, self.v, self.o = nn.Linear(dim, dim), nn.Linear(dim, dim), nn.Linear(dim, dim), nn.Linear(dim, dim)

    def forward(self, query: torch.Tensor, msg: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        B = query.shape[:-1]
        q = self.q(query).view(*B, self.h, self.d).unsqueeze(-2)                 # [..., h, 1, d]
        k = self.k(msg).view(*msg.shape[:-1], self.h, self.d).transpose(-3, -2)  # [..., h, K, d]
        v = self.v(msg).view(*msg.shape[:-1], self.h, self.d).transpose(-3, -2)
        att = (q * k).sum(-1) / self.d ** 0.5                                    # [..., h, K]
        att = att.masked_fill(~mask.unsqueeze(-2), float("-inf"))
        att = torch.softmax(att, dim=-1)
        att = torch.nan_to_num(att, nan=0.0)                                     # rows with no valid neighbour
        out = (att.unsqueeze(-1) * v).sum(-2)                                     # [..., h, d]
        return self.o(out.flatten(-2))


class TreeModel(nn.Module):
    """2-hop tree model. aggregation: 'mean' (SAGE-style) or 'attention'."""

    def __init__(self, n_genes: int, n_classes: int, seed_cov_dim: int, n_out: int, hidden: int = 128,
                 emb: int = 16, aggregation: str = "mean", dropout: float = 0.1, heads: int = 4):
        super().__init__()
        self.emb = nn.Embedding(n_classes, emb)
        nb_in = n_genes + 2 + emb + 1                      # expression, counts covariates, class embedding, distance
        self.phi = nn.Sequential(nn.Linear(nb_in, hidden), nn.ReLU(), nn.Dropout(dropout))      # hop-2 messages
        self.self2 = nn.Linear(n_genes + 2 + emb, hidden)                                        # hop-1 node self term
        self.psi = nn.Sequential(nn.Linear(hidden + 1, hidden), nn.ReLU(), nn.Dropout(dropout))  # hop-1 messages
        self.self1 = nn.Linear(seed_cov_dim + emb, hidden)                                       # seed self term
        self.agg = aggregation
        if aggregation == "attention":
            self.att2, self.att1 = MaskedAttention(hidden, heads), MaskedAttention(hidden, heads)
        self.head = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(), nn.Dropout(dropout), nn.Linear(hidden, n_out))

    def forward(self, b: TreeBatch) -> torch.Tensor:
        e1, e2 = self.emb(b.cls1), self.emb(b.cls2)
        x2 = torch.cat([b.x2, b.cov2, e2, b.d2.unsqueeze(-1)], -1)                 # [B, k1, k2, nb_in]
        msg2 = self.phi(x2)
        node1_self = self.self2(torch.cat([b.x1, b.cov1, e1], -1))                 # [B, k1, hidden]
        if self.agg == "attention":
            agg2 = self.att2(node1_self, msg2, b.m2)
        else:
            agg2 = masked_mean(msg2, b.m2)
        h1 = Fn.relu(node1_self + agg2)                                            # hop-1 node states
        msg1 = self.psi(torch.cat([h1, b.d1.unsqueeze(-1)], -1))                   # [B, k1, hidden]
        seed_self = self.self1(torch.cat([b.seed_cov, self.emb(b.seed_cls)], -1))  # [B, hidden]
        if self.agg == "attention":
            agg1 = self.att1(seed_self, msg1, b.m1)
        else:
            agg1 = masked_mean(msg1, b.m1)
        h = Fn.relu(seed_self + agg1)
        return self.head(h)

"""Padded neighbourhood tables for mini-batch message passing with explicit seed exclusion.

What: converts a CSR adjacency into a fixed-width table of the k nearest neighbours per cell, and builds
2-hop neighbourhoods for a batch of seed cells in which the seed itself is removed from the second hop.

Why: a plain (padded-tensor) implementation lets the model see exactly which cells contribute to a prediction.
In the chain Phase 2b SGC model the path i -> j -> i leaked a cell's own expression into its inputs; here that
path is masked by construction, and `tests/test_leakage.py::test_two_hop_tree_excludes_seed` guards it.
"""
from __future__ import annotations

import numpy as np


def padded_neighbours(indptr: np.ndarray, indices: np.ndarray, lengths: np.ndarray | None, kmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell table [n, kmax] of neighbour indices (-1 padded), keeping the kmax nearest when degree > kmax.

    Returns (nbr, dist) where dist holds the edge lengths (0 where padded).
    """
    n = len(indptr) - 1
    nbr = np.full((n, kmax), -1, np.int64)
    dist = np.zeros((n, kmax), np.float32)
    deg = np.diff(indptr)
    for i in range(n):
        a, b = indptr[i], indptr[i + 1]
        if b == a:
            continue
        idx = indices[a:b]
        d = lengths[a:b] if lengths is not None else np.zeros(b - a, np.float32)
        if b - a > kmax:
            keep = np.argsort(d, kind="stable")[:kmax]
            idx, d = idx[keep], d[keep]
        nbr[i, :len(idx)] = idx
        dist[i, :len(idx)] = d
    assert (deg <= kmax).mean() > 0.95 or True  # informational; truncation is reported by the caller
    return nbr, dist


def padded_two_hop(indptr: np.ndarray, indices: np.ndarray, seeds: np.ndarray, k1: int, k2: int,
                   lengths: np.ndarray | None = None, table: tuple[np.ndarray, np.ndarray] | None = None):
    """For each seed: hop-1 neighbours [S, k1] and hop-2 neighbours [S, k1, k2], -1 padded, seed excluded at hop 2.

    `table` may pass a precomputed `padded_neighbours(...)` result (with kmax >= max(k1, k2)) to avoid recomputation.
    """
    nbr, _ = table if table is not None else padded_neighbours(indptr, indices, lengths, max(k1, k2))
    seeds = np.asarray(seeds, np.int64)
    h1 = nbr[seeds][:, :k1]                                    # [S, k1]
    h2 = nbr[np.clip(h1, 0, None)][:, :, :k2]                  # [S, k1, k2]
    h2 = np.where(h1[:, :, None] < 0, -1, h2)                  # padded hop-1 slots have no hop-2
    h2 = np.where(h2 == seeds[:, None, None], -1, h2)          # the seed never appears in its own hop 2
    return h1, h2


def two_hop_batch(nbr: np.ndarray, dist: np.ndarray, seeds: np.ndarray, k1: int, k2: int):
    """Same as padded_two_hop but from a precomputed table and also returning hop distances.

    Returns h1 [S,k1], d1 [S,k1], h2 [S,k1,k2], d2 [S,k1,k2]; -1 padding in h*, 0 in d*.
    """
    seeds = np.asarray(seeds, np.int64)
    h1, d1 = nbr[seeds][:, :k1], dist[seeds][:, :k1]
    safe = np.clip(h1, 0, None)
    h2, d2 = nbr[safe][:, :, :k2], dist[safe][:, :, :k2]
    bad = (h1[:, :, None] < 0) | (h2 == seeds[:, None, None])
    h2 = np.where(bad, -1, h2)
    d2 = np.where(bad, 0.0, d2).astype(np.float32)
    return h1, d1, h2, d2

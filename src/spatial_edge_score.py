"""Edge-based spatial ligand-receptor communication scoring.

THEORETICAL CONTEXT:
--------------------
Unlike cluster-marginal scoring (such as Squidpy ligrec / CellPhoneDB), which evaluates
(mean_L(sender) + mean_R(receiver)) / 2 over all cells irrespective of spatial proximity,
spatial-native cell-cell communication tools (e.g., COMMOT, Cang et al., Nature Methods 2023)
evaluate signaling potential directly over spatial neighborhood graphs.

This module computes the edge-based spatial LR interaction score:
  For directed edges (i -> j) in the spatial connectivity graph:
    edge_score = mean_{ (i->j) : cell_type[i] == sender, cell_type[j] == receiver } ( expr_L[i] * expr_R[j] )

It also tracks n_qualifying_edges (the count of spatial edges connecting sender to receiver cells).
If n_qualifying_edges == 0 for observed data, the interaction is flagged as 'insufficient_edges'.
"""
from __future__ import annotations

from typing import Sequence
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def get_spatial_graph_edges(adata: ad.AnnData) -> tuple[np.ndarray, np.ndarray]:
    """Extract directed edge pairs (row_indices, col_indices) from spatial graph.

    Uses adata.obsp['spatial_connectivities'] as-is without symmetrization.
    """
    if "spatial_connectivities" not in adata.obsp:
        raise KeyError("spatial_connectivities not found in adata.obsp")

    conn = adata.obsp["spatial_connectivities"]
    if sp.issparse(conn):
        rows, cols = conn.nonzero()
    else:
        rows, cols = np.nonzero(conn)

    return np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)


def precompute_edge_products(
    adata: ad.AnnData,
    pairs: Sequence[tuple[str, str]],
    row_indices: np.ndarray,
    col_indices: np.ndarray,
) -> dict[tuple[str, str], np.ndarray]:
    """Precompute edge expression products expr_L[row] * expr_R[col] for unique L-R pairs.

    Uses adata.raw.X if available, matching baseline pipeline behavior.
    """
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)

    unique_genes = sorted(set([g for pair in pairs for g in pair]))
    gene_indices = {g: var_names.index(g) for g in unique_genes}

    # Extract dense column per gene
    expr_cache: dict[str, np.ndarray] = {}
    for gene in unique_genes:
        idx = gene_indices[gene]
        col = raw_source.X[:, idx]
        if hasattr(col, "toarray"):
            expr_cache[gene] = col.toarray().ravel().astype(np.float64)
        else:
            expr_cache[gene] = np.asarray(col).ravel().astype(np.float64)

    edge_products: dict[tuple[str, str], np.ndarray] = {}
    for lig, rec in set(pairs):
        l_expr_sender = expr_cache[lig][row_indices]
        r_expr_receiver = expr_cache[rec][col_indices]
        edge_products[(lig, rec)] = l_expr_sender * r_expr_receiver

    return edge_products


def compute_edge_scores_for_labels(
    labels: np.ndarray,
    row_indices: np.ndarray,
    col_indices: np.ndarray,
    edge_products: dict[tuple[str, str], np.ndarray],
    interactions: list[tuple[str, str, str, str]],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute edge scores and qualifying edge counts for a given cell_type label assignment.

    Args:
        labels: 1D array of cell type labels (length = n_cells).
        row_indices: 1D array of sender cell indices for all spatial edges.
        col_indices: 1D array of receiver cell indices for all spatial edges.
        edge_products: Precomputed edge products dict keyed by (ligand, receptor).
        interactions: List of (sender, receiver, ligand, receptor) tuples.

    Returns:
        tuple (scores, edge_counts):
          scores: 1D float array of mean edge scores (np.nan if n_qualifying == 0).
          edge_counts: 1D int array of qualifying edge counts.
    """
    sender_types = labels[row_indices]
    receiver_types = labels[col_indices]

    n_interactions = len(interactions)
    scores = np.zeros(n_interactions, dtype=np.float64)
    edge_counts = np.zeros(n_interactions, dtype=np.int64)

    # Pre-filter masks per (sender_type, receiver_type) pair to avoid redundant boolean ops
    unique_type_pairs = set([(s, r) for s, r, _, _ in interactions])
    pair_masks: dict[tuple[str, str], np.ndarray] = {}
    for s, r in unique_type_pairs:
        pair_masks[(s, r)] = (sender_types == s) & (receiver_types == r)

    for i, (s, r, lig, rec) in enumerate(interactions):
        mask = pair_masks[(s, r)]
        n_edges = int(np.sum(mask))
        edge_counts[i] = n_edges

        if n_edges > 0:
            prod_arr = edge_products[(lig, rec)]
            scores[i] = float(np.mean(prod_arr[mask]))
        else:
            scores[i] = np.nan

    return scores, edge_counts

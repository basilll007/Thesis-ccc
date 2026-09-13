"""Spatial-null permutation testing for ligand-receptor cell-cell communication.

SPATIAL-NULL MODEL VARIANT AND THEORETICAL JUSTIFICATION:
---------------------------------------------------------
This module implements the Spatial Label-Permutation Null Model for spatial cell-cell
communication (CCC) inference, drawing directly on the methodology and evaluation
frameworks established in:
  1. CONCISE (Zhao et al., 2026, bioRxiv / PMC13320749: "Spatial co-expression and
     cell-cell communication inference from spatially resolved transcriptomics with CONCISE")
  2. SOAAR (Khatri et al., 2026, bioRxiv: "Spatial Autocorrelation Aware Resampling
     Improves Cell-Cell Interaction Inference in Spatial Transcriptomics Data")

In spatially resolved transcriptomics, naive ligand-receptor co-expression scores are
heavily susceptible to confounding by spatial autocorrelation (the tendency of proximate
cells and tissue regions to display correlated expression regardless of functional signaling)
and tissue-wide abundance differences across cell populations.

As shown by SOAAR (Khatri et al., 2026), testing spatial correlation without preserving
spatial structure leads to severe underestimation of null variance and inflated false-positive
rates. Similarly, CONCISE (Zhao et al., 2026) demonstrates via controlled permutation
scenarios that spatial graph structure and single-cell expression profiles must be held fixed
when constructing reference null distributions.

In our implementation:
  - The spatial coordinate embedding (adata.obsm['spatial']), generic 6-NN spatial graph
    (adata.obsp['spatial_connectivities']), and single-cell gene expression profiles (adata.raw.X)
    are held strictly constant across all permutations.
  - The cell-type label vector (adata.obs['cell_type']) is randomly permuted across the fixed cells.
  - For each permutation, sender and receiver cell-type mean expression values are recomputed to
    evaluate the identical interaction score metric as the baseline: (mean(L_sender) + mean(R_receiver)) / 2.

This null model specifically evaluates the null hypothesis that an observed high ligand-receptor
interaction score between a sender and receiver cell type is merely a reflection of global tissue
composition and non-specific spatial distribution, rather than genuine, cell-type-specific spatial
enrichment.
"""
from __future__ import annotations

from typing import Sequence
import numpy as np
import pandas as pd
import anndata as ad


def extract_relevant_expression(
    adata: ad.AnnData,
    genes: Sequence[str],
) -> tuple[np.ndarray, dict[str, int]]:
    """Extract dense expression matrix for the subset of unique interacting genes.

    Uses adata.raw.X if available (matching baseline_lr behavior), falling back to adata.X.
    """
    raw_source = adata.raw if adata.raw is not None else adata
    var_names = list(raw_source.var_names)
    gene_to_idx: dict[str, int] = {}
    col_indices: list[int] = []

    for gene in sorted(set(genes)):
        if gene in var_names:
            gene_to_idx[gene] = len(col_indices)
            col_indices.append(var_names.index(gene))
        else:
            raise KeyError(f"Gene '{gene}' not found in AnnData var_names")

    expr_sub = raw_source.X[:, col_indices]
    if hasattr(expr_sub, "toarray"):
        expr_dense = expr_sub.toarray()
    else:
        expr_dense = np.asarray(expr_sub)

    return expr_dense, gene_to_idx


def compute_scores_from_labels(
    labels: np.ndarray,
    expr_dense: np.ndarray,
    gene_to_idx: dict[str, int],
    interactions: list[tuple[str, str, str, str]],
    unique_cell_types: list[str],
) -> np.ndarray:
    """Compute interaction scores for an array of cell_type labels.

    Score definition matches Squidpy sq.gr.ligrec / baseline:
      score = (mean_sender(ligand) + mean_receiver(receptor)) / 2.0
    """
    n_interactions = len(interactions)
    n_cell_types = len(unique_cell_types)
    n_genes = len(gene_to_idx)

    # Precompute mean expression per cell type: shape (n_cell_types, n_genes)
    means = np.zeros((n_cell_types, n_genes), dtype=np.float64)
    ct_to_idx = {ct: i for i, ct in enumerate(unique_cell_types)}

    for ct, ct_idx in ct_to_idx.items():
        mask = (labels == ct)
        if np.any(mask):
            means[ct_idx, :] = expr_dense[mask, :].mean(axis=0)

    scores = np.zeros(n_interactions, dtype=np.float64)
    for i, (sender, receiver, ligand, receptor) in enumerate(interactions):
        s_idx = ct_to_idx[sender]
        r_idx = ct_to_idx[receiver]
        lig_col = gene_to_idx[ligand]
        rec_col = gene_to_idx[receptor]
        scores[i] = (means[s_idx, lig_col] + means[r_idx, rec_col]) / 2.0

    return scores


def generate_label_permutations(
    orig_labels: np.ndarray,
    n_perms: int,
    seed: int = 42,
) -> list[np.ndarray]:
    """Generate deterministic random permutations of cell_type labels.

    Args:
        orig_labels: 1D array of original cell_type labels.
        n_perms: Number of permutations to generate.
        seed: Random seed for deterministic reproducibility.

    Returns:
        List of 1D numpy arrays containing permuted labels.
    """
    rng = np.random.default_rng(seed)
    return [rng.permutation(orig_labels) for _ in range(n_perms)]


def run_spatial_null_permutations(
    adata: ad.AnnData,
    baseline_df: pd.DataFrame,
    n_perms: int,
    seed: int = 42,
) -> np.ndarray:
    """Run N spatial label permutations and return null distribution matrix.

    Args:
        adata: Processed AnnData containing cell_type and raw expression.
        baseline_df: DataFrame of baseline interactions with sender, receiver, ligand, receptor.
        n_perms: Number of permutations to execute.
        seed: Random seed for deterministic reproducibility.

    Returns:
        np.ndarray of shape (n_perms, len(baseline_df)) containing the null scores.
    """
    # Extract interactions list
    interactions: list[tuple[str, str, str, str]] = [
        (row["sender"], row["receiver"], row["ligand"], row["receptor"])
        for _, row in baseline_df.iterrows()
    ]

    all_genes = set()
    for _, _, lig, rec in interactions:
        all_genes.add(lig)
        all_genes.add(rec)

    expr_dense, gene_to_idx = extract_relevant_expression(adata, list(all_genes))
    orig_labels = np.asarray(adata.obs["cell_type"].astype(str).to_numpy())
    unique_cell_types = sorted(set(orig_labels))

    null_matrix = np.zeros((n_perms, len(interactions)), dtype=np.float64)
    permutations = generate_label_permutations(orig_labels, n_perms=n_perms, seed=seed)

    for p, permuted_labels in enumerate(permutations):
        null_matrix[p, :] = compute_scores_from_labels(
            labels=permuted_labels,
            expr_dense=expr_dense,
            gene_to_idx=gene_to_idx,
            interactions=interactions,
            unique_cell_types=unique_cell_types,
        )

    return null_matrix


def calculate_null_statistics(
    baseline_df: pd.DataFrame,
    null_matrix: np.ndarray,
) -> pd.DataFrame:
    """Compute null_mean, null_std, empirical_pvalue, z_score, and significance flag.

    Args:
        baseline_df: Original baseline DataFrame with sender, receiver, ligand, receptor, score.
        null_matrix: Array of shape (n_perms, len(baseline_df)) with permuted scores.

    Returns:
        pd.DataFrame with columns:
          [sender, receiver, ligand, receptor, observed_score, null_mean, null_std,
           empirical_pvalue, z_score, significant]
        sorted by observed_score descending.
    """
    results = []
    n_perms = null_matrix.shape[0]

    for i, row in baseline_df.iterrows():
        observed_score = float(row["score"])
        null_dist = null_matrix[:, i]
        null_mean = float(np.mean(null_dist))
        null_std = float(np.std(null_dist, ddof=1)) if n_perms > 1 else 0.0

        # Empirical p-value: fraction of null scores >= observed
        empirical_pvalue = float(np.mean(null_dist >= observed_score))

        # Z-score: (observed - null_mean) / null_std
        if null_std > 1e-12:
            z_score = float((observed_score - null_mean) / null_std)
        else:
            z_score = 0.0

        significant = bool(empirical_pvalue < 0.05)

        results.append({
            "sender": row["sender"],
            "receiver": row["receiver"],
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "observed_score": observed_score,
            "null_mean": null_mean,
            "null_std": null_std,
            "empirical_pvalue": empirical_pvalue,
            "z_score": z_score,
            "significant": significant,
        })

    out_df = pd.DataFrame(results)
    out_df = out_df.sort_values("observed_score", ascending=False).reset_index(drop=True)
    return out_df

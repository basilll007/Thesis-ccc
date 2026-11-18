"""Vectorised spatial cell-cell communication scores and dyadic baselines.

What this module does:
Computes the four canonical scores that dissect spatial communication hypotheses:
1. Compositional score C: (Lbar_s + Rbar_r) / 2 (Squidpy / CellPhoneDB style, non-spatial mean).
2. Edge score E: mean_{edges s->r} L_i * R_j (graph-restricted ligand-receptor product).
3. Matched dyadic baseline M: Lbar_s * Rbar_r (same product form as E, but evaluated without spatial restriction).
4. Proximity ratio Q: E / M (relative edge enrichment over the matched non-spatial expectation).

Why:
Standard spatial tools conflate two distinct phenomena:
(a) Score form (sum vs product): whether an interaction score reflects a multiplicative relation.
(b) Proximity (spatial restriction): whether cells in contact have higher co-expression than expected
    given cell-type abundance alone.
Comparing C, M, E, and Q isolates whether an apparent 'communication' signal is purely an abundance
artefact (C and M match E) or a true spatial contact enrichment (Q > 1 with empirical significance).

Conventions:
- In permutations with zero qualifying edges, E is set to 0.0 (hence Q = 0.0).
- Observed interactions with zero qualifying edges evaluate to NaN, status 'insufficient_edges'.
- Empirical p-values use the primary add-one convention: p = (k + 1) / (n + 1).
- Standardised z-scores use sample standard deviation (ddof=1).

Adapted from:
`src/pdac/ccc_audit.py` and `src/pdac/phase1_matched_baseline.py`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCORES: tuple[str, ...] = ("comp", "edge", "matched", "ratio")


@dataclass
class LRSetup:
    """Precomputed gene expression matrices and edge-wise ligand-receptor products."""
    lr_pairs: list[tuple[str, str]]
    genes: list[str]
    X: np.ndarray  # (n_cells, n_genes)
    lig_idx: np.ndarray
    rec_idx: np.ndarray
    prods: np.ndarray  # (n_pairs, n_edges)
    rows: np.ndarray  # directed edge sender indices
    cols: np.ndarray  # directed edge receiver indices


@dataclass
class ScoreResults:
    """Observed scores, null distributions, and summary statistics."""
    dataframe: pd.DataFrame
    null_matrices: dict[str, np.ndarray]


def prepare_lr_setup(
    adata: ad.AnnData,
    lr_pairs: Sequence[tuple[str, str]],
    rows: np.ndarray,
    cols: np.ndarray,
    use_raw: bool = True,
    dtype=np.float64,
) -> LRSetup:
    """Extract expression matrix and precompute edge-wise L_i * R_j products.

    Source file: `src/pdac/ccc_audit.py` (setup function).
    """
    pairs = list(lr_pairs)
    genes = sorted({g for p in pairs for g in p})

    # Extract expression
    src = adata.raw if (use_raw and adata.raw is not None) else adata
    var_names = list(src.var_names)
    for g in genes:
        if g not in var_names:
            raise ValueError(f"Gene '{g}' not found in AnnData var_names.")

    gene_indices = [var_names.index(g) for g in genes]
    X = src.X[:, gene_indices]
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=dtype)

    gi = {g: i for i, g in enumerate(genes)}
    lig_idx = np.array([gi[l] for l, _ in pairs], dtype=np.int64)
    rec_idx = np.array([gi[r] for _, r in pairs], dtype=np.int64)

    # Edge products
    n_pairs = len(pairs)
    n_edges = len(rows)
    prods = np.empty((n_pairs, n_edges), dtype=dtype)
    for p, (l, r) in enumerate(pairs):
        prods[p] = X[rows, gi[l]] * X[cols, gi[r]]

    return LRSetup(
        lr_pairs=pairs,
        genes=genes,
        X=X,
        lig_idx=lig_idx,
        rec_idx=rec_idx,
        prods=prods,
        rows=rows,
        cols=cols,
    )


def compute_scores_for_labels(
    labels: np.ndarray,
    n_types: int,
    setup: LRSetup,
    zero_fill_edge: bool = False,
) -> dict[str, np.ndarray]:
    """Compute C, E, M, Q and edge counts for all cell-type pairs under given labels.

    Source file: `src/pdac/ccc_audit.py` (scores_for_labels).
    """
    n_cells = labels.shape[0]
    onehot = sp.csr_matrix(
        (np.ones(n_cells), (np.arange(n_cells), labels)),
        shape=(n_cells, n_types),
    )
    counts = np.asarray(onehot.sum(axis=0)).ravel()
    sums = np.asarray(onehot.T @ setup.X)  # shape (n_types, n_genes)

    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / counts[:, None]
    means[counts == 0] = 0.0

    Lm = means[:, setup.lig_idx].T  # shape (n_pairs, n_types)
    Rm = means[:, setup.rec_idx].T  # shape (n_pairs, n_types)

    # C = (Lm + Rm) / 2
    comp = (Lm[:, :, None] + Rm[:, None, :]) / 2.0
    # M = Lm * Rm
    matched = Lm[:, :, None] * Rm[:, None, :]

    # Edge-based counts and sums
    if len(setup.rows) > 0:
        code = labels[setup.rows] * n_types + labels[setup.cols]
        ecount = np.bincount(code, minlength=n_types * n_types).reshape(n_types, n_types)
        n_pairs = len(setup.lr_pairs)
        esum = np.empty((n_pairs, n_types * n_types), dtype=float)
        for p in range(n_pairs):
            esum[p] = np.bincount(code, weights=setup.prods[p], minlength=n_types * n_types)
        esum = esum.reshape(n_pairs, n_types, n_types)
        with np.errstate(invalid="ignore", divide="ignore"):
            edge = esum / ecount[None, :, :]
    else:
        ecount = np.zeros((n_types, n_types), dtype=np.int64)
        edge = np.full((len(setup.lr_pairs), n_types, n_types), np.nan)

    if zero_fill_edge:
        edge = np.where(ecount[None] == 0, 0.0, edge)

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(matched > 0, edge / matched, np.nan)

    return {
        "comp": comp,
        "edge": edge,
        "matched": matched,
        "ratio": ratio,
        "edge_count": ecount,
    }


def compute_statistics(obs: np.ndarray, null: np.ndarray) -> dict[str, np.ndarray]:
    """Compute summary statistics (null mean, SD, z-score, add-one p-value) against null.

    Source file: `src/pdac/ccc_audit.py` (stats).
    """
    n_valid = np.sum(np.isfinite(null), axis=0)
    mean = np.nanmean(null, axis=0)
    sd = np.nanstd(null, axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 1e-12, (obs - mean) / sd, 0.0)

    ge = np.sum(np.nan_to_num(null, nan=-np.inf) >= obs[None, :], axis=0)
    p_plus1 = (ge + 1) / (np.maximum(n_valid, 1) + 1)
    p_legacy = ge / np.maximum(n_valid, 1)

    bad = ~np.isfinite(obs)
    for arr in (mean, sd, z, p_plus1, p_legacy):
        arr[bad] = np.nan

    return {
        "null_mean": mean,
        "null_sd": sd,
        "z": z,
        "p_plus1": p_plus1,
        "p_legacy": p_legacy,
        "n_null_valid": n_valid,
    }


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Vectorised Benjamini-Hochberg FDR correction.

    Source file: `src/pdac/ccc_audit.py` (bh).
    """
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    m = int(valid.sum())
    if m == 0:
        return q

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[valid] = out
    return q


def compute_ccc_scores(
    adata: ad.AnnData,
    lr_pairs: Sequence[tuple[str, str]],
    interactions: Sequence[tuple[str, str, str, str]],  # (sender, receiver, ligand, receptor)
    rows: np.ndarray,
    cols: np.ndarray,
    perms: Sequence[np.ndarray],
    use_raw: bool = True,
    alpha: float = 0.05,
) -> ScoreResults:
    """Compute full CCC audit scores under observed and permuted label vectors.

    Source file: `src/pdac/ccc_audit.py` (tidy, run_null).
    """
    cell_types = sorted(set(adata.obs["cell_type"].astype(str)))
    type_to_idx = {t: i for i, t in enumerate(cell_types)}
    n_types = len(cell_types)
    obs_labels = np.array([type_to_idx[x] for x in adata.obs["cell_type"].astype(str)], dtype=np.int64)

    setup = prepare_lr_setup(adata, lr_pairs, rows, cols, use_raw=use_raw)

    # Observed scores
    obs = compute_scores_for_labels(obs_labels, n_types, setup, zero_fill_edge=False)

    # Permuted scores
    n_perms = len(perms)
    n_pairs = len(setup.lr_pairs)
    null = {k: np.empty((n_perms, n_pairs, n_types, n_types), dtype=float) for k in SCORES}
    null_zero_edges = np.zeros((n_types, n_types), dtype=np.int64)

    for k, pl in enumerate(perms):
        res = compute_scores_for_labels(pl, n_types, setup, zero_fill_edge=True)
        for key in SCORES:
            null[key][k] = res[key]
        null_zero_edges += (res["edge_count"] == 0)

    # Tidy results for selected interactions
    pair_to_idx = {p: i for i, p in enumerate(setup.lr_pairs)}
    interaction_indices = [
        (pair_to_idx[(l, r)], type_to_idx[s], type_to_idx[rc])
        for s, rc, l, r in interactions
    ]

    stat_dfs = []
    null_out = {}
    for key in SCORES:
        obs_vals = np.array([obs[key][p, a, b] for p, a, b in interaction_indices], dtype=float)
        null_vals = np.stack([null[key][:, p, a, b] for p, a, b in interaction_indices], axis=1)
        null_out[key] = null_vals
        stats = compute_statistics(obs_vals, null_vals)
        stat_dfs.append(pd.DataFrame({
            f"{key}_obs": obs_vals,
            **{f"{key}_{k}": v for k, v in stats.items()}
        }))

    base = pd.DataFrame(interactions, columns=["sender", "receiver", "ligand", "receptor"])
    base["n_qualifying_edges"] = [int(obs["edge_count"][a, b]) for _, a, b in interaction_indices]
    base["n_perms_zero_edges"] = [int(null_zero_edges[a, b]) for _, a, b in interaction_indices]

    df = pd.concat([base] + stat_dfs, axis=1)
    df["edge_status"] = np.where(df["n_qualifying_edges"] == 0, "insufficient_edges", "evaluated")

    for key in SCORES:
        df[f"{key}_sig"] = df[f"{key}_p_plus1"] < alpha
        df[f"{key}_q_bh"] = benjamini_hochberg(df[f"{key}_p_plus1"].to_numpy())

    return ScoreResults(dataframe=df, null_matrices=null_out)

"""Vectorised CCC audit scores under a shared label-permutation null.

For a ligand L, receptor R, sender type s and receiver type r:

  compositional  C_sr = (Lbar_s + Rbar_r) / 2                 (Squidpy/CellPhoneDB-style; no graph)
  edge           E_sr = mean_{(i->j) in G, t(i)=s, t(j)=r} L_i * R_j   (graph-restricted product)
  matched        M_sr = Lbar_s * Rbar_r                        (non-spatial matched dyadic baseline:
                                                                same product form as E, no graph)
  ratio          Q_sr = E_sr / M_sr                            (edge enrichment over matched baseline)

C vs M differ only in score form (sum vs product); M vs E differ only in spatial restriction.
All four are evaluated on the same permuted label vectors, so their nulls are directly comparable.

Conventions (identical to src/run_spatial_edge_null.py so paper numbers are reproduced):
  * Expression is taken from adata.raw.X (log-normalised) when present.
  * In a permutation with zero qualifying edges, E is set to 0 (and therefore Q = 0).
  * Observed E with zero qualifying edges -> NaN, status 'insufficient_edges'.
  * null sd uses ddof=1; z = (obs - null_mean) / null_sd (0 if sd < 1e-12).
  * p_legacy = mean(null >= obs)  (as in existing files);  p_plus1 = (k + 1) / (n + 1)  (primary).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

SCORES = ("comp", "edge", "matched", "ratio")


def gene_matrix(adata, genes: list[str], use_raw: bool = True, dtype=np.float64) -> np.ndarray:
    src = adata.raw if (use_raw and adata.raw is not None) else adata
    names = list(src.var_names)
    idx = [names.index(g) for g in genes]
    X = src.X[:, idx]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return X.astype(dtype)


def graph_edges(conn) -> tuple[np.ndarray, np.ndarray]:
    conn = sp.csr_matrix(conn)
    rows, cols = conn.nonzero()
    return rows.astype(np.int64), cols.astype(np.int64)


@dataclass
class LRSetup:
    lr_pairs: list[tuple[str, str]]  # unique (ligand, receptor)
    genes: list[str]
    X: np.ndarray  # n_cells x n_genes
    lig_idx: np.ndarray
    rec_idx: np.ndarray
    prods: np.ndarray  # n_pairs x n_edges   L_i * R_j on each directed edge
    rows: np.ndarray
    cols: np.ndarray


def setup(adata, lr_pairs: list[tuple[str, str]], rows: np.ndarray, cols: np.ndarray,
          use_raw: bool = True, prod_dtype=np.float64) -> LRSetup:
    genes = sorted({g for p in lr_pairs for g in p})
    X = gene_matrix(adata, genes, use_raw=use_raw)
    gi = {g: k for k, g in enumerate(genes)}
    lig_idx = np.array([gi[l] for l, _ in lr_pairs])
    rec_idx = np.array([gi[r] for _, r in lr_pairs])
    prods = np.empty((len(lr_pairs), len(rows)), dtype=prod_dtype)
    for p, (l, r) in enumerate(lr_pairs):
        prods[p] = X[rows, gi[l]] * X[cols, gi[r]]
    return LRSetup(lr_pairs, genes, X, lig_idx, rec_idx, prods, rows, cols)


def scores_for_labels(lab: np.ndarray, T: int, S: LRSetup, zero_fill_edge: bool) -> dict[str, np.ndarray]:
    """All four scores for all (pair, sender, receiver): arrays of shape (P, T, T); plus edge counts (T, T)."""
    n = lab.shape[0]
    onehot = sp.csr_matrix((np.ones(n), (np.arange(n), lab)), shape=(n, T))
    counts = np.asarray(onehot.sum(axis=0)).ravel()
    sums = np.asarray(onehot.T @ S.X)  # T x G
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / counts[:, None]
    means[counts == 0] = 0.0
    Lm = means[:, S.lig_idx].T  # P x T
    Rm = means[:, S.rec_idx].T  # P x T
    comp = (Lm[:, :, None] + Rm[:, None, :]) / 2.0
    matched = Lm[:, :, None] * Rm[:, None, :]

    code = lab[S.rows] * T + lab[S.cols]
    ecount = np.bincount(code, minlength=T * T).reshape(T, T)
    P = len(S.lr_pairs)
    esum = np.empty((P, T * T))
    for p in range(P):
        esum[p] = np.bincount(code, weights=S.prods[p], minlength=T * T)
    esum = esum.reshape(P, T, T)
    with np.errstate(invalid="ignore", divide="ignore"):
        edge = esum / ecount[None, :, :]
    if zero_fill_edge:
        edge = np.where(ecount[None] == 0, 0.0, edge)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(matched > 0, edge / matched, np.nan)
    return {"comp": comp, "edge": edge, "matched": matched, "ratio": ratio, "edge_count": ecount}


def run_null(lab: np.ndarray, T: int, S: LRSetup, perms: list[np.ndarray]):
    obs = scores_for_labels(lab, T, S, zero_fill_edge=False)
    n = len(perms)
    P = len(S.lr_pairs)
    null = {k: np.empty((n, P, T, T)) for k in SCORES}
    null_zero_edges = np.zeros((T, T), dtype=np.int64)
    for k, pl in enumerate(perms):
        res = scores_for_labels(pl, T, S, zero_fill_edge=True)
        for key in SCORES:
            null[key][k] = res[key]
        null_zero_edges += (res["edge_count"] == 0)
    return obs, null, null_zero_edges


def stats(obs: np.ndarray, null: np.ndarray) -> dict[str, np.ndarray]:
    """obs: (K,), null: (n, K). NaN-safe: null NaNs are dropped per column (only ratio can produce them)."""
    n_valid = np.sum(np.isfinite(null), axis=0)
    mean = np.nanmean(null, axis=0)
    sd = np.nanstd(null, axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 1e-12, (obs - mean) / sd, 0.0)
    ge = np.sum(np.nan_to_num(null, nan=-np.inf) >= obs[None, :], axis=0)
    p_legacy = ge / np.maximum(n_valid, 1)
    p_plus1 = (ge + 1) / (n_valid + 1)
    bad = ~np.isfinite(obs)
    for arr in (mean, sd, z, p_legacy, p_plus1):
        arr[bad] = np.nan
    return {"null_mean": mean, "null_sd": sd, "z": z, "p_legacy": p_legacy, "p_plus1": p_plus1, "n_null_valid": n_valid}


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN entries are ignored and returned as NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def tidy(interactions: list[tuple[str, str, str, str]], types: list[str], S: LRSetup, obs, null, null_zero_edges,
         alpha: float = 0.05) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Flatten selected (sender, receiver, ligand, receptor) rows into one table + per-row null matrices."""
    ti = {t: k for k, t in enumerate(types)}
    pi = {p: k for k, p in enumerate(S.lr_pairs)}
    idx = [(pi[(l, r)], ti[s], ti[rc]) for s, rc, l, r in interactions]
    rows = []
    null_out: dict[str, np.ndarray] = {}
    for key in SCORES:
        o = np.array([obs[key][p, a, b] for p, a, b in idx], dtype=float)
        nm = np.stack([null[key][:, p, a, b] for p, a, b in idx], axis=1)
        null_out[key] = nm
        st = stats(o, nm)
        rows.append(pd.DataFrame({f"{key}_obs": o, **{f"{key}_{k}": v for k, v in st.items()}}))
    base = pd.DataFrame(interactions, columns=["sender", "receiver", "ligand", "receptor"])
    base["n_qualifying_edges"] = [int(obs["edge_count"][a, b]) for _, a, b in idx]
    base["n_perms_zero_edges"] = [int(null_zero_edges[a, b]) for _, a, b in idx]
    df = pd.concat([base] + rows, axis=1)
    df["edge_status"] = np.where(df["n_qualifying_edges"] == 0, "insufficient_edges", "evaluated")
    for key in SCORES:
        df[f"{key}_sig"] = df[f"{key}_p_plus1"] < alpha
        df[f"{key}_q_bh"] = bh(df[f"{key}_p_plus1"].to_numpy())
    return df, null_out

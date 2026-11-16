"""Shared feature/graph utilities for chain Phases 2-4 and 6 (cache readers, tile means, neighbour aggregates)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

TYPES = {"tumour": "epithelial_ductal_tumor", "CAF": "CAF_fibroblast", "T": "T_cell", "myeloid": "myeloid"}
TILE = 200.0


def programs_used(phase="phase1b") -> dict:
    return json.loads((cc.RESULTS / phase / "programs_used.json").read_text())["programs"]


def load_cache(section: str, sub="mendeley_v1b") -> dict:
    z = np.load(cc.CACHE / sub / f"{section}.npz")
    d = {k: z[k] for k in z.files}
    n = len(d["cell_id"])
    d["W"] = sp.csr_matrix((np.ones(len(d["knn_indices"]), np.float32), d["knn_indices"], d["knn_indptr"]), shape=(n, n))
    d["R"] = sp.csr_matrix((np.ones(len(d["r30_indices"]), np.float32), d["r30_indices"], d["r30_indptr"]), shape=(n, n))
    d["score_names"] = [str(x) for x in d["score_names"]]
    d["lineage_names"] = [str(x) for x in d["lineage_names"]]
    return d


def tile_codes(xy: np.ndarray, size: float = TILE, offset=(0.0, 0.0)) -> np.ndarray:
    t = np.floor((xy - xy.min(axis=0) + np.asarray(offset)) / size).astype(np.int64)
    return t[:, 0] * (t[:, 1].max() + 2) + t[:, 1]


def group_mean_loo(values: np.ndarray, group: np.ndarray, member: np.ndarray):
    """Leave-one-out mean of `values` over members sharing a group. Returns (loo_mean, has_other) for all cells.

    For non-members the mean over all members in the group (nothing to leave out) is returned.
    """
    g = np.unique(group, return_inverse=True)[1]
    ng = g.max() + 1
    v = np.where(member, values, 0.0)
    s = np.bincount(g, weights=v, minlength=ng)
    c = np.bincount(g, weights=member.astype(float), minlength=ng)
    s_i, c_i = s[g] - np.where(member, values, 0.0), c[g] - member
    has = c_i > 0
    out = np.where(has, s_i / np.maximum(c_i, 1), np.nan)
    return out, has


def group_mean(values: np.ndarray, group: np.ndarray, member: np.ndarray):
    g = np.unique(group, return_inverse=True)[1]
    ng = g.max() + 1
    s = np.bincount(g, weights=np.where(member, values, 0.0), minlength=ng)
    c = np.bincount(g, weights=member.astype(float), minlength=ng)
    out = np.where(c[g] > 0, s[g] / np.maximum(c[g], 1), np.nan)
    return out, c[g] > 0


def neighbour_mean(A: sp.csr_matrix, values: np.ndarray, member: np.ndarray):
    """Mean of `values` (n x p or n) over neighbours that are members. Returns (mean, count)."""
    m = member.astype(np.float32)
    cnt = A @ m
    v = values if values.ndim == 2 else values[:, None]
    s = A @ (v * m[:, None])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s / cnt[:, None]
    mean[cnt == 0] = np.nan
    return (mean if values.ndim == 2 else mean[:, 0]), cnt

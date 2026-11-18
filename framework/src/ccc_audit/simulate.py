"""Synthetic spatial tissues for unit tests and calibration (no real data needed).

`synthetic_tissue` draws cell positions uniformly, assigns three cell types from a smooth random field (so
type composition is spatially structured, like stroma vs tumour nests), and draws counts for a small panel
(ligand L, receptor R, 4 markers per type, 8 background genes) as independent negative-binomial variables
with type-specific means. Optional planted effects:
  * contact: R mean x f in receivers touching an L+ sender (true contact signalling);
  * spill: move a fraction of every sender's counts into its contact neighbours (segmentation error).
"""
from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors

TYPES = ("Sender", "Receiver", "Other")


def nb(rng: np.random.Generator, mean: np.ndarray, theta: float = 3.0) -> np.ndarray:
    lam = rng.gamma(theta, np.asarray(mean, float) / theta)
    return rng.poisson(lam)


def symmetric_knn(xy: np.ndarray, k: int = 6) -> tuple[np.ndarray, np.ndarray]:
    n = xy.shape[0]
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy)
    a = np.repeat(np.arange(n), k)
    b = idx[:, 1:].ravel()
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    key = np.unique(lo * n + hi)
    return key // n, key % n


def synthetic_tissue(seed: int, n_cells: int = 800, size_um: float = 600.0, contact_effect: float = 1.0,
                     spill_frac: float = 0.0, theta: float = 3.0):
    """Return (xy, labels, counts, genes)."""
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, size_um, size=(n_cells, 2))
    ph = rng.uniform(0, 2 * np.pi, 4)
    field = (np.sin(xy[:, 0] / 90 + ph[0]) + np.cos(xy[:, 1] / 110 + ph[1])
             + 0.5 * np.sin((xy[:, 0] + xy[:, 1]) / 70 + ph[2]) + rng.normal(0, 0.6, n_cells))
    q = np.quantile(field, [0.4, 0.7])
    lab = np.where(field < q[0], 0, np.where(field < q[1], 1, 2))
    labels = np.array(TYPES)[lab]
    genes = ["L", "R"] + [f"{t}_m{j}" for t in TYPES for j in range(4)] + [f"bg{j}" for j in range(8)]
    mu = np.zeros((3, len(genes)))
    mu[:, 0] = [6.0, 1.0, 1.0]       # ligand mostly in senders
    mu[:, 1] = [5.0, 2.5, 1.0]       # receptor in receivers, higher in senders (so spillover can mimic signal)
    for t in range(3):
        for j in range(4):
            col = 2 + 4 * t + j
            mu[:, col] = 0.2
            mu[t, col] = 8.0
    mu[:, 2 + 12:] = 3.0
    mean = mu[lab]
    sf = rng.lognormal(0, 0.3, n_cells)[:, None]
    mean = mean * sf
    pa, pb = symmetric_knn(xy)
    if contact_effect != 1.0:
        L0 = nb(rng, mean[:, 0], theta)
        s_pos = (lab == 0) & (L0 > 0)
        exp_ = np.zeros(n_cells, bool)
        f1 = s_pos[pa] & (lab[pb] == 1)
        f2 = s_pos[pb] & (lab[pa] == 1)
        exp_[pb[f1]] = True
        exp_[pa[f2]] = True
        mean[exp_, 1] *= contact_effect
        counts = nb(rng, mean, theta)
        counts[:, 0] = L0
    else:
        counts = nb(rng, mean, theta)
    if spill_frac > 0:
        counts = spill(counts, lab == 0, pa, pb, spill_frac, rng)
    return xy, labels, counts.astype(float), genes


def spill(counts: np.ndarray, sender_mask: np.ndarray, pa: np.ndarray, pb: np.ndarray, frac: float,
          rng: np.random.Generator) -> np.ndarray:
    """Move `frac` of each sender's transcripts (binomial thinning) to its contact neighbours (uniform split)."""
    counts = counts.astype(np.int64).copy()
    n = counts.shape[0]
    nbrs: list[list[int]] = [[] for _ in range(n)]
    for a, b in zip(pa.tolist(), pb.tolist()):
        nbrs[a].append(b)
        nbrs[b].append(a)
    add = np.zeros_like(counts)
    for i in np.flatnonzero(sender_mask):
        if not nbrs[i]:
            continue
        moved = rng.binomial(counts[i], frac)
        counts[i] -= moved
        nb_ = np.array(nbrs[i])
        for g in np.flatnonzero(moved):
            dest = rng.integers(0, len(nb_), size=moved[g])
            np.add.at(add[:, g], nb_[dest], 1)
    return counts + add

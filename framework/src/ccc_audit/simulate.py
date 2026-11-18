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


# ----------------------------------------------------------------------------------------------
# Semi-synthetic benchmark on REAL geometry and REAL cell types (plan: framework/plans/v2_phase4.md)
# ----------------------------------------------------------------------------------------------
class SemiSyntheticSection:
    """Real section geometry + real cell-type labels; panel counts simulated from real type-level means.

    Panel = ligand, receptor, top `n_markers` markers per cell type (types with >= 20 cells) and `n_bg`
    high-mean background genes. mu[type, gene] = real mean raw count of the gene in that type; each cell keeps
    its real library-size factor (total counts / mean total of its type). Counts ~ NB(sf * mu, theta).
    """

    def __init__(self, xy, labels, counts, gene_names, ligand, receptor, sender, receiver,
                 n_markers: int = 5, n_bg: int = 10, theta: float = 3.0, k: int = 6):
        import scipy.sparse as sp
        self.xy = np.asarray(xy, float)
        self.labels = np.asarray(labels).astype(str)
        self.sender, self.receiver, self.ligand, self.receptor = sender, receiver, ligand, receptor
        self.theta = theta
        C = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
        genes = [str(g) for g in gene_names]
        gi = {g: i for i, g in enumerate(genes)}
        tot = np.asarray(C.sum(1)).ravel()
        types = sorted(set(self.labels))
        lin = (sp.diags(1e4 / np.maximum(tot, 1)) @ C).tocsr()
        tmeans_lin = np.vstack([np.asarray(lin[self.labels == t].mean(0)).ravel() for t in types])
        tmeans_cnt = np.vstack([np.asarray(C[self.labels == t].mean(0)).ravel() for t in types])
        excl = {gi[ligand], gi[receptor]}
        panel = [gi[ligand], gi[receptor]]
        for ti, t in enumerate(types):
            if (self.labels == t).sum() < 20:
                continue
            other = np.asarray(lin[self.labels != t].mean(0)).ravel()
            ratio = tmeans_lin[ti] / np.maximum(other, 1e-3)
            ok = np.flatnonzero(tmeans_cnt[ti] >= 0.5)
            ok = [g for g in ok[np.argsort(-ratio[ok])] if g not in excl and g not in panel]
            panel += ok[:n_markers]
        overall = np.asarray(C.mean(0)).ravel()
        bg = [g for g in np.argsort(-overall) if g not in panel and g not in excl][:n_bg]
        panel += bg
        self.panel_idx = np.array(panel)
        self.genes = [genes[g] for g in panel]
        self.types = types
        self.lab = np.array([types.index(x) for x in self.labels])
        self.mu = tmeans_cnt[:, self.panel_idx]
        tmean_tot = np.array([tot[self.labels == t].mean() for t in types])
        self.sf = tot / tmean_tot[self.lab]
        self.pa, self.pb = symmetric_knn(self.xy, k)
        n = len(self.labels)
        deg = np.bincount(np.concatenate([self.pa, self.pb]), minlength=n)
        order = np.argsort(np.concatenate([self.pa, self.pb]), kind="stable")
        self.nbr = np.concatenate([self.pb, self.pa])[order]
        self.nbr_start = np.concatenate([[0], np.cumsum(deg)[:-1]])
        self.deg = deg
        self.s_mask = self.labels == sender
        self.r_mask = self.labels == receiver

    def subset(self, keep: np.ndarray) -> "SemiSyntheticSection":
        """Same means/panel; geometry restricted to `keep` cells (used for receiver subsampling)."""
        new = object.__new__(SemiSyntheticSection)
        new.__dict__.update(self.__dict__)
        new.xy = self.xy[keep]
        new.labels = self.labels[keep]
        new.lab = self.lab[keep]
        new.sf = self.sf[keep]
        new.pa, new.pb = symmetric_knn(new.xy)
        n = len(new.labels)
        deg = np.bincount(np.concatenate([new.pa, new.pb]), minlength=n)
        order = np.argsort(np.concatenate([new.pa, new.pb]), kind="stable")
        new.nbr = np.concatenate([new.pb, new.pa])[order]
        new.nbr_start = np.concatenate([[0], np.cumsum(deg)[:-1]])
        new.deg = deg
        new.s_mask = new.labels == self.sender
        new.r_mask = new.labels == self.receiver
        return new

    def _exposed(self, L_counts: np.ndarray) -> np.ndarray:
        s_pos = self.s_mask & (L_counts > 0)
        ex = np.zeros(len(self.labels), bool)
        f1 = s_pos[self.pa] & self.r_mask[self.pb]
        f2 = s_pos[self.pb] & self.r_mask[self.pa]
        ex[self.pb[f1]] = True
        ex[self.pa[f2]] = True
        return ex

    def draw(self, scenario: str, level: float | None, seed: int) -> np.ndarray:
        """Simulate panel counts for one replicate. Scenarios: S0, S1 (level = fold), S2 (level = fold),
        S3 (level = tile um; fold 2), S4 (level = fraction moved)."""
        rng = np.random.default_rng(seed)
        mean = self.mu[self.lab] * self.sf[:, None]
        if scenario == "S2":
            mean[self.s_mask, 0] *= level
            mean[self.r_mask, 1] *= level
        if scenario == "S3":
            off = rng.uniform(0, level, 2)
            g = np.floor((self.xy - self.xy.min(0) + off) / level).astype(np.int64)
            _, code = np.unique(g[:, 0] * (g[:, 1].max() + 1) + g[:, 1], return_inverse=True)
            high = rng.random(code.max() + 1) < 0.5
            h = high[code]
            mean[self.s_mask & h, 0] *= 2.0
            mean[self.r_mask & h, 1] *= 2.0
        if scenario == "S1":
            L = nb(rng, mean[:, 0], self.theta)
            ex = self._exposed(L)
            mean[ex, 1] *= level
            counts = nb(rng, mean, self.theta)
            counts[:, 0] = L
        else:
            counts = nb(rng, mean, self.theta)
        if scenario == "S4":
            counts = self.spill_fast(counts, level, rng)
        return counts

    def spill_fast(self, counts: np.ndarray, frac: float, rng: np.random.Generator) -> np.ndarray:
        """Move `frac` of each sender cell's transcripts to its contact neighbours (each moved transcript goes to
        a uniformly chosen neighbour). Conserves total transcripts."""
        counts = counts.astype(np.int64)
        S = np.flatnonzero(self.s_mask & (self.deg > 0))
        sub = counts[S]
        moved = rng.binomial(sub, frac)
        counts[S] = sub - moved
        ci, gi = np.nonzero(moved)
        reps = moved[ci, gi]
        cell = np.repeat(S[ci], reps)
        gene = np.repeat(gi, reps)
        pick = (rng.random(len(cell)) * self.deg[cell]).astype(np.int64)
        dest = self.nbr[self.nbr_start[cell] + pick]
        G = counts.shape[1]
        add = np.bincount(dest * G + gene, minlength=counts.shape[0] * G).reshape(counts.shape)
        return counts + add

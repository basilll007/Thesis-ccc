"""Feature construction with provenance, and the leakage guard.

What: builds per-section feature blocks (own covariates, neighbourhood composition, neighbour program means,
neighbour gene means) for receiver cells, each column tagged with where its values come from (`self_covariate`,
`tile`, `neighbour`) and which genes/programs it involves. `assert_no_target_leak` refuses any feature matrix in
which a self-sourced column involves a target gene/program, or is numerically a copy of a target.

Why: the earlier SGC bug (chain Phase 2b) leaked the cell's own score into its inputs through a return path.
Provenance tags plus a numeric probe make that class of error fail loudly instead of silently inflating R2.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

TYPED = {"tumour": "epithelial_ductal_tumor", "CAF": "CAF_fibroblast", "T": "T_cell", "myeloid": "myeloid"}


class LeakError(RuntimeError):
    """Raised when a receiver cell's own target information is found among its inputs."""


@dataclass
class Provenance:
    source: str                      # 'self_covariate' | 'tile' | 'neighbour'
    genes: frozenset = frozenset()   # gene symbols involved (empty for count/density covariates)
    programs: frozenset = frozenset()


@dataclass
class FeatureSet:
    """Column-aligned feature blocks for one set of receiver rows."""
    rows: np.ndarray                                   # indices into the section (receiver cells)
    blocks: dict[str, np.ndarray] = field(default_factory=dict)
    names: dict[str, list[str]] = field(default_factory=dict)
    prov: dict[str, list[Provenance]] = field(default_factory=dict)

    def add(self, block: str, X: np.ndarray, names: list[str], prov: list[Provenance]) -> None:
        if X.shape != (len(self.rows), len(names)) or len(prov) != len(names):
            raise ValueError(f"block {block}: shape {X.shape} vs {len(self.rows)} rows x {len(names)} names")
        self.blocks[block], self.names[block], self.prov[block] = X.astype(np.float32), names, prov

    def matrix(self, blocks: list[str] | None = None) -> np.ndarray:
        blocks = blocks or list(self.blocks)
        return np.column_stack([self.blocks[b] for b in blocks]) if blocks else np.zeros((len(self.rows), 0), np.float32)

    def column_names(self, blocks: list[str] | None = None) -> list[str]:
        return [n for b in (blocks or list(self.blocks)) for n in self.names[b]]

    def provenance(self, blocks: list[str] | None = None) -> list[Provenance]:
        return [p for b in (blocks or list(self.blocks)) for p in self.prov[b]]


def assert_no_target_leak(fs: FeatureSet, blocks: list[str] | None, target: np.ndarray, target_genes: set[str],
                          target_programs: set[str], r_max: float = 0.999) -> None:
    """Refuse a feature matrix if the receiver's own target information is in it.

    Rule 1 (provenance): any self-sourced column (source == 'self_covariate') whose genes or programs intersect the
    target's is a leak. Rule 2 (numeric probe): any self-sourced column with |r| > r_max against any target column
    is a leak, whatever its tag says. Tile and neighbour columns are other cells' values and are allowed.
    """
    names, prov = fs.column_names(blocks), fs.provenance(blocks)
    X = fs.matrix(blocks)
    T = target.reshape(len(target), -1).astype(np.float64)
    bad = []
    for j, (n, p) in enumerate(zip(names, prov)):
        if p.source != "self_covariate":
            continue
        if p.genes & target_genes or p.programs & target_programs:
            bad.append(f"{n}: self-sourced and involves target genes/programs")
            continue
        x = X[:, j].astype(np.float64)
        if x.std() == 0:
            continue
        xz = (x - x.mean()) / x.std()
        for k in range(T.shape[1]):
            t = T[:, k]
            if t.std() == 0:
                continue
            r = float(np.abs(np.mean(xz * (t - t.mean()) / t.std())))
            if r > r_max:
                bad.append(f"{n}: self-sourced and |r|={r:.4f} with target column {k}")
                break
    if bad:
        raise LeakError("target leakage detected in inputs:\n  " + "\n  ".join(bad))


# ----------------------------------------------------------------------------- per-section helpers
def csr(d: dict, key: str, n: int) -> sp.csr_matrix:
    return sp.csr_matrix((np.ones(len(d[f"{key}_indices"]), np.float32), d[f"{key}_indices"], d[f"{key}_indptr"]),
                         shape=(n, n))


def neighbour_mean(A: sp.csr_matrix, values: np.ndarray, member: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean of `values` (n x p) over neighbours that are members; NaN where none. Returns (mean, count)."""
    m = member.astype(np.float32)
    cnt = np.asarray(A @ m).ravel()
    s = A @ (values.astype(np.float32) * m[:, None])
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.asarray(s) / cnt[:, None]
    out[cnt == 0] = np.nan
    return out, cnt


def tile_loo_mean(values: np.ndarray, tile: np.ndarray, member: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Leave-one-out mean of `values` (n x p) over member cells sharing a tile; NaN where no other member."""
    g = np.unique(tile, return_inverse=True)[1]
    ng = g.max() + 1
    v = values.astype(np.float64) * member[:, None]
    s = np.vstack([np.bincount(g, weights=v[:, j], minlength=ng) for j in range(v.shape[1])]).T[g] - v
    c = np.bincount(g, weights=member.astype(float), minlength=ng)[g] - member
    with np.errstate(invalid="ignore", divide="ignore"):
        out = s / c[:, None]
    out[c <= 0] = np.nan
    return out, c > 0


def build_section_features(d: dict, receiver_class: str, target_programs: list[str], target_genes: list[str],
                           classes: list[str], include_gene_block: bool = True) -> tuple[FeatureSet, np.ndarray]:
    """All feature blocks for the receiver cells of one section; also returns the target matrix (programs + genes).

    Targets are returned separately and never enter a block. Missing neighbour means are filled with the section
    mean of that class (absence is carried by composition), and missing tile means with NaN (imputed from the
    training set by the caller) with an indicator column.
    """
    lab = d["label"].astype(str)
    n = len(lab)
    rows = np.where(lab == receiver_class)[0]
    fs = FeatureSet(rows=rows)
    W, R = csr(d, "knn", n), csr(d, "r30", n)
    S, snames = d["scores"].astype(np.float32), [str(x) for x in d["score_names"]]
    genes = [str(g) for g in d["genes"]]
    gidx = {g: i for i, g in enumerate(genes)}
    Xlog = d["X_log"]
    # targets
    tcols = [S[rows, snames.index(p)] for p in target_programs] + [Xlog[rows, gidx[g]].astype(np.float32) for g in target_genes]
    target = np.column_stack(tcols)
    tnames = target_programs + target_genes
    # own covariates
    deg = np.diff(d["knn_indptr"])
    mean_len = np.bincount(np.repeat(np.arange(n), deg), weights=d["knn_length"], minlength=n) / np.maximum(deg, 1)
    dens = np.diff(d["r30_indptr"]).astype(np.float32)
    own = np.column_stack([np.log1p(d["total_counts"][rows]), np.log1p(d["n_genes_detected"][rows]), dens[rows], mean_len[rows]])
    fs.add("own", own, ["log1p_total_counts", "log1p_n_genes", "density_r30", "mean_knn_edge_len"],
           [Provenance("self_covariate")] * 4)
    # zonation control: LOO tile mean of each target over same-type cells
    member = lab == receiver_class
    tm, has = tile_loo_mean(target_full(d, snames, gidx, target_programs, target_genes), d["tile"], member)
    tile_block = np.column_stack([tm[rows], (~has[rows]).astype(np.float32)[:, None]])
    fs.add("tile", tile_block, [f"tile_loo_mean__{t}" for t in tnames] + ["tile_empty"],
           [Provenance("tile", genes=frozenset(g for g in [t] if g in gidx),
                       programs=frozenset(p for p in [t] if p in snames)) for t in tnames] + [Provenance("tile")])
    # composition
    onehot = np.zeros((n, len(classes)), np.float32)
    onehot[np.arange(n), [classes.index(x) for x in lab]] = 1
    comp6 = np.asarray(W @ onehot) / np.maximum(deg, 1)[:, None]
    compr = np.asarray(R @ onehot)
    fs.add("composition", np.column_stack([comp6[rows], compr[rows]]),
           [f"knn_frac__{c}" for c in classes] + [f"r30_count__{c}" for c in classes],
           [Provenance("neighbour")] * (2 * len(classes)))
    # neighbour programs per typed class (knn and r30)
    cols, names, prov = [], [], []
    for t, cls in TYPED.items():
        mem = lab == cls
        fill = S[mem].mean(axis=0) if mem.any() else np.zeros(S.shape[1], np.float32)
        for gname, A in (("knn", W), ("r30", R)):
            m, _ = neighbour_mean(A, S, mem)
            m = np.where(np.isnan(m), fill[None, :], m)
            cols.append(m[rows])
            names += [f"nb_{gname}_{t}__{s}" for s in snames]
            prov += [Provenance("neighbour", programs=frozenset([s])) for s in snames]
    fs.add("nb_programs", np.column_stack(cols), names, prov)
    if include_gene_block:
        cols, names, prov = [], [], []
        typed_mask = np.isin(lab, list(TYPED.values()))
        for t, mem in list((t, lab == c) for t, c in TYPED.items()) + [("other", ~typed_mask)]:
            fill = Xlog[mem].astype(np.float32).mean(axis=0) if mem.any() else np.zeros(len(genes), np.float32)
            m, _ = neighbour_mean(W, Xlog.astype(np.float32), mem)
            m = np.where(np.isnan(m), fill[None, :], m)
            cols.append(m[rows])
            names += [f"nb_knn_{t}__{g}" for g in genes]
            prov += [Provenance("neighbour", genes=frozenset([g])) for g in genes]
        fs.add("nb_genes", np.column_stack(cols), names, prov)
    return fs, target


def target_full(d: dict, snames: list[str], gidx: dict, programs: list[str], genes: list[str]) -> np.ndarray:
    """Target columns for ALL cells of a section (needed for the tile LOO means)."""
    S = d["scores"]
    return np.column_stack([S[:, snames.index(p)] for p in programs] + [d["X_log"][:, gidx[g]].astype(np.float32) for g in genes])


def fit_scaler(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column mean/sd from TRAINING rows only (constant columns get sd 1); NaNs ignored in the fit."""
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    mu[~np.isfinite(mu)] = 0.0
    return mu.astype(np.float32), sd.astype(np.float32)


def apply_scaler(X: np.ndarray, mu: np.ndarray, sd: np.ndarray) -> np.ndarray:
    """Standardise; NaN (empty tile) is imputed with the training mean, i.e. 0 after scaling."""
    Z = (X - mu) / sd
    return np.nan_to_num(Z, nan=0.0).astype(np.float32)

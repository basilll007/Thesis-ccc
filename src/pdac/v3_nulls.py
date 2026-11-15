"""v3 shared helpers: spatially constrained label-permutation nulls and section loading.

Nulls (all return a list of int label arrays, same length as `lab`; global type counts preserved exactly):
  global  : v2 reference null (spatial_null.generate_label_permutations, seed 42) -> identical draws to Phase B3.
  block<s>: labels permuted independently within square tiles of side s um anchored at the section min x/y.
  dist    : labels permuted within signed-distance bands from the tumour/stroma boundary of the same annotation.
Declared in docs/v3_plans/phase1_contact_vs_zonation.md before any result was computed.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, RESULTS_PDAC, SEED  # noqa: E402

import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import knn_edges, qc  # noqa: E402
from spatial_null import generate_label_permutations  # noqa: E402

TILE_SIZES = (50.0, 100.0, 200.0)
DIST_BANDS = np.array([0.0, 25.0, 50.0, 100.0, 200.0, 400.0, np.inf])
TUMOUR_CLASS = "epithelial_ductal_tumor"
MIN_COMPARTMENT = 50
NULL_ORDER = ["global", "block200", "block100", "block50", "dist"]
ANNOTATIONS = ["consensus_T", "marker", "author_mapped"]


def load_section(section: str, log: logging.Logger):
    """Load + QC + normalise one Mendeley PDAC-SMI section; attach v2 label columns; return (adata, meta, labels_df, pairs, det)."""
    adata, meta = pdac_loaders.load(f"mendeley_pdac_smi:{section}", DATA_DIR, log)
    adata, _, _ = qc(adata, log)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    lab_dir = RESULTS_PDAC / "phaseB" / "sections"
    labs = pd.read_csv(lab_dir / f"{section}_labels.csv.gz", usecols=["cell_id", "marker", "author_mapped"])
    cons = pd.read_csv(lab_dir / f"{section}_consensus_labels.csv.gz", usecols=["cell_id", "consensus_T"])
    ids = adata.obs_names.astype(str).to_numpy()
    if not (labs.cell_id.astype(str).to_numpy() == ids).all() or not (cons.cell_id.astype(str).to_numpy() == ids).all():
        raise ValueError(f"cell order mismatch in section {section}")
    labs["consensus_T"] = cons["consensus_T"].to_numpy()
    vn = list(adata.var_names)
    pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
    det = {g: int((adata.layers["counts"][:, vn.index(g)].toarray().ravel() > 0).sum()) for g in {x for p in pairs for x in p}}
    return adata, meta, labs, pairs, det


def encode(labels_str: np.ndarray):
    types = sorted(set(labels_str))
    ti = {t: k for k, t in enumerate(types)}
    return types, np.array([ti[x] for x in labels_str], dtype=np.int64)


def global_perms(labels_str: np.ndarray, types: list[str], n_perms: int, seed: int = SEED) -> list[np.ndarray]:
    ti = {t: k for k, t in enumerate(types)}
    return [np.array([ti[x] for x in p], dtype=np.int64) for p in generate_label_permutations(labels_str, n_perms=n_perms, seed=seed)]


def _within_group_perms(lab: np.ndarray, group: np.ndarray, n_perms: int, seed: int) -> list[np.ndarray]:
    """Permute labels independently within each group code (vectorised: lexsort by (group, random key))."""
    rng = np.random.default_rng(seed)
    pos = np.argsort(group, kind="stable")  # cell positions grouped by group id
    out = []
    for _ in range(n_perms):
        order = np.lexsort((rng.random(lab.shape[0]), group))  # grouped by group, random order inside
        new = np.empty_like(lab)
        new[pos] = lab[order]
        out.append(new)
    return out


def tile_codes(xy: np.ndarray, size: float) -> np.ndarray:
    t = np.floor((xy - xy.min(axis=0)) / size).astype(np.int64)
    return t[:, 0] * (t[:, 1].max() + 1) + t[:, 1]


def block_perms(lab: np.ndarray, xy: np.ndarray, size: float, n_perms: int, seed: int = SEED):
    codes = tile_codes(xy, size)
    _, cnt = np.unique(codes, return_counts=True)
    info = {"tile_um": size, "n_tiles": int(len(cnt)), "median_cells_per_tile": float(np.median(cnt)),
            "n_singleton_tiles": int((cnt == 1).sum())}
    return _within_group_perms(lab, codes, n_perms, seed), info


def signed_boundary_distance(xy: np.ndarray, is_tumour: np.ndarray) -> np.ndarray:
    d = np.empty(xy.shape[0])
    nn_t = NearestNeighbors(n_neighbors=1).fit(xy[is_tumour])
    nn_s = NearestNeighbors(n_neighbors=1).fit(xy[~is_tumour])
    d[~is_tumour] = nn_t.kneighbors(xy[~is_tumour])[0].ravel()
    d[is_tumour] = -nn_s.kneighbors(xy[is_tumour])[0].ravel()
    return d


def dist_perms(lab: np.ndarray, types: list[str], xy: np.ndarray, n_perms: int, seed: int = SEED):
    """Distance-stratified null; returns (perms or None, info). None when the tumour compartment is too small."""
    if TUMOUR_CLASS not in types:
        return None, {"skipped": "no tumour class"}
    is_t = lab == types.index(TUMOUR_CLASS)
    nt, ns = int(is_t.sum()), int((~is_t).sum())
    if nt < MIN_COMPARTMENT or ns < MIN_COMPARTMENT:
        return None, {"skipped": f"compartment too small (tumour {nt}, non-tumour {ns})"}
    d = signed_boundary_distance(xy, is_t)
    band = np.searchsorted(DIST_BANDS, np.abs(d), side="right") - 1  # 0..5
    codes = band + 10 * is_t.astype(np.int64)
    _, cnt = np.unique(codes, return_counts=True)
    info = {"n_tumour": nt, "n_non_tumour": ns, "n_bands_used": int(len(cnt)), "min_band_cells": int(cnt.min()),
            "bands_um": DIST_BANDS[:-1].tolist(), "median_abs_distance_um": float(np.median(np.abs(d)))}
    return _within_group_perms(lab, codes, n_perms, seed), info


def make_perms(null: str, labels_str: np.ndarray, types: list[str], lab: np.ndarray, xy: np.ndarray, n_perms: int,
               seed: int = SEED):
    if null == "global":
        return global_perms(labels_str, types, n_perms, seed), {"null": "global"}
    if null.startswith("block"):
        return block_perms(lab, xy, float(null[5:]), n_perms, seed)
    if null == "dist":
        return dist_perms(lab, types, xy, n_perms, seed)
    raise ValueError(null)


def section_to_tumour() -> dict[str, str]:
    info = pd.read_csv(RESULTS_PDAC / "phase3" / "cross_section_section_info.csv")
    info["tumour"] = info.donor.astype(str).str.replace(r"(_2nd|_B)$", "", regex=True)
    return dict(zip(info.section, info.tumour))


def nonbg_axes() -> list[tuple[str, str]]:
    bgf = pd.read_csv(RESULTS_PDAC / "phase3" / "background_flag_per_axis.csv")
    return [tuple(x) for x in bgf.loc[~bgf.background_level_gene_diagnostic, ["ligand", "receptor"]].to_numpy()]

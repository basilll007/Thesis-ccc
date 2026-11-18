"""Spatially constrained and unconstrained null models for CCC audits.

What this module does:
Generates permutation null draws for hypothesis testing:
1. Global label permutation: shuffles cell-type labels across all cells in the section,
   preserving global cell-type abundance exactly.
2. Block label permutation: shuffles cell-type labels strictly within local spatial tiles
   (e.g., 50 um, 100 um, 200 um) with random coordinate offset, preserving local composition
   and coarse-scale zonation.
3. Expression-level within-tile null: preserves cell labels and locations, permuting gene
   expression within local tiles to test whether expression within cell type is spatially coordinated.

Why:
A global null rejects whenever cell types are non-randomly positioned, conflating large-scale
tissue compartmentalisation (tumour nests, stroma) with specific juxtacrine contact signaling.
Spatially constrained block nulls enforce that senders and receivers must be in the same local
microenvironment under the null, testing whether direct contact adds signal beyond regional co-occurrence.

Adapted from:
`src/pdac/v3_nulls.py` and `src/spatial_null.py`.
"""

from __future__ import annotations

from typing import Any
import numpy as np

from ccc_audit.graphs import assign_tiles


def global_null(
    labels: np.ndarray,
    n_perms: int = 1000,
    seed: int = 42,
) -> list[np.ndarray]:
    """Generate global label permutations across all cells.

    Source file: `src/spatial_null.py` (generate_label_permutations).
    """
    labels = np.asarray(labels)
    n_cells = labels.shape[0]
    rng = np.random.default_rng(seed)
    perms = []
    for _ in range(n_perms):
        perms.append(labels[rng.permutation(n_cells)])
    return perms


def _within_group_perms(
    labels: np.ndarray,
    groups: np.ndarray,
    n_perms: int = 1000,
    seed: int = 42,
) -> list[np.ndarray]:
    """Vectorised independent permutation within discrete group codes.

    Source file: `src/pdac/v3_nulls.py` (_within_group_perms).
    """
    labels = np.asarray(labels)
    groups = np.asarray(groups)
    n_cells = labels.shape[0]

    rng = np.random.default_rng(seed)
    pos = np.argsort(groups, kind="stable")
    perms = []

    for _ in range(n_perms):
        # lexsort on random float inside group maintains group grouping while shuffling within
        order = np.lexsort((rng.random(n_cells), groups))
        new_labels = np.empty_like(labels)
        new_labels[pos] = labels[order]
        perms.append(new_labels)

    return perms


def block_null(
    labels: np.ndarray,
    xy: np.ndarray,
    tile_size_um: float,
    n_perms: int = 1000,
    random_offset: bool = True,
    seed: int = 42,
) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Generate block label permutations within square spatial tiles.

    Source file: `src/pdac/v3_nulls.py` (block_perms).
    """
    codes, info = assign_tiles(xy, tile_size_um, random_offset=random_offset, seed=seed)
    perms = _within_group_perms(labels, codes, n_perms=n_perms, seed=seed)
    return perms, info


def expression_tile_null(
    expression_matrix: np.ndarray,
    xy: np.ndarray,
    tile_size_um: float,
    n_perms: int = 1000,
    random_offset: bool = True,
    seed: int = 42,
) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Permute expression values within spatial tiles keeping cell-type labels fixed.

    Source file: `src/pdac/v3_phase1_constrained_nulls.py`.
    """
    X = np.asarray(expression_matrix)
    codes, info = assign_tiles(xy, tile_size_um, random_offset=random_offset, seed=seed)
    n_cells = X.shape[0]
    pos = np.argsort(codes, kind="stable")
    rng = np.random.default_rng(seed)

    perm_matrices = []
    for _ in range(n_perms):
        order = np.lexsort((rng.random(n_cells), codes))
        new_X = np.empty_like(X)
        new_X[pos] = X[order]
        perm_matrices.append(new_X)

    return perm_matrices, info


# Distance bands (um) from the tumour boundary, as in the thesis (src/pdac/v3_nulls.py::DIST_BANDS).
DIST_BANDS = np.array([0.0, 25.0, 50.0, 100.0, 200.0, 400.0, np.inf])


def distance_band_null(
    labels: np.ndarray,
    xy: np.ndarray,
    tumour_label: int,
    n_perms: int = 1000,
    seed: int = 42,
    min_compartment: int = 50,
) -> tuple[list[np.ndarray] | None, dict[str, Any]]:
    """Distance-stratified label null: permute labels within (compartment x distance-band) strata.

    Port of `src/pdac/v3_nulls.py::dist_perms`. Returns (None, info) when a compartment has fewer than
    `min_compartment` cells, exactly as the thesis did.
    """
    from sklearn.neighbors import NearestNeighbors

    labels = np.asarray(labels)
    is_t = labels == tumour_label
    nt, ns = int(is_t.sum()), int((~is_t).sum())
    if nt < min_compartment or ns < min_compartment:
        return None, {"skipped": f"compartment too small (tumour {nt}, non-tumour {ns})"}
    d = np.empty(xy.shape[0])
    nn_t = NearestNeighbors(n_neighbors=1).fit(xy[is_t])
    nn_s = NearestNeighbors(n_neighbors=1).fit(xy[~is_t])
    d[~is_t] = nn_t.kneighbors(xy[~is_t])[0].ravel()
    d[is_t] = -nn_s.kneighbors(xy[is_t])[0].ravel()
    band = np.searchsorted(DIST_BANDS, np.abs(d), side="right") - 1
    codes = band + 10 * is_t.astype(np.int64)
    _, cnt = np.unique(codes, return_counts=True)
    info = {"n_tumour": nt, "n_non_tumour": ns, "n_bands_used": int(len(cnt)), "min_band_cells": int(cnt.min())}
    return _within_group_perms(labels, codes, n_perms=n_perms, seed=seed), info

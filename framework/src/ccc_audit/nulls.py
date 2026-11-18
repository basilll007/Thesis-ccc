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
    n_perms: int = 500,
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
    n_perms: int = 500,
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
    n_perms: int = 500,
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
    n_perms: int = 500,
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

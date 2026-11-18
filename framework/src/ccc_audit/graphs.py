"""Spatial graph construction and spatial tiling assignment.

What this module does:
Builds directed cell-cell contact graphs (k-nearest neighbours and fixed-radius Euclidean graphs)
from physical 2D coordinates and assigns cells to spatial grid tiles with random offsets.

Why:
Cell-cell contact hypotheses require defining physical neighbours in tissue. Using physical
micrometre distances ensures scale independence across imaging platforms. Furthermore,
spatial tile assignment with random offset prevents artificial boundary pinning artefacts
when computing spatially constrained block nulls or controlling for regional expression zonation.

Adapted from:
`src/pdac/v3_nulls.py` and `src/pdac/phase3_pdac_audit.py`.
"""

from __future__ import annotations

from typing import Any
import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors


def knn_graph(
    xy: np.ndarray,
    k: int = 6,
    include_self: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct a directed k-nearest neighbours graph from 2D coordinates.

    Parameters
    ----------
    xy : np.ndarray
        Array of shape (n_cells, 2) containing physical coordinates in micrometres.
    k : int, default=6
        Number of nearest neighbours per cell.
    include_self : bool, default=False
        Whether to include self-loops (i -> i).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        rows, cols: 1D arrays of integer indices defining directed edges (rows[e] -> cols[e]).
    """
    xy = np.asarray(xy, dtype=float)
    n_cells = xy.shape[0]
    if n_cells <= 1:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    actual_k = min(k + (1 if not include_self else 0), n_cells)
    nn = NearestNeighbors(n_neighbors=actual_k, metric="euclidean").fit(xy)
    # v2 fix: query with the points themselves so column 0 is the cell itself. v1 called kneighbors() without X,
    # which already excludes self, and then dropped column 0 -> it discarded each cell's true nearest neighbour.
    indices = nn.kneighbors(xy, return_distance=False)

    if not include_self:
        # Exclude the 0th neighbour (the cell itself)
        nbr_indices = indices[:, 1:]
    else:
        nbr_indices = indices

    actual_nbr_count = nbr_indices.shape[1]
    rows = np.repeat(np.arange(n_cells, dtype=np.int64), actual_nbr_count)
    cols = nbr_indices.ravel().astype(np.int64)

    return rows, cols


def radius_graph(
    xy: np.ndarray,
    radius_um: float,
    include_self: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct a directed radius graph where edges connect cells within radius_um.

    Parameters
    ----------
    xy : np.ndarray
        Array of shape (n_cells, 2) containing physical coordinates in micrometres.
    radius_um : float
        Maximum distance in micrometres for an edge.
    include_self : bool, default=False
        Whether to include self-loops.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        rows, cols: 1D arrays of integer indices defining directed edges.
    """
    xy = np.asarray(xy, dtype=float)
    nn = NearestNeighbors(radius=radius_um, metric="euclidean").fit(xy)
    adj = nn.radius_neighbors_graph(mode="connectivity")

    if not include_self:
        adj.setdiag(0)
        adj.eliminate_zeros()

    adj_coo = adj.tocoo()
    return adj_coo.row.astype(np.int64), adj_coo.col.astype(np.int64)


def assign_tiles(
    xy: np.ndarray,
    tile_size_um: float,
    random_offset: bool = True,
    seed: int | None = 42,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Assign cells to regular 2D square spatial tiles with optional random offset.

    Parameters
    ----------
    xy : np.ndarray
        Array of shape (n_cells, 2) containing coordinates in micrometres.
    tile_size_um : float
        Side length of square tiles in micrometres (e.g. 50, 100, 200).
    random_offset : bool, default=True
        Whether to apply a uniform random offset in [0, tile_size_um) before tiling.
    seed : int, optional
        Random seed for the coordinate offset.

    Returns
    -------
    tuple[np.ndarray, dict[str, Any]]
        codes : 1D array of unique integer tile identifiers for each cell.
        info : summary dictionary containing offset, tile count, cells per tile.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.shape[0] == 0:
        return np.array([], dtype=np.int64), {}

    min_xy = xy.min(axis=0)
    offset = np.zeros(2, dtype=float)
    if random_offset:
        rng = np.random.default_rng(seed)
        offset = rng.uniform(0.0, tile_size_um, size=2)

    shifted = xy - min_xy + offset
    grid = np.floor(shifted / tile_size_um).astype(np.int64)
    max_y = grid[:, 1].max() + 1
    codes = grid[:, 0] * max_y + grid[:, 1]

    # Remap codes to compact 0..K-1
    _, compact_codes = np.unique(codes, return_inverse=True)
    _, counts = np.unique(compact_codes, return_counts=True)

    info = {
        "tile_size_um": float(tile_size_um),
        "random_offset": random_offset,
        "offset_um": offset.tolist(),
        "n_tiles": int(len(counts)),
        "median_cells_per_tile": float(np.median(counts)),
        "min_cells_per_tile": int(counts.min()),
        "max_cells_per_tile": int(counts.max()),
        "n_singleton_tiles": int((counts == 1).sum()),
    }
    return compact_codes, info

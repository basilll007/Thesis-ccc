"""Loaders for candidate PDAC spatial datasets.

Each loader returns (AnnData with raw counts in .X and centroids in .obsm['spatial'] in micrometres, metadata dict).
Downloaded files live under data/pdac/<dataset>/ (git-ignored) and are treated as untrusted data:
they are only parsed as tables/HDF5, never executed.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _filemeta(paths):
    return [{"path": Path(p).as_posix().split("data/", 1)[-1], "bytes": Path(p).stat().st_size, "sha256": sha256(Path(p))}
            for p in paths]


def load_xenium_dir(d: Path, log):
    """10x Xenium output bundle: cell_feature_matrix.h5 + cells.csv.gz / cells.parquet."""
    h5 = next(iter(sorted(d.rglob("cell_feature_matrix.h5"))))
    adata = sc.read_10x_h5(h5, gex_only=True)
    adata.var_names_make_unique()
    cells_p = next(iter(sorted(list(d.rglob("cells.parquet")) + list(d.rglob("cells.csv.gz")))))
    cells = pd.read_parquet(cells_p) if cells_p.suffix == ".parquet" else pd.read_csv(cells_p)
    cells["cell_id"] = cells["cell_id"].astype(str).str.replace("^b'|'$", "", regex=True)
    cells = cells.set_index("cell_id").reindex(adata.obs_names.astype(str))
    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    for c in ("transcript_counts", "cell_area", "nucleus_area"):
        if c in cells:
            adata.obs[c] = cells[c].to_numpy()
    log.info("Xenium loaded: %d cells x %d genes from %s", adata.n_obs, adata.n_vars, h5.parent.name)
    return adata, _filemeta([h5, cells_p])


DATASETS = {}


def load(key: str, data_dir: Path, log):
    if key not in DATASETS:
        raise KeyError(f"Unknown dataset key {key}; known: {list(DATASETS)}")
    return DATASETS[key](data_dir / "pdac" / key, log)

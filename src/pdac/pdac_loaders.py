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


COSMX_UM_PER_PX = 0.18  # Shiau et al. bioRxiv 2023 (doi 10.1101/2023.06.28.546848) Methods: "60 pixels (10.8 um)"


def load_mendeley_pdac_smi(d: Path, log, sample: str | None = None):
    """Mendeley PDAC-SMI (CosMx). One section (`sample_id`) per call; Base probes only.

    Spike-in probes (NaN where not measured) and negative probes are dropped. Coordinates are converted
    from global pixels to micrometres with COSMX_UM_PER_PX.
    """
    f = d / "raw_meta_data_final.h5ad"
    full = ad.read_h5ad(f, backed="r")
    base = (full.var["probe_type"].astype(str) == "Base").to_numpy()
    mask = np.ones(full.n_obs, bool) if sample is None else (full.obs["sample_id"].astype(str) == sample).to_numpy()
    if mask.sum() == 0:
        raise ValueError(f"sample {sample} not found")
    sub = full[mask, base].to_memory()
    X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    if np.isnan(X.data).any():
        raise ValueError("NaN found in Base probes; unexpected")
    adata = ad.AnnData(X=X.astype(np.float32), obs=sub.obs.copy(), var=sub.var.copy())
    adata.obsm["spatial"] = adata.obs[["CenterX_global_px", "CenterY_global_px"]].to_numpy(dtype=float) * COSMX_UM_PER_PX
    adata.obs["author_majortype"] = adata.obs["annotation_majortypes"].astype(str)
    log.info("Mendeley PDAC-SMI section %s: %d cells x %d Base genes (dropped %d spike-in/negative probes); "
             "um/px=%.2f", sample, adata.n_obs, adata.n_vars, int((~base).sum()), COSMX_UM_PER_PX)
    meta = {"dataset": "Mendeley PDAC-SMI (doi 10.17632/kx6b69n3cb.1, CC BY 4.0)", "section": sample,
            "files": _filemeta([f]), "n_cells_section": int(adata.n_obs), "n_base_genes": int(adata.n_vars),
            "dropped_features": full.var_names[~base].tolist(), "um_per_px": COSMX_UM_PER_PX,
            "treatment_status": str(adata.obs["treatment_status"].iloc[0]), "donor": str(adata.obs["donor"].iloc[0]),
            "n_fov": int(adata.obs["fov"].nunique())}
    return adata, meta


DATASETS = {"mendeley_pdac_smi": load_mendeley_pdac_smi}


def load(key: str, data_dir: Path, log):
    """key = '<dataset>' or '<dataset>:<section>' (e.g. 'mendeley_pdac_smi:U5-a')."""
    name, _, section = key.partition(":")
    if name not in DATASETS:
        raise KeyError(f"Unknown dataset key {name}; known: {list(DATASETS)}")
    return DATASETS[name](data_dir / "pdac" / name, log, section or None)

"""Loaders that turn the thesis' original input files into package-ready AnnData objects.

What this module does:
Re-implements, inside the package, the input preparation the thesis used for the Mendeley PDAC-SMI
CosMx data (`src/pdac/pdac_loaders.py::load_mendeley_pdac_smi`, `src/pdac/phase3_pdac_audit.py::qc`,
`src/pdac/v3_nulls.py::load_section`) and for GEO GSE310352 slides:
  1. Read the raw author h5ad (backed), keep Base probes, keep one section (`sample_id`).
  2. Coordinates: CenterX/Y_global_px * 0.18 um/px.
  3. QC: cells with >= 10 counts, genes detected in >= 3 cells.
  4. counts -> layers['counts']; X = log1p(normalize_total(1e4)).
  5. Attach the thesis cell annotations from `results/pdac/phaseB/sections/<section>_*.csv.gz`.

Why:
Reproduction tests must recompute thesis numbers from the original inputs, not reload result files.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

COSMX_UM_PER_PX = 0.18  # thesis constant (src/pdac/pdac_loaders.py)
MENDELEY_SECTIONS = ["T1-b", "T2-b", "T3-b", "T4-a", "T4-b", "T5-a", "T6-b", "U1-b", "U2-b", "U3-b",
                     "U4-a", "U4-b", "U5-a", "U6-a_1", "U6-a_2", "U7-a"]


def _normalise_log1p(counts: sp.csr_matrix, target: float = 1e4) -> sp.csr_matrix:
    tot = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.where(tot > 0, target / np.maximum(tot, 1e-12), 0.0)
    X = sp.diags(scale) @ counts
    X = X.tocsr()
    X.data = np.log1p(X.data)
    return X


def load_mendeley_section(
    h5ad_path: str | Path,
    section: str,
    labels_dir: str | Path | None = None,
    annotations: tuple[str, ...] = ("consensus_T",),
    tumour_map: dict[str, str] | None = None,
) -> ad.AnnData:
    """Load one Mendeley PDAC-SMI section exactly as the thesis v3 pipeline did.

    Returns AnnData with X = log1p CP10K, layers['counts'] = raw counts, obs columns x_um, y_um,
    section_id, patient_id (tumour/donor), and one column per requested annotation.
    """
    full = ad.read_h5ad(h5ad_path, backed="r")
    base = (full.var["probe_type"].astype(str) == "Base").to_numpy()
    mask = (full.obs["sample_id"].astype(str) == section).to_numpy()
    if mask.sum() == 0:
        raise ValueError(f"section {section} not found")
    sub = full[mask, base].to_memory()
    full.file.close()
    X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    X = X.astype(np.float32)
    obs = sub.obs.copy()
    var = sub.var.copy()
    # QC identical to thesis: filter_cells(min_counts=10) then filter_genes(min_cells=3)
    tot = np.asarray(X.sum(axis=1)).ravel()
    keep_c = tot >= 10
    X = X[keep_c]
    obs = obs.loc[keep_c]
    ncell = np.asarray((X > 0).sum(axis=0)).ravel()
    keep_g = ncell >= 3
    X = X[:, keep_g]
    var = var.loc[keep_g]
    adata = ad.AnnData(X=_normalise_log1p(X), obs=obs, var=var)
    adata.layers["counts"] = X
    xy = obs[["CenterX_global_px", "CenterY_global_px"]].to_numpy(dtype=float) * COSMX_UM_PER_PX
    adata.obs["x_um"] = xy[:, 0]
    adata.obs["y_um"] = xy[:, 1]
    adata.obs["section_id"] = section
    donor = str(obs["donor"].iloc[0]) if "donor" in obs else section
    adata.obs["patient_id"] = (tumour_map or {}).get(section, donor)
    if labels_dir is not None:
        labels_dir = Path(labels_dir)
        ids = adata.obs_names.astype(str).to_numpy()
        for a in annotations:
            if a == "consensus_T":
                df = pd.read_csv(labels_dir / f"{section}_consensus_labels.csv.gz", usecols=["cell_id", a])
            else:
                df = pd.read_csv(labels_dir / f"{section}_labels.csv.gz", usecols=["cell_id", a])
            if not (df.cell_id.astype(str).to_numpy() == ids).all():
                raise ValueError(f"cell order mismatch for {section} / {a}")
            adata.obs[a] = df[a].astype(str).to_numpy()
    return adata


def mendeley_section_to_tumour(section_info_csv: str | Path) -> dict[str, str]:
    """Thesis rule (src/pdac/v3_nulls.py::section_to_tumour): tumour = donor without _2nd/_B suffix."""
    info = pd.read_csv(section_info_csv)
    info["tumour"] = info.donor.astype(str).str.replace(r"(_2nd|_B)$", "", regex=True)
    return dict(zip(info.section, info.tumour))

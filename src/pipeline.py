"""Core AnnData, QC, clustering, spatial graph, and baseline CCC steps."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from config import SEED

LR_PAIRS = [
    ("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"), ("EGF", "EGFR"),
    ("FGF2", "FGFR1"), ("VEGFA", "KDR"), ("CXCL12", "CXCR4"),
    ("CCL2", "CCR2"), ("CCL5", "CCR5"), ("CXCL9", "CXCR3"),
    ("CXCL10", "CXCR3"), ("TNF", "TNFRSF1A"), ("IFNG", "IFNGR1"),
    ("IL6", "IL6R"), ("IL1B", "IL1R1"), ("HGF", "MET"),
    ("PDGFB", "PDGFRA"), ("PDGFB", "PDGFRB"), ("JAG1", "NOTCH1"),
    ("CD274", "PDCD1"), ("MIF", "CD74"),
]
MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "immune": ["PTPRC", "CD74", "HLA-DRA", "LYZ"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
}


def locate_h5(raw_dir: Path) -> Path:
    candidates = list(raw_dir.rglob("cell_feature_matrix.h5")) + list(raw_dir.rglob("*.h5"))
    if not candidates:
        raise FileNotFoundError("No Xenium/10x HDF5 matrix found under data/raw")
    return candidates[0]


def locate_cells(raw_dir: Path) -> Path | None:
    candidates = list(raw_dir.rglob("cells.csv.gz")) + list(raw_dir.rglob("cells.csv"))
    return candidates[0] if candidates else None


def load_data(raw_dir: Path) -> ad.AnnData:
    adata = sc.read_10x_h5(locate_h5(raw_dir))
    adata.var_names_make_unique()
    cells_path = locate_cells(raw_dir)
    if cells_path:
        cells = pd.read_csv(cells_path)
        cells.index = cells.get("cell_id", cells.index).astype(str)
        cells = cells.reindex(adata.obs_names.astype(str))
        for column in ("x_centroid", "y_centroid", "x", "y"):
            if column in cells:
                adata.obs[column] = pd.to_numeric(cells[column], errors="coerce").to_numpy()
        if {"x_centroid", "y_centroid"}.issubset(cells.columns):
            adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()
    if "spatial" not in adata.obsm:
        raise ValueError("Xenium cells.csv did not provide x_centroid/y_centroid coordinates")
    return adata


def quality_control(adata: ad.AnnData, log: list[str]) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["n_counts"] = np.asarray(adata.X.sum(axis=0)).ravel()
    before = (adata.n_obs, adata.n_vars)
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=3)
    log.append(f"QC counts: {before[0]} cells/{before[1]} genes -> {adata.n_obs} cells/{adata.n_vars} genes")
    return adata


def annotate_and_cluster(adata: ad.AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars), flavor="cell_ranger", subset=False)
    sc.pp.pca(adata, n_comps=min(30, adata.n_obs - 1, adata.n_vars - 1), random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)), random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED, key_added="leiden")
    scores = {}
    for label, markers in MARKERS.items():
        present = [gene for gene in markers if gene in adata.var_names]
        scores[label] = np.asarray(adata[:, present].X.mean(axis=1)).ravel() if present else np.zeros(adata.n_obs)
    score_frame = pd.DataFrame(scores, index=adata.obs_names)
    adata.obs["cell_type"] = score_frame.idxmax(axis=1).astype("category")
    adata.uns["cell_type_marker_genes"] = MARKERS


def add_spatial_graph(adata: ad.AnnData) -> None:
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6, set_diag=False)
    adata.uns["spatial_graph_parameters"] = {"method": "squidpy.gr.spatial_neighbors", "coord_type": "generic", "n_neighs": 6, "set_diag": False}


def baseline_lr(adata: ad.AnnData) -> pd.DataFrame:
    present = [(ligand, receptor) for ligand, receptor in LR_PAIRS if ligand in adata.var_names and receptor in adata.var_names]
    if not present:
        return pd.DataFrame(columns=["sender", "receiver", "ligand", "receptor", "score", "pvalue"])
    # One permutation is sufficient for v0 smoke validation; calibrated inference is future work.
    result = sq.gr.ligrec(
        adata,
        cluster_key="cell_type",
        interactions=pd.DataFrame(present, columns=["source", "target"]),
        n_perms=1,
        seed=SEED,
        copy=True,
        key_added="baseline_ligrec",
    )
    means = result["means"]
    pvalues = result["pvalues"]
    rows = []
    for (ligand, receptor), mean_row in means.iterrows():
        for sender, receiver in means.columns:
            score = float(mean_row[(sender, receiver)])
            pvalue = float(pvalues.loc[(ligand, receptor), (sender, receiver)])
            if score > 0:
                rows.append({"sender": sender, "receiver": receiver, "ligand": ligand, "receptor": receptor, "score": score, "pvalue": pvalue})
    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)

"""Shared helpers for the v2 analyses (spillover control, cross-type share on aggregated baselines).

What: memory gate (log free RAM, refuse to start below 8 GB), v2 cache paths, per-cell metadata from the
Mendeley h5ad (cell area and protein channels), nearest-CAF distances, and column-group definitions that split the
aggregated neighbourhood feature blocks by neighbour cell type.

Why: v2 phases reuse the v1 datasets, folds, feature cache and leakage guard unchanged; everything new that more
than one v2 script needs lives here so the definitions cannot drift between Phase 1 and Phase 2.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

V2_CACHE = C.CACHE / "v2"
MIN_FREE_GB = 8.0
CAF = "CAF_fibroblast"
MYELOID = "myeloid"
TUMOUR = "epithelial_ductal_tumor"
UM_PER_PX = 0.18          # src/pdac/pdac_loaders.COSMX_UM_PER_PX (Mendeley CosMx)
TYPED_PREFIX = {"tumour": TUMOUR, "CAF": CAF, "T": "T_cell", "myeloid": MYELOID}


def free_ram_gb() -> float:
    import psutil
    return psutil.virtual_memory().available / 2**30


def memory_gate(log, wait_s: int = 60, max_wait_s: int = 3600) -> float:
    """Log free RAM; if below MIN_FREE_GB, log and wait (up to max_wait_s) instead of starting."""
    waited = 0
    while True:
        gb = free_ram_gb()
        log.info("free RAM %.2f GB (minimum to start %.1f GB)", gb, MIN_FREE_GB)
        if gb >= MIN_FREE_GB:
            return gb
        if waited >= max_wait_s:
            raise RuntimeError(f"free RAM {gb:.2f} GB stayed below {MIN_FREE_GB} GB for {max_wait_s} s; not starting")
        log.warning("free RAM below %.1f GB: waiting %d s", MIN_FREE_GB, wait_s)
        time.sleep(wait_s)
        waited += wait_s


def cell_metadata(log=None) -> pd.DataFrame:
    """Per-cell area (um^2), equivalent diameter (um) and protein channel means for every development cell in the
    Phase 0 datasets, matched by (section, cell_id). Cached in model_cache/v2/cell_meta.csv.gz (checksummed)."""
    path = V2_CACHE / "cell_meta.csv.gz"
    if path.exists():
        return pd.read_csv(path, dtype={"section": str, "group": str, "cell_id": str})
    import anndata as ad
    a = ad.read_h5ad(C.DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    cols = ["sample_id", "Area", "Mean.CD45", "Mean.PanCK", "Mean.CD3", "Mean.MembraneStain", "Mean.DAPI"]
    obs = a.obs[cols].copy()
    obs["cell_id"] = obs.index.astype(str)
    obs["sample_id"] = obs["sample_id"].astype(str)
    a.file.close()
    parts = []
    for f in sorted((C.CACHE / "datasets" / "mendeley_pdac_smi").glob("*.npz")):
        z = np.load(f, allow_pickle=False)
        sec = str(z["section"])
        d = pd.DataFrame({"section": sec, "row": np.arange(len(z["cell_id"])), "cell_id": z["cell_id"].astype(str),
                          "label": z["label"].astype(str), "group": str(z["group"])})
        m = d.merge(obs[obs.sample_id == sec].drop(columns="sample_id"), on="cell_id", how="left", validate="one_to_one")
        if m["Area"].isna().any():
            raise ValueError(f"{sec}: {int(m['Area'].isna().sum())} cells without h5ad metadata")
        parts.append(m)
    out = pd.concat(parts, ignore_index=True)
    out["area_um2"] = out["Area"] * UM_PER_PX ** 2
    out["eq_diam_um"] = 2 * np.sqrt(out["area_um2"] / np.pi)
    V2_CACHE.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    register_v2_artefact(path, "v2_phase1", "per-cell area and protein channel means (Mendeley h5ad obs) matched to Phase 0 rows")
    if log:
        log.info("cell metadata: %d cells, saved %s", len(out), C.rel(path))
    return out


def nearest_type_distance(xy: np.ndarray, label: np.ndarray, cls: str, k: int = 1):
    """Distance (um) and index of the k nearest cells of class `cls` for every cell (excluding the cell itself)."""
    from scipy.spatial import cKDTree
    idx_c = np.where(label == cls)[0]
    tree = cKDTree(xy[idx_c])
    kk = k + 1
    dist, j = tree.query(xy, k=kk)
    j = idx_c[j]
    self_hit = j == np.arange(len(xy))[:, None]
    # drop the self match (only possible for cells of class cls) and keep the first k remaining
    dist = np.where(self_hit, np.inf, dist)
    order = np.argsort(dist, axis=1, kind="stable")[:, :k]
    return np.take_along_axis(dist, order, 1), np.take_along_axis(j, order, 1)


def column_type(name: str) -> str:
    """Neighbour cell type a feature column describes: 'tumour', 'CAF', 'T', 'myeloid', 'other', or 'none'
    (own covariates, tile means, distances). Used for block ablations by neighbour type."""
    if name.startswith(("knn_frac__", "r30_count__")):
        cls = name.split("__", 1)[1]
        for t, c in TYPED_PREFIX.items():
            if cls == c:
                return t
        return "other"
    if name.startswith(("nb_knn_", "nb_r30_")):
        t = name.split("_", 3)[2].split("__")[0]
        return t if t in TYPED_PREFIX or t == "other" else "other"
    if name.startswith("ncaf_"):
        return "CAF"
    return "none"


def column_kind(name: str) -> str:
    """'composition' (who is nearby: fractions, counts, distances) or 'expression' (what they express)."""
    if name.startswith(("knn_frac__", "r30_count__", "ncaf_dist")):
        return "composition"
    if name.startswith(("nb_knn_", "nb_r30_", "ncaf_expr__")):
        return "expression"
    return "own"


def register_v2_artefact(path: Path, phase: str, description: str) -> None:
    """Record a large v2 artefact (outside git) with its checksum in model_dev/results/v2_artefacts.csv.
    v2 keeps its own registry so the v1 registry (results/phase0/artefacts.csv) is never modified."""
    from datetime import datetime
    reg = C.RESULTS / "v2_artefacts.csv"
    row = {"path": C.rel(path), "phase": phase, "description": description, "bytes": Path(path).stat().st_size,
           "sha256": C.sha256(path), "registered_at": datetime.now().isoformat(timespec="seconds")}
    df = pd.read_csv(reg) if reg.exists() else pd.DataFrame(columns=list(row))
    df = df[df.path != row["path"]]
    pd.concat([df, pd.DataFrame([row])], ignore_index=True).to_csv(reg, index=False)

"""Phase 4a (env thesis-spatial-ccc): our C/E/M/Q scores and Squidpy gr.ligrec on the benchmark inputs.

Benchmark inputs (pre-declared before any benchmark result was seen):
  * crop : the CROP_N cells nearest to the cell closest to the section's median coordinate (deterministic,
           contiguous). Needed because COMMOT builds a dense n x n distance matrix (infeasible for 74k cells).
  * full : the whole primary section (Squidpy ligrec and our scores; SpatialDM runs in phase4b).
Labels: cell_type_marker rebuilt here exactly as in Phase 3 (load -> QC -> normalise/log1p -> per-cell argmax of
mean marker expression); marker scoring does not depend on PCA/Leiden, so the labels are identical to Phase 3.
(The Phase 3 processed AnnData was not written because the host stopped that run for low memory.)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, SEED, finish_run, save_csv, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import squidpy as sq  # noqa: E402

from lr_catalog import catalogue  # noqa: E402
import pdac_loaders  # noqa: E402
import scanpy as sc  # noqa: E402
from phase3_pdac_audit import PDAC_MARKERS, knn_edges, qc, run_scores  # noqa: E402

PHASE = "phase4"
TAG = "mendeley_pdac_smi__U5-a"
CROP_N = 8000
N_PERMS = 500


def crop_ids(adata) -> np.ndarray:
    xy = np.asarray(adata.obsm["spatial"], float)
    med = np.median(xy, axis=0)
    c = np.argmin(np.linalg.norm(xy - med, axis=1))
    d = np.linalg.norm(xy - xy[c], axis=1)
    return np.sort(np.argsort(d, kind="stable")[:CROP_N])


def ligrec_table(adata, pairs, log, tag):
    t0 = time.time()
    res = sq.gr.ligrec(adata, cluster_key="cell_type_marker", interactions=pd.DataFrame(pairs, columns=["source", "target"]),
                       n_perms=N_PERMS, seed=SEED, copy=True, use_raw=True, threshold=0.01, n_jobs=4, show_progress_bar=False)
    means, pvals = res["means"], res["pvalues"]
    rows = []
    for (l, r), mrow in means.iterrows():
        for (s, rc) in means.columns:
            rows.append({"ligand": l, "receptor": r, "sender": s, "receiver": rc,
                         "ligrec_mean": float(mrow[(s, rc)]), "ligrec_pvalue": float(pvals.loc[(l, r), (s, rc)])})
    df = pd.DataFrame(rows)
    log.info("[%s] squidpy ligrec: %d rows, %d LR pairs kept by ligrec, %.1fs", tag, len(df), means.shape[0], time.time() - t0)
    return df


def main() -> None:
    log = start_run("phase4a_ours_and_squidpy")
    full, _ = pdac_loaders.load("mendeley_pdac_smi:U5-a", DATA_DIR, log)
    full, _, _ = qc(full, log)
    full.layers["counts"] = full.X.copy()
    sc.pp.normalize_total(full, target_sum=1e4)
    sc.pp.log1p(full)
    full.raw = full
    sf = pd.DataFrame({t_: np.asarray(full[:, [g for g in gs if g in full.var_names]].X.mean(axis=1)).ravel()
                       for t_, gs in PDAC_MARKERS.items() if sum(g in full.var_names for g in gs) >= 2}, index=full.obs_names)
    full.obs["cell_type_marker"] = sf.idxmax(axis=1).astype("category")
    log.info("marker label counts: %s", full.obs["cell_type_marker"].value_counts().to_dict())
    full.write_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_processed.h5ad")
    vn = list(full.var_names)
    pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
    idx = crop_ids(full)
    crop = full[idx].copy()
    crop_df = pd.DataFrame({"cell_id": crop.obs_names, "x_um": crop.obsm["spatial"][:, 0], "y_um": crop.obsm["spatial"][:, 1],
                            "cell_type_marker": crop.obs["cell_type_marker"].astype(str).to_numpy()})
    save_csv(crop_df, PHASE, "benchmark_crop_cells.csv", f"The {CROP_N} cells of the deterministic benchmark crop (U5-a)", log)
    ext = np.ptp(crop.obsm["spatial"], axis=0)
    save_json({"crop_n": CROP_N, "rule": "CROP_N nearest cells to the cell closest to the median coordinate",
               "extent_um": ext.tolist(), "cell_type_counts": crop_df["cell_type_marker"].value_counts().to_dict(),
               "n_pairs": len(pairs)}, PHASE, "benchmark_crop_info.json", "Benchmark crop definition and composition", log)
    log.info("crop: %d cells, extent %s um", crop.n_obs, ext)

    for tag, a in (("crop", crop), ("full", full)):
        labels = np.asarray(a.obs["cell_type_marker"].astype(str).to_numpy())
        rows, cols = knn_edges(np.asarray(a.obsm["spatial"], float), 6)
        ours, _ = run_scores(a, labels, rows, cols, pairs, log, f"ours/{tag}")
        save_csv(ours, PHASE, f"ours_CEMQ_6NN_{tag}.csv", f"Our C/E/M/Q scores, 6-NN, 500 perms, marker labels ({tag})", log)
        lr = ligrec_table(a, pairs, log, tag)
        save_csv(lr, PHASE, f"squidpy_ligrec_{tag}.csv",
                 f"squidpy.gr.ligrec means/p-values (500 perms, seed 42, threshold 0.01) ({tag})", log)
    crop.write_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_benchmark_crop.h5ad")
    log.info("crop AnnData written to data/processed/pdac_%s_benchmark_crop.h5ad (git-ignored) for phase4b", TAG)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

"""Phase 0: reproducibility check of numbers already reported for the breast run.

Re-derives, WITHOUT touching any existing output file:
  (1) PDCD1 / CD274 detection counts from data/processed/xenium_breast_baseline.h5ad and
      compares them with results/cd274_pdcd1_detection_diagnostic.csv;
  (2) the 15 um contact-restricted (radius-graph) multi-axis results, including
      CD86->CTLA4 (reported z=3.56, p=0.035), using the identical algorithm and RNG stream as
      src/run_multi_axis_expansion.py (RandomState(42), 200 permutations), and compares them with
      results/multi_axis_expansion/multi_axis_summary.csv.
Outputs go to results/pdac/phase0/.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO_ROOT, SEED, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

PHASE = "phase0"
H5AD = REPO_ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad"
AXES = [  # identical to src/run_multi_axis_expansion.py
    ("CXCL12 -> CXCR4", "CXCL12", "CXCR4", "fibroblast", "immune"),
    ("PTN -> SDC4", "PTN", "SDC4", "fibroblast", "epithelial"),
    ("CD86 -> CTLA4", "CD86", "CTLA4", "immune", "immune"),
    ("CD274 -> PDCD1", "CD274", "PDCD1", "epithelial", "immune"),
]


def main() -> None:
    log = start_run("phase0_repro_check")
    adata = ad.read_h5ad(H5AD)
    raw = adata.raw if adata.raw is not None else adata
    names = list(raw.var_names)

    def gene(g):
        col = raw.X[:, names.index(g)]
        return np.asarray(col.toarray() if hasattr(col, "toarray") else col, dtype=np.float32).ravel()

    # (1) detection counts
    diag = pd.read_csv(REPO_ROOT / "results" / "cd274_pdcd1_detection_diagnostic.csv")
    det_rows = []
    for g in ("CD274", "PDCD1"):
        x = gene(g)
        n_det = int((x > 0).sum())
        rep = int(diag.query("gene == @g and cell_type == 'ALL'")["n_detected"].iloc[0])
        det_rows.append({"gene": g, "n_cells": adata.n_obs, "n_detected_recomputed": n_det,
                         "n_detected_in_existing_csv": rep, "match": n_det == rep})
        log.info("%s detected in %d/%d cells (existing CSV: %d) match=%s", g, n_det, adata.n_obs, rep, n_det == rep)
    det_df = pd.DataFrame(det_rows)
    save_csv(det_df, PHASE, "pdcd1_cd274_detection_recheck.csv",
             "Recomputed CD274/PDCD1 detection counts vs results/cd274_pdcd1_detection_diagnostic.csv", log)

    # (2) 15 um radius-graph multi-axis recheck (same algorithm/RNG as run_multi_axis_expansion.py)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    cell_types = adata.obs["cell_type"].astype(str).values
    unique_types, types_int = np.unique(cell_types, return_inverse=True)
    type_to_id = {t: int(np.where(unique_types == t)[0][0]) for t in unique_types}
    nn = NearestNeighbors(metric="euclidean", n_jobs=-1).fit(coords)
    n_perms = 200
    rng = np.random.RandomState(SEED)
    perm_list = [rng.permutation(types_int) for _ in range(n_perms)]
    log.info("RNG: numpy.random.RandomState(%d), %d label permutations (legacy multi-axis stream)", SEED, n_perms)
    r = 15.0
    adj = nn.radius_neighbors(coords, radius=r, return_distance=False)
    rows_l, cols_l = [], []
    for i, neighs in enumerate(adj):
        m = neighs != i
        rows_l.extend([i] * int(np.sum(m)))
        cols_l.extend(neighs[m])
    rows = np.asarray(rows_l, dtype=np.int64)
    cols = np.asarray(cols_l, dtype=np.int64)
    log.info("15 um radius graph: %d directed edges", len(rows))

    prior = pd.read_csv(REPO_ROOT / "results" / "multi_axis_expansion" / "multi_axis_summary.csv").query("radius_um == 15.0")
    out_rows, null_store = [], {}
    for name, lig, rec, s, rcv in AXES:
        prod = gene(lig)[rows] * gene(rec)[cols]
        s_id, r_id = type_to_id[s], type_to_id[rcv]
        mask = (types_int[rows] == s_id) & (types_int[cols] == r_id)
        obs = float(prod[mask].mean()) if mask.any() else 0.0
        null = np.zeros(n_perms)
        for k, p in enumerate(perm_list):
            pm = (p[rows] == s_id) & (p[cols] == r_id)
            if pm.any():
                null[k] = prod[pm].mean()
        sd = float(np.std(null, ddof=1))
        z = (obs - null.mean()) / sd if sd > 1e-9 else 0.0
        p_raw = float(np.mean(null >= obs))
        p_plus1 = float((np.sum(null >= obs) + 1) / (n_perms + 1))
        pr = prior[prior["axis_name"] == name].iloc[0]
        match = bool(round(z, 2) == float(pr["z_score"]) and p_raw == float(pr["p_value"]))
        out_rows.append({"axis": name, "sender": s, "receiver": rcv, "radius_um": r, "qualifying_edges": int(mask.sum()),
                         "obs_score": obs, "null_mean": float(null.mean()), "null_std": sd, "z_score": z,
                         "p_value_legacy_mean_ge": p_raw, "p_value_plus1": p_plus1,
                         "reported_z": float(pr["z_score"]), "reported_p": float(pr["p_value"]), "exact_match": match})
        null_store[name] = null
        log.info("%s 15um: z=%.2f p_legacy=%.3f (reported z=%.2f p=%.3f) match=%s; p_plus1=%.4f",
                 name, z, p_raw, pr["z_score"], pr["p_value"], match, p_plus1)
    rc = pd.DataFrame(out_rows)
    save_csv(rc, PHASE, "multi_axis_15um_recheck.csv",
             "Recomputed 15 um contact-restricted multi-axis z/p vs results/multi_axis_expansion/multi_axis_summary.csv", log)

    fig, ax = plt.subplots(1, len(AXES), figsize=(16, 3.4))
    for a, (name, *_r) in zip(ax, AXES):
        a.hist(null_store[name], bins=30, color="#9ecae1")
        a.axvline(rc.loc[rc.axis == name, "obs_score"].iloc[0], color="#de2d26", lw=2)
        zz = rc.loc[rc.axis == name, "z_score"].iloc[0]
        a.set_title(f"{name} (15 µm)\nz={zz:.2f}", fontsize=10)
        a.set_xlabel("edge score under label permutation")
    fig.suptitle("Phase 0 recheck: 15 µm contact-restricted null (200 perms, RandomState(42)); red = observed", fontsize=11)
    save_figure(fig, PHASE, "phase0_multi_axis_15um_recheck_nulls",
                "Null histograms for 15 um contact-restricted recheck of the four breast axes", log)

    summary = {"pdcd1_detected": int(det_df.loc[det_df.gene == "PDCD1", "n_detected_recomputed"].iloc[0]),
               "all_detection_counts_match": bool(det_df["match"].all()),
               "multi_axis_15um_all_exact_match": bool(rc["exact_match"].all()),
               "cd86_ctla4_15um_z": float(rc.loc[rc.axis == "CD86 -> CTLA4", "z_score"].iloc[0]),
               "cd86_ctla4_15um_p_legacy": float(rc.loc[rc.axis == "CD86 -> CTLA4", "p_value_legacy_mean_ge"].iloc[0])}
    save_json(summary, PHASE, "phase0_repro_summary.json", "Phase 0 reproducibility-check summary", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

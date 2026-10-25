"""Phase 4b (env pdac-bench: commot 0.0.3, SpatialDM 0.3.1): COMMOT and SpatialDM on the benchmark inputs.

COMMOT  : crop only (dense n x n distance matrix). L-R list = our testable catalogue pairs (non-heteromeric),
          dis_thr = 50 um (pre-declared; > typical 6-NN distance, < the 75 um radius-sweep step), other params default.
          Cluster-level summary via commot.tl.cluster_communication (marker labels, 100 permutations = package
          default, random_seed = 42).
SpatialDM: crop and full section. weight_matrix(l=25 um, cutoff=0.2, single_cell=True) (pre-declared: weights fall
          below the cutoff at ~45 um, comparable to COMMOT's 50 um); extract_lr(species='human', datahost='package',
          min_cell=3) = built-in CellChatDB; spatialdm_global(method='z-score'). Output is per L-R pair (global),
          not per cell-type pair.
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

PHASE = "phase4"
TAG = "mendeley_pdac_smi__U5-a"
DIS_THR = 50.0
SDM_L = 25.0


def catalogue_pairs(vn):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from lr_catalog import catalogue
    return [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]


def run_commot(adata, pairs, log):
    import commot as ct
    a = adata.copy()
    a.X = a.raw.X.copy() if a.raw is not None else a.X  # log-normalised expression
    df_lr = pd.DataFrame([[l, r, f"{l}-{r}"] for l, r in pairs], columns=["ligand", "receptor", "pathway"])
    t0 = time.time()
    ct.tl.spatial_communication(a, database_name="catalog", df_ligrec=df_lr, dis_thr=DIS_THR, heteromeric=False,
                                pathway_sum=False)
    log.info("COMMOT spatial_communication: %d pairs, dis_thr=%.0f um, %.1fs", len(pairs), DIS_THR, time.time() - t0)
    rows, glob = [], []
    for l, r in pairs:
        key = f"commot-catalog-{l}-{r}"
        if key not in a.obsp:
            log.warning("COMMOT produced no matrix for %s-%s", l, r)
            continue
        S = a.obsp[key]
        glob.append({"ligand": l, "receptor": r, "commot_total_signal": float(S.sum()), "commot_nnz": int(S.nnz)})
        ct.tl.cluster_communication(a, database_name="catalog", lr_pair=(l, r), clustering="cell_type_marker",
                                    n_permutations=100, random_seed=SEED)
        res = a.uns[f"commot_cluster-cell_type_marker-catalog-{l}-{r}"]
        cm, cp = res["communication_matrix"], res["communication_pvalue"]
        for s in cm.index:
            for rc in cm.columns:
                rows.append({"ligand": l, "receptor": r, "sender": s, "receiver": rc,
                             "commot_cluster_score": float(cm.loc[s, rc]), "commot_cluster_pvalue": float(cp.loc[s, rc])})
    log.info("COMMOT cluster_communication done (%d rows), total %.1fs", len(rows), time.time() - t0)
    return pd.DataFrame(rows), pd.DataFrame(glob)


def run_spatialdm(adata, log, tag):
    import spatialdm as sdm
    a = ad.AnnData(X=adata.raw.X.copy(), obs=adata.obs[["cell_type_marker"]].copy(), var=pd.DataFrame(index=adata.raw.var_names))
    a.obsm["spatial"] = np.asarray(adata.obsm["spatial"], float)
    a.raw = a
    t0 = time.time()
    sdm.weight_matrix(a, l=SDM_L, cutoff=0.2, single_cell=True)
    sdm.extract_lr(a, "human", min_cell=3, datahost="package")
    sdm.spatialdm_global(a, n_perm=1000, method="z-score", nproc=1)
    res = a.uns["global_res"].copy()
    gi = a.uns["geneInter"]
    res["interaction_name"] = gi.loc[res.index, "interaction_name"].values if "interaction_name" in gi else res.index
    res["pathway_name"] = gi.loc[res.index, "pathway_name"].values if "pathway_name" in gi else ""
    lig_cols = [c for c in res.columns if c.startswith("Ligand")]
    rec_cols = [c for c in res.columns if c.startswith("Receptor")]
    res["ligand_genes"] = res[lig_cols].apply(lambda x: ";".join(sorted(v for v in x if isinstance(v, str))), axis=1)
    res["receptor_genes"] = res[rec_cols].apply(lambda x: ";".join(sorted(v for v in x if isinstance(v, str))), axis=1)
    res = res.drop(columns=lig_cols + rec_cols).reset_index().rename(columns={"index": "cellchat_id"})
    log.info("[%s] SpatialDM global: %d CellChatDB pairs tested in %.1fs", tag, len(res), time.time() - t0)
    return res


def main() -> None:
    log = start_run("phase4b_commot_spatialdm", extra_packages=["commot", "spatialdm"])
    status = {}
    crop = ad.read_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_benchmark_crop.h5ad")
    pairs = catalogue_pairs(list(crop.var_names))
    try:
        cl, gl = run_commot(crop, pairs, log)
        save_csv(cl, PHASE, "commot_cluster_crop.csv", "COMMOT cluster-level score/p per (L,R,sender,receiver), crop", log)
        save_csv(gl, PHASE, "commot_global_crop.csv", "COMMOT total signal per L-R pair, crop", log)
        status["commot"] = "DONE"
    except Exception as exc:  # record failure instead of faking a comparison
        log.exception("COMMOT failed: %s", exc)
        status["commot"] = f"FAILED: {type(exc).__name__}: {exc}"
    for tag, path in (("crop", DATA_DIR / "processed" / f"pdac_{TAG}_benchmark_crop.h5ad"),
                      ("full", DATA_DIR / "processed" / f"pdac_{TAG}_processed.h5ad")):
        try:
            a = crop if tag == "crop" else ad.read_h5ad(path)
            res = run_spatialdm(a, log, tag)
            save_csv(res, PHASE, f"spatialdm_global_{tag}.csv", f"SpatialDM global z / p per CellChatDB pair ({tag})", log)
            status[f"spatialdm_{tag}"] = "DONE"
        except Exception as exc:
            log.exception("SpatialDM (%s) failed: %s", tag, exc)
            status[f"spatialdm_{tag}"] = f"FAILED: {type(exc).__name__}: {exc}"
    save_json(status, PHASE, "phase4b_status.json", "Run status of COMMOT / SpatialDM (DONE or FAILED with error)", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

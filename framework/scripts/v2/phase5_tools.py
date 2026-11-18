"""v2 Phase 5 step 2 (pdac-bench env): run existing tools on the exported Phase 4 replicates.

Tools: squidpy gr.ligrec, SpatialDM, COMMOT (smallest section only), LIANA+ bivariate.
One resumable chunk per (section, setting, tool). Plan: framework/plans/v2_phase5.md.
Usage: python -I phase5_tools.py [--workers 4] [--tools squidpy spatialdm liana commot] [--max-reps N]
"""
from __future__ import annotations

import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[v] = "1"

import argparse
import json
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
REP = Path("F:/Thesis/data/v2_cache/phase5_replicates")
OUT = ROOT / "framework/results/v2_phase5/chunks"
SENDER, RECEIVER, LIG, REC = "CAF_fibroblast", "myeloid", "FN1", "ITGB1"


def make_adata(counts, xy, labels, genes):
    import anndata as ad
    import scipy.sparse as sp
    counts = counts.astype(np.float64)
    tot = counts.sum(1, keepdims=True)
    X = np.log1p(counts / np.maximum(tot, 1) * 1e4).astype(np.float32)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"cell_type": pd.Categorical(labels)},
                                                        index=[f"c{i}" for i in range(len(labels))]),
                   var=pd.DataFrame(index=list(genes)))
    a.obsm["spatial"] = np.asarray(xy, float)
    a.raw = a
    return a


def run_squidpy(a):
    """Amendment 1: at the default threshold (0.01) squidpy 1.6.5 in this env returns NaN for 167/169 cluster pairs,
    including CAF_fibroblast->myeloid, although FN1 and ITGB1 are detected in far more than 1 % of those cells. We run
    threshold=0 and apply the default 1 % expression rule ourselves for the claim pair, which is what a working default
    would do for this pair."""
    import squidpy as sq
    fl = float((a[a.obs.cell_type == SENDER, [LIG]].X.toarray() > 0).mean())
    fr = float((a[a.obs.cell_type == RECEIVER, [REC]].X.toarray() > 0).mean())
    res = sq.gr.ligrec(a, cluster_key="cell_type", interactions=pd.DataFrame([(LIG, REC)], columns=["source", "target"]),
                       copy=True, use_raw=True, seed=0, n_jobs=1, numba_parallel=False, show_progress_bar=False,
                       threshold=0.0)
    p = float(res["pvalues"].loc[(LIG, REC), (SENDER, RECEIVER)])
    if fl < 0.01 or fr < 0.01:
        p = np.nan
    return {"p": p, "significant": bool(np.isfinite(p) and p < 0.05),
            "stat": float(res["means"].loc[(LIG, REC), (SENDER, RECEIVER)]),
            "note": f"frac_L_in_sender={fl:.3f}; frac_R_in_receiver={fr:.3f}; threshold=0 workaround"}


def run_spatialdm(a):
    import anndata as ad
    import spatialdm as sdm
    b = ad.AnnData(X=a.X.copy(), obs=a.obs.copy(), var=a.var.copy())
    b.obsm["spatial"] = a.obsm["spatial"].copy()
    b.raw = b
    sdm.weight_matrix(b, l=25, single_cell=True)
    sdm.extract_lr(b, "human", datahost="package")
    lig, rec = b.uns["ligand"], b.uns["receptor"]
    ok = [i for i in lig.index if set(lig.loc[i].dropna()) == {LIG} and set(rec.loc[i].dropna()) == {REC}]
    if not ok:
        return {"p": np.nan, "significant": False, "note": "FN1-ITGB1 not in CellChatDB after subunit filtering"}
    i = ok[0]
    b.uns["ligand"], b.uns["receptor"] = lig.loc[[i]], rec.loc[[i]]
    b.uns["geneInter"] = b.uns["geneInter"].loc[[i]]
    b.uns["num_pairs"] = 1
    sdm.spatialdm_global(b)
    sdm.sig_pairs(b)
    r = b.uns["global_res"].iloc[0]
    return {"p": float(r["z_pval"]), "significant": bool(r["selected"]), "stat": float(r.get("z", np.nan)),
            "note": f"cellchat_id={i}; n_matching_entries={len(ok)}"}


def run_commot(a):
    import commot as ct
    df_lr = pd.DataFrame([[LIG, REC, f"{LIG}-{REC}"]], columns=["ligand", "receptor", "pathway"])
    ct.tl.spatial_communication(a, database_name="bench", df_ligrec=df_lr, dis_thr=50.0, heteromeric=False,
                                pathway_sum=False)
    ct.tl.cluster_communication(a, database_name="bench", lr_pair=(LIG, REC), clustering="cell_type")
    res = a.uns[f"commot_cluster-cell_type-bench-{LIG}-{REC}"]
    p = float(res["communication_pvalue"].loc[SENDER, RECEIVER])
    return {"p": p, "significant": bool(p < 0.05), "stat": float(res["communication_matrix"].loc[SENDER, RECEIVER])}


def run_liana(a):
    import liana as li
    li.ut.spatial_neighbors(a, bandwidth=25)
    res = li.mt.bivariate(a, resource=pd.DataFrame({"ligand": [LIG], "receptor": [REC]}), local_name="cosine",
                          global_name="morans", n_perms=100, use_raw=False)
    var = res.var
    row = var.iloc[0]
    p = float(row["morans_pvals"]) if "morans_pvals" in var else np.nan
    return {"p": p, "significant": bool(np.isfinite(p) and p < 0.05), "stat": float(row.get("morans", np.nan))}


TOOLS = {"squidpy": run_squidpy, "spatialdm": run_spatialdm, "commot": run_commot, "liana": run_liana}


def run_task(sec, scen, level, tool, reps):
    out = OUT / f"{sec}__{scen}__{level}__{tool}.csv"
    if out.exists():
        return f"skip {out.name}"
    g = np.load(REP / f"{sec}__geometry.npz", allow_pickle=False)
    rows = []
    t0 = time.time()
    for rep in reps:
        counts = np.load(REP / f"{sec}__{scen}__{level}__{rep}.npz")["counts"]
        a = make_adata(counts, g["xy"], g["labels"], g["genes"])
        t1 = time.time()
        try:
            r = TOOLS[tool](a)
            status = "ok"
        except Exception as e:
            r = {"p": np.nan, "significant": False, "note": f"{type(e).__name__}: {e}"[:300]}
            status = "error"
        rows.append({"section": sec, "scenario": scen, "level": level, "rep": rep, "tool": tool, "status": status,
                     "seconds": time.time() - t1, **r})
        del a
    pd.DataFrame(rows).to_csv(out, index=False)
    return f"done {out.name} in {time.time() - t0:.0f}s"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tools", nargs="*", default=["squidpy", "spatialdm", "liana", "commot"])
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--settings", nargs="*", default=None)
    ap.add_argument("--max-reps", type=int, default=None)
    ap.add_argument("--commot-reps", type=int, default=5)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(REP / "manifest.csv")
    man["level"] = man["level"].astype(str).replace({"nan": "None"})
    sel = json.loads((ROOT / "framework/results/v2_phase4/section_selection.json").read_text())
    tasks = []
    for (sec, scen, level), g in man.groupby(["section", "scenario", "level"], sort=False):
        if args.sections and sec not in args.sections:
            continue
        if args.settings and f"{scen}:{level}" not in args.settings:
            continue
        reps = sorted(g.rep.tolist())[: args.max_reps] if args.max_reps else sorted(g.rep.tolist())
        for tool in args.tools:
            if tool == "commot":
                if sec != sel["roles"]["smallest"]:
                    continue
                tasks.append((sec, scen, level, tool, reps[: args.commot_reps]))
                continue
            tasks.append((sec, scen, level, tool, reps))
    # COMMOT tasks are memory-heavy: run them serially after the others
    light = [t for t in tasks if t[3] != "commot"]
    heavy = [t for t in tasks if t[3] == "commot"]
    if light:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(run_task, *t): t for t in light}
            for f in as_completed(futs):
                try:
                    print(f.result(), flush=True)
                except Exception:
                    print("ERROR", futs[f][:4], traceback.format_exc()[-500:], flush=True)
    for t in heavy:
        try:
            print(run_task(*t), flush=True)
        except Exception:
            print("ERROR", t[:4], traceback.format_exc()[-500:], flush=True)


if __name__ == "__main__":
    main()

"""Phase 3: run the breast audit pipeline on a PDAC single-cell spatial dataset.

Steps (each writes tables + a figure under results/pdac/phase3/ and figures/pdac/phase3/):
  3a  load + QC (same thresholds as breast: min_counts=10 per cell, min_cells=3 per gene)
  3b  normalise/log1p, PCA, kNN, Leiden (res 0.5), two annotations:
        - marker scoring (per-cell argmax of mean marker expression)
        - Leiden consensus (cluster-mean marker score argmax)
  3c  panel coverage table for the curated L-R catalogue + per-gene detection counts
  3d  6-NN graph: compositional, edge, matched baseline, ratio; 500 label permutations (seed 42)
      for both annotations; uncorrected p + BH q across all tested rows
  3e  sensitivity: k in {4, 6, 8, 12}; radius in {15..300} um (sensitivity only); annotation rank shifts
Usage: python src/pdac/phase3_pdac_audit.py --dataset <key>
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA_DIR, SEED, finish_run, register, save_csv, save_figure, save_json,  # noqa: E402
                    save_npz, start_run)

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

import ccc_audit as ca  # noqa: E402
import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from spatial_null import generate_label_permutations  # noqa: E402

PHASE = "phase3"
N_PERMS = 500
SPARSE_MIN = 20
K_LIST = [4, 6, 8, 12]
RADII = [15.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0]
MAX_EDGES_RADIUS = 12_000_000  # cap on directed edges per radius graph (receiver subsampling beyond this)

# Marker sets declared before looking at the PDAC data; only genes present on the panel are used,
# a type needs >= 2 present markers to be scored.
PDAC_MARKERS = {
    "epithelial_ductal_tumor": ["KRT19", "KRT7", "EPCAM", "KRT8", "KRT18", "KRT17", "MUC1", "CEACAM6", "TFF1", "SOX9"],
    "acinar": ["PRSS1", "CPA1", "CELA3A", "CTRB2", "PNLIP", "CPB1", "REG1A"],
    "endocrine": ["INS", "GCG", "SST", "CHGA", "PPY"],
    "CAF_fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "PDGFRB", "SFRP2", "POSTN"],
    "T_cell": ["CD3D", "CD3E", "CD2", "CD3G", "TRAC", "CD8A", "CD4", "IL7R", "CD247"],
    "B_plasma": ["MS4A1", "CD79A", "CD19", "JCHAIN", "MZB1", "IGKC"],
    "myeloid": ["CD68", "CD163", "CD14", "LYZ", "CSF1R", "C1QA", "C1QB", "ITGAM", "FCGR3A"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "PLVAP", "KDR"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "pericyte_smc": ["RGS5", "ACTA2", "MYH11", "NOTCH3"],
}
HEADLINE = [("CXCL12", "CXCR4", "CAF_fibroblast", "T_cell"),
            ("CXCL12", "CXCR4", "CAF_fibroblast", "epithelial_ductal_tumor")]


# ----------------------------------------------------------------------------------------------
def qc(adata, log):
    before = (adata.n_obs, adata.n_vars)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    qc_before = adata.obs[["total_counts", "n_genes"]].describe().T
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=3)
    log.info("QC (breast thresholds min_counts=10, min_cells=3): %s -> %s", before, (adata.n_obs, adata.n_vars))
    return adata, before, qc_before


def preprocess(adata, log):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1), random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    log.info("Leiden res 0.5: %d clusters", adata.obs["leiden"].nunique())


def annotate(adata, log):
    present = {}
    for t, genes in PDAC_MARKERS.items():
        g = [x for x in genes if x in adata.var_names]
        if len(g) >= 2:
            present[t] = g
    dropped = sorted(set(PDAC_MARKERS) - set(present))
    log.info("Marker types used: %s; dropped (<2 markers on panel): %s", {k: len(v) for k, v in present.items()}, dropped)
    scores = {}
    for t, g in present.items():
        x = adata[:, g].X
        scores[t] = np.asarray(x.mean(axis=1)).ravel()
    sf = pd.DataFrame(scores, index=adata.obs_names)
    adata.obs["cell_type_marker"] = sf.idxmax(axis=1).astype("category")
    cons = sf.groupby(adata.obs["leiden"], observed=False).mean().idxmax(axis=1)
    adata.obs["cell_type_leiden"] = adata.obs["leiden"].map(cons).astype(str).astype("category")
    return present, dropped, sf


def knn_edges(coords, k):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    # neighbour i of query j -> directed edge (j -> neighbour) and the graph is used as Squidpy does:
    # row = cell, col = its neighbour. We keep (row=query cell, col=neighbour) exactly like obsp connectivities.
    rows = np.repeat(np.arange(coords.shape[0]), k)
    cols = idx[:, 1:].ravel()
    return rows.astype(np.int64), cols.astype(np.int64)


def radius_edges(coords, r, rng, log):
    n = coords.shape[0]
    nn = NearestNeighbors(radius=r).fit(coords)
    # estimate edge count from a sample of queries
    probe = rng.choice(n, size=min(n, 2000), replace=False)
    cnt = np.array([len(a) - 1 for a in nn.radius_neighbors(coords[probe], return_distance=False)])
    est = cnt.mean() * n
    frac = min(1.0, MAX_EDGES_RADIUS / max(est, 1))
    q = np.arange(n) if frac >= 1.0 else np.sort(rng.choice(n, size=int(frac * n), replace=False))
    adj = nn.radius_neighbors(coords[q], return_distance=False)
    rows_l, cols_l = [], []
    for j, neigh in zip(q, adj):
        neigh = neigh[neigh != j]
        rows_l.append(np.full(len(neigh), j))
        cols_l.append(neigh)
    rows = np.concatenate(rows_l).astype(np.int64)
    cols = np.concatenate(cols_l).astype(np.int64)
    log.info("radius %.0f um: est. %.3g edges; query-cell fraction %.3f; %d directed edges", r, est, frac, len(rows))
    return rows, cols, frac


def run_scores(adata, labels_str, rows, cols, pairs, log, tag, n_perms=N_PERMS, prod_dtype=np.float64,
               interactions=None):
    types = sorted(set(labels_str))
    ti = {t: k for k, t in enumerate(types)}
    lab = np.array([ti[x] for x in labels_str], dtype=np.int64)
    perms = [np.array([ti[x] for x in p], dtype=np.int64)
             for p in generate_label_permutations(labels_str, n_perms=n_perms, seed=SEED)]
    S = ca.setup(adata, pairs, rows, cols, prod_dtype=prod_dtype)
    t0 = time.time()
    obs, null, nze = ca.run_null(lab, len(types), S, perms)
    if interactions is None:
        interactions = [(s, r, l, rc) for (l, rc) in pairs for s in types for r in types]
    df, nulls = ca.tidy(interactions, types, S, obs, null, nze)
    log.info("[%s] %d interactions x %d perms in %.1fs (%d edges)", tag, len(df), n_perms, time.time() - t0, len(rows))
    return df, nulls


def add_flags(df, det):
    df["ligand_n_cells_detected"] = df["ligand"].map(det)
    df["receptor_n_cells_detected"] = df["receptor"].map(det)
    df["sparse_gene_diagnostic"] = df[["ligand_n_cells_detected", "receptor_n_cells_detected"]].min(axis=1) < SPARSE_MIN
    return df


# ----------------------------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--skip-radius", action="store_true")
    args = ap.parse_args()
    log = start_run(f"phase3_pdac_audit_{args.dataset}")
    rng = np.random.default_rng(SEED)
    log.info("Subsampling RNG for radius graphs: numpy default_rng(%d) (separate from the permutation stream)", SEED)

    # 3a load + QC -------------------------------------------------------------------------
    adata, meta = pdac_loaders.load(args.dataset, DATA_DIR, log)
    save_json(meta, PHASE, f"{args.dataset}_load_metadata.json", "Loader provenance: files, sizes, sample selection", log)
    adata, before, qc_before = qc(adata, log)
    qc_tab = qc_before.copy()
    qc_tab["stage"] = "before_qc"
    qa = adata.obs[["total_counts", "n_genes"]].describe().T
    qa["stage"] = "after_qc"
    qc_tab = pd.concat([qc_tab, qa])
    qc_tab.loc["shape_before", "count"] = f"{before[0]} cells x {before[1]} genes"
    qc_tab.loc["shape_after", "count"] = f"{adata.n_obs} cells x {adata.n_vars} genes"
    save_csv(qc_tab, PHASE, f"{args.dataset}_qc_table.csv", "QC summary before/after (breast thresholds)", log, index=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].hist(np.log10(adata.obs["total_counts"] + 1), bins=60, color="#6baed6")
    ax[0].set_xlabel("log10(total counts + 1) after QC")
    ax[1].hist(adata.obs["n_genes"], bins=60, color="#74c476")
    ax[1].set_xlabel("genes detected per cell after QC")
    xy = adata.obsm["spatial"]
    ax[2].scatter(xy[:, 0], xy[:, 1], s=0.2, c="k", rasterized=True)
    ax[2].set_aspect("equal")
    ax[2].set_title("cell centroids (µm)")
    fig.suptitle(f"{args.dataset}: QC ({adata.n_obs} cells x {adata.n_vars} genes after min_counts=10, min_cells=3)")
    save_figure(fig, PHASE, f"{args.dataset}_3a_qc", "QC distributions and centroid map", log)

    # 3b preprocess + annotate ------------------------------------------------------------
    preprocess(adata, log)
    present, dropped, sf = annotate(adata, log)
    ct = pd.crosstab(adata.obs["cell_type_marker"], adata.obs["cell_type_leiden"], margins=True)
    save_csv(ct, PHASE, f"{args.dataset}_annotation_crosstab.csv", "Marker-scoring vs Leiden-consensus labels", log, index=True)
    save_json({"markers_used": present, "types_dropped": dropped,
               "marker_counts": adata.obs["cell_type_marker"].value_counts().to_dict(),
               "leiden_consensus_counts": adata.obs["cell_type_leiden"].value_counts().to_dict(),
               "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
               "agreement_fraction": float((adata.obs["cell_type_marker"].astype(str)
                                            == adata.obs["cell_type_leiden"].astype(str)).mean())},
              PHASE, f"{args.dataset}_annotation_summary.json", "Annotation summary for both strategies", log)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for a, key in zip(axes, ["cell_type_marker", "cell_type_leiden"]):
        cats = adata.obs[key].astype(str)
        for k, c in enumerate(sorted(cats.unique())):
            m = (cats == c).to_numpy()
            a.scatter(xy[m, 0], xy[m, 1], s=0.3, label=f"{c} ({m.sum()})", color=plt.cm.tab10(k % 10), rasterized=True)
        a.set_aspect("equal")
        a.set_title(key)
        a.legend(markerscale=12, fontsize=7, loc="upper right")
    save_figure(fig, PHASE, f"{args.dataset}_3b_annotations_spatial", "Spatial maps of both annotations", log)

    # 3c coverage + detection --------------------------------------------------------------
    cat = catalogue()
    raw_counts = adata.layers["counts"]
    vn = list(adata.var_names)

    def ncells(g):
        if g not in vn:
            return 0
        col = raw_counts[:, vn.index(g)]
        return int((col.toarray() if sp.issparse(col) else np.asarray(col)).ravel().__gt__(0).sum())

    cov = pd.DataFrame([{"ligand": l, "receptor": r, "source_list": s, "ligand_on_panel": l in vn,
                         "receptor_on_panel": r in vn, "testable": (l in vn) and (r in vn),
                         "ligand_n_cells": ncells(l), "receptor_n_cells": ncells(r)} for l, r, s in cat])
    cov["sparse_gene_diagnostic"] = cov["testable"] & (cov[["ligand_n_cells", "receptor_n_cells"]].min(axis=1) < SPARSE_MIN)
    save_csv(cov, PHASE, f"{args.dataset}_lr_coverage.csv", "Present/absent coverage of curated L-R catalogue + detection counts", log)
    genes = sorted({g for l, r, _ in cat for g in (l, r)})
    det_rows = []
    for g in genes:
        row = {"gene": g, "on_panel": g in vn, "n_cells_detected": ncells(g), "n_cells_total": adata.n_obs}
        if g in vn:
            col = raw_counts[:, vn.index(g)]
            col = (col.toarray() if sp.issparse(col) else np.asarray(col)).ravel()
            for t in sorted(adata.obs["cell_type_marker"].astype(str).unique()):
                m = (adata.obs["cell_type_marker"].astype(str) == t).to_numpy()
                row[f"n_detected_{t}"] = int((col[m] > 0).sum())
        det_rows.append(row)
    det_df = pd.DataFrame(det_rows)
    det_df["role"] = det_df["gene"].map(lambda g: "+".join(sorted({("ligand" if g == l else "receptor")
                                                                    for l, r, _ in cat if g in (l, r)})))
    save_csv(det_df, PHASE, f"{args.dataset}_gene_detection_counts.csv",
             "Detection counts (raw counts > 0) per catalogue gene, overall and per marker cell type", log)
    det = dict(zip(det_df["gene"], det_df["n_cells_detected"]))
    pairs = [(l, r) for l, r, _ in cat if l in vn and r in vn]
    log.info("Testable L-R pairs: %d of %d catalogue pairs: %s", len(pairs), len(cat), pairs)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.25 * len(cov))))
    yy = np.arange(len(cov))
    ax.barh(yy - 0.2, cov["ligand_n_cells"], height=0.4, color="#fd8d3c", label="ligand cells detected")
    ax.barh(yy + 0.2, cov["receptor_n_cells"], height=0.4, color="#6baed6", label="receptor cells detected")
    ax.axvline(SPARSE_MIN, color="r", ls="--", lw=0.8, label=f"sparse threshold ({SPARSE_MIN} cells)")
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{l}-{r}{'' if t else '  (absent)'}" for l, r, t in cov[["ligand", "receptor", "testable"]].itertuples(index=False)], fontsize=7)
    ax.legend(fontsize=8)
    ax.set_title(f"{args.dataset}: curated L-R coverage and detection")
    save_figure(fig, PHASE, f"{args.dataset}_3c_lr_coverage_detection", "Coverage/detection of curated L-R pairs", log)

    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    # 3d main analysis: 6-NN, both annotations ----------------------------------------------
    rows6, cols6 = knn_edges(coords, 6)
    main_tabs = {}
    for ann in ("cell_type_marker", "cell_type_leiden"):
        labels = np.asarray(adata.obs[ann].astype(str).to_numpy())
        df, nulls = run_scores(adata, labels, rows6, cols6, pairs, log, f"6NN/{ann}")
        df = add_flags(df, det)
        df.insert(0, "annotation", ann)
        main_tabs[ann] = df
        save_npz(PHASE, f"{args.dataset}_null_6NN_{ann}.npz",
                 f"Null distributions ({N_PERMS} perms) for comp/edge/matched/ratio, 6-NN, {ann}; column order = CSV row order",
                 log, **{f"{k}_null": v.astype(np.float32) for k, v in nulls.items()})
    allmain = pd.concat(main_tabs.values(), ignore_index=True)
    # BH across ALL rows tested (both annotations pooled, per score) + per annotation
    for key in ca.SCORES:
        allmain[f"{key}_q_bh_all_rows"] = ca.bh(allmain[f"{key}_p_plus1"].to_numpy())
        allmain[f"{key}_q_bh_within_annotation"] = allmain.groupby("annotation")[f"{key}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    allmain = allmain.drop(columns=[f"{k}_q_bh" for k in ca.SCORES])
    save_csv(allmain, PHASE, f"{args.dataset}_6NN_scores_both_annotations.csv",
             "Main table: C/E/M/Q obs, null stats, uncorrected p, BH q (all rows; within annotation), flags; 6-NN", log)
    n_tested = {ann: int(np.isfinite(d["edge_p_plus1"]).sum()) for ann, d in main_tabs.items()}
    n_tested_comp = {ann: int(np.isfinite(d["comp_p_plus1"]).sum()) for ann, d in main_tabs.items()}
    log.info("Rows tested: compositional %s; edge (>=1 qualifying edge) %s", n_tested_comp, n_tested)

    # headline extract
    hl = allmain[allmain.apply(lambda r: (r.ligand, r.receptor, r.sender, r.receiver) in HEADLINE, axis=1)]
    save_csv(hl, PHASE, f"{args.dataset}_headline_CXCL12_CXCR4.csv", "Headline CXCL12-CXCR4 CAF->T and CAF->tumour rows", log)

    # annotation sensitivity: rank shifts on the shared (L,R,sender,receiver) rows
    k4 = ["sender", "receiver", "ligand", "receptor"]
    m = main_tabs["cell_type_marker"].merge(main_tabs["cell_type_leiden"], on=k4, suffixes=("_marker", "_leiden"))
    ann_rows = []
    for key in ca.SCORES:
        a, b = m[f"{key}_z_marker"], m[f"{key}_z_leiden"]
        ok = np.isfinite(a) & np.isfinite(b)
        m[f"{key}_rank_marker"] = a.rank(ascending=False)
        m[f"{key}_rank_leiden"] = b.rank(ascending=False)
        m[f"{key}_rank_shift"] = m[f"{key}_rank_leiden"] - m[f"{key}_rank_marker"]
        ann_rows.append({"score": key, "n_shared_rows": int(ok.sum()),
                         "spearman_z": float(spearmanr(a[ok], b[ok]).statistic) if ok.sum() > 2 else np.nan,
                         "median_abs_rank_shift": float(m.loc[ok, f"{key}_rank_shift"].abs().median()) if ok.any() else np.nan,
                         "n_verdict_changes": int((m[f"{key}_sig_marker"] != m[f"{key}_sig_leiden"]).sum())})
    save_csv(m, PHASE, f"{args.dataset}_annotation_sensitivity_rows.csv", "Per-row z/rank under both annotations (shared rows)", log)
    ann_sum = pd.DataFrame(ann_rows)
    save_csv(ann_sum, PHASE, f"{args.dataset}_annotation_sensitivity_summary.csv", "Spearman / rank shifts / verdict changes between annotations", log)

    # 3e k sensitivity (marker annotation, all testable rows) -------------------------------
    labels_m = np.asarray(adata.obs["cell_type_marker"].astype(str).to_numpy())
    k_tabs = []
    for k in K_LIST:
        rk, ck = knn_edges(coords, k)
        dfk, nk = run_scores(adata, labels_m, rk, ck, pairs, log, f"k={k}")
        dfk.insert(0, "k", k)
        k_tabs.append(add_flags(dfk, det))
        save_npz(PHASE, f"{args.dataset}_null_k{k}_marker.npz", f"Null distributions, k={k}-NN, marker labels", log,
                 **{f"{kk}_null": v.astype(np.float32) for kk, v in nk.items() if kk in ("edge", "ratio")})
    ksens = pd.concat(k_tabs, ignore_index=True)
    save_csv(ksens, PHASE, f"{args.dataset}_k_sensitivity.csv", "Edge/ratio scores for k in {4,6,8,12} (marker labels)", log)

    # 3e radius sensitivity (marker labels; testable rows) -----------------------------------
    rad_tabs = []
    if not args.skip_radius:
        for r in RADII:
            rr, cr, frac = radius_edges(coords, r, rng, log)
            dfr, nr = run_scores(adata, labels_m, rr, cr, pairs, log, f"r={r}", prod_dtype=np.float32)
            dfr.insert(0, "radius_um", r)
            dfr.insert(1, "query_cell_fraction", frac)
            rad_tabs.append(add_flags(dfr, det))
            save_npz(PHASE, f"{args.dataset}_null_radius{int(r)}_marker.npz",
                     f"Null distributions, radius {r} um (SENSITIVITY ONLY), marker labels", log,
                     **{f"{kk}_null": v.astype(np.float32) for kk, v in nr.items() if kk in ("edge", "ratio")})
        rsens = pd.concat(rad_tabs, ignore_index=True)
        save_csv(rsens, PHASE, f"{args.dataset}_radius_sensitivity.csv",
                 "SENSITIVITY ONLY: edge/ratio scores for radius graphs 15-300 um (marker labels)", log)

    # figures -------------------------------------------------------------------------------
    make_main_figure(allmain, args.dataset, log)
    make_sensitivity_figure(ksens, rad_tabs, args.dataset, log)
    adata.write_h5ad(DATA_DIR / "processed" / f"pdac_{args.dataset}_processed.h5ad")
    log.info("Processed AnnData written to data/processed/pdac_%s_processed.h5ad (git-ignored)", args.dataset)
    finish_run(log, PHASE)


def make_main_figure(allmain, ds, log):
    for ann, d in allmain.groupby("annotation"):
        d = d[np.isfinite(d["comp_z"])].copy()
        top = d.reindex(d["comp_z"].abs().sort_values(ascending=False).index).head(40).sort_values("comp_z")
        y = np.arange(len(top))
        fig, ax = plt.subplots(figsize=(12, 0.3 * len(top) + 2))
        for key, c, mk in (("comp", "#7f7f7f", "o"), ("matched", "#1f77b4", "D"), ("edge", "#d62728", "s"), ("ratio", "#2ca02c", "^")):
            ax.scatter(top[f"{key}_z"], y, c=c, marker=mk, s=30, label=key)
        for i, r in enumerate(top.itertuples()):
            if r.sparse_gene_diagnostic:
                ax.text(1.01, i, "sparse", transform=ax.get_yaxis_transform(), fontsize=6, color="#b15928", va="center")
        ax.set_xscale("symlog", linthresh=2)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels((top["ligand"] + "→" + top["receptor"] + "  " + top["sender"] + "→" + top["receiver"]).tolist(), fontsize=7)
        ax.legend(fontsize=8)
        ax.set_title(f"{ds} 6-NN, {ann}: top 40 rows by |compositional z| (500 perms)")
        fig.tight_layout()
        save_figure(fig, PHASE, f"{ds}_3d_z_comparison_6NN_{ann}", f"z of C/M/E/Q for top rows, {ann}", log)


def make_sensitivity_figure(ksens, rad_tabs, ds, log):
    sel = [(l, r, s, rc) for l, r, s, rc in HEADLINE]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for (l, r, s, rc) in sel:
        d = ksens[(ksens.ligand == l) & (ksens.receptor == r) & (ksens.sender == s) & (ksens.receiver == rc)]
        if len(d):
            axes[0].plot(d["k"], d["edge_z"], "o-", label=f"{l}→{r} {s}→{rc} edge z")
            axes[0].plot(d["k"], d["ratio_z"], "^--", label=f"{l}→{r} {s}→{rc} ratio z")
        if rad_tabs:
            rs = pd.concat(rad_tabs)
            d = rs[(rs.ligand == l) & (rs.receptor == r) & (rs.sender == s) & (rs.receiver == rc)]
            axes[1].plot(d["radius_um"], d["edge_z"], "o-", label=f"{s}→{rc} edge z")
            axes[1].plot(d["radius_um"], d["ratio_z"], "^--", label=f"{s}→{rc} ratio z")
    axes[0].set_xlabel("k (kNN graph)")
    axes[1].set_xlabel("radius (µm) — SENSITIVITY ONLY")
    for a in axes:
        a.axhline(1.645, color="#aaa", ls=":")
        a.set_ylabel("z")
        a.legend(fontsize=7)
    fig.suptitle(f"{ds}: headline CXCL12→CXCR4 sensitivity to k and radius (marker labels)")
    save_figure(fig, PHASE, f"{ds}_3e_headline_k_radius_sensitivity", "Headline rows vs k and radius", log)


if __name__ == "__main__":
    main()

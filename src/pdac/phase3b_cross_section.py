"""Phase 3b: cross-section consistency check on all 16 PDAC-SMI sections.

For every section independently (labels are permuted within the section only):
  load -> QC (breast thresholds) -> normalise/PCA/Leiden (breast calls) -> marker + Leiden-consensus labels
  -> 6-NN graph -> C/E/M/Q with 500 label permutations (seed 42).
Summarises, per (L, R, sender, receiver, annotation), in how many sections each score is significant,
and the median z. Sections are separate tumours/slides; this is NOT a pooled analysis.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import ccc_audit as ca  # noqa: E402
import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import HEADLINE, SPARSE_MIN, add_flags, annotate, knn_edges, preprocess, qc, run_scores  # noqa: E402

PHASE = "phase3"


def main() -> None:
    log = start_run("phase3b_cross_section")
    full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    sections = sorted(full.obs["sample_id"].astype(str).unique())
    del full
    log.info("Sections: %s", sections)
    tabs, sec_info = [], []
    for s in sections:
        t0 = time.time()
        adata, meta = pdac_loaders.load(f"mendeley_pdac_smi:{s}", DATA_DIR, log)
        adata, before, _ = qc(adata, log)
        preprocess(adata, log)
        present, dropped, _ = annotate(adata, log)
        vn = list(adata.var_names)
        pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
        cnt = adata.layers["counts"]
        det = {g: int((cnt[:, vn.index(g)].toarray().ravel() > 0).sum()) for g in sorted({x for p in pairs for x in p})}
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        rows, cols = knn_edges(coords, 6)
        for ann in ("cell_type_marker", "cell_type_leiden"):
            labels = np.asarray(adata.obs[ann].astype(str).to_numpy())
            df, _ = run_scores(adata, labels, rows, cols, pairs, log, f"{s}/{ann}")
            df = add_flags(df, det)
            df.insert(0, "annotation", ann)
            df.insert(0, "section", s)
            df.insert(1, "treatment_status", meta["treatment_status"])
            tabs.append(df)
        info = {"section": s, "treatment_status": meta["treatment_status"], "donor": meta["donor"],
                "n_cells_before_qc": before[0], "n_cells_after_qc": int(adata.n_obs), "n_genes_after_qc": int(adata.n_vars),
                "n_testable_pairs": len(pairs), "types_dropped": ";".join(dropped),
                **{f"n_{t}": int(v) for t, v in adata.obs["cell_type_marker"].value_counts().items()},
                "runtime_s": round(time.time() - t0, 1)}
        sec_info.append(info)
        log.info("section %s done: %s", s, info)
    allt = pd.concat(tabs, ignore_index=True)
    for key in ca.SCORES:
        allt[f"{key}_q_bh_within_section_annotation"] = allt.groupby(["section", "annotation"])[f"{key}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    allt = allt.drop(columns=[f"{k}_q_bh" for k in ca.SCORES])
    save_csv(allt, PHASE, "cross_section_6NN_all_rows.csv",
             "All 16 sections x 2 annotations: C/E/M/Q obs, z, p (+1), BH q within section/annotation; 6-NN, 500 perms", log)
    save_csv(pd.DataFrame(sec_info), PHASE, "cross_section_section_info.csv", "Per-section QC and cell-type counts", log)

    k = ["annotation", "ligand", "receptor", "sender", "receiver"]
    g = allt.groupby(k)
    summ = g.agg(n_sections=("section", "nunique"),
                 n_edge_evaluable=("edge_status", lambda x: int((x == "evaluated").sum())),
                 n_comp_sig=("comp_sig", "sum"), n_matched_sig=("matched_sig", "sum"),
                 n_edge_sig=("edge_sig", "sum"), n_ratio_sig=("ratio_sig", "sum"),
                 n_edge_sig_bh=("edge_q_bh_within_section_annotation", lambda x: int((x < 0.05).sum())),
                 n_ratio_sig_bh=("ratio_q_bh_within_section_annotation", lambda x: int((x < 0.05).sum())),
                 median_comp_z=("comp_z", "median"), median_matched_z=("matched_z", "median"),
                 median_edge_z=("edge_z", "median"), median_ratio_z=("ratio_z", "median"),
                 median_ratio_obs=("ratio_obs", "median"),
                 n_sections_sparse_gene=("sparse_gene_diagnostic", "sum")).reset_index()
    save_csv(summ, PHASE, "cross_section_summary_by_interaction.csv",
             "Per interaction & annotation: number of sections significant under C/M/E/Q and median z", log)
    hl = allt[allt.apply(lambda r: (r.ligand, r.receptor, r.sender, r.receiver) in HEADLINE, axis=1)]
    save_csv(hl, PHASE, "cross_section_headline_CXCL12_CXCR4.csv", "Headline CXCL12-CXCR4 rows in every section", log)
    hs = summ[summ.apply(lambda r: (r.ligand, r.receptor, r.sender, r.receiver) in HEADLINE, axis=1)]
    save_json(hs.to_dict(orient="records"), PHASE, "cross_section_headline_summary.json",
              "Headline rows: number of sections significant under each score", log)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    for ax, (l, r, s_, rc) in zip(axes, HEADLINE):
        d = hl[(hl.sender == s_) & (hl.receiver == rc) & (hl.annotation == "cell_type_marker")].sort_values("section")
        x = np.arange(len(d))
        for key, c, mk, off in (("comp", "#7f7f7f", "o", -0.27), ("matched", "#1f77b4", "D", -0.09),
                                ("edge", "#d62728", "s", 0.09), ("ratio", "#2ca02c", "^", 0.27)):
            ax.scatter(x + off, d[f"{key}_z"], c=c, marker=mk, s=36, label=key)
        ax.axhline(1.645, color="#aaa", ls=":")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{a}\n({b})" for a, b in zip(d["section"], d["treatment_status"])], fontsize=7, rotation=90)
        ax.set_yscale("symlog", linthresh=2)
        ax.set_title(f"{l}→{r}  {s_}→{rc} (marker labels, 6-NN)", fontsize=10)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("z vs within-section label permutation (500)")
    fig.suptitle("Cross-section consistency check: headline axis in 16 PDAC-SMI sections")
    fig.tight_layout()
    save_figure(fig, PHASE, "cross_section_headline_CXCL12_CXCR4_z_by_section", "Headline z per section and score", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

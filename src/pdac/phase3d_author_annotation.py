"""Phase 3d: annotation sensitivity against the dataset authors' cell-type labels (reference annotation).

Motivation (from phase3 crosstab, U5-a): our marker-scoring 'T_cell' class overlaps poorly with the authors'
T-cell calls, so the headline CAF->T row may depend on label purity. Here the same 6-NN audit
(C/E/M/Q, 500 label permutations, seed 42, within section) is run with the authors' `annotation_majortypes`,
mapped to our coarse type names by the fixed table below (declared before running).
Low memory: no PCA/Leiden/radius graphs.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

import ccc_audit as ca  # noqa: E402
import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import HEADLINE, add_flags, knn_edges, qc, run_scores  # noqa: E402

PHASE = "phase3"
AUTHOR_MAP = {
    "CAF": "CAF_fibroblast", "Malignant": "epithelial_ductal_tumor", "ductal": "epithelial_ductal_tumor",
    "CD8+ T": "T_cell", "CD4+ T": "T_cell", "Treg": "T_cell", "Natural killer": "NK",
    "Macrophage": "myeloid", "Dendritic (activated)": "myeloid", "Dendritic (conventional)": "myeloid",
    "Dendritic (plasmacytoid)": "myeloid", "Neutrophil": "myeloid",
    "B": "B_plasma", "Plasma": "B_plasma", "Endothelial": "endothelial", "Pericyte": "pericyte_smc",
    "Smooth muscle": "pericyte_smc", "Mast": "mast", "Endocrine": "endocrine", "Schwann/nerve": "schwann_nerve",
    "Noise": "unassigned", "uninformative_EVEN": "unassigned",
}


def main() -> None:
    log = start_run("phase3d_author_annotation")
    full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    sections = sorted(full.obs["sample_id"].astype(str).unique())
    unmapped = sorted(set(full.obs["annotation_majortypes"].astype(str)) - set(AUTHOR_MAP))
    del full
    if unmapped:
        raise ValueError(f"unmapped author labels: {unmapped}")
    tabs = []
    for s in sections:
        adata, meta = pdac_loaders.load(f"mendeley_pdac_smi:{s}", DATA_DIR, log)
        adata, _, _ = qc(adata, log)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        labels = adata.obs["author_majortype"].map(AUTHOR_MAP).astype(str).to_numpy()
        vn = list(adata.var_names)
        pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
        det = {g: int((adata.layers["counts"][:, vn.index(g)].toarray().ravel() > 0).sum()) for g in {x for p in pairs for x in p}}
        rows, cols = knn_edges(np.asarray(adata.obsm["spatial"], float), 6)
        df, _ = run_scores(adata, labels, rows, cols, pairs, log, f"{s}/author")
        df = add_flags(df, det)
        df.insert(0, "annotation", "author_mapped")
        df.insert(0, "section", s)
        df.insert(1, "treatment_status", meta["treatment_status"])
        tabs.append(df)
        log.info("section %s author-mapped counts: %s", s, pd.Series(labels).value_counts().to_dict())
    allt = pd.concat(tabs, ignore_index=True)
    for key in ca.SCORES:
        allt[f"{key}_q_bh_within_section"] = allt.groupby("section")[f"{key}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    allt = allt.drop(columns=[f"{k}_q_bh" for k in ca.SCORES])
    save_csv(allt, PHASE, "author_annotation_6NN_all_sections.csv",
             "6-NN C/E/M/Q (500 perms) under authors' mapped annotation, all 16 sections", log)
    hl = allt[allt.apply(lambda r: (r.ligand, r.receptor, r.sender, r.receiver) in HEADLINE, axis=1)]
    save_csv(hl, PHASE, "author_annotation_headline_CXCL12_CXCR4.csv", "Headline rows per section, author labels", log)
    summ = hl.groupby(["sender", "receiver"]).agg(
        n_sections=("section", "nunique"), n_comp_sig=("comp_sig", "sum"), n_matched_sig=("matched_sig", "sum"),
        n_edge_sig=("edge_sig", "sum"), n_ratio_sig=("ratio_sig", "sum"),
        n_edge_sig_bh=("edge_q_bh_within_section", lambda x: int((x < 0.05).sum())),
        n_ratio_sig_bh=("ratio_q_bh_within_section", lambda x: int((x < 0.05).sum())),
        median_edge_z=("edge_z", "median"), median_ratio_z=("ratio_z", "median"),
        median_ratio_obs=("ratio_obs", "median")).reset_index()
    save_csv(summ, PHASE, "author_annotation_headline_summary.csv", "Headline: sections significant under author labels", log)
    ev = allt[allt.edge_status == "evaluated"]
    save_json({"map": AUTHOR_MAP, "n_section_rows": int(len(allt)), "n_edge_sig": int(ev.edge_sig.sum()),
               "n_edge_sig_and_ratio_sig": int((ev.edge_sig & ev.ratio_sig).sum()),
               "headline": summ.to_dict(orient="records")}, PHASE, "author_annotation_summary.json",
              "Author-annotation sensitivity summary", log)

    # compare headline across annotations (marker / leiden from phase3b vs author)
    xs = pd.read_csv(Path(__file__).resolve().parents[2] / "results" / "pdac" / "phase3" / "cross_section_headline_CXCL12_CXCR4.csv")
    comb = pd.concat([xs[["section", "annotation", "sender", "receiver", "edge_z", "ratio_z", "edge_sig", "ratio_sig"]],
                      hl[["section", "annotation", "sender", "receiver", "edge_z", "ratio_z", "edge_sig", "ratio_sig"]]])
    save_csv(comb, PHASE, "headline_by_annotation_all_sections.csv", "Headline edge/ratio z per section under marker, Leiden, author labels", log)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), sharey=False)
    for ax, (l, r, s_, rc) in zip(axes, HEADLINE):
        d = comb[(comb.sender == s_) & (comb.receiver == rc)]
        secs = sorted(d.section.unique())
        x = {s: i for i, s in enumerate(secs)}
        for ann, c, off in (("cell_type_marker", "#d62728", -0.2), ("cell_type_leiden", "#9467bd", 0.0), ("author_mapped", "#1f77b4", 0.2)):
            dd = d[d.annotation == ann]
            ax.scatter([x[s] + off for s in dd.section], dd.edge_z, c=c, marker="s", s=30, label=f"{ann} edge z")
            ax.scatter([x[s] + off for s in dd.section], dd.ratio_z, c=c, marker="^", s=30, alpha=0.6, label=f"{ann} ratio z")
        ax.axhline(1.645, color="#aaa", ls=":")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_yscale("symlog", linthresh=2)
        ax.set_xticks(range(len(secs)))
        ax.set_xticklabels(secs, rotation=90, fontsize=7)
        ax.set_title(f"{l}→{r} {s_}→{rc}", fontsize=10)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle("Headline axis under three annotations, 16 sections (6-NN, 500 perms); Leiden has no T-cell class")
    fig.tight_layout()
    save_figure(fig, PHASE, "headline_edge_ratio_z_by_annotation_all_sections", "Headline z by annotation and section", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

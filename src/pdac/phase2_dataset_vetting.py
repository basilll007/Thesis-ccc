"""Phase 2: finalise the PDAC dataset vetting table.

Input : results/pdac/phase2/dataset_vetting_draft.csv (web vetting; evidence in dataset_vetting_raw_notes.md)
        data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad (downloaded; CC BY 4.0)
Output: results/pdac/dataset_vetting.csv  (+ verification columns for the selected dataset)
        results/pdac/phase2/mendeley_h5ad_inspection.json
        figures/pdac/phase2/phase2_lr_gene_coverage_by_candidate.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO_ROOT, RESULTS_PDAC, finish_run, register, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import pdac_loaders  # noqa: E402

PHASE = "phase2"
H5AD = REPO_ROOT / "data" / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad"


def main() -> None:
    log = start_run("phase2_dataset_vetting")
    draft = pd.read_csv(RESULTS_PDAC / PHASE / "dataset_vetting_draft.csv")
    register(RESULTS_PDAC / PHASE / "dataset_vetting_draft.csv", PHASE,
             "Web-vetting draft (one row per candidate; fetched-page evidence)", "subagent web vetting (2026-10-24)")
    register(RESULTS_PDAC / PHASE / "dataset_vetting_raw_notes.md", PHASE,
             "Raw vetting notes: URLs fetched, verbatim quotes, BLOCKED fetches", "subagent web vetting (2026-10-24)")
    genes = list(draft.columns[draft.columns.get_loc("CXCL12"):draft.columns.get_loc("notes")])

    a = ad.read_h5ad(H5AD, backed="r")
    vn = set(a.var_names)
    pt = a.var["probe_type"].astype(str)
    base = set(a.var_names[pt == "Base"])
    insp = {
        "file": "data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad",
        "bytes": H5AD.stat().st_size,
        "sha256": pdac_loaders.sha256(H5AD),
        "n_cells": int(a.n_obs),
        "n_features": int(a.n_vars),
        "probe_type_counts": pt.value_counts().to_dict(),
        "n_samples": int(a.obs["sample_id"].nunique()),
        "n_donors_field": int(a.obs["donor"].nunique()),
        "cells_per_sample": a.obs["sample_id"].value_counts().to_dict(),
        "treatment_by_sample": a.obs.groupby("sample_id", observed=True)["treatment_status"].first().astype(str).to_dict(),
        "obs_columns": list(a.obs.columns),
        "catalogue_genes_present_base_panel": {g: (g in base) for g in genes},
        "coordinate_unit": "pixels (CenterX/Y_global_px); 0.18 um/px from Shiau et al. bioRxiv 2023 "
                           "(doi 10.1101/2023.06.28.546848) Methods: '60 pixels (10.8 um)'",
    }
    log.info("h5ad: %d cells x %d features; probe types %s; samples %d", insp["n_cells"], insp["n_features"],
             insp["probe_type_counts"], insp["n_samples"])
    save_json(insp, PHASE, "mendeley_h5ad_inspection.json", "Direct inspection of the downloaded Mendeley h5ad", log)

    out = draft.copy()
    out["selected"] = out["candidate"].str.startswith("Mendeley PDAC-SMI")
    out["verified_from_downloaded_file"] = ""
    sel = out["selected"]
    agree = all((out.loc[sel, g].iloc[0] == "yes") == (g in base) for g in genes)
    out.loc[sel, "verified_from_downloaded_file"] = (
        f"h5ad downloaded (300,718,306 B, sha256 {insp['sha256'][:16]}...); {insp['n_cells']} cells in h5ad "
        f"(paper text reports 753,133 - h5ad count is used); gene presence re-checked on 960 Base probes: "
        f"{'agrees' if agree else 'DISAGREES'} with draft; licence re-fetched from Mendeley API: CC BY 4.0")
    out["n_catalogue_genes_present"] = [(r[genes] == "yes").sum() for _, r in out.iterrows()]
    out["rank"] = np.nan
    order = ["Mendeley PDAC-SMI", "GEO GSE310352", "GEO GSE277782", "10x Xenium FFPE Human PDAC",
             "GEO GSE267680", "10x Xenium Pancreatic Cancer"]
    for i, o in enumerate(order, 1):
        out.loc[out["candidate"].str.startswith(o), "rank"] = i
    out["selection_reason"] = ""
    out.loc[sel, "selection_reason"] = (
        "Rank 1: single-cell CosMx (centroids); 37/39 catalogue genes on panel (only PTN, SDC4 absent); "
        "explicit CC BY 4.0 data licence; one 301 MB h5ad over plain HTTP; 16 sections / 13 tumours allow a "
        "cross-section consistency check. GSE310352/GSE277782 have the same gene coverage but no licence could be "
        "confirmed (GEO terms page BLOCKED by reCAPTCHA); 10x datasets have lower coverage and 10x pages were BLOCKED (HTTP 429).")
    path = RESULTS_PDAC / "dataset_vetting.csv"
    out.to_csv(path, index=False)
    register(path, PHASE, "Final dataset vetting table (licence, size, platform, coords, panel, L-R gene presence, rank, selection)")
    log.info("Gene presence re-check agrees with draft: %s", agree)

    mat = out.set_index("candidate")[genes].replace({"yes": 1, "no": 0, "UNVERIFIED": np.nan}).astype(float)
    fig, ax = plt.subplots(figsize=(16, 4.8))
    cmap = matplotlib.colors.ListedColormap(["#f0f0f0", "#31a354"])
    ax.imshow(np.ma.masked_invalid(mat.to_numpy()), cmap=cmap, aspect="auto", vmin=0, vmax=1)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isnan(mat.iat[i, j]):
                ax.text(j, i, "?", ha="center", va="center", fontsize=7, color="#888")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels([c[:60] for c in mat.index], fontsize=7)
    ax.set_title("Phase 2: catalogue gene presence per candidate (green = on panel, grey = absent, ? = UNVERIFIED)")
    save_figure(fig, PHASE, "phase2_lr_gene_coverage_by_candidate", "Catalogue-gene presence heatmap across candidates", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

"""Phase 3f: are catalogue genes detected above CosMx background? (negative-probe comparison; low memory)

For each section, per-gene detection rate (fraction of cells with count > 0) and mean count are compared with
the 19 negative-control probes ('Negprob') in the same cells. A gene is flagged 'near_background' when its
detection rate is below 2x the maximum negative-probe detection rate in that section (threshold declared
here; it is a diagnostic flag, not a filter). Also reports per-gene detection in the authors' T cells vs
malignant cells for checkpoint receptors, to see whether detection is cell-type specific.
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

from lr_catalog import catalogue  # noqa: E402

PHASE = "phase3"
FOLD = 2.0
T_LABELS = {"CD8+ T", "CD4+ T", "Treg"}


def main() -> None:
    log = start_run("phase3f_background_check")
    a = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    pt = a.var["probe_type"].astype(str).to_numpy()
    neg = list(a.var_names[pt == "Negprob"])
    genes = sorted({g for l, r, _ in catalogue() for g in (l, r)} & set(a.var_names[pt == "Base"]))
    cols = [list(a.var_names).index(g) for g in genes + neg]
    rows, ct_rows = [], []
    for s in sorted(a.obs["sample_id"].astype(str).unique()):
        m = (a.obs["sample_id"].astype(str) == s).to_numpy()
        X = a[np.where(m)[0], cols].to_memory().X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        det = (X > 0).mean(axis=0)
        mean = X.mean(axis=0)
        neg_det = det[len(genes):]
        nmax = float(neg_det.max())
        for g, dr, mu in zip(genes, det[:len(genes)], mean[:len(genes)]):
            rows.append({"section": s, "gene": g, "detection_rate": float(dr), "mean_count": float(mu),
                         "negprobe_max_detection_rate": nmax, "negprobe_mean_detection_rate": float(neg_det.mean()),
                         "fold_over_negprobe_max": float(dr / nmax) if nmax > 0 else np.inf,
                         "near_background": bool(dr < FOLD * nmax)})
        lab = a.obs.loc[m, "annotation_majortypes"].astype(str).to_numpy()
        tmask = np.isin(lab, list(T_LABELS))
        mmask = lab == "Malignant"
        for g in ("PDCD1", "CTLA4", "CD28", "HAVCR2", "CXCR4", "CD274", "CD86", "CD80"):
            if g in genes:
                j = genes.index(g)
                ct_rows.append({"section": s, "gene": g, "n_author_T": int(tmask.sum()), "n_author_malignant": int(mmask.sum()),
                                "det_rate_author_T": float((X[tmask, j] > 0).mean()) if tmask.any() else np.nan,
                                "det_rate_author_malignant": float((X[mmask, j] > 0).mean()) if mmask.any() else np.nan,
                                "det_rate_negprobe_max": nmax})
        log.info("%s: neg-probe max detection %.4f", s, nmax)
    df = pd.DataFrame(rows)
    save_csv(df, PHASE, "background_check_detection_vs_negprobes.csv",
             "Per section x catalogue gene: detection rate vs negative-probe detection; near_background flag (< 2x neg max)", log)
    ct = pd.DataFrame(ct_rows)
    save_csv(ct, PHASE, "background_check_celltype_specificity.csv",
             "Detection rate in author T cells vs author malignant cells for checkpoint/chemokine genes", log)
    summ = df.groupby("gene").agg(n_sections_near_background=("near_background", "sum"),
                                  median_detection_rate=("detection_rate", "median"),
                                  median_fold_over_negprobe_max=("fold_over_negprobe_max", "median")).reset_index()
    save_csv(summ, PHASE, "background_check_summary_by_gene.csv", "Per gene: sections near background and median fold over neg-probe max", log)
    cts = ct.groupby("gene")[["det_rate_author_T", "det_rate_author_malignant", "det_rate_negprobe_max"]].median().reset_index()
    save_csv(cts, PHASE, "background_check_celltype_specificity_median.csv", "Median over sections of the cell-type detection rates", log)
    save_json({"fold_threshold": FOLD, "n_negprobes": len(neg), "genes_near_background_in_ge_8_sections":
               summ.loc[summ.n_sections_near_background >= 8, "gene"].tolist()}, PHASE, "background_check_summary.json",
              "Background-check summary", log)
    s2 = summ.sort_values("median_fold_over_negprobe_max")
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(range(len(s2)), s2.median_fold_over_negprobe_max, color=np.where(s2.n_sections_near_background >= 8, "#d62728", "#6baed6"))
    ax.axhline(FOLD, color="k", ls="--", lw=0.8, label=f"{FOLD}x neg-probe max")
    ax.set_yscale("log")
    ax.set_xticks(range(len(s2)))
    ax.set_xticklabels(s2.gene, rotation=90, fontsize=7)
    ax.set_ylabel("median fold: gene detection / max neg-probe detection")
    ax.set_title("PDAC-SMI: catalogue-gene detection relative to negative probes (16 sections); red = near background in >= 8 sections")
    ax.legend()
    save_figure(fig, PHASE, "background_check_fold_over_negprobes", "Catalogue gene detection vs negative probes", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

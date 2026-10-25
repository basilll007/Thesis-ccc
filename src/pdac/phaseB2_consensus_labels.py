"""Phase B2 (v2): POST-HOC consensus_T annotation (addendum in docs/PDAC_V2_PROTEIN_RULE.md) + protein triangulation.

consensus_T = author_mapped labels, except: T_cell only where author T (CD8+ T / CD4+ T / Treg) AND v1 marker T_cell;
author T without marker confirmation -> 'T_cell_unconfirmed'.
Triangulation (independent third source): CD3 / CD45 protein in consensus-T vs author-only T vs marker-only T vs rest,
as within-section AUC (consensus-T vs all non-T cells) and median log1p intensity.
Reads the per-cell tables written by phaseB1 (results/pdac/phaseB/sections/*_labels.csv.gz).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, out_dir, register, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

PHASE = "phaseB"
AUTHOR_T = {"CD8+ T", "CD4+ T", "Treg"}


def consensus(c):
    at = c.author_majortype.isin(AUTHOR_T).to_numpy()
    mt = (c.marker == "T_cell").to_numpy()
    lab = c.author_mapped.astype(str).to_numpy().copy()
    lab[at & mt] = "T_cell"
    lab[at & ~mt] = "T_cell_unconfirmed"
    grp = np.select([at & mt, at & ~mt, ~at & mt], ["consensus_T", "author_only_T", "marker_only_T"], "other")
    return lab, grp


def main() -> None:
    log = start_run("phaseB2_consensus_labels")
    files = sorted((RESULTS_PDAC / "phaseB" / "sections").glob("*_labels.csv.gz"))
    tri, counts = [], []
    d = out_dir("phaseB/sections")
    allc = []
    for f in files:
        c = pd.read_csv(f)
        s = c.section.iloc[0]
        lab, grp = consensus(c)
        c["consensus_T"] = lab
        c["T_group"] = grp
        out = c[["section", "cell_id", "consensus_T", "T_group"]]
        p = d / f"{s}_consensus_labels.csv.gz"
        out.to_csv(p, index=False)
        register(p, PHASE, f"Per-cell POST-HOC consensus_T labels and T-group (section {s})")
        for g in ("consensus_T", "author_only_T", "marker_only_T", "other"):
            m = grp == g
            row = {"section": s, "group": g, "n": int(m.sum())}
            for ch in ("CD3", "CD45", "PanCK"):
                x = c[f"log1p_{ch}"].to_numpy()
                row[f"median_log1p_{ch}"] = float(np.median(x[m])) if m.any() else np.nan
                if g != "other" and m.any():
                    y = np.r_[np.ones(m.sum()), np.zeros((grp == "other").sum())]
                    row[f"AUC_{ch}_vs_other"] = float(roc_auc_score(y, np.r_[x[m], x[grp == "other"]]))
            tri.append(row)
        counts.append(pd.Series(lab).value_counts().rename("n").reset_index().rename(columns={"index": "label"}).assign(section=s))
        allc.append(c[["section", "T_group", "log1p_CD3", "log1p_CD45"]])
        log.info("section %s: consensus T=%d, author-only T=%d, marker-only T=%d", s, int((grp == "consensus_T").sum()),
                 int((grp == "author_only_T").sum()), int((grp == "marker_only_T").sum()))
    tri = pd.DataFrame(tri)
    save_csv(tri, PHASE, "consensus_T_protein_triangulation.csv",
             "Per section x T-group: n, median log1p CD3/CD45/PanCK, AUC vs non-T cells (independent protein check)", log)
    save_csv(pd.concat(counts, ignore_index=True), PHASE, "consensus_T_label_counts.csv", "consensus_T label counts per section", log)
    med = tri.groupby("group")[[c for c in tri.columns if c.startswith("AUC_") or c.startswith("median_")]].median()
    tot = tri.groupby("group").n.sum()
    med.insert(0, "n_total", tot)
    save_csv(med, PHASE, "consensus_T_protein_triangulation_median.csv", "Median over 16 sections of triangulation stats", log, index=True)
    log.info("triangulation medians over sections:\n%s", med.round(3).to_string())
    u5 = tri[tri.section == "U5-a"].set_index("group")
    save_json({"n_by_group_all_sections": tot.astype(int).to_dict(), "U5-a_n_by_group": u5.n.astype(int).to_dict(),
               "median_over_sections": med.round(4).to_dict(orient="index")}, PHASE, "consensus_T_summary.json",
              "consensus_T counts and protein triangulation summary", log)

    cc = pd.concat(allc, ignore_index=True)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    cols = {"consensus_T": "#d62728", "author_only_T": "#ff7f0e", "marker_only_T": "#9467bd", "other": "#9e9e9e"}
    for ax, ch in zip(axes[:2], ("CD3", "CD45")):
        # within-section centring so slides with different staining are comparable
        x = cc[f"log1p_{ch}"] - cc.groupby("section")[f"log1p_{ch}"].transform("median")
        for g, col in cols.items():
            ax.hist(x[cc.T_group == g], bins=120, range=(-4, 5), density=True, histtype="step", color=col, lw=1.4,
                    label=f"{g} (n={int((cc.T_group == g).sum())})")
        ax.set_xlabel(f"log1p(Mean.{ch}) minus section median")
        ax.legend(fontsize=7)
        ax.set_title(f"{ch} protein by T-group, 16 sections pooled", fontsize=10)
    ax = axes[2]
    for k, g in enumerate(("consensus_T", "author_only_T", "marker_only_T")):
        v = tri[tri.group == g]
        ax.scatter(v["AUC_CD3_vs_other"], v["AUC_CD45_vs_other"], c=cols[g], label=g, s=28)
    ax.axhline(0.5, color="k", lw=0.5)
    ax.axvline(0.5, color="k", lw=0.5)
    ax.set_xlabel("AUC CD3 (group vs non-T cells)")
    ax.set_ylabel("AUC CD45 (group vs non-T cells)")
    ax.legend(fontsize=7)
    ax.set_title("Per-section protein AUC (one point per section)", fontsize=10)
    fig.suptitle("Phase B2 (post-hoc): consensus T cells (author ∧ marker) — independent protein check")
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseB2_consensus_T_protein_triangulation", "Protein intensity by T-group; per-section AUC", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

"""Phase D2 (v2): tumour-level view of the CXCL12->CXCR4 CAF->T_cell headline (method declared before running).

Unit = tumour, not section. 16 sections -> 13 tumours: a section's tumour is its donor id with the suffixes
'_2nd' / '_B' removed (T4-a/T4-b -> 2762, U4-a/U4-b -> 011, U6-a_1/U6-a_2 -> 2626). This grouping is an assumption
consistent with the dataset's stated 13 tumours (results/pdac/dataset_vetting.csv); it is not otherwise verified.

Per annotation (marker, author_mapped, consensus_T; protein_gmm / protein_mad reported but they failed validation)
and per score (edge E, matched M, ratio Q):
  * tumour p = Fisher combination of the section add-one p-values within the tumour (sections are independent
    label-permutation tests); single-section tumours keep their p. Discrete add-one p makes Fisher conservative.
  * BH across tumours (within annotation x score); tumour significant if q < 0.05.
  * tumour mean z and tumour mean observed Q.
  * across tumours (exploratory): one-sided Wilcoxon signed-rank test of tumour-mean z > 0, and the number of
    tumours with mean observed Q > 1.
Sections where the annotation has no T_cell class (or 0 edges) are skipped; tumours with no evaluable section drop out.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import combine_pvalues, wilcoxon  # noqa: E402

import ccc_audit as ca  # noqa: E402

PHASE = "phaseD"
ANN = ["marker", "author_mapped", "consensus_T", "protein_gmm", "protein_mad"]
SCORES = {"edge": "E", "matched": "M", "ratio": "Q"}


def main() -> None:
    log = start_run("phaseD2_tumour_level")
    info = pd.read_csv(RESULTS_PDAC / "phase3" / "cross_section_section_info.csv")
    info["tumour"] = info.donor.astype(str).str.replace(r"(_2nd|_B)$", "", regex=True)
    tmap = dict(zip(info.section, info.tumour))
    log.info("section -> tumour: %s; %d tumours", tmap, info.tumour.nunique())
    hl = pd.read_csv(RESULTS_PDAC / "phaseB" / "headline_CAF_T_by_annotation_section.csv")
    hl = hl[hl.evaluable].copy()
    hl["tumour"] = hl.section.map(tmap)
    rows = []
    for a in ANN:
        d = hl[hl.annotation == a]
        for sc, lab in SCORES.items():
            for t, g in d.groupby("tumour"):
                p = g[f"{sc}_p"].to_numpy()
                pt = float(combine_pvalues(p, method="fisher").pvalue) if len(p) > 1 else float(p[0])
                rows.append({"annotation": a, "score": lab, "tumour": t, "n_sections": len(g),
                             "sections": ";".join(g.section), "treatment_status": g.treatment_status.iloc[0],
                             "tumour_p_fisher": pt, "tumour_mean_z": float(g[f"{sc}_z"].mean()),
                             "tumour_mean_obs": float(g[f"{sc}_obs"].mean())})
    tt = pd.DataFrame(rows)
    tt["tumour_q_bh"] = tt.groupby(["annotation", "score"]).tumour_p_fisher.transform(
        lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    tt["sig_uncorrected"] = tt.tumour_p_fisher < 0.05
    tt["sig_bh"] = tt.tumour_q_bh < 0.05
    save_csv(tt, PHASE, "tumour_level_headline_CAF_T.csv",
             "CAF->T CXCL12-CXCR4 per tumour (13): Fisher-combined section p, BH across tumours, mean z / obs", log)
    summ = []
    for (a, s), g in tt.groupby(["annotation", "score"]):
        z = g.tumour_mean_z.to_numpy()
        w = wilcoxon(z, alternative="greater") if len(z) >= 5 else None
        summ.append({"annotation": a, "score": s, "n_tumours": len(g), "n_sig_uncorrected": int(g.sig_uncorrected.sum()),
                     "n_sig_bh": int(g.sig_bh.sum()), "min_q_bh": float(g.tumour_q_bh.min()),
                     "n_mean_z_pos": int((z > 0).sum()), "median_tumour_mean_z": float(np.median(z)),
                     "wilcoxon_z_gt0_p": float(w.pvalue) if w is not None else np.nan,
                     "n_mean_obs_gt1": int((g.tumour_mean_obs > 1).sum()) if s == "Q" else np.nan})
    summ = pd.DataFrame(summ)
    summ["annotation"] = pd.Categorical(summ.annotation, ANN, ordered=True)
    summ = summ.sort_values(["annotation", "score"]).reset_index(drop=True)
    save_csv(summ, PHASE, "tumour_level_headline_summary.csv",
             "Tumour-level summary per annotation x score: tumours significant (Fisher, BH), Wilcoxon of tumour z > 0", log)
    for r in summ.itertuples():
        log.info("TUMOUR %s %s: %d tumours; sig %d (BH %d, min q %.3f); z>0 in %d; median z %.2f; Wilcoxon p %.4f; obs>1 %s",
                 r.annotation, r.score, r.n_tumours, r.n_sig_uncorrected, r.n_sig_bh, r.min_q_bh, r.n_mean_z_pos,
                 r.median_tumour_mean_z, r.wilcoxon_z_gt0_p, r.n_mean_obs_gt1)
    save_json({"section_to_tumour": tmap, "summary": summ.astype({"annotation": str}).to_dict(orient="records")}, PHASE,
              "tumour_level_summary.json", "Tumour-level headline summary", log)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=False)
    cols = {"marker": "#7f7f7f", "author_mapped": "#1f77b4", "consensus_T": "#d62728", "protein_gmm": "#ff7f0e",
            "protein_mad": "#ffbb78"}
    tumours = sorted(tt.tumour.unique())
    xi = {t: i for i, t in enumerate(tumours)}
    for ax, s in zip(axes, ("E", "M", "Q")):
        for k, a in enumerate(ANN):
            g = tt[(tt.annotation == a) & (tt.score == s)]
            ax.scatter([xi[t] + (k - 2) * 0.13 for t in g.tumour], g.tumour_mean_z, c=cols[a], s=np.where(g.sig_bh, 55, 16),
                       edgecolors=np.where(g.sig_bh, "k", "none"), label=a)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(tumours)))
        ax.set_xticklabels(tumours, rotation=90, fontsize=7)
        ax.set_title(f"{s}: tumour-mean z (large/black edge = Fisher+BH q<.05 across tumours)", fontsize=8)
        if s != "Q":
            ax.set_yscale("symlog", linthresh=2)
    axes[0].legend(fontsize=7)
    axes[0].set_ylabel("tumour-mean z")
    fig.suptitle("CXCL12→CXCR4 CAF→T_cell at tumour level (13 tumours; audit/diagnostic only)")
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseD2_tumour_level_headline", "Tumour-level E/M/Q for the CAF->T headline, 5 annotations", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

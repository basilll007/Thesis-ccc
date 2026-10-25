"""Phase A (v2): spot-check the five headline numbers of docs/PDAC_SUMMARY.md against row-level result files.

Each number is RECOMPUTED from the most granular saved table (not copied from a summary file), then compared
with the value claimed in the summary. Output: results/pdac/phaseA/spotcheck.csv (+ figure, run log).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, save_csv, save_figure, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

PHASE = "phaseA"
ALPHA = 0.05
K4 = ["sender", "receiver", "ligand", "receptor"]


def main() -> None:
    log = start_run("phaseA_spotcheck")
    rows = []

    # (1) breast gate: 3/9 C/E disagreements explained by score form -----------------------------
    src1 = RESULTS_PDAC / "phase1" / "breast_matched_baseline_paper22.csv"
    p = pd.read_csv(src1)
    ev = p[p.edge_status == "evaluated"]
    c, e, m = ev.comp_p_plus1 < ALPHA, ev.edge_p_plus1 < ALPHA, ev.matched_p_plus1 < ALPHA
    dis = c != e
    n_form = int((dis & (m == e)).sum())
    rec1 = f"{n_form}/{int(dis.sum())}"
    rows.append({"item": 1, "claim": "breast gate: matched baseline explains C/E disagreements", "claimed_value": "3/9",
                 "recomputed_value": rec1, "source_file": src1.relative_to(RESULTS_PDAC.parents[1]).as_posix(),
                 "method": "rows with edge evaluated; disagreement = (comp_p_plus1<.05) != (edge_p_plus1<.05); "
                           "score-form if matched verdict == edge verdict"})

    # (2)/(3) PDAC CAF->T CXCL12->CXCR4, 16 sections, marker and author labels ----------------------
    src_m = RESULTS_PDAC / "phase3" / "cross_section_6NN_all_rows.csv"
    src_a = RESULTS_PDAC / "phase3" / "author_annotation_6NN_all_sections.csv"
    xm = pd.read_csv(src_m)
    xa = pd.read_csv(src_a)
    sel = lambda d: d[(d.ligand == "CXCL12") & (d.receptor == "CXCR4") & (d.sender == "CAF_fibroblast") & (d.receiver == "T_cell")]
    hm = sel(xm[xm.annotation == "cell_type_marker"])
    ha = sel(xa)
    n_m, n_a = hm.section.nunique(), ha.section.nunique()
    e_m, e_a = int((hm.edge_p_plus1 < ALPHA).sum()), int((ha.edge_p_plus1 < ALPHA).sum())
    q_m, q_a = int((hm.ratio_p_plus1 < ALPHA).sum()), int((ha.ratio_p_plus1 < ALPHA).sum())
    qb_m = int((hm.ratio_q_bh_within_section_annotation < ALPHA).sum())
    qb_a = int((ha.ratio_q_bh_within_section < ALPHA).sum())
    rel = lambda f: f.relative_to(RESULTS_PDAC.parents[1]).as_posix()
    rows.append({"item": 2, "claim": "PDAC CAF->T CXCL12-CXCR4 edge-significant sections (marker / author)",
                 "claimed_value": "16/16 / 16/16", "recomputed_value": f"{e_m}/{n_m} / {e_a}/{n_a}",
                 "source_file": f"{rel(src_m)}; {rel(src_a)}", "method": "count sections with edge_p_plus1 < 0.05"})
    rows.append({"item": 3, "claim": "PDAC CAF->T Q-significant sections (marker / author); after BH within section",
                 "claimed_value": "4/16 / 2/16; 0/16 after BH (both annotations)",
                 "recomputed_value": f"{q_m}/{n_m} / {q_a}/{n_a}; {qb_m}/{n_m} / {qb_a}/{n_a} after BH",
                 "source_file": f"{rel(src_m)}; {rel(src_a)}",
                 "method": "ratio_p_plus1 < 0.05; BH q (within section[/annotation]) < 0.05"})

    # (4) background flag: 24/31 axes -----------------------------------------------------------
    src4 = RESULTS_PDAC / "phase3" / "background_check_detection_vs_negprobes.csv"
    src4b = RESULTS_PDAC / "phase3" / "background_flag_per_axis.csv"
    bg = pd.read_csv(src4)
    nsec = bg.groupby("gene").near_background.sum()
    axes = pd.read_csv(src4b)
    axes = axes[axes.testable]
    flag = axes.apply(lambda r: max(nsec.get(r.ligand, 0), nsec.get(r.receptor, 0)) >= 8, axis=1)
    rows.append({"item": 4, "claim": "PDAC axes flagged background-level (either gene near background in >= 8/16 sections)",
                 "claimed_value": "24/31", "recomputed_value": f"{int(flag.sum())}/{len(axes)}",
                 "source_file": f"{rel(src4)} (per section x gene); axis list from {rel(src4b)}",
                 "method": "recount near_background per gene across sections from the per-section table; flag axis if max >= 8"})
    agree_flag = bool((flag.to_numpy() == axes.background_level_gene_diagnostic.to_numpy()).all())
    log.info("item 4: recomputed per-axis flags identical to stored column: %s", agree_flag)

    # (5) Spearman edge z vs ligrec mean, crop ---------------------------------------------------
    src5a = RESULTS_PDAC / "phase4" / "ours_CEMQ_6NN_crop.csv"
    src5b = RESULTS_PDAC / "phase4" / "squidpy_ligrec_crop.csv"
    o = pd.read_csv(src5a)
    lr = pd.read_csv(src5b)
    mm = o[K4 + ["edge_z"]].merge(lr, on=K4, how="left")
    ok = np.isfinite(mm.edge_z) & np.isfinite(mm.ligrec_mean)
    rho = float(spearmanr(mm.edge_z[ok], mm.ligrec_mean[ok]).statistic)
    rows.append({"item": 5, "claim": "Spearman(edge z, ligrec mean), cluster level, crop", "claimed_value": "0.02",
                 "recomputed_value": f"{rho:.3f} (n={int(ok.sum())})", "source_file": f"{rel(src5a)}; {rel(src5b)}",
                 "method": "merge on (sender, receiver, ligand, receptor); Spearman on finite rows"})

    df = pd.DataFrame(rows)
    # numeric comparison: claimed and recomputed numbers (in order of appearance) must be identical
    # (item 5: equal after rounding to the 2 decimals reported in the summary)
    claimed_nums = {1: [3, 9], 2: [16, 16, 16, 16], 3: [4, 16, 2, 16, 0, 16, 0, 16], 4: [24, 31], 5: [0.02]}
    recomputed_nums = {1: [n_form, int(dis.sum())], 2: [e_m, n_m, e_a, n_a], 3: [q_m, n_m, q_a, n_a, qb_m, n_m, qb_a, n_a],
                       4: [int(flag.sum()), len(axes)], 5: [round(rho, 2)]}
    df["match"] = ["yes" if claimed_nums[i] == recomputed_nums[i] else "no" for i in df["item"]]
    for r in df.itertuples():
        log.info("spot-check %d: claimed=%s recomputed=%s match=%s", r.item, r.claimed_value, r.recomputed_value, r.match)
    save_csv(df, PHASE, "spotcheck.csv", "Five headline numbers of PDAC_SUMMARY.md: claimed vs recomputed, source, match", log)

    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.axis("off")
    tab = ax.table(cellText=df[["item", "claimed_value", "recomputed_value", "match"]].values,
                   colLabels=["item", "claimed", "recomputed", "match"], loc="center", cellLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1, 1.5)
    for i, mt in enumerate(df["match"], start=1):
        tab[i, 3].set_facecolor("#c7e9c0" if mt == "yes" else "#fcbba1")
    ax.set_title("Phase A: headline-number spot-check (PDAC_SUMMARY.md vs recomputation from row-level files)", fontsize=9)
    save_figure(fig, PHASE, "phaseA_spotcheck_table", "Spot-check table of the five headline numbers", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

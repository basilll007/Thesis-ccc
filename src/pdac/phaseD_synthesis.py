"""Phase D (v2): synthesis table + figure for docs/PDAC_SUMMARY_v2.md.

Collects every number quoted in the v2 summary from the Phase A/B result files into one table
(results/pdac/phaseD/v2_key_numbers.csv: item, value, source file) and draws one synthesis figure:
  left  : CAF->T CXCL12-CXCR4, sections significant under E / Q / Q-after-BH per annotation
  middle: per-section z_M vs z_E for consensus_T (is E carried by the cell-type product?)
  right : share of edge-significant rows that are M-significant / Q-significant / Q-significant after BH
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, save_csv, save_figure, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phaseD"
B = RESULTS_PDAC / "phaseB"
A = RESULTS_PDAC / "phaseA"
ANN = ["marker", "author_mapped", "protein_gmm", "protein_mad", "consensus_T"]


def r(p):
    return p.relative_to(RESULTS_PDAC.parents[1]).as_posix()


def main() -> None:
    log = start_run("phaseD_synthesis")
    rows = []

    def add(item, value, src):
        rows.append({"item": item, "value": value, "source_file": r(src)})

    sp = pd.read_csv(A / "spotcheck.csv")
    add("spot-check: headline numbers matching", f"{int((sp.match == 'yes').sum())}/{len(sp)}", A / "spotcheck.csv")
    fs = pd.read_csv(A / "pvalue_conventions_flip_summary.csv")
    br = fs[fs.dataset.str.startswith("breast")]
    pd_ = fs[fs.dataset.str.startswith("pdac")]
    add("breast: uncorrected verdict flips legacy vs add-one", f"{int(br.n_flip_uncorrected.sum())} of {int(br.n_rows.sum())} rows",
        A / "pvalue_conventions_flip_summary.csv")
    add("breast: BH verdict flips legacy vs add-one", int(br.n_flip_bh.sum()), A / "pvalue_conventions_flip_summary.csv")
    add("PDAC: uncorrected verdict flips (all ratio, n_null<500)", f"{int(pd_.n_flip_uncorrected.sum())} of {int(pd_.n_rows.sum())}",
        A / "pvalue_conventions_flip_summary.csv")
    add("PDAC: BH verdict flips", int(pd_.n_flip_bh.sum()), A / "pvalue_conventions_flip_summary.csv")
    pj = json.loads((B / "protein_labels_summary.json").read_text())
    for sc in ("U5-a", "all_16_sections"):
        for a in ("protein_gmm", "marker"):
            v = pj["T_cell_vs_author_T"][sc][a]
            add(f"T_cell precision/recall vs author T, {a}, {sc}", f"{v['precision']:.3f} / {v['recall']:.3f} (tp {v['tp']}/{v['n_pred']})",
                B / "protein_labels_summary.json")
    add("GMM accepted (Ashman D>=2) sections per channel", json.dumps(pj["gmm_accepted_by_channel"]), B / "protein_labels_summary.json")
    tri = pd.read_csv(B / "consensus_T_protein_triangulation_median.csv", index_col=0)
    for g in ("consensus_T", "author_only_T", "marker_only_T"):
        add(f"protein AUC vs non-T (median over sections), {g}: CD3 / CD45",
            f"{tri.loc[g, 'AUC_CD3_vs_other']:.3f} / {tri.loc[g, 'AUC_CD45_vs_other']:.3f} (n={int(tri.loc[g, 'n_total'])})",
            B / "consensus_T_protein_triangulation_median.csv")
    hs = pd.read_csv(B / "headline_CAF_T_summary_by_annotation.csv").set_index("annotation")
    for a in ANN:
        h = hs.loc[a]
        add(f"CAF->T headline, {a}: evaluable; E sig (BH); Q sig (BH); median z_E, z_Q, Q",
            f"{int(h.n_sections_evaluable)}/16; {int(h.n_edge_sig)} ({int(h.n_edge_sig_bh)}); {int(h.n_ratio_sig)} ({int(h.n_ratio_sig_bh)}); "
            f"{h.median_edge_z:.2f}, {h.median_ratio_z:.2f}, {h.median_ratio_obs:.3f}", B / "headline_CAF_T_summary_by_annotation.csv")
    hl = pd.read_csv(B / "headline_CAF_T_by_annotation_section.csv")
    c = hl[hl.annotation == "consensus_T"]
    add("consensus_T: sections with z_M > z_E", f"{int((c.matched_z > c.edge_z).sum())}/{len(c)}", B / "headline_CAF_T_by_annotation_section.csv")
    add("consensus_T: sections with observed Q > 1", f"{int((c.ratio_obs > 1).sum())}/{len(c)}", B / "headline_CAF_T_by_annotation_section.csv")
    add("consensus_T: smallest Q BH q over sections", f"{c.ratio_q.min():.3f}", B / "headline_CAF_T_by_annotation_section.csv")
    mi = pd.read_csv(B / "edge_vs_product_mirroring_by_annotation.csv")
    for a in ("marker", "author_mapped", "consensus_T"):
        for ax_ in ("all_axes", "unflagged_axes_only"):
            m = mi[(mi.annotation == a) & (mi["axes"] == ax_)].iloc[0]
            add(f"mirroring {a} {ax_}: rho(zE,zM); frac E-sig also M-sig / Q-sig / Q-sig BH",
                f"{m.spearman_edge_z_vs_matched_z:.3f}; {m.frac_edge_sig_also_matched_sig:.3f} / {m.frac_edge_sig_also_ratio_sig:.3f} / "
                f"{m.frac_edge_sig_also_ratio_sig_bh:.3f} (n E-sig {int(m.n_edge_sig)})", B / "edge_vs_product_mirroring_by_annotation.csv")
    kn = pd.DataFrame(rows)
    save_csv(kn, PHASE, "v2_key_numbers.csv", "Every number quoted in docs/PDAC_SUMMARY_v2.md with its source file", log)
    for x in rows:
        log.info("%s = %s  [%s]", x["item"], x["value"], x["source_file"])

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw={"width_ratios": [1.5, 1, 1.2]})
    ax = axes[0]
    x = np.arange(len(ANN))
    for k, (col, lab, colr) in enumerate((("n_edge_sig", "E sig", "#d62728"), ("n_ratio_sig", "Q sig", "#2ca02c"),
                                          ("n_ratio_sig_bh", "Q sig after BH", "#006d2c"))):
        ax.bar(x + (k - 1) * 0.27, hs.loc[ANN, col], width=0.27, color=colr, label=lab)
    ax.plot(x, hs.loc[ANN, "n_sections_evaluable"], "k_", ms=30, label="sections evaluable")
    ax.set_xticks(x)
    ax.set_xticklabels([a + ("\n(failed validation)" if a.startswith("protein") else "\n(post-hoc)" if a == "consensus_T" else "")
                        for a in ANN], fontsize=8)
    ax.set_ylabel("sections (of 16)")
    ax.set_ylim(0, 17)
    ax.legend(fontsize=8)
    ax.set_title("CXCL12→CXCR4 CAF→T_cell: sections significant (add-one p<.05; BH within section)", fontsize=9)
    ax = axes[1]
    ax.scatter(c.matched_z, c.edge_z, c="#d62728")
    for s, xm, ye in zip(c.section, c.matched_z, c.edge_z):
        ax.annotate(s, (xm, ye), fontsize=6)
    lim = [0, max(c.matched_z.max(), c.edge_z.max()) * 1.05]
    ax.plot(lim, lim, "k--", lw=0.7)
    ax.set_xlabel("z_M (matched product L̄s·R̄r, no graph)")
    ax.set_ylabel("z_E (6-NN edge score)")
    ax.set_title("consensus_T, 16 sections: E vs matched product", fontsize=9)
    ax = axes[2]
    mm = mi[mi.annotation.isin(["marker", "author_mapped", "consensus_T"])]
    lab = [f"{a}\n{t.replace('_axes', '').replace('_only', '')}" for a, t in zip(mm.annotation, mm["axes"])]
    y = np.arange(len(mm))
    ax.barh(y - 0.27, mm.frac_edge_sig_also_matched_sig, height=0.27, color="#1f77b4", label="also M-sig")
    ax.barh(y, mm.frac_edge_sig_also_ratio_sig, height=0.27, color="#2ca02c", label="also Q-sig")
    ax.barh(y + 0.27, mm.frac_edge_sig_also_ratio_sig_bh, height=0.27, color="#006d2c", label="also Q-sig after BH")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of edge-significant section-rows")
    ax.legend(fontsize=7)
    ax.set_title("Full catalogue, 16 sections", fontsize=9)
    fig.suptitle("PDAC v2 synthesis (audit/diagnostic only)")
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseD_v2_synthesis", "v2 synthesis: headline per annotation; E vs M; E-sig overlap with M/Q", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

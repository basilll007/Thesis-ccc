"""Phase 3c: summarise the saved Phase 3 tables (no new permutations; low memory).

Reads the CSVs written by phase3_pdac_audit.py (primary section U5-a; the run was stopped by the host at the
200 um radius step, after the 6-NN, annotation and k analyses were saved) and phase3b_cross_section.py.
Produces the headline table, multiple-testing counts, C/E decomposition (score-form vs proximity),
sparse-gene list, k-sensitivity summary, cross-section counts, and the figures the stopped run did not write.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import finish_run, out_dir, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

PHASE = "phase3"
TAG = "mendeley_pdac_smi__U5-a"
HEAD = [("CXCL12", "CXCR4", "CAF_fibroblast", "T_cell"), ("CXCL12", "CXCR4", "CAF_fibroblast", "epithelial_ductal_tumor")]
CHECKPOINT = {("CD274", "PDCD1"), ("PDCD1LG2", "PDCD1"), ("CD86", "CTLA4"), ("CD80", "CTLA4"), ("CD80", "CD28"),
              ("CD86", "CD28"), ("LGALS9", "HAVCR2")}


def decomp(r):
    if r["edge_status"] != "evaluated":
        return "edge not evaluable"
    if r["comp_sig"] == r["edge_sig"]:
        return "C and E agree"
    return "score-form" if r["matched_sig"] == r["edge_sig"] else "proximity"


def is_head(r):
    return (r["ligand"], r["receptor"], r["sender"], r["receiver"]) in HEAD


def main() -> None:
    log = start_run("phase3c_summarise")
    d = out_dir(PHASE)
    main_t = pd.read_csv(d / f"{TAG}_6NN_scores_both_annotations.csv")
    cov = pd.read_csv(d / f"{TAG}_lr_coverage.csv")
    det = pd.read_csv(d / f"{TAG}_gene_detection_counts.csv")
    ks = pd.read_csv(d / f"{TAG}_k_sensitivity.csv")
    ann = pd.read_csv(d / f"{TAG}_annotation_sensitivity_summary.csv")
    xs = pd.read_csv(d / "cross_section_summary_by_interaction.csv")
    xall = pd.read_csv(d / "cross_section_6NN_all_rows.csv",
                       usecols=["section", "annotation", "ligand", "receptor", "sender", "receiver", "edge_status",
                                "comp_sig", "matched_sig", "edge_sig", "ratio_sig", "sparse_gene_diagnostic"])

    main_t["decomposition"] = main_t.apply(decomp, axis=1)
    main_t["checkpoint_axis"] = [(l, r) in CHECKPOINT for l, r in zip(main_t.ligand, main_t.receptor)]
    out = {"primary_section": "U5-a (largest untreated section; pre-declared)"}
    for a, g in main_t.groupby("annotation"):
        ev = g[g.edge_status == "evaluated"]
        dis = ev[ev.comp_sig != ev.edge_sig]
        es = ev[ev.edge_sig]
        out[a] = {
            "n_rows_tested": int(len(g)), "n_lr_pairs": int(g[["ligand", "receptor"]].drop_duplicates().shape[0]),
            "n_cell_types": int(g["sender"].nunique()),
            "n_edge_evaluable": int(len(ev)),
            **{f"n_{k}_sig_uncorrected": int(g[f"{k}_sig"].sum()) for k in ("comp", "matched", "edge", "ratio")},
            **{f"n_{k}_q_lt_0.05_bh_within_annotation": int((g[f"{k}_q_bh_within_annotation"] < 0.05).sum())
               for k in ("comp", "matched", "edge", "ratio")},
            **{f"n_{k}_q_lt_0.05_bh_all_rows": int((g[f"{k}_q_bh_all_rows"] < 0.05).sum()) for k in ("comp", "matched", "edge", "ratio")},
            "n_C_E_disagree": int(len(dis)), "n_disagree_score_form": int((dis.decomposition == "score-form").sum()),
            "n_disagree_proximity": int((dis.decomposition == "proximity").sum()),
            "n_disagree_sparse_gene": int(dis.sparse_gene_diagnostic.sum()),
            "n_edge_sig": int(len(es)), "n_edge_sig_and_matched_sig": int(es.matched_sig.sum()),
            "n_edge_sig_and_ratio_sig": int(es.ratio_sig.sum()),
            "n_rows_sparse_gene": int(g.sparse_gene_diagnostic.sum()),
            "spearman_z_comp_vs_edge": float(spearmanr(ev.comp_z, ev.edge_z).statistic),
            "spearman_z_matched_vs_edge": float(spearmanr(ev.matched_z, ev.edge_z).statistic),
            "checkpoint_rows": int(g.checkpoint_axis.sum()),
            "checkpoint_rows_edge_sig": int((g.checkpoint_axis & g.edge_sig).sum()),
            "checkpoint_rows_comp_sig": int((g.checkpoint_axis & g.comp_sig).sum()),
        }
        log.info("[%s] %s", a, out[a])
    out["n_rows_tested_total_both_annotations"] = int(len(main_t))
    head = main_t[main_t.apply(is_head, axis=1)]
    save_csv(head, PHASE, f"{TAG}_headline_summary_rows.csv", "Headline CXCL12-CXCR4 rows (primary section) with decomposition", log)
    out["headline"] = head[["annotation", "sender", "receiver", "n_qualifying_edges", "comp_z", "matched_z", "edge_z", "ratio_obs",
                            "ratio_z", "comp_p_plus1", "edge_p_plus1", "ratio_p_plus1", "edge_q_bh_all_rows", "ratio_q_bh_all_rows",
                            "comp_sig", "matched_sig", "edge_sig", "ratio_sig", "decomposition"]].to_dict(orient="records")

    # checkpoint axes table
    cp = main_t[main_t.checkpoint_axis & (main_t.annotation == "cell_type_marker")]
    save_csv(cp, PHASE, f"{TAG}_checkpoint_axes_marker.csv", "Checkpoint-axis rows (marker labels, 6-NN) incl. sparse flags", log)

    # sparse genes (receptors and ligands < 20 cells)
    tested_genes = set(cov[cov.testable].ligand) | set(cov[cov.testable].receptor)
    sp_ = det[det.gene.isin(tested_genes) & (det.n_cells_detected < 20)][["gene", "role", "n_cells_detected"]]
    save_csv(sp_, PHASE, f"{TAG}_sparse_genes.csv", "Tested genes detected in fewer than 20 cells (sparse-gene diagnostic)", log)
    rec_det = det[det.role.str.contains("receptor")][["gene", "on_panel", "n_cells_detected", "n_cells_total"]]
    save_csv(rec_det, PHASE, f"{TAG}_receptor_detection_counts.csv", "Detection counts for every catalogue receptor", log)
    out["coverage"] = {"n_catalogue_pairs": int(len(cov)), "n_testable": int(cov.testable.sum()),
                       "absent_pairs": cov.loc[~cov.testable, ["ligand", "receptor"]].apply(lambda r: f"{r.ligand}-{r.receptor}", axis=1).tolist(),
                       "testable_pairs_sparse": cov.loc[cov.sparse_gene_diagnostic, ["ligand", "receptor"]].apply(lambda r: f"{r.ligand}-{r.receptor}", axis=1).tolist(),
                       "sparse_genes": dict(zip(sp_.gene, sp_.n_cells_detected.astype(int)))}

    # k sensitivity
    kk = []
    base = ks[ks.k == 6].set_index(["ligand", "receptor", "sender", "receiver"])
    for k, g in ks.groupby("k"):
        g = g.set_index(["ligand", "receptor", "sender", "receiver"])
        j = base.join(g, rsuffix="_k", how="inner")
        ok = np.isfinite(j.edge_z) & np.isfinite(j.edge_z_k)
        kk.append({"k": int(k), "n_edge_sig": int(g.edge_sig.sum()), "n_ratio_sig": int(g.ratio_sig.sum()),
                   "spearman_edge_z_vs_k6": float(spearmanr(j.edge_z[ok], j.edge_z_k[ok]).statistic),
                   "n_edge_verdict_changes_vs_k6": int((j.edge_sig != j.edge_sig_k).sum()),
                   "n_ratio_verdict_changes_vs_k6": int((j.ratio_sig != j.ratio_sig_k).sum())})
    ksum = pd.DataFrame(kk)
    save_csv(ksum, PHASE, f"{TAG}_k_sensitivity_summary.csv", "k in {4,6,8,12}: significant counts and agreement with k=6", log)
    khead = ks[ks.apply(is_head, axis=1)][["k", "sender", "receiver", "n_qualifying_edges", "edge_z", "ratio_z", "edge_sig", "ratio_sig"]]
    save_csv(khead, PHASE, f"{TAG}_k_sensitivity_headline.csv", "Headline rows across k", log)
    out["k_sensitivity"] = ksum.to_dict(orient="records")
    out["annotation_sensitivity"] = ann.to_dict(orient="records")

    # cross-section
    xm = xs[xs.annotation == "cell_type_marker"]
    xh = xs[xs.apply(is_head, axis=1)]
    out["cross_section_headline"] = xh.to_dict(orient="records")
    xe = xall[(xall.annotation == "cell_type_marker") & (xall.edge_status == "evaluated")]
    out["cross_section_marker_all_rows"] = {
        "n_section_rows": int(len(xe)), "n_edge_sig": int(xe.edge_sig.sum()),
        "n_edge_sig_and_matched_sig": int((xe.edge_sig & xe.matched_sig).sum()),
        "n_edge_sig_and_ratio_sig": int((xe.edge_sig & xe.ratio_sig).sum()),
        "n_C_E_disagree": int((xe.comp_sig != xe.edge_sig).sum()),
        "n_disagree_score_form": int(((xe.comp_sig != xe.edge_sig) & (xe.matched_sig == xe.edge_sig)).sum()),
        "n_disagree_proximity": int(((xe.comp_sig != xe.edge_sig) & (xe.matched_sig == xe.comp_sig)).sum()),
        "n_interactions_edge_sig_in_all_16": int((xm.n_edge_sig == 16).sum()),
        "n_interactions_ratio_sig_in_ge_12_of_16": int((xm.n_ratio_sig >= 12).sum()),
    }
    rob = xm[(xm.n_ratio_sig >= 12)].sort_values("n_ratio_sig", ascending=False)
    save_csv(rob, PHASE, "cross_section_ratio_sig_in_ge12_sections_marker.csv",
             "Interactions whose edge/matched ratio is significant in >= 12 of 16 sections (marker labels)", log)
    save_json(out, PHASE, f"{TAG}_phase3_summary.json", "Phase 3 key numbers (primary section + cross-section)", log)

    # figures --------------------------------------------------------------------------------
    for a, g in main_t.groupby("annotation"):
        g = g[np.isfinite(g.comp_z)]
        top = g.reindex(g.comp_z.abs().sort_values(ascending=False).index).head(40).sort_values("comp_z")
        y = np.arange(len(top))
        fig, ax = plt.subplots(figsize=(12, 0.3 * len(top) + 2))
        for key, c, mk in (("comp", "#7f7f7f", "o"), ("matched", "#1f77b4", "D"), ("edge", "#d62728", "s"), ("ratio", "#2ca02c", "^")):
            ax.scatter(top[f"{key}_z"], y, c=c, marker=mk, s=30, label=key)
        ax.set_xscale("symlog", linthresh=2)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels((top.ligand + "→" + top.receptor + "  " + top.sender + "→" + top.receiver).tolist(), fontsize=7)
        ax.legend(fontsize=8)
        ax.set_title(f"U5-a 6-NN, {a}: top 40 rows by |compositional z| (500 perms)")
        fig.tight_layout()
        save_figure(fig, PHASE, f"{TAG}_3d_z_comparison_6NN_{a}", f"z of C/M/E/Q, top 40 rows, {a}", log)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, (l, r, s, rc) in zip(axes, HEAD):
        h = khead[(khead.sender == s) & (khead.receiver == rc)]
        ax.plot(h.k, h.edge_z, "s-", c="#d62728", label="edge z")
        ax.plot(h.k, h.ratio_z, "^--", c="#2ca02c", label="ratio z")
        ax.axhline(1.645, color="#aaa", ls=":")
        ax.set_xlabel("k (kNN)")
        ax.set_title(f"{l}→{r} {s}→{rc}", fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("U5-a: headline rows vs k (marker labels, 500 perms)")
    save_figure(fig, PHASE, f"{TAG}_3e_headline_k_sensitivity", "Headline edge/ratio z across k", log)

    fig, ax = plt.subplots(figsize=(9, 5))
    sub = xm[xm.n_edge_evaluable == 16]
    ax.scatter(sub.n_edge_sig, sub.n_ratio_sig, s=8, alpha=0.4, c="#555")
    for (l, r, s, rc) in HEAD:
        h = sub[(sub.ligand == l) & (sub.receptor == r) & (sub.sender == s) & (sub.receiver == rc)]
        ax.scatter(h.n_edge_sig, h.n_ratio_sig, s=80, c="#d62728", marker="*")
        for row in h.itertuples():
            ax.annotate(f"{l}→{r} {s}→{rc}", (row.n_edge_sig, row.n_ratio_sig), fontsize=7)
    ax.set_xlabel("sections (of 16) with edge score significant")
    ax.set_ylabel("sections (of 16) with edge/matched ratio significant")
    ax.set_title("Cross-section: edge significance is common, ratio significance is rare (marker labels, 6-NN)")
    save_figure(fig, PHASE, "cross_section_edge_vs_ratio_sig_counts", "Per-interaction counts of sections with edge vs ratio significance", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

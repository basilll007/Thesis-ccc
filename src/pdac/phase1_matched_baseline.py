"""Phase 1: matched dyadic baseline on the breast Xenium data.

Question: when the compositional score (Lbar_s + Rbar_r)/2 and the edge-restricted score
mean_{edges s->r} L_i R_j disagree, is the disagreement caused by the score FORM (sum vs product)
or by PROXIMITY (restricting to graph edges)?

Decomposition used here:
  C (compositional) vs M (matched = Lbar_s * Rbar_r): same cells, no graph, different form  -> score-form effect
  M vs E (edge):                                          same form, graph restriction added   -> proximity effect
  Q = E / M: edge enrichment over the matched non-spatial expectation, tested under the same null.

All four scores use the SAME 500 label permutations as the paper's 6-NN nulls
(spatial_null.generate_label_permutations, numpy default_rng(42)), on the stored Squidpy 6-NN graph.

Families:
  paper  : the 22 interactions of results/baseline_lr_ranked.csv (CXCL12-CXCR4, CD274-PDCD1; per-cell labels)
  extra  : CD86-CTLA4 and PTN-SDC4, all 16 sender-receiver combinations (per-cell labels)
  leiden : CXCL12-CXCR4, CD274-PDCD1 under Leiden-consensus labels (annotation sensitivity),
           combinations with compositional score > 0 (same rule as src/run_leiden_robustness.py)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO_ROOT, SEED, finish_run, register, save_csv, save_figure, save_json, save_npz, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import ccc_audit as ca  # noqa: E402
from pipeline import MARKERS  # noqa: E402
from spatial_null import generate_label_permutations  # noqa: E402

PHASE = "phase1"
N_PERMS = 500
SPARSE_GENE_MIN_CELLS = 20
PAPER_PAIRS = [("CXCL12", "CXCR4"), ("CD274", "PDCD1")]
EXTRA_PAIRS = [("CD86", "CTLA4"), ("PTN", "SDC4")]


def classify(row) -> str:
    c, e, m = row["comp_sig"], row["edge_sig"], row["matched_sig"]
    if row["edge_status"] != "evaluated":
        return "edge not evaluable (0 edges)"
    if c == e:
        return "C and E agree"
    if m == e:
        return "C/E disagree: SCORE-FORM (matched baseline reproduces edge verdict)"
    return "C/E disagree: PROXIMITY (matched baseline reproduces compositional verdict)"


def run_family(name, adata, labels_str, types, rows, cols, lr_pairs, interactions, log):
    t0 = time.time()
    ti = {t: k for k, t in enumerate(types)}
    lab = np.array([ti[x] for x in labels_str], dtype=np.int64)
    perms_str = generate_label_permutations(labels_str, n_perms=N_PERMS, seed=SEED)
    perms = [np.array([ti[x] for x in p], dtype=np.int64) for p in perms_str]
    log.info("[%s] %d label permutations via spatial_null.generate_label_permutations (default_rng(%d)); "
             "%d types; %d LR pairs; %d interactions", name, N_PERMS, SEED, len(types), len(lr_pairs), len(interactions))
    S = ca.setup(adata, lr_pairs, rows, cols)
    obs, null, nze = ca.run_null(lab, len(types), S, perms)
    df, null_rows = ca.tidy(interactions, types, S, obs, null, nze)
    # detection counts (raw > 0) for sparse-gene diagnostic
    det = {g: int((S.X[:, k] > 0).sum()) for k, g in enumerate(S.genes)}
    df["ligand_n_cells_detected"] = df["ligand"].map(det)
    df["receptor_n_cells_detected"] = df["receptor"].map(det)
    df["sparse_gene_diagnostic"] = (df[["ligand_n_cells_detected", "receptor_n_cells_detected"]].min(axis=1)
                                    < SPARSE_GENE_MIN_CELLS)
    df["decomposition"] = df.apply(classify, axis=1)
    df.insert(0, "family", name)
    log.info("[%s] done in %.1fs", name, time.time() - t0)
    return df, null_rows


def main() -> None:
    log = start_run("phase1_matched_baseline")
    adata = ad.read_h5ad(REPO_ROOT / "data" / "processed" / "xenium_breast_baseline.h5ad")
    rows, cols = ca.graph_edges(adata.obsp["spatial_connectivities"])
    log.info("Breast AnnData: %d cells; 6-NN graph (stored adata.obsp['spatial_connectivities']): %d directed edges",
             adata.n_obs, len(rows))

    # ---------------- paper family (per-cell labels) ----------------
    labels = np.asarray(adata.obs["cell_type"].astype(str).to_numpy())
    types = sorted(set(labels))
    base = pd.read_csv(REPO_ROOT / "results" / "baseline_lr_ranked.csv")
    paper_inter = [(r.sender, r.receiver, r.ligand, r.receptor) for r in base.itertuples()]
    paper, paper_null = run_family("paper", adata, labels, types, rows, cols, PAPER_PAIRS, paper_inter, log)

    # reproduction check vs existing paper files
    comp_prev = pd.read_csv(REPO_ROOT / "results" / "spatial_null_corrected.csv")
    edge_prev = pd.read_csv(REPO_ROOT / "results" / "spatial_edge_null.csv")
    k = ["sender", "receiver", "ligand", "receptor"]
    chk = paper.merge(comp_prev[k + ["observed_score", "z_score", "empirical_pvalue"]], on=k) \
               .merge(edge_prev[k + ["edge_observed_score", "edge_z_score", "edge_empirical_pvalue"]], on=k)
    repro = {
        "n_rows_matched": int(len(chk)),
        "max_abs_diff_comp_obs": float(np.max(np.abs(chk["comp_obs"] - chk["observed_score"]))),
        "max_abs_diff_comp_z": float(np.max(np.abs(chk["comp_z"] - chk["z_score"]))),
        "max_abs_diff_comp_p_legacy": float(np.max(np.abs(chk["comp_p_legacy"] - chk["empirical_pvalue"]))),
        "max_abs_diff_edge_obs": float(np.nanmax(np.abs(chk["edge_obs"] - chk["edge_observed_score"]))),
        "max_abs_diff_edge_z": float(np.nanmax(np.abs(chk["edge_z"] - chk["edge_z_score"]))),
        "max_abs_diff_edge_p_legacy": float(np.nanmax(np.abs(chk["edge_p_legacy"] - chk["edge_empirical_pvalue"]))),
    }
    log.info("Reproduction vs results/spatial_null_corrected.csv + spatial_edge_null.csv: %s", repro)

    # ---------------- extra family ----------------
    extra_inter = [(s, r, l, rc) for (l, rc) in EXTRA_PAIRS for s in types for r in types]
    extra, extra_null = run_family("extra", adata, labels, types, rows, cols, EXTRA_PAIRS, extra_inter, log)

    # ---------------- leiden-consensus family ----------------
    scores = {}
    for lbl, markers in MARKERS.items():
        present = [g for g in markers if g in adata.var_names]
        scores[lbl] = np.asarray(adata[:, present].X.mean(axis=1)).ravel()
    sf = pd.DataFrame(scores, index=adata.obs_names)
    consensus = sf.groupby(adata.obs["leiden"], observed=False).mean().idxmax(axis=1)
    lab_leiden = np.asarray(adata.obs["leiden"].map(consensus).astype(str).to_numpy())
    types_l = sorted(set(lab_leiden))
    log.info("Leiden-consensus types: %s", pd.Series(lab_leiden).value_counts().to_dict())
    # candidate rows: compositional score > 0 (rule of run_leiden_robustness.py)
    S_tmp = ca.setup(adata, PAPER_PAIRS, rows[:1], cols[:1])
    til = {t: i for i, t in enumerate(types_l)}
    o_tmp = ca.scores_for_labels(np.array([til[x] for x in lab_leiden]), len(types_l), S_tmp, zero_fill_edge=False)
    leiden_inter = [(s, r, l, rc) for pi_, (l, rc) in enumerate(PAPER_PAIRS) for s in types_l for r in types_l
                    if o_tmp["comp"][pi_, til[s], til[r]] > 0]
    leiden, leiden_null = run_family("leiden", adata, lab_leiden, types_l, rows, cols, PAPER_PAIRS, leiden_inter, log)

    allres = pd.concat([paper, extra, leiden], ignore_index=True)
    save_csv(paper, PHASE, "breast_matched_baseline_paper22.csv",
             "Per-axis C/E/M/Q scores, null stats, BH q, decomposition for the 22 paper interactions (per-cell labels, 6-NN)", log)
    save_csv(extra, PHASE, "breast_matched_baseline_extra_CD86CTLA4_PTNSDC4.csv",
             "Same table for CD86-CTLA4 and PTN-SDC4, all 16 sender-receiver combinations (per-cell labels, 6-NN)", log)
    save_csv(leiden, PHASE, "breast_matched_baseline_leiden.csv",
             "Same table for paper pairs under Leiden-consensus labels (annotation sensitivity)", log)
    save_csv(allres, PHASE, "breast_matched_baseline_all_families.csv", "All three families concatenated", log)
    save_json(repro, PHASE, "reproduction_check_vs_paper_files.json",
              "Max abs differences of recomputed C/E scores, z and p vs results/spatial_null_corrected.csv and spatial_edge_null.csv", log)
    for nm, nl, df in (("paper", paper_null, paper), ("extra", extra_null, extra), ("leiden", leiden_null, leiden)):
        save_npz(PHASE, f"null_distributions_{nm}.npz",
                 f"Null distributions (500 perms x interactions) for comp/edge/matched/ratio, family={nm}; row order = CSV order",
                 log, **{f"{k}_null": v for k, v in nl.items()},
                 interaction_labels=np.array([f"{a}|{b}|{c}|{d}" for a, b, c, d in
                                              df[["sender", "receiver", "ligand", "receptor"]].itertuples(index=False)]))

    # ---------------- decomposition summary and gate ----------------
    summ_rows = []
    for fam, df in (("paper", paper), ("extra", extra), ("leiden", leiden)):
        ev = df[df["edge_status"] == "evaluated"]
        dis = ev[ev["comp_sig"] != ev["edge_sig"]]
        sf_ = dis["decomposition"].str.contains("SCORE-FORM")
        summ_rows.append({
            "family": fam, "n_interactions": len(df), "n_edge_evaluable": len(ev),
            "n_comp_sig": int(df["comp_sig"].sum()), "n_edge_sig": int(df["edge_sig"].sum()),
            "n_matched_sig": int(df["matched_sig"].sum()), "n_ratio_sig": int(df["ratio_sig"].sum()),
            "n_C_E_disagree": len(dis), "n_disagree_score_form": int(sf_.sum()),
            "n_disagree_proximity": int((~sf_).sum()),
            "n_disagree_sparse_gene": int(dis["sparse_gene_diagnostic"].sum()),
            "n_C_M_disagree": int((df["comp_sig"] != df["matched_sig"]).sum()),
            "n_M_E_disagree": int((ev["matched_sig"] != ev["edge_sig"]).sum()),
            "spearman_z_comp_vs_edge": float(spearmanr(ev["comp_z"], ev["edge_z"]).statistic) if len(ev) > 2 else np.nan,
            "spearman_z_comp_vs_matched": float(spearmanr(df["comp_z"], df["matched_z"]).statistic) if len(df) > 2 else np.nan,
            "spearman_z_matched_vs_edge": float(spearmanr(ev["matched_z"], ev["edge_z"]).statistic) if len(ev) > 2 else np.nan,
            "n_comp_sig_bh": int((df["comp_q_bh"] < 0.05).sum()), "n_edge_sig_bh": int((df["edge_q_bh"] < 0.05).sum()),
            "n_matched_sig_bh": int((df["matched_q_bh"] < 0.05).sum()), "n_ratio_sig_bh": int((df["ratio_q_bh"] < 0.05).sum()),
        })
    summ = pd.DataFrame(summ_rows)
    save_csv(summ, PHASE, "decomposition_summary.csv",
             "Per-family counts: significant under each score, C/E disagreements split into score-form vs proximity", log)
    for r in summ.itertuples():
        log.info("[%s] %s", r.family, r._asdict())

    p = summ.set_index("family").loc["paper"]
    frac = p["n_disagree_score_form"] / p["n_C_E_disagree"] if p["n_C_E_disagree"] else np.nan
    gate = {
        "family": "paper (22 interactions, per-cell labels, 6-NN, 500 perms)",
        "n_C_E_disagreements": int(p["n_C_E_disagree"]),
        "n_explained_by_score_form": int(p["n_disagree_score_form"]),
        "n_attributed_to_proximity": int(p["n_disagree_proximity"]),
        "fraction_score_form": frac,
        "matched_baseline_explains_disagreement": bool(frac > 0.5) if np.isfinite(frac) else None,
        "rule": "A C/E disagreement is 'score-form' if the matched baseline M (no graph, product form) gives the same "
                "verdict as E; 'proximity' if M gives the same verdict as C. Verdict = p_plus1 < 0.05 (500 perms). "
                "'Explains' = more than half of disagreements are score-form.",
    }
    save_json(gate, PHASE, "phase1_gate.json", "Phase 1 gate: does the matched baseline explain the C-vs-E disagreement?", log)
    log.info("GATE: %s", gate)

    # ---------------- figures ----------------
    make_z_figure(paper, "paper", "22 paper interactions (per-cell labels, 6-NN, 500 perms)", log)
    make_z_figure(extra, "extra", "CD86–CTLA4 and PTN–SDC4, all 16 combinations (per-cell, 6-NN)", log)
    make_z_figure(leiden, "leiden", "Paper pairs under Leiden-consensus labels (6-NN)", log)
    make_null_figure(paper, paper_null, log)
    make_ratio_figure(paper, extra, log)
    finish_run(log, PHASE)


def _lab(df):
    return (df["ligand"] + "→" + df["receptor"] + "  " + df["sender"] + "→" + df["receiver"]).tolist()


def make_z_figure(df, tag, title, log):
    d = df.sort_values("comp_z").reset_index(drop=True)
    y = np.arange(len(d))
    fig, axes = plt.subplots(1, 2, figsize=(15, max(5, 0.32 * len(d) + 1.5)), sharey=True,
                             gridspec_kw={"width_ratios": [3, 1.3]})
    ax = axes[0]
    style = {"comp": ("#7f7f7f", "o", "compositional C=(L̄s+R̄r)/2"),
             "matched": ("#1f77b4", "D", "matched M=L̄s·R̄r (no graph)"),
             "edge": ("#d62728", "s", "edge E (6-NN)"),
             "ratio": ("#2ca02c", "^", "ratio Q=E/M")}
    for key, (c, m, lab) in style.items():
        ax.scatter(d[f"{key}_z"], y, c=c, marker=m, s=38, label=lab, zorder=3, alpha=0.9)
    ax.axvline(0, color="k", lw=0.6)
    ax.axvline(1.645, color="#999", ls="--", lw=0.8)
    ax.set_xscale("symlog", linthresh=2)
    ax.set_yticks(y)
    ax.set_yticklabels(_lab(d), fontsize=7.5)
    ax.set_xlabel("z vs label-permutation null (symlog axis; dashed: z=1.645)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(title, fontsize=10)
    ax2 = axes[1]
    for j, key in enumerate(["comp", "matched", "edge", "ratio"]):
        col = ["#2ca25f" if s else "#dddddd" for s in d[f"{key}_sig"]]
        ax2.scatter(np.full(len(d), j), y, c=col, s=60, marker="s", edgecolors="k", linewidths=0.4)
    for i, r in d.iterrows():
        if r["sparse_gene_diagnostic"]:
            ax2.text(4.0, i, "sparse", fontsize=6.5, va="center", color="#b15928")
        if "SCORE-FORM" in r["decomposition"]:
            ax2.text(4.9, i, "form", fontsize=6.5, va="center", color="#1f77b4")
        elif "PROXIMITY" in r["decomposition"]:
            ax2.text(4.9, i, "prox", fontsize=6.5, va="center", color="#d62728")
    ax2.set_xticks(range(4))
    ax2.set_xticklabels(["C", "M", "E", "Q"])
    ax2.set_xlim(-0.6, 5.8)
    ax2.set_title("p+1 < 0.05 (green)\nC/E split: form / prox", fontsize=9)
    fig.tight_layout()
    save_figure(fig, PHASE, f"phase1_z_compositional_vs_edge_vs_matched_{tag}",
                f"z-scores of compositional, matched baseline, edge and ratio scores; family={tag}", log)


def make_null_figure(df, nulls, log):
    pick = df[(df["receiver"].isin(["immune", "epithelial"])) & (df["sender"] == "fibroblast")
              & (df["ligand"] == "CXCL12")].index.tolist()
    pick += df[(df["ligand"] == "CD274")].index.tolist()[:2]
    fig, axes = plt.subplots(len(pick), 4, figsize=(15, 2.6 * len(pick)))
    for r_i, idx in enumerate(pick):
        for c_i, key in enumerate(["comp", "matched", "edge", "ratio"]):
            ax = axes[r_i, c_i]
            v = nulls[key][:, idx]
            ax.hist(v[np.isfinite(v)], bins=30, color="#9ecae1")
            o = df.loc[idx, f"{key}_obs"]
            if np.isfinite(o):
                ax.axvline(o, color="#de2d26", lw=2)
            ax.set_title(f"{_lab(df.loc[[idx]])[0]}\n{key}: z={df.loc[idx, f'{key}_z']:.2f}, "
                         f"p+1={df.loc[idx, f'{key}_p_plus1']:.3f}", fontsize=8)
    fig.suptitle("Phase 1: null distributions (500 label permutations, default_rng(42)); red line = observed", fontsize=11)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase1_null_distributions_selected_interactions",
                "Null histograms of C/M/E/Q for selected paper interactions", log)


def make_ratio_figure(paper, extra, log):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, df, t in ((axes[0], paper, "paper 22"), (axes[1], extra, "CD86–CTLA4 / PTN–SDC4")):
        ev = df[df["edge_status"] == "evaluated"]
        ax.scatter(ev["matched_z"], ev["edge_z"], c=np.where(ev["ratio_sig"], "#2ca02c", "#999999"), s=30)
        lim = [min(ev["matched_z"].min(), ev["edge_z"].min()) - 1, max(ev["matched_z"].max(), ev["edge_z"].max()) + 1]
        ax.plot(lim, lim, "k--", lw=0.7)
        ax.axhline(1.645, color="#ccc", lw=0.7)
        ax.axvline(1.645, color="#ccc", lw=0.7)
        ax.set_xlabel("matched baseline z (no graph)")
        ax.set_ylabel("edge z (6-NN)")
        ax.set_title(f"{t}: green = ratio Q=E/M significant (p+1<0.05)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase1_matched_vs_edge_z_scatter",
                "Matched-baseline z vs edge z; colour marks significant edge/matched ratio", log)


if __name__ == "__main__":
    main()

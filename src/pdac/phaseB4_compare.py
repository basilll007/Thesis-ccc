"""Phase B4 (v2): compare the CCC audit across annotations, 16 PDAC-SMI sections.

Annotations: v1 marker (results/pdac/phase3/cross_section_6NN_all_rows.csv), v1 author_mapped
(results/pdac/phase3/author_annotation_6NN_all_sections.csv), v2 protein_gmm / protein_mad (pre-registered; failed
validation) and consensus_T (post-hoc) (results/pdac/phaseB/v2_ccc_6NN_all_rows.csv.gz).
All tables use add-one p; BH within section x annotation (v1 author file: within section = same thing).

Outputs:
  headline (CXCL12->CXCR4, CAF_fibroblast->T_cell): per section x annotation z for C/M/E/Q, verdicts, BH; counts
  pairwise annotation comparison on shared rows: Spearman of z per score, verdict changes (U5-a and pooled)
  'edge mirrors product?' diagnostics per annotation: Spearman(E z, M z); share of E-significant rows that are
  M-significant / Q-significant; with and without background-flagged axes
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import ccc_audit as ca  # noqa: E402

PHASE = "phaseB"
ALPHA = 0.05
KEY = ["section", "sender", "receiver", "ligand", "receptor"]
ANN_ORDER = ["marker", "author_mapped", "protein_gmm", "protein_mad", "consensus_T"]
ANN_NOTE = {"marker": "v1 marker scoring (T precision 0.35 pooled)", "author_mapped": "v1 author labels (reference)",
            "protein_gmm": "v2 pre-registered protein (FAILED validation)", "protein_mad": "v2 protein MAD (FAILED validation)",
            "consensus_T": "v2 POST-HOC author AND marker T"}


def harmonise(df, ann, qsuffix):
    out = df[KEY + ["treatment_status", "n_qualifying_edges", "edge_status"]].copy()
    out["annotation"] = ann
    for k in ca.SCORES:
        out[f"{k}_obs"] = df[f"{k}_obs"]
        out[f"{k}_z"] = df[f"{k}_z"]
        out[f"{k}_p"] = df[f"{k}_p_plus1"]
        out[f"{k}_q"] = df[f"{k}_q_bh_{qsuffix}"]
        out[f"{k}_sig"] = df[f"{k}_p_plus1"] < ALPHA
        out[f"{k}_sig_bh"] = df[f"{k}_q_bh_{qsuffix}"] < ALPHA
    return out


def load_all(log):
    xs = pd.read_csv(RESULTS_PDAC / "phase3" / "cross_section_6NN_all_rows.csv")
    mk = harmonise(xs[xs.annotation == "cell_type_marker"], "marker", "within_section_annotation")
    au = harmonise(pd.read_csv(RESULTS_PDAC / "phase3" / "author_annotation_6NN_all_sections.csv"), "author_mapped",
                   "within_section")
    v2 = pd.read_csv(RESULTS_PDAC / "phaseB" / "v2_ccc_6NN_all_rows.csv.gz")
    v2h = [harmonise(v2[v2.annotation == a], a, "within_section_annotation") for a in ("protein_gmm", "protein_mad", "consensus_T")]
    allt = pd.concat([mk, au] + v2h, ignore_index=True)
    log.info("rows per annotation: %s", allt.annotation.value_counts().to_dict())
    return allt


def headline(allt, log):
    hl = allt[(allt.ligand == "CXCL12") & (allt.receptor == "CXCR4") & (allt.sender == "CAF_fibroblast") & (allt.receiver == "T_cell")]
    sections = sorted(allt.section.unique())
    full = pd.MultiIndex.from_product([ANN_ORDER, sections], names=["annotation", "section"]).to_frame(index=False)
    hl = full.merge(hl, on=["annotation", "section"], how="left")
    hl["evaluable"] = hl.edge_status.eq("evaluated")
    save_csv(hl, PHASE, "headline_CAF_T_by_annotation_section.csv",
             "CXCL12-CXCR4 CAF->T_cell per section x annotation (NaN = no T_cell class or 0 edges)", log)
    summ = hl.groupby("annotation").agg(
        n_sections_evaluable=("evaluable", "sum"), median_n_edges=("n_qualifying_edges", "median"),
        n_comp_sig=("comp_sig", lambda x: int(x.fillna(False).sum())),
        n_matched_sig=("matched_sig", lambda x: int(x.fillna(False).sum())),
        n_edge_sig=("edge_sig", lambda x: int(x.fillna(False).sum())),
        n_edge_sig_bh=("edge_sig_bh", lambda x: int(x.fillna(False).sum())),
        n_ratio_sig=("ratio_sig", lambda x: int(x.fillna(False).sum())),
        n_ratio_sig_bh=("ratio_sig_bh", lambda x: int(x.fillna(False).sum())),
        n_ratio_z_pos=("ratio_z", lambda x: int((x > 0).sum())),
        median_edge_z=("edge_z", "median"), median_matched_z=("matched_z", "median"),
        median_ratio_z=("ratio_z", "median"), median_ratio_obs=("ratio_obs", "median")).reindex(ANN_ORDER).reset_index()
    summ["note"] = summ.annotation.map(ANN_NOTE)
    save_csv(summ, PHASE, "headline_CAF_T_summary_by_annotation.csv",
             "Headline: sections significant per score, uncorrected and BH within section, per annotation", log)
    for r in summ.itertuples():
        log.info("HEADLINE %s: evaluable %d/16; E sig %d (BH %d); Q sig %d (BH %d); Q z>0 in %d; median zE %.2f zQ %.2f Q %.3f",
                 r.annotation, r.n_sections_evaluable, r.n_edge_sig, r.n_edge_sig_bh, r.n_ratio_sig, r.n_ratio_sig_bh,
                 r.n_ratio_z_pos, r.median_edge_z, r.median_ratio_z, r.median_ratio_obs)
    return hl, summ


def pairwise(allt, log):
    rows = []
    for scope in ("U5-a", "all_16_sections"):
        sub = allt if scope == "all_16_sections" else allt[allt.section == scope]
        for a, b in combinations(ANN_ORDER, 2):
            A = sub[sub.annotation == a].set_index(KEY)
            B = sub[sub.annotation == b].set_index(KEY)
            j = A.join(B, how="inner", lsuffix="_a", rsuffix="_b")
            for k in ca.SCORES:
                x, y = j[f"{k}_z_a"], j[f"{k}_z_b"]
                ok = np.isfinite(x) & np.isfinite(y)
                rows.append({"scope": scope, "annotation_a": a, "annotation_b": b, "score": k, "n_shared_rows": int(ok.sum()),
                             "spearman_z": float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() > 2 else np.nan,
                             "n_verdict_changes": int((j.loc[ok, f"{k}_sig_a"] != j.loc[ok, f"{k}_sig_b"]).sum()),
                             "n_sig_a": int(j.loc[ok, f"{k}_sig_a"].sum()), "n_sig_b": int(j.loc[ok, f"{k}_sig_b"].sum())})
    pw = pd.DataFrame(rows)
    save_csv(pw, PHASE, "annotation_pairwise_comparison.csv",
             "Shared rows between annotations: Spearman of z and verdict changes per score (U5-a, pooled 16 sections)", log)
    return pw


def mirroring(allt, log):
    bgf = pd.read_csv(RESULTS_PDAC / "phase3" / "background_flag_per_axis.csv")
    flagged = set(map(tuple, bgf.loc[bgf.background_level_gene_diagnostic, ["ligand", "receptor"]].to_numpy()))
    rows = []
    for ann in ANN_ORDER:
        d0 = allt[(allt.annotation == ann) & (allt.edge_status == "evaluated")]
        for subset in ("all_axes", "unflagged_axes_only"):
            d = d0 if subset == "all_axes" else d0[[(l, r) not in flagged for l, r in zip(d0.ligand, d0.receptor)]]
            ok = np.isfinite(d.edge_z) & np.isfinite(d.matched_z)
            es = d.edge_sig
            rows.append({"annotation": ann, "axes": subset, "n_rows": len(d),
                         "spearman_edge_z_vs_matched_z": float(spearmanr(d.edge_z[ok], d.matched_z[ok]).statistic),
                         "spearman_edge_z_vs_ratio_z": float(spearmanr(d.edge_z[np.isfinite(d.ratio_z)], d.ratio_z[np.isfinite(d.ratio_z)]).statistic),
                         "n_edge_sig": int(es.sum()), "n_edge_sig_and_matched_sig": int((es & d.matched_sig).sum()),
                         "n_edge_sig_and_ratio_sig": int((es & d.ratio_sig).sum()),
                         "n_edge_sig_and_ratio_sig_bh": int((es & d.ratio_sig_bh).sum()),
                         "frac_edge_sig_also_matched_sig": float((es & d.matched_sig).sum() / max(es.sum(), 1)),
                         "frac_edge_sig_also_ratio_sig": float((es & d.ratio_sig).sum() / max(es.sum(), 1)),
                         "frac_edge_sig_also_ratio_sig_bh": float((es & d.ratio_sig_bh).sum() / max(es.sum(), 1))})
    mi = pd.DataFrame(rows)
    save_csv(mi, PHASE, "edge_vs_product_mirroring_by_annotation.csv",
             "Does the edge score mirror the matched cell-type product? Spearman(E z, M z) and overlap of verdicts", log)
    for r in mi.itertuples():
        log.info("MIRROR %s/%s: rho(Ez,Mz)=%.3f; E-sig %d, of which M-sig %.3f, Q-sig %.3f, Q-sig BH %.3f", r.annotation, r.axes,
                 r.spearman_edge_z_vs_matched_z, r.n_edge_sig, r.frac_edge_sig_also_matched_sig, r.frac_edge_sig_also_ratio_sig,
                 r.frac_edge_sig_also_ratio_sig_bh)
    return mi


def figure(hl, pw, log):
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5), gridspec_kw={"width_ratios": [2.2, 2.2, 1.4]})
    sections = sorted(hl.section.unique())
    xi = {s: i for i, s in enumerate(sections)}
    cols = {"marker": "#7f7f7f", "author_mapped": "#1f77b4", "protein_gmm": "#ff7f0e", "protein_mad": "#ffbb78",
            "consensus_T": "#d62728"}
    for ax, score, ttl in ((axes[0], "edge_z", "edge E z"), (axes[1], "ratio_z", "ratio Q = E/M z")):
        for k, a in enumerate(ANN_ORDER):
            d = hl[hl.annotation == a]
            off = (k - 2) * 0.14
            sig = d[score.replace("_z", "_sig")].fillna(False).astype(bool)
            ax.scatter([xi[s] + off for s in d.section], d[score], c=cols[a], s=np.where(sig, 40, 14),
                       edgecolors=np.where(d[score.replace("_z", "_sig_bh")].fillna(False).astype(bool), "k", "none"),
                       label=a)
        ax.axhline(1.645, color="#999", ls=":", lw=0.8)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(sections)))
        ax.set_xticklabels(sections, rotation=90, fontsize=7)
        ax.set_title(f"CXCL12→CXCR4 CAF→T_cell: {ttl} (large = p+1<.05; black edge = BH<.05)", fontsize=9)
        if score == "edge_z":
            ax.set_yscale("symlog", linthresh=2)
        ax.legend(fontsize=7)
    ax = axes[2]
    p = pw[(pw.scope == "all_16_sections") & (pw.score.isin(["edge", "ratio"]))]
    lab = [f"{a} vs {b}" for a, b in zip(p.annotation_a, p.annotation_b)]
    piv = p.assign(pair=lab).pivot(index="pair", columns="score", values="spearman_z")
    im = ax.imshow(piv.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels(piv.index, fontsize=7)
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.iat[i, j]:.2f}", ha="center", va="center", fontsize=7, color="w")
    fig.colorbar(im, ax=ax, fraction=0.05)
    ax.set_title("Spearman of z on shared rows, 16 sections", fontsize=9)
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseB4_headline_and_annotation_agreement",
                "Headline edge/ratio z per section under 5 annotations; Spearman of z between annotations", log)


def main() -> None:
    log = start_run("phaseB4_compare")
    allt = load_all(log)
    hl, summ = headline(allt, log)
    pw = pairwise(allt, log)
    mi = mirroring(allt, log)
    u5 = hl[hl.section == "U5-a"][["annotation", "n_qualifying_edges", "comp_z", "matched_z", "edge_z", "ratio_z", "ratio_obs",
                                   "ratio_p", "ratio_q", "edge_sig", "ratio_sig", "ratio_sig_bh"]]
    save_csv(u5, PHASE, "headline_CAF_T_U5-a_by_annotation.csv", "U5-a headline row under each annotation", log)
    save_json({"headline_summary": summ.to_dict(orient="records"), "U5-a": u5.to_dict(orient="records"),
               "mirroring": mi.to_dict(orient="records")}, PHASE, "phaseB4_summary.json", "Phase B4 key numbers", log)
    figure(hl, pw, log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

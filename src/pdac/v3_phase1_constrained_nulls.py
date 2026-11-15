"""v3 Phase 1: contact vs zonation. C/E/M/Q under the global null and spatially constrained nulls
(block 200/100/50 um, distance-stratified), 500 perms, seed 42, 6-NN, 16 sections, 3 annotations.
Plan: docs/v3_plans/phase1_contact_vs_zonation.md (committed before this script ran).
Usage: python src/pdac/v3_phase1_constrained_nulls.py [--sections A B] [--aggregate-only]
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, finish_run, mem_status, out_dir, register, save_csv, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import combine_pvalues, spearmanr, wilcoxon  # noqa: E402

import ccc_audit as ca  # noqa: E402
import v3_nulls as vn  # noqa: E402
from phase3_pdac_audit import add_flags  # noqa: E402

PHASE = "v3_phase1"
N_PERMS = 500
HEAD = ("CAF_fibroblast", "T_cell", "CXCL12", "CXCR4")
SCORE_LAB = {"edge": "E", "matched": "M", "ratio": "Q"}
P_NOTE = "add-one (k+1)/(n+1)"


def run_section(section, log, od):
    todo = [(a, n) for a in vn.ANNOTATIONS for n in vn.NULL_ORDER if not (od / f"{section}__{a}__{n}.csv.gz").exists()]
    if not todo:
        log.info("section %s: all files present (resume)", section)
        return
    t0 = time.time()
    adata, meta, labs, pairs, det = vn.load_section(section, log)
    xy = np.asarray(adata.obsm["spatial"], float)
    rows, cols = vn.knn_edges(xy, 6)
    S = ca.setup(adata, pairs, rows, cols)
    null_info = {}
    for a in vn.ANNOTATIONS:
        labels_str = labs[a].astype(str).to_numpy()
        types, lab = vn.encode(labels_str)
        for n in vn.NULL_ORDER:
            if (a, n) not in todo:
                continue
            perms, info = vn.make_perms(n, labels_str, types, lab, xy, N_PERMS)
            null_info[f"{a}/{n}"] = info
            p = od / f"{section}__{a}__{n}.csv.gz"
            if perms is None:
                log.warning("section %s / %s / %s SKIPPED: %s", section, a, n, info)
                pd.DataFrame({"section": [section], "annotation": [a], "null": [n], "skipped": [str(info)]}).to_csv(p, index=False)
                register(p, PHASE, f"SKIPPED marker: {section} {a} {n} ({info.get('skipped')})")
                continue
            t1 = time.time()
            obs, null, nze = ca.run_null(lab, len(types), S, perms)
            inter = [(s, r, l, rc) for (l, rc) in pairs for s in types for r in types]
            df, _ = ca.tidy(inter, types, S, obs, null, nze)
            df = df.drop(columns=[c for c in df.columns if c.endswith("_sig") or c.endswith("_q_bh")])
            df = add_flags(df, det)
            df.insert(0, "null", n)
            df.insert(0, "annotation", a)
            df.insert(0, "section", section)
            df.insert(1, "treatment_status", meta["treatment_status"])
            df.to_csv(p, index=False)
            register(p, PHASE, f"6-NN C/E/M/Q ({N_PERMS} perms) section {section}, annotation {a}, null {n}")
            del null, obs, perms
            gc.collect()
            log.info("section %s / %s / %s: %d rows in %.1fs; %s; %s", section, a, n, len(df), time.time() - t1, info, mem_status())
    save_json({"section": section, "n_cells": int(adata.n_obs), "n_edges": int(len(rows)), "null_info": null_info},
              f"{PHASE}/sections", f"{section}_null_info.json", f"Null construction details (tiles, bands) for section {section}", log)
    del adata, S
    gc.collect()
    log.info("section %s done in %.1fs; %s", section, time.time() - t0, mem_status())


def aggregate(log, od):
    tmap = vn.section_to_tumour()
    nonbg = set(vn.nonbg_axes())
    log.info("non-background axes (%d): %s", len(nonbg), sorted(nonbg))
    files = sorted(od.glob("*.csv.gz"))
    frames, skipped = [], []
    for f in files:
        d = pd.read_csv(f)
        if "skipped" in d.columns:
            skipped.append(d)
            continue
        d = d[[c for c in d.columns if not c.startswith("comp_")]]
        d = d[[tuple(x) in nonbg for x in d[["ligand", "receptor"]].to_numpy()]]
        frames.append(d)
    allt = pd.concat(frames, ignore_index=True)
    skipped = pd.concat(skipped, ignore_index=True) if skipped else pd.DataFrame(columns=["section", "annotation", "null", "skipped"])
    save_csv(skipped, PHASE, "skipped_runs.csv",
             "Section x annotation x null runs skipped (distance null needs >=50 tumour and non-tumour cells)", log)
    allt["tumour"] = allt.section.map(tmap)
    for k in ("edge", "matched", "ratio"):
        allt[f"{k}_q_bh"] = allt.groupby(["section", "annotation", "null"])[f"{k}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    p = out_dir(PHASE) / "nonbg_catalogue_by_null.csv.gz"
    allt.to_csv(p, index=False)
    register(p, PHASE, "7 non-background axes x all sender/receiver x 16 sections x 3 annotations x 5 nulls: "
                       "E/M/Q obs, z, add-one p, BH within section/annotation/null")
    log.info("non-bg catalogue rows: %d (%d files, %d skipped)", len(allt), len(files), len(skipped))

    # headline ----------------------------------------------------------------------------------------
    hl = allt[(allt.sender == HEAD[0]) & (allt.receiver == HEAD[1]) & (allt.ligand == HEAD[2]) & (allt.receptor == HEAD[3])].copy()
    hl["null"] = pd.Categorical(hl["null"], vn.NULL_ORDER, ordered=True)
    hl = hl.sort_values(["annotation", "null", "section"])
    hl["p_value_convention"] = P_NOTE
    save_csv(hl, PHASE, "headline_by_null.csv", "CXCL12->CXCR4 CAF->T_cell per section x annotation x null (p = add-one)", log)
    hs = []
    for (a, n), g in hl.groupby(["annotation", "null"], observed=True):
        ev = g[g.edge_status == "evaluated"]
        rec = {"annotation": a, "null": n, "n_sections_evaluable": len(ev)}
        for k, lab in SCORE_LAB.items():
            rec[f"n_{lab}_p05"] = int((ev[f"{k}_p_plus1"] < 0.05).sum())
            rec[f"n_{lab}_bh05"] = int((ev[f"{k}_q_bh"] < 0.05).sum())
            rec[f"median_z_{lab}"] = float(ev[f"{k}_z"].median())
        rec.update({"n_zQ_pos": int((ev.ratio_z > 0).sum()), "median_Q_obs": float(ev.ratio_obs.median()),
                    "median_null_mean_Q": float(ev.ratio_null_mean.median()), "median_null_sd_Q": float(ev.ratio_null_sd.median()),
                    "p_value_convention": P_NOTE})
        hs.append(rec)
    hs = pd.DataFrame(hs)
    save_csv(hs, PHASE, "headline_summary_by_null.csv",
             "Headline per annotation x null: sections with p<0.05 / BH q<0.05 per score, median z, median Q", log)
    for r in hs.itertuples():
        log.info("HEADLINE %s %s: eval %d; Q p<.05 %d (BH %d); zQ>0 %d; median zE %.2f zM %.2f zQ %.2f; Q %.3f (null mean %.3f sd %.3f)",
                 r.annotation, r.null, r.n_sections_evaluable, r.n_Q_p05, r.n_Q_bh05, r.n_zQ_pos, r.median_z_E, r.median_z_M,
                 r.median_z_Q, r.median_Q_obs, r.median_null_mean_Q, r.median_null_sd_Q)

    # tumour level (Phase D2 method) -----------------------------------------------------------------
    rows = []
    for (a, n), d in hl[hl.edge_status == "evaluated"].groupby(["annotation", "null"], observed=True):
        for k, lab in SCORE_LAB.items():
            for t, g in d.groupby("tumour"):
                pv = g[f"{k}_p_plus1"].to_numpy()
                pt = float(combine_pvalues(pv, method="fisher").pvalue) if len(pv) > 1 else float(pv[0])
                rows.append({"annotation": a, "null": n, "score": lab, "tumour": t, "n_sections": len(g),
                             "sections": ";".join(g.section), "tumour_p_fisher": pt, "tumour_mean_z": float(g[f"{k}_z"].mean()),
                             "tumour_mean_obs": float(g[f"{k}_obs"].mean()), "n_edges_total": int(g.n_qualifying_edges.sum())})
    tt = pd.DataFrame(rows)
    tt["tumour_q_bh"] = tt.groupby(["annotation", "null", "score"]).tumour_p_fisher.transform(
        lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    tt["sig_bh"] = tt.tumour_q_bh < 0.05
    tt["p_value_convention"] = "section p = add-one; tumour p = Fisher of section p; BH across tumours"
    save_csv(tt, PHASE, "tumour_level_headline_by_null.csv",
             "Headline per tumour x annotation x null x score: Fisher p, BH q, mean z, mean obs", log)
    ts = []
    for (a, n, s), g in tt.groupby(["annotation", "null", "score"]):
        z = g.tumour_mean_z.to_numpy()
        w = wilcoxon(z, alternative="greater") if len(z) >= 5 and np.any(z != 0) else None
        ts.append({"annotation": a, "null": n, "score": s, "n_tumours": len(g),
                   "n_sig_uncorrected": int((g.tumour_p_fisher < 0.05).sum()), "n_sig_bh": int(g.sig_bh.sum()),
                   "min_q_bh": float(g.tumour_q_bh.min()), "n_mean_z_pos": int((z > 0).sum()),
                   "median_tumour_mean_z": float(np.median(z)), "mean_tumour_mean_z": float(np.mean(z)),
                   "wilcoxon_z_gt0_p": float(w.pvalue) if w is not None else np.nan,
                   "n_mean_obs_gt1": int((g.tumour_mean_obs > 1).sum()) if s == "Q" else np.nan})
    ts = pd.DataFrame(ts)
    ts["null"] = pd.Categorical(ts["null"], vn.NULL_ORDER, ordered=True)
    ts = ts.sort_values(["annotation", "score", "null"]).reset_index(drop=True)
    ts["p_value_convention"] = "section p = add-one; tumour p = Fisher; BH across tumours; Wilcoxon exploratory"
    save_csv(ts, PHASE, "tumour_level_summary_by_null.csv",
             "Tumour-level headline summary per annotation x null x score (Fisher+BH, z>0 count, Wilcoxon exploratory, Q>1 count)", log)
    for r in ts[ts.score == "Q"].itertuples():
        log.info("TUMOUR Q %s %s: n %d; sig %d (BH %d, min q %.3f); z>0 %d; median z %.3f; Wilcoxon p %.4f; Q>1 %s",
                 r.annotation, r.null, r.n_tumours, r.n_sig_uncorrected, r.n_sig_bh, r.min_q_bh, r.n_mean_z_pos,
                 r.median_tumour_mean_z, r.wilcoxon_z_gt0_p, r.n_mean_obs_gt1)
    verdict = {}
    for a in vn.ANNOTATIONS:
        q = ts[(ts.annotation == a) & (ts.score == "Q")].copy()
        q["null"] = q["null"].astype(str)
        q = q.set_index("null")
        if "global" in q.index and "block50" in q.index:
            g0, b50 = q.loc["global"], q.loc["block50"]
            shrinks = abs(b50.median_tumour_mean_z) < abs(g0.median_tumour_mean_z) and b50.n_mean_z_pos <= g0.n_mean_z_pos
            verdict[a] = {"global_median_z": float(g0.median_tumour_mean_z), "block200_median_z": float(q.loc["block200"].median_tumour_mean_z),
                          "block100_median_z": float(q.loc["block100"].median_tumour_mean_z),
                          "block50_median_z": float(b50.median_tumour_mean_z),
                          "dist_median_z": float(q.loc["dist"].median_tumour_mean_z) if "dist" in q.index else None,
                          "global_n_zpos": int(g0.n_mean_z_pos), "block50_n_zpos": int(b50.n_mean_z_pos),
                          "shrinks_toward_zero": bool(shrinks)}
    log.info("SHRINK VERDICT: %s", verdict)

    # non-background catalogue summary ---------------------------------------------------------------
    ev = allt[allt.edge_status == "evaluated"]
    key = ["section", "sender", "receiver", "ligand", "receptor"]
    cs = []
    for (a, n), g in ev.groupby(["annotation", "null"]):
        ref = ev[(ev.annotation == a) & (ev["null"] == "global")].set_index(key)
        j = g.set_index(key).join(ref, how="inner", rsuffix="_glob")
        rec = {"annotation": a, "null": n, "n_rows_evaluated": len(g), "n_sections": g.section.nunique()}
        for k, lab in SCORE_LAB.items():
            rec[f"n_{lab}_p05"] = int((g[f"{k}_p_plus1"] < 0.05).sum())
            rec[f"frac_{lab}_p05"] = float((g[f"{k}_p_plus1"] < 0.05).mean())
            rec[f"n_{lab}_bh05"] = int((g[f"{k}_q_bh"] < 0.05).sum())
            x, y = j[f"{k}_z"], j[f"{k}_z_glob"]
            ok = np.isfinite(x) & np.isfinite(y)
            rec[f"spearman_z{lab}_vs_global"] = float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() > 2 else np.nan
            rec[f"median_z_{lab}"] = float(g[f"{k}_z"].median())
            rec[f"n_z{lab}_pos"] = int((g[f"{k}_z"] > 0).sum())
        rec["p_value_convention"] = P_NOTE + "; BH within section x annotation x null"
        cs.append(rec)
    cs = pd.DataFrame(cs)
    cs["null"] = pd.Categorical(cs["null"], vn.NULL_ORDER, ordered=True)
    cs = cs.sort_values(["annotation", "null"]).reset_index(drop=True)
    save_csv(cs, PHASE, "nonbg_summary_by_null.csv",
             "Non-background catalogue per annotation x null: rows with p<0.05 / BH q<0.05 per score, Spearman of z vs global null", log)
    for r in cs.itertuples():
        log.info("CATALOGUE %s %s: rows %d; Q p<.05 %d (%.3f) BH %d; E p<.05 %d; rho(zQ,global) %.3f rho(zE,global) %.3f",
                 r.annotation, r.null, r.n_rows_evaluated, r.n_Q_p05, r.frac_Q_p05, r.n_Q_bh05, r.n_E_p05,
                 r.spearman_zQ_vs_global, r.spearman_zE_vs_global)
    save_json({"shrink_verdict": verdict, "nulls": vn.NULL_ORDER, "n_perms": N_PERMS, "annotations": vn.ANNOTATIONS,
               "nonbg_axes": sorted(nonbg), "section_to_tumour": tmap, "n_skipped_runs": int(len(skipped))},
              PHASE, "phase1_summary.json", "Phase 1 summary: shrink verdict per annotation, settings", log)

    # figure --------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    xs = {n: i for i, n in enumerate(vn.NULL_ORDER)}
    cols = {"consensus_T": "#d62728", "marker": "#7f7f7f", "author_mapped": "#1f77b4"}
    ax = axes[0]
    h = hl[(hl.annotation == "consensus_T") & (hl.edge_status == "evaluated")].copy()
    h["null"] = h["null"].astype(str)
    for s, g in h.groupby("section"):
        g = g.iloc[np.argsort([xs[n] for n in g["null"]])]
        ax.plot([xs[n] for n in g["null"]], g.ratio_z, color="#9e9e9e", lw=0.7, alpha=0.7)
    tq = tt[(tt.annotation == "consensus_T") & (tt.score == "Q")]
    jit = np.random.default_rng(0)
    for n, g in tq.groupby("null"):
        ax.scatter(np.full(len(g), xs[n]) + jit.uniform(-0.12, 0.12, len(g)), g.tumour_mean_z, c="#d62728", s=18, zorder=3)
        ax.scatter([xs[n]], [g.tumour_mean_z.median()], marker="_", s=600, c="k", zorder=4)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(vn.NULL_ORDER)))
    ax.set_xticklabels(vn.NULL_ORDER)
    ax.set_ylabel("z of Q = E/M")
    ax.set_title("CXCL12→CXCR4 CAF→T, consensus_T\ngrey: sections; red: tumour means; black: median", fontsize=9)
    ax = axes[1]
    for a in vn.ANNOTATIONS:
        q = ts[(ts.annotation == a) & (ts.score == "Q")].sort_values("null")
        xq = [xs[str(n)] for n in q["null"]]
        ax.plot(xq, q.median_tumour_mean_z, "o-", color=cols[a], label=f"{a}: median tumour z_Q")
        ax.plot(xq, q.n_mean_z_pos / q.n_tumours, "s--", color=cols[a], alpha=0.5, label=f"{a}: frac tumours z_Q>0")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(vn.NULL_ORDER)))
    ax.set_xticklabels(vn.NULL_ORDER)
    ax.set_title("Headline at tumour level vs null constraint", fontsize=9)
    ax.legend(fontsize=6)
    ax = axes[2]
    for a in vn.ANNOTATIONS:
        c = cs[cs.annotation == a].sort_values("null")
        xc = [xs[str(n)] for n in c["null"]]
        ax.plot(xc, c.frac_Q_p05, "o-", color=cols[a], label=f"{a}: Q")
        ax.plot(xc, c.frac_E_p05, "^:", color=cols[a], alpha=0.6, label=f"{a}: E")
    ax.axhline(0.05, color="k", lw=0.5, ls="--")
    ax.set_xticks(range(len(vn.NULL_ORDER)))
    ax.set_xticklabels(vn.NULL_ORDER)
    ax.set_ylabel("fraction of rows with add-one p < 0.05")
    ax.set_title("7 non-background axes, all sender→receiver rows", fontsize=9)
    ax.legend(fontsize=6, ncol=2)
    fig.suptitle("v3 Phase 1: proximity signal under spatially constrained nulls (audit/diagnostic only; 500 perms, seed 42)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase1_Qz_vs_null",
                "Q z vs null constraint: headline sections/tumours, tumour-level summary, non-background catalogue", log)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--aggregate-only", action="store_true")
    ap.add_argument("--no-aggregate", action="store_true")
    args = ap.parse_args()
    log = start_run("v3_phase1_constrained_nulls")
    log.info("nulls=%s annotations=%s n_perms=%d; %s", vn.NULL_ORDER, vn.ANNOTATIONS, N_PERMS, mem_status())
    od = out_dir(f"{PHASE}/sections")
    if not args.aggregate_only:
        full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
        sections = sorted(full.obs["sample_id"].astype(str).unique())
        del full
        if args.sections:
            sections = [s for s in sections if s in args.sections]
        for s in sections:
            run_section(s, log, od)
    if not args.no_aggregate:
        aggregate(log, od)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

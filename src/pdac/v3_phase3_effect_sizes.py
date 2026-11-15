"""v3 Phase 3: effect sizes for the headline CXCL12->CXCR4 CAF->T_cell axis (plan: docs/v3_plans/phase3_effect_sizes_mixed_model.md).
3a  Q with spatial block-bootstrap 95% CI per section (100 um tiles, B=1000, seed 42), edge counts per estimate.
3b  Tumour-level mixed model log Q ~ 1 + (1|tumour), statsmodels MixedLM (REML); exploratory.
3c  z vs edge-count inflation: raw z next to z/sqrt(n_edges), from the Phase 1 global-null files.
Usage: python src/pdac/v3_phase3_effect_sizes.py [--sections ...] [--aggregate-only]
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, RESULTS_PDAC, SEED, finish_run, mem_status, out_dir, register, save_csv, save_figure, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import stats as sps  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import ccc_audit as ca  # noqa: E402
import v3_nulls as vn  # noqa: E402

PHASE = "v3_phase3"
TILE_UM = 100.0
B = 1000
SENDER, RECEIVER = "CAF_fibroblast", "T_cell"
HEAD = ("CXCL12", "CXCR4")


def block_bootstrap(xy, lab, types, rows, cols, X, gi, axes_list, rng):
    """Return list of dicts (one per axis) with Q_obs and bootstrap CI; tile draws shared across axes."""
    tiles = vn.tile_codes(xy, TILE_UM)
    _, tiles = np.unique(tiles, return_inverse=True)
    n_tiles = tiles.max() + 1
    W = rng.multinomial(n_tiles, np.full(n_tiles, 1.0 / n_tiles), size=B).astype(float)  # B x n_tiles
    iS, iR = types.index(SENDER), types.index(RECEIVER)
    s_cells, r_cells = np.flatnonzero(lab == iS), np.flatnonzero(lab == iR)
    e = np.flatnonzero((lab[rows] == iS) & (lab[cols] == iR))
    den_s = W @ np.bincount(tiles[s_cells], minlength=n_tiles)
    den_r = W @ np.bincount(tiles[r_cells], minlength=n_tiles)
    den_e = W @ np.bincount(tiles[rows[e]], minlength=n_tiles)
    out = []
    for lig, rec in axes_list:
        L, R = X[:, gi[lig]], X[:, gi[rec]]
        prods = L[rows[e]] * R[cols[e]]
        Ls = (W @ np.bincount(tiles[s_cells], weights=L[s_cells], minlength=n_tiles)) / den_s
        Rr = (W @ np.bincount(tiles[r_cells], weights=R[r_cells], minlength=n_tiles)) / den_r
        with np.errstate(invalid="ignore", divide="ignore"):
            E = (W @ np.bincount(tiles[rows[e]], weights=prods, minlength=n_tiles)) / den_e
            Q = E / (Ls * Rr)
        ok = np.isfinite(Q) & (den_e > 0)
        Qb = Q[ok]
        E_obs = prods.mean() if len(e) else np.nan
        M_obs = L[s_cells].mean() * R[r_cells].mean()
        Q_obs = E_obs / M_obs if M_obs > 0 else np.nan
        lo, hi = (np.percentile(Qb, [2.5, 97.5]) if len(Qb) > 10 else (np.nan, np.nan))
        out.append({"ligand": lig, "receptor": rec, "sender": SENDER, "receiver": RECEIVER, "E_obs": E_obs, "M_obs": M_obs, "Q_obs": Q_obs,
                    "Q_ci_low": lo, "Q_ci_high": hi, "Q_boot_median": float(np.median(Qb)) if len(Qb) else np.nan,
                    "logQ_boot_se": float(np.std(np.log(Qb[Qb > 0]), ddof=1)) if (Qb > 0).sum() > 2 else np.nan,
                    "n_edges": int(len(e)), "n_sender": int(len(s_cells)), "n_receiver": int(len(r_cells)), "n_tiles": int(n_tiles),
                    "n_boot": B, "n_boot_dropped": int(B - ok.sum()), "tile_um": TILE_UM, "ci_excludes_1": bool(np.isfinite(lo) and (lo > 1 or hi < 1))})
    return out


def run_section(section, log, od):
    outp = od / f"{section}_block_bootstrap.csv"
    if outp.exists():
        log.info("section %s: exists (resume)", section)
        return
    t0 = time.time()
    adata, meta, labs, pairs, det = vn.load_section(section, log)
    xy = np.asarray(adata.obsm["spatial"], float)
    rows, cols = vn.knn_edges(xy, 6)
    nonbg = vn.nonbg_axes()
    genes = sorted({g for p in nonbg for g in p})
    X = ca.gene_matrix(adata, genes)
    gi = {g: k for k, g in enumerate(genes)}
    del adata
    gc.collect()
    rng = np.random.default_rng(SEED)
    recs = []
    for a in vn.ANNOTATIONS:
        types, lab = vn.encode(labs[a].astype(str).to_numpy())
        if SENDER not in types or RECEIVER not in types:
            log.warning("section %s %s: missing %s/%s", section, a, SENDER, RECEIVER)
            continue
        axes_list = nonbg if a == "consensus_T" else [HEAD]
        for r in block_bootstrap(xy, lab, types, rows, cols, X, gi, axes_list, np.random.default_rng(rng.integers(1 << 31))):
            recs.append({"section": section, "treatment_status": meta["treatment_status"], "annotation": a, **r})
    df = pd.DataFrame(recs)
    df.to_csv(outp, index=False)
    register(outp, PHASE, f"Block-bootstrap (100 um, B=1000) Q CI for section {section}: headline x 3 annotations + non-bg axes (consensus_T)")
    h = df[(df.annotation == "consensus_T") & (df.ligand == HEAD[0]) & (df.receptor == HEAD[1])].iloc[0]
    log.info("section %s: headline consensus_T Q=%.3f [%.3f, %.3f] edges %d (dropped %d) in %.1fs; %s", section, h.Q_obs, h.Q_ci_low, h.Q_ci_high,
             h.n_edges, h.n_boot_dropped, time.time() - t0, mem_status())


def mixed_models(hl, tmap, log):
    rows = []
    for a in vn.ANNOTATIONS:
        d = hl[(hl.annotation == a) & (hl.Q_obs > 0) & np.isfinite(hl.Q_obs)].copy()
        d["tumour"] = d.section.map(tmap)
        y = np.log(d.Q_obs.to_numpy())
        rec = {"annotation": a, "n_sections": len(d), "n_tumours": d.tumour.nunique(),
               "n_singleton_tumours": int((d.tumour.value_counts() == 1).sum()), "response": "log Q_obs (headline)"}
        t = sps.ttest_1samp(y, 0.0)
        ci = t.confidence_interval(0.95)
        rec.update({"ols_mean_logQ": float(y.mean()), "ols_ci_low": float(ci.low), "ols_ci_high": float(ci.high), "ols_p_two_sided": float(t.pvalue)})
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                md = sm.MixedLM(y, np.ones((len(y), 1)), groups=d.tumour.to_numpy())
                fit = md.fit(reml=True)
                warn = "; ".join(sorted({str(x.message)[:80] for x in w}))
            b, se = float(fit.params[0]), float(fit.bse[0])
            rec.update({"mm_intercept_logQ": b, "mm_se": se, "mm_ci_low": b - 1.96 * se, "mm_ci_high": b + 1.96 * se,
                        "mm_geomean_Q": float(np.exp(b)), "mm_geomean_Q_ci_low": float(np.exp(b - 1.96 * se)), "mm_geomean_Q_ci_high": float(np.exp(b + 1.96 * se)),
                        "mm_tumour_var": float(np.asarray(fit.cov_re).ravel()[0]), "mm_resid_var": float(fit.scale),
                        "mm_converged": bool(fit.converged), "mm_warnings": warn, "mm_p_two_sided": float(fit.pvalues[0])})
        except Exception as exc:  # pragma: no cover
            rec.update({"mm_intercept_logQ": np.nan, "mm_warnings": f"FAILED: {exc}"})
        rec["label"] = "exploratory (16 sections, 13 tumours)"
        rows.append(rec)
        log.info("MIXED %s: n=%d sections/%d tumours; mm geomean Q %.3f [%.3f, %.3f] (converged %s; tumour var %.4f; %s); t-interval mean logQ %.3f [%.3f, %.3f]",
                 a, rec["n_sections"], rec["n_tumours"], rec.get("mm_geomean_Q", np.nan), rec.get("mm_geomean_Q_ci_low", np.nan),
                 rec.get("mm_geomean_Q_ci_high", np.nan), rec.get("mm_converged"), rec.get("mm_tumour_var", np.nan), rec.get("mm_warnings", ""),
                 rec["ols_mean_logQ"], rec["ols_ci_low"], rec["ols_ci_high"])
    return pd.DataFrame(rows)


def inflation(log):
    allt = pd.read_csv(RESULTS_PDAC / "v3_phase1" / "nonbg_catalogue_by_null.csv.gz")
    g = allt[(allt["null"] == "global") & (allt.edge_status == "evaluated") & (allt.n_qualifying_edges > 0)].copy()
    for k in ("edge", "matched", "ratio"):
        g[f"{k}_z_per_sqrt_edges"] = g[f"{k}_z"] / np.sqrt(g.n_qualifying_edges)
    g["is_headline"] = (g.sender == SENDER) & (g.receiver == RECEIVER) & (g.ligand == HEAD[0]) & (g.receptor == HEAD[1])
    keep = ["section", "annotation", "sender", "receiver", "ligand", "receptor", "n_qualifying_edges", "is_headline"] + \
           [f"{k}_{s}" for k in ("edge", "matched", "ratio") for s in ("obs", "z", "z_per_sqrt_edges", "p_plus1")]
    tab = g[keep].sort_values(["annotation", "section", "ligand", "receptor", "sender", "receiver"])
    save_csv(tab, PHASE, "z_vs_edge_count.csv", "Global-null rows (7 non-bg axes, all sender/receiver, 3 annotations): raw z next to z/sqrt(n_edges); headline flagged", log)
    rows = []
    for a, d in g.groupby("annotation"):
        rec = {"annotation": a, "n_rows": len(d), "scope": "7 non-background axes, all sender->receiver, 16 sections"}
        for k in ("edge", "matched", "ratio"):
            ok = np.isfinite(d[f"{k}_z"])
            rec[f"spearman_{k}_z_vs_edges"] = float(spearmanr(d.loc[ok, f"{k}_z"], d.loc[ok, "n_qualifying_edges"]).statistic)
            rec[f"spearman_{k}_z_per_sqrt_vs_edges"] = float(spearmanr(d.loc[ok, f"{k}_z_per_sqrt_edges"], d.loc[ok, "n_qualifying_edges"]).statistic)
            rec[f"spearman_abs_{k}_z_vs_edges"] = float(spearmanr(d.loc[ok, f"{k}_z"].abs(), d.loc[ok, "n_qualifying_edges"]).statistic)
            sl = np.polyfit(np.log10(d.loc[ok, "n_qualifying_edges"]), d.loc[ok, f"{k}_z"].abs(), 1)[0]
            rec[f"slope_abs_{k}_z_per_log10_edges"] = float(sl)
        h = d[d.is_headline]
        rec["headline_n_sections"] = len(h)
        rec["headline_median_edge_z"] = float(h.edge_z.median())
        rec["headline_median_edge_z_per_sqrt"] = float(h.edge_z_per_sqrt_edges.median())
        rec["headline_median_ratio_z"] = float(h.ratio_z.median())
        rec["headline_median_ratio_z_per_sqrt"] = float(h.ratio_z_per_sqrt_edges.median())
        rec["headline_spearman_edge_z_vs_edges"] = float(spearmanr(h.edge_z, h.n_qualifying_edges).statistic) if len(h) > 3 else np.nan
        rec["headline_spearman_ratio_z_vs_edges"] = float(spearmanr(h.ratio_z, h.n_qualifying_edges).statistic) if len(h) > 3 else np.nan
        rows.append(rec)
        log.info("INFLATION %s: rho(zE, edges) %.3f -> normalised %.3f; rho(|zE|, edges) %.3f; rho(zQ, edges) %.3f -> %.3f; headline rho(zE,edges) %.3f rho(zQ,edges) %.3f",
                 a, rec["spearman_edge_z_vs_edges"], rec["spearman_edge_z_per_sqrt_vs_edges"], rec["spearman_abs_edge_z_vs_edges"],
                 rec["spearman_ratio_z_vs_edges"], rec["spearman_ratio_z_per_sqrt_vs_edges"], rec["headline_spearman_edge_z_vs_edges"], rec["headline_spearman_ratio_z_vs_edges"])
    summ = pd.DataFrame(rows)
    save_csv(summ, PHASE, "z_vs_edge_count_summary.csv", "Spearman of z (raw and /sqrt(n_edges)) with edge count per score and annotation; headline medians", log)
    return g, summ


def aggregate(log, od):
    tmap = vn.section_to_tumour()
    bb = pd.concat([pd.read_csv(f) for f in sorted(od.glob("*_block_bootstrap.csv"))], ignore_index=True)
    bb["tumour"] = bb.section.map(tmap)
    hl = bb[(bb.ligand == HEAD[0]) & (bb.receptor == HEAD[1])].copy()
    save_csv(hl.sort_values(["annotation", "section"]), PHASE, "Q_block_bootstrap_by_section.csv",
             "Headline Q per section x annotation with 100 um block-bootstrap 95% CI (B=1000), edge and cell counts", log)
    nb = bb[bb.annotation == "consensus_T"].copy()
    save_csv(nb.sort_values(["ligand", "receptor", "section"]), PHASE, "Q_block_bootstrap_nonbg_axes.csv",
             "7 non-background axes, CAF->T, consensus_T: Q per section with block-bootstrap CI", log)
    summ = hl.groupby("annotation").agg(n_sections=("section", "count"), n_ci_excl_1=("ci_excludes_1", "sum"),
                                        n_ci_above_1=("Q_ci_low", lambda x: int((x > 1).sum())), n_ci_below_1=("Q_ci_high", lambda x: int((x < 1).sum())),
                                        median_Q=("Q_obs", "median"), median_ci_width=("Q_ci_high", lambda x: np.nan),
                                        median_edges=("n_edges", "median"), total_dropped=("n_boot_dropped", "sum")).reset_index()
    summ["median_ci_width"] = [float((hl[hl.annotation == a].Q_ci_high - hl[hl.annotation == a].Q_ci_low).median()) for a in summ.annotation]
    nbs = nb.groupby(["ligand", "receptor"]).agg(n_sections=("section", "count"), n_ci_above_1=("Q_ci_low", lambda x: int((x > 1).sum())),
                                                 n_ci_below_1=("Q_ci_high", lambda x: int((x < 1).sum())), median_Q=("Q_obs", "median"),
                                                 median_edges=("n_edges", "median")).reset_index()
    save_csv(pd.concat([summ.assign(ligand=HEAD[0], receptor=HEAD[1], scope="headline by annotation"),
                        nbs.assign(annotation="consensus_T", scope="non-bg axes CAF->T")], ignore_index=True), PHASE, "Q_block_bootstrap_summary.csv",
             "Sections whose block-bootstrap 95% CI for Q lies above / below 1, per annotation (headline) and per non-bg axis (consensus_T)", log)
    for r in summ.itertuples():
        log.info("BOOT %s: %d sections; CI above 1 in %d, below 1 in %d; median Q %.3f; median CI width %.3f; median edges %.0f",
                 r.annotation, r.n_sections, r.n_ci_above_1, r.n_ci_below_1, r.median_Q, r.median_ci_width, r.median_edges)
    for r in nbs.itertuples():
        log.info("BOOT nonbg %s-%s CAF->T consensus_T: CI above 1 in %d/%d, below 1 in %d; median Q %.3f", r.ligand, r.receptor, r.n_ci_above_1, r.n_sections, r.n_ci_below_1, r.median_Q)
    mm = mixed_models(hl, tmap, log)
    save_csv(mm, PHASE, "mixed_model_logQ.csv", "Exploratory tumour-level mixed model of log Q (random intercept per tumour) and one-sample t-interval, per annotation", log)
    g, infl = inflation(log)
    save_json({"headline": HEAD, "tile_um": TILE_UM, "B": B, "seed": SEED, "section_to_tumour": tmap,
               "bootstrap_summary": summ.to_dict(orient="records"), "mixed_model": mm.to_dict(orient="records"),
               "inflation": infl.to_dict(orient="records")}, PHASE, "phase3_summary.json", "Phase 3 summary (bootstrap, mixed model, inflation)", log)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    ax = axes[0]
    h = hl[hl.annotation == "consensus_T"].sort_values("Q_obs").reset_index(drop=True)
    ax.errorbar(range(len(h)), h.Q_obs, yerr=[h.Q_obs - h.Q_ci_low, h.Q_ci_high - h.Q_obs], fmt="none", ecolor="#999999", lw=1)
    ax.scatter(range(len(h)), h.Q_obs, s=10 + 60 * h.n_edges / h.n_edges.max(), c="#d62728", zorder=3)
    ax.axhline(1, color="k", ls="--", lw=0.7)
    ax.set_xticks(range(len(h)))
    ax.set_xticklabels([f"{s} (n={n})" for s, n in zip(h.section, h.n_edges)], rotation=90, fontsize=7)
    ax.set_ylabel("Q = E/M with 95% block-bootstrap CI (100 um)")
    ax.set_title("Headline CXCL12→CXCR4 CAF→T, consensus_T\n(point size ∝ CAF→T edge count)", fontsize=9)
    ax = axes[1]
    for k, r in enumerate(mm.itertuples()):
        if np.isfinite(getattr(r, "mm_geomean_Q", np.nan)):
            ax.errorbar([k - 0.1], [r.mm_geomean_Q], yerr=[[r.mm_geomean_Q - r.mm_geomean_Q_ci_low], [r.mm_geomean_Q_ci_high - r.mm_geomean_Q]], fmt="o", color="#d62728",
                        label="mixed model (1|tumour)" if k == 0 else None)
        ax.errorbar([k + 0.1], [np.exp(r.ols_mean_logQ)], yerr=[[np.exp(r.ols_mean_logQ) - np.exp(r.ols_ci_low)], [np.exp(r.ols_ci_high) - np.exp(r.ols_mean_logQ)]], fmt="s",
                    color="#7f7f7f", label="t-interval over sections" if k == 0 else None)
    ax.axhline(1, color="k", ls="--", lw=0.7)
    ax.set_xticks(range(len(mm)))
    ax.set_xticklabels(mm.annotation, fontsize=8)
    ax.set_ylabel("geometric-mean Q (exp of log-Q intercept) with 95% CI")
    ax.set_title("Tumour-level model of log Q (exploratory; 16 sections, 13 tumours)", fontsize=9)
    ax.legend(fontsize=7)
    ax = axes[2]
    d = g[g.annotation == "consensus_T"]
    ax.scatter(d.n_qualifying_edges, d.edge_z, s=4, c="#1f77b4", alpha=0.4, label="z_E (raw)")
    ax.scatter(d.n_qualifying_edges, d.edge_z_per_sqrt_edges, s=4, c="#ff7f0e", alpha=0.4, label="z_E / sqrt(n_edges)")
    hh = d[d.is_headline]
    ax.scatter(hh.n_qualifying_edges, hh.edge_z, s=30, facecolors="none", edgecolors="k", label="headline sections (raw)")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("n qualifying edges (sender→receiver)")
    ax.set_ylabel("z (global null)")
    ax.set_title("z vs edge count, 7 non-bg axes, consensus_T\n" + f"rho(z_E, n) = {infl[infl.annotation == 'consensus_T'].spearman_edge_z_vs_edges.iloc[0]:.2f}; "
                 f"normalised {infl[infl.annotation == 'consensus_T'].spearman_edge_z_per_sqrt_vs_edges.iloc[0]:.2f}", fontsize=9)
    ax.legend(fontsize=7)
    fig.suptitle("v3 Phase 3: effect sizes instead of z (audit/diagnostic only)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase3_effect_sizes", "Headline Q with block-bootstrap CI per section; mixed-model geometric-mean Q; z vs edge count", log)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()
    log = start_run("v3_phase3_effect_sizes")
    log.info("tile_um=%.0f B=%d seed=%d; %s", TILE_UM, B, SEED, mem_status())
    od = out_dir(f"{PHASE}/sections")
    if not args.aggregate_only:
        full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
        sections = sorted(full.obs["sample_id"].astype(str).unique())
        del full
        if args.sections:
            sections = [s for s in sections if s in args.sections]
        for s in sections:
            run_section(s, log, od)
    aggregate(log, od)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

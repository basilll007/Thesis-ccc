"""v3 Phase 6: synthesis. Collects every number quoted in docs/PDAC_SUMMARY_v3.md into results/pdac/v3_phase6/v3_key_numbers.csv
(value + source file + how it was read), writes the claims register, and runs the manifest completeness check for v3 files.
Nothing here recomputes a statistic; it only reads registered Phase 1-5 outputs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, REPO_ROOT, finish_run, rel, save_csv, save_json, start_run  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "v3_phase6"
P1, P2, P3 = RESULTS_PDAC / "v3_phase1", RESULTS_PDAC / "v3_phase2", RESULTS_PDAC / "v3_phase3"
NULLS = ["global", "block200", "block100", "block50", "dist"]


def main() -> None:
    log = start_run("v3_phase6_synthesis")
    K: list[dict] = []

    def add(name, value, source, how):
        K.append({"key": name, "value": value, "source_file": rel(source), "how": how})

    # ---- Phase 1 ------------------------------------------------------------------------------------------
    hs = pd.read_csv(P1 / "headline_summary_by_null.csv")
    ts = pd.read_csv(P1 / "tumour_level_summary_by_null.csv")
    cs = pd.read_csv(P1 / "nonbg_summary_by_null.csv")
    p1 = json.loads((P1 / "phase1_summary.json").read_text())
    for a in ("consensus_T", "marker", "author_mapped"):
        for n in NULLS:
            r = hs[(hs.annotation == a) & (hs["null"] == n)]
            if len(r):
                r = r.iloc[0]
                add(f"p1.headline.{a}.{n}.n_Q_p05", int(r.n_Q_p05), P1 / "headline_summary_by_null.csv", "n_Q_p05")
                add(f"p1.headline.{a}.{n}.n_Q_bh05", int(r.n_Q_bh05), P1 / "headline_summary_by_null.csv", "n_Q_bh05")
                add(f"p1.headline.{a}.{n}.median_zQ", round(float(r.median_z_Q), 3), P1 / "headline_summary_by_null.csv", "median_z_Q")
                add(f"p1.headline.{a}.{n}.median_zE", round(float(r.median_z_E), 2), P1 / "headline_summary_by_null.csv", "median_z_E")
                add(f"p1.headline.{a}.{n}.n_zQ_pos", int(r.n_zQ_pos), P1 / "headline_summary_by_null.csv", "n_zQ_pos (of n_sections_evaluable)")
                add(f"p1.headline.{a}.{n}.n_sections", int(r.n_sections_evaluable), P1 / "headline_summary_by_null.csv", "n_sections_evaluable")
                add(f"p1.headline.{a}.{n}.median_null_mean_Q", round(float(r.median_null_mean_Q), 3), P1 / "headline_summary_by_null.csv", "median_null_mean_Q")
            t = ts[(ts.annotation == a) & (ts["null"] == n) & (ts.score == "Q")]
            if len(t):
                t = t.iloc[0]
                add(f"p1.tumour.{a}.{n}.n_tumours", int(t.n_tumours), P1 / "tumour_level_summary_by_null.csv", "n_tumours")
                add(f"p1.tumour.{a}.{n}.n_sig_bh", int(t.n_sig_bh), P1 / "tumour_level_summary_by_null.csv", "n_sig_bh (Q)")
                add(f"p1.tumour.{a}.{n}.n_sig_uncorrected", int(t.n_sig_uncorrected), P1 / "tumour_level_summary_by_null.csv", "n_sig_uncorrected (Q)")
                add(f"p1.tumour.{a}.{n}.n_mean_z_pos", int(t.n_mean_z_pos), P1 / "tumour_level_summary_by_null.csv", "n_mean_z_pos (Q)")
                add(f"p1.tumour.{a}.{n}.median_tumour_mean_zQ", round(float(t.median_tumour_mean_z), 3), P1 / "tumour_level_summary_by_null.csv", "median_tumour_mean_z (Q)")
                add(f"p1.tumour.{a}.{n}.wilcoxon_p", round(float(t.wilcoxon_z_gt0_p), 4), P1 / "tumour_level_summary_by_null.csv", "wilcoxon_z_gt0_p (Q)")
                add(f"p1.tumour.{a}.{n}.n_Q_gt1", int(t.n_mean_obs_gt1), P1 / "tumour_level_summary_by_null.csv", "n_mean_obs_gt1 (Q)")
                add(f"p1.tumour.{a}.{n}.min_q_bh", round(float(t.min_q_bh), 3), P1 / "tumour_level_summary_by_null.csv", "min_q_bh (Q)")
            c = cs[(cs.annotation == a) & (cs["null"] == n)]
            if len(c):
                c = c.iloc[0]
                add(f"p1.catalogue.{a}.{n}.n_rows", int(c.n_rows_evaluated), P1 / "nonbg_summary_by_null.csv", "n_rows_evaluated")
                add(f"p1.catalogue.{a}.{n}.frac_Q_p05", round(float(c.frac_Q_p05), 3), P1 / "nonbg_summary_by_null.csv", "frac_Q_p05")
                add(f"p1.catalogue.{a}.{n}.frac_E_p05", round(float(c.frac_E_p05), 3), P1 / "nonbg_summary_by_null.csv", "frac_E_p05")
                add(f"p1.catalogue.{a}.{n}.n_Q_bh05", int(c.n_Q_bh05), P1 / "nonbg_summary_by_null.csv", "n_Q_bh05")
                add(f"p1.catalogue.{a}.{n}.rho_zQ_vs_global", round(float(c.spearman_zQ_vs_global), 3), P1 / "nonbg_summary_by_null.csv", "spearman_zQ_vs_global")
                add(f"p1.catalogue.{a}.{n}.rho_zE_vs_global", round(float(c.spearman_zE_vs_global), 3), P1 / "nonbg_summary_by_null.csv", "spearman_zE_vs_global")
        v = p1["shrink_verdict"].get(a)
        if v:
            add(f"p1.shrink_verdict.{a}", v["shrinks_toward_zero"], P1 / "phase1_summary.json", "shrink_verdict.shrinks_toward_zero (rule in plan)")
    add("p1.n_skipped_runs", p1["n_skipped_runs"], P1 / "phase1_summary.json", "n_skipped_runs")

    # ---- Phase 2 ------------------------------------------------------------------------------------------
    cal = pd.read_csv(P2 / "calibration_deliverable.csv")
    tab = pd.read_csv(P2 / "calibration_table.csv")
    for r in cal.to_dict(orient="records"):
        for col in cal.columns:
            if col.startswith("proximity_f") or col.startswith("zonation_f"):
                add(f"p2.pooled.{r['score']}.{r['null']}.{col}", round(float(r[col]), 2), P2 / "calibration_deliverable.csv", col)
        add(f"p2.pooled.{r['score']}.{r['null']}.min_effect_power_ge_0.8", r["min_effect_power_ge_0.8"], P2 / "calibration_deliverable.csv", "min_effect_power_ge_0.8")
        add(f"p2.pooled.{r['score']}.{r['null']}.fooled_by_zonation_f1.5", bool(r["fooled_by_zonation_f1.5"]), P2 / "calibration_deliverable.csv", "fooled_by_zonation_f1.5")
    per = tab[(tab.section != "pooled_3_sections") & (tab.score == "Q")]
    for r in per.itertuples():
        add(f"p2.section.{r.section}.Q.{r.null}.{r.scenario}.f{r.effect}.rate", round(float(r.rate_p05), 2), P2 / "calibration_table.csv", "rate_p05")
    for f in sorted(P2.glob("sections/*_setup.json")):
        s = json.loads(f.read_text())
        add(f"p2.setup.{s['section']}.n_T", s["n_T"], f, "n_T")
        add(f"p2.setup.{s['section']}.n_CAF_T_edges", s["n_CAF_T_edges"], f, "n_CAF_T_edges")
        add(f"p2.setup.{s['section']}.n_targets", s["n_proximity_targets"], f, "n_proximity_targets")
        add(f"p2.setup.{s['section']}.mean_CXCR4_count_targets", round(s["mean_CXCR4_count_targets"], 3), f, "mean_CXCR4_count_targets")
        add(f"p2.setup.{s['section']}.n_sims", s["n_sims"], f, "n_sims")

    # ---- Phase 3 ------------------------------------------------------------------------------------------
    bs = pd.read_csv(P3 / "Q_block_bootstrap_summary.csv")
    for r in bs.itertuples():
        tag = f"{r.annotation}" if r.scope.startswith("headline") else f"{r.ligand}_{r.receptor}"
        add(f"p3.boot.{r.scope.split()[0]}.{tag}.n_ci_above_1", int(r.n_ci_above_1), P3 / "Q_block_bootstrap_summary.csv", "n_ci_above_1")
        add(f"p3.boot.{r.scope.split()[0]}.{tag}.n_ci_below_1", int(r.n_ci_below_1), P3 / "Q_block_bootstrap_summary.csv", "n_ci_below_1")
        add(f"p3.boot.{r.scope.split()[0]}.{tag}.n_sections", int(r.n_sections), P3 / "Q_block_bootstrap_summary.csv", "n_sections")
        add(f"p3.boot.{r.scope.split()[0]}.{tag}.median_Q", round(float(r.median_Q), 3), P3 / "Q_block_bootstrap_summary.csv", "median_Q")
        if np.isfinite(r.median_ci_width):
            add(f"p3.boot.{r.scope.split()[0]}.{tag}.median_ci_width", round(float(r.median_ci_width), 3), P3 / "Q_block_bootstrap_summary.csv", "median_ci_width")
    hb = pd.read_csv(P3 / "Q_block_bootstrap_by_section.csv")
    for r in hb[hb.annotation == "consensus_T"].itertuples():
        add(f"p3.boot.section.{r.section}.Q", f"{r.Q_obs:.3f} [{r.Q_ci_low:.3f}, {r.Q_ci_high:.3f}] edges={r.n_edges}", P3 / "Q_block_bootstrap_by_section.csv", "Q_obs, Q_ci_low, Q_ci_high, n_edges (consensus_T)")
    mm = pd.read_csv(P3 / "mixed_model_logQ.csv")
    for r in mm.itertuples():
        add(f"p3.mm.{r.annotation}.geomean_Q", round(float(r.mm_geomean_Q), 3), P3 / "mixed_model_logQ.csv", "mm_geomean_Q")
        add(f"p3.mm.{r.annotation}.ci", f"[{r.mm_geomean_Q_ci_low:.3f}, {r.mm_geomean_Q_ci_high:.3f}]", P3 / "mixed_model_logQ.csv", "mm_geomean_Q_ci_low/high")
        add(f"p3.mm.{r.annotation}.n_tumours", int(r.n_tumours), P3 / "mixed_model_logQ.csv", "n_tumours")
        add(f"p3.mm.{r.annotation}.n_sections", int(r.n_sections), P3 / "mixed_model_logQ.csv", "n_sections")
        add(f"p3.mm.{r.annotation}.tumour_var", round(float(r.mm_tumour_var), 4), P3 / "mixed_model_logQ.csv", "mm_tumour_var")
        add(f"p3.mm.{r.annotation}.converged", bool(r.mm_converged), P3 / "mixed_model_logQ.csv", "mm_converged")
        add(f"p3.mm.{r.annotation}.warnings", str(r.mm_warnings), P3 / "mixed_model_logQ.csv", "mm_warnings")
        add(f"p3.t.{r.annotation}.geomean_Q_ci", f"{np.exp(r.ols_mean_logQ):.3f} [{np.exp(r.ols_ci_low):.3f}, {np.exp(r.ols_ci_high):.3f}]", P3 / "mixed_model_logQ.csv", "exp(ols_mean_logQ, ols_ci_low, ols_ci_high)")
    inf = pd.read_csv(P3 / "z_vs_edge_count_summary.csv")
    for r in inf.itertuples():
        add(f"p3.infl.{r.annotation}.rho_zE_edges", round(float(r.spearman_edge_z_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "spearman_edge_z_vs_edges")
        add(f"p3.infl.{r.annotation}.rho_zE_norm_edges", round(float(r.spearman_edge_z_per_sqrt_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "spearman_edge_z_per_sqrt_vs_edges")
        add(f"p3.infl.{r.annotation}.rho_absZE_edges", round(float(r.spearman_abs_edge_z_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "spearman_abs_edge_z_vs_edges")
        add(f"p3.infl.{r.annotation}.rho_zQ_edges", round(float(r.spearman_ratio_z_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "spearman_ratio_z_vs_edges")
        add(f"p3.infl.{r.annotation}.rho_zQ_norm_edges", round(float(r.spearman_ratio_z_per_sqrt_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "spearman_ratio_z_per_sqrt_vs_edges")
        add(f"p3.infl.{r.annotation}.n_rows", int(r.n_rows), P3 / "z_vs_edge_count_summary.csv", "n_rows")
        add(f"p3.infl.{r.annotation}.headline_median_zE", round(float(r.headline_median_edge_z), 2), P3 / "z_vs_edge_count_summary.csv", "headline_median_edge_z")
        add(f"p3.infl.{r.annotation}.headline_median_zE_norm", round(float(r.headline_median_edge_z_per_sqrt), 3), P3 / "z_vs_edge_count_summary.csv", "headline_median_edge_z_per_sqrt")
        add(f"p3.infl.{r.annotation}.headline_rho_zE_edges", round(float(r.headline_spearman_edge_z_vs_edges), 3), P3 / "z_vs_edge_count_summary.csv", "headline_spearman_edge_z_vs_edges")

    # ---- gates ---------------------------------------------------------------------------------------------
    for ph in ("v3_phase4", "v3_phase5"):
        f = RESULTS_PDAC / ph / f"{ph.replace('v3_', '')}_gate_status.json"
        if f.exists():
            s = json.loads(f.read_text())
            add(f"{ph}.status", s["status"], f, "status")
            add(f"{ph}.reason", s["reason"], f, "reason")

    kn = pd.DataFrame(K)
    save_csv(kn, PHASE, "v3_key_numbers.csv", "Every number quoted in docs/PDAC_SUMMARY_v3.md with its source file and column", log)
    log.info("%d key numbers collected", len(kn))

    # ---- manifest check for v3 files -------------------------------------------------------------------------
    man = pd.read_csv(RESULTS_PDAC / "MANIFEST.csv")
    issues = []
    for p in man.path:
        if not (REPO_ROOT / p).exists():
            issues.append({"path": p, "issue": "registered but missing on disk"})
    reg = set(man.path)
    for d in ("results/pdac/v3_phase1", "results/pdac/v3_phase2", "results/pdac/v3_phase3", "results/pdac/v3_phase4", "results/pdac/v3_phase5",
              "results/pdac/v3_phase6", "figures/pdac/v3_phase1", "figures/pdac/v3_phase2", "figures/pdac/v3_phase3"):
        for f in sorted((REPO_ROOT / d).rglob("*")):
            if f.is_file() and rel(f) not in reg and f.name != "v3_manifest_check.csv":
                issues.append({"path": rel(f), "issue": "on disk but not registered"})
    for f in sorted((REPO_ROOT / "logs").glob("v3_*.log")):
        if rel(f) not in reg:
            issues.append({"path": rel(f), "issue": "run log not registered"})
    chk = pd.DataFrame(issues, columns=["path", "issue"])
    save_csv(chk, PHASE, "v3_manifest_check.csv", "Manifest completeness check for v3 outputs (missing files, unregistered files)", log)
    log.info("manifest check: %d v3-related issues; manifest rows %d", len(chk), len(man))
    save_json({"n_key_numbers": len(kn), "n_manifest_issues": len(chk), "manifest_rows": int(len(man))}, PHASE, "phase6_summary.json", "Phase 6 synthesis summary", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

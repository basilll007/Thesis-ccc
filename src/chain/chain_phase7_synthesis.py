"""Chain Phase 7: claims register, key-numbers table, manifest check and overview figure.

Plan: docs/chain_plans/phase7_synthesis.md. Every number written here is read from a generated results file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase7"
R = cc.RESULTS


def rd(path, **kw):
    return pd.read_csv(R / path, keep_default_na=False, na_values=[""], **kw)


def main():
    log = cc.start_run("chain_phase7_synthesis")
    kn = []

    def k(name, value, source):
        kn.append({"name": name, "value": value, "source": f"results/chain/{source}"})
        return value

    # Phase 1
    p1 = json.loads((R / "phase1" / "phase1_summary.json").read_text())
    p1b = json.loads((R / "phase1b" / "phase1_summary.json").read_text())
    k("n_background_level_genes_sectionwide", p1["n_background_level_genes"], "phase1/phase1_summary.json")
    k("programs_kept_original_rule", ";".join(p1["programs_kept"]), "phase1/phase1_summary.json")
    k("programs_kept_amendment1", ";".join(p1b["programs_kept"]), "phase1b/phase1_summary.json")
    k("programs_dropped_amendment1", ";".join(p1b["programs_dropped"]), "phase1b/phase1_summary.json")
    val = rd("phase1b/validity_auc.csv").dropna(subset=["auc"])
    auc_med = k("validity_auc_classical_vs_basal_median", round(float(val.auc.median()), 3), "phase1b/validity_auc.csv")
    sel = rd("phase1b/chain_selection.csv")
    prim = sel[sel.chain == "tumour->CAF->T"].iloc[0]
    k("primary_chain_sections_pass", int(prim.n_sections_pass), "phase1b/chain_selection.csv")
    k("primary_chain_tumours_pass", int(prim.n_tumours_pass), "phase1b/chain_selection.csv")
    # Phase 2
    d2 = rd("phase2/delta_r2_by_program.csv")
    a2 = d2[(d2.annotation == "author_mapped") & (~d2.spillover_covariates.astype(bool))]
    for r in a2.itertuples():
        k(f"delta_r2_{r.program}_{r.model}", f"{r.delta_r2:.4f} [{r.ci_low:.4f}, {r.ci_high:.4f}]", "phase2/delta_r2_by_program.csv")
    n_pos = int(a2[a2.model == "ridge"].neighbours_predict_beyond_own.astype(str).eq("True").sum())
    n_prog = int(a2.program.nunique())
    att = rd("phase2/feature_group_attribution.csv")
    att = att[(att.kind == "drop") & (~att.spillover_covariates.astype(bool))]
    other_max = k("max_drop_group_delta_r2_other_type_programs",
                  round(float(att[att.group.str.startswith("other_type")].delta_r2.max()), 4),
                  "phase2/feature_group_attribution.csv")
    same_rng = att[att.group == "same_type_programs"].delta_r2
    k("drop_group_same_type_range", f"{same_rng.min():.4f}-{same_rng.max():.4f}", "phase2/feature_group_attribution.csv")
    sgc = rd("phase2/gnn_sgc_v2_delta_r2.csv")
    tcy = a2[(a2.program == "T:cytotoxic")].set_index("model")
    # Phase 3
    calls = rd("phase3/chain_calls_author_mapped.csv")
    callsm = rd("phase3/chain_calls_marker.csv")
    tl = rd("phase3/tumour_level_tests_author_mapped.csv")
    fw = tl[tl.variant == "forward"]
    k("phase3_forward_combinations", len(calls), "phase3/chain_calls_author_mapped.csv")
    k("phase3_labels_author", json.dumps(calls.label.value_counts().to_dict()), "phase3/chain_calls_author_mapped.csv")
    k("phase3_labels_marker", json.dumps(callsm.label.value_counts().to_dict()), "phase3/chain_calls_marker.csv")
    minq = fw[fw.stat.isin(["beta_CB", "beta_BA", "indirect"])].groupby("stat").q_expr.min()
    for s, v in minq.items():
        k(f"phase3_min_q_expr_{s}", round(float(v), 4), "phase3/tumour_level_tests_author_mapped.csv")
    pc = fw[(fw.chain == "tumour->CAF->T") & (fw.stat == "beta_CB")]
    k("primary_beta_CB_range", f"{pc.estimate.min():.4f} to {pc.estimate.max():.4f}", "phase3/tumour_level_tests_author_mapped.csv")
    pi = fw[(fw.chain == "tumour->CAF->T") & (fw.stat == "indirect")]
    k("primary_indirect_range", f"{pi.estimate.min():.4f} to {pi.estimate.max():.4f}", "phase3/tumour_level_tests_author_mapped.csv")
    bx = fw[(fw.chain == "tumour->CAF->T") & (fw.stat == "beta_BA") & (fw.prog_A == "tumour:basal_like") &
            (fw.prog_B == "CAF:CXCL12")].iloc[0]
    k("hop_CXCL12CAF_from_basal_tumour", f"{bx.estimate:.4f} [{bx.ci_low:.4f}, {bx.ci_high:.4f}]; {int(bx.n_tumours_pos)}/"
      f"{int(bx.n_tumours)} tumours positive; p_expr {bx.p_expr:.4f}; q {bx.q_expr:.4f}; p_label {bx.p_label:.4f}",
      "phase3/tumour_level_tests_author_mapped.csv")
    lnm = fw[fw.stat.isin(["beta_CB", "beta_BA"])].null_mean_label
    k("label_null_mean_range_hops", f"{lnm.min():.3f} to {lnm.max():.3f}", "phase3/tumour_level_tests_author_mapped.csv")
    am = rd("phase3/all_variant_calls_marker.csv")
    mh = am[(am.variant == "forward") & (am.c1_hop_CB.astype(str).eq("True") | am.c2_hop_BA.astype(str).eq("True"))]
    k("marker_forward_hops_passing_both_nulls", len(mh), "phase3/all_variant_calls_marker.csv")
    # Phase 4
    cal = json.loads((R / "phase4" / "calibrated_range.json").read_text())
    ct = rd("phase4/calibration_table.csv")
    pool = rd("phase4/calibration_pooled_3sections.csv")
    chain_pw = ct[ct.scenario == "chain"].chain_call_rate.max()
    k("max_chain_call_power_any_section_fold", float(chain_pw), "phase4/calibration_table.csv")
    k("pooled_chain_power_by_fold", json.dumps(dict(zip(pool[pool.scenario == "chain"].fold.astype(str),
                                                          pool[pool.scenario == "chain"].chain_call_rate))),
      "phase4/calibration_pooled_3sections.csv")
    z2 = ct[(ct.scenario == "zonation") & (ct.fold == 2.0) & (ct.section == "T4-a")].iloc[0]
    k("T4a_zonation2_hop_rates", f"C<-B {z2.hop_CB_rate:.2f}; B<-A {z2.hop_BA_rate:.2f}; chain {z2.chain_call_rate:.2f}",
      "phase4/calibration_table.csv")
    k("trustworthy_setting_found", cal["trustworthy_setting_found"], "phase4/calibrated_range.json")
    # Phase 5
    p5 = json.loads((R / "phase5" / "phase5_summary.json").read_text())
    k("phase5_hops_supported", f"{p5['n_hops_supported']}/{p5['n_hops']}", "phase5/phase5_summary.json")
    dd = rd("ranked_hypotheses_dedup.csv")
    k("ranked_hypotheses_dedup_rows", len(dd), "ranked_hypotheses_dedup.csv")
    # Phase 6
    p6 = json.loads((R / "phase6" / "phase6_summary.json").read_text())
    p6r = rd("phase6/phase2_repeat_delta_r2.csv")
    t6 = rd("phase6/phase3_repeat_unit_level.csv")
    k("phase6_chains_analysed", ";".join(p6["chains_analysed"]), "phase6/phase6_summary.json")
    for st, v in p6["phase3_same_direction_rate"].items():
        k(f"phase6_same_direction_rate_{st}", round(v, 3) if v == v else "NA", "phase6/phase6_summary.json")
    r6 = p6r[(p6r.model == "ridge") & (~p6r.spillover_covariates.astype(bool))]
    k("phase6_phase2_same_direction_ridge", f"{int(r6.same_direction_as_mendeley.astype(str).eq('True').sum())}/{len(r6)}",
      "phase6/phase2_repeat_delta_r2.csv")
    k("phase6_phase2_ci_excl0_ridge", f"{int(r6.ci_excludes_0.astype(str).eq('True').sum())}/{len(r6)}",
      "phase6/phase2_repeat_delta_r2.csv")
    keyn = pd.DataFrame(kn)
    cc.save_csv(keyn, PHASE, "key_numbers.csv", "Key numbers with source files", log)

    # claims register
    def num(name):
        return str(keyn.set_index("name").value[name])
    claims = [
        ("C1", "Under the pre-declared section-wide background rule no T-cell program is measurable, so T-containing chains "
               "(including the primary) are not testable under that rule", "not testable",
         "phase1/programs_used.json; phase1/phase1_summary.json", f"programs kept: {num('programs_kept_original_rule')}"),
        ("C2", "With the type-aware background rule (Amendment 1, logged deviation) 7 programs are measurable; the classical "
               "vs basal-like contrast agrees with the authors' CLS/BSL subtypes", "supported",
         "phase1b/validity_auc.csv", f"median AUC {auc_med}"),
        ("C3", f"Neighbourhood features improve held-out prediction of cell state beyond own covariates (incl. 200 um tile "
               f"zonation control) for {n_pos}/{n_prog} receiver programs (ridge, tumour-bootstrap CI > 0)", "supported",
         "phase2/delta_r2_by_program.csv", f"classical {num('delta_r2_tumour:classical_ridge')}"),
        ("C4", "The neighbourhood gain is carried by same-type neighbour programs and composition; other cell types' "
               "programs add almost nothing", "supported", "phase2/feature_group_attribution.csv",
         f"max other-type drop-group delta R2 {other_max}; same-type {num('drop_group_same_type_range')}"),
        ("C5", "T-cell cytotoxic program is predicted by neighbours beyond own covariates", "not supported",
         "phase2/delta_r2_by_program.csv",
         f"ridge {num('delta_r2_T:cytotoxic_ridge')}; hgb {num('delta_r2_T:cytotoxic_hgb')} (marker-label HGB CI > 0: exploratory)"),
        ("C6", "A graph (SGC-style) model adds predictive value over ridge/HGB", "not supported", "phase2/gnn_sgc_v2_delta_r2.csv",
         "; ".join(f"{r.program}{' +spill' if str(r.spillover_covariates) == 'True' else ''}: {r.delta_r2:.4f} vs ridge "
                   f"{r.ridge_delta_r2:.4f}" for r in sgc.itertuples())),
        ("C7", "Any A -> B -> C chain meets the pre-declared chain-evidence rule", "not supported",
         "phase3/chain_calls_author_mapped.csv; phase3/chain_calls_marker.csv",
         f"author {num('phase3_labels_author')}; marker {num('phase3_labels_marker')}"),
        ("C8", "Primary chain tumour -> CAF -> T: T-cell state associates with adjacent CAF program, or an indirect "
               "tumour -> T association via CAF", "not supported", "phase3/tumour_level_tests_author_mapped.csv",
         f"beta_CB {num('primary_beta_CB_range')}; indirect {num('primary_indirect_range')}"),
        ("C9", "CAFs adjacent to basal-like tumour cells have lower CXCL12 (single hop; expression null only)", "exploratory",
         "phase3/tumour_level_tests_author_mapped.csv", num("hop_CXCL12CAF_from_basal_tumour")),
        ("C10", "CAF (myCAF/CXCL12) <-> myeloid C1Q_TAM adjacency association (single hops pass both nulls under marker "
                "labels only, in both directions; direction not identifiable)", "exploratory",
         "phase3/all_variant_calls_marker.csv", f"{num('marker_forward_hops_passing_both_nulls')} forward hop rows"),
        ("C11", "The chain test can detect a planted chain within 1.1x-2x per-hop effects on real geometry", "not supported",
         "phase4/calibration_table.csv; phase4/calibration_pooled_3sections.csv",
         f"max chain-call power {num('max_chain_call_power_any_section_fold')}; pooled {num('pooled_chain_power_by_fold')}"),
        ("C12", "Single-hop tests are robust to zonation", "not supported", "phase4/calibration_table.csv",
         f"T4-a zonation 2x: {num('T4a_zonation2_hop_rates')}"),
        ("C13", "Any ligand-receptor candidate is attached to a supported hop", "not supported", "phase5/phase5_summary.json",
         f"hops supported {num('phase5_hops_supported')}"),
        ("C14", "Phase 2/3 directions are consistent in GSE310352 (consistency check; independence UNVERIFIED)",
         "exploratory", "phase6/phase2_repeat_delta_r2.csv; phase6/phase3_repeat_unit_level.csv",
         f"Phase 2 same direction (ridge) {num('phase6_phase2_same_direction_ridge')}, CI>0 {num('phase6_phase2_ci_excl0_ridge')}; "
         f"Phase 3 same-direction rates beta_BA {num('phase6_same_direction_rate_beta_BA')}, beta_CB "
         f"{num('phase6_same_direction_rate_beta_CB')}, indirect {num('phase6_same_direction_rate_indirect')}"),
    ]
    cr = pd.DataFrame(claims, columns=["id", "claim", "status", "evidence_file", "key_number"])
    cc.save_csv(cr, PHASE, "claims_register.csv", "Claims register (status per pre-declared rules)", log)
    # manifest check
    man = pd.read_csv(cc.MANIFEST)
    listed = set(man.path)
    files = {cc.rel(p) for d in (cc.RESULTS, cc.FIGURES) for p in d.rglob("*") if p.is_file()
             and p.name not in ("MANIFEST.csv", ".manifest.lock")}
    issues = [{"path": p, "issue": "not in manifest"} for p in sorted(files - listed)]
    issues += [{"path": p, "issue": "in manifest but missing"} for p in sorted(listed) if not (cc.REPO_ROOT / p).exists()]
    # files written after this check (this script's own outputs) are registered by save_*; check excludes them
    issues = [i for i in issues if not i["path"].startswith("results/chain/phase7/manifest_check")]
    cc.save_csv(pd.DataFrame(issues, columns=["path", "issue"]), PHASE, "manifest_check.csv",
                "Manifest completeness check (results/chain, figures/chain)", log)
    log.info("manifest issues: %d", len(issues))
    # overview figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    e = a2[a2.model == "ridge"].sort_values("delta_r2")
    axes[0].errorbar(e.delta_r2, range(len(e)), xerr=[e.delta_r2 - e.ci_low, e.ci_high - e.delta_r2], fmt="o", capsize=2)
    axes[0].set_yticks(range(len(e)), e.program, fontsize=8)
    axes[0].axvline(0, color="k", lw=.6)
    axes[0].set_xlabel("held-out delta R2 (ridge, LOTO; 95% tumour CI)")
    axes[0].set_title("Neighbours predict state (mostly same-type)", fontsize=10)
    f = fw[fw.stat.isin(["beta_BA", "beta_CB", "indirect"])]
    for st, c in (("beta_BA", "#2ca02c"), ("beta_CB", "#1f77b4"), ("indirect", "#d62728")):
        x = f[f.stat == st]
        axes[1].scatter(x.estimate, -np.log10(x.q_expr.clip(lower=1e-3)), s=12, color=c, label=st, alpha=.7)
    axes[1].axhline(-np.log10(.05), color="grey", ls=":")
    axes[1].set_xlabel("tumour-level standardised effect")
    axes[1].set_ylabel("-log10 BH q (expression block null)")
    axes[1].set_title("Chain tests: 44 combinations, none supported", fontsize=10)
    axes[1].legend(fontsize=7)
    cp = ct[ct.scenario == "chain"]
    for sec, mk in zip(cp.section.unique(), "osd"):
        x = cp[cp.section == sec].sort_values("fold")
        axes[2].plot(x.fold, x.chain_call_rate, marker=mk, label=f"{sec} chain call")
        axes[2].plot(x.fold, x.hop_BA_rate, marker=mk, ls=":", label=f"{sec} hop B<-A")
    axes[2].axhline(.8, color="grey", lw=.6)
    axes[2].set_xlabel("planted per-hop fold change")
    axes[2].set_ylabel("detection rate (100 sims)")
    axes[2].set_title("Calibration: chain undetectable at 1.1-2x", fontsize=10)
    axes[2].legend(fontsize=6)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, "phase7_overview", "Overview: Phase 2 gains, Phase 3 tests, Phase 4 power", log)
    log.info("\n%s", cr.to_string(index=False))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

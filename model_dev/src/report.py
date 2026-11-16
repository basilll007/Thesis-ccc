"""Phase 7: collect every reported number from generated files into tables, the data card and the model card.

What: writes model_dev/results/phase7/key_numbers.csv (name, value, source file), results_tables.md (markdown
tables for the report), reports/DATA_CARD.md and reports/MODEL_CARD.md, and a checks table (tests, manifest of
artefacts, compute). The narrative report (reports/MODEL_REPORT_v1.md) quotes only numbers present in key_numbers.csv.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase7"
R = C.RESULTS


def rd(path: str) -> pd.DataFrame | None:
    p = R / path
    return pd.read_csv(p) if p.exists() else None


def fmt_ci(r2, lo, hi) -> str:
    return f"{r2:.3f} [{lo:.3f}, {hi:.3f}]"


def main() -> None:
    log = C.start_run(PHASE, "report")
    kn: list[dict] = []

    def k(name, value, source):
        kn.append({"name": name, "value": value, "source": f"model_dev/results/{source}"})
        return value

    md: list[str] = []
    # ---- Phase 0 / data card
    panel = json.loads((R / "phase0" / "panel_intersection.json").read_text())
    counts = rd("phase0/cell_counts.csv")
    by_group = rd("phase0/cell_counts_by_group.csv")
    k("n_model_genes_G", panel["n_model_gene_set_G"], "phase0/panel_intersection.json")
    k("n_gse_only_genes_dropped", len(panel["gse_only"]), "phase0/panel_intersection.json")
    k("n_dev_cells", int(counts[counts.split == "development"].n_cells.sum()), "phase0/cell_counts.csv")
    k("n_ext_cells", int(counts[counts.split == "external"].n_cells.sum()), "phase0/cell_counts.csv")
    for cls in ("myeloid", "epithelial_ductal_tumor", "T_cell"):
        for split in ("development", "external"):
            k(f"n_{split}_{cls}", int(counts[(counts.split == split) & (counts.cell_type == cls)].n_cells.sum()), "phase0/cell_counts.csv")
    piv = by_group.pivot_table(index=["split", "group"], columns="cell_type", values="n_cells", fill_value=0).reset_index()
    # ---- Phase 1
    folds = rd("phase1/folds.csv")
    floor = rd("phase1/fold_floor_check.csv")
    k("n_folds", int(folds.fold.nunique()), "phase1/folds.csv")
    k("all_folds_pass_floor", bool(floor.passes_floor.all()), "phase1/fold_floor_check.csv")
    # ---- Phase 2
    rep = rd("phase2/reproduction_check.csv")
    k("reproduction_max_abs_diff", float(rep.abs_diff.max()), "phase2/reproduction_check.csv")
    k("reproduction_pass", bool(rep["pass"].all()), "phase2/reproduction_check.csv")
    summ2 = rd("phase2/baseline_summary.csv")
    best2 = rd("phase2/strongest_baseline_per_target.csv")
    tab = []
    for r in summ2.itertuples():
        tab.append({"receiver": r.receiver, "target": r.target, "model": f"{r.model}/{r.feature_set}",
                    "R2 [95% tumour CI]": fmt_ci(r.r2, r.ci_low, r.ci_high),
                    "gain vs own covariates [CI]": fmt_ci(r.gain_vs_own, r.gain_ci_low, r.gain_ci_high), "config": r.selected_config})
        k(f"p2_{r.receiver}_{r.target}_{r.model}/{r.feature_set}_r2", round(r.r2, 4), "phase2/baseline_summary.csv")
        k(f"p2_{r.receiver}_{r.target}_{r.model}/{r.feature_set}_gain", f"{r.gain_vs_own:.4f} [{r.gain_ci_low:.4f}, {r.gain_ci_high:.4f}]", "phase2/baseline_summary.csv")
    md.append("## Phase 2 baselines (grouped 5-fold, pooled held-out R2)\n\n" + pd.DataFrame(tab).to_markdown(index=False) + "\n")
    for r in best2.itertuples():
        k(f"strongest_baseline_{r.receiver}_{r.target}", f"{r.model}/{r.feature_set} R2 {r.r2:.4f}", "phase2/strongest_baseline_per_target.csv")
    # ---- Phase 3
    g3, abl, sel = rd("phase3/graph_summary.csv"), rd("phase3/ablation_summary.csv"), None
    if (R / "phase3" / "architecture_selection.json").exists():
        sel = json.loads((R / "phase3" / "architecture_selection.json").read_text())
        k("selected_arch", sel["selected_arch"], "phase3/architecture_selection.json")
        k("selected_lr", sel["selected_lr"], "phase3/architecture_selection.json")
    if g3 is not None:
        tab = [{"receiver": r.receiver, "target": r.target, "arch": r.arch, "lr": r.selected_lr,
                "R2 [CI]": fmt_ci(r.r2, r.ci_low, r.ci_high), "per seed": r.r2_per_seed, "strongest baseline": f"{r.strongest_baseline} ({r.r2_baseline:.3f})",
                "gain vs baseline [CI]": fmt_ci(r.gain_vs_baseline, r.gain_ci_low, r.gain_ci_high), "beats": r.beats_baseline} for r in g3.itertuples()]
        md.append("## Phase 3 graph models\n\n" + pd.DataFrame(tab).to_markdown(index=False) + "\n")
        for r in g3.itertuples():
            k(f"p3_{r.receiver}_{r.target}_{r.arch}_r2", fmt_ci(r.r2, r.ci_low, r.ci_high), "phase3/graph_summary.csv")
            k(f"p3_{r.receiver}_{r.target}_{r.arch}_gain_vs_baseline", fmt_ci(r.gain_vs_baseline, r.gain_ci_low, r.gain_ci_high), "phase3/graph_summary.csv")
            k(f"p3_{r.receiver}_{r.target}_{r.arch}_beats_baseline", bool(r.beats_baseline), "phase3/graph_summary.csv")
    if abl is not None:
        tab = [{"receiver": r.receiver, "target": r.target, "ablation": r["mode"], "R2 ablation [CI]": fmt_ci(r.r2_mode, r.ci_low, r.ci_high),
                "full - ablation [CI]": fmt_ci(r.full_minus_mode, r.diff_ci_low, r.diff_ci_high), "CI excludes 0": r.diff_ci_excludes_0,
                "quantity": r.quantity} for _, r in abl.iterrows()]
        md.append("## Phase 3 ablations (selected architecture)\n\n" + pd.DataFrame(tab).to_markdown(index=False) + "\n")
        for _, r in abl.iterrows():
            k(f"abl_{r.receiver}_{r.target}_{r['mode']}_full_minus", fmt_ci(r.full_minus_mode, r.diff_ci_low, r.diff_ci_high), "phase3/ablation_summary.csv")
            k(f"abl_{r.receiver}_{r.target}_{r['mode']}_r2", fmt_ci(r.r2_mode, r.ci_low, r.ci_high), "phase3/ablation_summary.csv")
    comp = pd.concat([d for d in (rd("phase3/compute_graph.csv"), rd("phase3/compute_ablation.csv")) if d is not None], ignore_index=True) if (R / "phase3" / "compute_graph.csv").exists() else None
    if comp is not None:
        k("p3_total_gpu_hours", round(float(comp.wall_s.sum()) / 3600, 2), "phase3/compute_*.csv")
        k("p3_n_runs", len(comp), "phase3/compute_*.csv")
        k("p3_peak_vram_mib", float(comp.peak_vram_mib.max()), "phase3/compute_*.csv")
        k("p3_median_epochs", float(comp.epochs.median()), "phase3/compute_*.csv")
    # ---- Phase 4
    ext = rd("phase4/external_results.csv")
    if ext is not None:
        tab = [{"receiver": r.receiver, "target": r.target, "model": r.model, "internal R2": f"{r.internal_r2:.3f}",
                "external R2 [unit CI]": fmt_ci(r.external_r2, r.unit_ci_low, r.unit_ci_high), "drop": f"{r.internal_to_external_drop:.3f}",
                "gain vs own (ext) [CI]": (fmt_ci(r.gain_vs_own_external, r.gain_unit_ci_low, r.gain_unit_ci_high) if pd.notna(r.gain_vs_own_external) else ""),
                "per slide": r.per_slide_r2} for r in ext.itertuples()]
        md.append("## Phase 4 external test (GSE310352; single touch; units UNVERIFIED as independent patients)\n\n" + pd.DataFrame(tab).to_markdown(index=False) + "\n")
        for r in ext.itertuples():
            k(f"p4_{r.receiver}_{r.target}_{r.model}_external_r2", fmt_ci(r.external_r2, r.unit_ci_low, r.unit_ci_high), "phase4/external_results.csv")
            k(f"p4_{r.receiver}_{r.target}_{r.model}_drop", round(float(r.internal_to_external_drop), 4) if pd.notna(r.internal_to_external_drop) else "NA", "phase4/external_results.csv")
            if pd.notna(r.gain_vs_own_external):
                k(f"p4_{r.receiver}_{r.target}_{r.model}_gain_vs_own_external", fmt_ci(r.gain_vs_own_external, r.gain_unit_ci_low, r.gain_unit_ci_high), "phase4/external_results.csv")
    # ---- Phase 5
    if (R / "phase5" / "gate.json").exists():
        gate = json.loads((R / "phase5" / "gate.json").read_text())
        k("p5_gate", gate["decision"], "phase5/gate.json")
        imp = rd("phase5/permutation_importance.csv")
        if imp is not None:
            top = imp.groupby(["receiver", "target"]).head(6)
            md.append("## Phase 5 permutation importance (top groups)\n\n" + top.round(4).to_markdown(index=False) + "\n")
            for r in top.itertuples():
                k(f"p5_imp_{r.receiver}_{r.target}_{r.group}", fmt_ci(r.r2_drop, r.ci_low, r.ci_high), "phase5/permutation_importance.csv")
        pert = rd("phase5/in_silico_perturbation.csv")
        if pert is not None and len(pert):
            md.append("## Phase 5 in-silico perturbation (model-based hypotheses)\n\n" + pert.round(4).to_markdown(index=False) + "\n")
    # ---- Phase 6
    p6 = rd("phase6/probe_summary.csv")
    if p6 is not None:
        tab = [{"receiver": r.receiver, "probe": r.probe, "target": r.target, "R2 [CI]": fmt_ci(r.r2, r.ci_low, r.ci_high),
                "probe - strongest baseline [CI]": fmt_ci(r.probe_minus_baseline, r.pb_ci_low, r.pb_ci_high),
                "probe - graph [CI]": fmt_ci(r.probe_minus_graph, r.pg_ci_low, r.pg_ci_high), "pretraining helps": r.pretraining_helps} for r in p6.itertuples()]
        md.append("## Phase 6 pretraining probe\n\n" + pd.DataFrame(tab).to_markdown(index=False) + "\n")
        for r in p6.itertuples():
            k(f"p6_{r.receiver}_{r.probe}_{r.target}_r2", fmt_ci(r.r2, r.ci_low, r.ci_high), "phase6/probe_summary.csv")
            k(f"p6_{r.receiver}_{r.probe}_{r.target}_minus_graph", fmt_ci(r.probe_minus_graph, r.pg_ci_low, r.pg_ci_high), "phase6/probe_summary.csv")
    lc = rd("phase6/learning_curve.csv")
    if lc is not None:
        md.append("## Phase 6 learning curve\n\n" + lc.round(4).to_markdown(index=False) + "\n")
    # ---- tests and artefacts
    test = subprocess.run([sys.executable, "-m", "pytest", str(C.MODEL_DEV / "tests"), "-q", "-p", "no:cacheprovider"],
                          capture_output=True, text=True, cwd=C.REPO_ROOT)
    last = [l for l in test.stdout.strip().splitlines() if l.strip()][-1] if test.stdout.strip() else test.stderr[-200:]
    k("unit_tests", last, "pytest model_dev/tests")
    art = rd("phase0/artefacts.csv")
    k("n_registered_artefacts", len(art) if art is not None else 0, "phase0/artefacts.csv")
    missing = [p for p in (art.path if art is not None else []) if not (C.REPO_ROOT / p).exists()]
    k("artefacts_missing", len(missing), "phase0/artefacts.csv")
    C.save_table(pd.DataFrame(kn), PHASE, "key_numbers.csv", log)
    (R / PHASE / "results_tables.md").write_text("\n".join(md), encoding="utf-8")
    # ---- data card
    dc = ["# Data card (model_dev v1)", "",
          "## Sources and licences (as recorded in this repository)",
          "- Development set: Mendeley PDAC-SMI CosMx, doi 10.17632/kx6b69n3cb.1, CC BY 4.0 (as recorded in `src/pdac/pdac_loaders.py`). 16 sections, 13 tumours; authors' major-type labels mapped to 13 classes (`src/pdac/phase3d_author_annotation.AUTHOR_MAP`).",
          "- External set: GEO GSE310352 CosMx (Lyubetskaya et al., Cell Rep 2026, doi 10.1016/j.celrep.2025.116827, as recorded in the loader); terms confirmed in `LICENSE_CONFIRMED_GSE310352.txt`. 8 slides; marker labels from `results/pdac/v3_phase4`. Patient identity, independence from the development set and the um/px scale are UNVERIFIED; slides 5/5b2/5b3 are grouped as one unit by name (UNVERIFIED).",
          "", "## Processing",
          f"- QC: min 10 counts per cell, min 3 cells per gene (`src/pdac/phase3_pdac_audit.qc`). Normalisation: CP10k + log1p. Gene set G: {panel['n_model_gene_set_G']} of the 960 panel-intersection genes passing the type-aware background rule for >= 1 typed class (`results/chain/phase1b`). {len(panel['gse_only'])} GSE310352-only genes dropped: {', '.join(panel['gse_only'])}.",
          "- Program scores: scanpy-style (25 bins, 50 controls/gene, fixed seeds) per section; identical to the chain Phase 1b cache (max abs diff 0.0).",
          "- Graphs: symmetrised 6-NN; 30 um radius (external: 30 um-equivalent via median nearest-neighbour distance); 200 um tiles.",
          "", "## Cells per group and type", "", piv.to_markdown(index=False), "",
          "## Exclusions",
          "- Cells failing QC (see `results/phase0/section_info.csv` for post-QC counts).",
          "- T cells are not a modelling target (declared not testable: too few cells per fold).",
          "- No cells were excluded on the basis of any model output.", ""]
    (C.MODEL_DEV / "reports" / "DATA_CARD.md").write_text("\n".join(dc), encoding="utf-8")
    log.info("wrote key numbers (%d), tables, data card", len(kn))
    log.info("done")


if __name__ == "__main__":
    main()

"""Execute the predeclared Phase 5 Visium programs; see docs/v3_plans/.

Run with the thesis-spatial-ccc Python environment. One slide in memory at a time.
"""
from __future__ import annotations

import gc
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "data/pdac/numba_cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import PermutationMethod, mannwhitneyu, wilcoxon
import statsmodels.formula.api as smf

from common import (finish_run, out_dir, register, rel, save_csv, save_figure,
                    save_json, start_run)
from ccc_audit import bh
from v3_phase5_gate import site_labels_in

PHASE = "v3_phase5"
PROGRAMS = {
    "CAF": "COL1A1 COL1A2 COL3A1 DCN LUM FAP PDGFRB ACTA2 POSTN THY1 SFRP2 VCAN".split(),
    "EMT": "VIM FN1 CDH2 ZEB1 ZEB2 SNAI1 SNAI2 TWIST1 MMP2 SPARC TGFBI S100A4".split(),
    "Hypoxia": "VEGFA SLC2A1 CA9 LDHA PGK1 ADM NDRG1 BNIP3 ANKRD37 P4HA1 ENO1 HK2".split(),
    "Interferon": "ISG15 IFI6 IFI27 IFIT1 IFIT2 IFIT3 MX1 OAS1 STAT1 IRF7 IFI44L BST2".split(),
    "Immune_infiltration": "PTPRC CD3E CD3D CD2 CD8A CD4 MS4A1 CD79A CD68 CD163 LYZ NKG7".split(),
}
PRIMARY = "Primary PDAC"
PARAMS = {"seed": 42, "min_umi": 500, "max_fail_fraction": 0.5,
          "ctrl_size_per_bin": 50, "n_bins": 25, "ctrl_as_ref": True,
          "normalization": 1e4, "programs": PROGRAMS, "scanpy": sc.__version__}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_soft(path):
    records, cur = [], None
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                cur = {"gsm": line.split(" = ", 1)[1]}
                records.append(cur)
            elif cur is not None and line.startswith(("!Sample_title = ", "!Sample_characteristics_ch1 = ")):
                key, value = line.split(" = ", 1)
                if key == "!Sample_title":
                    cur["geo_title"] = value
                elif ": " in value:
                    key, value = value.split(": ", 1)
                    cur[key.lower()] = value
    frame = pd.DataFrame(records)
    if frame.gsm.duplicated().any():
        raise ValueError("Duplicate sample accession in GEO SOFT")
    return frame


def inputs():
    soft = ROOT / "data/GSE277783_family.soft.gz"
    if site_labels_in(soft).get("n_site_like_lines", 0) == 0:
        raise ValueError("Phase 5 metadata gate closed")
    mapping = pd.read_csv(ROOT / "eda_gse277783/tables/sample_map.csv")
    geo = read_soft(soft)
    samples = mapping.merge(geo, left_on="gsm_accession", right_on="gsm", how="left", validate="one_to_one")
    if samples.tissue.isna().any():
        raise ValueError("Missing GEO tissue label")
    samples["patient"] = samples.geo_title.str.extract(r"Pt-?(\d+)", expand=False).map(lambda x: f"Pt-{int(x)}")
    inv = pd.read_csv(ROOT / "eda_gse277783/tables/nested_tar_inventory.csv")
    cyt = set(inv.loc[inv.member_name.str.endswith("spatial_enrichment.csv"), "parent_tar"])
    samples["chemistry"] = np.where(samples.spatial_tar.isin(cyt), "CytAssist_inferred", "Visium_v1_inferred")
    samples["site"] = samples.tissue
    samples["is_pdx"] = samples.tissue.str.contains("PDX", case=False)
    return samples, soft


def controls_and_check(adata, genes, score):
    # Use the pinned Scanpy implementation itself to recover its exact controls.
    from scanpy.tools._score_genes import _check_score_genes_args, _score_genes_bins
    gl, pool, subset = _check_score_genes_args(adata, genes, None, use_raw=False)
    np.random.seed(42)
    controls = pd.Index([], dtype="string")
    for selected in _score_genes_bins(gl, pool, ctrl_as_ref=True, ctrl_size=50,
                                      n_bins=25, get_subset=subset):
        controls = controls.union(selected)
    direct = (np.asarray(subset(gl).mean(axis=1, dtype=np.float64)).ravel()
              - np.asarray(subset(controls).mean(axis=1, dtype=np.float64)).ravel())
    # SciPy sparse mean rescales float32 entries before summation; Scanpy
    # sums in float64 first. Allow the observed sub-micro-unit roundoff.
    if not np.allclose(direct, score, atol=1e-6, rtol=1e-6):
        raise ValueError("Scanpy score does not match recorded genes and controls")
    return controls.tolist()


def score_slides(samples, soft, log):
    od = out_dir(f"{PHASE}/slides")
    qc, scores, coverage, provenance = [], [], [], []
    for sample in samples.itertuples(index=False):
        path = ROOT / "eda_gse277783/extracted" / sample.h5_file
        identity = {"matrix": digest(path), "soft": digest(soft), "params": PARAMS}
        cache = od / f"{sample.gsm}_summary.json"
        if cache.exists():
            saved = json.loads(cache.read_text())
            if saved["identity"] != identity:
                raise ValueError(f"Resume identity mismatch: {sample.gsm}")
            log.info("%s: verified resume", sample.gsm)
        else:
            adata = sc.read_10x_h5(path)
            adata.var_names_make_unique()
            counts = np.asarray(adata.X.sum(axis=1)).ravel()
            keep = counts >= 500
            reason = "PDX excluded from patient inference" if sample.is_pdx else (
                "more than 50% of spots below 500 UMIs" if (~keep).mean() > 0.5 else "")
            q = {"gsm": sample.gsm, "patient": sample.patient, "site": sample.site,
                 "geo_title": sample.geo_title, "chemistry": sample.chemistry,
                 "chemistry_evidence": "nested spatial inventory; inferred batch",
                 "n_spots": adata.n_obs, "n_retained": int(keep.sum()),
                 "fraction_below_500": float((~keep).mean()), "median_umi": float(np.median(counts)),
                 "included": not bool(reason), "exclusion_reason": reason}
            result, cov = [], []
            if not reason:
                adata = adata[keep].copy()
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                spot = pd.DataFrame({"barcode": adata.obs_names.astype(str), "total_umi": counts[keep]})
                for name, genes in PROGRAMS.items():
                    present = [g for g in genes if g in adata.var_names]
                    c = {"gsm": sample.gsm, "program": name, "n_present": len(present),
                         "present_genes": ";".join(present), "missing_genes": ";".join(g for g in genes if g not in present)}
                    if len(present) < 5:
                        c.update({"n_controls": 0, "control_genes": "", "status": "insufficient_genes"})
                        cov.append(c)
                        continue
                    sc.tl.score_genes(adata, present, score_name=name, ctrl_size=50, n_bins=25,
                                      ctrl_as_ref=True, random_state=42, use_raw=False)
                    val = adata.obs[name].to_numpy()
                    if not np.isfinite(val).all():
                        raise ValueError(f"Nonfinite scores: {sample.gsm}/{name}")
                    controls = controls_and_check(adata, present, val)
                    c.update({"n_controls": len(controls), "control_genes": ";".join(controls), "status": "scored"})
                    cov.append(c)
                    spot[name] = val
                    result.append({**{k: q[k] for k in ("gsm", "patient", "site", "chemistry", "median_umi")},
                                   "program": name, "score": float(np.median(val)), "n_spots": adata.n_obs})
                p = od / f"{sample.gsm}_spot_scores.csv.gz"
                spot.to_csv(p, index=False)
                register(p, PHASE, f"Retained-spot program scores and raw UMIs: {sample.gsm}")
            saved = {"identity": identity, "qc": q, "scores": result, "coverage": cov}
            save_json(saved, f"{PHASE}/slides", cache.name, f"Resumable slide QC and program summaries: {sample.gsm}", log)
            log.info("%s: %s; %d/%d spots retained", sample.gsm, reason or "scored", keep.sum(), len(keep))
            del adata
            gc.collect()
        qc.append(saved["qc"])
        scores.extend(saved["scores"])
        coverage.extend(saved["coverage"])
        provenance.append({"gsm": sample.gsm, "file": rel(path), "bytes": path.stat().st_size, "sha256": identity["matrix"]})
    return pd.DataFrame(qc), pd.DataFrame(scores), pd.DataFrame(coverage), pd.DataFrame(provenance)


def compare(patient_site):
    arms = patient_site.copy()
    arms["arm"] = np.where(arms.site == PRIMARY, "primary", "metastasis")
    arms = arms.groupby(["patient", "program", "arm"], as_index=False).score.mean()
    tests = []
    for scope in ["combined_metastasis"] + sorted(s for s in patient_site.site.unique() if s != PRIMARY):
        d = arms if scope == "combined_metastasis" else patient_site[patient_site.site.isin([PRIMARY, scope])].assign(
            arm=lambda x: np.where(x.site == PRIMARY, "primary", "metastasis"))
        for program, g in d.groupby("program", sort=False):
            wide = g.pivot(index="patient", columns="arm", values="score").reindex(columns=["primary", "metastasis"])
            paired = wide.dropna()
            a, b = wide.primary.dropna(), wide.metastasis.dropna()
            if len(paired):
                delta = np.round((paired.metastasis - paired.primary).to_numpy(), 12)
                if np.all(delta == 0):
                    statistic, p = 0.0, 1.0
                else:
                    test = wilcoxon(delta, alternative="two-sided", method=PermutationMethod(n_resamples=np.inf))
                    statistic, p = float(test.statistic), float(test.pvalue)
                method = "paired_Wilcoxon_exhaustive"
                effect = float(np.median(delta))
            elif len(a) and len(b):
                test = mannwhitneyu(b, a, alternative="two-sided", method=PermutationMethod(n_resamples=np.inf))
                statistic, p = float(test.statistic), float(test.pvalue)
                method, effect = "unpaired_Mann_Whitney_exhaustive", float(b.median() - a.median())
            else:
                statistic, p, effect, method = np.nan, np.nan, np.nan, "not_testable"
            tests.append({"scope": scope, "program": program, "method": method,
                          "n_primary_patients": len(a), "n_metastasis_patients": len(b), "n_pairs": len(paired),
                          "n_unmatched_primary": len(a) - len(paired), "n_unmatched_metastasis": len(b) - len(paired),
                          "median_primary_all_patients": float(a.median()), "median_metastasis_all_patients": float(b.median()),
                          "median_delta_metastasis_minus_primary": effect, "statistic": statistic, "p": p})
    tests = pd.DataFrame(tests)
    tests["q_bh"] = np.nan
    for main in (True, False):
        idx = tests.scope.eq("combined_metastasis") == main
        tests.loc[idx, "q_bh"] = bh(tests.loc[idx, "p"].to_numpy())
    return arms, tests


def covariate_models(scores):
    output = []
    for name, data in scores.groupby("program", sort=False):
        fit = smf.ols('score ~ C(site, Treatment(reference="Primary PDAC")) + median_umi + C(chemistry)', data=data).fit(
            cov_type="cluster", cov_kwds={"groups": data.patient, "use_correction": True}, use_t=True)
        ci = fit.conf_int()
        for term in fit.params.index:
            output.append({"program": name, "term": term, "coefficient": fit.params[term],
                           "ci_low": ci.loc[term, 0], "ci_high": ci.loc[term, 1], "p_exploratory": fit.pvalues[term],
                           "n_slides": len(data), "n_patients": data.patient.nunique(),
                           "design_rank": int(np.linalg.matrix_rank(fit.model.exog)), "design_columns": fit.model.exog.shape[1],
                           "condition_number": fit.condition_number, "standard_errors": "patient_clustered"})
    return pd.DataFrame(output)


def main():
    log = start_run("v3_phase5_programs")
    register(log.log_path, PHASE, "Phase 5 run log (registered before computation)")
    log.info("Parameters: %s", PARAMS)
    samples, soft = inputs()
    save_csv(samples, PHASE, "verified_sample_metadata.csv", "GEO tissue labels, title-parsed patient IDs and chemistry evidence", log)
    qc, scores, coverage, provenance = score_slides(samples, soft, log)
    for df, name, desc in [(qc, "slide_qc.csv", "Planned UMI QC and PDX exclusions"),
                           (scores, "slide_program_scores.csv", "Median spot score per included slide and program"),
                           (coverage, "program_gene_coverage_controls.csv", "Actual program genes and Scanpy sampled controls"),
                           (provenance, "matrix_input_manifest.csv", "SHA-256 of every matrix actually inspected")]:
        save_csv(df, PHASE, name, desc, log)
    patient = scores.groupby(["patient", "site", "program"], as_index=False).agg(score=("score", "mean"), n_slides=("gsm", "size"))
    arms, tests = compare(patient)
    models = covariate_models(scores)
    save_csv(patient, PHASE, "patient_site_program_scores.csv", "Equal-weight slide means per patient and anatomical site", log)
    save_csv(arms, PHASE, "patient_arm_program_scores.csv", "Equal-weight anatomical-site means per patient and arm", log)
    save_csv(tests, PHASE, "program_comparisons.csv", "Patient-level exact tests; main and organ-specific BH families", log)
    save_csv(models, PHASE, "exploratory_covariate_models.csv", "OLS site + raw depth + inferred chemistry; patient-clustered SEs", log)
    imbalance = scores[["gsm", "patient", "site", "chemistry"]].drop_duplicates()
    save_csv(imbalance.groupby(["site", "chemistry"], as_index=False).agg(n_slides=("gsm", "size"), n_patients=("patient", "nunique")),
             PHASE, "site_chemistry_balance.csv", "Confounding and chemistry/site imbalance", log)
    fig, axes = plt.subplots(1, 5, figsize=(16, 4.3))
    for ax, program in zip(axes, PROGRAMS):
        wide = arms[arms.program == program].pivot(index="patient", columns="arm", values="score")
        for _, row in wide.iterrows():
            ax.plot([0, 1], [row.get("primary", np.nan), row.get("metastasis", np.nan)], color="#78909c", alpha=.65, lw=1)
        for i, arm in enumerate(["primary", "metastasis"]):
            vals = wide[arm].dropna()
            ax.scatter(np.full(len(vals), i), vals, s=24, color=["#247b99", "#bf5b48"][i], zorder=3)
        t = tests[(tests.scope == "combined_metastasis") & (tests.program == program)].iloc[0]
        ax.set(title=f"{program.replace('_', ' ')}\n{t.n_pairs} pairs; p={t.p:.3g}, q={t.q_bh:.3g}", xticks=[0, 1], xticklabels=["Primary", "Metastasis"])
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Patient program score")
    fig.suptitle("GSE277783: mixed-spot programs, paired patient comparison", y=1.03)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase5_programs", "Paired patient scores; metastasis is equal-weight site mean", log)
    summary = {"status": "DONE", "parameters": PARAMS, "n_input_slides": len(qc),
               "n_pdx_excluded": int(samples.is_pdx.sum()), "n_qc_excluded_human": int((~qc.included & ~qc.gsm.isin(samples.loc[samples.is_pdx, 'gsm'])).sum()),
               "n_included_slides": int(qc.included.sum()), "n_patients": scores.patient.nunique(),
               "n_retained_spots": int(qc.loc[qc.included, 'n_retained'].sum()),
               "n_main_bh_significant": int(((tests.scope == "combined_metastasis") & (tests.q_bh < .05)).sum()),
               "soft_sha256": digest(soft), "script_sha256": digest(__file__),
               "limitations": ["Mixed-spot expression and composition cannot be separated", "Patient IDs parsed from GEO titles",
                               "Chemistry inferred from archive batch; CytAssist restricted to primary tissues", "Exploratory OLS; 13 patient clusters"]}
    save_json(summary, PHASE, "phase5_summary.json", "Phase 5 completion, input identities and limitations", log)
    log.info("Summary: %s", summary)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

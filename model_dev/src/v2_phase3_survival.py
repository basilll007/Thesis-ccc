"""v2 Phase 3: the prepared translational survival analysis, executed per docs/translational_plans/phase3.md with
the amendments in model_dev/plans/v2_phase3_survival.md (three macrophage signatures, 9-test BH family, CPTAC CI,
power statement).

What: reuses src/translational/survival.py (scoring, fixed Cox PHReg models with BH, validated rank-time PH
diagnostic, Kaplan-Meier steps, Spearman correlations) and src/translational/phase3.py (verbatim v3 gene sets,
SHA-256-checked expression loading, cohort configs) without editing them. Writes to
results/translational/phase3/<UTC stamp>/ with a run-local MANIFEST_v2.csv.
Why: the survival models were planned and coded in translational v1 but never run.
"""
from __future__ import annotations

import csv
import hashlib
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "translational"))
import phase3 as T3  # noqa: E402  (imports translational common: caps threads at 4)
import survival as S  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

BASE6 = ["Hypoxia", "CAF", "EMT", "Interferon", "Immune_infiltration", "Immune_exclusion_proxy"]
NEW3 = ["Macrophage_C1Q_TAM", "Macrophage_SPP1_TAM", "Macrophage_C1Q_minus_SPP1"]
SIGS = BASE6 + NEW3
FAMILIES = ["unadjusted", "clinical", "clinical_purity"]


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        self.out = ROOT / "results" / "translational" / "phase3" / stamp
        self.out.mkdir(parents=True)
        self.log = logging.getLogger("v2_phase3" + stamp)
        self.log.setLevel(logging.INFO)
        for h in (logging.FileHandler(self.out / "run.log", encoding="utf-8"), logging.StreamHandler(sys.stdout)):
            h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            self.log.addHandler(h)
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        import importlib.metadata as md
        vers = {n: md.version(n) for n in ["numpy", "pandas", "scipy", "statsmodels", "matplotlib"]}
        self.log.info("v2 Phase 3 survival; commit=%s seed=42 max_threads=4 versions=%s", commit, vers)
        self.log.info("plan: docs/translational_plans/phase3.md + model_dev/plans/v2_phase3_survival.md")
        self.manifest = []

    def csv(self, df: pd.DataFrame, name: str, desc: str):
        p = self.out / name
        if p.exists():
            raise FileExistsError(p)
        df.to_csv(p, index=False)
        self._reg(p, desc)

    def json(self, obj, name: str, desc: str):
        p = self.out / name
        if p.exists():
            raise FileExistsError(p)
        p.write_text(json.dumps(obj, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)), encoding="utf-8")
        self._reg(p, desc)

    def _reg(self, p: Path, desc: str):
        self.manifest.append({"path": p.relative_to(ROOT).as_posix(), "description": desc, "bytes": p.stat().st_size, "sha256": sha(p)})

    def close(self):
        for h in self.log.handlers:
            h.flush()
        self._reg(self.out / "run.log", "run log")
        with (self.out / "MANIFEST_v2.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["path", "description", "bytes", "sha256"])
            w.writeheader(); w.writerows(self.manifest)


def gene_sets() -> dict:
    sets = T3.gene_sets()   # verbatim v3 programs + fixed T-cell proxy
    prog = json.loads((ROOT / "results" / "chain" / "phase1b" / "programs_used.json").read_text())["programs"]
    sets["Macrophage_C1Q_TAM"] = list(prog["myeloid:C1Q_TAM"])
    sets["Macrophage_SPP1_TAM"] = list(prog["myeloid:SPP1_TAM"])
    return sets


def add_difference(expression: pd.DataFrame, sets: dict, scores: pd.DataFrame, cov: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Macrophage_C1Q_minus_SPP1 from the two component mean-z scores (same construction as the plan's proxy)."""
    raw = {}
    for name in ("Macrophage_C1Q_TAM", "Macrophage_SPP1_TAM"):
        row = cov[cov.signature == name].iloc[0]
        if row.status != "SCORABLE":
            return scores, cov
        genes = row.genes_used.split(";")
        v = expression[genes]
        raw[name] = ((v - v.mean()) / v.std(ddof=1)).mean(axis=1)
    d = raw["Macrophage_C1Q_TAM"] - raw["Macrophage_SPP1_TAM"]
    scores = scores.copy()
    scores["Macrophage_C1Q_minus_SPP1"] = (d - d.mean()) / d.std(ddof=1)
    extra = {"signature": "Macrophage_C1Q_minus_SPP1", "n_declared": 11, "n_present_finite_variable": np.nan,
             "coverage_fraction": np.nan, "genes_used": "C1Q mean gene z minus SPP1 mean gene z",
             "missing_or_nonfinite_or_constant": "See component coverage", "status": "SCORABLE"}
    return scores, pd.concat([cov, pd.DataFrame([extra])], ignore_index=True)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def power_table(clinical: pd.DataFrame, scores: pd.DataFrame, models: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """Minimum detectable HR per SD at 80 % power (Schoenfeld), with variance inflation from the adjustment."""
    d = clinical.set_index("PATIENT_ID").join(scores)
    rows = []
    for _, m in models.iterrows():
        if m.signature not in d or not np.isfinite(m.get("events", np.nan)):
            continue
        req = [m.signature, "time_months", "event"] + (["age", "sex", "stage"] if m.model != "unadjusted" else []) + (["purity"] if m.model == "clinical_purity" else [])
        g = d.dropna(subset=req)
        r2 = 0.0
        if m.model != "unadjusted":
            X = pd.DataFrame({"age": g.age / 10, "male": g.sex.eq("Male").astype(float), "s2": g.stage.eq("II").astype(float),
                              "s34": g.stage.eq("III_IV").astype(float)})
            if m.model == "clinical_purity":
                X["purity"] = g.purity
            X = np.column_stack([np.ones(len(X)), X.to_numpy(float)])
            y = g[m.signature].to_numpy(float)
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            r2 = 1 - ((y - X @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        sd = g[m.signature].std(ddof=1)
        D = int(m.events)
        for alpha, label in ((0.05, "0.05"), (0.05 / 9, "0.05/9")):
            z = norm.ppf(1 - alpha / 2) + norm.ppf(0.8)
            rows.append({"cohort": cohort, "signature": m.signature, "model": m.model, "events": D, "alpha": label,
                         "signature_sd_in_model_subset": sd, "r2_on_adjusters": r2,
                         "min_detectable_hr_per_sd": float(np.exp(z / np.sqrt(D * (1 - r2) * sd ** 2)))})
    return pd.DataFrame(rows)


def main() -> None:
    run = Run()
    log = run.log
    try:
        validation, receipt = S.validate_ph()
        run.csv(validation, "PH_synthetic_validation.csv", "PH diagnostic validation: 100 proportional + 100 time-varying synthetic examples, gradient check")
        run.json(receipt, "PH_validation_receipt.json", "PH diagnostic validation receipt")
        log.info("PH validation: %s", receipt)
        if receipt["status"] != "PASS":
            raise RuntimeError("PH diagnostic validation failed; no outcome model fitted")
        sets = gene_sets()
        run.json(sets, "declared_gene_sets.json", "verbatim v3 programs, T-cell proxy, and the two macrophage programs (programs_used.json)")
        models, coefs, ph, resid, corr, km, risk, cuts, cov_all, power = [], [], [], [], [], [], [], [], [], []
        scored = {}
        for conf in T3.cohort_configs():
            name = conf["cohort"]
            clinical = pd.read_csv(ROOT / conf["clinical_table"])
            expression = T3.load_expression(conf, clinical)
            scores, cov = S.score_signatures(expression, sets)
            scores, cov = add_difference(expression, sets, scores, cov)
            cov["cohort"] = name
            cov_all.append(cov)
            run.csv(scores.rename_axis("PATIENT_ID").reset_index(), f"{name}_patient_scores.csv", "per-patient signature scores (SD 1)")
            m, c, p, r, _ = S.fit_models(clinical, scores, SIGS, name)       # BH over all 9 signatures (A2)
            m["q_bh_original_six"] = np.nan
            for fam, ix in m[m.signature.isin(BASE6)].groupby("model").groups.items():
                pv = m.loc[ix, "p"]
                q = np.full(len(pv), np.nan)
                ok = pv.notna().to_numpy()
                if ok.any():
                    q[ok] = bh(pv[ok].to_numpy())
                m.loc[ix, "q_bh_original_six"] = q
            models.append(m); coefs.append(c); ph.append(p); resid.append(r)
            co = S.correlations(clinical, scores, SIGS, name); corr.append(co)
            cu, ar, ct = S.km_data(clinical, scores, SIGS, name); km.append(cu); risk.append(ar); cuts.append(ct)
            power.append(power_table(clinical, scores, m, name))
            if conf["strict"]:
                scored[conf["dataset"]] = {"cohort": name, "models": m, "corr": co, "curves": cu}
            log.info("%s: n=%d; model status %s", name, len(clinical), m.status.value_counts().to_dict())

        def cat(fr):
            fr = [f for f in fr if len(f)]
            return pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()
        M = cat(models)
        for df, nm, desc in [(M, "cox_models.csv", "9 signatures x 3 fixed Cox families per cohort; q_bh = BH over 9 (primary, A2); q_bh_original_six = plan's six-signature BH"),
                             (cat(coefs), "cox_all_coefficients.csv", "all fitted coefficients"),
                             (cat(ph), "PH_diagnostics.csv", "rank-time scaled Schoenfeld diagnostic; BH across all PH tests per cohort/family"),
                             (cat(resid), "scaled_schoenfeld_residuals.csv", "event-level scaled Schoenfeld residuals"),
                             (cat(corr), "signature_correlations.csv", "Spearman: signatures vs purity and each other; purity-proxy flag |rho|>=0.7"),
                             (cat(cov_all), "gene_coverage.csv", "declared/present/finite/variable genes per signature and cohort"),
                             (cat(km), "KM_curves.csv", "median-split descriptive KM steps"), (cat(risk), "KM_at_risk.csv", "at-risk counts"),
                             (cat(cuts), "KM_cutoffs.csv", "per-cohort median cut-points (not optimised)"),
                             (cat(power), "power_min_detectable_hr.csv", "Schoenfeld minimum detectable HR per SD at 80% power (A4)")]:
            run.csv(df, nm, desc)
        # validation in CPTAC (A3)
        rows = []
        tc, cp = M[M.cohort == "TCGA_confirmed_PDAC"], M[M.cohort == "CPTAC_confirmed_PDAC"]
        for sgn in SIGS:
            for fam in FAMILIES:
                a = tc[(tc.signature == sgn) & (tc.model == fam)]
                b = cp[(cp.signature == sgn) & (cp.model == fam)]
                a = a.iloc[0] if len(a) else None
                b = b.iloc[0] if len(b) else None
                ok = a is not None and b is not None and a.status == "DONE" and b.status == "DONE"
                same = ("yes" if np.sign(a.beta) == np.sign(b.beta) else "no") if ok else "UNTESTABLE"
                rows.append({"signature": sgn, "model": fam, "TCGA_status": None if a is None else a.status,
                             "TCGA_HR": a.get("hr", np.nan) if a is not None else np.nan, "TCGA_ci_low": a.get("ci_low", np.nan) if a is not None else np.nan,
                             "TCGA_ci_high": a.get("ci_high", np.nan) if a is not None else np.nan, "TCGA_q_bh9": a.get("q_bh", np.nan) if a is not None else np.nan,
                             "CPTAC_status": None if b is None else b.status, "CPTAC_HR": b.get("hr", np.nan) if b is not None else np.nan,
                             "CPTAC_ci_low": b.get("ci_low", np.nan) if b is not None else np.nan, "CPTAC_ci_high": b.get("ci_high", np.nan) if b is not None else np.nan,
                             "CPTAC_p": b.get("p", np.nan) if b is not None else np.nan, "same_direction": same,
                             "CPTAC_ci_excludes_1": bool(ok and (b.ci_low > 1 or b.ci_high < 1)),
                             "holds_in_CPTAC_A3": bool(ok and same == "yes" and (b.ci_low > 1 or b.ci_high < 1))})
        V = pd.DataFrame(rows)
        run.csv(V, "external_validation_cptac.csv", "unchanged models in CPTAC: HR, CI, same direction, holds (A3)")
        primary = tc[tc.model == "clinical_purity"]
        passed = primary[(primary.status == "DONE") & (primary.q_bh < 0.05)].signature.tolist()
        pht = cat(ph)
        ph_sig = pht[(pht.cohort == "TCGA_confirmed_PDAC") & (pht.model == "clinical_purity") & (pht.term == "signature")]
        summary = {"status": "DONE", "n_signatures": len(SIGS), "bh_family_size": 9, "n_primary_models_strict_tcga": int(len(tc)),
                   "primary_adjusted_q_lt_0.05": passed,
                   "primary_adjusted": primary[["signature", "status", "reason", "n", "events", "hr", "ci_low", "ci_high", "p", "q_bh", "q_bh_original_six"]].to_dict(orient="records"),
                   "ph_signature_tests_primary": ph_sig[["signature", "p", "q_bh_PH"]].to_dict(orient="records") if len(ph_sig) else [],
                   "cptac_holds": V[V.holds_in_CPTAC_A3].to_dict(orient="records"),
                   "PH_validation": receipt}
        run.json(summary, "summary.json", "Phase 3 outcome summary")
        figure(run, scored, M)
        log.info("primary adjusted q<0.05: %s", passed)
    except Exception as e:
        log.exception("Phase 3 failure")
        run.json({"status": "BLOCKED", "reason": str(e)}, "blocked.json", "Phase 3 blocker")
        raise
    finally:
        run.close()


def figure(run: Run, scored: dict, M: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(24, 13))
    gs = fig.add_gridspec(3, 9, height_ratios=[1, 1, 1.3])
    for row, ds in enumerate(["TCGA", "CPTAC"]):
        info = scored.get(ds)
        if info is None:
            continue
        for col, sgn in enumerate(SIGS):
            ax = fig.add_subplot(gs[row, col])
            g = info["curves"][info["curves"].signature == sgn]
            for lab, color in (("low", "#2a6f97"), ("high", "#c0504d")):
                v = g[g.group == lab]
                if len(v):
                    ax.step(v.time, v.survival, where="post", color=color, lw=1.2, label=f"{lab} (n={int(v.n_group.iloc[0])})")
            ax.set_title(f"{ds}: {sgn.replace('_', ' ')}", fontsize=8)
            ax.set_ylim(0, 1.05); ax.set_xlabel("Months", fontsize=7); ax.tick_params(labelsize=6)
            ax.legend(fontsize=5, frameon=False); ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylabel("OS probability (descriptive)", fontsize=7)
    ax = fig.add_subplot(gs[2, :5])
    y = 0; ticks = []
    for sgn in SIGS:
        for ds, color, off in (("TCGA", "#2a6f97", 0.0), ("CPTAC", "#c0504d", 0.3)):
            s = M[(M.cohort == f"{ds}_confirmed_PDAC") & (M.signature == sgn) & (M.model == "clinical_purity")]
            if len(s) and s.status.iloc[0] == "DONE":
                r = s.iloc[0]
                ax.errorbar(r.hr, y + off, xerr=[[r.hr - r.ci_low], [r.ci_high - r.hr]], fmt="o", color=color, ms=4, lw=1)
            ticks.append((y + off, f"{ds} {sgn}"))
        y += 1
    ax.set_yticks([t[0] for t in ticks], [t[1] for t in ticks], fontsize=6)
    ax.axvline(1, color="black", ls="--", lw=0.6); ax.set_xscale("log")
    ax.set_xlabel("HR per SD (age, sex, stage, purity adjusted)"); ax.set_title("Strict ductal cohorts; fixed models", fontsize=9)
    ax = fig.add_subplot(gs[2, 6:])
    co = scored["TCGA"]["corr"]; names = SIGS + ["purity"]
    mat = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for r in co.itertuples():
        mat.loc[r.signature_a, r.signature_b] = mat.loc[r.signature_b, r.signature_a] = r.spearman_rho
    im = ax.imshow(mat.to_numpy(float), vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(names)), names, rotation=90, fontsize=6); ax.set_yticks(range(len(names)), names, fontsize=6)
    ax.set_title("TCGA Spearman correlations (signatures, purity)", fontsize=9); fig.colorbar(im, ax=ax, fraction=0.04)
    fig.suptitle("Spatial-program and macrophage signatures in bulk PDAC cohorts: prognostic expression associations (not causal)", y=1.0)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = run.out / f"signature_survival_v2.{ext}"
        fig.savefig(p, dpi=150 if ext == "png" else None, bbox_inches="tight")
        run._reg(p, "Phase 3 composite figure")
    plt.close(fig)


if __name__ == "__main__":
    main()

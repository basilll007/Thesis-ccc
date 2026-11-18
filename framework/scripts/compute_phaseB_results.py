import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import norm

from run_phaseB_survival import df_gse62, df_gse28, PROBE_MAP

def score_cohort(df):
    genes = list(PROBE_MAP.keys())
    # Verify gene coverage
    covered = [g for g in genes if g in df.columns and not df[g].isna().all()]
    coverage_pct = len(covered) / len(genes) * 100.0
    print(f"Gene coverage: {len(covered)}/{len(genes)} ({coverage_pct:.1f}%)")
    
    # z-score each gene
    Z = pd.DataFrame(index=df.index)
    for g in covered:
        vals = df[g].to_numpy(dtype=float)
        mean, std = np.nanmean(vals), np.nanstd(vals, ddof=1)
        Z[g] = (vals - mean) / (std if std > 1e-12 else 1.0)
        
    raw_score = Z.mean(axis=1)
    # Standardise signature
    std_score = (raw_score - raw_score.mean()) / raw_score.std(ddof=1)
    df["hypoxia_score"] = std_score
    return df, coverage_pct

print("=== Scoring GSE62452 ===")
df_gse62, cov62 = score_cohort(df_gse62)

print("\n=== Scoring GSE28735 ===")
df_gse28, cov28 = score_cohort(df_gse28)

# 1. Cox model GSE62452 (unadjusted)
cph62 = CoxPHFitter()
cph62.fit(df_gse62[["os_months", "event", "hypoxia_score"]], duration_col="os_months", event_col="event")
s62 = cph62.summary.loc["hypoxia_score"]
b62, se62, p62 = float(s62["coef"]), float(s62["se(coef)"]), float(s62["p"])
hr62, ci_low62, ci_high62 = float(np.exp(b62)), float(np.exp(b62 - 1.96 * se62)), float(np.exp(b62 + 1.96 * se62))
print(f"GSE62452 (unadjusted): HR = {hr62:.3f} [{ci_low62:.3f}, {ci_high62:.3f}], beta = {b62:.4f}, SE = {se62:.4f}, p = {p62:.4f}")

# GSE62452 stage-adjusted
df62_adj = df_gse62[["os_months", "event", "hypoxia_score", "stage"]].copy()
# Simplify stage to early (I/IIA) vs advanced (IIB/III/IV)
df62_adj["is_advanced_stage"] = df62_adj["stage"].astype(str).str.contains("IIB|III|IV").astype(int)
cph62_adj = CoxPHFitter()
cph62_adj.fit(df62_adj[["os_months", "event", "hypoxia_score", "is_advanced_stage"]], duration_col="os_months", event_col="event")
s62_adj = cph62_adj.summary.loc["hypoxia_score"]
b62_adj, se62_adj, p62_adj = float(s62_adj["coef"]), float(s62_adj["se(coef)"]), float(s62_adj["p"])
hr62_adj = float(np.exp(b62_adj))
ci_low62_adj = float(np.exp(b62_adj - 1.96 * se62_adj))
ci_high62_adj = float(np.exp(b62_adj + 1.96 * se62_adj))
print(f"GSE62452 (stage-adjusted): HR = {hr62_adj:.3f} [{ci_low62_adj:.3f}, {ci_high62_adj:.3f}], beta = {b62_adj:.4f}, SE = {se62_adj:.4f}, p = {p62_adj:.4f}")

# 2. Cox model GSE28735
cph28 = CoxPHFitter()
cph28.fit(df_gse28[["os_months", "event", "hypoxia_score"]], duration_col="os_months", event_col="event")
s28 = cph28.summary.loc["hypoxia_score"]
b28, se28, p28 = float(s28["coef"]), float(s28["se(coef)"]), float(s28["p"])
hr28, ci_low28, ci_high28 = float(np.exp(b28)), float(np.exp(b28 - 1.96 * se28)), float(np.exp(b28 + 1.96 * se28))
print(f"GSE28735 (unadjusted): HR = {hr28:.3f} [{ci_low28:.3f}, {ci_high28:.3f}], beta = {b28:.4f}, SE = {se28:.4f}, p = {p28:.4f}")

def random_effects_meta_analysis(betas, ses):
    betas = np.array(betas, dtype=float)
    ses = np.array(ses, dtype=float)
    k = len(betas)
    w_fe = 1.0 / (ses ** 2)
    beta_fe = np.sum(w_fe * betas) / np.sum(w_fe)
    Q = np.sum(w_fe * (betas - beta_fe) ** 2)
    df = k - 1
    C = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_re = 1.0 / (ses ** 2 + tau2)
    beta_re = np.sum(w_re * betas) / np.sum(w_re)
    se_re = np.sqrt(1.0 / np.sum(w_re))
    z = beta_re / se_re
    p_one_sided = 1.0 - norm.cdf(z)
    p_two_sided = 2.0 * (1.0 - norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    return {
        "beta_pooled": float(beta_re),
        "se_pooled": float(se_re),
        "hr_pooled": float(np.exp(beta_re)),
        "ci_low": float(np.exp(beta_re - 1.96 * se_re)),
        "ci_high": float(np.exp(beta_re + 1.96 * se_re)),
        "z": float(z),
        "p_one_sided": float(p_one_sided),
        "p_two_sided": float(p_two_sided),
        "Q": float(Q),
        "I2_pct": float(I2),
        "tau2": float(tau2),
    }

# Primary Meta-analysis: NEW cohorts only (using stage-adjusted for GSE62452)
meta_new = random_effects_meta_analysis([b62_adj, b28], [se62_adj, se28])
print("\n=== PRIMARY META-ANALYSIS (NEW INDEPENDENT COHORTS ONLY) ===")
print(f"Pooled HR = {meta_new['hr_pooled']:.3f} [{meta_new['ci_low']:.3f}, {meta_new['ci_high']:.3f}]")
print(f"Pooled beta = {meta_new['beta_pooled']:.4f}, SE = {meta_new['se_pooled']:.4f}")
print(f"One-sided p (beta > 0) = {meta_new['p_one_sided']:.4f}")
print(f"Two-sided p = {meta_new['p_two_sided']:.4f}")
print(f"Heterogeneity: Q = {meta_new['Q']:.3f}, I2 = {meta_new['I2_pct']:.1f}%, tau2 = {meta_new['tau2']:.4f}")

# Secondary Meta-analysis: including TCGA and CPTAC
b_tcga = float(np.log(1.32))
se_tcga = float((np.log(1.62) - np.log(1.07)) / (2 * 1.96))
b_cptac = float(np.log(1.43))
se_cptac = float((np.log(1.86) - np.log(1.10)) / (2 * 1.96))

meta_all = random_effects_meta_analysis([b62_adj, b28, b_tcga, b_cptac], [se62_adj, se28, se_tcga, se_cptac])
print("\n=== SECONDARY META-ANALYSIS (INCLUDING DISCOVERY TCGA & CPTAC) ===")
print(f"Pooled HR = {meta_all['hr_pooled']:.3f} [{meta_all['ci_low']:.3f}, {meta_all['ci_high']:.3f}]")
print(f"Pooled beta = {meta_all['beta_pooled']:.4f}, SE = {meta_all['se_pooled']:.4f}")
print(f"One-sided p (beta > 0) = {meta_all['p_one_sided']:.4f}")
print(f"Two-sided p = {meta_all['p_two_sided']:.4f}")
print(f"Heterogeneity: Q = {meta_all['Q']:.3f}, I2 = {meta_all['I2_pct']:.1f}%, tau2 = {meta_all['tau2']:.4f}")

# Power & MDHR
pooled_events = 49 + 29
mdhr = float(np.exp((1.645 + 0.842) / np.sqrt(pooled_events)))
print(f"\nTotal pooled independent events: {pooled_events}")
print(f"Minimum detectable HR at 80% power (alpha=0.05, one-sided): {mdhr:.3f}")

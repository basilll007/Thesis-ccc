import gzip
import io
import urllib.request
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import norm

from compute_phaseB_results import random_effects_meta_analysis, df_gse62, df_gse28

# Load all declared gene sets from translational phase 3
with open("results/translational/phase3/20261117T024638285521Z/declared_gene_sets.json", "r") as f:
    ALL_SIGNATURES = json.load(f)

# Load GPL6244 annotation map
print("Loading GPL6244 gene to probe map...")
url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    content = resp.read()
gz = gzip.GzipFile(fileobj=io.BytesIO(content))

table_started = False
all_target_genes = {g for sig in ALL_SIGNATURES.values() for g in sig}
gene_to_probe = {}

for line in gz:
    l = line.decode("utf-8", errors="ignore")
    if l.startswith("!platform_table_begin"):
        table_started = True
        header = next(gz).decode("utf-8", errors="ignore").strip().split("\t")
        id_col = header.index("ID")
        sym_col = header.index("Gene symbol")
        continue
    if table_started:
        if l.startswith("!platform_table_end"):
            break
        parts = l.strip().split("\t")
        if len(parts) > sym_col:
            sym = parts[sym_col].strip()
            for s in sym.split("///"):
                s = s.strip()
                if s in all_target_genes and s not in gene_to_probe:
                    gene_to_probe[s] = parts[id_col]

print(f"Mapped {len(gene_to_probe)} of {len(all_target_genes)} genes in GPL6244.")

def extract_cohort_matrix(acc, target_probes):
    prefix = acc[:5] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        content = resp.read()
    gz = gzip.GzipFile(fileobj=io.BytesIO(content))
    table_started = False
    sample_ids = []
    matrix = {}
    for line in gz:
        l = line.decode("utf-8", errors="ignore")
        if not table_started:
            if l.startswith("!Sample_geo_accession"):
                sample_ids = [p.strip('"\n ') for p in l.split("\t")[1:]]
            if l.startswith("!series_matrix_table_begin"):
                table_started = True
            continue
        if l.startswith("!series_matrix_table_end"):
            break
        parts = l.strip().split("\t")
        pid = parts[0].strip('"\n ')
        if pid in target_probes:
            gene = target_probes[pid]
            matrix[gene] = [float(x) if x != "null" and x != "" else np.nan for x in parts[1:]]
    return pd.DataFrame(matrix, index=sample_ids)

probe_to_gene = {p: g for g, p in gene_to_probe.items()}
print("Extracting expression for GSE62452...")
exp62 = extract_cohort_matrix("GSE62452", probe_to_gene)
print("Extracting expression for GSE28735...")
exp28 = extract_cohort_matrix("GSE28735", probe_to_gene)

# Match with clinical tables
exp62_tumour = exp62.reindex(df_gse62["sample_id"]).dropna(axis=1, how="all")
exp28_tumour = exp28.reindex(df_gse28["sample_id"]).dropna(axis=1, how="all")

exploratory_results = []

def score_sig(exp_df, gene_list):
    present = [g for g in gene_list if g in exp_df.columns]
    coverage = len(present) / len(gene_list)
    if coverage < 0.5:
        return None, coverage
    # z-score
    Z = pd.DataFrame(index=exp_df.index)
    for g in present:
        vals = exp_df[g].to_numpy(dtype=float)
        m, s = np.nanmean(vals), np.nanstd(vals, ddof=1)
        Z[g] = (vals - m) / (s if s > 1e-12 else 1.0)
    raw = Z.mean(axis=1)
    std = (raw - raw.mean()) / raw.std(ddof=1)
    return std, coverage

for sig_name, genes in ALL_SIGNATURES.items():
    if sig_name == "Hypoxia":
        continue # Primary
        
    score62, cov62 = score_sig(exp62_tumour, genes)
    score28, cov28 = score_sig(exp28_tumour, genes)
    
    if score62 is None or score28 is None:
        continue
        
    # Fit Cox GSE62452
    d62 = df_gse62[["os_months", "event"]].copy()
    d62["score"] = score62.to_numpy()
    c62 = CoxPHFitter().fit(d62, duration_col="os_months", event_col="event")
    b62 = float(c62.summary.loc["score", "coef"])
    se62 = float(c62.summary.loc["score", "se(coef)"])
    
    # Fit Cox GSE28735
    d28 = df_gse28[["os_months", "event"]].copy()
    d28["score"] = score28.to_numpy()
    c28 = CoxPHFitter().fit(d28, duration_col="os_months", event_col="event")
    b28 = float(c28.summary.loc["score", "coef"])
    se28 = float(c28.summary.loc["score", "se(coef)"])
    
    meta = random_effects_meta_analysis([b62, b28], [se62, se28])
    exploratory_results.append({
        "signature": sig_name,
        "n_genes_declared": len(genes),
        "coverage_gse62": cov62,
        "coverage_gse28": cov28,
        "gse62_hr": float(np.exp(b62)),
        "gse62_p": float(c62.summary.loc["score", "p"]),
        "gse28_hr": float(np.exp(b28)),
        "gse28_p": float(c28.summary.loc["score", "p"]),
        "pooled_hr": meta["hr_pooled"],
        "ci_low": meta["ci_low"],
        "ci_high": meta["ci_high"],
        "pooled_beta": meta["beta_pooled"],
        "pooled_se": meta["se_pooled"],
        "p_two_sided": meta["p_two_sided"],
    })

# Add C1Q minus SPP1
score62_c1q, _ = score_sig(exp62_tumour, ALL_SIGNATURES["Macrophage_C1Q_TAM"])
score62_spp1, _ = score_sig(exp62_tumour, ALL_SIGNATURES["Macrophage_SPP1_TAM"])
diff62 = score62_c1q - score62_spp1
diff62_std = (diff62 - diff62.mean()) / diff62.std(ddof=1)

score28_c1q, _ = score_sig(exp28_tumour, ALL_SIGNATURES["Macrophage_C1Q_TAM"])
score28_spp1, _ = score_sig(exp28_tumour, ALL_SIGNATURES["Macrophage_SPP1_TAM"])
diff28 = score28_c1q - score28_spp1
diff28_std = (diff28 - diff28.mean()) / diff28.std(ddof=1)

d62 = df_gse62[["os_months", "event"]].copy()
d62["score"] = diff62_std.to_numpy()
c62 = CoxPHFitter().fit(d62, duration_col="os_months", event_col="event")
b62 = float(c62.summary.loc["score", "coef"])
se62 = float(c62.summary.loc["score", "se(coef)"])

d28 = df_gse28[["os_months", "event"]].copy()
d28["score"] = diff28_std.to_numpy()
c28 = CoxPHFitter().fit(d28, duration_col="os_months", event_col="event")
b28 = float(c28.summary.loc["score", "coef"])
se28 = float(c28.summary.loc["score", "se(coef)"])

meta_diff = random_effects_meta_analysis([b62, b28], [se62, se28])
exploratory_results.append({
    "signature": "C1Q_minus_SPP1",
    "n_genes_declared": len(ALL_SIGNATURES["Macrophage_C1Q_TAM"]) + len(ALL_SIGNATURES["Macrophage_SPP1_TAM"]),
    "coverage_gse62": 1.0,
    "coverage_gse28": 1.0,
    "gse62_hr": float(np.exp(b62)),
    "gse62_p": float(c62.summary.loc["score", "p"]),
    "gse28_hr": float(np.exp(b28)),
    "gse28_p": float(c28.summary.loc["score", "p"]),
    "pooled_hr": meta_diff["hr_pooled"],
    "ci_low": meta_diff["ci_low"],
    "ci_high": meta_diff["ci_high"],
    "pooled_beta": meta_diff["beta_pooled"],
    "pooled_se": meta_diff["se_pooled"],
    "p_two_sided": meta_diff["p_two_sided"],
})

df_exp = pd.DataFrame(exploratory_results)
# BH correction across the exploratory signatures
p_vals = df_exp["p_two_sided"].to_numpy()
m = len(p_vals)
order = np.argsort(p_vals)
ranked = p_vals[order] * m / np.arange(1, m + 1)
ranked = np.minimum.accumulate(ranked[::-1])[::-1]
q_vals = np.empty(m)
q_vals[order] = np.minimum(ranked, 1.0)
df_exp["bh_q"] = q_vals

print("\n=== EXPLORATORY SIGNATURES META-ANALYSIS ===")
print(df_exp[["signature", "pooled_hr", "ci_low", "ci_high", "p_two_sided", "bh_q"]].to_string(index=False))
df_exp.to_csv("results/translational/phaseB/exploratory_signatures.csv", index=False)

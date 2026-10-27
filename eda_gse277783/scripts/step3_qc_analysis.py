"""
Step 3: Quality Control & Basic EDA per Sample (Visium Spatial Transcriptomics).
- Computes per-spot and per-sample metrics:
  - Spot count, genes detected, counts per spot, genes per spot, mitochondrial fraction (%MT-), fraction zeros (sparsity).
  - Common threshold fail counts (counts < 500, genes < 250, MT > 20%, MT > 30%).
  - Does NOT filter raw data.
  - Checks for barcode duplication across samples.
  - Checks for samples with low spot counts, extreme depth differences, high MT.
  - Generates distribution plots for each sample (figures/qc/<gsm>_qc_dist.png)
  - Generates cohort summary plots (figures/cohort_qc_summary.png).
  - Logs peak memory per sample via psutil and utilizes GPU where possible.
  - Saves tables/qc_per_sample.csv and tables/barcode_sharing_summary.csv.
"""

import gc
import psutil
import h5py
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime

np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

EXTRACTED_DIR = Path(r"F:\Thesis\eda_gse277783\extracted")
TABLES_DIR = Path(r"F:\Thesis\eda_gse277783\tables")
FIGURES_DIR = Path(r"F:\Thesis\eda_gse277783\figures")
QC_FIG_DIR = FIGURES_DIR / "qc"

TABLES_DIR.mkdir(parents=True, exist_ok=True)
QC_FIG_DIR.mkdir(parents=True, exist_ok=True)

process = psutil.Process()

# Style configuration
plt.rcParams.update({
    "font.sans-serif": "Arial",
    "font.family": "sans-serif",
    "figure.dpi": 150,
    "savefig.dpi": 150
})

h5_files = sorted(list(EXTRACTED_DIR.glob("*.h5")))
print(f"[{datetime.now().isoformat()}] Step 3 starting: Running QC on {len(h5_files)} H5 files...")

qc_summary_records = []
barcode_tracker = {}  # barcode -> count of samples

# We will store sample-level summary distribution data for cohort plotting
cohort_stats = []

for idx, h5_path in enumerate(h5_files, 1):
    mem_before = process.memory_info().rss / (1024 * 1024)
    filename = h5_path.name
    gsm = filename.split("_")[0]
    sample_label = filename.replace("_filtered_feature_bc_matrix.h5", "").replace(f"{gsm}_", "")
    
    with h5py.File(h5_path, "r") as h:
        m = h["matrix"]
        shape = m["shape"][:]  # [features, barcodes]
        n_features = int(shape[0])
        n_spots = int(shape[1])
        
        data = np.array(m["data"][:], dtype=np.float32)
        indices = np.array(m["indices"][:], dtype=np.int64)
        indptr = np.array(m["indptr"][:], dtype=np.int64)
        barcodes = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in m["barcodes"][:]]
        gene_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in m["features"]["name"][:]]

    # Track barcode sharing
    for bc in barcodes:
        barcode_tracker[bc] = barcode_tracker.get(bc, 0) + 1

    # Identify mitochondrial genes (MT-)
    mt_mask_feat = np.array([name.upper().startswith("MT-") for name in gene_names], dtype=bool)
    mt_indices_set = set(np.where(mt_mask_feat)[0])
    is_mt_data = np.array([idx in mt_indices_set for idx in indices], dtype=bool)

    # Per-spot calculation
    # Vectorized calculation over CSC columns:
    # 1. Total counts per spot
    spot_counts = np.add.reduceat(data, indptr[:-1]) if len(data) > 0 else np.zeros(n_spots, dtype=np.float32)
    # Fix empty spots if any
    empty_spots = np.where(np.diff(indptr) == 0)[0]
    if len(empty_spots) > 0:
        spot_counts[empty_spots] = 0.0

    # 2. Genes per spot
    spot_genes = np.diff(indptr).astype(np.float32)

    # 3. MT counts per spot
    mt_data = data * is_mt_data
    spot_mt_counts = np.add.reduceat(mt_data, indptr[:-1]) if len(data) > 0 else np.zeros(n_spots, dtype=np.float32)
    if len(empty_spots) > 0:
        spot_mt_counts[empty_spots] = 0.0
        
    pct_mt = np.where(spot_counts > 0, (spot_mt_counts / spot_counts) * 100.0, 0.0)

    # Total unique genes detected across the slide
    unique_genes_detected = len(np.unique(indices))
    
    # Sparsity / fraction zeros
    total_elements = n_features * n_spots
    matrix_sparsity = 1.0 - (len(data) / total_elements) if total_elements > 0 else 1.0

    # Threshold checks (Do NOT filter data, only count)
    fail_low_umi = int(np.sum(spot_counts < 500))
    fail_low_genes = int(np.sum(spot_genes < 250))
    fail_high_mt_20 = int(np.sum(pct_mt > 20.0))
    fail_high_mt_30 = int(np.sum(pct_mt > 30.0))
    fail_any = int(np.sum((spot_counts < 500) | (spot_genes < 250) | (pct_mt > 20.0)))
    pct_fail_any = round((fail_any / n_spots) * 100.0, 2) if n_spots > 0 else 0.0

    # Transfer spot metrics to PyTorch GPU for rapid quantiles & acceleration
    counts_t = torch.from_numpy(spot_counts).cuda() if torch.cuda.is_available() else torch.from_numpy(spot_counts)
    genes_t = torch.from_numpy(spot_genes).cuda() if torch.cuda.is_available() else torch.from_numpy(spot_genes)
    mt_t = torch.from_numpy(pct_mt).cuda() if torch.cuda.is_available() else torch.from_numpy(pct_mt)

    # Quantiles
    q_levels = torch.tensor([0.25, 0.5, 0.75], device=counts_t.device)
    c_q = torch.quantile(counts_t, q_levels).cpu().numpy()
    g_q = torch.quantile(genes_t, q_levels).cpu().numpy()
    mt_q = torch.quantile(mt_t, q_levels).cpu().numpy()

    # QC Flags
    flags = []
    if n_spots < 500:
        flags.append("LOW_SPOT_COUNT")
    if c_q[1] < 1000:
        flags.append("LOW_MEDIAN_DEPTH")
    if c_q[1] > 20000:
        flags.append("HIGH_MEDIAN_DEPTH")
    if mt_q[1] > 15.0:
        flags.append("HIGH_MEDIAN_MT")
    qc_flag = "; ".join(flags) if flags else "OK"

    gc.collect()
    mem_after = process.memory_info().rss / (1024 * 1024)
    peak_mem = max(mem_before, mem_after)

    rec = {
        "gsm_accession": gsm,
        "sample_label": sample_label,
        "n_spots": n_spots,
        "n_genes_total": n_features,
        "n_genes_detected": unique_genes_detected,
        "total_umi_counts": int(np.sum(spot_counts)),
        "umi_min": round(float(torch.min(counts_t).item()), 1),
        "umi_q25": round(float(c_q[0]), 1),
        "umi_median": round(float(c_q[1]), 1),
        "umi_mean": round(float(torch.mean(counts_t).item()), 1),
        "umi_q75": round(float(c_q[2]), 1),
        "umi_max": round(float(torch.max(counts_t).item()), 1),
        "genes_per_spot_min": int(torch.min(genes_t).item()),
        "genes_per_spot_q25": round(float(g_q[0]), 1),
        "genes_per_spot_median": round(float(g_q[1]), 1),
        "genes_per_spot_mean": round(float(torch.mean(genes_t).item()), 1),
        "genes_per_spot_q75": round(float(g_q[2]), 1),
        "genes_per_spot_max": int(torch.max(genes_t).item()),
        "pct_mt_min": round(float(torch.min(mt_t).item()), 2),
        "pct_mt_q25": round(float(mt_q[0]), 2),
        "pct_mt_median": round(float(mt_q[1]), 2),
        "pct_mt_mean": round(float(torch.mean(mt_t).item()), 2),
        "pct_mt_q75": round(float(mt_q[2]), 2),
        "pct_mt_max": round(float(torch.max(mt_t).item()), 2),
        "matrix_sparsity": round(float(matrix_sparsity), 4),
        "fail_low_umi_500": fail_low_umi,
        "fail_low_genes_250": fail_low_genes,
        "fail_high_mt_20": fail_high_mt_20,
        "fail_high_mt_30": fail_high_mt_30,
        "fail_any_qc": fail_any,
        "pct_fail_any_qc": pct_fail_any,
        "qc_flag": qc_flag,
        "peak_mem_mb": round(peak_mem, 1)
    }
    qc_summary_records.append(rec)
    cohort_stats.append({
        "gsm": gsm,
        "sample": sample_label,
        "n_spots": n_spots,
        "umi_median": float(c_q[1]),
        "genes_median": float(g_q[1]),
        "mt_median": float(mt_q[1]),
        "fail_pct": pct_fail_any
    })

    # Generate individual sample QC plot (3-panel violin)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    fig.suptitle(f"{gsm} ({sample_label}) - QC Distributions (N = {n_spots} spots)", fontsize=11, fontweight="bold")
    
    # Counts
    sns.violinplot(y=spot_counts, ax=axes[0], color="#4C72B0", inner="quartile")
    axes[0].set_title("Total UMI Counts / Spot", fontsize=10)
    axes[0].set_ylabel("Counts")
    axes[0].axhline(500, color="red", linestyle="--", alpha=0.7, label="Threshold (500)")
    axes[0].legend(loc="upper right", fontsize=8)
    
    # Genes
    sns.violinplot(y=spot_genes, ax=axes[1], color="#55A868", inner="quartile")
    axes[1].set_title("Detected Genes / Spot", fontsize=10)
    axes[1].set_ylabel("Genes")
    axes[1].axhline(250, color="red", linestyle="--", alpha=0.7, label="Threshold (250)")
    axes[1].legend(loc="upper right", fontsize=8)
    
    # MT%
    sns.violinplot(y=pct_mt, ax=axes[2], color="#C44E52", inner="quartile")
    axes[2].set_title("Mitochondrial %", fontsize=10)
    axes[2].set_ylabel("% MT Reads")
    axes[2].axhline(20, color="red", linestyle="--", alpha=0.7, label="Threshold (20%)")
    axes[2].legend(loc="upper right", fontsize=8)
    
    plt.tight_layout()
    fig_path = QC_FIG_DIR / f"{gsm}_qc_dist.png"
    fig.savefig(fig_path)
    plt.close(fig)

    if idx % 10 == 0 or idx == len(h5_files):
        print(f"  Processed QC for {idx}/{len(h5_files)}: {gsm} ({sample_label}) - Spots: {n_spots}, Median UMI: {c_q[1]:.0f}, Fail%: {pct_fail_any}%")

# Save tables/qc_per_sample.csv
df_qc = pd.DataFrame(qc_summary_records)
qc_csv_path = TABLES_DIR / "qc_per_sample.csv"
df_qc.to_csv(qc_csv_path, index=False)
print(f"Saved {qc_csv_path}")

# Barcode sharing summary
bc_counts = list(barcode_tracker.values())
bc_df = pd.DataFrame([
    {"barcode": bc, "samples_sharing_count": count}
    for bc, count in barcode_tracker.items()
])
bc_csv_path = TABLES_DIR / "barcode_sharing_summary.csv"
bc_df.to_csv(bc_csv_path, index=False)
print(f"Saved {bc_csv_path} (Total unique barcodes across cohort: {len(bc_df)})")

# Cohort Summary Plots
cohort_df = pd.DataFrame(cohort_stats)
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle(f"Cohort QC Overview Across 57 Visium Samples (GSE277783)", fontsize=13, fontweight="bold")

# 1. Spot counts per sample
axes[0, 0].bar(range(len(cohort_df)), cohort_df["n_spots"], color="#3470a3", edgecolor="none")
axes[0, 0].set_title("Spot Counts per Sample", fontsize=11)
axes[0, 0].set_ylabel("Number of Spots")
axes[0, 0].set_xlabel("Samples (1-57)")
axes[0, 0].axhline(cohort_df["n_spots"].median(), color="black", linestyle="--", label=f"Cohort Median: {cohort_df['n_spots'].median():.0f}")
axes[0, 0].legend()

# 2. Median UMI counts
axes[0, 1].bar(range(len(cohort_df)), cohort_df["umi_median"], color="#3b9a67", edgecolor="none")
axes[0, 1].set_title("Median UMI Counts per Spot", fontsize=11)
axes[0, 1].set_ylabel("Median UMIs")
axes[0, 1].set_xlabel("Samples (1-57)")
axes[0, 1].axhline(cohort_df["umi_median"].median(), color="black", linestyle="--", label=f"Cohort Median: {cohort_df['umi_median'].median():.0f}")
axes[0, 1].legend()

# 3. Median Genes detected
axes[1, 0].bar(range(len(cohort_df)), cohort_df["genes_median"], color="#e07a5f", edgecolor="none")
axes[1, 0].set_title("Median Detected Genes per Spot", fontsize=11)
axes[1, 0].set_ylabel("Median Genes")
axes[1, 0].set_xlabel("Samples (1-57)")
axes[1, 0].axhline(cohort_df["genes_median"].median(), color="black", linestyle="--", label=f"Cohort Median: {cohort_df['genes_median'].median():.0f}")
axes[1, 0].legend()

# 4. Spots failing QC (%)
axes[1, 1].bar(range(len(cohort_df)), cohort_df["fail_pct"], color="#9b5de5", edgecolor="none")
axes[1, 1].set_title("Fraction of Spots Failing QC Thresholds (<500 UMIs, <250 Genes, >20% MT)", fontsize=11)
axes[1, 1].set_ylabel("% Spots Failing")
axes[1, 1].set_xlabel("Samples (1-57)")
axes[1, 1].axhline(cohort_df["fail_pct"].median(), color="black", linestyle="--", label=f"Cohort Median: {cohort_df['fail_pct'].median():.1f}%")
axes[1, 1].legend()

plt.tight_layout()
cohort_fig_path = FIGURES_DIR / "cohort_qc_summary.png"
cohort_fig_pdf = FIGURES_DIR / "cohort_qc_summary.pdf"
fig.savefig(cohort_fig_path)
fig.savefig(cohort_fig_pdf)
plt.close(fig)
print(f"Saved cohort QC plots: {cohort_fig_path} and {cohort_fig_pdf}")

print(f"[{datetime.now().isoformat()}] Step 3 completed successfully.")

"""
Step 5: Metadata Structure, Patient Cohort Enumeration, and Confounds Analysis.
- Builds patient x site x condition x modality table from verified/parsed filename labels.
- Evaluates patient-level sample sizes for primary vs metastasis comparisons.
- Systematically documents and quantifies confound variables:
  - Pseudoreplication (spots vs patients)
  - Anatomic site label uncertainty (unverified letter suffixes A-E)
  - Extreme sequencing depth variation (114-fold disparity)
  - Slide chemistry / version differences (standard Visium vs CytAssist Visium)
  - PDX xenograft reference disparity (18k vs 36k genes, mouse stroma)
  - Unbalanced sample distribution per patient (3 to 8 samples/patient)
  - Modality restriction (Visium only; CosMx missing from RAW tar)
- Produces:
  - tables/patient_cohort_summary.csv
  - tables/metadata_confounds_matrix.csv
  - figures/cohort_patient_sample_distribution.png
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

OUT_ROOT = Path(r"F:\Thesis\eda_gse277783")
TABLES_DIR = OUT_ROOT / "tables"
FIGURES_DIR = OUT_ROOT / "figures"

TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sample_map = pd.read_csv(TABLES_DIR / "sample_map.csv")
qc_df = pd.read_csv(TABLES_DIR / "qc_per_sample.csv")
spatial_df = pd.read_csv(TABLES_DIR / "spatial_slide_summary.csv")
h5_df = pd.read_csv(TABLES_DIR / "h5_structure_summary.csv")

# Merge sample metadata
merged = sample_map.merge(qc_df[["gsm_accession", "n_spots", "umi_median", "pct_mt_median", "pct_fail_any_qc", "qc_flag"]], on="gsm_accession")
merged = merged.merge(h5_df[["gsm_accession", "n_genes", "file_size_mb"]], on="gsm_accession")
merged = merged.merge(spatial_df[["gsm_accession", "spot_diameter_fullres"]], on="gsm_accession")

# Extract clean patient ID
def clean_pt(label):
    lbl = str(label)
    if "Pt-10" in lbl or "Pt10" in lbl:
        return "Pt-10"
    if "Pt-11" in lbl or "Pt11" in lbl:
        return "Pt-11"
    if "Pt-12" in lbl or "Pt12" in lbl:
        return "Pt-12"
    if "Pt-13" in lbl or "Pt13" in lbl:
        return "Pt-13"
    for i in range(1, 10):
        if f"Pt-{i}" in lbl or f"Pt{i}" in lbl:
            return f"Pt-{i}"
    return "UNKNOWN"

merged["patient_id"] = merged["sample_label_raw"].apply(clean_pt)

# Classify sample condition from filename
def classify_condition(row):
    lbl = str(row["sample_label_raw"])
    if "PDX" in lbl:
        return "PDX Model (Mouse Xenograft)"
    elif "Pri" in lbl:
        return "Primary Tumor (Explicit 'Pri' in filename)"
    else:
        # Suffix letter A, B, C, D, E
        return "Unverified Specimen (Letter suffix A-E; site UNKNOWN)"

merged["specimen_type"] = merged.apply(classify_condition, axis=1)

# Group by patient
patient_groups = []
for pt, grp in merged.groupby("patient_id"):
    n_samples = len(grp)
    total_spots = grp["n_spots"].sum()
    has_explicit_pri = any(grp["specimen_type"] == "Primary Tumor (Explicit 'Pri' in filename)")
    has_pdx = any(grp["specimen_type"] == "PDX Model (Mouse Xenograft)")
    unverified_letters = [lbl.split("_")[-1].replace(" (from filename, unverified)", "") for lbl in grp[grp["specimen_type"] == "Unverified Specimen (Letter suffix A-E; site UNKNOWN)"]["sample_label_raw"]]
    
    patient_groups.append({
        "patient_id": pt,
        "total_visium_slides": n_samples,
        "total_spots": total_spots,
        "has_explicit_primary_in_filename": has_explicit_pri,
        "has_explicit_met_in_filename": False,  # 0 in filename
        "has_pdx_model": has_pdx,
        "specimen_labels": ", ".join(sorted(unverified_letters)),
        "explicit_pri_count": sum(grp["specimen_type"] == "Primary Tumor (Explicit 'Pri' in filename)"),
        "unverified_letter_count": len(unverified_letters),
        "pdx_count": sum(grp["specimen_type"] == "PDX Model (Mouse Xenograft)"),
        "median_umi_range": f"{grp['umi_median'].min():.0f} - {grp['umi_median'].max():.0f}",
        "spots_per_slide_range": f"{grp['n_spots'].min()} - {grp['n_spots'].max()}"
    })

df_patient = pd.DataFrame(patient_groups)
patient_summary_path = TABLES_DIR / "patient_cohort_summary.csv"
df_patient.to_csv(patient_summary_path, index=False)
print(f"Saved {patient_summary_path}")

# Build Confound Matrix Table
confounds_records = [
    {
        "confound_id": "CONF-01",
        "confound_category": "Sample Size & Pseudoreplication",
        "description": "Treating spots (N = 149,088) as independent units inflated N by >11,000-fold. True patient cohort size is N = 13 human patients.",
        "severity": "CRITICAL",
        "impact_on_analysis": "Severe type I error inflation if pseudoreplication is unaddressed. Linear mixed models (LMM) with patient as random effect are mandatory.",
        "affected_samples": "All 57 samples (100%)"
    },
    {
        "confound_id": "CONF-02",
        "confound_category": "Anatomic Site Label Uncertainty",
        "description": "49 of 57 samples (86%) only contain letter codes (A-E) in filenames. No primary vs metastatic organ metadata exists in the RAW tar archive.",
        "severity": "CRITICAL",
        "impact_on_analysis": "Primary vs metastasis comparison is IMPOSSIBLE without external GEO Series Matrix / sample metadata sheet.",
        "affected_samples": "49 samples (Pt-1A through Pt-13E)"
    },
    {
        "confound_id": "CONF-03",
        "confound_category": "Sequencing Depth Heterogeneity",
        "description": "114.2-fold difference in median UMI depth across slides (min = 130 UMIs/spot in GSM8452848; max = 14,841 UMIs/spot in GSM8452857).",
        "severity": "HIGH",
        "impact_on_analysis": "Differential expression will capture library size artifacts rather than biological differences unless size-factor normalized.",
        "affected_samples": "6 flagged low-depth slides (<1000 UMIs), highly variable across all 57"
    },
    {
        "confound_id": "CONF-04",
        "confound_category": "Slide Chemistry & Resolution Batch Effect",
        "description": "8 samples (GSM8714177-GSM8714184) utilized Visium CytAssist / v2 architecture with 81px/280px spot scalefactors and spatial enrichment files, whereas 49 samples utilized standard Visium v1 (~179px spot diameter).",
        "severity": "HIGH",
        "impact_on_analysis": "Slide capture chemistry and optical resolution differ systematically between the replicate set and main patient set.",
        "affected_samples": "8 samples (GSM8714177 - GSM8714184)"
    },
    {
        "confound_id": "CONF-05",
        "confound_category": "Species & Reference Disparity (PDX Xenografts)",
        "description": "GSM8714183 (PDXPt10) and GSM8714184 (PDXPt13) are mouse xenograft models aligned against an 18,085 gene human reference, while patient samples have 36,601 genes. Mouse stromal contamination and human probe depletion exist.",
        "severity": "HIGH",
        "impact_on_analysis": "Cannot be directly pooled with primary patient tissue without gene-space subsetting and xenograft-specific stromal deconvolution.",
        "affected_samples": "2 samples (GSM8714183, GSM8714184)"
    },
    {
        "confound_id": "CONF-06",
        "confound_category": "Imbalanced Patient Representation",
        "description": "Slide count varies widely per patient from 3 slides (Pt-1, Pt-2, Pt-5) to 8 slides (Pt-13). Pt-13 represents 14% of all spatial spots.",
        "severity": "MEDIUM",
        "impact_on_analysis": "Pooled analyses will be heavily biased towards patient Pt-13's idiosyncratic tumor biology.",
        "affected_samples": "Pt-13 (8 slides, 19,850 spots)"
    },
    {
        "confound_id": "CONF-07",
        "confound_category": "Modality Incompleteness in Download",
        "description": "NanoString CosMx SMI data is completely absent from GSE277783_RAW.tar (0 CosMx files). Only 10x Visium is present.",
        "severity": "HIGH",
        "impact_on_analysis": "Single-cell spatial resolution and CosMx-specific validation cannot be performed from this archive alone.",
        "affected_samples": "Entire series (requires separate download of CosMx subseries e.g. GSE277782 / GSE310352)"
    }
]

df_confounds = pd.DataFrame(confounds_records)
confounds_path = TABLES_DIR / "metadata_confounds_matrix.csv"
df_confounds.to_csv(confounds_path, index=False)
print(f"Saved {confounds_path}")

# Plot Patient Sample Distribution
plt.figure(figsize=(10, 5))
patients = df_patient["patient_id"]
x = np.arange(len(patients))
width = 0.45

p1 = plt.bar(x, df_patient["unverified_letter_count"], width, label="Unverified Letter Specimen (A-E)", color="#4575b4")
p2 = plt.bar(x, df_patient["explicit_pri_count"], width, bottom=df_patient["unverified_letter_count"], label="Explicit Primary ('Pri2/3')", color="#d73027")
p3 = plt.bar(x, df_patient["pdx_count"], width, bottom=df_patient["unverified_letter_count"] + df_patient["explicit_pri_count"], label="PDX Xenograft", color="#fdae61")

plt.xlabel("Patient ID", fontsize=11, fontweight="bold")
plt.ylabel("Number of Visium Slides", fontsize=11, fontweight="bold")
plt.title("Sample Distribution per Patient across GSE277783 (Total = 57 Slides)", fontsize=12, fontweight="bold")
plt.xticks(x, patients, rotation=45, ha="right")
plt.legend(loc="upper left", fontsize=9)
plt.tight_layout()

patient_fig_path = FIGURES_DIR / "cohort_patient_sample_distribution.png"
patient_fig_pdf = FIGURES_DIR / "cohort_patient_sample_distribution.pdf"
plt.savefig(patient_fig_path)
plt.savefig(patient_fig_pdf)
plt.close()
print(f"Saved patient sample distribution plot: {patient_fig_path} and {patient_fig_pdf}")

print("Step 5 completed successfully.")

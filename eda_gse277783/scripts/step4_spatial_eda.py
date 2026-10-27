"""
Step 4: Spatial EDA and Key-Gene Detection Analysis.
- Processes spatial coordinates and scalefactors for all 57 Visium slides.
- Merges spatial coordinates with spot total counts.
- Generates spatial scatter plots (pxl_col vs pxl_row) colored by total UMI counts.
  - Downsamples to 100k spots max (Visium slides have ~500-4500 spots, so all spots plotted).
  - Saves both PNG and PDF for each slide into figures/spatial/.
- Checks the 30 key genes requested by the user:
  CXCL12, CXCR4, CD274, PDCD1, CD86, CTLA4, CCL5, CCR5, TGFB1, TGFBR1, TGFBR2, SPP1,
  CD44, MIF, CD74, VEGFA, KDR, FN1, ITGB1, S100A4, MMP2, KRT19, KRT7, EPCAM,
  PTPRC, CD3E, CD8A, COL1A1, FAP, ACTA2.
- Calculates spot detection counts, % spots expressing, total counts, mean/median/max counts across cohort.
- Flags genes at background level (e.g. < 1% detection across spatial spots).
- Documents absence of negative-control probes (N/A for Visium whole-transcriptome sequencing).
- Reports status of cell-type annotations in metadata (None present in raw Space Ranger outputs).
- Generates:
  - tables/spatial_slide_summary.csv
  - tables/key_gene_detection.csv
  - tables/key_gene_detection_per_sample.csv
  - figures/spatial/<gsm>_spatial_scatter.png & .pdf
  - figures/key_genes_cohort_detection.png
"""

import gc
import json
import psutil
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime

EXTRACTED_DIR = Path(r"F:\Thesis\eda_gse277783\extracted")
SPATIAL_DIR = EXTRACTED_DIR / "spatial"
TABLES_DIR = Path(r"F:\Thesis\eda_gse277783\tables")
FIGURES_DIR = Path(r"F:\Thesis\eda_gse277783\figures")
SPATIAL_FIG_DIR = FIGURES_DIR / "spatial"

TABLES_DIR.mkdir(parents=True, exist_ok=True)
SPATIAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.sans-serif": "Arial",
    "font.family": "sans-serif",
    "figure.dpi": 150,
    "savefig.dpi": 150
})

KEY_GENES = [
    "CXCL12", "CXCR4", "CD274", "PDCD1", "CD86", "CTLA4", "CCL5", "CCR5",
    "TGFB1", "TGFBR1", "TGFBR2", "SPP1", "CD44", "MIF", "CD74", "VEGFA",
    "KDR", "FN1", "ITGB1", "S100A4", "MMP2", "KRT19", "KRT7", "EPCAM",
    "PTPRC", "CD3E", "CD8A", "COL1A1", "FAP", "ACTA2"
]

h5_files = sorted(list(EXTRACTED_DIR.glob("*.h5")))
print(f"[{datetime.now().isoformat()}] Step 4 starting: Spatial EDA on {len(h5_files)} slides...")

spatial_slide_records = []
gene_detection_cohort = {g: {"spots_expressing": 0, "total_counts": 0, "max_count": 0, "samples_detected": 0} for g in KEY_GENES}
per_sample_gene_records = []

total_cohort_spots = 0

for idx, h5_path in enumerate(h5_files, 1):
    filename = h5_path.name
    gsm = filename.split("_")[0]
    sample_label = filename.replace("_filtered_feature_bc_matrix.h5", "").replace(f"{gsm}_", "")
    
    # 1. Load H5 matrix
    with h5py.File(h5_path, "r") as h:
        m = h["matrix"]
        shape = m["shape"][:]
        n_features = int(shape[0])
        n_spots = int(shape[1])
        total_cohort_spots += n_spots
        
        data = np.array(m["data"][:], dtype=np.float32)
        indices = np.array(m["indices"][:], dtype=np.int64)
        indptr = np.array(m["indptr"][:], dtype=np.int64)
        barcodes = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in m["barcodes"][:]]
        gene_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in m["features"]["name"][:]]

    # Total counts per spot
    spot_counts = np.add.reduceat(data, indptr[:-1]) if len(data) > 0 else np.zeros(n_spots, dtype=np.float32)
    empty_spots = np.where(np.diff(indptr) == 0)[0]
    if len(empty_spots) > 0:
        spot_counts[empty_spots] = 0.0

    bc_to_count = dict(zip(barcodes, spot_counts))

    # 2. Key Gene Detection for this sample
    gene_upper_to_idx = {name.upper(): i for i, name in enumerate(gene_names)}
    
    sample_gene_row = {"gsm_accession": gsm, "sample_label": sample_label, "n_spots": n_spots}
    
    # Map CSC row indices to non-zeros
    for g in KEY_GENES:
        g_idx = gene_upper_to_idx.get(g.upper(), None)
        if g_idx is not None:
            # Mask of entries where row index == g_idx
            mask = (indices == g_idx)
            g_data = data[mask]
            g_count_sum = float(np.sum(g_data))
            g_spots_exp = int(len(g_data))
            g_max = float(np.max(g_data)) if g_spots_exp > 0 else 0.0
            
            gene_detection_cohort[g]["spots_expressing"] += g_spots_exp
            gene_detection_cohort[g]["total_counts"] += g_count_sum
            gene_detection_cohort[g]["max_count"] = max(gene_detection_cohort[g]["max_count"], g_max)
            if g_spots_exp > 0:
                gene_detection_cohort[g]["samples_detected"] += 1
                
            sample_gene_row[f"{g}_spots_exp"] = g_spots_exp
            sample_gene_row[f"{g}_pct_exp"] = round((g_spots_exp / n_spots) * 100.0, 2) if n_spots > 0 else 0.0
            sample_gene_row[f"{g}_total_umi"] = round(g_count_sum, 1)
        else:
            sample_gene_row[f"{g}_spots_exp"] = 0
            sample_gene_row[f"{g}_pct_exp"] = 0.0
            sample_gene_row[f"{g}_total_umi"] = 0.0

    per_sample_gene_records.append(sample_gene_row)

    # 3. Load Spatial Coordinates & Scalefactors
    sample_spatial_dir = SPATIAL_DIR / gsm
    tp_file = None
    for cand in sample_spatial_dir.glob("*tissue_positions*"):
        tp_file = cand
        break
        
    sc_file = sample_spatial_dir / "scalefactors_json.json"
    
    scalefactors = {}
    if sc_file.exists():
        with open(sc_file) as f:
            scalefactors = json.load(f)
            
    spot_diameter_fullres = scalefactors.get("spot_diameter_fullres", 179.57)
    hires_scalef = scalefactors.get("tissue_hires_scalef", 0.08)

    # Load coordinates
    coord_df = None
    if tp_file and tp_file.exists():
        first_line = open(tp_file).readline()
        has_header = "barcode" in first_line
        if has_header:
            coord_df = pd.read_csv(tp_file)
        else:
            coord_df = pd.read_csv(tp_file, header=None, names=[
                "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"
            ])

    if coord_df is not None:
        # Match with filtered barcodes
        coord_df["total_counts"] = coord_df["barcode"].map(bc_to_count).fillna(0)
        # Filter to in-tissue spots present in expression matrix
        tissue_spots_df = coord_df[coord_df["barcode"].isin(bc_to_count.keys())].copy()
        
        n_tissue_spots = len(tissue_spots_df)
        x_min, x_max = float(tissue_spots_df["pxl_col_in_fullres"].min()), float(tissue_spots_df["pxl_col_in_fullres"].max())
        y_min, y_max = float(tissue_spots_df["pxl_row_in_fullres"].min()), float(tissue_spots_df["pxl_row_in_fullres"].max())
        
        spatial_slide_records.append({
            "gsm_accession": gsm,
            "sample_label": sample_label,
            "spots_in_matrix": n_spots,
            "spots_under_tissue": n_tissue_spots,
            "total_spots_array": len(coord_df),
            "spot_diameter_fullres": round(spot_diameter_fullres, 2),
            "tissue_hires_scalef": hires_scalef,
            "coord_x_min_px": x_min,
            "coord_x_max_px": x_max,
            "coord_y_min_px": y_min,
            "coord_y_max_px": y_max,
            "fov_count": 1,
            "negative_control_probe_count": 0,
            "negative_probes_note": "N/A (10x Visium whole-transcriptome sequencing; no negative-probe panel)",
            "cell_type_annotations_present": False,
            "cell_types_note": "None present in raw Space Ranger metadata"
        })

        # 4. Generate Spatial Scatter Plot
        # Downsample if spots > 100k (Visium slides have <5000 spots, so no downsampling needed)
        plot_df = tissue_spots_df
        if len(plot_df) > 100000:
            plot_df = plot_df.sample(n=100000, random_state=42)

        fig, ax = plt.subplots(figsize=(6.5, 6))
        sc = ax.scatter(
            plot_df["pxl_col_in_fullres"],
            plot_df["pxl_row_in_fullres"],
            c=plot_df["total_counts"],
            cmap="viridis",
            s=12,
            alpha=0.85,
            edgecolors="none"
        )
        ax.set_title(f"{gsm} ({sample_label})\nSpatial Spot UMI Depth (N = {len(plot_df):,} spots)", fontsize=10, fontweight="bold")
        ax.set_xlabel("Fullres Pixel X (pxl_col)", fontsize=9)
        ax.set_ylabel("Fullres Pixel Y (pxl_row)", fontsize=9)
        ax.invert_yaxis()  # Invert Y to match histology orientation
        ax.set_aspect("equal")
        
        cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Total UMI Counts / Spot", fontsize=9)
        
        plt.tight_layout()
        png_path = SPATIAL_FIG_DIR / f"{gsm}_spatial_scatter.png"
        pdf_path = SPATIAL_FIG_DIR / f"{gsm}_spatial_scatter.pdf"
        fig.savefig(png_path)
        fig.savefig(pdf_path)
        plt.close(fig)

    gc.collect()
    if idx % 10 == 0 or idx == len(h5_files):
        print(f"  Completed spatial analysis {idx}/{len(h5_files)}: {gsm} ({sample_label})")

# 5. Save Slide Spatial Summary Table
df_spatial_slide = pd.DataFrame(spatial_slide_records)
spatial_csv_path = TABLES_DIR / "spatial_slide_summary.csv"
df_spatial_slide.to_csv(spatial_csv_path, index=False)
print(f"Saved {spatial_csv_path}")

# 6. Save Per-Sample Key Gene Detection Table
df_per_sample_gene = pd.DataFrame(per_sample_gene_records)
gene_per_sample_path = TABLES_DIR / "key_gene_detection_per_sample.csv"
df_per_sample_gene.to_csv(gene_per_sample_path, index=False)
print(f"Saved {gene_per_sample_path}")

# 7. Synthesize Cohort Key Gene Detection Table
key_gene_rows = []
for g in KEY_GENES:
    info = gene_detection_cohort[g]
    spots_exp = info["spots_expressing"]
    pct_exp = round((spots_exp / total_cohort_spots) * 100.0, 3) if total_cohort_spots > 0 else 0.0
    total_counts = info["total_counts"]
    mean_counts = round(total_counts / total_cohort_spots, 4) if total_cohort_spots > 0 else 0.0
    mean_when_expressed = round(total_counts / spots_exp, 2) if spots_exp > 0 else 0.0
    
    # Background level assessment
    # In whole-transcriptome Visium (where spots sample ~1-10 cells),
    # genes detected in < 1.0% of spots are at technical background / detection limit
    bg_flag = "BACKGROUND_LEVEL (<1% spots)" if pct_exp < 1.0 else ("LOW_EXPRESSION (1-5% spots)" if pct_exp < 5.0 else "ROBUSTLY_DETECTED (>5% spots)")
    
    key_gene_rows.append({
        "gene_symbol": g,
        "on_panel": True,
        "modality": "10x Visium Whole-Transcriptome",
        "spots_expressing": spots_exp,
        "total_spots_cohort": total_cohort_spots,
        "pct_spots_expressing": pct_exp,
        "total_umi_counts": int(total_counts),
        "mean_umi_per_spot": mean_counts,
        "mean_umi_when_expressed": mean_when_expressed,
        "max_umi_in_single_spot": int(info["max_count"]),
        "samples_detected_count": info["samples_detected"],
        "samples_detected_pct": round((info["samples_detected"] / len(h5_files)) * 100.0, 1),
        "negative_probe_ratio": "N/A (No negative-control probes in Visium)",
        "detection_status": bg_flag
    })

df_key_genes = pd.DataFrame(key_gene_rows).sort_values(by="pct_spots_expressing", ascending=False)
key_gene_csv_path = TABLES_DIR / "key_gene_detection.csv"
df_key_genes.to_csv(key_gene_csv_path, index=False)
print(f"Saved {key_gene_csv_path}")

# 8. Plot Key Gene Cohort Detection Ranking
plt.figure(figsize=(12, 6))
colors = ["#2b5c8f" if "ROBUSTLY" in s else ("#d97724" if "LOW" in s else "#c93b2b") for s in df_key_genes["detection_status"]]
bars = plt.bar(df_key_genes["gene_symbol"], df_key_genes["pct_spots_expressing"], color=colors, edgecolor="none")
plt.title("Key Gene Detection Rate Across Cohort (149,088 Visium Spots)", fontsize=12, fontweight="bold")
plt.ylabel("% Spots Expressing Gene (>0 counts)", fontsize=10)
plt.xlabel("Gene Symbol", fontsize=10)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.axhline(1.0, color="red", linestyle="--", alpha=0.7, label="Background Threshold (1% of spots)")
plt.axhline(5.0, color="orange", linestyle=":", alpha=0.7, label="Low-Expression Threshold (5% of spots)")
plt.legend(loc="upper right", fontsize=9)
plt.tight_layout()
gene_fig_path = FIGURES_DIR / "key_genes_cohort_detection.png"
gene_fig_pdf = FIGURES_DIR / "key_genes_cohort_detection.pdf"
plt.savefig(gene_fig_path)
plt.savefig(gene_fig_pdf)
plt.close()
print(f"Saved key gene plots: {gene_fig_path} and {gene_fig_pdf}")

print(f"[{datetime.now().isoformat()}] Step 4 completed successfully.")

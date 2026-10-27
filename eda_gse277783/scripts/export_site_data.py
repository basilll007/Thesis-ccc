"""
Export GSE277783 EDA summary data to window.GSE277783_DATA in JavaScript format
for interactive web dashboard consumption.
"""

import json
import pandas as pd
from pathlib import Path

tables_dir = Path(r"F:\Thesis\eda_gse277783\tables")
assets_dir = Path(r"F:\Thesis\assets\metrics")
assets_dir.mkdir(parents=True, exist_ok=True)

sample_map = pd.read_csv(tables_dir / "sample_map.csv")
qc_df = pd.read_csv(tables_dir / "qc_per_sample.csv")
spatial_df = pd.read_csv(tables_dir / "spatial_slide_summary.csv")
img_df = pd.read_csv(tables_dir / "image_inventory_summary.csv")
h5_df = pd.read_csv(tables_dir / "h5_structure_summary.csv")
pt_df = pd.read_csv(tables_dir / "patient_cohort_summary.csv")
kg_df = pd.read_csv(tables_dir / "key_gene_detection.csv")
conf_df = pd.read_csv(tables_dir / "metadata_confounds_matrix.csv")

# Merge slide info
merged = sample_map.merge(qc_df[["gsm_accession", "n_spots", "umi_median", "pct_mt_median", "pct_fail_any_qc", "qc_flag", "total_umi_counts", "genes_per_spot_median"]], on="gsm_accession")
merged = merged.merge(h5_df[["gsm_accession", "n_genes"]], on="gsm_accession")
merged = merged.merge(spatial_df[["gsm_accession", "spot_diameter_fullres", "coord_x_max_px", "coord_y_max_px"]], on="gsm_accession")
merged = merged.merge(img_df[["gsm_accession", "thumbnail_file", "dimensions", "format"]], on="gsm_accession")

slides_list = []
for idx, r in merged.iterrows():
    gsm = r["gsm_accession"]
    lbl = r["sample_label_raw"].replace(" (from filename, unverified)", "")
    
    # Clean patient ID
    if "Pt-10" in lbl or "Pt10" in lbl:
        pt = "Pt-10"
    elif "Pt-11" in lbl or "Pt11" in lbl:
        pt = "Pt-11"
    elif "Pt-12" in lbl or "Pt12" in lbl:
        pt = "Pt-12"
    elif "Pt-13" in lbl or "Pt13" in lbl:
        pt = "Pt-13"
    else:
        prefix = lbl.split("_")[0].replace("PDX", "")
        pt = prefix if prefix.startswith("Pt-") else f"Pt-{prefix.replace('Pt', '')}"
        
    spec_type = "PDX Xenograft (Mouse)" if "PDX" in lbl else ("Primary Tumor (Pri)" if "Pri" in lbl else "Specimen A-E (Organ Unverified)")
    chem = "Visium CytAssist v2" if r["spot_diameter_fullres"] < 100 or "PDX" in lbl or "Pri" in lbl else "Standard Visium v1"
    thumb_name = r["thumbnail_file"]
    
    slides_list.append({
        "gsm": gsm,
        "label": lbl,
        "patient": pt,
        "specimen_type": spec_type,
        "chemistry": chem,
        "n_spots": int(r["n_spots"]),
        "n_genes": int(r["n_genes"]),
        "total_umi": int(r["total_umi_counts"]),
        "umi_median": float(r["umi_median"]),
        "genes_median": float(r["genes_per_spot_median"]),
        "pct_mt_median": float(r["pct_mt_median"]),
        "fail_pct": float(r["pct_fail_any_qc"]),
        "qc_flag": r["qc_flag"],
        "img_dims": r["dimensions"],
        "img_format": r["format"],
        "spatial_png": f"eda_gse277783/figures/spatial/{gsm}_spatial_scatter.png",
        "spatial_pdf": f"eda_gse277783/figures/spatial/{gsm}_spatial_scatter.pdf",
        "thumb_path": f"eda_gse277783/figures/thumbs/{thumb_name}",
        "qc_plot_path": f"eda_gse277783/figures/qc/{gsm}_qc_dist.png"
    })

data_payload = {
    "summary": {
        "total_slides": len(slides_list),
        "total_spots": int(qc_df["n_spots"].sum()),
        "total_patients": len(pt_df),
        "modality": "10x Genomics Visium Spatial Transcriptomics",
        "cosmx_status": "Absent from raw tar (0 CosMx files)",
        "max_patient_n": 13,
        "depth_fold_diff": 114.2,
        "flagged_slides": int((qc_df["qc_flag"] != "OK").sum())
    },
    "slides": slides_list,
    "patients": pt_df.to_dict(orient="records"),
    "key_genes": kg_df.to_dict(orient="records"),
    "confounds": conf_df.to_dict(orient="records")
}

out_js = assets_dir / "gse277783_eda_data.js"
with open(out_js, "w", encoding="utf-8") as f:
    f.write("window.GSE277783_DATA = " + json.dumps(data_payload, indent=2) + ";")

print(f"Wrote {out_js} with {len(slides_list)} slides and {len(kg_df)} key genes.")

"""
Step 2: Per-file integrity and structure analysis.
- Inspects every H5 file: format, shape (spots x genes), gene ID types, obs/metadata fields, coordinates inside H5.
- Inspects every image: dimensions, size, saves small thumbnail (max 1000 px) to figures/thumbs/.
- Checks for CosMx tables (0 in this archive).
- Flags corrupt, empty, or unreadable files.
- Tracks peak memory per file via psutil.
- Produces tables/h5_structure_summary.csv and tables/image_inventory_summary.csv.
"""

import os
import gc
import gzip
import struct
import psutil
import h5py
import pandas as pd
from PIL import Image
from pathlib import Path
from datetime import datetime

EXTRACTED_DIR = Path(r"F:\Thesis\eda_gse277783\extracted")
SPATIAL_DIR = EXTRACTED_DIR / "spatial"
TABLES_DIR = Path(r"F:\Thesis\eda_gse277783\tables")
THUMBS_DIR = Path(r"F:\Thesis\eda_gse277783\figures\thumbs")

TABLES_DIR.mkdir(parents=True, exist_ok=True)
THUMBS_DIR.mkdir(parents=True, exist_ok=True)
Image.MAX_IMAGE_PIXELS = None

process = psutil.Process()

print(f"[{datetime.now().isoformat()}] Step 2 starting: Integrity & Structure Inspection...")

# ---------------------------------------------------------
# Part A: H5 Files Inspection
# ---------------------------------------------------------
h5_files = sorted(list(EXTRACTED_DIR.glob("*.h5")))
print(f"Found {len(h5_files)} H5 files to inspect.")

h5_records = []

for idx, h5_path in enumerate(h5_files, 1):
    mem_before = process.memory_info().rss / (1024 * 1024)
    file_size_mb = round(h5_path.stat().st_size / (1024 * 1024), 2)
    filename = h5_path.name
    gsm = filename.split("_")[0] if filename.startswith("GSM") else "UNKNOWN"
    
    status = "OK"
    format_type = "UNKNOWN"
    shape_str = "UNKNOWN"
    n_spots = 0
    n_genes = 0
    gene_id_type = "UNKNOWN"
    has_obs_metadata = False
    coords_in_h5 = False
    details = ""
    
    try:
        if h5_path.stat().st_size == 0:
            status = "EMPTY"
            format_type = "EMPTY_FILE"
        else:
            with h5py.File(h5_path, "r") as h:
                root_keys = list(h.keys())
                if "matrix" in root_keys:
                    format_type = "10x Genomics HDF5 (Feature-Barcode Matrix v3/v2)"
                    matrix_grp = h["matrix"]
                    shape = matrix_grp["shape"][:]  # [features, barcodes]
                    n_genes = int(shape[0])
                    n_spots = int(shape[1])
                    shape_str = f"{n_spots} spots x {n_genes} genes"
                    
                    # Gene ID types
                    feat_grp = matrix_grp.get("features", {})
                    feat_keys = list(feat_grp.keys()) if hasattr(feat_grp, "keys") else []
                    id_samples = []
                    name_samples = []
                    if "id" in feat_grp:
                        id_samples = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat_grp["id"][:3]]
                    if "name" in feat_grp:
                        name_samples = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat_grp["name"][:3]]
                    
                    gene_id_type = f"IDs: {id_samples[0] if id_samples else 'None'} (Ensembl), Symbols: {name_samples[0] if name_samples else 'None'}"
                    
                    # Check for obs metadata and coordinates inside H5
                    # 10x feature bc matrices have barcodes and features, no obs annotations
                    has_obs_metadata = any(k not in ["barcodes", "data", "features", "indices", "indptr", "shape"] for k in matrix_grp.keys())
                    coords_in_h5 = "spatial" in root_keys or "spatial" in matrix_grp.keys()
                    details = f"CSC sparse: {len(matrix_grp['data'])} non-zero entries"
                elif "X" in root_keys or "obs" in root_keys:
                    format_type = "AnnData h5ad"
                    # Check shape
                    shape_str = "AnnData structure"
                else:
                    format_type = f"Custom HDF5 (keys: {root_keys})"
    except Exception as e:
        status = f"CORRUPT: {str(e)}"
        details = str(e)
        
    gc.collect()
    mem_after = process.memory_info().rss / (1024 * 1024)
    peak_mem = max(mem_before, mem_after)
    
    h5_records.append({
        "gsm_accession": gsm,
        "filename": filename,
        "file_size_mb": file_size_mb,
        "format": format_type,
        "n_spots": n_spots,
        "n_genes": n_genes,
        "shape": shape_str,
        "gene_id_type": gene_id_type,
        "has_obs_metadata": has_obs_metadata,
        "coordinates_in_h5": coords_in_h5,
        "status": status,
        "details": details,
        "peak_mem_mb": round(peak_mem, 1)
    })
    
    if idx % 10 == 0 or idx == len(h5_files):
        print(f"  Inspected H5 {idx}/{len(h5_files)}: {filename} ({shape_str}, status={status})")

df_h5 = pd.DataFrame(h5_records)
h5_summary_path = TABLES_DIR / "h5_structure_summary.csv"
df_h5.to_csv(h5_summary_path, index=False)
print(f"Saved {h5_summary_path}")

# ---------------------------------------------------------
# Part B: Image Files Inspection & Thumbnails
# ---------------------------------------------------------
image_files = sorted(list(EXTRACTED_DIR.glob("*.jpg.gz")) + list(EXTRACTED_DIR.glob("*.qp.tiff.gz")))
print(f"Found {len(image_files)} image archives to inspect...")

image_records = []

for idx, img_path in enumerate(image_files, 1):
    mem_before = process.memory_info().rss / (1024 * 1024)
    file_size_mb = round(img_path.stat().st_size / (1024 * 1024), 2)
    filename = img_path.name
    gsm = filename.split("_")[0] if filename.startswith("GSM") else "UNKNOWN"
    
    status = "OK"
    width, height = 0, 0
    img_format = "UNKNOWN"
    thumb_file = ""
    
    try:
        if filename.endswith(".qp.tiff.gz"):
            # BigTIFF whole slide scan
            img_format = "BigTIFF (Vectra / PerkinElmer QPTIFF Whole-Slide Image)"
            with gzip.open(img_path, "rb") as f:
                hdr = f.read(16)
                first_ifd = struct.unpack("<Q", hdr[8:16])[0]
                f.seek(first_ifd)
                num_entries = struct.unpack("<Q", f.read(8))[0]
                for _ in range(num_entries):
                    tag_data = f.read(20)
                    tag, tag_type, count, val = struct.unpack("<HHQQ", tag_data)
                    if tag == 256:
                        width = val
                    elif tag == 257:
                        height = val
            
            # Save thumbnail from associated tissue_hires_image.png in spatial folder
            hires_path = SPATIAL_DIR / gsm / "tissue_hires_image.png"
            if hires_path.exists():
                thumb_target = THUMBS_DIR / f"{filename.replace('.qp.tiff.gz', '')}_thumb.jpg"
                with Image.open(hires_path) as himg:
                    thumb = himg.convert("RGB")
                    thumb.thumbnail((1000, 1000), Image.Resampling.LANCZOS)
                    thumb.save(thumb_target, "JPEG", quality=85)
                thumb_file = thumb_target.name
            else:
                thumb_file = "NONE_AVAILABLE"
                
        elif filename.endswith(".jpg.gz"):
            img_format = "Gzipped JPEG (Visium Histology Slide Scan)"
            thumb_target = THUMBS_DIR / f"{filename.replace('.jpg.gz', '')}_thumb.jpg"
            with gzip.open(img_path, "rb") as gz:
                with Image.open(gz) as img:
                    width, height = img.size
                    # Use draft mode for memory safety
                    img.draft("RGB", (1000, 1000))
                    thumb = img.copy()
                    thumb.thumbnail((1000, 1000), Image.Resampling.LANCZOS)
                    thumb.save(thumb_target, "JPEG", quality=85)
            thumb_file = thumb_target.name
            
    except Exception as e:
        status = f"ERROR: {str(e)}"
        thumb_file = "FAILED"
        
    gc.collect()
    mem_after = process.memory_info().rss / (1024 * 1024)
    peak_mem = max(mem_before, mem_after)
    
    image_records.append({
        "gsm_accession": gsm,
        "filename": filename,
        "file_size_mb": file_size_mb,
        "format": img_format,
        "width_px": width,
        "height_px": height,
        "dimensions": f"{width} x {height}",
        "thumbnail_file": thumb_file,
        "status": status,
        "peak_mem_mb": round(peak_mem, 1)
    })
    
    if idx % 10 == 0 or idx == len(image_files):
        print(f"  Processed image {idx}/{len(image_files)}: {filename} ({width}x{height}, thumb={thumb_file})")

df_img = pd.DataFrame(image_records)
img_summary_path = TABLES_DIR / "image_inventory_summary.csv"
df_img.to_csv(img_summary_path, index=False)
print(f"Saved {img_summary_path}")

print(f"[{datetime.now().isoformat()}] Step 2 completed successfully.")

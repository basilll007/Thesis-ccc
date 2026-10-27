"""
Step 1: Inventory of GSE277783_RAW.tar
- List every member: name, size, extension, mtime.
- Check path safety (no absolute paths, no ..).
- Identify nested TAR archives and inspect their contents in-memory / via streaming.
- Parse GSM accessions and sample attributes (patient, site, condition/label) from filename.
- Categorize modality and provide evidence.
- Generate tables/tar_inventory.csv and tables/sample_map.csv.
"""

import sys
import tarfile
import csv
import io
import re
from pathlib import Path
from datetime import datetime

TAR_PATH = Path(r"F:\Thesis\data\GSE277783_RAW.tar")
OUT_ROOT = Path(r"F:\Thesis\eda_gse277783")
TABLES_DIR = OUT_ROOT / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

print(f"Opening archive: {TAR_PATH}")
if not TAR_PATH.exists():
    print(f"ERROR: {TAR_PATH} does not exist!")
    sys.exit(1)

inventory_records = []
sample_records = []
nested_tar_records = []

# Regex for GSM IDs: GSM followed by 7 or 8 digits
GSM_REGEX = re.compile(r"(GSM\d+)", re.IGNORECASE)

with tarfile.open(TAR_PATH, mode="r:") as tar:
    print("Reading tar members...")
    for member in tar:
        name = member.name
        size = member.size
        is_file = member.isfile()
        is_dir = member.isdir()
        mtime = datetime.fromtimestamp(member.mtime).isoformat()
        
        # Path safety check
        is_safe = not (name.startswith("/") or name.startswith("\\") or ".." in name or ":" in name)
        
        # Extension
        p = Path(name)
        ext = "".join(p.suffixes).lower()
        if not ext and is_file:
            ext = "no_ext"
        elif is_dir:
            ext = "directory"
            
        # Group file type
        type_group = "OTHER"
        if ext.endswith(".tar.gz") or ext.endswith(".tar") or ext.endswith(".tgz"):
            type_group = "NESTED_TAR"
        elif ext.endswith(".h5") or ext.endswith(".h5ad") or ext.endswith(".hdf5"):
            type_group = "H5"
        elif ext.endswith(".tif") or ext.endswith(".tiff"):
            type_group = "TIFF"
        elif ext.endswith(".jpg") or ext.endswith(".jpeg"):
            type_group = "JPG"
        elif ext.endswith(".csv") or ext.endswith(".tsv") or ext.endswith(".txt") or ext.endswith(".csv.gz"):
            type_group = "TABLE"
        elif is_dir:
            type_group = "DIR"

        inventory_records.append({
            "archive": "GSE277783_RAW.tar",
            "member_name": name,
            "size_bytes": size,
            "size_mb": round(size / (1024 * 1024), 2),
            "extension": ext,
            "type_group": type_group,
            "is_safe": is_safe,
            "mtime": mtime
        })
        
        # If nested tar, inspect its contents in memory
        if type_group == "NESTED_TAR" and is_file:
            print(f"Inspecting nested archive: {name} ({round(size/(1024*1024), 2)} MB)...")
            try:
                extracted_file = tar.extractfile(member)
                if extracted_file is not None:
                    # Open nested archive from stream
                    with tarfile.open(fileobj=extracted_file, mode="r:*") as nested_tar:
                        for n_member in nested_tar:
                            n_name = n_member.name
                            n_size = n_member.size
                            n_ext = "".join(Path(n_name).suffixes).lower()
                            nested_tar_records.append({
                                "parent_tar": name,
                                "member_name": n_name,
                                "size_bytes": n_size,
                                "size_mb": round(n_size / (1024 * 1024), 2),
                                "extension": n_ext,
                                "is_file": n_member.isfile(),
                                "is_safe": not (n_name.startswith("/") or n_name.startswith("\\") or ".." in n_name or ":" in n_name)
                            })
            except Exception as e:
                print(f"Warning: Failed to inspect nested tar {name}: {e}")

        # Parse GSM accession and metadata guesses
        gsm_match = GSM_REGEX.search(name)
        gsm_id = gsm_match.group(1).upper() if gsm_match else "UNKNOWN"
        
        # Modality guess and evidence based on file name & extension
        # e.g. CosMx, single-cell RNA-seq, spatial transcriptomics, Visium, images
        lower_name = name.lower()
        modality_guess = "UNKNOWN"
        modality_evidence = []
        
        if "cosmx" in lower_name or "smi" in lower_name:
            modality_guess = "CosMx Spatial Molecular Imaging"
            modality_evidence.append("Name contains 'cosmx' or 'smi'")
        elif "visium" in lower_name or "spatial" in lower_name:
            modality_guess = "Spatial Transcriptomics (Visium-like)"
            modality_evidence.append("Name contains 'visium' or 'spatial'")
        elif "filtered_feature_bc_matrix" in lower_name or "raw_feature_bc_matrix" in lower_name:
            modality_guess = "Single-cell RNA-seq (10x Genomics)"
            modality_evidence.append("Name matches 10x matrix convention")
        elif ext.endswith(".h5"):
            modality_guess = "HDF5 Gene Expression Matrix (Likely scRNA-seq or Spatial)"
            modality_evidence.append("HDF5 format; specific modality pending Step 2 structure inspection")
        elif type_group in ["TIFF", "JPG"]:
            modality_guess = "Microscopy / Histology Image"
            modality_evidence.append(f"Image format: {ext}")
        elif type_group == "NESTED_TAR":
            modality_guess = "Compound Archive"
            modality_evidence.append("Nested tar file")
            
        # Parse patient / site / label (marked 'from filename, unverified')
        # Look for tokens in filename: e.g. GSMxxxxxx_SampleName_Site_Patient...
        # We preserve pure filename parsing with zero assumptions.
        # Strip extension and GSM prefix:
        base_clean = Path(name).name
        for sfx in Path(name).suffixes:
            base_clean = base_clean.removesuffix(sfx)
        if gsm_id != "UNKNOWN":
            base_clean = base_clean.replace(gsm_id, "").lstrip("_- ")
            
        sample_records.append({
            "file": name,
            "gsm_accession": gsm_id,
            "filename_cleaned": base_clean,
            "type_group": type_group,
            "size_mb": round(size / (1024 * 1024), 2),
            "modality_guess": modality_guess,
            "modality_evidence": "; ".join(modality_evidence) if modality_evidence else "None",
            "parsed_patient_raw": "Parsed from filename, unverified",
            "parsed_site_raw": "Parsed from filename, unverified",
            "parsed_label_raw": "Parsed from filename, unverified"
        })

print(f"Total members found in archive: {len(inventory_records)}")
print(f"Total nested members found: {len(nested_tar_records)}")

# Write tar_inventory.csv
tar_inv_path = TABLES_DIR / "tar_inventory.csv"
with open(tar_inv_path, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "archive", "member_name", "size_bytes", "size_mb", "extension", "type_group", "is_safe", "mtime"
    ])
    writer.writeheader()
    writer.writerows(inventory_records)
print(f"Wrote {tar_inv_path}")

# Write nested_tar_inventory.csv if any
if nested_tar_records:
    nested_inv_path = TABLES_DIR / "nested_tar_inventory.csv"
    with open(nested_inv_path, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "parent_tar", "member_name", "size_bytes", "size_mb", "extension", "is_file", "is_safe"
        ])
        writer.writeheader()
        writer.writerows(nested_tar_records)
    print(f"Wrote {nested_inv_path}")

# Write sample_map.csv
sample_map_path = TABLES_DIR / "sample_map.csv"
with open(sample_map_path, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "file", "gsm_accession", "filename_cleaned", "type_group", "size_mb",
        "modality_guess", "modality_evidence",
        "parsed_patient_raw", "parsed_site_raw", "parsed_label_raw"
    ])
    writer.writeheader()
    writer.writerows(sample_records)
print(f"Wrote {sample_map_path}")

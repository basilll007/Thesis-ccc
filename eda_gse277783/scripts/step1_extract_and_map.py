"""
Step 1: Extract GSE277783_RAW.tar and build refined inventory and sample map.
- Validate path safety for every member (no absolute paths, no .. traversal).
- Extract tar into F:\\Thesis\\eda_gse277783\\extracted\\ in streaming mode.
- Build tables/tar_inventory.csv (name, size, extension, type_group).
- Build tables/sample_map.csv (file, GSM, parsed patient/site/label, modality guess, evidence).
- Log decisions, findings, and modality assessment.
"""

import os
import re
import csv
import sys
import tarfile
from pathlib import Path
from datetime import datetime

TAR_PATH = Path(r"F:\Thesis\data\GSE277783_RAW.tar")
OUT_ROOT = Path(r"F:\Thesis\eda_gse277783")
EXTRACTED_DIR = OUT_ROOT / "extracted"
TABLES_DIR = OUT_ROOT / "tables"
LOGS_DIR = OUT_ROOT / "logs"

EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

GSM_PATTERN = re.compile(r"(GSM\d+)")

print(f"[{datetime.now().isoformat()}] Step 1 starting: Streaming extraction & mapping...")

tar_inventory = []
gsm_to_files = {}

# 1. Read tar and extract in streaming mode with path safety checks
with tarfile.open(TAR_PATH, mode="r:") as tar:
    members = tar.getmembers()
    print(f"Total archive members: {len(members)}")
    
    for m in members:
        # Path safety validation
        target_path = (EXTRACTED_DIR / m.name).resolve()
        if not str(target_path).startswith(str(EXTRACTED_DIR.resolve())):
            raise ValueError(f"SECURITY ALERT: Path traversal attempt in member: {m.name}")
        if m.name.startswith("/") or m.name.startswith("\\") or ".." in m.name:
            raise ValueError(f"SECURITY ALERT: Insecure member name: {m.name}")

        # Extension and grouping
        p = Path(m.name)
        suffixes = p.suffixes
        ext = "".join(suffixes).lower()
        
        type_group = "OTHER"
        if ext.endswith(".tar.gz") or ext.endswith(".tar"):
            type_group = "NESTED_TAR"
        elif ext.endswith(".h5") or ext.endswith(".h5ad"):
            type_group = "H5"
        elif ext.endswith(".tif") or ext.endswith(".tiff") or ext.endswith(".qp.tiff.gz"):
            type_group = "TIFF"
        elif ext.endswith(".jpg") or ext.endswith(".jpeg") or ext.endswith(".jpg.gz"):
            type_group = "JPG"
        elif ext.endswith(".csv") or ext.endswith(".tsv"):
            type_group = "TABLE"

        tar_inventory.append({
            "name": m.name,
            "size": m.size,
            "size_mb": round(m.size / (1024 * 1024), 2),
            "extension": ext,
            "type_group": type_group
        })
        
        # Associate with GSM
        gsm_match = GSM_PATTERN.search(m.name)
        gsm = gsm_match.group(1).upper() if gsm_match else "NO_GSM"
        if gsm not in gsm_to_files:
            gsm_to_files[gsm] = []
        gsm_to_files[gsm].append(m)

    # Now extract files in streaming mode if not already extracted
    print(f"Extracting {len(members)} files to {EXTRACTED_DIR}...")
    for idx, m in enumerate(members, 1):
        target_file = EXTRACTED_DIR / m.name
        if target_file.exists() and target_file.stat().st_size == m.size:
            continue  # Already extracted identically
        tar.extract(m, path=EXTRACTED_DIR, filter="data" if hasattr(tarfile, "data_filter") else None)
        if idx % 20 == 0 or idx == len(members):
            print(f"  Extracted {idx}/{len(members)} files...")

print(f"Extraction completed safely.")

# 2. Write tables/tar_inventory.csv
inv_csv = TABLES_DIR / "tar_inventory.csv"
with open(inv_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["name", "size", "size_mb", "extension", "type_group"])
    writer.writeheader()
    writer.writerows(tar_inventory)
print(f"Saved {inv_csv}")

# 3. Build refined sample map
# Expected GSM ranges:
# GSM8452847 - GSM8452895 (49 samples)
# GSM8529618 - GSM8529624 (7 samples) -> Check if present
# GSM8714177 - GSM8714184 (8 samples)
sample_map_rows = []

for gsm in sorted(gsm_to_files.keys()):
    files = [m.name for m in gsm_to_files[gsm]]
    h5_file = next((f for f in files if f.endswith(".h5")), "UNKNOWN")
    spatial_tar = next((f for f in files if "spatial.tar" in f), "UNKNOWN")
    img_file = next((f for f in files if (".jpg" in f or ".tiff" in f or ".tif" in f)), "UNKNOWN")
    
    # Parse patient, site, label from filename
    # E.g. GSM8452847_Pt-1A_filtered_feature_bc_matrix.h5 -> Pt-1A
    # E.g. GSM8714177_Pt3_Pri2_filtered_feature_bc_matrix.h5 -> Pt3_Pri2
    # E.g. GSM8714183_PDXPt10_filtered_feature_bc_matrix.h5 -> PDXPt10
    parsed_patient = "UNKNOWN"
    parsed_site = "UNKNOWN (not in filename)"
    parsed_label = "UNKNOWN"
    
    # Extract representative sample name by removing GSM prefix and standard suffixes
    rep_name = h5_file if h5_file != "UNKNOWN" else files[0]
    clean_name = re.sub(r"^GSM\d+_", "", rep_name)
    clean_name = clean_name.replace("_filtered_feature_bc_matrix.h5", "")
    clean_name = clean_name.replace("_spatial.tar.gz", "")
    clean_name = re.sub(r"\.(jpg|jpeg|tif|tiff)(\.gz)?$", "", clean_name)
    parsed_label = f"{clean_name} (from filename, unverified)"
    
    # Patient extraction
    if clean_name.startswith("Pt-") or clean_name.startswith("Pt"):
        m_pt = re.match(r"(Pt-?\d+)", clean_name, re.IGNORECASE)
        if m_pt:
            parsed_patient = f"{m_pt.group(1)} (from filename, unverified)"
    elif clean_name.startswith("PDXPt"):
        m_pt = re.match(r"PDXPt(\d+)", clean_name, re.IGNORECASE)
        if m_pt:
            parsed_patient = f"Pt{m_pt.group(1)} (PDX model) (from filename, unverified)"
            
    # Site / Condition extraction from filename
    if "Pri" in clean_name:
        parsed_site = "Primary tumor (parsed from 'Pri' in filename, unverified)"
    elif "PDX" in clean_name:
        parsed_site = "PDX xenograft (parsed from 'PDX' in filename, unverified)"
    else:
        # Letters A, B, C, D, E without site indication
        parsed_site = "UNKNOWN (site not indicated in filename; specimen labeled with letter suffix)"

    # Modality guess & evidence
    modality_guess = "10x Genomics Visium Spatial Transcriptomics"
    evidence_parts = []
    if h5_file != "UNKNOWN":
        evidence_parts.append(f"10x HDF5 expression matrix ({h5_file})")
    if spatial_tar != "UNKNOWN":
        evidence_parts.append(f"Space Ranger spatial folder ({spatial_tar})")
    if img_file != "UNKNOWN":
        evidence_parts.append(f"Histology image ({img_file})")
    modality_evidence = "; ".join(evidence_parts)

    sample_map_rows.append({
        "gsm_accession": gsm,
        "sample_label_raw": parsed_label,
        "parsed_patient": parsed_patient,
        "parsed_site": parsed_site,
        "data_modality": modality_guess,
        "modality_evidence": modality_evidence,
        "h5_file": h5_file,
        "spatial_tar": spatial_tar,
        "image_file": img_file,
        "total_files_for_sample": len(files)
    })

sample_map_csv = TABLES_DIR / "sample_map.csv"
with open(sample_map_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "gsm_accession", "sample_label_raw", "parsed_patient", "parsed_site",
        "data_modality", "modality_evidence", "h5_file", "spatial_tar", "image_file", "total_files_for_sample"
    ])
    writer.writeheader()
    writer.writerows(sample_map_rows)
print(f"Saved {sample_map_csv} with {len(sample_map_rows)} samples.")

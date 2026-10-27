# Exploratory Data Analysis Log: GEO GSE277783

- **Dataset**: GEO Series GSE277783 ("Spatial mapping of transcriptomic and lineage plasticity in metastatic pancreatic cancer")
- **Input File**: `F:\Thesis\data\GSE277783_RAW.tar`
- **Output Directory**: `F:\Thesis\eda_gse277783\`
- **Virtual Environment**: `F:\Thesis\.venv_eda\`
- **Random Seed**: 42 (fixed across all scripts)
- **Start Time**: 2026-10-27T16:55:00

---

## Initialization & Environment Setup
- **Timestamp**: 2026-10-27T16:55:00
- **Action**: Check system resources and setup isolated Python virtual environment.
- **System Findings**:
  - GPU: NVIDIA GeForce RTX 5060 Laptop GPU (8151 MiB total VRAM, ~7040 MiB free)
  - Memory: 32 GB RAM (~9.5 GB free physical memory)
  - Disk: Drive F has ~382 GB free space
  - Input tar file verified at `F:\Thesis\data\GSE277783_RAW.tar` (6,738,790,400 bytes, ~6.28 GB). Untouched and preserved read-only.
  - No GEO SOFT or series-matrix files found in `F:\Thesis\data\`.
- **Environment**: Created virtual environment at `F:\Thesis\.venv_eda\` using Python 3.12.13. Installed core dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `pillow`, `h5py`, `anndata`, `psutil`, `torch` (CUDA 12.4).
- **Rule Adherence**: Memory safety mode, streaming extraction with path bounds validation, lazy/backed H5 reads, no permanent filtering, strict provenance tracking for sample metadata.

---

## Step 1: Inventory & Modality Assessment (STOP Point for Judgment)
- **Timestamp**: 2026-10-27T17:01:29.758393
- **Action**: Comprehensive archive inventory, member path safety verification, streaming extraction, and modality classification.
- **Archive Contents**:
  - Total members: 171
  - Safe paths: 171 (100% relative, no directory traversal or absolute paths)
  - Nested archives: 57 *_spatial.tar.gz files containing 412 internal files (extracted safely to xtracted/spatial/<gsm>/).
  - Total GSM samples: 57
  - GSM Ranges:
    - GSM8452847 - GSM8452895 (49 samples: Pt-1A to Pt-13E)
    - GSM8714177 - GSM8714184 (8 samples: Pt3_Pri2, Pt6_Pri2, Pt6_Pri3, Pt8_Pri2, Pt13_Pri2, Pt13_Pri3, PDXPt10, PDXPt13)
    - **Missing Expected Range**: GSM8529618 - GSM8529624 is NOT present in GSE277783_RAW.tar.
- **Modality Discovery & Judgment**:
  - **Archive Modality**: 100% of samples (57/57) are **10x Genomics Visium Spatial Transcriptomics**.
  - Components per sample:
    1. Gene expression matrix: 10x HDF5 format (*_filtered_feature_bc_matrix.h5), 36,601 features (whole transcriptome).
    2. Spatial layout: Space Ranger nested archive (*_spatial.tar.gz) containing spot coordinates (	issue_positions.csv), scaling factors (scalefactors_json.json), and alignment images.
    3. Histology images: 55 .jpg.gz and 2 .qp.tiff.gz microscopy slide images.
  - **CosMx Modality Status**: No NanoString CosMx SMI files (e.g. cell_by_gene_counts.csv, FOV tables, negative-probe tables) are present inside GSE277783_RAW.tar.
- **Operational Plan**:
  1. Step 2: Lazy header inspection of 57 H5 files (shapes, features, barcodes); PIL dimension extraction and thumbnail generation (max 1000px) for 57 histology images into igures/thumbs/.
  2. Step 3: High-throughput single-sample QC across all 57 Visium slides (spot counts, UMI depth, gene counts, %MT, %zero), logging peak memory per file; generate 	ables/qc_per_sample.csv and distribution plots in igures/.
  3. Step 4: Visium spatial EDA (spot counts under tissue, coordinate ranges, spatial scatter plots colored by UMI depth); evaluate detection of the 30 user-specified key genes against whole-transcriptome background (	ables/key_gene_detection.csv).
  4. Step 5: Metadata structuring (patient x site x replicate), confound analysis, and patient-level sample size enumeration.
  5. Step 6: Comprehensive EDA_REPORT.md synthesis referencing all generated tables.

## Step 2: Per-file Integrity & Structure Results
- **Timestamp**: 2026-10-27T17:04:32.065645
- **Action**: Lazy inspection of H5 structures, image dimension parsing, thumbnail generation, and integrity validation.
- **Key Findings**:
  - **H5 Integrity**: 57/57 files verified intact (0 corrupt, 0 empty).
  - **H5 Format**: All 57 are 10x Genomics HDF5 format (matrix group containing CSC sparse data, barcodes, features).
  - **Feature Dimensions**:
    - 55 files have 36,601 features (human GRCh38 whole transcriptome).
    - 2 files (GSM8714183_PDXPt10 and GSM8714184_PDXPt13) have 18,085 features (human GRCh38 curated/protein-coding subset for xenograft tumor mapping).
  - **Coordinates & Metadata inside H5**: Spatial coordinates and cell-level annotations are absent from inside the H5 files; coordinates are maintained in the companion spatial/tissue_positions.csv files.
  - **CosMx Tables**: 0 CosMx expression or metadata CSV tables found in this archive.
  - **Image Integrity**: 57/57 images intact (0 corrupt).
    - 55 gzipped JPEG whole-slide histology scans (dimensions ~22k x 22k pixels).
    - 2 BigTIFF multi-resolution whole-slide scans (GSM8714183_Pt-10_PDX.qp.tiff.gz 69120x87840 px, and GSM8714184_Pt-13_PDX.qp.tiff.gz 78720x115200 px).
    - 57 thumbnails (max 1000px) safely generated and saved to igures/thumbs/.
  - **Peak Memory**: Monitored per file; peak process memory remained <120 MB throughout due to streaming and JPEG draft IDCT downsampling.

## Step 3: QC & Basic EDA Results
- **Timestamp**: 2026-10-27T17:07:15.379748
- **Action**: Per-sample spot-level quality control across all 57 Visium slides without permanent filtering.
- **Key Metrics & Distributions**:
  - Total spots in dataset: 149,088 spots.
  - Spots per slide: Min = 419 (GSM8714181), Median = 2,696, Max = 4,339 (GSM8452895).
  - Sequencing depth per spot: Min median depth = 130 UMIs (GSM8452848), Median = 3,895 UMIs, Max median depth = 14,841 UMIs (GSM8452857) -> **114.2-fold sequencing depth disparity across the cohort**.
  - Genes per spot: Min median = 93 genes (GSM8452848), Median = 1,878 genes, Max median = 4,382 genes.
  - Mitochondrial reads (%MT): Generally low/healthy throughout (Cohort median = 1.83%, max sample median = 10.8%).
  - Matrix Sparsity: 80.36% to 99.73% zeros (Cohort median = 94.1%).
- **Identified Problems & Flagged Samples (6/57 flagged)**:
  1. GSM8452848_Pt-1B: Severe technical drop-out / degradation (median 130 UMIs, 93.36% spots fail QC <500 UMIs).
  2. GSM8452866_Pt-6D: Very low depth (median 748 UMIs, 29.62% spots fail QC).
  3. GSM8452878_Pt-9D: Low depth (median 735 UMIs, 33.44% spots fail QC).
  4. GSM8714177_Pt3_Pri2: Very small tissue section (444 spots).
  5. GSM8714180_Pt8_Pri2: Low depth (median 948 UMIs, 18.08% spots fail QC).
  6. GSM8714181_Pt13_Pri2: Both low spot count and low depth (419 spots, median 673 UMIs, 39.14% fail QC).
- **Barcode Whitelist Recurrence**:
  - 9,541 unique barcodes observed.
  - A single barcode is shared by up to 54 slides, which confirms standard 10x Visium array coordinate design (not cell droplet duplication).
- **Outputs Produced**: 	ables/qc_per_sample.csv, 	ables/barcode_sharing_summary.csv, 57 individual distribution plots in igures/qc/, and cohort overview plots igures/cohort_qc_summary.png and .pdf.

## Step 4: Spatial EDA & Key-Gene Detection Results
- **Timestamp**: 2026-10-27T17:09:22.561667
- **Action**: Spatial coordinate validation, slide-by-slide spatial plotting, and key gene detection analysis across 149,088 spots.
- **Key Spatial Findings**:
  - All 57 slides contain valid spatial coordinates from Space Ranger (	issue_positions.csv or 	issue_positions_list.csv).
  - Spots under tissue per slide range from 419 to 4,339 (median 2,696).
  - Individual spatial scatter plots generated in both PNG and PDF formats (57 pairs, 114 files in igures/spatial/).
  - Negative-control probes: 0 (N/A for 10x Visium whole-transcriptome sequencing).
  - Cell-type annotations: False (none present in raw Space Ranger outputs; per rule 4, no annotations imputed).
- **Key Gene Panel Detection (30 genes tested across 149,088 spots)**:
  - Robustly Detected (>5% spots, 22/30 genes):
    - Epithelial / Ductal: MIF (88.1%, 2.04M UMIs), KRT19 (79.1%, 2.12M UMIs), KRT7 (63.1%), EPCAM (49.9%).
    - Stromal / Desmoplastic: FN1 (74.2%, 1.39M UMIs), COL1A1 (71.9%, 1.09M UMIs), ITGB1 (68.2%), S100A4 (62.9%), MMP2 (37.6%), ACTA2 (34.3%), FAP (10.6%).
    - Immune / Signaling: CD74 (70.1%), CD44 (56.9%), VEGFA (55.5%), SPP1 (49.2%), TGFBR2 (38.6%), TGFB1 (37.9%), CXCR4 (24.4%), CXCL12 (16.8%), TGFBR1 (15.2%), CCL5 (11.2%), PTPRC (8.8%).
  - Low Expression (1-5% spots, 6/30 genes):
    - CD86 (4.19%), KDR (3.73%), CD274 (PD-L1, 3.23%), CD3E (3.01%), PDCD1 (PD-1, 2.72%), CD8A (2.65%).
  - Background Floor (<1% spots, 2/30 genes):
    - CCR5: 0.98% spots (1,541 total UMIs).
    - CTLA4: 0.53% spots (826 total UMIs).
  - Both CCR5 and CTLA4 are effectively at the technical background detection floor in these spatial transcriptomics profiles.

## Step 5: Metadata Structure & Confounds Results
- **Timestamp**: 2026-10-27T17:10:21.748493
- **Action**: Patient-level aggregation, confound quantification, and sample size enumeration.
- **Patient Cohort Summary**:
  - Exactly 13 unique human patients identified: Pt-1 through Pt-13.
  - Total Visium slides: 57 slides (3 to 8 slides per patient).
  - Explicit Primary in filename: 4 patients (Pt-3, Pt-6, Pt-8, Pt-13 via Pri2/Pri3).
  - Explicit Metastasis in filename: 0 patients (all metastatic and primary site identities in Pt-1 through Pt-13 are masked by letter suffixes A-E).
  - PDX models: 2 patients (Pt-10, Pt-13).
  - Maximum patient-level sample size for primary vs metastasis comparison: N <= 13 patients (NOT 149,088 spots).
- **Major Confounds Quantified**:
  1. Pseudoreplication (149k spots clustered across only 13 patients).
  2. Missing anatomic site metadata (unverified letter suffixes A-E).
  3. Sequencing depth variation (114-fold disparity across slides).
  4. Slide architecture batch effects (standard Visium v1 vs CytAssist v2).
  5. PDX xenograft reference disparity (18k vs 36k genes, mouse stroma).
  6. Patient sample imbalance (Pt-13 has 8 slides; Pt-1, Pt-2, Pt-5 have 3 slides).
  7. Modality limitation (CosMx missing from this raw archive).

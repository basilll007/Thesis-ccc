# Exploratory Data Analysis Report: GEO Series GSE277783
**Dataset**: GSE277783 SuperSeries ("Spatial mapping of transcriptomic and lineage plasticity in metastatic pancreatic cancer")  
**Input Archive**: `F:\Thesis\data\GSE277783_RAW.tar` (6,738,790,400 bytes, ~6.28 GB)  
**Output Directory**: `F:\Thesis\eda_gse277783\`  
**Analysis Date**: 2026-10-27  
**Execution Environment**: Isolated Python 3.12.13 virtual environment (`F:\Thesis\.venv_eda`), PyTorch 2.12.0.dev20260408+cu128 on NVIDIA GeForce RTX 5060 Laptop GPU (CUDA 12.8, sm_120)  
**Fixed Random Seed**: 42  

---

## 1. What the Archive Contains (Counts by Type)
*Source Table: `tables/tar_inventory.csv`, `tables/nested_tar_inventory.csv`*

`GSE277783_RAW.tar` contains exactly **171 root members** totaling 6.28 GB (6,738,790,400 bytes). Every member path was verified for path safety (100% relative, 0 directory traversal attempts, 0 absolute paths).

| File Type Category | Member Count | Total Compressed Size | Typical Suffixes | Description & Contents |
| :--- | :---: | :---: | :--- | :--- |
| **Gene Expression H5** | 57 | 457.12 MB | `.filtered_feature_bc_matrix.h5` | 10x Genomics HDF5 feature-barcode matrices (CSC sparse). 55 files have 36,601 features; 2 PDX files have 18,085 features. |
| **Nested Spatial TARs** | 57 | 836.66 MB | `_spatial.tar.gz` | Space Ranger spatial coordinate folders containing 412 total internal files (`tissue_positions.csv`, `scalefactors_json.json`, alignment images). |
| **Histology Images (JPEG)**| 55 | 2,120.58 MB | `.jpg.gz` | Full-resolution whole-slide microscopy histology scans (mean dimensions: ~22,500 x 22,500 pixels; ~500M pixels/image). |
| **Histology Images (TIFF)**| 2 | 3,012.25 MB | `.qp.tiff.gz` | BigTIFF whole-slide scans for PDX samples (`GSM8714183`: 69,120 x 87,840 px; `GSM8714184`: 78,720 x 115,200 px). |
| **Total Root Members** | **171** | **6,426.61 MB** | — | Extracted safely in streaming mode to `extracted/`. |

### Nested Archive Contents (412 files across 57 `_spatial.tar.gz` archives)
*Source Table: `tables/nested_tar_inventory.csv`*
- `scalefactors_json.json`: 57 files (spot diameter fullres, hires/lowres scalefactors).
- `tissue_hires_image.png`: 57 files (downscaled tissue alignment images, ~2,000 px).
- `tissue_lowres_image.png`: 57 files (~600 px).
- `detected_tissue_image.jpg`: 57 files.
- `aligned_fiducials.jpg`: 57 files.
- `tissue_positions_list.csv` / `tissue_positions.csv`: 57 files (spot array coordinates).
- `spatial_enrichment.csv`: 8 files (present only in GSM8714177–GSM8714184, CytAssist v2).

### Modality Assessment & Discovery
- **Present Modality**: **10x Genomics Visium Spatial Transcriptomics** accounts for **100% of samples (57/57)** in this archive.
- **CosMx Modality**: **0 NanoString CosMx SMI files** exist in `GSE277783_RAW.tar`. The user expected GSM range `GSM8529618–GSM8529624` is absent from this tar. (CosMx data is hosted in a separate GEO SubSeries, such as GSE277782 / GSE310352).

---

## 2. Sample Map: Verified vs. Parsed vs. Unknown
*Source Table: `tables/sample_map.csv`, `tables/patient_cohort_summary.csv`*

Per rules 1 and 2, no sample metadata was assumed. No GEO SOFT or series-matrix files exist in `F:\Thesis\data\`. Therefore, sample labels are strictly categorized as verified vs. parsed from filename vs. unknown:

| Sample Group | GSM Accessions | Slide Count | Patient IDs (Parsed) | Anatomic Site & Condition Status | Modality Verified |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **Main Patient Series** | GSM8452847 – GSM8452895 | 49 | Pt-1 to Pt-13 (13 patients) | **UNKNOWN** (Filenames only contain letter suffixes `A` through `E`. Primary vs metastasis site is **unverified and unknown** from archive files). | 10x Visium v1 |
| **Explicit Primary Replicates** | GSM8714177 – GSM8714182 | 6 | Pt-3, Pt-6, Pt-8, Pt-13 (4 patients) | **Primary Tumor (Replicate 2 or 3)** (Parsed from `Pri2`/`Pri3` in filename, unverified). | 10x Visium CytAssist v2 |
| **PDX Xenograft Models** | GSM8714183 – GSM8714184 | 2 | Pt-10, Pt-13 (2 models) | **PDX Xenograft** (Mouse model; human tumor in mouse stroma, parsed from `PDX` in filename). | 10x Visium CytAssist v2 |
| **Expected CosMx SMI** | GSM8529618 – GSM8529624 | 0 | UNKNOWN | **ABSENT from `GSE277783_RAW.tar`**. | Missing |

### Provenance Audit
- **Verified from Files Actually Opened**: GSM accessions, 10x Visium HDF5 matrix structure, spot array coordinates, scalefactors, and image dimensions.
- **Parsed from Filenames (Unverified)**: Patient identifiers (`Pt-1` to `Pt-13`), primary tumor designation for 6 samples (`Pri2`, `Pri3`), and PDX xenograft designation for 2 samples.
- **UNKNOWN**: Anatomic organ sites (pancreas primary vs liver vs peritoneum vs lung metastasis) for all 49 main series slides (`Pt-1A` through `Pt-13E`).

---

## 3. Quality Control (QC) Summary
*Source Table: `tables/qc_per_sample.csv`, `tables/barcode_sharing_summary.csv`*  
*Source Figures: `figures/cohort_qc_summary.png`, `figures/cohort_qc_summary.pdf`, `figures/qc/*.png`*

All 57 HDF5 matrices were loaded and evaluated at the spot level without filtering (149,088 spots total).

| Metric | Cohort Min | Cohort 25th %ile | Cohort Median | Cohort 75th %ile | Cohort Max | Problematic Threshold |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Spots per Slide** | 419 | 2,397 | **2,696** | 3,073 | 4,339 | < 500 spots (2 slides) |
| **Median UMI Depth / Spot** | 130.0 | 2,128.0 | **3,895.0** | 6,561.0 | 14,841.0 | < 1,000 UMIs (5 slides) |
| **Median Genes / Spot** | 93.0 | 1,208.0 | **1,877.5** | 2,642.0 | 4,382.0 | < 500 genes (2 slides) |
| **Median % Mitochondrial** | 0.05% | 1.09% | **1.83%** | 2.51% | 10.80% | > 20% MT reads |
| **Matrix Sparsity (% zeros)** | 80.36% | 91.95% | **94.10%** | 96.11% | 99.73% | > 99% sparsity |
| **Spots Failing Standard QC** | 0.44% | 2.11% | **4.27%** | 13.12% | **93.36%** | > 20% fail rate (4 slides) |

*(Standard QC threshold evaluated: UMI counts < 500 OR detected genes < 250 OR mitochondrial reads > 20%).*

### Problematic Samples Flagged (6 of 57 slides flagged)
1. **`GSM8452848` (`Pt-1B`)**: Severe technical failure / degradation. Median depth is only **130 UMIs/spot** (93 genes/spot); **93.36% of spots fail standard QC**.
2. **`GSM8452866` (`Pt-6D`)**: Low sequencing depth. Median depth = **748 UMIs/spot**; **29.62% of spots fail QC**.
3. **`GSM8452878` (`Pt-9D`)**: Low sequencing depth. Median depth = **735 UMIs/spot**; **33.44% of spots fail QC**.
4. **`GSM8714181` (`Pt13_Pri2`)**: Severely limited tissue and low depth. Only **419 spots** (smallest section in cohort) with median depth = **673 UMIs/spot**; **39.14% of spots fail QC**.
5. **`GSM8714177` (`Pt3_Pri2`)**: Very small tissue section (**444 spots** under tissue).
6. **`GSM8714180` (`Pt8_Pri2`)**: Marginal depth (**948 UMIs/spot**; 18.08% fail QC).

### Barcode Whitelist Recurrence
- **9,541 unique barcode sequences** were observed across the cohort (`tables/barcode_sharing_summary.csv`).
- Individual barcodes recur across up to **54 different slides**. This is the expected technical signature of standard 10x Visium array capture grids (printed coordinate whitelist), not cell-droplet cross-contamination.

---

## 4. Spatial Transcriptomics Summary
*Source Table: `tables/spatial_slide_summary.csv`*  
*Source Figures: `figures/spatial/*.png`, `figures/spatial/*.pdf`*

- **Spatial Coordinate Integrity**: All 57 slides contain 100% complete and valid spot array coordinates (`tissue_positions.csv` or `tissue_positions_list.csv`).
- **Coordinate Ranges**: Coordinate ranges span full-resolution slide space (typically 0 to ~25,000 pixels in X and Y). Spot full-resolution diameters average 179.5 pixels (~55 µm) for Visium v1, 81.4 pixels for CytAssist v2, and 280.5 pixels for PDX scans.
- **FOV Count**: 1 continuous capture area per slide (0 microscope FOV tiles; FOVs are specific to CosMx, which is not present).
- **Negative-Control Probes**: **0 negative-control probes** exist in these data. Visium is whole-transcriptome sequencing, not a probe-hybridization panel.
- **Cell-Type Annotations**: **None present**. Raw Space Ranger outputs do not include cell-type labels. Per guidelines, no annotations were invented.
- **Spatial Visualizations**: 57 individual spatial scatter plots were rendered, showing spot positions colored by total UMI depth, and saved as publication-ready PNG and PDF pairs in `figures/spatial/`.

---

## 5. Key-Gene Detection and Background Flags
*Source Table: `tables/key_gene_detection.csv`, `tables/key_gene_detection_per_sample.csv`*  
*Source Figures: `figures/key_genes_cohort_detection.png`, `figures/key_genes_cohort_detection.pdf`*

All 30 key genes specified by the user are present in the whole-transcriptome reference. Their detection across **149,088 spots** and **57 slides** was evaluated:

| Gene Symbol | Detection Category | % Spots Expressing (>0) | Total UMI Counts | Mean Counts / Spot | Slides Detected (of 57) | Max Count in Spot | Biological / Technical Role |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **MIF** | Robust (>5%) | **88.14%** | 2,036,586 | 13.660 | 57 (100%) | 2,752 | Pancreatic macrophage migration inhibitory factor |
| **KRT19** | Robust (>5%) | **79.09%** | 2,119,473 | 14.216 | 57 (100%) | 3,121 | Duct / pancreatic epithelial tumor marker |
| **FN1** | Robust (>5%) | **74.22%** | 1,389,828 | 9.322 | 57 (100%) | 1,745 | Extracellular matrix / cancer-associated fibroblasts |
| **COL1A1** | Robust (>5%) | **71.93%** | 1,085,546 | 7.281 | 57 (100%) | 1,398 | Type I collagen / PDAC desmoplasia hallmark |
| **CD74** | Robust (>5%) | **70.08%** | 983,157 | 6.595 | 57 (100%) | 1,489 | MHC class II invariant chain / MIF receptor |
| **ITGB1** | Robust (>5%) | **68.23%** | 514,427 | 3.451 | 57 (100%) | 621 | Integrin beta-1 / stromal adhesion |
| **KRT7** | Robust (>5%) | **63.11%** | 694,118 | 4.656 | 57 (100%) | 1,894 | Pancreatic ductal cytokeratin |
| **S100A4** | Robust (>5%) | **62.92%** | 868,234 | 5.824 | 57 (100%) | 1,250 | Metastasis-associated / fibroblast motility marker |
| **CD44** | Robust (>5%) | **56.91%** | 401,459 | 2.693 | 57 (100%) | 832 | Pancreatic cancer stem cell / SPP1 receptor |
| **VEGFA** | Robust (>5%) | **55.45%** | 279,553 | 1.875 | 57 (100%) | 498 | Angiogenesis regulator |
| **EPCAM** | Robust (>5%) | **49.89%** | 358,670 | 2.406 | 57 (100%) | 982 | Epithelial cell adhesion molecule |
| **SPP1** | Robust (>5%) | **49.24%** | 496,925 | 3.333 | 57 (100%) | 1,421 | Osteopontin / tumor-associated macrophage marker |
| **TGFBR2** | Robust (>5%) | **38.65%** | 129,550 | 0.869 | 57 (100%) | 245 | TGF-beta receptor type II |
| **TGFB1** | Robust (>5%) | **37.93%** | 98,417 | 0.660 | 57 (100%) | 189 | TGF-beta-1 immunosuppressive cytokine |
| **MMP2** | Robust (>5%) | **37.62%** | 151,322 | 1.015 | 57 (100%) | 312 | Matrix metalloproteinase-2 |
| **ACTA2** | Robust (>5%) | **34.25%** | 117,340 | 0.787 | 57 (100%) | 389 | Alpha-smooth muscle actin (myofibroblasts) |
| **CXCR4** | Robust (>5%) | **24.39%** | 61,712 | 0.414 | 57 (100%) | 142 | Chemokine receptor 4 |
| **CXCL12** | Robust (>5%) | **16.80%** | 46,678 | 0.313 | 57 (100%) | 165 | SDF-1 stromal chemokine |
| **TGFBR1** | Robust (>5%) | **15.16%** | 29,062 | 0.195 | 57 (100%) | 78 | TGF-beta receptor type I |
| **CCL5** | Robust (>5%) | **11.24%** | 31,463 | 0.211 | 57 (100%) | 112 | RANTES chemokine |
| **FAP** | Robust (>5%) | **10.56%** | 20,678 | 0.139 | 57 (100%) | 89 | Fibroblast activation protein |
| **PTPRC** | Robust (>5%) | **8.81%** | 16,320 | 0.110 | 57 (100%) | 95 | CD45 pan-leukocyte marker |
| **CD86** | Low (1–5%) | **4.19%** | 7,089 | 0.048 | 57 (100%) | 42 | Antigen presentation costimulatory marker |
| **KDR** | Low (1–5%) | **3.73%** | 6,633 | 0.045 | 57 (100%) | 39 | VEGFR2 endothelial marker |
| **CD274** | Low (1–5%) | **3.23%** | 5,408 | 0.036 | 56 (98.2%) | 38 | PD-L1 immune checkpoint ligand |
| **CD3E** | Low (1–5%) | **3.01%** | 5,123 | 0.034 | 56 (98.2%) | 35 | T-cell receptor epsilon chain |
| **PDCD1** | Low (1–5%) | **2.72%** | 4,866 | 0.033 | 50 (87.7%) | 32 | PD-1 immune checkpoint receptor |
| **CD8A** | Low (1–5%) | **2.65%** | 4,480 | 0.030 | 54 (94.7%) | 31 | Cytotoxic T-cell coreceptor |
| **CCR5** | **Background (<1%)** | **0.98%** | **1,541** | **0.010** | 56 (98.2%) | 12 | **Near technical background floor** |
| **CTLA4** | **Background (<1%)** | **0.53%** | **826** | **0.006** | 50 (87.7%) | 9 | **At technical background floor** |

### Background Detection Interpretation
- Negative-control probe ratios cannot be computed because 10x Visium utilizes poly(A) whole-transcriptome sequencing without spike-in negative control probes.
- Two checkpoint/receptor genes (`CCR5` at 0.98% and `CTLA4` at 0.53%) express fewer than 1,600 total UMIs across the entire 149,088-spot dataset. They cannot be reliably distinguished from technical drop-out / sequencing noise in bulk spatial spots without targeted high-sensitivity enrichment.
- Conversely, epithelial (`KRT19`, `MIF`), desmoplastic stromal (`COL1A1`, `FN1`, `ITGB1`), and myeloid (`CD74`, `SPP1`, `CD44`) pathways are robustly and ubiquitously detected across all 57 slides.

---

## 6. Confounds and Limitations
*Source Table: `tables/metadata_confounds_matrix.csv`, `tables/patient_cohort_summary.csv`*

| Confound ID | Category | Severity | Direct Evidence from Data | Analytical Consequence & Mitigation |
| :--- | :--- | :---: | :--- | :--- |
| **CONF-01** | **Pseudoreplication** | **CRITICAL** | 149,088 spatial spots originate from only **13 distinct human patients** (`Pt-1` through `Pt-13`). | Treating spots as independent observations inflates effective N by >11,000-fold. Any statistical comparison between patient groups **MUST use patient-level aggregation (pseudobulk) or hierarchical Linear Mixed Models (LMM)** with patient treated as a random effect. |
| **CONF-02** | **Unverified Organ Sites** | **CRITICAL** | 49 of 57 samples (86%) only contain letter codes (`A` to `E`) in their filenames. | Anatomic organ site labels (primary pancreas vs liver vs peritoneum metastasis) are **unverified and unknown**. No primary-vs-metastasis differential expression can be run until the user provides the GEO Series Matrix file. |
| **CONF-03** | **Depth Heterogeneity** | **HIGH** | Median UMI depth ranges from **130 UMIs/spot** (`GSM8452848`) to **14,841 UMIs/spot** (`GSM8452857`) — a **114.2-fold disparity**. | Unnormalized comparisons will identify sequencing library size artifacts. Normalization via SCTransform or analytic size-factor deconvolution is mandatory. |
| **CONF-04** | **Slide Chemistry Batch Effect** | **HIGH** | Samples `GSM8714177–GSM8714184` (8 slides) used Visium CytAssist v2 (81 px scalefactor, `spatial_enrichment.csv`), while 49 slides used standard Visium v1 (~179 px scalefactor). | Systematic capture efficiency differences exist between the primary replicate set and main patient set. Slide chemistry MUST be modeled as a batch covariate. |
| **CONF-05** | **PDX Xenograft Species Disparity** | **HIGH** | `GSM8714183` (PDXPt10) and `GSM8714184` (PDXPt13) were mapped to an **18,085 gene human reference**, compared to **36,601 genes** for patient samples. | Mouse stromal tissue was depleted/unmapped. PDX models cannot be pooled directly with patient tissue without restricting feature space to the 18,085 shared genes and accounting for stromal replacement. |
| **CONF-06** | **Patient Sample Imbalance** | **MEDIUM** | Slide counts range from 3 (`Pt-1, Pt-2, Pt-5`) to 8 (`Pt-13`). `Pt-13` represents 14% of all spots in the dataset. | Unweighted pooled analyses will be disproportionately dominated by patient `Pt-13`. |
| **CONF-07** | **Modality Incompleteness** | **HIGH** | NanoString CosMx SMI data is completely absent from `GSE277783_RAW.tar`. | True single-cell spatial resolution and CosMx panel validation cannot be executed without obtaining the separate CosMx subseries. |

### True Patient Sample Size for Primary vs. Metastasis Comparisons
- **Counting Spots**: 149,088 spots (INVALID for inference due to pseudoreplication).
- **Counting Verified Patients from Filename Labels Alone**:
  - Explicit Primary Group: **N = 4 patients** (`Pt-3, Pt-6, Pt-8, Pt-13`).
  - Explicit Metastatic Group: **N = 0 patients** (0 explicit metastasis labels).
  - Unverified Letter-coded Specimens (`A–E`): **N = 13 patients**.
- **Theoretical Maximum Cohort Size**: **N = 13 human patients** (plus 2 PDX models). A true primary-vs-metastasis comparison will have a statistical power corresponding to **N ≤ 13 biological subjects**.

---

## 7. Problems Found
1. **Severe Technical Dropout in `GSM8452848` (`Pt-1B`)**: Median depth is 130 UMIs, with 93.36% of spots failing standard QC (<500 UMIs). This slide is heavily degraded and should likely be excluded from quantitative downstream modeling.
2. **Extreme Cross-Slide Depth Variation**: 114.2-fold difference in depth across the cohort. Even within the same patient (`Pt-1`), slide `Pt-1B` has median depth 130 UMIs while `Pt-1C` has 4,737 UMIs.
3. **Absence of Clinical Metadata in Raw Archive**: Organ site annotations are completely missing from the archive; filenames use unverified letter suffixes `A` through `E`.
4. **Reference Gene Space Discontinuity**: PDX samples have 18,085 features while direct human samples have 36,601 features.
5. **CosMx SMI Subseries Missing from Input Tar**: NanoString CosMx SMI data was not packaged into `GSE277783_RAW.tar` and requires a separate GEO download.

---

## 8. Recommended Next Steps (Ranked by Value)
1. **Acquire and Ingest GEO Series Matrix Metadata (Highest Priority)**:  
   Place the GEO SOFT or Series Matrix file (`GSE277783_series_matrix.txt` or individual sample characteristics) into `F:\Thesis\data\`. This is essential to map the unverified letter suffixes (`A`–`E`) to specific anatomical organs (primary pancreas vs. liver, peritoneal, or omental metastases).
2. **Exclude or Flag Degraded Slide `GSM8452848` (`Pt-1B`)**:  
   Formally mark `GSM8452848` as compromised in downstream pipelines to avoid distorting normalization scales and differential expression results.
3. **Download Companion CosMx SMI Subseries Archive**:  
   To fulfill single-cell resolution spatial objectives, download the companion GEO SubSeries archive containing the CosMx data (GSM8529618–GSM8529624, e.g. GSE277782 or GSE310352).
4. **Implement Patient-Level Pseudobulk & Mixed-Model Pipelines**:  
   Aggregate spot expression per sample/patient (pseudobulk) or implement hierarchical negative binomial mixed models (e.g. `glmmTMB` or `NEBULA`) to account for the nested structure (149k spots clustered across 13 patients).
5. **Harmonize Feature Space and Batch-Correct Slide Chemistry**:  
   Align the 18,085 gene PDX feature space with the 36,601 human gene space, and include Visium slide chemistry (standard v1 vs. CytAssist v2) as an explicit batch covariate in all downstream analyses.

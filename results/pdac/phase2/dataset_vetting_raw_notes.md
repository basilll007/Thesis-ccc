# PDAC single-cell spatial dataset vetting: raw notes

Session date: 2026-10-24. Every value below comes from a page or endpoint fetched in this session. Anything not confirmed is marked UNVERIFIED. Fetch failures are marked BLOCKED, with the error text.

Method note on gene panels: no large file was downloaded. Gene lists were read with HTTP Range requests that fetched only small byte ranges, under a hard 5 MB cap per run:
- zip central directory plus the first ~200 KB of each deflated `exprMat_file.csv` member (Mendeley, ~3.7 MB total);
- the first 200 KB of gzipped GEO CSVs;
- HDF5 metadata plus the `var` index of `.h5ad` files, via h5py on a range-reading file object (~1.1-1.2 MB per file).

Saved panel/header files are under `F:\Thesis\data\pdac_panels\<candidate>\`.

---

## 1. Mendeley Data "PDAC-SMI and co-culture dataset" (kx6b69n3cb)

URLs fetched
- https://data.mendeley.com/public-api/datasets/kx6b69n3cb : HTTP 200 (JSON)
- https://data.mendeley.com/public-api/datasets/kx6b69n3cb/files?folder_id=fbd2fed3-...&version=1 : HTTP 200
- https://data.mendeley.com/public-api/datasets/kx6b69n3cb/folders/1 : HTTP 200 (folders "co-culture", "SMI")
- https://data.mendeley.com/datasets/kx6b69n3cb (WebFetch) : 200. The page showed the title, description and "CC BY 4.0". The Files section was empty in the rendered capture.
- HEAD https://data.mendeley.com/public-files/datasets/kx6b69n3cb/files/ba959853-2391-459a-afca-73ad8db8b63c/file_downloaded : 302 to S3, then 200, `Content-Disposition: attachment; filename="SMI_raw.zip"`, `Accept-Ranges: bytes`, `Content-Length: 5519589653`
- Range reads of SMI_raw.zip (central directory plus member headers) and raw_meta_data_final.h5ad (HDF5 metadata and `var`)
- Paper: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10327107/fullTextXML : 200 (bioRxiv preprint "Therapy-associated remodeling of pancreatic cancer revealed by single-cell spatial transcriptomics and optimal transport analysis")
- Zenodo: https://zenodo.org/api/records/7963531 : 200

Verbatim evidence
- API: `"name":"PDAC-SMI and co-culture dataset"`; `"The spatial molecular imaging (SMI) data for 13 PDAC tumors that either received neoadjuvant chemo/radiation therapy or were treatment naive. Raw data files (images, cell segmentation, expression data, metadata) are included in the SMI_raw folder. The count matrix, metadata, and cell type annotations are included in the PDAC_raw_meta_data.h5ad file."`
- API licence: `"short_name":"CC BY 4.0","full_name":"Creative Commons Attribution 4.0 International"`, url http://creativecommons.org/licenses/by/4.0
- API dataset size: `"size":6405533496`; published 2024-07-11; DOI 10.17632/kx6b69n3cb.1
- Files listed by the API:
  - `Co-culture_seurat.rds`: 585,225,537 B
  - `raw_meta_data_final.h5ad`: 300,718,306 B. Download: https://data.mendeley.com/public-files/datasets/kx6b69n3cb/files/26999a73-4b05-4f47-9cab-023cd0de80a9/file_downloaded
  - `SMI_raw.zip`: 5,519,589,653 B. Download: https://data.mendeley.com/public-files/datasets/kx6b69n3cb/files/ba959853-2391-459a-afca-73ad8db8b63c/file_downloaded
- Paper: "dissected the spatial transcriptomic landscape of primary resected human PDAC tumors with ( n = 6) or without ( n = 7) neoadjuvant chemotherapy and radiotherapy (CRT) at single-cell resolution using spatial molecular imaging (SMI) ... and a 960-plex gene panel, with ( n = 6) and without ( n = 9) 30 custom probes"
- Paper: "across 320 FOVs ... This yielded 753,133 high confidence single-cell profiles."
- Paper: "U, untreated; T, treated; b, base 960-plex panel; a, augmented 990-plex panel."
- Paper data availability: "The raw and processed SMI and co-culture snRNA-seq data are available at Mendeley Data ( https://doi.org/10.17632/kx6b69n3cb.1 ) and Zenodo ( https://doi.org/10.5281/zenodo.7963531 )."
- Paper: the authors ran ligand-receptor CCC analysis (SCOTIA). Supplementary Table S3 is "Ligand receptor pairs used for SCOTIA analysis" (media-4.xlsx, not fetched).
- Zenodo 7963531: "PDAC-SMI transcript locations", 16 `*_tx_file.csv.zip` files totalling 9,041,394,811 B, `access_right: open`. No licence field in the API response, so the Zenodo licence is UNVERIFIED.

What the zip central directory shows (saved as `SMI_raw_zip_central_directory_listing.tsv`, 3,167 entries)
- Contents: 16 `SMI_images/<sample>/<sample>_exprMat_file.csv` files, 185,107,582 B compressed (1,916,247,861 B uncompressed). Also 15 `_metadata_file.csv` files (36,585,022 B compressed) and 16 `_fov_positions_file.csv` files, plus image folders (CellLabels, CompartmentLabels, CellOverlay, CellComposite).
- Missing file: U4-a has no `_metadata_file.csv` in the zip.
- Sample IDs: U1-b, U2-b, U3-b, U4-a, U4-b, U5-a, U6-a_1, U6-a_2, U7-a, T1-b, T2-b, T3-b, T4-a, T4-b, T5-a, T6-b.
- Metadata header (U1-b): `"fov","cell_ID","Area","AspectRatio","CenterX_local_px","CenterY_local_px","CenterX_global_px","CenterY_global_px",...`. Single-cell centroids are present, in pixels. The µm-per-pixel factor is not stated, so it is UNVERIFIED.
- exprMat headers:
  - "-b" files: 981 columns = fov, cell_ID, 960 genes, 19 NegPrb.
  - "-a" files: 1011 columns = fov, cell_ID, 990 genes, 19 NegPrb.
  - The 30 genes in the augmented panel only: ANXA10, BTNL8, CDH17, CELA3B, CFTR, CNTN4, CPA1, CPA2, CTNND2, DCDC2, FLRT2, GALR1, GATA6, GCK, GLIS3, HNF1B, KCNJ16, LGALS4, NFIB, NKX2-2, NR5A2, PLA2G10, RELN, RLN1, SEMA3E, SEMA5A, SYTL2, TFF1, TP63, TSPAN8.
- h5ad structure, read by range: top-level keys `X, obs, var`. `var` has 1009 entries (index `gene`, column `probe_type`). `obs` columns include `CenterX_global_px, CenterY_global_px, annotation_majortypes, annotation_subtypes, donor, fov, sample_id, treatment_status`.
- Whether `X` holds raw counts is UNVERIFIED.
- Gene check: 37 of 39 genes are present in both panels and in the h5ad var list. PTN and SDC4 are absent.

Saved files
- `F:\Thesis\data\pdac_panels\mendeley_pdac_smi\<sample>_exprMat_header.txt` (16 files)
- `raw_meta_data_final_h5ad_var_names.txt`
- `SMI_raw_zip_central_directory_listing.tsv`

---

## 2. GEO GSE267680

URLs fetched
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267680 (HTML): BLOCKED. curl got a Google reCAPTCHA challenge page instead of the GEO page (`<base href="https://www.google.com/recaptcha/challengepage/">`).
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267680&targ=self&form=text&view=brief : 200 (SOFT text)
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267680/suppl/filelist.txt : 200
- GSM8273019 SOFT text: 200. GPL33762 SOFT text: 200.
- Range read of the first 1.5 MB of https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8273nnn/GSM8273021/suppl/GSM8273021_PanIN1132.tar.gz (tar member order only)

Verbatim evidence
- "!Series_title = PanIN and CAF Transitions in Pancreatic Carcinogenesis Revealed with Spatial Data Integration [Xenium]"
- "!Series_overall_design = 5 tissue sections from surgical sections of resected primary pancreatic ductal adenocarcinoma from 3 unique human subjects, which are Xenium data inputs to recreate results of validation experiments conducted as an extension of analyses conducted with GSE254829."
- "!Series_platform_id = GPL33762". The platform record reads "!Platform_title = Xenium In Situ Analyzer: Homo sapiens".
- Samples: GSM8273019 to GSM8273023. "!Sample_source_name_ch1 = Pancreatic Intraductal Neoplasm"; "subject: PanIN-LG1".
- Panel: "!Sample_extract_protocol_ch1 = The capture probe panel was comprised of the 10x Genomic 280 gene breast cancer panel with the addition of 100 additional custom probes."
- Per-sample contents listed in GEO: "cells.csv.gz: compressed cell locations", "cell_feature_matrix.h5", "gene_panel.json: gene panel list", "cell_feature_matrix/features.tsv.gz".
- filelist.txt:
  - GSE267680_RAW.tar: 29,470,064,640 B
  - GSM8273019_PanIN1131_S1A.tar.gz: 4,742,128,544 B
  - GSM8273020_PanIN1131_S1C.tar.gz: 4,230,082,018 B
  - GSM8273021_PanIN1132.tar.gz: 5,554,325,014 B
  - GSM8273022_PanIN1134.tar.gz: 6,289,287,592 B
  - GSM8273023_PanIN1144.tar.gz: 8,654,235,621 B
- The tar.gz peek showed `PanIN1132/PanIN1132_de_gene_groups.csv` (1,493 B) first, then `analysis.h5` (14,301,892 B). `gene_panel.json` was not reachable within the 5 MB cap.
- Gene presence: all UNVERIFIED.
- Licence: not stated on fetched pages. The GEO disclaimer page was BLOCKED (see section 6).

---

## 3. bioRxiv 10.1101/2025.02.11.637627

URLs fetched
- https://api.biorxiv.org/details/biorxiv/10.1101/2025.02.11.637627 : 200
- https://www.biorxiv.org/content/10.1101/2025.02.11.637627v1 (curl): 200. The page contained the abstract only.
- https://www.biorxiv.org/content/early/2025/02/15/2025.02.11.637627.source.xml : BLOCKED, "HTTP 429" (Too Many Requests)
- https://www.biorxiv.org/content/10.1101/2025.02.11.637627v1.full (WebFetch, and a later curl retry): BLOCKED, "The server returned HTTP 429 Too Many Requests. Retry-After: 114"
- https://link.springer.com/article/10.1007/s13402-025-01100-6 : BLOCKED, 303 redirect to idp.springer.com/authorize (login/cookie). The redirect was not followed.
- Europe PMC search plus https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12528271/fullTextXML : 200 (published version)

Verbatim evidence
- API: "title":"Integrating multi-modal transcriptomics identifies cellular subtypes with distinct roles in PDAC progression", "license":"cc_no", "published":"10.1007\/s13402-025-01100-6"
- Abstract: "... and 2 Xenium high-resolution ST-seq samples with 307,679 cells"
- PMC12528271: "One sample was derived from the FFPE Human PDAC dataset available in the 10x Genomics database, specifically the Human Immuno-Oncology Profiling Panel ( https://www.10xgenomics.com/datasets/ffpe-human-ductal-adenocarcinoma-data-with-human-immuno-oncology-profiling-panel-1-standard ). The second sample was obtained from the GEO database, consisting of a PDAC surgical specimen associated with high-grade pancreatic intraepithelial neoplasia (PanIN) ( GSE267680 , GSM8273021 )"
- Data availability: "The data and accession numbers used in this study can be found in the article and supplementary tables. The codes used for bioinformatics analysis are available at https://github.com/JunekureWu/multi-modal_PDAC ."
- Article licence: "licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License"

Conclusion: no new single-cell spatial data. The paper reuses candidate 4a (10x PDAC IO) and GSE267680 GSM8273021.

---

## 4. 10x Genomics public datasets

URLs fetched, all BLOCKED
- https://www.10xgenomics.com/datasets/ffpe-human-ductal-adenocarcinoma-data-with-human-immuno-oncology-profiling-panel-1-standard : BLOCKED. curl got HTTP 429 with body "Vercel Security Checkpoint". WebFetch got "The server returned HTTP 429 Too Many Requests."
- /jp/ variant of the same page (WebFetch): BLOCKED, HTTP 429
- .../pancreatic-cancer-with-xenium-human-multi-tissue-and-cancer-panel-1-standard : BLOCKED, HTTP 429 (Vercel Security Checkpoint)
- .../ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard : BLOCKED, HTTP 429
- .../human-pancreas-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard : BLOCKED, HTTP 429
- https://www.10xgenomics.com/support/software/xenium-explorer/latest/resources/xenium-explorer-demos : BLOCKED, HTTP 429
- Two retries later in the session also returned 429.
- No cf.10xgenomics.com download URL was guessed or probed. The 10x bundle sizes, download URLs and licences are therefore all UNVERIFIED.

Search results (WebSearch summaries; the pages themselves were not fetched, so this is not evidence)
- Dataset pages exist for:
  - "FFPE Human Pancreatic Ductal Adenocarcinoma Data with Human Immuno-Oncology Profiling Panel"
  - "Pancreatic Cancer with Xenium Human Multi-Tissue and Cancer Panel"
  - "FFPE Human Pancreas with Xenium Multimodal Cell Segmentation" (normal)
  - "Human Pancreas Preview Data" (non-diseased)
- No Xenium Prime 5K human pancreas dataset page was found. Its existence is UNVERIFIED.

Other sources fetched
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12528271/fullTextXML. Quoted in section 3; it names the 10x PDAC IO dataset URL.
- https://yusuke-imoto-lab.github.io/RECODE/Tutorials/Tutorial_SpatialTranscriptome_Xenium.html : 200
  - "We use spatial transcriptome data of Xenium 10X, FFPE Human Pancreatic Ductal Adenocarcinoma Data with Human Immuno-Oncology Profiling Panel."
  - Code: `INPUT_DIR = 'data/Xenium_V1_Human_Ductal_Adenocarcinoma_FFPE_outs' INPUT_FILE = "cell_feature_matrix.h5"`, which reads `cells.csv.gz`.
  - Output: "AnnData object with n_obs × n_vars = 234857 × 380 obs: 'cell_id', 'x_centroid', 'y_centroid', ...". These counts are after filtering out zero-count cells and genes.
- https://api.github.com/search/repositories?q=xenium_v1_human_pancreas : 200. It found broadinstitute/Xenium_V1_human_Pancreas_FFPE_outs: "This repo contains a test dataset for Celldega. The data was obtained from 10X Genomics Inc". This is the normal-pancreas dataset, not PDAC.
- https://science.ecosyste.ms/projects/161584 and /161568 : BLOCKED, HTTP 402
- Hugging Face mirror (prism-oncology/novae):
  - https://huggingface.co/api/datasets/prism-oncology/novae : 200. Licence tag `license:bsd-3-clause`. This is the mirror repo's licence; the upstream 10x data licence is UNVERIFIED.
  - https://huggingface.co/datasets/prism-oncology/novae/resolve/main/metadata.csv : 200. Rows:
    - `Xenium_V1_Human_Ductal_Adenocarcinoma_FFPE_outs,235099,0,380,0,xenium,pancreas,human`
    - `Xenium_V1_hPancreas_Cancer_Add_on_FFPE_outs,190965,0,474,0,xenium,pancreas,human`
    - `Xenium_V1_human_Pancreas_FFPE_outs,140702,0,377,...`
    - `Xenium_V1_hPancreas_nondiseased_section_outs,103901,0,377,...`
    - `cosmx_pancreas,48944,0,21569,...`
    - The disease status of `cosmx_pancreas` is UNVERIFIED.
  - (`metadata_all.csv` at /resolve/main returned 404 "Entry not found".)
  - https://huggingface.co/api/datasets/prism-oncology/novae/tree/main/human/pancreas : 200. File sizes:
    - `Xenium_V1_Human_Ductal_Adenocarcinoma_FFPE_outs.h5ad`: 333,536,825 B
    - `Xenium_V1_hPancreas_Cancer_Add_on_FFPE_outs.h5ad`: 359,122,103 B
    - `Xenium_V1_human_Pancreas_FFPE_outs.h5ad`: 169,400,146 B
    - `Xenium_V1_hPancreas_nondiseased_section_outs.h5ad`: 165,879,647 B
    - `cosmx_pancreas.h5ad`: 634,518,168 B
  - h5ad range reads (PDAC IO and Cancer Add-on files):
    - Keys `X, layers, obs, obsm, obsp, uns, var, varm, varp`; `obsm` has `X_pca, X_umap, spatial`.
    - var: 380 genes (PDAC IO) and 474 genes (Cancer Add-on).
    - These are Novae-preprocessed objects. Whether `X` holds raw counts is UNVERIFIED.

Gene coverage (from the mirror var lists)
- PDAC IO: 25 of 39 present. Missing: TGFBR1, TGFBR2, SPP1, MIF, KDR, PTN, SDC4, HGF, IL6R, TNFRSF1A, JAG1, PDGFB, PDGFRB, ITGB1.
- Pancreas Cancer Add-on: 13 of 39 present.

Saved files
- `F:\Thesis\data\pdac_panels\10x_xenium_pdac_io\novae_mirror_var_names.txt`
- `F:\Thesis\data\pdac_panels\10x_xenium_pancreas_cancer_addon\novae_mirror_var_names.txt`

---

## 5. Extra candidates (4 checked)

Discovery: an NCBI E-utilities esearch on db=gds for "(xenium OR cosmx OR merfish OR merscope) AND (pancreatic OR pancreas OR PDAC) AND gse[ETYP]" returned 30 hits, followed by esummary. Other human PDAC hits that were not vetted because of the 4-candidate cap: GSE274673 (Xenium PDAC apCAF niches, 6 samples), GSE280634 (Xenium breast panel, treatment-naive PDAC, 8 samples), GSE310433 (Xenium PDA, 8 samples), GSE343070, GSE313662.

### 5a. GSE310352: CosMx, treatment-naive PDAC (SubSeries of GSE310353; PMID 41533516)
- Fetched (all 200): SOFT text for GSE310352, GSE310353 and GSM9294399; suppl directory listing and filelist.txt; first 200 KB of each of the 8 `cell_by_gene_counts.csv.gz` files; first 30 KB of slide1 `cell_metadata.csv.gz`.
- "!Series_overall_design = FFPE sections of treatment naïve PDAC tumor resections were characterized using the CosMx tecnology through the Nanostring technology access program (TAP, 1k universal RNA panel)."
- "characterize primary untreated tumors from 39 patients. We capture 340K spatial low-bulk and 530K spatial single-cell transcriptomes". This covers the whole SuperSeries. The CosMx patient count is UNVERIFIED.
- Platform GPL33484, "NanoString CosMx SMI (human)". 8 samples, GSM9294399 to GSM9294406, "tissue: PDAC tumor resection".
- Data processing note: "Cell-by-gene counts are provided per slide as <slideN>_cell_by_gene_counts.csv.gz (from CosMx exprMat_file.csv). Cell metadata are provided per slide as <slideN>_cell_metadata.csv.gz (from CosMx metadata_file.csv)."
- Sizes (filelist.txt):
  - RAW.tar: 6,331,320,320 B
  - counts files: 96,587,832 B total
  - metadata files: 27,939,772 B total
  - minimal set (counts + metadata): 124,527,604 B
- Counts header: 989 targets = 970 genes + 19 NegPrb, identical in all 8 slides. Metadata header includes CenterX_global_px and CenterY_global_px.
- Gene check: 37 of 39 present. PTN and SDC4 are absent.
- Licence: not stated on fetched pages.
- Saved: `F:\Thesis\data\pdac_panels\geo_gse310352_cosmx\` (8 header files)

### 5b. GSE277782: CosMx, metastatic PDAC (Nature 2025, "Spatial mapping of transcriptomic plasticity in metastatic pancreatic cancer")
- Fetched (all 200): SOFT text for GSE277782 and GSM8529618 to GSM8529624; suppl directory listing; HEAD on minimal files; first 200 KB of both exprMat files; first 30 KB of the metadata and annotation files.
- "!Series_overall_design = CosMx SMI profiling (Nanostring) was performed using a 1000-plex gene panel on fresh frozen tissues. A total of 90 Fields of View (FOVs, each 0.51 x 0.51 mm) ... across seven tissue samples"
- "we performed CosMX SMI profiling (Nanostring) on 7 samples from 3 patients, including primary and/or liver metastasis."
- Samples:
  - Pt-8: Pri, LiM-1, LiM-2
  - Pt-5: Pri, LiM-1, LiM-2
  - Pt-13: Pri
  - Treatments (e.g. "FOLFIRINOX, FOLFIRI, Abraxane/Gemcitabine").
- Sizes from Content-Length:
  - exprMat: Slide1 5,590,303 B; Slide2 7,741,395 B
  - metadata: Slide1 1,348,988 B; Slide2 2,003,027 B
  - Meta.data_CoxMx (annotations): 569,141 B
  - minimal total: 17,254,766 B
- exprMat header: 1000 genes + 10 NegPrb. The annotation file has an `Annotation_main` column.
- Gene check: 37 of 39 present. PTN and SDC4 are absent.
- https://www.nature.com/articles/s41586-025-08927-x : BLOCKED, 303 redirect to idp.nature.com/authorize (login)
- Europe PMC fullTextXML for PMC13242257 : BLOCKED/failed, HTTP 500
- Saved: `F:\Thesis\data\pdac_panels\geo_gse277782_cosmx\` (2 header files)

### 5c. GSE274103: Visium (excluded, not single-cell)
- SOFT text and filelist fetched (200).
- "we performed Visium spatial transcriptomics (10X Genomics) on FFPE samples from five treatment-naïve PDAC patient tissues."
- RAW.tar: 186,603,520 B

### 5d. Cancer Letters 2026, Xenium 5K PDAC primary tumours + hepatic metastases (doi 10.1016/j.canlet.2026.218613)
- https://www.sciencedirect.com/science/article/pii/S0304383526003769 : BLOCKED, "The server returned HTTP 403 Forbidden."
- Europe PMC search: 200. PMID 42214806, no PMC full text, licence "cc by-nc".
- Abstract: "we performed high-definition single-cell spatial transcriptomics using Xenium In Situ (5K) on formalin-fixed paraffin-embedded (FFPE) matched PDAC PTs and HMs, generating a single-cell spatial atlas encompassing more than 900,000 resolved cells."
- Data accession: none found (WebSearch). Availability is UNVERIFIED.

---

## 6. Other BLOCKED fetches
- https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html : curl got a reCAPTCHA challenge page. WebFetch reported "only a browser verification notice". The GEO data-use terms are therefore UNVERIFIED for all GEO candidates.

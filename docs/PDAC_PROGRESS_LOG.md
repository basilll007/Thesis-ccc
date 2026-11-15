# PDAC Extension — Progress Log

Append-only log. Every number here comes from a file listed in `results/pdac/MANIFEST.csv`. Interpretation is limited to audit and diagnostic statements. Nothing here supports a biological or clinical claim.

## Task checklist

- [x] Phase 0: Repo hygiene (PDCD1 caption, F:\ paths, settings note, "replication" renamed)
- [x] Phase 1: Matched dyadic baseline (breast) — gate: matched baseline does NOT explain most C/E disagreements (3/9), but it carries all edge-significant calls (0/8 ratio-significant)
- [x] Phase 2: Find and vet a PDAC spatial dataset — selected Mendeley PDAC-SMI (CosMx, CC BY 4.0, 37/39 genes)
- [x] Phase 3: Run the audit on PDAC — DONE except full-catalogue radius sweep (host stopped run for low memory; headline-only radius sweep done; full sweep BLOCKED pending user go-ahead)
- [x] Phase 4: Benchmark — COMMOT DONE (crop), SpatialDM DONE, Squidpy ligrec DONE (crop) / full-section p-values BLOCKED (all NaN, 2 attempts)
- [x] Phase 5: Synthesis (PDAC_SUMMARY.md, README section)

## Conventions used throughout (decided 2026-10-24)

- Python: conda env `thesis-spatial-ccc` (Python 3.11.10, scanpy 1.10.3, squidpy 1.6.5, anndata 0.10.9, numpy 1.26.4, pandas 2.2.3, scipy 1.14.1). Exact versions are written to each run log.
- New scripts live in `src/pdac/`. They use `src/pdac/common.py`, which writes `logs/<script>_<timestamp>.log` (package versions and git commit), registers outputs in `results/pdac/MANIFEST.csv`, and saves figures as PNG + PDF without overwriting.
- Seed 42 everywhere. New label permutations use `numpy.random.default_rng(42)` through the existing `spatial_null.generate_label_permutations`, so they are the same stream as the paper's 6-NN nulls.
- p-values: new tables report both the legacy convention `mean(null >= obs)` (to compare with existing files) and the add-one convention `(k+1)/(n+1)`. Add-one is primary for new claims and for BH. Reason: the legacy convention can return p = 0, which is not a valid permutation p-value.

---

## 2026-10-24 21:05 — Phase 0 start: branch and environment

- **Done**: Created branch `pdac-extension` from `main` @ 3cde438. Confirmed the project env `thesis-spatial-ccc` imports scanpy/squidpy/anndata. The base Anaconda env is broken for this code (`np.unicode_` removed under NumPy 2), so it is not used.
- **Command**: `git checkout -b pdac-extension`
- **Next**: Phase 0 hygiene items.

## 2026-10-24 21:10 — Phase 0.1: PDCD1 caption mismatch

- **Done**: `ablation.html` (Figure 4(a) caption) said "12 out of 7,163 cells (0.17%)". `results/cd274_pdcd1_detection_diagnostic.csv` says PDCD1 `ALL` n_detected = 5 of 7,163. Recomputed from `data/processed/xenium_breast_baseline.h5ad` (`adata.raw.X > 0`): **5/7,163**. The figure itself (`src/figures.py::make_biology_spatial_expression`) computes its count dynamically, so only the HTML caption was wrong. Caption changed to "5 out of 7,163 cells (0.07%)" with the source CSV cited.
- `docs/PROGRESS_LOG.md` line ~254 (2026-09 entry) also says 12/7,163. That historical entry is kept, and a correction note is appended to PROGRESS_LOG.md.
- **Output**: `results/pdac/phase0/pdcd1_cd274_detection_recheck.csv` (CD274 53/7,163 and PDCD1 5/7,163; both match the existing CSV).

## 2026-10-24 21:12 — Phase 0.2: hard-coded absolute paths

- **Done**: Replaced all `r"F:\Thesis\..."` literals (20 occurrences in 5 scripts: `scgpt_step0_overlap_check.py`, `scgpt_step1_embedding_gate.py`, `scgpt_step2_embedding_comm_score.py`, `gnn_step1_link_pred_gate.py`, `gnn_step2_pca_comm_score.py`) with `REPO_ROOT / ...`, where `REPO_ROOT = Path(__file__).resolve().parents[1]`. Every `src/*.py` passes `py_compile`.
- Old `F:\` strings in historical log files (`logs/spatial_null.log`, `logs/spatial_edge_null.log`) are records of past runs and were left as they are.
- The scGPT/GNN scripts were not re-run because they need torch / scGPT envs and their outputs feed the submitted paper. Only the path strings changed.

## 2026-10-24 21:14 — Phase 0.3: separate the 15 µm results from the radius-sweep peaks

- **Done**: Wrote `docs/RESULTS_SETTINGS_NOTE.md`, which separates (a) contact-restricted 15 µm results (paper), (b) radius-sweep peaks (sensitivity only), and (c) 6-NN null comparison, with the source file for each number.
- Added a qualifier to the `enhancements.html` badge: "z = +19.96, 200 µm radius-sweep peak; 15 µm contact: z = +3.56".
- **Reproducibility check** (`src/pdac/phase0_repro_check.py`, log `logs/phase0_repro_check_20261024T211538.log`): recomputed the 15 µm radius graph (40,618 directed edges) with the identical algorithm and RNG (`RandomState(42)`, 200 permutations). All four axes match `results/multi_axis_expansion/multi_axis_summary.csv` exactly:
  - CD86→CTLA4 immune→immune: z = 3.56, p_legacy = 0.035 (p_add-one = 0.0398)
  - CXCL12→CXCR4 fibroblast→immune: z = 36.20, p_legacy = 0.000 (p_add-one = 0.0050)
  - PTN→SDC4 fibroblast→epithelial: z = −1.36, p = 0.920
  - CD274→PDCD1 epithelial→immune: z = 0.00, p = 1.000
- **Outputs**: `results/pdac/phase0/multi_axis_15um_recheck.csv`, `results/pdac/phase0/phase0_repro_summary.json`, `figures/pdac/phase0/phase0_multi_axis_15um_recheck_nulls.{png,pdf}`.
- **Decision**: Paper numbers are reproducible from the current code and data. No existing result file was modified.

## 2026-10-24 21:18 — Phase 0.4: "replication" renamed to "internal regional consistency check"

- **Done**: Renamed the margin/core split wording in `enhancements.html` (section badge, title, intro, figure card, overlay, alt text, and the item (7) summary), `index.html`, `methods.html`, `docs/PROGRESS_LOG.md` (2026-09-24 entry), `src/run_cross_fov_replication.py` (docstring, JSON dataset label, print), and `src/visualize_pcst_and_replication.py` (docstrings and figure titles).
- **Deliberately not renamed**: file and directory names (`results/cross_fov_replication/`, `src/run_cross_fov_replication.py`, `figures/cross_fov_replication_analysis.*`) and the JSON key `kinetics_replication`, because renaming them would break links and the existing outputs. "Replication" in `docs/aaai_student_abstract_draft.md` (the submitted paper; "pending replication" and "no cross-patient replicate") refers to future external data, not the margin/core split, so it is out of scope and not edited. `ablation.html:120` uses "replicates" for a GNN failure mode, which is unrelated.
- **Next**: Commit Phase 0, then Phase 1.

## 2026-10-24 21:20 — Phase 1: matched dyadic baseline (breast, 6-NN, 500 perms)

- **What was done**: Added `src/pdac/ccc_audit.py`, a vectorised scorer. It computes four scores for every (L–R pair, sender, receiver) on the same permuted label vectors:
  - compositional C = (L̄_s + R̄_r)/2
  - edge E = mean over s→r graph edges of L_i·R_j
  - matched dyadic baseline M = L̄_s·R̄_r (product form, no graph)
  - ratio Q = E/M

  C vs M differ only in score form; M vs E differ only in spatial restriction. `src/pdac/phase1_matched_baseline.py` runs three families on the stored 6-NN graph (42,978 directed edges):
  - **paper**: the 22 rows of `results/baseline_lr_ranked.csv`
  - **extra**: CD86–CTLA4 and PTN–SDC4, all 16 sender–receiver combinations each (32 rows)
  - **leiden**: paper pairs under Leiden-consensus labels (18 rows, compositional score > 0)

  Permutations: `spatial_null.generate_label_permutations(..., 500, seed=42)`, the same stream as the paper.
- **Commands**: `python src/pdac/phase1_matched_baseline.py` (log `logs/phase1_matched_baseline_20261024T212000.log`); `python src/pdac/phase1b_gate_details.py` (log `logs/phase1b_gate_details_*.log`).
- **Reproduction check** (`results/pdac/phase1/reproduction_check_vs_paper_files.json`): all 22 rows matched. Edge obs, z and p equal the paper files to within 1e-13 (p exact). Compositional obs and p are exact. Compositional z differs by at most 2.7e-4, probably a dtype difference in the original run; no verdict changes.
- **Key numbers, paper family** (`results/pdac/phase1/decomposition_summary.csv`, `phase1_gate_details.json`):
  - Significant (add-one p < 0.05): C 17/22, M 14/22, E 8/22, Q 2/22. BH q < 0.05 within the family of 22 gives the same counts.
  - C/E verdict disagreements: 9 (the paper's 9 flips). Under the pre-registered rule, 3 are score-form (M agrees with E) and 6 are proximity (M agrees with C).
  - 4 of the 9 are CD274–PDCD1, a sparse-gene diagnostic (PDCD1 detected in 5 cells). Among the 5 CXCL12–CXCR4 disagreements, 1 is score-form (epithelial→immune) and 4 are proximity. 2 of those 4 have fewer than 30 qualifying edges (immune→endothelial 17, endothelial→endothelial 27). For them, "proximity" cannot be told apart from the loss of power that comes from averaging over few edges.
  - **8 of 8 edge-significant interactions are also matched-significant, and 0 of 8 have a significant ratio Q.** Q is significant only in 2 rows, both with a strongly negative matched z (epithelial→immune, z_M = −7.23; epithelial→fibroblast, z_M = −21.12). That is the pattern of a ratio with a small denominator.
- **Extra family**: E significant in 4/32, all 4 matched-significant, 0/4 ratio-significant. 5 C/E disagreements, all CD86–CTLA4 (1 score-form, 4 proximity). Q significant in 4 PTN–SDC4 rows, all with negative matched z.
- **Leiden family**: 6 C/E disagreements, all CD274–PDCD1 (sparse-gene; 3 form, 3 proximity). CXCL12–CXCR4 has no C/E disagreements under Leiden labels. 1 of 7 edge-significant rows is ratio-significant (CXCL12 epithelial→immune, z_M = −10.03).
- **Gate (stated plainly)**: By the pre-registered rule, the matched baseline does **not** explain most C-vs-E verdict disagreements: 3 of 9 (`results/pdac/phase1/phase1_gate.json`). Two qualifications are just as important:
  - (i) Almost half of the disagreements (4/9) are uninformative sparse-gene cases, and 2 of the remaining "proximity" cases rest on fewer than 30 edges.
  - (ii) Where the edge score **is** significant, the matched non-spatial baseline accounts for that significance in every case (8/8 paper, 4/4 extra). Restricting to graph edges adds no detectable enrichment beyond the cell-type mean product (Q stays within its null). The edge score's positive calls therefore mainly reflect composition in product form, not proximity-specific coupling.
  - On this dataset, the edge score's value as an audit tool lies in **removing** compositional calls (especially zero-contact or low-edge ones). There is no evidence here that it **adds** proximity-specific positives.
- **Caveat**: Q = E/M is unstable when M is small. All 8 Q-significant rows across the families have negative matched z, so Q significance should not be read as evidence of proximity coupling.
- **Outputs**: `results/pdac/phase1/` (per-family CSVs, `null_distributions_{paper,extra,leiden}.npz`, `decomposition_summary.csv`, `disagreements_detail.csv`, gate JSONs); `figures/pdac/phase1/` (5 figures from phase1 plus 1 from phase1b, PNG + PDF).
- **Next**: Phase 2. Dataset vetting is running in parallel; then Phase 3.

## 2026-10-24 21:35 — Phase 2: PDAC dataset vetting

- **What was done**: A delegated web-vetting pass fetched the candidate pages and saved verbatim evidence (`results/pdac/phase2/dataset_vetting_raw_notes.md`, draft table `results/pdac/phase2/dataset_vetting_draft.csv`). Small panel and header files were stored git-ignored under `data/pdac_panels/`. I then independently re-verified the selected dataset:
  - Mendeley public API (`https://data.mendeley.com/public-api/datasets/kx6b69n3cb`): licence **CC BY 4.0**; files `raw_meta_data_final.h5ad` 300,718,306 B, `SMI_raw.zip` 5,519,589,653 B, `Co-culture_seurat.rds` 585,225,537 B.
  - Downloaded `raw_meta_data_final.h5ad` (plain HTTPS, size matches exactly, sha256 `7c65c071...`) to `data/pdac/mendeley_pdac_smi/` (git-ignored).
  - Inspected it (`results/pdac/phase2/mendeley_h5ad_inspection.json`): 717,493 cells × 1,009 features (960 Base, 30 Spike-in, 19 Negprob); 16 sections (`sample_id`); integer counts (max 521); author annotations (`annotation_majortypes` / `annotation_subtypes`); centroids `CenterX/Y_global_px`. The 30 spike-in genes are NaN in every cell of the "-b" sections (284,492 cells), i.e. not measured there. They are dropped from all analyses, together with the 19 negative probes.
  - Gene presence on the 960 Base probes, re-checked from the file, agrees with the draft: 37/39 catalogue genes present; PTN and SDC4 absent.
  - Source paper, fetched via the Europe PMC REST API (`data/pdac_panels/mendeley_pdac_smi/PMC10327107_fulltext.xml`): Shiau, Cao, Gregory, Gong et al., "Therapy-associated remodeling of pancreatic cancer revealed by single-cell spatial transcriptomics and optimal transport analysis", bioRxiv 2023, doi 10.1101/2023.06.28.546848. It states "960-plex gene panel, with (n = 6) and without (n = 9) 30 custom probes" and "753,133 high confidence single-cell profiles" (the h5ad holds 717,493; the difference is unexplained and the h5ad count is used). The data-availability statement points to this Mendeley DOI.
  - **Pixel scale**: the paper's Methods say "60 pixels (10.8 μm)", i.e. 0.18 µm/px, used here to convert `CenterX/Y_global_px` to µm. Local coordinates reach 5,460 × 3,640 px, which fits a 5,472 × 3,648 px FOV. Caveat: the sentence is about a downsampled PanCK image, so the scale is inferred, not explicitly stated for the centroid columns. Only the radius sweep depends on it; kNN graphs are scale-free.
- **BLOCKED**:
  - PMC HTML (reCAPTCHA), so Europe PMC's public API was used instead.
  - GEO usage-terms page (reCAPTCHA), so licences for GSE310352, GSE277782 and GSE267680 are not stated.
  - All 10x Genomics pages (HTTP 429), so 10x licences, bundle sizes and URLs are UNVERIFIED.
  - Cancer Letters 2026 Xenium 5K article (HTTP 403).

  None of these blocks was worked around.
- **Excluded**: bioRxiv 10.1101/2025.02.11.637627 (no new spatial data; re-uses the 10x PDAC IO Xenium set plus one GSE267680 sample); GSE274103 (Visium, spot-level).
- **Ranking** (single-cell > L–R coverage > licence > downloadability):
  1. Mendeley PDAC-SMI
  2. GSE310352 (CosMx, 37/39, licence not stated)
  3. GSE277782 (CosMx, 37/39, mostly liver metastases, licence not stated)
  4. 10x Xenium PDAC IO (25/39, 10x page blocked)
  5. GSE267680 (PanIN, 5.5 GB per sample, panel UNVERIFIED)
  6. 10x Xenium pancreas cancer add-on (13/39)
- **Decision**: Mendeley PDAC-SMI. It is single-cell, has 37/39 genes, an explicit CC BY 4.0 licence, a 301 MB file, and 16 sections from 13 tumours, which allows a cross-section consistency check.
- **Outputs**: `results/pdac/dataset_vetting.csv` (final), `results/pdac/phase2/*`, `figures/pdac/phase2/phase2_lr_gene_coverage_by_candidate.{png,pdf}`, log `logs/phase2_dataset_vetting_20261024T213443.log`.
- **Phase 4 prep (same session)**: Cloned env `thesis-spatial-ccc` → `pdac-bench`, so benchmark installs cannot alter the paper env. `pip install commot` installed commot 0.0.3 but upgraded numpy to 2.2.6, and import then failed (`AttributeError: np.Inf was removed in the NumPy 2.0 release`). Attempt 2 pinned `numpy<2` (1.26.4): commot imports (pip warns that tobler 0.13.0 wants numpy>=2; tobler is not used). `pip install SpatialDM` installed SpatialDM 0.3.1, which imports. Logs: `logs/phase4_install_commot_*.log`, `logs/phase4_install_spatialdm_*.log`.
- **Next**: Phase 3. Pre-declared before loading any PDAC expression into the pipeline:
  - Primary section = the largest untreated section (by the rule: U5-a, 74,793 cells).
  - All 16 sections get the 6-NN main analysis as a cross-section consistency check.
  - QC thresholds are the breast ones (min_counts=10, min_cells=3). This deviates from the authors' ≥60-transcript filter, deliberately, to keep the pipeline identical.
  - Marker sets and the L–R catalogue are fixed in `src/pdac/phase3_pdac_audit.py` and `src/pdac/lr_catalog.py`.

## 2026-10-24 21:37–22:20 — Phase 3: audit on PDAC (Mendeley PDAC-SMI, CosMx)

- **Scripts**:
  - `src/pdac/phase3_pdac_audit.py --dataset mendeley_pdac_smi:U5-a` (primary section; log `logs/phase3_pdac_audit_mendeley_pdac_smi__U5-a_*.log`)
  - `src/pdac/phase3b_cross_section.py` (all 16 sections)
  - `src/pdac/phase3c_summarise.py` (summaries and figures from saved tables)
  - `src/pdac/phase3d_author_annotation.py` (authors' labels as a reference annotation)
  - Shared: `src/pdac/lr_catalog.py`, `src/pdac/pdac_loaders.py`
- **Data handling**: one section at a time; 960 Base genes (30 spike-ins and 19 negative probes dropped); coordinates px × 0.18 → µm. QC with breast thresholds (min_counts = 10, min_cells = 3) removed **0 cells and 0 genes** in every section (`results/pdac/phase3/cross_section_section_info.csv`), because the authors had already filtered cells (≥ 60 transcripts in their Methods). This is not a deviation in thresholds, but QC was effectively inert here.
- **Annotations**:
  - (1) Marker scoring: per-cell argmax of mean marker expression over 10 PDAC types, markers declared before the data were loaded; none dropped.
  - (2) Leiden consensus: res 0.5, same calls as breast. In U5-a this gives only 4 types (CAF_fibroblast, endocrine, epithelial_ductal_tumor, myeloid). **No T-cell class**, so the CAF→T headline cannot be tested under Leiden labels.
  - (3) Added during the run: the authors' `annotation_majortypes`, mapped to our coarse types (`AUTHOR_MAP` in `phase3d_author_annotation.py`) as a reference. Reason: the crosstab `mendeley_pdac_smi__U5-a_author_vs_marker_crosstab.csv` shows our marker "T_cell" class (945 cells) contains 179 CD8+ T, 20 CD4+ T, 62 Treg and 29 NK author calls (290 T/NK in total), plus 254 "Malignant" and 221 "Noise". Of the authors' 879 T cells (CD4+ 50, CD8+ 476, Treg 353), our marker scoring labels 467 as CAF_fibroblast. The marker T-cell label is impure.
- **Coverage** (`mendeley_pdac_smi__U5-a_lr_coverage.csv`): 31/33 catalogue pairs testable; absent: PTN–SDC4, CXCL13–CXCR5. **No tested gene is detected in fewer than 20 cells** in U5-a (`mendeley_pdac_smi__U5-a_sparse_genes.csv` is empty), so no axis carries the sparse-gene flag. Receptor detection counts are in `mendeley_pdac_smi__U5-a_receptor_detection_counts.csv`.
- **Multiple testing (U5-a, 6-NN, 500 perms; `mendeley_pdac_smi__U5-a_phase3_summary.json`)**: **3,596 rows tested** = 3,100 marker (31 pairs × 10 × 10 types) + 496 Leiden (31 × 4 × 4). Every row had ≥ 1 qualifying edge.

  | Annotation | Significant (uncorrected p < 0.05, add-one p): C / M / E / Q | BH q < 0.05 within annotation: C / M / E / Q | BH across all 3,596 rows: E / Q |
  |---|---|---|---|
  | marker | 724 / 645 / 385 / 405 | 643 / 544 / 214 / 186 | 239 / 207 |
  | Leiden | 164 / 137 / 155 / 147 | — | — |

- **C/E decomposition (U5-a, marker)**: 495 C/E disagreements; 106 score-form and 389 proximity. 0 involve a sparse gene. Of 385 edge-significant rows, 303 are also matched-significant and **141 are also ratio-significant**. Unlike breast (0/8), PDAC shows proximity enrichment beyond the matched baseline in a sizeable minority of edge-significant rows.
- **Headline CXCL12→CXCR4** (U5-a):
  - CAF→T_cell (marker): 2,101 edges; C z = 33.17, M z = 36.05, E z = 23.38, Q = 1.62 with z_Q = 2.81 (p_add-one 0.0080). After BH across all 3,596 rows: E q = 0.022 (significant), Q q = 0.082 (not significant).
  - CAF→epithelial_ductal_tumor: depleted under both marker labels (E z = −15.22, Q z = −8.85; 30,053 edges) and Leiden labels (E z = −18.65); not significant.
- **Cross-section consistency check (16 sections; `cross_section_summary_by_interaction.csv`, `author_annotation_headline_summary.csv`)**:
  - CAF→T_cell, **marker** labels: C 16/16, M 16/16, E 16/16 (13/16 after within-section BH), **Q 4/16 (0/16 after BH)**; median E z = 20.87, median Q = 1.35 (median z_Q = 0.75).
  - CAF→T_cell, **author** labels: C 16/16, M 16/16, E 16/16 (14/16 BH), **Q 2/16 (0/16 BH)**; median E z = 21.63, median Q = 1.35.
  - CAF→tumour: E 0/16 under marker, Leiden and author labels. C 3/16 (marker), 3/16 (Leiden), 2/16 (author).
  - All sections, marker labels: 49,569 section-rows; 6,058 edge-significant, of which 4,676 matched-significant and 2,127 ratio-significant. 22 interactions are edge-significant in all 16 sections; 21 interactions are ratio-significant in ≥ 12/16 (`cross_section_ratio_sig_in_ge12_sections_marker.csv`).
- **Annotation sensitivity (U5-a, marker vs Leiden, 496 shared rows; `mendeley_pdac_smi__U5-a_annotation_sensitivity_summary.csv`)**:

  | Score | Spearman z | median abs rank shift | verdict changes |
  |---|---|---|---|
  | C | 0.92 | 18 | 49 |
  | M | 0.89 | 22 | 50 |
  | E | 0.83 | 29 | 77 |
  | Q | 0.65 (484 rows) | 48 | 106 |

  The ratio score is the most annotation-sensitive.
- **k sensitivity (U5-a, marker; `mendeley_pdac_smi__U5-a_k_sensitivity_summary.csv`)**: edge-significant rows at k = 4 / 6 / 8 / 12 = 378 / 385 / 413 / 433. Spearman of edge z vs k = 6: 0.955 / 1 / 0.974 / 0.945. Edge verdict changes vs k = 6: 105 / 0 / 80 / 134; ratio verdict changes: 195 / 0 / 161 / 239. The headline CAF→T rows across k are in `mendeley_pdac_smi__U5-a_k_sensitivity_headline.csv`.
- **Radius sensitivity: PARTIAL.** The primary run was **stopped by the host (Claude Code) because the system was critically low on memory** at the 200 µm radius step (about 12M edges after receiver subsampling; the 150 and 200 µm graphs used query-cell fractions of 0.538 and 0.321). Null distributions for 15–150 µm were saved (`*_null_radius{15..150}_marker.npz`). The per-radius observed table, the figures and the processed AnnData were held in memory and lost. Per the host's instruction, the full run was **not restarted automatically**. A low-memory, headline-only radius sweep is done separately (next entry). The full-catalogue radius sweep needs to be re-run when memory allows, and that requires the user's go-ahead.
- **Edge lengths** (`mendeley_pdac_smi__U5-a_6NN_edge_length_summary.json`): 6-NN max edge 74.4 µm; 3.6e-5 of edges are longer than 50 µm. FOV-gap bridging is negligible.
- **Interpretation (audit level)**: The headline CXCL12→CXCR4 CAF→T row is significant under C, M and E in every section and under both annotations that have a T-cell class. However, its edge significance is largely matched by the non-spatial matched baseline (M), and proximity enrichment beyond that baseline (Q) is significant in only 4/16 (marker) or 2/16 (author) sections, with 0/16 surviving BH. The marker T-cell label is substantially impure. Nothing here supports a biological claim about CAF–T-cell signalling.

## 2026-10-24 22:43–22:49 — Phase 3e: headline-only radius sensitivity (low memory; SENSITIVITY ONLY)

- **Done**: `src/pdac/phase3e_radius_headline.py` (log `logs/phase3e_radius_headline_20261024T224325.log`) ran the 15–300 µm radius sweep for CXCL12–CXCR4 only (all 100 sender/receiver combinations, marker labels, 500 perms), using the same receiver-subsampling RNG sequence as the stopped run. Above 12M edges, receivers are subsampled: query-cell fraction 0.538 at 150 µm, 0.321 at 200 µm, 0.163 at 300 µm.
- **Reproducibility**: for 15–150 µm, the edge nulls are identical (max abs diff 0) to the nulls saved by the stopped full run (`mendeley_pdac_smi__U5-a_radius_headline_repro_check.json`).
- **Key numbers** (`mendeley_pdac_smi__U5-a_radius_sensitivity_headline_pair.csv`):
  - CAF→T_cell is edge- and ratio-significant at every radius. Q falls from 1.742 (15 µm) to 1.494–1.532 (150–300 µm), while z_Q rises from 2.09 to 6.84 (200 µm) as qualifying edges grow from 1,455 to about 59,000.
  - CAF→tumour is depleted at every radius (E z −9.7 to −20.7).
  - k sweep for the headline (`mendeley_pdac_smi__U5-a_k_sensitivity_headline.csv`): CAF→T z_Q = 2.54 / 2.81 / 3.24 / 4.11 at k = 4 / 6 / 8 / 12; edges 1,429 → 4,274.
- **Interpretation**: As in breast, larger neighbourhoods inflate z through edge count, while the effect size Q does not grow. Radius and k z-values must not be read as stronger communication.
- **Still not done**: full-catalogue radius sweep (all 31 pairs) beyond what the stopped run saved. Status: **BLOCKED pending a user go-ahead to re-run** (host memory reaping). Saved null files: `*_null_radius{15..150}_marker.npz`.

## 2026-10-24 22:21–22:43 — Phase 4: benchmark vs Squidpy ligrec, COMMOT, SpatialDM

- **Inputs (pre-declared)**:
  - Primary section U5-a with marker labels. Labels were rebuilt by `phase4a` without PCA/Leiden and are identical to Phase 3 (same counts).
  - **crop**: the 8,000 cells nearest to the cell closest to the median coordinate (extent 1,800 × 2,739 µm, because FOVs are discontiguous; `benchmark_crop_info.json`). The crop exists because COMMOT builds a dense n × n distance matrix (`commot/tools/_spatial_communication.py`: `distance_matrix(adata.obsm["spatial"], ...)`), which is infeasible for 74,793 cells.
  - **full**: the whole section, for our scores, ligrec and SpatialDM.
- **Scripts and logs**:
  - `src/pdac/phase4a_ours_and_squidpy.py` (env thesis-spatial-ccc)
  - `src/pdac/phase4a2_squidpy_serial.py`
  - `src/pdac/phase4b_commot_spatialdm.py` (env pdac-bench; the log shows commot `__version__` as "unknown" because the package does not define it; `pip show` = 0.0.3; SpatialDM 0.3.1)
  - `src/pdac/phase4c_compare.py`
- **Parameters**:
  - Squidpy ligrec: 500 perms, seed 42, threshold 0.01, use_raw.
  - COMMOT: our 31 pairs, `dis_thr` = 50 µm, defaults otherwise; `cluster_communication` with 100 perms (package default), random_seed 42; 544 s.
  - SpatialDM: `weight_matrix(l=25 µm, cutoff=0.2, single_cell=True)`, built-in CellChatDB (`datahost='package'`, no network), `spatialdm_global(method='z-score')`; 565 CellChatDB pairs. 26 of our 31 pairs map exactly, 5 via receptor-in-complex (`spatialdm_mapping_full.csv`).
- **Status** (`phase4b_status.json`, `squidpy_ligrec_attempts.json`):
  - COMMOT DONE (crop).
  - SpatialDM DONE (crop and full).
  - Squidpy ligrec DONE on crop (1,868/3,100 p-values NaN by design, from the 1% expression mask; 203 rows p < 0.05).
  - **Squidpy ligrec on the full section: BLOCKED.** Attempt 1 (default, `n_jobs=4`, numba-parallel kernel because n ≥ 2^15) and attempt 2 (`numba_parallel=False`, `n_jobs=1`) both returned NaN for all 3,100 p-values, with identical means. The input has no NaN/inf and unique cell names. Cause not identified. Per the two-failure rule, no third attempt. Full-section ligrec **means** are used for rank correlation; its significance calls are excluded.
  - Crop control: the serial re-run reproduces the NaN mask exactly. p-values differ by up to 0.072 because Squidpy's per-job seeds depend on `n_jobs`.
- **Agreement, crop, cluster level (3,100 rows; `comparison_spearman_crop.csv`, `comparison_top20_jaccard_crop.csv`, `comparison_significance_agreement_crop.csv`)**:
  - Spearman of our edge z with: ligrec mean 0.022, COMMOT cluster score 0.118. Our compositional z with ligrec mean: 0.258. ligrec vs COMMOT: 0.556.
  - Top-20 Jaccard between any of our scores and ligrec or COMMOT is ≤ 0.053; ligrec vs COMMOT is 0.026.
  - Significance kappa:

    | Pair | kappa |
    |---|---|
    | ligrec vs our compositional | 0.449 |
    | ligrec vs our edge | 0.004 |
    | COMMOT vs our edge | 0.443 |
    | COMMOT vs our compositional | 0.387 |
    | COMMOT vs ligrec | 0.125 |

    So ligrec agrees with the compositional score and COMMOT agrees with the edge score, as expected from their designs.
- **Agreement, full section, cluster level (`comparison_spearman_full.csv`)**: Spearman of ligrec mean with our compositional z 0.098, edge z −0.081, ratio z −0.197.
- **L–R level (31 pairs; `comparison_lr_level_spearman_{crop,full}.csv`)**:
  - Full: SpatialDM global z vs our max edge z 0.260, vs ligrec max mean 0.542, vs our max compositional z 0.038.
  - Crop: COMMOT total signal vs our max edge z 0.603 and vs our max ratio z 0.641; SpatialDM vs COMMOT 0.294.
- **Interpretation (audit level)**: The four method families do not agree on the top interactions (top-20 Jaccard ≤ 0.053 for every cross-family pair). Agreement follows score design: non-spatial cluster-mean methods (ligrec) agree with our compositional score; contact or transport methods (COMMOT) agree with our edge score. The rankings depend on the method choice, which is itself an audit-relevant result.

## 2026-10-24 22:51–22:53 — Phase 3f: background-level detection check (added diagnostic)

- **Why**: In U5-a, the checkpoint genes had implausibly high detection for T-cell-restricted genes (`breast_vs_pdac_axes.csv`: PDCD1 in 7,560 cells, CTLA4 in 7,427, CD274 in 13,426 of 74,793). The breast "fewer than 20 cells" sparse-gene rule cannot catch the inverse problem, where non-specific background counts make a gene look abundant.
- **Done**: `src/pdac/phase3f_background_check.py` compared each catalogue gene's detection rate with the 19 CosMx negative-control probes per section. Flag `near_background` = detection < 2 × the maximum negative-probe detection; the threshold was declared in the script before running. This is a diagnostic flag, not a filter.
- **Key numbers** (`background_check_summary_by_gene.csv`, `background_check_celltype_specificity_median.csv`, `background_check_summary.json`):
  - 33 catalogue genes are near background in ≥ 8/16 sections. They include **every checkpoint gene**, with median fold over the negative-probe max of: PDCD1 0.96, CTLA4 1.00, PDCD1LG2 0.99, CD86 1.08, HAVCR2 1.12, CD28 1.17, CD274 1.37, CD80 1.42.
  - CXCL12 (fold 2.32; near background in 4/16 sections) and CXCR4 (2.13; 6/16) are above background in most sections.
  - Median detection in author T cells vs author malignant cells vs negative-probe max: PDCD1 0.115 / 0.073 / 0.072; CTLA4 0.239 / 0.068 / 0.072; CXCR4 0.560 / 0.076 / 0.072. Detection in malignant cells sits at the negative-probe level, so any T-cell enrichment is modest for PDCD1 and clearer for CTLA4 and CXCR4.
- **Decision**: All PDAC checkpoint-axis results (CD274–PDCD1, PDCD1LG2–PDCD1, CD86/CD80–CTLA4, CD80/CD86–CD28, LGALS9–HAVCR2) carry a **"background-level gene diagnostic"** label. Their significance calls are not interpretable as signal. The rule, and the per-axis result, is: an axis is flagged if either gene is near background in ≥ 8/16 sections. Applied to the 31 testable pairs (`results/pdac/phase3/background_flag_per_axis.csv`): **24/31 flagged**. 7 unflagged: TGFB1–TGFBR1, TGFB1–TGFBR2, CXCL12–CXCR4, JAG1–NOTCH1, MIF–CD74, SPP1–CD44, FN1–ITGB1.

## 2026-10-24 22:55 — Phase 5: synthesis

- **Done**:
  - `src/pdac/phase5_synthesis_tables.py` → `results/pdac/phase5/breast_vs_pdac_axes.csv` and `figures/pdac/phase5/phase5_breast_vs_pdac_axes_z.{png,pdf}`. Axis mapping is declared in the script docstring.
  - Wrote `docs/PDAC_SUMMARY.md` (results table, matched-baseline verdicts, sensitivities, failures and blocks, limitations, 5 next steps).
  - Added the "PDAC extension" section to `README.md`.
- **Final status**:
  - Phases 0, 1, 2 and 5: DONE.
  - Phase 3: DONE, except the full-catalogue radius sweep (BLOCKED: the host stopped the run for low memory; a re-run needs the user's go-ahead; the headline-only sweep is done).
  - Phase 4: DONE for COMMOT, SpatialDM and crop-level ligrec. Full-section ligrec p-values BLOCKED (all NaN in 2 attempts).

## 2026-10-25 — Branch consolidated into `main`

- **Done**: At the user's request, no separate branch is kept. `main` was fast-forwarded to the PDAC work (3cde438 → 40c9707; `git merge --ff-only pdac-extension`) and the local branch `pdac-extension` was deleted (`git branch -d pdac-extension`). All six phase commits are now on `main` with unchanged history. Nothing was pushed. Mentions of the branch name in the README, the summary, the log titles and the module docstrings were updated. The 2026-10-24 entries above that describe creating the branch are kept as a historical record.

---

# PDAC v2 (branch `pdac-v2`, from `main` @ df3867c)

## v2 task checklist

- [x] Phase A: Housekeeping — tag pushed; 5/5 headline numbers match; no uncorrected verdict flips (BH-level flips documented)
- [x] Phase B: DONE — pre-registered protein rule FAILED validation (T precision 0.022 in U5-a); post-hoc consensus_T: CAF→T Q 2/16 uncorrected, 0/16 after BH (does not survive)
- [ ] Phase C: BLOCKED — GSE310352 licence not confirmed (GEO terms reCAPTCHA, publisher 403); author labels only in seurat.RDS (no R on host). Gated script ready; items needed listed in the Phase C entry
- [x] Phase D: DONE — docs/PDAC_SUMMARY_v2.md, synthesis figure, README v2 links (push of commits after 161418a BLOCKED: TLS read errors on large upload)

## 2026-10-25 18:15 — v2 start: reading of the brief

- **Brief**: the pasted v2 brief is cut off in the middle of Phase C, after "write what you need from me (exact items to paste f". It resumes with the end of a synthesis question ("…validated labels and an independent dataset, or is the edge score mostly mirroring cell-type products?") and "Update the README PDAC section with links to v2."
- **Reading taken** (per the "do not stop to ask" rule):
  - Phase C: try the public file listing and the paper or data repository through normal fetches. If the licence or format cannot be confirmed, mark it BLOCKED and write down exactly what the user needs to paste.
  - Phase D: write `docs/PDAC_SUMMARY_v2.md`. It answers "Does the CAF→T proximity-specific (Q) signal survive validated labels and an independent dataset, or is the edge score mostly mirroring cell-type products?" and updates the README with v2 links.
- **Environment**:
  - conda env `thesis-spatial-ccc` (`D:/Anaconda/envs/thesis-spatial-ccc/python.exe`), the same as v1.
  - Free RAM at start: 6.6 GB of 32.6 GB, so the work is chunked per section.
- **Command**: `git checkout -b pdac-v2`

## 2026-10-25 18:18 — Phase A.1: tag `aaai27-submitted`

- **Done**:
  - `git ls-remote --tags origin` listed no tags (exit 0), so the tag did not exist on the remote.
  - Ran `git tag aaai27-submitted 3cde438` (3cde438ec565…, "feat: S-PCST signaling backbone optimization…").
  - Pushed **only the tag**: `git push origin refs/tags/aaai27-submitted`. The remote now lists `refs/tags/aaai27-submitted` → 3cde438. No branch was pushed.

## 2026-10-25 18:20 — Phase A.2: headline-number spot-check

- **Command**: `python src/pdac/phaseA_spotcheck.py` (log `logs/phaseA_spotcheck_20261025T182020.log`)
- **Method**: each number is recomputed from the row-level table, not copied from a summary file.
- **Output**:
  - `results/pdac/phaseA/spotcheck.csv`
  - `figures/pdac/phaseA/phaseA_spotcheck_table.{png,pdf}`
- **Key numbers** (claimed → recomputed; **all 5 match**):
  1. Breast gate 3/9 → 3/9 (`results/pdac/phase1/breast_matched_baseline_paper22.csv`)
  2. CAF→T edge-significant sections 16/16 (marker) and 16/16 (author) → 16/16 and 16/16
  3. Q-significant 4/16 (marker) and 2/16 (author); 0/16 after BH → 4/16 and 2/16; 0/16 and 0/16 after BH within section
  4. Background-flagged axes 24/31 → 24/31. The recount from the per-section × gene table equals the stored per-axis flag column.
  5. Spearman(edge z, ligrec mean), crop: 0.02 → 0.022 (n = 2,821 rows with finite values)
- **Decision**: no correction to `docs/PDAC_SUMMARY.md` is needed.

## 2026-10-25 18:22 — Phase A.3: p-value conventions side by side

- **Command**: `python src/pdac/phaseA_pvalue_conventions.py` (log `logs/phaseA_pvalue_conventions_20261025T182230.log`)
- **Run note**: the first run wrote the PDAC long table as a 101.6 MB CSV, which is over GitHub's 100 MB file limit. That file and the run's other outputs, log and MANIFEST rows were all from minutes earlier and uncommitted, so they were deleted. The script was changed to write `pvalue_conventions_pdac.csv.gz` (11.7 MB) and re-run. The numbers are identical.
- **Recompute vs stored**: p was recomputed from the saved null `.npz` files for breast paper/extra/leiden and for PDAC U5-a marker/Leiden. The maximum absolute difference from the stored `*_p_legacy` and `*_p_plus1` columns is 1.1e-16, with 0 verdict differences. The PDAC 16-section tables have no saved nulls, so their stored columns were used (both were written by the same `ccc_audit.stats`).
- **Key numbers** (`results/pdac/phaseA/pvalue_conventions_flip_summary.csv`, `pvalue_conventions_flipped_rows.csv`, `pvalue_conventions_headline_CAF_T.csv`):
  - **Breast**: 0 uncorrected flips in 292 rows. The Phase 1 gate is 3/9 under both conventions. 3 BH flips in the `extra` family, all significant only under legacy:
    - CD86→CTLA4 fibroblast→immune, edge: q 0.032 → 0.064
    - PTN→SDC4 epithelial→epithelial, edge: q 0 → 0.064
    - CD86→CTLA4 fibroblast→endothelial, compositional: q 0.0498 → 0.057
  - **PDAC**:
    - 9 uncorrected flips in 532,952 PDAC rows [corrected in Phase D from 533,244, the breast + PDAC total]. All are ratio Q with n_valid nulls < 500 and p_legacy 0.048–0.0495.
    - 4,530 BH flips, all significant only under legacy (legacy p = 0 lowers q). [Corrected in Phase D: first written as 4,533, which is the breast + PDAC total from `pvalue_conventions_summary.json`.]
    - CAF→T headline, uncorrected: identical under both conventions.
    - CAF→T headline, BH within section, edge: marker 16 → 13/16, author 16 → 14/16.
    - CAF→T headline, BH within section, Q: marker 0 → 0/16, author 1 → 0/16.
- **Answer**: no uncorrected verdict flips for any breast row or any headline row. The flips that do occur are BH-level, and the legacy convention is always the more liberal one.
- **Doc**: added a "p-value convention" section to `docs/RESULTS_SETTINGS_NOTE.md`.
- **Next**: commit Phase A, then Phase B.

## 2026-10-25 18:25 — Phase B.0: protein-guided rule pre-registered

- **Done**: wrote `docs/PDAC_V2_PROTEIN_RULE.md` (channels, per-section GMM thresholds with an Ashman-D acceptance check, median + 3·MAD fallback, gating order, RNA argmax within each compartment, the `protein_mad` sensitivity variant, validation plan and CCC plan). It was committed alone (5b6089e) **before** any label or CCC result existed.
- **Protein channels present**: `Mean./Max.` for CD3, CD45, PanCK, MembraneStain and DAPI. No other lineage protein channels. MembraneStain and DAPI are not used for gating.
- **Helper**: added `common.mem_status()` (Windows ctypes; system available RAM and process working set). psutil is not in the env, and the env was not modified.

## 2026-10-25 18:27–18:31 — Phase B.1: protein labels and validation — **pre-registered rule FAILED validation**

- **Command**: `python src/pdac/phaseB1_protein_labels.py` (log `logs/phaseB1_protein_labels_20261025T183119.log`)
- **Run note**: the first two runs crashed in the figure code (duplicate U5-a scope; then a pivot shape mismatch). Their logs were deleted while re-running. To keep the provenance complete, all per-section outputs were then deleted and the script was re-run from scratch. Per-cell label files are byte-identical (md5 of the decompressed content) to the earlier run.
- **Outputs**:
  - `results/pdac/phaseB/sections/<section>_labels.csv.gz` and `_thresholds.csv` (16 each)
  - `protein_threshold_diagnostics.csv`, `label_counts_by_section.csv`, `validation_precision_recall.csv`, `confusion_{protein_gmm,protein_mad,marker}_vs_author_{U5-a,all_16_sections}.csv`, `protein_labels_summary.json`
  - `figures/pdac/phaseB/phaseB1_protein_gating_validation.{png,pdf}`
- **Memory**: peak process working set 2.21 GB; minimum system available 5.86 GB.
- **Key numbers**:
  - GMM acceptance (Ashman D ≥ 2): CD3 3/16 sections, CD45 2/16, PanCK 15/16. The other section × channel combinations used the MAD fallback.
  - U5-a: CD3 Ashman D = 0.43; T_cell count is 1,467 under `protein_gmm` and 0 under `protein_mad`.
  - **T_cell vs author T (CD8+ T, CD4+ T, Treg), precision / recall**:

    | Labels | U5-a | Pooled 16 sections |
    |---|---|---|
    | `protein_gmm` | 0.022 / 0.036 (32 of 1,467 predicted; 879 author T) | 0.077 / 0.049 (563 / 7,273; 11,557 author T) |
    | `protein_mad` | 0 cells predicted | 0.031 / 0.005 |
    | v1 marker | 0.276 / 0.297 (261 / 945) | 0.353 / 0.318 (3,676 / 10,427) |

  - Diagnostic (computed inline, not saved as a file; recomputable from the per-cell tables): per-section AUC for author T vs all other cells was CD3 0.451–0.921 and CD45 0.502–0.887. U5-a: CD3 0.703, CD45 0.549. Author malignant vs rest, PanCK: 0.744–0.993.
- **Where the new labels are poor**: T cells (above). Also, pooled per-class precision for `protein_gmm` is below marker for every class with predictions (figure, bottom right; acinar has no precision bar for either annotation). PanCK gating of tumour works, but the CD45 gate is weak, so myeloid, B_plasma and mast precision drop.
- **Decision**: the pre-registered rule **does not produce validated T-cell labels**. Its CCC re-run is still carried out as pre-registered and reported as "labels failed validation". A post-hoc consensus annotation was added (next entry). It is labelled post-hoc everywhere and was committed (161418a) before any v2 CCC run.

## 2026-10-25 18:34 — Phase B.2: POST-HOC `consensus_T` labels and protein triangulation

- **Command**: `python src/pdac/phaseB2_consensus_labels.py` (log `logs/phaseB2_consensus_labels_20261025T183400.log`)
- **Rule**: `T_cell` = author T **and** v1 marker T_cell. Author T without marker support → `T_cell_unconfirmed`. Every other cell keeps its author-mapped label.
- **Counts**:

  | Group | All 16 sections | U5-a |
  |---|---|---|
  | consensus T | 3,676 | 261 |
  | author-only T | 7,881 | 618 |
  | marker-only T | 6,751 | 684 |

- **Independent protein check** (median over sections, AUC vs all non-T cells; `consensus_T_protein_triangulation_median.csv`):

  | Group | CD3 AUC | CD45 AUC |
  |---|---|---|
  | consensus T | 0.598 | 0.769 |
  | author-only T | 0.580 | 0.739 |
  | marker-only T | 0.588 | 0.684 |

  Consensus T cells have the highest protein support, but only slightly. **The protein channels cannot independently confirm any T-cell set in this dataset.** So "validated" here means two RNA annotations agree; it does not mean protein-confirmed.
- **Figure**: `figures/pdac/phaseB/phaseB2_consensus_T_protein_triangulation.{png,pdf}`

## 2026-10-25 18:35–18:47 — Phase B.3: CCC re-run under v2 annotations (16 sections, full catalogue)

- **Commands**:
  - `python src/pdac/phaseB3_ccc_rerun.py T6-b`: smoke test, 16.5 s (log `logs/phaseB3_ccc_rerun_20261025T183455.log`)
  - `python src/pdac/phaseB3_ccc_rerun.py`: all 16 sections (log `logs/phaseB3_ccc_rerun_20261025T183529.log`). T6-b was re-used through resume; the aggregate files were rewritten by this run.
- **Settings**: C/E/M/Q; 500 perms (seed 42, `generate_label_permutations`); 6-NN; within-section permutations; 31 pairs × all sender-receiver combinations; add-one p; BH within section × annotation.
- **Outputs**:
  - `results/pdac/phaseB/ccc_sections/<section>__<annotation>.csv.gz` (48 files, each written immediately)
  - `v2_ccc_6NN_all_rows.csv.gz` (212,970 rows)
  - `v2_headline_CAF_T_rows.csv`
- **Memory**: peak working set 2.11 GB; minimum system available 5.83 GB.
- **Manifest race**: while this run was registering files, the Phase C self-test cleanup rewrote `MANIFEST.csv` from a separate process. The completeness check in Phase D re-registers any rows lost this way.

## 2026-10-25 18:47 — Phase B.4: comparison across annotations — **the key question**

- **Command**: `python src/pdac/phaseB4_compare.py` (log `logs/phaseB4_compare_20261025T184712.log`)
- **Outputs**:
  - `headline_CAF_T_by_annotation_section.csv`, `headline_CAF_T_summary_by_annotation.csv`, `headline_CAF_T_U5-a_by_annotation.csv`
  - `annotation_pairwise_comparison.csv`, `edge_vs_product_mirroring_by_annotation.csv`, `phaseB4_summary.json`
  - `figures/pdac/phaseB/phaseB4_headline_and_annotation_agreement.{png,pdf}`
- **Headline CXCL12→CXCR4 CAF→T_cell, 16 sections** (add-one p < 0.05; BH within section):

  | Annotation | Evaluable | E sig (BH) | Q sig (BH) | Q z > 0 | median z_E | median z_Q | median Q |
  |---|---|---|---|---|---|---|---|
  | marker (v1) | 16 | 16 (13) | 4 (**0**) | 14 | 20.87 | 0.75 | 1.347 |
  | author (v1) | 16 | 16 (14) | 2 (**0**) | 15 | 21.63 | 0.74 | 1.346 |
  | protein_gmm (failed validation) | 14 | 10 (6) | 4 (**0**) | 11 | 5.00 | 1.17 | 1.623 |
  | protein_mad (failed validation) | 13 | 9 (4) | 5 (**0**) | 10 | 4.69 | 1.56 | 1.799 |
  | **consensus_T (post-hoc)** | 16 | 15 (11) | 2 (**0**) | 11 | 14.68 | 0.36 | 1.320 |

  Not evaluable (no `T_cell` class): `protein_gmm` T3-b and U7-a; `protein_mad` T5-a, U5-a and U7-a.
- **U5-a**:

  | Annotation | z_C / z_M / z_E / z_Q | Q | Q p | Q BH q |
  |---|---|---|---|---|
  | marker | 33.2 / 36.1 / 23.4 / 2.81 | 1.62 | 0.008 | 0.108 |
  | author | 31.8 / 35.9 / 23.0 / 2.66 | 1.75 | 0.006 | 0.074 |
  | protein_gmm | 30.4 / 31.3 / 13.8 / 0.99 | 1.28 | 0.166 | 0.732 |
  | consensus_T | 21.6 / 24.6 / 17.5 / 1.94 | 1.96 | 0.038 | 0.314 (448 edges) |

- **Answer to the key question**: **No.** With the best-validated T-cell labels available (`consensus_T`), CAF→T Q is significant uncorrected in 2/16 sections (T1-b p = 0.044; U5-a p = 0.038) and in **0/16 after BH within section**. That is the same as the marker and author labels (0/16 after BH).
  - The edge result persists (15/16; 11/16 after BH).
  - It is carried by the matched product M: consensus z_M > z_E in all 16 sections (z_M 6.1–54.3 vs z_E 0.24–40.0).
  - The observed Q is above 1 in 14/16 consensus sections (median 1.32; below 1 in T4-a and U6-a_2), but that does not exceed the permutation null.
  - The two protein annotations failed validation and give no Q result surviving BH either.
- **Edge score vs cell-type product, full catalogue** (`edge_vs_product_mirroring_by_annotation.csv`):
  - Spearman(z_E, z_M) across evaluated rows: 0.66–0.68 for marker, author and consensus; 0.80–0.81 on the 7 axes without background flags.
  - 77% of edge-significant rows are also M-significant, 35–37% are Q-significant uncorrected, and **≈ 10% are Q-significant after BH** (marker 0.098, author 0.104, consensus 0.096).
- **Annotation agreement** (pooled 16 sections, Spearman of z on shared rows; `annotation_pairwise_comparison.csv`):
  - author vs consensus: E 0.970, Q 0.951 (they differ only in the T class).
  - marker vs author: E 0.655, Q 0.394.
  - author vs protein_gmm: E 0.348, Q 0.188.
  - Q remains the most annotation-sensitive score.
- **Next**: Phase C.

## 2026-10-25 18:36–18:42 — Phase C: independent replication on GSE310352 — **BLOCKED (licence not confirmed; validated labels unreadable on this host)**

- **Fetch attempts (normal fetches only; nothing worked around)**:
  - `https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html` (GEO terms page): **reCAPTCHA**, "Checking your browser - reCAPTCHA". Not bypassed.
  - Europe PMC REST record for PMID 41533516 (200). Title "In situ multi-modal characterization of pancreatic cancer reveals tumor cell identity as a defining factor…"; Cell Reports 45(1), 2026 Jan; doi 10.1016/j.celrep.2025.116827; first authors Lyubetskaya A, Rabe B, Kavran A, Bai Y, Fisher A. Fields `license: "cc by-nc-nd"`, `isOpenAccess: "N"`, `inPMC: "N"`. This licence is the **article's**; it is not a data licence.
  - `https://doi.org/10.1016/j.celrep.2025.116827` → 302 to linkinghub.elsevier.com, which returned only "Redirecting." `https://www.cell.com/cell-reports/fulltext/S2211-1247(25)01599-2`: **HTTP 403**. The paper's "Data and code availability" section could not be read.
  - `https://www.ncbi.nlm.nih.gov/home/about/policies/` (200). Its "Molecular Data Usage" section covers gene-expression data. Quotes: "NCBI itself places no restrictions on the use or distribution of the data contained therein"; "Nor do we accept data when the submitter has requested restrictions on reuse or redistribution"; "some submitters of the original data ... may claim patent, copyright, or other intellectual property rights in all or a portion of the data"; "NCBI cannot provide comment or unrestricted permission concerning the use, copying, or distribution". This is a general NCBI policy. It is **not a dataset licence**, and it explicitly leaves submitter rights open.
  - `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310352/suppl/filelist.txt` (200). Per slide: `cell_by_gene_counts.csv.gz` (e.g. slide1 13,671,200 B), `cell_metadata.csv.gz` (slide1 5,260,222 B), `raw_core.tar.gz` (157,983,978–676,615,407 B) and `seurat.RDS` (121,711,464–526,674,009 B). One range request: the first 40 KB of `GSM9294399_slide1_cell_metadata.csv.gz` (HTTP 206), saved to `data/pdac/geo_gse310352/probe/` (git-ignored).
- **Format: CONFIRMED** (v1 and this check):
  - counts CSV = `fov`, `cell_ID` + 970 genes + 19 NegPrb
  - metadata CSV = `fov`, `cell_ID`, `Area`, `AspectRatio`, `CenterX/Y_local_px`, `CenterX/Y_global_px`, `Width`, `Height`, `Mean./Max.` for MembraneStain, PanCK, CD45, CD3 and DAPI
  - **No cell-type labels in the CSVs.** Author labels exist only inside the `seurat.RDS` objects (column names UNVERIFIED).
- **Licence: NOT CONFIRMED.** Per the v2 rule ("If you cannot confirm the license … log BLOCKED"), nothing beyond the 40 KB header probe was downloaded and no analysis was run on GSE310352.
- **Second blocker: validated labels.** Phase B shows that the CosMx protein channels (same panel here) do not give valid T-cell labels. A validated T-cell receiver in GSE310352 therefore needs the authors' labels from `seurat.RDS`. This host has no R (`Rscript` not on PATH, no `C:/Program Files/R`) and no Python RDS reader (`rdata` and `pyreadr` not installed; neither reads Seurat S4 objects reliably). The project env was not modified.
- **Prepared so it can run as soon as it is unblocked**:
  - `pdac_loaders.load_geo_gse310352` (slide → GSM map taken from `filelist.txt`)
  - `src/pdac/phaseC_gse310352_audit.py`, gated by `--licence-confirmed`. Without the flag it prints BLOCKED and exits with code 2 (tested). It downloads only the 16 small CSVs (124,527,604 B), runs one slide at a time with resumable outputs, and covers marker labels, the pre-registered protein labels and optional user-supplied labels (`--labels`).
  - **Self-test on SYNTHETIC random data** (8 fake slides × 1,500 cells, real header layout) ran end to end, exit 0. All synthetic outputs and their MANIFEST rows were deleted. The run log was kept as `logs/phaseC_gse310352_audit_20261025T184035_SYNTHETIC_SELFTEST.log`. No real GSE310352 numbers exist anywhere in the repository.
- **What I need from you (any one of 1a–1c, plus 2 for validated labels)**:
  1. Licence or terms, one of:
     - (a) paste the "Data and code availability" section (and any data-licence sentence) of Lyubetskaya et al., Cell Reports 2026, doi 10.1016/j.celrep.2025.116827;
     - (b) paste the text of https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html;
     - (c) state that NCBI's Molecular Data Usage policy (quoted above) is sufficient for this non-redistributing audit use.
  2. Validated labels: either install R + SeuratObject (then I export the metadata myself), or run this in R and give me the CSV:
     ```r
     library(SeuratObject)
     for (s in c("slide1","slide2","slide3","slide4","slide5","slide5b2","slide5b3","slide6")) {
       f <- list.files(".", pattern = paste0("_", s, "_seurat.RDS$")); o <- readRDS(f)
       write.csv(cbind(slide = s, o@meta.data), paste0(s, "_seurat_metadata.csv"), row.names = TRUE)
     }
     ```
     Then tell me which column holds the cell type, and which columns hold `fov` and `cell_ID` (or how the row names encode them).
  3. Optional: the µm/px for this CosMx run. It is only needed for radius analyses; the 6-NN audit is scale-free.
- **Next**: Phase D (synthesis) with Phase C marked BLOCKED.

## 2026-10-25 18:44–19:00 — Push of `pdac-v2`: **BLOCKED** for commits after 161418a

- `git push origin pdac-v2` (Phase B + C commits, about 86 MB of new results) failed:
  - `error: RPC failed; curl 56 schannel: failed to read data from server: SEC_E_MESSAGE_ALTERED (0x8009030f)`
  - `send-pack: unexpected disconnect while reading sideband packet`
- Pushing commit by commit with `-c http.postBuffer=524288000`: 5b6089e and 161418a went through, and the remote is at 161418a. 5166445 (Phase B results) failed again with the same schannel error.
- SSH: `git@github.com: Permission denied (publickey)` (no GitHub key on this host).
- One attempt with `-c http.sslBackend=openssl -c http.sslVerify=true` (certificate verification **on**): "remote end hung up".
- **Decision**: the same failure happened twice. Per the stop rule, push is BLOCKED and no further attempts were made. TLS was not disabled. All commits are on the local `pdac-v2` branch. Nothing was pushed to `main`, and nothing was force-pushed.
- **To finish from a stable network**: `git push origin pdac-v2`

## 2026-10-25 18:52–19:00 — Phase D: synthesis

- **Commands**:
  - `python src/pdac/phaseD_manifest_check.py` (logs `logs/phaseD_manifest_check_*.log`)
  - `python src/pdac/phaseD_synthesis.py` (log `logs/phaseD_synthesis_*.log`). The first run crashed: the `axes` column collided with `DataFrame.axes`. Its log was deleted, the script was fixed, and it was re-run.
- **Outputs**:
  - `results/pdac/phaseD/v2_key_numbers.csv`: every number quoted in the v2 summary, with its source file.
  - `results/pdac/phaseD/manifest_check.csv` (final run: 0 issues). The first run's file was overwritten by the second run under the same name; it was restored from commit c95d2fa as `manifest_check_first_run.csv`.
  - `figures/pdac/phaseD/phaseD_v2_synthesis.{png,pdf}`
  - `docs/PDAC_SUMMARY_v2.md`
  - README "PDAC extension" section: v2 links added.
  - `docs/PDAC_SUMMARY.md`: a pointer to v2 at the top. The v1 text is otherwise unchanged.
- **Manifest check**:
  - 2 files had lost their rows to the concurrent manifest write (`phaseB/ccc_sections/U2-b__consensus_T.csv.gz` and the Phase C synthetic self-test log). Both were re-registered.
  - 1 stale row (deleted log `phaseB1_protein_labels_20261025T183049.log`) was removed.
- **Correction found while building the key-numbers table**: the Phase A statements "4,533 PDAC BH flips" and "533,244 PDAC rows" were breast + PDAC totals. Corrected in place to 4,530 and 532,952, in this log and in `docs/RESULTS_SETTINGS_NOTE.md`. The breast share is 3 BH flips and 292 rows.
- **Final status**: A DONE; B DONE (negative result); C BLOCKED (licence, labels); D DONE; push BLOCKED after 161418a.

## 2026-10-25 19:10 — Phase D2: tumour-level view, and summary restructured (the remaining part of the brief)

- **Brief**: the user supplied the cut-off end of the v2 brief. The summary must lead with "Does a proximity-specific signal exist at tumour level, or is the edge score mainly mirroring cell-type products?", followed by results, limitations, BLOCKED items and at most 5 next steps, linked from the README.
- **Why a new analysis**: all v2 results so far were per section. Answering "at tumour level" needs tumours as the unit.
- **Command**: `python src/pdac/phaseD2_tumour_level.py` (log `logs/phaseD2_tumour_level_*.log`)
- **Method** (declared in the script docstring before running):
  - 16 sections → 13 tumours (donor id without `_2nd`/`_B`; an assumption consistent with the stated 13 tumours).
  - Tumour p = Fisher combination of the section add-one p; BH across tumours.
  - Exploratory one-sided Wilcoxon of tumour-mean z > 0.
- **Outputs**:
  - `results/pdac/phaseD/tumour_level_headline_CAF_T.csv`, `tumour_level_headline_summary.csv`, `tumour_level_summary.json`
  - `figures/pdac/phaseD/phaseD2_tumour_level_headline.{png,pdf}`
- **Key numbers**:
  - CAF→T Q, tumours significant after BH: **0/13** for marker, author and consensus_T; 0/11 for protein_gmm; 0/10 for protein_mad.
  - Uncorrected: marker 3, author 2, consensus 2.
  - E and M: 13/13 tumours significant after BH for marker, author and consensus.
  - Tumour-mean z_Q > 0 and Wilcoxon p:
    - marker: 12/13, p = 0.0085
    - author: 12/13, p = 0.0004
    - consensus: 10/13, p = 0.016 (smallest consensus BH q = 0.285)
  - Mean observed Q > 1: marker 12/13, author 13/13, consensus 11/13.
- **Interpretation (audit level)**: no tumour shows a proximity-specific excess after correction. The direction is consistent across tumours, which is exploratory and confounded with tissue-scale co-zonation under the label-permutation null. The edge score mainly mirrors the cell-type product.
- **Docs**:
  - `docs/PDAC_SUMMARY_v2.md` rewritten to the requested structure (lead question → results → limitations → BLOCKED → 5 next steps).
  - README v2 entry updated with the lead question and answer.

# PDAC v3 (branch `thesis-v3`) — contact vs zonation, calibration, effect sizes

## 2026-11-15 15:40 — Setup

- **Brief**: v3 brief (Phases 1–6; Phases 4 and 5 gated). Autonomous run; questions are resolved by the most reasonable reading and logged here.
- **Branch**: `thesis-v3` created from `origin/pdac-v2` (1173e52); `origin/main` (ed86917, GSE277783 Visium EDA in `eda_gse277783/`) merged in. The merge was a fast-forward (pdac-v2 is an ancestor of main), so no conflicts arose. Upstream tracking removed so nothing points at pdac-v2.
- **Environment**: conda env `thesis-spatial-ccc` (`D:/Anaconda/envs/thesis-spatial-ccc/python.exe`, Python 3.11.10, numpy 1.26.4, pandas 2.2.3, scipy 1.14.1, scanpy 1.10.3, statsmodels 0.15.0), the same as v1/v2. The base Anaconda Python is unusable (NumPy 2.5.2 breaks compiled modules) and was not used.
- **Gates checked at setup** (re-checked by scripts later):
  - Phase 4: `LICENSE_CONFIRMED_GSE310352.txt` absent → will be BLOCKED (licence unconfirmed).
  - Phase 5: no GEO SOFT/series-matrix file for GSE277783 under `data/`; `eda_gse277783/EDA_REPORT.md` records the organ site as UNKNOWN for 49/57 slides → will be BLOCKED. Gene lists were still pre-declared in `docs/v3_plans/phase5_visium_programs_gse277783.md`.
- **Plans committed before computing**: `docs/v3_plans/phase1_contact_vs_zonation.md` (a73da4d), `phase2_calibration_simulation.md` (580987b), `phase3_effect_sizes_mixed_model.md`, `phase4_replication_gse310352.md`, `phase5_visium_programs_gse277783.md` (aef6b71).
- **Conventions**: add-one p = (k+1)/(n+1) everywhere (stated in a `p_value_convention` column of every new table); seed 42; 500 permutations; relative paths; audit/diagnostic framing only.
- **Shared code**: `src/pdac/v3_nulls.py` (block and distance-stratified permutation generators; section loader re-using the v2 label files `results/pdac/phaseB/sections/*`).

## 2026-11-15 15:50 — Phase 1: implementation checks (before the full run)

- `python src/pdac/v3_phase1_constrained_nulls.py --sections T6-b --no-aggregate` (test on the smallest section, 58.7 s).
- **Reproduction check**: the `global` null re-creates the Phase B3 file `results/pdac/phaseB/ccc_sections/T6-b__consensus_T.csv.gz` exactly (max |diff| = 0 for E obs, E z, Q obs, Q z, Q add-one p over all 5,239 rows). The test output was kept (it is the first section of the full run; the script resumes from existing files).
- **Tile occupancy, T6-b** (from the run log): 50 um → 2,140 tiles, median 5 cells/tile, 264 singleton tiles; 100 um → 635 tiles, median 16; 200 um → 205 tiles, median 51. Distance null: 8 of 12 bands used, smallest band 6 cells.
- **Phase 2 fast evaluator check** (ad hoc, not saved): the per-row bincount evaluation (`PermStore.null_scores`) gives the same z for C/M/E/Q and the same add-one p as the Phase 1 files for T6-b under `global`, `block50` and `dist` (identical to 6 decimals).
- **Housekeeping**: a 3-simulation test of `v3_phase2_calibration.py` on T6-b wrote `results/pdac/v3_phase2/sections/sim_rows_T6-b.csv.gz`, `T6-b_setup.json` and a run log; these were deleted and their 3 manifest rows removed so that the real 100-simulation run is the only Phase 2 output. The manifest was edited while Phase 1 was running; the manifest check at the end of v3 will recover any row lost to the concurrent write (as happened in v2). No further manifest-writing jobs were started while Phase 1 ran.
- **Full run started 15:51**: `python src/pdac/v3_phase1_constrained_nulls.py` (all 16 sections, 3 annotations, 5 nulls; stdout in `logs/v3_phase1_stdout.txt`).

## 2026-11-15 15:51–16:40 — Phase 1: contact vs zonation — DONE

- **Command**: `python src/pdac/v3_phase1_constrained_nulls.py` (log `logs/v3_phase1_constrained_nulls_20261115T155116.log`; stdout `logs/v3_phase1_stdout.txt`). 16 sections x 3 annotations x 5 nulls = 240 section files under `results/pdac/v3_phase1/sections/`; 0 runs skipped (every section has >= 50 tumour and >= 50 non-tumour cells under every annotation, so the distance null always ran). Peak process RSS 1.7–2.4 GB; per-section wall time 59–420 s.
- **Outputs** (all registered): `headline_by_null.csv`, `headline_summary_by_null.csv`, `tumour_level_headline_by_null.csv`, `tumour_level_summary_by_null.csv`, `nonbg_catalogue_by_null.csv.gz` (227,920 rows), `nonbg_summary_by_null.csv`, `skipped_runs.csv` (empty), `phase1_summary.json`, per-section `*_null_info.json` (tile and band occupancy); figure `figures/pdac/v3_phase1/phase1_Qz_vs_null.{png,pdf}`.
- **Key numbers, headline CXCL12→CXCR4 CAF→T, consensus_T** (`headline_summary_by_null.csv`; 16 sections; add-one p):
  - median z_Q: global 0.36 → block200 0.84 → block100 0.88 → block50 1.35 → dist 0.38
  - sections with Q p<0.05 (BH within section): global 2 (0), block200 3 (1), block100 6 (0), block50 6 (2), dist 3 (0)
  - median null mean of Q: global 1.067, block200 1.009, block100 0.954, block50 1.056, dist 1.136; median null sd: 0.433, 0.293, 0.257, 0.236, 0.331. Observed median Q 1.320 is unchanged by construction.
  - median z_E drops from 14.7 (global) to 8.3 (block50): the constrained nulls remove part of the composition/zonation component of E, as intended.
- **Tumour level, Q, consensus_T** (`tumour_level_summary_by_null.csv`; 13 tumours; Fisher + BH across tumours):
  - median tumour-mean z_Q: global 0.43, block200 0.86, block100 0.81, block50 1.30, dist 0.58
  - tumours with z_Q>0: 10, 11, 12, 11, 10 of 13; tumours significant after BH: 0, 2, 0, 3, 0; Wilcoxon p (exploratory): 0.016, 0.0006, 0.0002, 0.0009, 0.020
  - marker: median z_Q 0.90 → 1.00 → 1.37 → 1.23 → 1.08; BH-significant tumours 0, 1, 4, 3, 0. author_mapped: 0.78 → 0.97 → 1.70 → 1.49 → 0.82; BH 0, 3, 3, 4, 0.
- **Answer to the key question (rule declared in the plan)**: the tumour-level positive trend does **not** shrink toward zero as the null becomes more spatially constrained; `shrinks_toward_zero = False` for all three annotations (`phase1_summary.json`). The distance-stratified null leaves it about where the global null put it. The mechanism is visible in the null moments: the null mean of Q stays near 1 under every null (so zonation at the 50–200 um scale does not reproduce the observed Q ≈ 1.32), while the null sd falls, so z rises. Interpretation stays at audit level: the Q excess for this axis is not explained by label zonation at >= 50 um; whether the audit has the power to see a real contact effect, and whether zonation of *expression* (rather than labels) can fake it, is the Phase 2 question.
- **Full non-background catalogue** (7 axes x all sender→receiver, `nonbg_summary_by_null.csv`): the fraction of rows with Q p<0.05 is 0.138 (global) → 0.131 → 0.126 → 0.122 (block50) → 0.110 (dist) for consensus_T, far above 0.05 under every null; BH-significant Q rows fall from 1,234 (global) to 565 (block50). Spearman of z_Q against the global null z_Q: 0.79 (block200), 0.71 (block100), 0.57 (block50), 0.86 (dist); for z_E 0.90 / 0.86 / 0.81 / 0.81. The constrained nulls therefore re-rank the catalogue substantially at the 50 um scale but do not remove the catalogue-wide excess of small p-values.
- **Caveat logged**: block permutation within 50 um tiles (median 5 cells per tile, 12% singleton tiles in T6-b) still exchanges labels among immediate neighbours, so it is a contact-breaking null; it is not a test of whether the labels themselves are spatially accurate.
- **Next**: Phase 2 (calibration), Phase 3 (effect sizes) and the Phase 4/5 gates were started as one sequential job at 16:41 so that only one process writes the manifest at a time.

## 2026-11-15 16:41–16:58 — Phase 2: calibration by simulation — DONE (one crash, fixed)

- **Command**: `python src/pdac/v3_phase2_calibration.py` (sequential job with Phase 3 and the gates). Sections by the fixed rule: T6-b (smallest, 12,473 cells), T2-b (lower median, rank 8 of 16, 43,658), U5-a (largest, 74,793). 100 simulations x 500 permutations per setting, 5 nulls, scores C/M/E/Q; add-one p.
- **Setup numbers** (`results/pdac/v3_phase2/sections/<section>_setup.json`): T cells / CAF→T edges / proximity targets = T6-b 26 / 83 / 13; T2-b 23 / 74 / 10; U5-a 261 / 448 / 85. Mean CXCR4 count in targets 1.69, 0.70, 1.55. The renormalisation check against `adata.raw` passed (max |diff| 3e-7).
- **Crash**: the first run wrote all simulation rows and both tables, then failed in the summary JSON (`AttributeError: 'Pandas' object has no attribute '_14'`; positional itertuples fields). The log of that run is kept and registered (`logs/v3_phase2_calibration_20261115T164130.log`). The bug (also present in the unrun Phase 6 script) was fixed to use named columns and the aggregation re-run with `--aggregate-only` (no simulation was recomputed; `logs/v3_phase2_stdout_rerun.txt`). This counts as one failure, not two.
- **Outputs**: `sim_rows.csv.gz` (per-section files under `sections/`), `calibration_table.csv` (720 rows: per section and pooled), `calibration_deliverable.csv` (20 rows: score x null), `calibration_summary.json`; figure `figures/pdac/v3_phase2/phase2_calibration.{png,pdf}`.
- **Results, Q (the proximity-specific score)** — detection rate = fraction of 100 sims with add-one p < 0.05 (`calibration_table.csv`):
  - T6-b and T2-b (23–26 T cells): 0.00 for every null and every planted fold change up to f = 2, in both scenarios. Zero power and zero false positives: the row is simply not testable at this cell count.
  - U5-a (261 T cells): power at f = 1.1 / 1.25 / 1.5 / 2 — global 0.12 / 0.34 / 0.78 / 1.00; block200 0.00 / 0.00 / 0.02 / 0.49; block100 0.18 / 0.47 / 0.85 / 1.00; block50 0.01 / 0.11 / 0.52 / 1.00; dist 0.39 / 0.73 / 0.92 / 1.00. Control (f = 1): 0.02 / 0.00 / 0.06 / 0.00 / 0.18.
  - U5-a, zonation-only false-positive rate at f = 1.1 / 1.25 / 1.5 / 2: global 0.42 / 0.44 / 0.56 / 0.54; block200 0.10 / 0.06 / 0.06 / 0.08; block100 0.56 / 0.60 / 0.73 / 0.74; block50 0.35 / 0.32 / 0.46 / 0.43; dist 0.74 / 0.76 / 0.88 / 0.86.
  - Pooled over the 3 sections (`calibration_deliverable.csv`): no score/null combination reaches power 0.8 at any f (max pooled Q power 0.33, i.e. U5-a alone). Fooled-by-zonation flag (pooled FPR > 0.10 at f = 1.5): True for Q under global (0.19), block100 (0.24), block50 (0.15), dist (0.29); False only for block200 (0.02).
- **Results, C, M, E**: detection rate 1.00 for every setting including the f = 1 control and every zonation setting (E: 0.75–0.96 at f = 1 depending on the null). These scores are saturated by the real composition signal on these sections; they cannot distinguish a planted contact effect from zonation or from nothing.
- **Reading (audit level)**: the Q audit has no power in sections with ~25 T cells; in the largest section it detects a 1.5x contact effect with the global/block100/dist nulls, but those nulls also call a 200 um expression zonation positive in 42–88% of simulations. The only null that is not fooled (block200, matching the planted zonation scale) detects only a 2x effect, and only half the time. Note that the planted zonation scale (200 um) matches block200 by design; zonation at other scales was not simulated.
- The planted perturbation is a computational stand-in for a stimulation experiment; nothing experimental is claimed.

## 2026-11-15 16:54 — Phase 3: effect sizes and tumour-level model — DONE

- **Command**: `python src/pdac/v3_phase3_effect_sizes.py` (log `logs/v3_phase3_effect_sizes_*.log`). 100 um tiles, B = 1000, seed 42.
- **Outputs**: `results/pdac/v3_phase3/Q_block_bootstrap_by_section.csv`, `Q_block_bootstrap_nonbg_axes.csv`, `Q_block_bootstrap_summary.csv`, `mixed_model_logQ.csv`, `z_vs_edge_count.csv`, `z_vs_edge_count_summary.csv`, `phase3_summary.json`; figure `figures/pdac/v3_phase3/phase3_effect_sizes.{png,pdf}`.
- **Block-bootstrap CI, headline Q** (`Q_block_bootstrap_summary.csv`): consensus_T — CI above 1 in 4/16 sections, below 1 in 1/16, median Q 1.320, median CI width 0.767, median 448 edges. marker 8/16 above (1 below), median width 0.480. author_mapped 7/16 above (0 below), width 0.531. Non-background axes CAF→T (consensus_T): CXCL12-CXCR4 is the only axis with CIs above 1 (4/16); MIF-CD74 is below 1 in 6/16, SPP1-CD44 in 3/16.
- **Mixed model log Q ~ 1 + (1|tumour)** (`mixed_model_logQ.csv`; exploratory; 16 sections, 13 tumours, 10 singleton tumours): geometric-mean Q = 1.242 [1.013, 1.522] consensus_T (tumour variance 0.007; statsmodels warning "MLE may be on the boundary"), 1.315 [1.176, 1.471] marker, 1.348 [1.219, 1.491] author_mapped. All converged. The section-level t-interval is nearly identical (consensus_T 1.236 [1.007, 1.517]), so the tumour grouping changes little because most tumours have one section.
- **z vs edge count** (`z_vs_edge_count_summary.csv`, 7 non-bg axes, all rows, global null): Spearman(z_E, n_edges) is about −0.11 but Spearman(|z_E|, n_edges) is 0.52–0.56, i.e. the magnitude of z_E grows with edge count; after dividing by sqrt(n_edges) the signed correlation is 0.19–0.26. For the headline, Spearman(z_E, n_edges) across the 16 sections is 0.71 (consensus_T), while z_Q shows no edge-count dependence (−0.08). Normalised headline z_E: median 14.68 → 0.69 per sqrt(edge).

## 2026-11-15 16:55 — Phases 4 and 5: BLOCKED (gates closed)

- `python src/pdac/v3_phase4_gate.py` → `results/pdac/v3_phase4/phase4_gate_status.json`: **BLOCKED**, licence unconfirmed (`LICENSE_CONFIRMED_GSE310352.txt` absent). `data/GSE310352_RAW.tar` (6.33 GB) is on disk; nothing was extracted or analysed. Also outstanding: um/px unverified, no author labels, patient overlap with the Mendeley set unverified.
- `python src/pdac/v3_phase5_gate.py` → `results/pdac/v3_phase5/phase5_gate_status.json`: **BLOCKED**, no GEO SOFT/series-matrix file for GSE277783 under `data/`; organ site unknown for 49/57 slides (EDA report). The program gene lists stay pre-declared in `docs/v3_plans/phase5_visium_programs_gse277783.md` for a later run.

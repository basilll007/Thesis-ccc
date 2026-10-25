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

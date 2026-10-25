# PDAC Extension — Progress Log (branch `pdac-extension`)

Append-only log. Every number here comes from a file listed in `results/pdac/MANIFEST.csv`. Interpretation is limited to audit and diagnostic statements. Nothing here supports a biological or clinical claim.

## Task checklist

- [x] Phase 0: Repo hygiene (PDCD1 caption, F:\ paths, settings note, "replication" renamed)
- [x] Phase 1: Matched dyadic baseline (breast) — gate: matched baseline does NOT explain most C/E disagreements (3/9), but it carries all edge-significant calls (0/8 ratio-significant)
- [x] Phase 2: Find and vet a PDAC spatial dataset — selected Mendeley PDAC-SMI (CosMx, CC BY 4.0, 37/39 genes)
- [ ] Phase 3: Run the audit on PDAC
- [ ] Phase 4: Benchmark against COMMOT / SpatialDM / Squidpy ligrec
- [ ] Phase 5: Synthesis (PDAC_SUMMARY.md, README section)

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

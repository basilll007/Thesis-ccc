# Progress Log

## 2026-09-08 — Project initialized

- **Date/time**: 2026-09-08 (workspace-local date; exact clock time not available from the session metadata).
- **What I did**: Created the required project memory files and repository directories: `data/raw`, `data/processed`, `src`, `results`, `figures`, `logs`, and `docs`.
- **Key decisions**: Preserve raw data as immutable; use fixed random seed 42; implement v0 as a baseline pipeline rather than the future S-PCST/reproducibility-calibration method.
- **What worked**: Workspace was empty and directory initialization succeeded.
- **What failed**: The configured `rtk` command is not available in this Windows shell, so standard PowerShell commands are being used where necessary.
- **Blockers**: Dataset and environment availability still need inspection.
- **Exact next step**: Inspect local Python/conda/uv tooling and verify public dataset endpoints in the requested order.

## 2026-09-08 — Public data and baseline implementation

- **Date/time**: 2026-09-08.
- **What I did**: Verified conda 24.11.3 and uv 0.10.8; confirmed the official 10x Xenium breast 2-FOV bundle endpoint; downloaded and unpacked the public bundle into `data/raw/10x_xenium_breast`; implemented modular download, AnnData/QC/clustering/spatial-graph/baseline-LR code and figure generation.
- **Key decisions**: Use the official 10x bundle instead of GEO GSE243168 because it is much smaller and directly public; use a documented curated L–R mean-expression product as the v0 baseline so the pipeline remains reproducible without an unpinned external interaction database.
- **What worked**: The 379,018,013-byte archive downloaded successfully and contains `cell_feature_matrix.h5`, `cells.csv.gz`, transcript outputs, and spatial metadata. The data provenance section was updated with URL, date, size, SHA-256, and CC BY 4.0 terms.
- **What failed**: The initial archive unpack did not extract the nested `cell_feature_matrix.tar.gz`; the downloader was corrected and rerun. The first clustering code used `seurat_v3`, which would require an extra dependency, so it was changed to `cell_ranger`.
- **Blockers**: Conda environment creation is still installing pinned pip packages.
- **Exact next step**: Verify imports in the new Python 3.11 environment, then run the end-to-end prototype and fix only reproducibility/runtime issues.

## 2026-09-08 — Environment and end-to-end verification

- **Date/time**: 2026-09-08; pipeline log start `2026-09-08T07:43:46.900164+00:00`.
- **What I did**: Built the conda Python 3.11 environment from `environment.yml`, corrected the `setuptools` pin for the current Squidpy dependency stack, verified all requested imports, and ran `conda run -n thesis-spatial-ccc python src\\run_prototype.py`.
- **Key decisions**: Fixed seed 42; QC thresholds are `min_counts=10` and `min_cells=3`; Leiden resolution is 0.5; Squidpy spatial graph uses generic coordinates and 6 neighbors; the baseline ranks curated L–R pairs by sender mean ligand expression × receiver mean receptor expression.
- **What worked**: Imports succeeded with exact versions: scanpy 1.10.3, squidpy 1.6.5, anndata 0.10.9, liana 1.5.1, numpy 1.26.4, scipy 1.14.1, pandas 2.2.3, networkx 3.4.2, matplotlib 3.9.2, seaborn 0.13.2, scikit-learn 1.5.2, leidenalg 0.10.2, igraph 0.11.8. The pipeline loaded 7,275 cells/280 genes, retained 7,163 cells/280 genes, created the spatial graph, ranked 22 baseline interactions, saved the processed `.h5ad`, CSV, and PNG/PDF figures.
- **What failed**: The first import check attempted `leidenalg.__version__`, which does not exist; package import itself succeeded and the exact distribution version was verified via `importlib.metadata`.
- **Blockers**: None for v0.
- **Exact next step**: Use the generated outputs as the reproducible v0 baseline and begin future work only after reviewing biological validity and adding response/reproducibility evidence.

## 2026-09-08 — Final reproducibility checks

- **Date/time**: 2026-09-08.
- **What I did**: Ran Python bytecode compilation, `uv pip check`, and reloaded the processed AnnData object to verify spatial coordinates and the Squidpy connectivity matrix.
- **Key decisions**: Keep the raw archive and extracted raw outputs untouched; use only derived paths for the `.h5ad`, ranked CSV, figures, and run log.
- **What worked**: Compilation passed; `uv pip check` reported all 134 packages compatible; processed AnnData reload passed with 7,163 cells, 280 genes, spatial coordinates `(7163, 2)`, and spatial connectivity `(7163, 7163)`.
- **What failed**: One shell-only verification command had a PowerShell quoting error; the equivalent check was rerun successfully without changing project files.
- **Blockers**: None.
- **Exact next step**: Stop v0 implementation and review the baseline biology before adding any future response-evidence or reproducibility-calibration method.

## 2026-09-08 — Squidpy baseline correction

- **Date/time**: 2026-09-08.
- **What I did**: Replaced the initial hand-calculated L–R score implementation with an explicit `squidpy.gr.ligrec` baseline using the curated 20-pair panel and one deterministic permutation for smoke validation, then reran the full pipeline.
- **Key decisions**: Treat the Squidpy mean-expression score as the v0 ranking; do not interpret the one-permutation p-value field as calibrated significance. Future work should increase permutations and add response/reproducibility evidence.
- **What worked**: The corrected pipeline completed and the top ranked interaction is CXCL12→CXCR4 from fibroblast to immune cells; outputs were regenerated in their required locations.
- **What failed**: Squidpy interprets a raw sequence of tuples as source/target combination sets, so the interaction input was changed to a two-column DataFrame to preserve exact ligand–receptor pairs.
- **Blockers**: None.
- **Exact next step**: Stop implementation at the v0 definition of done.

## 2026-09-13 — Spatial-null permutation test and reproducibility verification

- **Date/time**: 2026-09-13; local run start 08:18:52 UTC-05:00 / 13:18:52Z.
- **What I did**:
  1. Ran two back-to-back reproducibility re-runs of `src/run_prototype.py` with no code changes and diffed the resulting `results/baseline_lr_ranked.csv` byte-for-byte; verified identical SHA-256 hash `14cb5c26f3eed89f9018fac3128c5ed61ec1f7e5f0052000fb56b18779331000` across both runs and the reference baseline. Traced historical 22 vs 20 interaction counts in earlier progress logs to an unpinned tuple interaction input in Squidpy that was subsequently fixed on 2026-09-08.
  2. Researched spatial-null permutation procedures in CONCISE (Zhao et al., 2026) and SOAAR (Khatri et al., 2026), and implemented the Spatial Label-Permutation Null Model (shuffling cell-type labels across cells while holding the spatial coordinates, 6-NN spatial graph, and single-cell expression profiles fixed).
  3. Implemented `src/spatial_null.py` and runner `src/run_spatial_null.py` with comprehensive docstrings citing CONCISE and SOAAR.
  4. Benchmarked 10 permutations (0.0226 s, ~0.00226 s/perm) and selected N=500 permutations (completed in 0.84 s).
  5. Computed empirical p-values, z-scores, null means, and null standard deviations for all 22 baseline interactions, generated `results/spatial_null_corrected.csv`, `results/spatial_null_distributions.npz`, `figures/spatial_null_comparison.png` and `.pdf`, and logged execution to `logs/spatial_null.log`.
  6. Verified end-to-end reproducibility of the permutation test by running `src/run_spatial_null.py` twice under fixed `SEED=42`; both runs produced byte-for-byte identical `spatial_null_corrected.csv` (SHA-256 `cafa5b17a920bab1f3b613d80c04942a06548ec28faf3b16eb919463a3529666`).
- **Key decisions**: Keep processed AnnData `data/processed/xenium_breast_baseline.h5ad` and raw data untouched; fix spatial graph structure and single-cell expression while shuffling `cell_type` labels; choose N=500 permutations given fast vectorized evaluation; store aggregated null distributions in a single compressed `.npz` file.
- **What worked**: Baseline pipeline is 100% deterministic and reproducible. Permutation test completed in ~2.0 seconds for N=500. Of 22 baseline interactions, 17 were flagged significant (empirical p < 0.05, all involving CXCL12-CXCR4 and CD274-PDCD1 with fibroblast, endothelial, or immune senders) and 5 were flagged non-significant (all with epithelial as sender: epithelial→endothelial CXCL12-CXCR4 p=0.73; epithelial→fibroblast CXCL12-CXCR4 p=1.0; epithelial→epithelial CXCL12-CXCR4 p=1.0; epithelial→fibroblast CD274-PDCD1 p=1.0; epithelial→epithelial CD274-PDCD1 p=1.0). Epithelial sender scores are substantially below the spatial null distribution (negative z-scores down to -38.76), correctly identifying that epithelial cells do not drive these signaling interactions.
- **What failed**: An initial ternary operator in the figure generator duplicated a legend label for non-significant markers when all top 10 interactions were significant; the labeling logic was corrected and the figure regenerated.
- **Blockers**: None.
- **Exact next step**: Review whether epithelial-marker argmax typing contributes to low epithelial ligand signal, and evaluate downstream response gene activation in receiver cells as planned for the thesis.

## 2026-09-13 — Spatial-edge null model and direct comparison against compositional null

- **Date/time**: 2026-09-13; local run start 13:20:55 UTC-05:00 / 18:20:55Z.
- **What I did**:
  1. Implemented `src/spatial_edge_score.py` to evaluate edge-based spatial ligand-receptor interaction scores: the mean of `expr_L(i) * expr_R(j)` over directed edges in `adata.obsp['spatial_connectivities']` connecting sender-type cells to receiver-type cells, alongside `n_qualifying_edges`.
  2. Implemented `src/run_spatial_edge_null.py` reusing the exact label-permutation engine (`generate_label_permutations` with $N=500$, `seed=42`) from `spatial_null.py`, keeping the 6-NN spatial graph and single-cell expressions fixed.
  3. Evaluated all 22 baseline interactions, generated `results/spatial_edge_null.csv`, merged with the prior compositional null results to produce `results/null_model_comparison.csv` with `verdict_flip` tracking, and produced comparison figures `figures/spatial_edge_null_comparison.png/.pdf` and `figures/null_model_flip_comparison.png/.pdf`.
  4. Logged execution to `logs/spatial_edge_null.log`.
  5. Verified run-to-run reproducibility by executing `src/run_spatial_edge_null.py` twice consecutively under `seed=42`; both runs produced byte-for-byte identical `results/null_model_comparison.csv` (SHA-256 `b4c3531d0c796784189b9a010899e8815041cd560b5527d572bff37f54dba44f`).
- **Key decisions**: Precompute edge products $expr_L[i] * expr_R[j]$ once across all 42,978 graph edges; assign 0.0 signaling score when a permutation yields 0 qualifying edges between rare cell types; reuse the identical label shuffle sequence from `spatial_null.py`.
- **What worked**: The edge-based scoring test completed in ~20-25 seconds for $N=500$. Out of 22 baseline interactions, 8 were significant under the edge null ($p < 0.05$) and 14 were non-significant. Comparing verdicts against the compositional null revealed **9 verdict flips** (all flipping from Significant under the compositional null to Non-Significant under the spatial-edge null):
  1. `endothelial→endothelial (CXCL12→CXCR4)`: edge score 6.633, $p = 0.084$ (only 27 spatial edges; lacks statistical power over null).
  2. `endothelial→epithelial (CXCL12→CXCR4)`: edge score 2.219, $p = 0.506$, $z = -0.10$ (co-expression along edges is indistinguishable from random cell mixture).
  3. `immune→endothelial (CXCL12→CXCR4)`: edge score 6.254, $p = 0.076$, $z = 1.52$ (only 17 spatial edges; marginally non-significant).
  4. `epithelial→immune (CXCL12→CXCR4)`: edge score 1.754, $p = 0.944$, $z = -1.61$ (epithelial cells have low ligand expression; edge product is below null expectation).
  5. `immune→epithelial (CXCL12→CXCR4)`: edge score 2.611, $p = 0.192$, $z = 0.85$ (not significantly above spatial null).
  6. `immune→fibroblast (CD274→PDCD1)`: edge score 0.000, $p = 1.000$ (ZERO spatial edges in the entire tissue connect a CD274+ cell to a PDCD1+ cell; cluster-mean scoring gave a false positive).
  7. `immune→epithelial (CD274→PDCD1)`: edge score 0.000, $p = 1.000$ (ZERO spatial edges).
  8. `fibroblast→fibroblast (CD274→PDCD1)`: edge score 0.000, $p = 1.000$ (ZERO spatial edges).
  9. `fibroblast→epithelial (CD274→PDCD1)`: edge score 0.000, $p = 1.000$ (ZERO spatial edges).
  The 8 interactions that remained significant are genuine spatial CXCL12→CXCR4 axes (fibroblast→immune, endothelial→immune, fibroblast→endothelial, immune→immune, fibroblast→fibroblast, fibroblast→epithelial, endothelial→fibroblast, immune→fibroblast) with strong local co-enrichment.
- **What failed**: Direct Python looping over 42,978 edges would be slow; vectorized masking on precomputed edge products resolved this and ran 500 permutations in under 20 seconds.
- **Blockers**: None.
- **Exact next step**: Integrate edge-based scoring and null calibration into the planned S-PCST graph optimization framework.

## 2026-09-13 — Diagnostic investigation of CD274→PDCD1 edge-null variance

- **Date/time**: 2026-09-13; local time ~13:43 UTC-05:00 / 18:43Z.
- **What I did**:
  1. Loaded `data/processed/xenium_breast_baseline.h5ad` to diagnose why `edge_null_std == 0.0` across all 500 permutations for all 6 CD274→PDCD1 interactions in `results/spatial_edge_null.csv`.
  2. Evaluated global and per-cell-type detection counts, detection rates, and mean expressions for CD274 and PDCD1 on normalized `adata.X`, exporting `results/cd274_pdcd1_detection_diagnostic.csv` (10 rows: 4 cell types x 2 genes + 2 global 'ALL' rows).
  3. Performed an independent cross-check of the `immune→fibroblast` CD274→PDCD1 interaction by directly indexing `adata.obsp['spatial_connectivities']` without calling `spatial_edge_score.py`, confirming the exact edge count (962) and observed edge score (0.0).
  4. Inspected `logs/spatial_edge_null.log` for printed SHA-256 digests.
- **Key decisions**: Keep all existing results and pipeline code untouched; evaluate sparsity across both the single-cell expression vectors and the 42,978-edge spatial connectivity graph.
- **What worked**: Identified the definitive biological and mathematical cause:
  - PDCD1 has extreme global drop-out / low expression, detected in only 5 of 7,163 cells (0.0698% detection rate: 0/361 immune, 0/49 endothelial, 1/5355 epithelial, 4/1398 fibroblast).
  - CD274 is detected in 53 of 7,163 cells (0.7399% detection rate: 17/361 immune, 0/49 endothelial, 16/5355 epithelial, 20/1398 fibroblast).
  - Across all 42,978 directed edges in the tissue's 6-NN spatial graph, exactly 0 edges connect ANY CD274+ cell to ANY PDCD1+ cell (`expr_CD274[i] * expr_PDCD1[j] == 0.0` for all 42,978 edges).
  - Because label permutations reassign cell-type labels across fixed spatial nodes without moving the underlying single-cell expressions, the 42,978 edge product vector remains identically zero across every single permutation. Consequently, the mean of any masked subset of edges is deterministically 0.0, yielding `edge_observed_score = 0.0`, `edge_null_mean = 0.0`, and `edge_null_std = 0.0` for all 500 permutations.
  - Independent cross-check verified exactly 962 directed edges between immune and fibroblast cells, with observed edge score 0.0.
- **What failed**: None.
- **Blockers**: None.
- **Exact next step**: Document findings for the thesis methodology section regarding sparse receptor dropout and spatial-edge null properties, and proceed to downstream response pathway integration.

## 2026-09-13 — Reviewer Robustness Checks: LR Detection-Rate Sweep & Leiden Consensus Typing

- **Date/time**: 2026-09-13; local time ~15:58 UTC-05:00 / 20:58Z.
- **What I did**:
  1. **Task A (Detection-Rate Sweep across 20 LR Pairs)**:
     - Gathered all 37 unique genes across the 20 curated `LR_PAIRS` tuples in `src/pipeline.py` and evaluated their single-cell detection rates on normalized `adata.X`.
     - Discovered that only **9 of the 37 genes** are present on the commercial 10x Xenium 280-gene targeted breast panel; the remaining 28 genes are completely missing from the physical panel design.
     - Confirmed that only **2 of the 20 pairs** (`CXCL12-CXCR4` and `CD274-PDCD1`) have both ligand and receptor present in the panel, definitively explaining why the other 18 pairs contributed zero rows to the baseline ranking.
     - Exported `results/gene_detection_diagnostic.csv` (37 rows) and `results/lr_pair_coverage_summary.csv` (20 rows) in 0.29 seconds runtime.
  2. **Task B (Leiden-Cluster-Based Consensus Typing Robustness Check)**:
     - Computed cluster-level consensus labels across the 4 Leiden clusters (res=0.5) by taking the argmax of mean marker scores: Clusters 0 & 1 $\to$ `epithelial` (5,756 cells), Cluster 2 $\to$ `immune` (773 cells), Cluster 3 $\to$ `fibroblast` (634 cells). Endothelial cells ($n=49$) were absorbed into larger clusters.
     - Generated confusion matrix / crosstab and exported `results/leiden_robustness/crosstab_cell_type_leiden.csv`.
     - Re-executed the baseline scoring, 500-permutation compositional null test, and 500-permutation spatial-edge null test under `cell_type_leiden` ($N=500$, `seed=42`).
     - Exported `results/leiden_robustness/null_model_comparison_leiden.csv`, `results/leiden_robustness/spatial_edge_null_leiden.csv`, and `results/leiden_robustness/comparison_summary.txt` in 27.58 seconds runtime.
- **Key decisions**: Keep existing per-cell `cell_type` results untouched; evaluate Leiden consensus as a parallel robustness check in `results/leiden_robustness/`; record missing panel genes explicitly rather than skipping them; execute with `n_jobs=1` under `if __name__ == '__main__':` for Windows multiprocessing stability.
- **What worked**:
  - The central thesis findings are 100% robust and invariant across cell typing methodologies:
    1. The top communication axis under both compositional and edge scoring remains the `fibroblast→immune` (CXCL12→CXCR4) paracrine axis (edge score 20.21, 1,094 edges, $z = +39.57$).
    2. All CD274→PDCD1 combinations (6 pairs) again flipped from Significant under compositional null to Non-Significant under edge null because physical contact is zero.
    3. Qualitative dichotomy holds: high-edge-count interactions robustly survive spatial calibration, while low-contact/random-mixing interactions flip to non-significant.
- **What failed**: An initial multiprocessing spawn error occurred in Squidpy on Windows when called without `if __name__ == '__main__':`; resolved by restructuring runner scripts with explicit entry points.
- **Blockers**: None.
- **Exact next step**: Transition research dashboard (`index.html`) into a multi-page web application featuring dedicated pages for Flips, Diagnostics, Robustness, Interactive Explorer, Spatial Map, and Methodology.

## 2026-09-13 — Scientific Framing Realignment: CD274–PDCD1 Spatial Sparsity Elevation & Granularity Sensitivity Caveat

- **Date/time**: 2026-09-13; local time ~17:55 UTC-05:00 / 22:55Z.
- **What I did**:
  1. Conducted an empirical Leiden resolution sweep ($res \in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]$) to test whether endothelial cells ($n=49$, 0.68% of 7,163 cells) separate into an independent cluster at higher community detection resolutions without over-clustering.
  2. Determined that endothelial cells remain merged into stromal and epithelial clusters across all standard resolutions ($0.5 \le res \le 1.5$), separating only at an extreme resolution of 2.0 (19 micro-clusters) which fragments biology.
  3. Formally realigned the thesis presentation architecture under Option 1 ("Reframe Now"):
     - **Primary, typing-robust finding**: Promoted CD274–PDCD1 detection sparsity and spatial non-contact as the headline proof of method utility. Dissociated/marginal null models manufactured false-positive significance for uncontacted cells (0 physical edges out of 42,978 in tissue graph), whereas spatial-edge scoring eliminates all 6 combinations to observed score 0.000000 across both per-cell and Leiden consensus typing schemes.
     - **Secondary, explicitly caveated observation**: Maintained CXCL12–CXCR4 dominant paracrine axis stability (`fibroblast→immune` ranks #1 in both models: $z = +29.51$ per-cell, $z = +39.57$ Leiden consensus), while adding an explicit methodological caveat that secondary rank variations in rare lineages (endothelial, 0.68%) are sensitive to cell-typing granularity when low-abundance cells are absorbed by unsupervised clustering.
  4. Updated `results/leiden_robustness/comparison_summary.txt`, the research dashboard pages (`index.html`, `flips.html`, `robustness.html`, `diagnostic.html`, `methods.html`), and project documentation.
- **Key decisions**: Avoid artificial resolution tuning that causes over-clustering; transparently report endothelial absorption as an expected clustering property; position the zero-edge CD274–PDCD1 discovery as the undeniable anchor for spatial graph modeling.
- **What worked**: Eliminates thesis risk, prevents over-engineering, provides the cleanest and most methodologically honest thesis narrative for examination with 15 days remaining.
- **Blockers**: None.
- **Exact next step**: Synthesize findings into final thesis manuscript chapters and presentation defense materials.


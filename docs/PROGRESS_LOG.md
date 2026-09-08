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

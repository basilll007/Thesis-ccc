# PDAC Extension — Progress Log (branch `pdac-extension`)

Append-only log. Every number here comes from a file listed in `results/pdac/MANIFEST.csv`. Interpretation is limited to audit and diagnostic statements. Nothing here supports a biological or clinical claim.

## Task checklist

- [x] Phase 0: Repo hygiene (PDCD1 caption, F:\ paths, settings note, "replication" renamed)
- [ ] Phase 1: Matched dyadic baseline (breast)
- [ ] Phase 2: Find and vet a PDAC spatial dataset
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

# model_dev v2 progress log

Branch `model-v2`, created 2026-11-16 from origin/main (71ebd43).

**Rules:**
- Never push to main and never force-push; push only `model-v2`.
- Plans are committed before each phase runs.
- Every number comes from a generated file.
- These are expression data (not mutation), and results are associations, not causal claims.
- The tumour is the unit for every CI.
- Phases run one at a time on the GPU, with CPU workers capped at 4 and free RAM logged before each phase (minimum 8 GB).

## 2026-11-16 Setup
- **Working tree.** F:/Thesis was clean and on main, and no git operation was in progress. Another agent's worktree is F:/Thesis-ccc-translational (branch translational-v1).
  - The branch was created in place, as in v1. The data and model caches are git-ignored paths relative to this tree, and a new worktree would not see them.
- **Environment.** `thesis-model` env (CUDA torch 2.11, xgboost). Free RAM at start: 10.1 GB. GPU: RTX 5060 Laptop, 910 MiB in use.
- **Parquet.** No parquet engine is installed, so the per-cell metadata cache is written as `model_cache/v2/cell_meta.csv.gz` (no new dependency).
- **Pre-plan descriptives** (no target values read): `results/v2_phase1/preplan/`.
  - Median cell diameter: 12.14 um.
  - Distance bins: touching 42 %, near 39 %, separated 19 %.
  - No CD68 protein channel; PanCK does not separate macrophages from CAFs.
  - CD45 rule: ~62 % of macrophages are above the same-section CAF 75th percentile.
- **XGBoost timing** (one fit, 60,000 x 3,300 features, 500 rounds): ~22-25 s, plus a refit of ~13 s. So Phase 1 uses the v1-selected XGBoost configuration (rounds still tuned on inner validation), and Phase 2 uses the full v1 grid (declared in the plans).
- **Plans** `plans/v2_phase1_spillover.md` and `plans/v2_phase2_crosstype.md` were committed before any Phase 1 or Phase 2 computation.

## 2026-11-16 Phase 1: spillover control (DONE)
- **prepare** (free RAM 9.97 GB):
  - The contamination-marker rule needed the pre-declared fallback. Only 2 genes (COL6A3, COL3A1) reached ratio_CAF/iso >= 10, so the ratio-5 threshold was used: 20 markers, all collagens and stromal genes.
  - S (5 genes): HLA-DRB1, HLA-DRA, SPP1, CD163, ITGAX.
  - Target genes removed (CAF/iso >= 0.5): SELENOP, CD9, LGALS3.
- **fit** (GPU, about 50 min, one process): groups ref, sc, cd45_gt_caf75 and cd45_gt_caf90, all 5 folds.
- **Decision JSON fix.** The first decision JSON serialised numpy booleans as strings. The decision was re-run from the saved summary table (`decide` step; no refit, no re-summarise). The verdict was unchanged: **explained by spillover**, with checks a, c and d failing.
- **Post-hoc table** (`f_posthoc_mixing_expectation.csv`, labelled exploratory and not in the plan): the expected level of each target gene in touching macrophages if they were pure mixtures of isolated-macrophage and CAF transcripts at the estimated 23 % rate.

## 2026-11-16 Phase 2: first attempt failed (failure 1)
- The single-process myeloid run completed fold 0 (950 s), then failed in fold 1 with an XGBoost CUDA out-of-memory error ("Requested 1.2 GB, free 4.5 GB"). Only desktop applications share the GPU, so the memory build-up was inside the process across folds.
- **Fix:**
  - one outer fold per process (`--fold-budget 1`, chained in `model_cache/v2/run_phase2_chain.sh`), as in v1;
  - XGBoost training retried once after releasing GPU caches (identical inputs).
- Fold 0 is kept (resume).

## 2026-11-16 Phase 3: survival (DONE, run between Phase 2 runs, not in parallel)
- **Amendment** `plans/v2_phase3_survival.md` was committed before any outcome model.
- **Run:** `results/translational/phase3/20261117T024638285521Z/`.
  - PH validation: PASS.
  - 108 of 108 models were estimable (4 cohorts x 9 signatures x 3 families).
- **Primary result:** no signature has adjusted BH q < 0.05 in strict TCGA.
- **Coverage:** SELENOP is absent from the TCGA expression matrix, so the TCGA C1Q signature uses 5 of 6 genes (the plan's coverage rule is met); CPTAC has 6 of 6. TRAC is absent in both cohorts (T-cell proxy 8 of 9).

## 2026-11-16 Phase 2: attempt 2 BLOCKED; environment lost; attempt 3
- **Attempt 2 (one fold per process) failed twice with the same error class:**
  - XGBoost CUDA out-of-memory errors at the start of fold 1 in fresh processes, and one host allocation failure.
  - Windows commit charge was ~46 of 58 GB, mostly held by other applications (as in v1).
  - Per the stop rule this configuration is **BLOCKED**. The tumour and adjusted runs then also failed fast.
  - The adjusted run completed fold 0 for both of its groups before the environment disappeared.
- **Environment lost at ~21:43.** The D: volume holding Anaconda and the `thesis-model` env vanished from the system: no D: volume is listed, and the last chain process returned rc=127 because python.exe was not found. Nothing was changed on any disk to recover it.
- **Attempt 3 (materially different, declared before any Phase 2 result was summarised):**
  - The pinned v1 environment (`model_dev/environment.yml`, plus xgboost 3.2.x and psutil) is rebuilt with uv at `model_cache/envs/thesis-model-v2` (git-ignored; caches and temp on F:).
  - The XGBoost loop now keeps at most two DMatrices alive and builds each from a GPU column slice. Each XGBoost DMatrix holds a dense host copy of about 1.4 GB, and attempt 2 had four alive at once. The inputs and seeds of every fit are unchanged.
  - To avoid mixing environments, **all** Phase 2 folds are recomputed in the new environment. The 3 old-environment fold-0 files are archived, unused, in `model_cache/v2/phase2_attempt1_oldenv/`.
  - Phase 1 and Phase 3 results (old environment) are unchanged. The reproduction check against v1 Phase 2 tests the new environment.

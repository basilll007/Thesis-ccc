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

# v3 Phase 3 plan — effect sizes and tumour-level model

Written 2026-11-15 BEFORE any Phase 3 result was computed. Branch `thesis-v3`.

## Goal
Replace z-only reporting with effect sizes and uncertainty for the headline CXCL12->CXCR4 CAF_fibroblast->T_cell axis,
and report how z scales with edge count.

## 3a. Q with spatial block-bootstrap 95% CI per section
- Data: as Phase 1 (`v3_nulls.load_section`; 6-NN graph; log1p(normalize_total 1e4) expression).
- Annotations: consensus_T (primary), marker, author_mapped. Axes: the headline for all three annotations;
  additionally the 7 non-background axes for CAF->T under consensus_T.
- Blocks: square tiles of 100 um anchored at the section min x/y (same construction as `v3_nulls.tile_codes`).
- Resample: B = 1000 replicates, numpy default_rng(42). Each replicate draws n_tiles tiles with replacement; tile
  multiplicity w_t becomes the weight of every cell in the tile (w_i).
  E* = sum over CAF->T edges of w_row(i) L_i R_j / sum w_row(i) (an edge belongs to the block of its sender cell);
  M* = (sum_CAF w_i L_i / sum_CAF w_i) x (sum_T w_j R_j / sum_T w_j); Q* = E*/M*.
  Replicates with zero CAF->T edge weight are dropped and counted.
- Report per section: Q_obs, bootstrap percentile 2.5 / 97.5 (and log-scale symmetric SE), n CAF->T edges,
  n CAF, n T, n tiles, n dropped replicates; CI excludes 1 (yes/no).
- Output `results/pdac/v3_phase3/Q_block_bootstrap_by_section.csv` (+ `_nonbg_axes.csv`).

## 3b. Tumour-level mixed model (exploratory)
- Observations: the 16 sections (one per row), response y = log Q_obs for the headline; grouping = tumour
  (16 sections -> 13 tumours, Phase D2 mapping). Model: y ~ 1 + (1 | tumour), statsmodels MixedLM, REML.
- Report: intercept (= population mean log Q), 95% Wald CI, exp() of both (geometric-mean Q), tumour variance,
  n tumours, n sections, convergence flag; per annotation. Also the fixed-effects-only comparison (log Q ~ 1, i.e. a one-sample
  t-interval over sections) to show what the random intercept changes. Weighted variant: none (declared).
- Label exploratory: 16 observations, 13 groups, 10 singleton tumours.
- Output `results/pdac/v3_phase3/mixed_model_logQ.csv`.

## 3c. z-vs-edge-count inflation
- From the Phase 1 global-null files (`results/pdac/v3_phase1/sections/*__global.csv.gz`): headline rows per section and the
  7 non-background axes (all sender/receiver rows), 3 annotations.
- Report raw z (E, M, Q) next to z / sqrt(n_qualifying_edges); Spearman(z, n_qualifying_edges) per score and annotation;
  the same for the normalised z. Also the slope of |z_E| on log10(n_edges).
- Output `results/pdac/v3_phase3/z_vs_edge_count.csv`, `z_vs_edge_count_summary.csv`.

## Figure: figures/pdac/v3_phase3/phase3_effect_sizes.{png,pdf}
A: headline Q with 95% block-bootstrap CI per section (consensus_T), sorted by Q, dashed line at 1, point size by edge count;
B: mixed-model intercept and CI per annotation (geometric-mean Q scale);
C: z_E and z_E/sqrt(n_edges) against n_edges for the non-background catalogue (consensus_T).

## Interpretation rules (declared now)
- Effect sizes are reported as audit quantities. No biological or clinical claim.
- "Supported" would require: mixed-model CI for geometric-mean Q excluding 1 AND the Phase 1 constrained nulls not shrinking
  the signal; otherwise "exploratory" or "not supported".

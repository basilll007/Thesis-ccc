# v2 Phase 2 plan: cross-type share on the aggregated baselines

Committed before any Phase 2 result is computed (written together with the Phase 1 plan, before Phase 1 ran). Development set only; GSE310352 is not used. Folds, inner-validation rule, the v1 feature caches (`development__myeloid.npz`, `development__tumour.npz`, tumour receivers capped at 8,000 per section as in v1), the feature registry and the leakage guard are reused unchanged.

## Question
How large is the other-cell-type share of the neighbourhood gain when it is measured on the best (aggregated) models rather than inside the graph model? v1 graph estimate for myeloid C1Q_TAM: cross-type contribution +0.019 [0.004, 0.028].

## Models (same protocol as v1 Phase 2)
- **XGBoost-CUDA on the `genes` blocks** (own, tile, composition, nb_programs, nb_genes); program targets only, one model per target.
  - Full v1 grid: learning_rate {0.05, 0.1} x max_leaves {15, 31} (lossguide, min_child_weight 50).
  - Rounds (<= 500) chosen per fold on the inner-validation RMSE. Refit on all four training folds with <= 60,000 seeded rows.
  - The configuration is selected per (receiver, target, variant) by the lowest inner-validation MSE averaged over folds, exactly as v1 `summarise`.
- **Ridge (GPU closed form) on the `programs` blocks** (own, tile, composition, nb_programs); alpha {0.1, 1, 10, 100, 1000} on inner validation, then refit on all training folds.
  - Ridge covers the programs and the gene vectors (multi-output). XGBoost is not fitted to gene vectors (v1 scope, logged).
- **Reference:** own = ridge on own + tile (v1).

## Variants (feature columns split by the neighbour cell type they describe; `v2common.column_type`)
Types: tumour, CAF, T, myeloid, other (all untyped classes, including the `other` gene block). Own and tile columns are always kept.
1. **all**: every neighbour column (v1 model);
2. **same**: same-type columns only (composition, programs and genes of the receiver's own type);
3. **other_only**: all other-type columns, no same-type columns;
4. **minus_X**: all columns except type X, for each other type X (myeloid receiver: CAF, tumour, T, other; tumour receiver: CAF, myeloid, T, other);
5. **same_plus_other_comp**: same-type columns plus the composition columns (fractions and counts) of the other types, without their expression.

## Quantities (pooled held-out R2; paired tumour bootstrap, 2,000 draws)
- **Neighbour gain** = R2(all) - R2(own).
- **Cross-type contribution** = R2(all) - R2(same).
- **Cross-type share** = (R2(all) - R2(same)) / (R2(all) - R2(own)), bootstrapped as a ratio of summed SSE differences per resample. Reported only if the neighbour-gain CI excludes 0; otherwise "not defined".
- **Same-type contribution** = R2(all) - R2(other_only).
- **Unique contribution of type X** = R2(all) - R2(minus_X).
- **"Who is nearby"** = R2(same_plus_other_comp) - R2(same).
- **"What they express"** = R2(all) - R2(same_plus_other_comp).

## Targets
- **Myeloid:** C1Q_TAM and SPP1_TAM (XGBoost and ridge); 11-gene vector (ridge).
- **Tumour:** classical and basal_like (XGBoost and ridge); 18-gene vector (ridge, reported but not requested).
- **Spillover-adjusted re-run (myeloid):**
  - targets from the Phase 1 rule: C1Q_TAM_adj and SPP1_TAM_adj (XGBoost and ridge), gene_vector_adj (ridge);
  - with the Phase 1 contamination score s_c added to the own block in every variant, including the own reference;
  - plus the unadjusted gene-mean targets C1Q_TAM_mean and SPP1_TAM_mean (XGBoost and ridge, no s_c) as the like-for-like reference;
  - targets that Phase 1 marks not estimable are skipped, with the reason.

## Agreement with the graph-model estimate
- For myeloid C1Q_TAM, the XGBoost cross-type contribution **agrees** with the graph estimate (+0.019 [0.004, 0.028]) if the XGBoost CI contains 0.019 and the graph CI contains the XGBoost point estimate.
- Otherwise it is reported as larger or smaller.
- The same comparison is made for SPP1_TAM (graph +0.015 [0.004, 0.028]), the gene vector (ridge vs graph +0.014 [0.009, 0.018]) and tumour classical (graph +0.007 [-0.004, 0.024]).

## Reproduction check
The `all` variant re-fits the v1 models.
- Pooled R2 should match v1 (`results/phase2/baseline_summary.csv`) within 0.005 for XGBoost (GPU histogram sums are not bit-deterministic) and 0.001 for ridge.
- A larger difference is logged and investigated before reporting.

## Outputs, compute, safety
- **Outputs:** `model_dev/results/v2_phase2/` (per-variant summary, quantities table, agreement table, figure, logs). Per-(task, fold) predictions go to `model_cache/v2/phase2/` (resumable).
- **Compute:**
  - GPU only; CPU workers capped at 4;
  - the run order is myeloid, then tumour, then the adjusted re-run, one at a time;
  - free RAM is logged before each run, and a run does not start below 8 GB.

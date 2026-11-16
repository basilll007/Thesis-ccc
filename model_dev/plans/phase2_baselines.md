# Phase 2 plan: baselines

Committed before `baselines.py` is run. Folds: `model_dev/results/phase1/folds.csv` (grouped 5-fold); targets: `model_dev/results/phase1/targets.json`. Receivers: myeloid (primary; programs C1Q_TAM, SPP1_TAM; 11-gene vector) and tumour (positive control; classical, basal_like; 18-gene vector).

## Models, on identical folds
1. **Mean predictor:** the training-fold mean of each target. Pooled R2 is 0 by construction; it anchors the scale.
2. **Own-covariates model:** ridge on the `own` block (log1p total counts, log1p detected genes, density within 30 um, mean 6-NN edge length) plus the `tile` block (leave-one-out 200 um tile means of the targets over same-type cells, with an empty-tile indicator). No neighbour information.
3. **Aggregated-neighbourhood models:**
   - ridge on feature set `programs` = own + tile + composition + neighbour program means (the chain Phase 2 set);
   - ridge on feature set `genes` = `programs` + per-type neighbour gene means over 6-NN (5 x 539 columns);
   - gradient boosting (sklearn `HistGradientBoostingRegressor`) on `programs` and on `genes`, program targets only (one model per target column; the gene-vector target would need 11-18 models per grid point per fold, so it is restricted to ridge and the MLP - logged scope).
4. **MLP** (torch) on `genes`: hidden (256, 128), ReLU, dropout 0.1, AdamW, batch 1024, up to 100 epochs, early stopping on the inner validation fold (patience 10), multi-output over all target columns of the receiver. 3 seeds (0, 1, 2); the reported number is the mean over seeds, with per-seed values in the table.

## Hyperparameter grid (pre-declared; chosen on the inner validation fold only; every configuration is reported)
- ridge alpha in {0.1, 1, 10, 100, 1000} (features standardised with training-fold statistics; NaN tile means imputed with the training mean).
- HGB: learning_rate in {0.05, 0.1} x max_leaf_nodes in {15, 31}; max_iter 500, with the iteration count chosen by the inner-validation loss along `staged_predict` (one fit per grid point); min_samples_leaf 50; at most 60,000 training rows (random subset, seed 42) for the `genes` set.
- MLP: learning rate in {1e-3, 3e-4}; weight decay 1e-4.
- After selection on the inner fold, ridge and HGB are refitted on all four training folds; the MLP keeps its early-stopped weights (trained on 3 folds).

## Metric and uncertainty
- Pooled held-out R2 with the training-fold mean in SST (Phase 1 plan). Per-tumour (SSE, SST) pairs are saved; 95 % CIs resample tumours (2,000 draws). Gains versus the own-covariates model use the paired tumour bootstrap.
- Primary endpoint: pooled R2 for `myeloid:C1Q_TAM` and its gain over own covariates.

## Reproduction check of chain Phase 2 (run first; stop if it fails)
- Re-run the chain Phase 2 ridge analysis from the Phase 0 datasets with the chain definitions: leave-one-tumour-out (13 folds), within-section demeaning of targets and features, the chain own covariates (log total counts, density, tile LOO mean of the program, empty flag), composition (6-NN fractions and 30 um counts over all classes), neighbour program means per typed class (6-NN and 30 um), a 20,000-cell cap per section with the chain seeds (`[42, section_index, type_index]`), `RidgeCV` with alphas 10^-2..10^4 (13 values), author labels, no spillover covariates.
- Compare delta R2 (neighbour - baseline) with `results/chain/phase2/delta_r2_by_program.csv` (author_mapped, no spillover, ridge) for the four receiver programs. **Tolerance: |difference| <= 0.01.** If any program fails, Phase 2 stops and the cause is found before continuing.

## Outputs
`model_dev/results/phase2/`: `reproduction_check.csv`; `baseline_results.csv` (every model x feature set x configuration x fold with inner-validation and test metrics); `baseline_summary.csv` (selected configuration per model with pooled R2, CI, gain over own covariates with CI); `per_tumour_sse.csv`; a figure. The strongest baseline per target (by pooled R2) is recorded for Phase 3.

## Success criterion
The reproduction check passes; every model runs on all five folds for both receivers; the strongest baseline per target is identified with a CI.

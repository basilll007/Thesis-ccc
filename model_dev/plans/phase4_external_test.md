# Phase 4 plan: external test (one touch)

Committed before Phase 2 results are summarised and before the external set is used by any model. GSE310352 was loaded once in Phase 0 for the data card only (`model_dev/results/phase0/external_first_load.json`, commit 58ce3c9).

## What is frozen
- **Preprocessing:** the Phase 0 pipeline (QC, panel intersection, log1p CP10k, gene set G, per-section label-free program scoring with fixed seeds, graphs scaled by nearest-neighbour distance because um/px is UNVERIFIED). No gene selection, threshold or scaler is fitted on the external set; all scalers use development-set statistics (all 13 tumours).
- **Models**, each retrained once on the full development set with hyperparameters fixed from the internal cross-validation (no external tuning):
  1. mean predictor (development mean);
  2. own-covariates ridge (alpha = the alpha selected most often across folds; ties -> the larger alpha);
  3. the strongest Phase 2 baseline per target (same rule for its hyperparameters; HGB uses the median selected iteration count; the MLP uses the median early-stopped epoch count, 3 seeds averaged);
  4. the graph model selected in Phase 3 (architecture and learning rate from `architecture_selection.json`; epochs = median early-stopped epoch count over folds and seeds; 3 seeds averaged), whether or not it beat the baseline internally.
  Evaluating this fixed list in one run is the single touch of the external set.

## Receivers and labels
- Receivers: myeloid and tumour by **marker labels** (the external set has no author labels; the development models were trained with author labels). Mendeley marker-label counts are in the data card for comparison.
- Targets: the same program definitions and declared gene sets.

## Reported
- Pooled R2 per model and target over all external cells, with the development-set mean in SST (the frozen mean predictor), and a **unit bootstrap** CI (6 units: slides 5, 5b2, 5b3 grouped; the grouping is UNVERIFIED and the units are not verified independent patients - stated in every table).
- Internal-to-external drop: R2(internal CV, Phase 2/3 summary) - R2(external) per model and target.
- Per-slide breakdown (8 slides).
- The gain of each model over the frozen own-covariates model, paired unit bootstrap.

## Decision rules
- "Transfers": the external gain over own covariates has a unit-bootstrap CI excluding 0. Because 6 units give wide intervals, the point estimates are reported alongside and no claim of replication is made.
- No model is modified after the external numbers are seen. If anything is changed later, the external result is marked stale and not reported as a test.

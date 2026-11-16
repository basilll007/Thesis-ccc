# Chain Phase 2 plan: does the neighbourhood predict a cell's state?

Committed after Phase 1 (programs, counts and graphs only) and before any Phase 2 model is fitted. Seed 42. Annotation: `author_mapped` (primary) and `marker` (sensitivity). Unit of inference: the tumour (13 tumours).

## Receivers and outcomes
- Receiver types are tumour, CAF, T cell and myeloid. Outcomes are every program kept in Phase 1 Amendment 1 for that type (results/chain/phase1b/programs_used.json; cache data/chain_cache/mendeley_v1b), including the CAF CXCL12 single-gene readout, which stays labelled as such.
- Rows are cells of the receiver type in all 16 sections. Compute cap: at most 20,000 randomly chosen cells per section and receiver type (seed 42). This affects tumour and CAF in the large sections. All T and myeloid cells are kept wherever they number <= 20,000.

## Features
Each block is computed per section from the Phase 1 cache, with the symmetrised 6-NN graph and the 30 um radius graph.
- **Baseline** (own covariates only):
  - log1p total counts.
  - Local cell density: number of cells within 30 um.
  - Zonation control: leave-one-out mean of the outcome program over cells of the receiver type in the same 200 um tile (grid anchored at the section min x/y). Where the tile has no other cell of that type, the section mean is used and an indicator is set.
  - Section fixed effect: the outcome and every feature are demeaned within section, using the rows in the analysis. Held-out R2 is therefore within-section R2. Section dummies cannot be estimated for a held-out tumour, so this is the stated reading of "section" as a covariate.
- **Neighbour model** = baseline plus:
  - Composition: the fraction of each annotation class among the 6-NN neighbours, and the count of each class within 30 um.
  - Mean neighbour program scores by type: for each of the four typed types t and each kept program of t, the mean score over neighbours of type t, for both 6-NN and 30 um. When there is no neighbour of type t, the value is set to the section mean of that program among type t, so that it is 0 after demeaning. Absence is carried by the composition features.
- **Spillover sensitivity** (both models): add the receiver's own four lineage-marker scores (Phase 1). CosMx segmentation spillover can make a cell look like its neighbours. A neighbour gain that disappears under this covariate is flagged "spillover-compatible".

## Models and cross-validation
- Leave-one-tumour-out cross-validation: 13 folds, with all sections of a tumour held out together.
- Ridge: features standardised on the training fold; alpha chosen by `RidgeCV`, alphas = 10^[-2..4] (13 values), efficient LOO on the training fold.
- Gradient boosting: `HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, early_stopping=True, validation_fraction=0.1, random_state=42)`, trained on at most 60,000 random training rows per fold.
- Metric: pooled held-out R2 = 1 - sum(SSE) / sum(SST) over held-out cells, where SST is the sum of the squared within-section-demeaned outcome. delta R2 = R2(neighbour) - R2(baseline).
- **Tumour-level bootstrap CI**: resample the 13 tumours with replacement 2,000 times and recompute pooled delta R2 from the per-tumour SSE and SST. Report a 95 % percentile CI, the per-tumour delta R2, and the number of tumours with delta R2 > 0.

## Which features carry the gain (ridge only)
- Leave-one-tumour-out, drop-group: the delta R2 lost when one group is removed from the full neighbour model.
- Add-group: the delta R2 gained when one group is added to the baseline.
- Groups: composition (6-NN and 30 um); same-type neighbour programs; and other-type neighbour programs, split by type (tumour, CAF, T, myeloid).

## Decision rules
- "Neighbours predict state beyond own covariates" for a receiver program: the ridge or gradient-boosting delta R2 has a tumour-bootstrap 95 % CI with lower bound > 0. Magnitude is always reported. A gain near zero is reported as such and is a valid result.
- **GNN decision.** A small graph model is added only if some receiver program has (a) a delta R2 lower CI > 0.01 and (b) a gradient-boosting delta R2 at least 0.01 above ridge, which indicates non-linearity that a more flexible neighbourhood model could use. torch is not installed. If the rule triggers, the model will be an SGC-style network: neighbour features propagated over 1 and 2 hops of the 6-NN graph, fitted with sklearn `MLPRegressor`. The decision and its inputs are logged in either case.
- Multiplicity: Phase 2 reports intervals, not tests; no BH is applied. The number of receiver-program x annotation x model combinations is stated.

## Outputs
results/chain/phase2/: per-fold predictions summary (per tumour SSE and SST), delta_r2_by_program.csv (ridge and gradient boosting, with and without spillover covariates, for both annotations), feature_group_attribution.csv, gnn_decision.json and phase2_summary.json. Figure figures/chain/phase2/: delta R2 with 95 % CI per receiver program.

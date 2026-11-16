# Phase 5 plan: interpretation (only for models that beat baselines)

Committed before Phase 3 results are summarised.

## Gate
Phase 5 runs for a receiver/target only if the Phase 3 selected graph model beats the strongest baseline (paired tumour-bootstrap CI of the gain excluding 0). Otherwise Phase 5 writes `results/phase5/gate.json` with the reason and stops. Attention weights are never presented as explanations.

## Tumour-level permutation importance (the gated analysis)
- Unit: the held-out tumour. For each outer fold and seed, with the trained Phase 3 model (checkpoint), compute the test R2 drop when a feature group is permuted across the held-out cells of that fold (within section, to keep section composition): groups = neighbour class (each annotation class, by replacing that class's neighbour nodes' expression with a random same-section cell's of the same class - i.e. breaking the link between the receiver and that class's state), and neighbour gene groups (the program gene lists and the lineage sets; all other G genes as one group).
- Importance = mean R2 drop over seeds; CI by resampling tumours of the per-tumour (SSE_perm - SSE_full) contributions (2,000 draws). Table `permutation_importance.csv`.

## In-silico perturbation (double gate)
Only if the cross-type contribution CI (Phase 3 ablation) excludes 0 for that target: for the neighbour program (among other-type classes) with the largest permutation importance, set the neighbour cells' expression of that program's genes to their 10th and 90th percentile (development distribution, within class), re-predict the held-out cells, and report the predicted change with a tumour-bootstrap CI. Every row is labelled "model-based hypothesis". No claim beyond the model.

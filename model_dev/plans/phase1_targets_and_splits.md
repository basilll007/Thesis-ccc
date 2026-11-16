# Phase 1 plan: targets and splits

Committed before `splits.py` is run and before any model is fitted. Inputs: Phase 0 datasets and `model_dev/results/phase0/cell_counts.csv`.

## Receivers
- **Primary receiver:** myeloid cells (`author_mapped == myeloid` in the development set; `marker == myeloid` in the external set). Programs: `myeloid:C1Q_TAM` (6 genes) and `myeloid:SPP1_TAM` (5 genes), as defined in `results/chain/phase1b/programs_used.json`.
- **Positive-control receiver:** tumour cells (`epithelial_ductal_tumor`). Programs `tumour:classical` (8 genes) and `tumour:basal_like` (10 genes). A neighbourhood gain is already known here (chain Phase 2 ridge delta R2 0.141 for classical).
- **T cells:** not a target. They are reported as "not testable" with their counts (11,557 author-labelled T cells in total; <= 3,333 per section; the chain Phase 2 gain was 0.0005 [-0.003, 0.005]).
- CAF cells are not a receiver in this baseline stage (scope), but their programs are neighbour features.

## Targets (per receiver cell)
- **(a) Program scores:** the Phase 0 scanpy-style scores (per-section, label-free control selection; 25 bins, 50 controls per gene, fixed seeds). One scalar target per program.
- **(b) Expression vector:** the receiver's log1p CP10k expression over the **declared gene set** = the union of its type's program genes (myeloid: C1QA, C1QB, C1QC, CD163, MRC1, SELENOP, SPP1, GPNMB, CD9, FABP5, LGALS3 = 11 genes; tumour: the 18 classical + basal_like genes). Multi-output; the metric pools SSE and SST over genes, and a per-gene table is also reported.
- **Own inputs never include the receiver's own expression.** Supervised models see only own covariates (below) and neighbourhood features. The feature registry carries provenance (`self_covariate`, `tile`, `neighbour`) and gene names; `features.assert_no_target_leak` raises if any self-sourced column involves a target gene or target program, and additionally raises if any self-sourced column has |Pearson r| > 0.999 with any target column (catches a planted copy even with wrong tags). Unit tests plant the target and require the error.

## Own covariates (baseline (2); "section-free")
- log1p total counts; log1p detected genes; local density = number of cells within 30 um; mean 6-NN edge length.
- Zonation control: leave-one-out mean of each target (programs and declared genes) over *same-type* cells in the cell's 200 um tile (grid anchored at the section origin), with an indicator for an empty tile (then the training-set mean is imputed). These are other cells' values, never the cell's own.

## Neighbourhood features (for aggregated-feature models)
- Composition: fraction of each annotation class among 6-NN neighbours; count of each class within 30 um.
- Neighbour programs: for each typed class (tumour, CAF, T, myeloid) the mean over its neighbours (6-NN, and 30 um) of the 7 program scores and 4 lineage scores; filled with the section mean of that class when absent (absence is carried by composition).
- Neighbour genes: for each of tumour, CAF, T, myeloid and "other", the mean over 6-NN neighbours of that class of the 539-gene log expression (G). This block is the one that can carry cross-type gene-level signal.
- Feature sets: `programs` = own covariates + composition + neighbour programs (the chain Phase 2 set, used for the reproduction check); `genes` = `programs` + neighbour genes.

## Metric
- **Pooled held-out R2** = 1 - sum SSE / sum SST over held-out cells, where SST uses the *training-fold mean* of the target (so the mean predictor scores exactly 0). Per-tumour (SSE, SST) pairs are stored; the 95 % CI resamples tumours (2,000 draws). Differences between models use a paired tumour bootstrap of the same pairs.
- Secondary: per-tumour R2, and Pearson r on held-out cells.

## Splits
- **Outer:** grouped 5-fold over the 13 tumours (all sections of a tumour in one fold). Assignment is deterministic: tumours sorted by myeloid cell count (descending, ties by name), each assigned greedily to the fold with the fewest myeloid cells so far. Written to `model_dev/results/phase1/folds.csv` and reused by every phase.
- **Inner validation:** for outer fold k, the inner validation set is training fold (k + 1) mod 5; inner training is the other three. Hyperparameters and early stopping use only the inner validation set. Baselines with a grid are then refitted on all four training folds with the chosen setting; neural models are not refitted (they keep the early-stopped weights), and this is stated with the results.
- **Reproduction check only:** leave-one-tumour-out (13 folds) is used once, in Phase 2, to reproduce the chain Phase 2 ridge numbers.
- **Cell-count floor:** each outer test fold must contain >= 500 receiver cells and >= 2 tumours for a receiver to be evaluated. T cells fail the primary-target criterion by declaration (above).

## Primary endpoint and success criteria
- **Primary endpoint:** pooled held-out R2 for `myeloid:C1Q_TAM` (program score) under the grouped 5-fold split, and the **neighbourhood gain** = R2(model) - R2(own-covariates baseline) with its tumour-bootstrap 95 % CI.
- "Neighbours predict beyond own covariates" is supported for a target when the gain CI excludes 0.
- A graph model "beats the strongest baseline" only when R2(graph) - R2(best baseline) has a CI excluding 0 (Phase 3).
- The **cross-type contribution** (Phase 3) = R2(full model) - R2(same-type-neighbours-only model), with CI; "other cell types carry signal" is supported only when that CI excludes 0.

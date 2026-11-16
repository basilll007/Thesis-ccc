# Phase 6 plan: pretraining feasibility (small, exploratory)

Committed before `pretrain.py` is run.

## Encoder and objective
- Self-supervised masked-expression objective on the development set: for each cell, mask a random 20 % of its G genes (set to 0 with a mask indicator), and predict the masked log expression from (a) the cell's remaining genes and covariates and (b) its 2-hop neighbourhood tree (same padded-tree encoder as Phase 3, SAGE-style mean aggregation, hidden 128). Loss: MSE on masked positions. All cell types are used.
- **Fold discipline:** the target in Phases 1-3 is expression, so pretraining on test-fold cells would leak their targets. The encoder is therefore pretrained separately inside each outer fold on the inner-training cells only (3 folds), with early stopping on the inner-validation masked loss. 1 seed per fold (exploratory; logged), 15 epochs maximum, batch 512.

## Probe
- Linear probe: ridge on the frozen 128-dim seed embedding computed with the receiver's **target genes masked** (its other G genes remain visible, as in the pretraining input). Alpha from the Phase 2 grid on the inner validation fold. Same targets, folds and metric as Phase 2.
- Comparison: pooled R2 of the probe versus (i) the own-covariates baseline, (ii) the strongest Phase 2 baseline, (iii) the Phase 3 graph model, with paired tumour-bootstrap CIs. Note that the probe sees the cell's own non-target genes, which the supervised models do not; a second probe with the cell's own genes fully masked (neighbourhood-only embedding) is reported for a like-for-like comparison.

## Data-scale estimate
- Learning curve: the supervised Phase 3 model and the probe trained with 3, 6, 9 and all 12-13 development tumours (nested subsets, fixed seed), evaluated on the same held-out folds; plot R2 against the number of training cells. The extrapolation is descriptive (log-linear fit) and labelled as a rough estimate.

## Decision
"Pretraining helps" only if the like-for-like (neighbourhood-only) probe beats the Phase 3 supervised model with a CI excluding 0. The result is not called a foundation model in any case.

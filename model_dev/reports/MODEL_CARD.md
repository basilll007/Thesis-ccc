# Model card: neighbourhood models of cell state (model_dev v1, baseline stage)

## Models
- **Selected graph model:** SAGE-style 2-hop padded-tree message passing in PyTorch.
  - Architecture: hidden 128, class embedding 16, dropout 0.1, AdamW lr 3e-4, weight decay 1e-4, batch 512, early stopping with patience 8. Neighbourhood: 12 + 12 neighbours of the symmetrised 6-NN graph.
  - Selection: chosen over an attention variant by inner-validation MSE on the primary receiver (`results/phase3/architecture_selection.json`).
- **Baselines:** mean predictor; ridge on own covariates; ridge and XGBoost-CUDA on aggregated neighbourhood features; MLP. See `plans/phase2_baselines.md`.
- **Checkpoints:** `model_cache/checkpoints/` (git-ignored; SHA-256 manifest in `results/phase7/checkpoint_manifest.csv`).

## Intended use
Research only. The models quantify how much of a cell's measured expression state (macrophage C1Q_TAM / SPP1_TAM programs, tumour classical / basal-like programs, or the declared gene vectors) can be predicted from its spatial neighbourhood in CosMx PDAC sections, and how much of that comes from other cell types.
- **Not for** clinical use, patient stratification, or inferring causal signalling.
- **Outputs are predictive associations.** In-silico perturbations are model-based hypotheses.

## Inputs and outputs
- **Inputs, receiver cell:**
  - its class;
  - log1p total counts, log1p detected genes, local density (cells within 30 um) and mean 6-NN edge length;
  - leave-one-out 200 um tile means of the targets over same-type cells.
  - **Never** its own expression.
- **Inputs, neighbour cells:** class, counts covariates, 539-gene log expression (standardised with training statistics), and edge lengths.
- **Outputs:** the receiver's program scores and declared-gene log1p CP10k expression.

## Training data
Mendeley PDAC-SMI CosMx (CC BY 4.0, as recorded in the repository): 16 sections, 13 tumours, author labels. Receivers: 95,941 myeloid cells, and tumour cells capped at 8,000 per section (102,546). See `reports/DATA_CARD.md`.

## Metrics (pooled held-out R2, 95 % tumour CI; `MODEL_REPORT_v1.md`)
- **Internal, grouped 5-fold:**
  - myeloid C1Q_TAM 0.336 [0.220, 0.432], not better than XGBoost on aggregated neighbour genes (+0.0025 [-0.013, 0.015]);
  - myeloid gene vector 0.289 [0.214, 0.348], better than the strongest baseline (+0.013 [0.005, 0.020]);
  - tumour classical 0.597 [0.444, 0.703], worse than the baseline (-0.097 [-0.172, -0.053]).
- **External (GSE310352, single touch, 6 units of unverified independence):**
  - myeloid C1Q_TAM 0.324 [0.239, 0.452];
  - myeloid gene vector 0.266 [0.203, 0.339];
  - tumour classical 0.466 [0.351, 0.593].

## Limitations and failure cases
- **Does not beat simple aggregated baselines** on the program targets. Prefer ridge or XGBoost on neighbour programs and genes for those.
- **Poor on tumour basal_like:** 0.281 internal, with unstable seeds.
- **Weakest on GSE310352 slide 5b2** (C1Q_TAM R2 0.06).
- **Sensitive to annotation shift** (author vs marker labels). The tumour own-covariate baseline collapses externally.
- **Spillover.** Neighbour co-variation may partly reflect transcript spillover between adjacent segmented cells.
- **Scope.** Trained on 13 tumours from one cohort, so generalisation to other platforms, panels or tissues is untested.
- **T cells** are not supported as receivers (too few per fold).

## Ethical and data notes
Public de-identified datasets only. No patient-level attributes are used or predicted.

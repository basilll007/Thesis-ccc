# v2 Phase 1 plan: segmentation-spillover control for the CAF -> macrophage signal

Committed before any Phase 1 result is computed. Branch `model-v2` (from origin/main 71ebd43). Development set only (Mendeley PDAC-SMI CosMx, 16 sections, 13 tumours, author labels). The external set GSE310352 is **not used** in Phase 1.

Reused unchanged from v1:
- the folds (`results/phase1/folds.csv`; grouped 5-fold; inner validation = fold (k+1) mod 5);
- the feature cache (`model_cache/features/development__myeloid.npz`) and the feature registry with provenance;
- the leakage guard (`features.assert_no_target_leak`).

Expression, not mutation; every result is a predictive association. The tumour is the unit for every CI (paired tumour bootstrap, 2,000 draws, seed 42).

## Question
Is the CAF -> macrophage co-variation explained by transcripts of neighbouring CAFs being assigned to the macrophage (segmentation spillover)? In v1, CAF neighbours were the largest other-type contributor to macrophage-state prediction.

## Pre-plan descriptives (no target values read; `results/v2_phase1/preplan/`)
- Median equivalent cell diameter over all development cells, from the segmentation area: **D = 12.14 um** (`preplan_summary.json`).
- Distance from each myeloid receiver centroid to the nearest CAF centroid: median 13.4 um, 90th percentile 35.6 um.
- Bins (fixed from D, not from results):
  - **touching**: d <= D (40,626 receivers, 42 %);
  - **near**: D < d <= 2D (37,235, 39 %);
  - **separated**: d > 2D = 24.3 um (18,080, 19 %).
- With centroid-to-edge about D/2, d > 2D means a gap of more than about one cell diameter. Every tumour has >= 194 separated receivers (`myeloid_distance_bins_by_tumour.csv`).
- Protein channels available per cell: Mean.CD45, Mean.PanCK, Mean.CD3, Mean.MembraneStain, Mean.DAPI, and cell area.
  - **No CD68 channel.**
  - PanCK does not separate macrophages from CAFs (both are PanCK-low), so the identity check uses CD45 only.

## Shared definitions
- **Expression units:** log1p CP10k from the Phase 0 datasets (`X_log`, 539-gene set G). Linear CP10k = expm1.
- **Isolated macrophages:** myeloid receivers in the separated bin (nearest CAF > 2D).
- **Gene statistics** (all genes in G):
  - per tumour, mean linear CP10k and detection fraction (counts > 0) in (i) all CAFs, (ii) isolated macrophages and (iii) touching macrophages;
  - the across-tumour median is the gene's summary value;
  - ratio_CAF/iso(g) = median CAF mean / median isolated-macrophage mean.
- **Contamination markers M** (part b):
  - genes g in G, excluding the 11 macrophage target genes and the myeloid lineage genes (CD68, CD14, LYZ, CSF1R, FCGR3A);
  - CAF detection fraction >= 0.25 and ratio_CAF/iso >= 10;
  - if more than 20 genes pass, keep the 20 with the highest CAF mean;
  - if fewer than 3 pass, the ratio threshold becomes 5 (declared fallback); if still fewer than 3, part (b) and (c-i) are BLOCKED.
- **Contamination score** s_c (per cell) = log1p(10^4 x sum of raw counts of M / total counts).
- **Spillover-free macrophage genes S** (part d):
  - genes in G with isolated-macrophage detection fraction >= 0.25 and ratio_iso/CAF >= 10 (isolated-macrophage mean at least 10x the CAF mean);
  - target genes are allowed in S;
  - keep the top 20 by isolated-macrophage mean if more pass;
  - fewer than 3 -> ratio-5 fallback; still fewer -> (d) BLOCKED.
- **Target-gene removal** (part c-ii):
  - a macrophage target gene is removed if ratio_CAF/iso >= 0.5, i.e. CAF expression is at least half the isolated-macrophage level;
  - adjusted targets: C1Q_TAM_adj and SPP1_TAM_adj = mean log1p CP10k over the retained genes of each program; gene_vector_adj = the retained genes among the 11;
  - a program with fewer than 2 retained genes is "not estimable".
- **Gene-mean reference targets:**
  - C1Q_TAM_mean and SPP1_TAM_mean = mean log1p CP10k over all genes of the program;
  - the adjusted targets use the same formula, so adjusted versus unadjusted is like for like;
  - the original v1 program scores (scanpy-style) are also kept as targets.
- **Set targets** (part d): S_mean and M_mean = mean log1p CP10k over the genes of S or M in the receiver.
- **Tile block for new targets:** leave-one-out 200 um tile mean over same-type (myeloid) cells of each new target, plus an empty indicator, built with the v1 function `features.tile_loo_mean`.
- **New neighbour block `ncaf`:**
  - mean log1p CP10k over G of the 3 nearest CAFs (any distance, same section), plus two distances (nearest CAF; mean of the 3 nearest);
  - needed because the v1 6-NN CAF block is empty (section-mean filled) for most separated receivers.
- **Phase 1 feature set "full"** = v1 `genes` set (own, tile, composition, nb_programs, nb_genes) + `ncaf`.
- **CAF expression columns** = nb_programs CAF (6-NN and 30 um), nb_genes CAF (6-NN) and ncaf expression.
  - "Minus CAF expression" removes these and keeps composition and CAF distances: who is nearby is kept; what CAFs express is removed.
- **CAF-expression contribution Q** = R2(full) - R2(full minus CAF expression).
- **Neighbour gain** = R2(full) - R2(own), where own = ridge on own + tile (v1 reference).

## Models and tuning
- **XGBoost-CUDA** (primary):
  - the v1 configuration selected for the myeloid genes set: learning_rate 0.05, max_leaves 15, lossguide, min_child_weight 50;
  - the number of rounds (<= 500) is chosen on the inner-validation RMSE, then the model is refitted on all four training folds;
  - <= 60,000 training rows (seeded subsample, as v1);
  - the full 4-point grid is not repeated in Phase 1, to fit the compute budget (each fit takes ~35 s on 3,370 columns). Phase 2 uses the full grid.
- **Ridge** (GPU closed form): alpha in {0.1, 1, 10, 100, 1000} chosen on inner validation, then refitted on all training folds. Reported alongside, and used for multi-gene vector targets.
- **Preprocessing and safety:**
  - features are standardised with training-fold statistics (as v1);
  - the leakage guard runs on every design matrix;
  - s_c is tagged `self_covariate` with genes = M (M is disjoint from every target by construction).
- **Outputs per fit:** held-out predictions per (task, fold) are saved (resumable) in `model_cache/v2/phase1/`. All metrics below are computed from them.

## Analyses
**a. Distance dependence.**
- Reference models (no adjustment). Targets: C1Q_TAM score (primary), SPP1_TAM score, S_mean and M_mean.
- Fit on all receivers; evaluate the test cells per bin.
- Report per bin: neighbour gain and Q, using pooled R2 within the bin (SST around the training-fold mean) with tumour-bootstrap CIs within the bin. Also report the ridge versions.

**b. Contamination-marker test.**
- Report M with its statistics.
- Among receivers, report mean marker CP10k (sum over M) by bin and by number of CAFs among the 6 nearest neighbours (0, 1, 2, >= 3), with tumour-bootstrap CIs.
- **Spillover-rate estimate** per bin = (marker CP10k in bin - marker CP10k in isolated macrophages) / (marker CP10k in CAFs - marker CP10k in isolated macrophages). CI from a tumour bootstrap of the ratio of pooled means.
- Correlation (Spearman, per tumour; median and range) between a touching receiver's s_c and the mean marker CP10k of its 3 nearest CAFs.

**c. Spillover-adjusted model.**
- Adjustments: (i) s_c added to the own block of every model (own and full); (ii) target-gene removal. (i) + (ii) combined is "the spillover-adjusted model".
- Report neighbour gain and Q for:
  - C1Q_TAM score with (i);
  - C1Q_TAM_mean without adjustment;
  - C1Q_TAM_adj with (ii) only;
  - C1Q_TAM_adj with (i) + (ii);
  - the same four for SPP1_TAM;
  - ridge gene_vector_adj with (i) + (ii).

**d. Negative-control targets.**
- Q for S_mean (spillover-free macrophage genes) and for M_mean (CAF-high genes in macrophages: the pure spillover expectation), using the reference features without s_c.
- Ridge per-gene Q for every gene in S and M (descriptive).

**e. Protein identity check.**
- Protein-confirmed macrophages: author-labelled myeloid cells with Mean.CD45 above the 75th percentile of Mean.CD45 among CAFs of the same section (primary; pre-plan check: ~62 % retained). Sensitivity: above the 90th percentile (~40 % retained).
- Reference models for C1Q_TAM are refitted with training and test restricted to confirmed cells; report neighbour gain and Q.
- Descriptive only: s_c by CD45 status, and s_c by cell-area quartile.

**f. Hypothesis genes.** For SELENOP, CD163, C1QA and SPP1 (and the other 7 target genes), report:
- mean CP10k and detection in CAFs, isolated macrophages and touching macrophages (across-tumour median and range);
- ratio_CAF/iso with a tumour-bootstrap CI.

## Decision rule (fixed now)
Primary quantity: Q for macrophage C1Q_TAM, XGBoost. "CI > 0" means the paired tumour-bootstrap 95 % CI lower bound is above 0.
- **Check (a)** passes if Q (reference model, C1Q_TAM score) has CI > 0 in the near bin or the separated bin, i.e. it is not confined to the touching bin. The separated-bin result alone is reported as the stricter reading.
- **Check (c)** passes if Q on the spillover-adjusted model ((i) + (ii); target C1Q_TAM_adj) has CI > 0. If C1Q_TAM_adj is not estimable, ridge gene_vector_adj is used; if that is not estimable either, (c) is "not estimable".
- **Check (d)** passes if Q for S_mean has CI > 0.
- **Verdict:**
  - **"survives spillover control"** if (a), (c) and (d) all pass;
  - **"explained by spillover"** if (c) or (d) fails AND the distance pattern matches spillover (touching-bin Q CI > 0 while the near- and separated-bin Q CIs include 0);
  - **"inconclusive"** otherwise, naming every failed check;
  - if the reference Q (all receivers) itself has a CI including 0, the verdict is "inconclusive: no reference CAF-expression signal in the aggregated model".
- Ridge, SPP1_TAM, the protein check and part (b) are reported but do not change the verdict.

## Outputs
- `model_dev/results/v2_phase1/`:
  - gene statistics; M and S lists; adjusted-target definitions;
  - per-bin table; contamination table; protein-check table; hypothesis-gene table;
  - per-task summary (R2, gains, Q with CIs);
  - decision JSON, figure, run logs.
- Large per-cell arrays go in `model_cache/v2/`, with checksums in `results/v2_artefacts.csv` (the v1 registry is not modified).

## Compute and safety
- GPU for all model fits; CPU workers capped at 4.
- Free RAM is logged before the phase, and the phase does not start below 8 GB.
- Outputs are resumable per (task, fold).
- If the phase fails twice with the same error, it is marked BLOCKED.

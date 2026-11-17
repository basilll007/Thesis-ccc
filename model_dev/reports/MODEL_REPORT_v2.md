# Model report v2: spillover control, cross-type share on the best model, and survival

Branch `model-v2` (from origin/main 71ebd43).

**Data and scope.**
- Development set: Mendeley PDAC-SMI CosMx (16 sections, 13 tumours).
- Survival cohorts: TCGA (strict histology, 141 patients / 81 events) and CPTAC (129 / 72).
- The external CosMx set GSE310352 was **not used** in v2.

**Provenance.**
- Every number below is read from a file generated in this work:
  - `model_dev/results/v2_phase1/`, `model_dev/results/v2_phase2/` and `results/translational/phase3/20261117T024638285521Z/`;
  - figure index with checksums: `model_dev/results/figures/v2_figure_index.csv`.
- Plans were committed before each phase: `model_dev/plans/v2_phase1_spillover.md`, `v2_phase2_crosstype.md` (+ Clarification 1) and `v2_phase3_survival.md`. Deviations and failures are in `reports/PROGRESS_LOG_v2.md`.

**How to read the numbers.**
- These are expression data, and every result is a predictive or prognostic association, not a causal, mechanistic or clinical claim.
- The tumour (patient) is the unit of every CI: paired tumour bootstrap with 2,000 draws.
- R2 values are pooled held-out R2 over grouped 5-fold cross-validation.

## Three answers
**1. Does the CAF -> macrophage signal survive spillover control? No: by the pre-declared rule it is "explained by spillover", with all three checks (a, c, d) failing.**
- **Reference signal.** The CAF-expression contribution Q (R2 lost when the CAF expression columns are removed, keeping who is nearby) for macrophage C1Q_TAM with XGBoost is +0.004 [0.0004, 0.007], about a tenth of the neighbour gain (+0.038 [0.017, 0.055]).
- **Check (a), distance bins.** Q is confined to macrophages that touch a CAF:
  - touching (nearest CAF <= 12.1 um, one median cell diameter): +0.008 [0.004, 0.013];
  - near (12.1-24.3 um): +0.004 [-0.001, 0.008];
  - separated (> 24.3 um): -0.004 [-0.010, -0.0003].
  - The touching-only pattern is what spillover predicts. FAIL.
- **Check (c), spillover-adjusted model.** Adding a per-cell contamination covariate and removing CAF-shared genes (SELENOP, CD9, LGALS3) from the target gives +0.005 [-0.001, 0.013]. The CI includes 0, a narrow fail. FAIL.
- **Check (d), negative controls.**
  - Spillover-free macrophage genes S (HLA-DRB1, HLA-DRA, SPP1, CD163, ITGAX; at least 10x higher in isolated macrophages than in CAFs): Q = +0.001 [-0.003, 0.005]. FAIL.
  - For comparison, CAF marker genes measured inside macrophages (the pure-spillover expectation): +0.009 [0.002, 0.012].
- **Contamination is large.** Macrophages touching a CAF carry an extra **23 % [17 %, 26 %]** of the CAF marker level over isolated macrophages (near: 16 % [11 %, 18 %]).
- **The protein identity check does not rescue the signal, and does not refute spillover either.** In CD45-confirmed macrophages Q stays small and positive (+0.005 [0.003, 0.007]; stricter threshold +0.007 [0.003, 0.010]). CD45 confirms identity, but it cannot remove transcripts carried into true macrophages.
- **What survives** (not part of the decision):
  - The neighbourhood gain itself survives at every distance (separated bin +0.052 [0.002, 0.074]).
  - Post hoc and exploratory: touching macrophages have *more* C1QA/B/C, CD163, MRC1 and SELENOP than pure mixing with average CAF RNA would produce.

**2. How large is the cross-type share on the best model? About zero for C1Q_TAM on XGBoost; a minority share on ridge that is carried by CAF expression.**
- **XGBoost (genes set), macrophage C1Q_TAM:**
  - neighbour gain +0.034 [0.019, 0.048], reproducing v1;
  - a same-type-only model reaches 0.339, versus 0.334 for all neighbours;
  - cross-type contribution -0.005 [-0.010, 0.0002]; share -15 % [-33 %, 1 %].
- **SPP1_TAM (XGBoost):** cross-type +0.003 [-0.002, 0.008].
- **Ridge (programs set):** cross-type +0.009 [0.006, 0.012] for C1Q_TAM (share 28 % [17 %, 39 %]), mostly from CAF expression (unique CAF contribution +0.006 [0.003, 0.009]).
- **Gene vector (ridge):** +0.002 [0.001, 0.004] (share 19 % [12 %, 30 %]).
- **Tumour cells:** no cross-type share (classical XGBoost -0.001 [-0.005, 0.003]).
- **Comparison with the graph model:** every aggregated estimate is **smaller** than the graph-model estimate (+0.019 [0.004, 0.028] for C1Q_TAM). By the pre-declared criterion none "agrees".
- **Spillover-adjusted re-run.** The C1Q_TAM share stays at about zero on XGBoost (-0.007 [-0.015, 0.007]). Ridge keeps +0.005 [0.002, 0.008], about the same as its unadjusted like-for-like reference (+0.006). For the SPP1 program, XGBoost finds a small share (+0.006 [0.004, 0.009], 14 %), mostly from *who is nearby* rather than from neighbour expression.

**3. Do any signatures relate to survival after adjustment, and does that hold in CPTAC? No signature passes the pre-declared rule (adjusted TCGA BH q < 0.05 over 9 signatures: 0 of 9).**
- **Closest: hypoxia.** HR per SD 1.32 [1.07, 1.62], p = 0.009, q = 0.078, PH test p = 0.14. It points the same way in CPTAC, HR 1.43 [1.10, 1.86], but it did not pass the primary rule in TCGA.
- **The two new macrophage signatures and their difference are null in both cohorts:**
  - C1Q-TAM: TCGA 0.93 [0.74, 1.17], CPTAC 1.00 [0.76, 1.30];
  - SPP1-TAM: 1.12 [0.85, 1.46] and 1.27 [0.98, 1.64];
  - C1Q minus SPP1: 0.88 [0.71, 1.11] and 0.82 [0.61, 1.10].
- **Power.** With 75 events (adjusted model), the minimum detectable HR per SD at 80 % power is about 1.38-1.40 (alpha 0.05), and about 1.52-1.55 at alpha 0.05/9. Effects of the size seen for hypoxia are below what this cohort can reliably detect.

## Phase 1 results: spillover control (development set; XGBoost unless stated)
Files: `results/v2_phase1/phase1_summary.csv` (all tasks, models, bins), `phase1_decision.json`, `b_*.csv`, `e_*.csv`, `f_*.csv`. Figure: `results/figures/v2_fig1_spillover.png`, with every estimate, its CI and the PASS/FAIL of each check printed on the figure.

**Gene sets, defined from the data by the pre-declared rules** (`gene_sets_and_adjusted_targets.json`):
- **Contamination markers M (20):**
  - COL1A1, COL3A1, COL1A2, IGFBP7, FN1, TAGLN, BGN, COL6A2, COL6A3, DCN, LUM, CALD1, COL6A1, TPM1, MMP2, THBS2, MGP, COL5A1, COL4A2, COL12A1.
  - Only 2 genes reached the declared 10x CAF/isolated-macrophage ratio, so the declared fallback ratio of 5 was used.
  - Even isolated macrophages carry about 10 % of CAF collagen levels (COL1A1: 49 vs 449 CP10k), so the 23 % estimate is an *excess* over a contamination floor.
- **Spillover-free genes S (5):** HLA-DRB1, HLA-DRA, SPP1, CD163, ITGAX.
- **Target genes removed (CAF/isolated >= 0.5):** SELENOP (0.81), LGALS3 (0.71), CD9 (0.61). Retained: C1Q_TAM_adj = C1QA, C1QB, C1QC, CD163, MRC1; SPP1_TAM_adj = SPP1, GPNMB, FABP5.

| Target (XGBoost; ridge own+tile as reference) | R2 own | R2 full | neighbour gain [CI] | Q, CAF-expression contribution [CI] |
|---|---|---|---|---|
| C1Q_TAM score (reference) | 0.299 | 0.337 | +0.038 [0.017, 0.055] | +0.004 [0.0004, 0.007] |
| C1Q_TAM score + contamination covariate s_c | 0.301 | 0.350 | +0.049 [0.031, 0.066] | +0.006 [0.001, 0.010] |
| C1Q_TAM gene mean (unadjusted, like for like) | 0.413 | 0.438 | +0.025 [0.002, 0.046] | +0.002 [-0.0001, 0.004] |
| C1Q_TAM_adj (genes removed) | 0.393 | 0.412 | +0.019 [0.001, 0.037] | +0.001 [-0.001, 0.004] |
| **C1Q_TAM_adj + s_c (spillover-adjusted model, check c)** | 0.395 | 0.424 | +0.029 [0.005, 0.054] | **+0.005 [-0.001, 0.013]** |
| SPP1_TAM score (reference) | 0.336 | 0.384 | +0.048 [0.032, 0.068] | -0.002 [-0.005, 0.002] |
| SPP1_TAM_adj + s_c | 0.501 | 0.545 | +0.044 [0.032, 0.059] | +0.003 [0.002, 0.004] |
| **S genes, spillover-free (check d)** | 0.308 | 0.336 | +0.028 [0.023, 0.035] | **+0.001 [-0.003, 0.005]** |
| M genes inside macrophages (pure spillover) | 0.434 | 0.547 | +0.113 [0.073, 0.142] | +0.009 [0.002, 0.012] |
| C1Q_TAM, CD45 > CAF 75th pct (64 % of cells) | 0.308 | 0.349 | +0.041 [0.016, 0.059] | +0.005 [0.003, 0.007] |
| C1Q_TAM, CD45 > CAF 90th pct (39 %) | 0.294 | 0.331 | +0.038 [0.005, 0.062] | +0.007 [0.003, 0.010] |

**By distance to the nearest CAF** (bins fixed before results from the median cell diameter, 12.14 um):

| Bin (receivers) | C1Q_TAM neighbour gain | C1Q_TAM Q | S genes Q | M genes Q | CAF-marker CP10k in macrophages |
|---|---|---|---|---|---|
| touching (40,626) | +0.032 [0.017, 0.048] | **+0.008 [0.004, 0.013]** | +0.007 [0.005, 0.010] | +0.016 [0.008, 0.020] | 845 [687, 937] |
| near (37,235) | +0.037 [0.019, 0.053] | +0.004 [-0.001, 0.008] | +0.004 [0.002, 0.005] | +0.006 [-0.0004, 0.009] | 687 [565, 754] |
| separated (18,080) | +0.052 [0.002, 0.074] | -0.004 [-0.010, -0.0003] | -0.013 [-0.021, -0.0003] | -0.003 [-0.006, -0.0001] | 334 [292, 392] |

CAFs themselves: 2,543 CP10k.

**Readings.**
- The CAF-expression contribution follows the contamination gradient: largest where cells touch, gone or negative beyond one cell diameter.
- It appears even for the spillover-free S genes in the touching bin (+0.007). This is consistent with **dilution**: extra CAF transcripts lower a macrophage's normalised expression of its own genes, which is also a segmentation artefact.
- **The ridge rows overfit on the 3,372-column Phase 1 set** (negative neighbour gains; v1 saw the same for ridge on the genes set). They are reported in the CSV but do not carry the reading, as declared.

**Hypothesis genes (part f; `f_target_genes_caf_vs_macrophage.csv`).** Median over tumours of mean CP10k in CAFs vs isolated macrophages, with the ratio CAF/isolated [tumour-bootstrap CI]:

| Gene | CAFs | isolated macrophages | ratio CAF/isolated |
|---|---|---|---|
| SELENOP | 6.7 | 8.3 | 0.81 [0.31, 1.36] |
| CD163 | 3.9 | 39.0 | 0.10 [0.07, 0.13] |
| C1QA | 6.4 | 47.7 | 0.13 [0.11, 0.19] |
| SPP1 | 8.8 | 96.8 | 0.09 [0.03, 0.24] |

- **SELENOP** is the only hypothesis gene for which carry-over of CAF RNA is a credible explanation of a macrophage-level shift.
- **Post hoc, exploratory** (`f_posthoc_mixing_expectation.csv`, not in the plan): if touching macrophages were pure mixtures (23 % CAF RNA), they would have C1QA 38, CD163 31 and SELENOP 7.9 CP10k. They have 46, 47 and 19.5.
  - These genes are *above* the mixing expectation, so the state difference near CAFs is not explained by average-CAF carry-over. Adjacent CAFs could differ from the average CAF, so this is not proof.
  - SPP1 (expected 76, observed 81) is within 10 % of the mixing expectation, so its lower level near CAFs is consistent with dilution.

## Phase 2 results: cross-type share on the aggregated baselines
Files: `results/v2_phase2/phase2_quantities.csv` (every quantity, model, target), `agreement_with_graph_estimate.csv`, `reproduction_check_vs_v1.csv`. Figure: `results/figures/v2_fig2_crosstype.png`.
- **Reproduction:** the `all` variant reproduces v1 Phase 2 in 10/10 comparisons (ridge |diff| < 1e-7; XGBoost |diff| <= 0.0008; tolerance 0.005), in the rebuilt environment (see Compute).

| Receiver / target | model | neighbour gain | cross-type (all - same) | cross-type share | who is nearby (other-type composition) | what they express (other-type expression) | unique CAF |
|---|---|---|---|---|---|---|---|
| macrophage C1Q_TAM | XGBoost | +0.034 [0.019, 0.048] | -0.005 [-0.010, 0.0002] | -15 % [-33 %, 1 %] | +0.0004 [-0.005, 0.003] | -0.006 [-0.010, 0.001] | -0.001 [-0.003, 0.003] |
| macrophage C1Q_TAM | ridge | +0.033 [0.023, 0.043] | +0.009 [0.006, 0.012] | 28 % [17 %, 39 %] | +0.002 [-0.002, 0.006] | +0.007 [0.004, 0.009] | +0.006 [0.003, 0.009] |
| macrophage SPP1_TAM | XGBoost | +0.050 [0.034, 0.069] | +0.003 [-0.002, 0.008] | 6 % [-4 %, 16 %] | +0.006 [0.004, 0.009] | -0.003 [-0.008, 0.003] | -0.001 [-0.005, 0.002] |
| macrophage SPP1_TAM | ridge | +0.048 [0.037, 0.057] | +0.007 [0.004, 0.012] | 15 % [9 %, 27 %] | +0.003 [0.001, 0.008] | +0.004 [0.001, 0.006] | +0.002 [-0.0002, 0.004] |
| macrophage 11-gene vector | ridge | +0.012 [0.008, 0.016] | +0.002 [0.001, 0.004] | 19 % [12 %, 30 %] | +0.001 [0.0004, 0.002] | +0.001 [0.0002, 0.002] | +0.001 [0.001, 0.001] |
| tumour classical | XGBoost | +0.075 [0.053, 0.106] | -0.001 [-0.005, 0.003] | -1 % [-8 %, 4 %] | +0.001 [0.0002, 0.001] | -0.001 [-0.005, 0.002] | +0.002 [0.001, 0.004] |
| tumour classical | ridge | +0.074 [0.052, 0.105] | +0.0002 [-0.001, 0.002] | 0 % [-1 %, 2 %] | +0.0001 [-0.0002, 0.001] | 0.000 [-0.001, 0.001] | +0.001 [0.0002, 0.001] |
| tumour basal_like | XGBoost | +0.046 [0.032, 0.061] | +0.003 [-0.001, 0.008] | 6 % [-2 %, 23 %] | -0.001 [-0.002, 0.0002] | +0.004 [-0.0003, 0.010] | +0.006 [0.001, 0.012] |
| tumour basal_like | ridge | +0.041 [0.029, 0.052] | -0.003 [-0.006, -0.001] | -8 % [-15 %, -2 %] | -0.0005 [-0.001, -0.0002] | -0.003 [-0.005, -0.0003] | +0.0004 [-0.0002, 0.001] |

**Agreement with the graph-model estimate** (pre-declared: the aggregated CI contains the graph point and vice versa): none agrees, and all are **smaller**.
- C1Q_TAM: XGBoost -0.005 vs graph +0.019 [0.004, 0.028];
- SPP1_TAM: XGBoost +0.003 vs graph +0.015;
- gene vector: ridge +0.002 vs graph +0.014;
- tumour classical: XGBoost -0.001 vs graph +0.007.

**Readings.**
- **The best model gets the macrophage neighbour gain from neighbouring macrophages.** A same-type-only XGBoost model reaches the full model's R2 for C1Q_TAM and SPP1_TAM.
- **Model dependence.** The linear ridge model, which can only use averaged programs, gains a minority share from other cell types, carried by CAF *expression*. That is the same kind of signal Phase 1 attributed to spillover. The flexible model does not need it.
- **"Who is nearby" vs "what they express":** both components are within about ±0.006 R2 for every target. The only robust positive composition term is SPP1_TAM (+0.006 [0.004, 0.009], XGBoost).

**Spillover-adjusted re-run** (myeloid; contamination covariate s_c in every variant; targets with CAF-shared genes removed):

| Target | model | neighbour gain | cross-type (all - same) | cross-type share | who is nearby | what they express | unique CAF |
|---|---|---|---|---|---|---|---|
| C1Q_TAM gene mean, unadjusted reference | XGBoost | +0.023 [0.003, 0.042] | -0.003 [-0.007, 0.004] | -15 % [-132 %, 14 %] | +0.002 [-0.002, 0.004] | -0.005 [-0.009, 0.004] | -0.001 [-0.005, 0.004] |
| **C1Q_TAM_adj + s_c** | XGBoost | +0.027 [0.004, 0.052] | **-0.007 [-0.015, 0.007]** | -24 % [-226 %, 18 %] | +0.002 [-0.001, 0.005] | -0.009 [-0.018, 0.005] | +0.002 [-0.003, 0.008] |
| C1Q_TAM gene mean, unadjusted reference | ridge | +0.027 [0.018, 0.036] | +0.006 [0.002, 0.010] | 23 % [7 %, 38 %] | +0.001 [-0.002, 0.005] | +0.005 [0.001, 0.007] | +0.005 [0.002, 0.008] |
| **C1Q_TAM_adj + s_c** | ridge | +0.026 [0.017, 0.035] | **+0.005 [0.002, 0.008]** | 20 % [8 %, 38 %] | +0.002 [-0.002, 0.005] | +0.004 [0.0001, 0.006] | +0.004 [0.002, 0.007] |
| SPP1_TAM gene mean, unadjusted reference | XGBoost | +0.044 [0.033, 0.057] | +0.005 [0.001, 0.012] | 12 % [3 %, 26 %] | +0.007 [0.004, 0.012] | -0.002 [-0.006, 0.002] | +0.002 [-0.0002, 0.004] |
| **SPP1_TAM_adj + s_c** | XGBoost | +0.044 [0.031, 0.058] | **+0.006 [0.004, 0.009]** | 14 % [9 %, 20 %] | +0.005 [0.004, 0.007] | +0.001 [-0.001, 0.003] | +0.003 [0.002, 0.005] |
| **SPP1_TAM_adj + s_c** | ridge | +0.034 [0.022, 0.045] | +0.012 [0.008, 0.014] | 33 % [27 %, 43 %] | +0.003 [0.001, 0.006] | +0.009 [0.004, 0.011] | +0.003 [0.001, 0.004] |
| gene_vector_adj + s_c (8 genes) | ridge | +0.012 [0.008, 0.017] | +0.002 [0.001, 0.003] | 19 % [13 %, 30 %] | +0.001 [0.0003, 0.002] | +0.001 [-0.00002, 0.003] | +0.001 [0.0004, 0.002] |

**Readings (adjusted re-run).**
- **C1Q_TAM.** On XGBoost the cross-type share stays at about zero after spillover adjustment. On ridge it is barely changed by the adjustment (+0.006 unadjusted gene-mean, +0.005 adjusted) and is still carried by CAF expression.
  - So the ridge share is **not removed** by the contamination covariate and gene removal. A linear model on averaged programs keeps using CAF information that a flexible model does not need.
  - This is a model-dependence result, not evidence of signalling: Phase 1 showed that the CAF-expression signal is confined to touching cells.
- **SPP1_TAM.** The gene-mean and adjusted targets show a small cross-type share on XGBoost, +0.005 to +0.006 (12-14 %).
  - It is carried mainly by *who is nearby* (other-type composition, +0.005 to +0.007), not by what the neighbours express.
  - Where macrophages sit (near CAFs, near tumour) is weakly informative of their SPP1 program beyond their macrophage neighbours.

## Phase 3 results: survival (prepared plan executed; amendment A1-A5)
Files: `results/translational/phase3/20261117T024638285521Z/`:
- `cox_models.csv` (108 models, all estimable);
- `PH_diagnostics.csv`;
- `external_validation_cptac.csv`;
- `signature_correlations.csv`;
- `gene_coverage.csv`;
- `power_min_detectable_hr.csv`;
- KM files;
- `MANIFEST_v2.csv`.

Figures: `model_dev/results/figures/v2_fig3_survival.png` (forest plot with HR [CI], q and PH p printed for every estimate, plus the CPTAC replication table) and the composite KM figure `signature_survival_v2.png` in the run folder.

**Setup.**
- **Number of tests:** 9 signatures x 3 model families per cohort. The primary family is the purity-adjusted model in strict TCGA (9 tests, BH). 27 models per cohort, 108 in total.
- **PH diagnostic** (validated before use): null rejection 0.04, time-varying detection 1.00.
- **Coverage:** SELENOP is absent from the TCGA matrix (C1Q signature 5 of 6 genes; rule met); TRAC is absent in both cohorts (T-cell proxy 8 of 9). All other signatures are complete.

| Signature (adjusted: age, sex, stage, purity) | TCGA HR per SD [95 % CI] | TCGA q (9 tests) | PH p | CPTAC HR [95 % CI] | same direction | holds in CPTAC* |
|---|---|---|---|---|---|---|
| Hypoxia | 1.32 [1.07, 1.62] | 0.078 | 0.14 | 1.43 [1.10, 1.86] | yes | yes |
| CAF | 0.98 [0.78, 1.23] | 0.850 | 0.31 | 1.24 [0.95, 1.60] | no | no |
| EMT | 1.08 [0.84, 1.39] | 0.614 | 0.93 | 1.23 [0.97, 1.57] | yes | no |
| Interferon | 1.21 [0.95, 1.52] | 0.351 | 0.98 | 1.03 [0.80, 1.32] | yes | no |
| Immune infiltration | 0.81 [0.63, 1.03] | 0.351 | 0.13 | 0.75 [0.57, 0.99] | yes | yes |
| Immune exclusion proxy | 1.10 [0.88, 1.38] | 0.614 | 0.88 | 1.42 [1.08, 1.85] | yes | yes |
| Macrophage C1Q-TAM | 0.93 [0.74, 1.17] | 0.614 | 0.47 | 1.00 [0.76, 1.30] | yes | no |
| Macrophage SPP1-TAM | 1.12 [0.85, 1.46] | 0.614 | 0.14 | 1.27 [0.98, 1.64] | yes | no |
| Macrophage C1Q minus SPP1 | 0.88 [0.71, 1.11] | 0.614 | 0.12 | 0.82 [0.61, 1.10] | yes | no |

\*holds = same direction AND the CPTAC CI excludes 1 (amendment A3). Because no signature passed in TCGA, "holds" is descriptive only.
- **Unadjusted TCGA:** hypoxia q = 0.034, and clinical-only adjustment gives q = 0.041. Both are below 0.05, but they are not the decision family.
- **Purity:** no signature is mostly a purity proxy (largest |Spearman rho| with purity 0.28 in TCGA and 0.38 in CPTAC).
- **Macrophage C1Q signature correlations (TCGA):** strongly with immune infiltration (rho 0.72) and moderately with CAF (0.53) and EMT (0.54). In bulk it tracks immune and stromal content, not a macrophage-specific state.
- **Power** (`power_min_detectable_hr.csv`; Schoenfeld): minimum detectable HR per SD at 80 % power is 1.37 (TCGA unadjusted, 81 events) to 1.40 (adjusted, 75 events), and 1.39-1.45 in CPTAC. At alpha 0.05/9 it is 1.49-1.55 (TCGA).

## Claims register (carried forward; status changed only with evidence from this work)

| # | Claim (source) | Status before v2 | Status after v2 | Evidence |
|---|---|---|---|---|
| 1 | Neighbourhood features improve macrophage C1Q_TAM prediction (+0.034) (model v1) | Supported association | **Supported, strengthened** | Reproduced (+0.034 [0.019, 0.048]). Present at every distance from CAFs, including separated cells (+0.052). |
| 2 | A graph model is not better than aggregated baselines (v1) | Supported | Unchanged | Not retested |
| 3 | Other cell types carry a small non-zero share for macrophages (graph +0.019) (v1) | Supported within the graph model | **Not supported on the best model for C1Q_TAM; model-dependent** | XGBoost cross-type -0.005 [-0.010, 0.0002] (adjusted -0.007 [-0.015, 0.007]). Ridge +0.009 (adjusted +0.005 [0.002, 0.008]), carried by CAF expression. SPP1 program: small XGBoost share (+0.006 [0.004, 0.009]), mostly composition. |
| 4 | CAF neighbours are the largest other-type contributor to macrophage prediction (v1) | Supported association (graph permutation importance) | **Explained by segmentation spillover** (pre-declared verdict) | Checks a, c, d fail. Q is confined to touching cells. Contamination is 23 % [17 %, 26 %]. |
| 5 | Model-based hypothesis: neighbouring CAF CXCL12 -> macrophage SELENOP/CD163/C1QA up, SPP1 down (v1) | Hypothesis | **Unsupported** (its source signal is spillover-explained); SELENOP and SPP1 shifts are compatible with carry-over or dilution | Phase 1 (f) and the post-hoc mixing table |
| 6 | Macrophages touching CAFs are more C1Q/CD163/MRC1-high than average-CAF mixing predicts (new) | - | **Exploratory, post hoc** | `f_posthoc_mixing_expectation.csv` |
| 7 | Tumour-cell state is predicted from nearby tumour cells, not other types (chain v1, model v1) | Supported | **Supported, extended to the best model** | Cross-type XGBoost -0.001 [-0.005, 0.003] |
| 8 | CAF <-> C1Q_TAM adjacency lead (chain C10; FINDINGS "only consistent cross-type lead") | Exploratory association | **Weakened**: the model-based part is spillover-explained; the adjacency of cell states (claim 6) remains exploratory | Phase 1 |
| 9 | Spatial program signatures relate to survival after purity/stage adjustment (translational v1) | Not assessed | **Not supported** (0 of 9 at adjusted q < 0.05; underpowered for HR < ~1.4 per SD) | Phase 3 |
| 10 | Hypoxia program higher in metastases (v3); not robust (translational v1) | Not robust | Unchanged. New, separate observation: hypoxia HR 1.32 in TCGA (q 0.078, fails the rule) and 1.43 in CPTAC. Exploratory. | Phase 3 |
| 11 | Expression data establish signalling, causation or clinical utility | Not supported | Not supported | Outside the design |

## Limitations
1. **Spillover is estimated, not measured.**
   - The contamination markers and the distance bins are 2-D proxies. Out-of-plane overlap or ambient RNA likely sets the ~10 % floor seen even in isolated macrophages.
   - The 23 % figure is the excess over that floor.
   - There is no CD68 or other macrophage-specific protein channel; CD45 confirms immune identity only.
2. **Narrow margins.**
   - Check (c) failed with CI [-0.001, 0.013].
   - All CAF-expression contributions are <= 0.008 R2.
   - A real but tiny paracrine component cannot be excluded; it is just not distinguishable from spillover here.
3. **Thirteen tumours.** CIs are wide, and every inference rests on 13 development tumours. The external set was not used.
4. **Model dependence of the cross-type share.** XGBoost (flexible) and ridge (linear on programs) disagree, and the graph model gave a third value. The share is not a model-free quantity.
5. **Bulk survival cohorts are small** (75-81 events) and underpowered for HR < ~1.4 per SD.
   - Purity measures differ between cohorts (TCGA ABSOLUTE vs CPTAC pathology cellularity).
   - Treatment is not modelled.
   - Bulk signatures mix cell types (the C1Q signature tracks immune infiltration, rho 0.72).
6. **Phase 1 tuning and environment.**
   - Phase 1 used the v1-selected XGBoost configuration (rounds tuned); Phase 2 used the full grid.
   - The Phase 2 environment was rebuilt after the original disappeared (see Compute). Reproduction of v1 was exact for ridge and within 0.001 for XGBoost.

## Next steps (at most 5)
1. Test the remaining exploratory lead (claim 6) experimentally without imaging: a transwell vs contact co-culture of macrophages with PDAC CAFs, with sorted-macrophage readouts (`WETLAB_BRIEF_v1.md`).
2. Re-segment with a membrane stain and a nucleus-shrunk mask (or use a platform with protein-validated macrophage masks), then repeat Phase 1: if spillover drives the signal, Q should fall with tighter masks.
3. Add a macrophage-specific protein channel (CD68 or CD163) to separate identity from contamination.
4. Pool survival analyses with additional PDAC cohorts (e.g. ICGC, PACA-AU) to reach about 200+ events before testing the hypoxia lead again.
5. Report cross-type shares as model-dependent ranges (ridge, XGBoost, graph) rather than single numbers in any downstream text.

## Compute and environment
- **GPU:** RTX 5060 Laptop (8 GB), 31 GB RAM. CPU workers capped at 4. Free RAM was logged before every phase and fold process (all >= 8 GB at start).
- **Phase 1:** one process, about 50 min.
- **Phase 2 attempts:**
  - Attempt 1 (one process): GPU out-of-memory error in fold 1.
  - Attempt 2 (one fold per process): out-of-memory error twice, so BLOCKED. Windows commit charge was about 46 of 58 GB, mostly other applications.
  - During attempt 2, the D: volume holding the conda environment disappeared, so the pinned environment was rebuilt on F: (`model_cache/envs/thesis-model-v2`, torch 2.11.0+cu128, xgboost 3.2.0).
  - Attempt 3 changed only the memory order (at most two XGBoost matrices alive) and recomputed every fold in the new environment, about 15-19 min per fold.
- **Phase 3:** CPU, about 15 s.
- **Tests:** all 20 tests pass (15 from v1 + 5 new); results in `results/v2_tests/pytest_junit.xml` and in figure `v2_fig0_tests_and_checks.png`.

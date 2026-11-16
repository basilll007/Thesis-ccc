# Model report v1: neighbourhood models of macrophage and tumour cell state (baseline stage)

Branch `model-v1`. Data: Mendeley PDAC-SMI CosMx (development; 16 sections, 13 tumours, author labels) and GSE310352 CosMx (external test; 8 slides, marker labels). All numbers are read from files in `model_dev/results/` (index: `results/phase7/key_numbers.csv`; tables: `results/phase7/results_tables.md`). Plans were committed before each phase (`model_dev/plans/`); amendments and failures are in `reports/PROGRESS_LOG.md`. These data measure **expression**; every result is a **predictive association**, not a causal or clinical claim. The tumour (patient) is the unit for every confidence interval.

## Lead question
**Does a graph model predict macrophage and tumour cell state from neighbours better than simple baselines, and how much of that comes from OTHER cell types?**

1. **Neighbours help, but a graph model is not better than simple baselines.**
   - For the primary endpoint (myeloid C1Q_TAM program, grouped 5-fold), own covariates give held-out R2 0.300 [0.188, 0.383].
   - Aggregated neighbourhood features add **+0.034 [0.018, 0.048]** (XGBoost on per-type neighbour gene means) and +0.033 [0.023, 0.043] (ridge on neighbour programs).
   - The SAGE graph model reaches 0.336 [0.220, 0.432], which is **+0.0025 [-0.013, 0.015] versus the strongest baseline: not better**.
   - The graph model beats the baselines only on the myeloid 11-gene expression vector (+0.013 [0.005, 0.020]). That edge does not transfer to the external set.
   - For tumour cells the graph models are clearly worse than the baselines (classical -0.097 [-0.172, -0.053]).
2. **Other cell types carry a small but non-zero share for macrophages, and none for tumour cells.**
   - Within the graph model, the cross-type contribution (full minus same-type-neighbours-only) is +0.019 [0.004, 0.028] for C1Q_TAM, +0.015 [0.004, 0.028] for SPP1_TAM and +0.014 [0.009, 0.018] for the gene vector.
   - For tumour classical it is +0.007 [-0.004, 0.024], CI including 0.
   - The CAF neighbourhood is the largest other-type contributor to the macrophage prediction: permutation importance 0.011 [0.009, 0.012] R2, versus 0.033 [0.022, 0.041] for other macrophages.
3. **On the external set (single touch), the aggregated-neighbourhood gain transfers.**
   - XGBoost +0.072 [0.031, 0.107] over own covariates for C1Q_TAM.
   - The graph model's C1Q_TAM gain does not transfer: +0.014 [-0.014, 0.046].
   - The external units are 6 slides or slide groups whose independence and patient identity are UNVERIFIED.

## Results: development set (grouped 5-fold over 13 tumours; pooled held-out R2, 95 % tumour-bootstrap CI)
Per-fold and per-configuration values: `results/phase2/baseline_results.csv` (5,250 rows, every grid point) and `results/phase3/graph_results.csv`. Selected configurations are listed in `results/phase2/baseline_summary.csv` and `results/phase3/graph_summary.csv`. Gains are versus own covariates (paired tumour bootstrap).

| Receiver / target | own covariates | ridge, neighbour programs | XGBoost, neighbour genes | MLP, neighbour genes | graph SAGE | graph attention |
|---|---|---|---|---|---|---|
| myeloid C1Q_TAM (primary) | 0.300 [0.188, 0.383] | 0.333 [0.223, 0.417]; +0.033 [0.023, 0.043] | 0.334 [0.222, 0.424]; +0.034 [0.018, 0.048] | 0.187 [0.008, 0.316]; -0.113 [-0.225, -0.045] | 0.336 [0.220, 0.432] | 0.335 [0.215, 0.431] |
| myeloid SPP1_TAM | 0.335 [0.254, 0.398] | 0.383 [0.299, 0.443] | 0.385 [0.312, 0.436]; best is XGBoost on programs, 0.387 [0.307, 0.446] | 0.334 [0.234, 0.405] | 0.377 [0.301, 0.428] | 0.381 [0.307, 0.432] |
| myeloid 11-gene vector | 0.265 [0.191, 0.320] | 0.277 [0.201, 0.332] | n/a (XGBoost fitted to programs only) | 0.183 [0.087, 0.265] | 0.289 [0.214, 0.348] | 0.286 [0.210, 0.344] |
| tumour classical (positive control) | 0.618 [0.502, 0.711] | 0.693 [0.597, 0.766] | 0.694 [0.601, 0.768]; +0.076 [0.054, 0.107] | 0.547 [0.416, 0.644] | 0.597 [0.444, 0.703] | 0.601 [0.448, 0.710] |
| tumour basal_like | 0.397 [0.289, 0.512] | 0.438 [0.325, 0.552] | 0.444 [0.342, 0.550]; +0.046 [0.033, 0.061] | 0.265 [0.126, 0.406] | 0.281 [0.196, 0.367] | 0.296 [0.234, 0.367] |
| tumour 18-gene vector | 0.425 [0.352, 0.514] | 0.435 [0.362, 0.524] | n/a | 0.252 [0.177, 0.353] | 0.436 [0.364, 0.523] | 0.432 [0.360, 0.519] |

**Graph model versus the strongest baseline** (paired tumour bootstrap; seed-mean of 3 seeds; per-seed values in `graph_summary.csv`, e.g. SAGE C1Q_TAM 0.330 / 0.340 / 0.339):

| Target | strongest baseline | SAGE minus baseline | attention minus baseline |
|---|---|---|---|
| myeloid C1Q_TAM | XGBoost/genes 0.334 | +0.0025 [-0.013, 0.015] | +0.0010 [-0.016, 0.014] |
| myeloid SPP1_TAM | XGBoost/programs 0.387 | -0.010 [-0.021, 0.007] | -0.006 [-0.017, 0.008] |
| myeloid gene vector | ridge/programs 0.277 | **+0.013 [0.005, 0.020]** | **+0.0095 [0.002, 0.018]** |
| tumour classical | XGBoost/genes 0.694 | -0.097 [-0.172, -0.053] | -0.093 [-0.165, -0.051] |
| tumour basal_like | XGBoost/genes 0.444 | -0.163 [-0.330, -0.045] | -0.148 [-0.279, -0.051] |
| tumour gene vector | ridge/programs 0.435 | +0.0005 [-0.005, 0.007] | -0.004 [-0.011, 0.003] |

**Reproduction check.** The chain Phase 2 ridge analysis reproduces exactly from the new datasets: delta R2 0.0299 / 0.0561 / 0.1411 / 0.0568, |diff| 0.0000 against a tolerance of 0.01 (`results/phase2/reproduction_check.csv`).

**GPU vs CPU check.** GPU ridge matches the archived CPU ridge on 1,050 fold x target results to within 4e-5 R2.

## Ablations (SAGE, lr 3e-4, 3 seeds; same folds; `results/phase3/ablation_summary.csv`)
| Quantity (R2 full - R2 ablated) | myeloid C1Q_TAM | myeloid SPP1_TAM | myeloid gene vector | tumour classical |
|---|---|---|---|---|
| **cross-type contribution** (minus same-type neighbours only) | **+0.019 [0.004, 0.028]** | **+0.015 [0.004, 0.028]** | **+0.014 [0.009, 0.018]** | +0.007 [-0.004, 0.024] |
| minus other-type neighbours only | +0.028 [0.013, 0.035] | +0.047 [0.038, 0.055] | +0.021 [0.013, 0.026] | +0.048 [0.027, 0.069] |
| minus composition only (no neighbour expression) | +0.028 [0.004, 0.048] | +0.031 [0.014, 0.050] | +0.019 [0.010, 0.028] | not run (declared scope) |
| minus tile-shuffled neighbour expression (zonation control) | +0.039 [0.026, 0.049] | +0.035 [0.025, 0.048] | +0.021 [0.015, 0.027] | not run |
| minus no graph | +0.029 [0.005, 0.051] | +0.052 [0.037, 0.071] | +0.022 [0.013, 0.033] | +0.031 [0.011, 0.052] |

Readings:
- **Neighbour expression matters beyond tile-level zonation.** Replacing each neighbour's expression with that of a random same-class cell from the same 200 um tile costs the macrophage prediction 0.035-0.039. The model therefore uses finer-scale, cell-to-cell co-variation, not only regional averages.
- **Cross-type information is a minority share.** For macrophages it is smaller than the contribution of same-type neighbours (full minus other-type-only, 0.028-0.047), but it is non-zero. For tumour cells it is indistinguishable from zero: tumour state is predicted from nearby tumour cells.
- **No-graph check.** The no-graph model (0.307 for C1Q_TAM) is close to the own-covariates ridge (0.300), as intended.
- **Tumour basal_like.** All tumour basal_like differences have very wide CIs. The full graph model is unstable for this target (seed values 0.266-0.294), and ablating the graph does not hurt it (no-graph -0.059 [-0.215, 0.045]).

## External test (GSE310352, single touch per receiver; `results/phase4/external_results.csv`)
**Setup.**
- Frozen preprocessing; models retrained once on the whole development set with internally selected hyperparameters.
- Receivers by marker labels: 68,094 myeloid cells; 64,000 tumour cells (capped at 8,000 per slide).
- CIs resample 6 units (slides 5/5b2/5b3 grouped). These units are **not verified independent patients**, and patient overlap with the development set is UNVERIFIED.

| Target | model | internal R2 | external R2 [unit CI] | gain vs own (external) [CI] |
|---|---|---|---|---|
| myeloid C1Q_TAM | own covariates | 0.300 | 0.310 [0.211, 0.463] | - |
| | XGBoost/genes | 0.334 | 0.382 [0.304, 0.500] | **+0.072 [0.031, 0.107]** |
| | XGBoost/programs | 0.331 | 0.391 [0.308, 0.520] | **+0.082 [0.051, 0.111]** |
| | ridge/programs | 0.333 | 0.359 [0.271, 0.496] | **+0.049 [0.029, 0.068]** |
| | graph SAGE | 0.336 | 0.324 [0.239, 0.452] | +0.014 [-0.014, 0.046] |
| myeloid SPP1_TAM | graph SAGE | 0.377 | 0.375 [0.272, 0.466] | **+0.050 [0.016, 0.086]** |
| | XGBoost/genes | 0.385 | 0.375 [0.270, 0.467] | **+0.050 [0.020, 0.076]** |
| myeloid gene vector | ridge/programs | 0.277 | 0.274 [0.212, 0.341] | **+0.023 [0.014, 0.033]** |
| | graph SAGE | 0.289 | 0.266 [0.203, 0.339] | +0.015 [-0.002, 0.029] |
| tumour classical | own covariates | 0.618 | **-0.074 [-0.346, 0.266]** | - |
| | XGBoost/genes | 0.694 | 0.587 [0.493, 0.707] | +0.661 [0.449, 0.842] |
| | graph SAGE | 0.597 | 0.466 [0.351, 0.593] | +0.540 [0.325, 0.724] |
| tumour basal_like | XGBoost/genes | 0.444 | 0.500 [0.419, 0.546] | +0.125 [0.100, 0.144] |
| | graph SAGE | 0.281 | 0.449 [0.358, 0.504] | +0.073 [0.024, 0.106] |

Readings:
- **Macrophages.** The aggregated-neighbourhood models keep, and even raise, their internal performance externally (internal-to-external drop -0.048 for XGBoost/genes). The graph model drops slightly (0.013 for C1Q_TAM), and its gene-vector advantage disappears (0.266 versus 0.274 for ridge/programs).
- **Per slide.** Slide 5b2 is hardest for every model, with C1Q_TAM R2 0.06-0.18 (`fig3_external_test`).
- **Tumour classical, own covariates.** This model collapses externally (-0.074). Its strongest input is the same-type tile mean of the target. With marker labels the "tumour" class differs from the author-labelled one, so this covariate shifts. The huge external "gains" for tumour classical mainly reflect that baseline's failure, not stronger neighbourhood signal.

## Interpretation (Phase 5; gated; `results/phase5/`)
- **Gate.** Only the myeloid gene vector passed (SAGE beats the strongest baseline, +0.013 [0.005, 0.020]). Program targets did not pass, so no interpretation is reported for them. Attention weights are not used.
- **Tumour-level permutation importance** (R2 lost; 13 tumours):
  - permuting **other myeloid neighbours** costs 0.033 [0.022, 0.041];
  - **CAF neighbours** 0.011 [0.009, 0.012];
  - tumour neighbours 0.004 [0.001, 0.008];
  - T-cell neighbours 0.0008 [0.0003, 0.0013].
  - Zeroing the declared CAF, tumour or T-cell program genes in neighbours costs <= 0.0005, with CIs including 0. So the CAF signal is carried by the broad neighbour profile, not by the curated program genes.
- **In-silico perturbation.** This ran because the cross-type CI excludes 0 for this target. **Every row is a model-based hypothesis.** Setting the neighbouring CAFs' CXCL12 readout to its 90th percentile changes predicted macrophage expression by:
  - SELENOP +0.074 [0.058, 0.094];
  - CD163 +0.044 [0.031, 0.060];
  - C1QA +0.024 [0.015, 0.035];
  - SPP1 -0.041 [-0.052, -0.033].
  - These are log1p CP10k units: small shifts in predicted expression, not cell-state switches.
  - Setting neighbour CAF myCAF genes high raises predicted GPNMB (+0.068 [0.052, 0.086]) and C1QC (+0.063 [0.051, 0.075]).
  - Setting basal-like genes of neighbouring tumour cells high raises predicted SPP1 (+0.048 [0.019, 0.092]).
  - The CXCL12-CAF -> C1Q/CD163/SELENOP macrophage direction agrees with the exploratory CAF <-> C1Q_TAM adjacency lead of the chain study (`docs/CHAIN_SUMMARY_v1.md`, claim C10). Both are associations in the same kind of data, so the agreement is not independent confirmation.

## Pretraining feasibility (Phase 6; exploratory; `results/phase6/`)
- A masked-expression encoder (pretrained inside each fold on training cells only) with a linear probe does **not** help like for like. The neighbourhood-only probe gets C1Q_TAM 0.192 [0.115, 0.237], which is -0.144 [-0.207, -0.080] versus the supervised graph model, and tumour classical 0.311 [0.217, 0.390].
- A probe that also sees the cell's own non-target genes scores higher (C1Q_TAM 0.485 [0.382, 0.544]). It is not comparable: the cell's other genes are correlated with its targets.
- **Learning curve** (1 seed, nested tumour subsets): R2 rises by about 0.064 per tenfold more training cells, for both macrophages and tumour cells. The rough log-linear extrapolation to 10x the current data gives about 0.36 for C1Q_TAM and about 0.59 for tumour classical. A larger pretraining effort would need on the order of 10x more tumours (around 100+ patients) to gain about 0.06 R2. This is a rough extrapolation, not a foundation-model claim.

## Leakage and quality tests (`model_dev/tests/`, 15 tests pass)
- **Target leakage.** A planted copy of the target in the inputs is detected even under an innocent provenance tag (|r| > 0.999 probe). Self-sourced target genes and programs are detected by provenance. Neighbour copies of target genes are allowed.
- **Leave-one-out tile means** exclude the cell itself.
- **Scalers** are fitted on training rows only.
- **2-hop neighbourhoods** never contain the seed cell (the chain-study SGC bug class). In real data, 1,861 return paths per 256 seeds were masked.
- **Folds:** a tumour never straddles folds, and each fold has >= 2 tumours. Inner validation is disjoint from test.
- **GPU estimators:** GPU ridge equals sklearn Ridge, and the GPU scaler equals the CPU scaler with NaN handling.
- **Metric sign:** the sign of paired R2 differences is pinned by a test. A sign-inversion bug in all callers was found and fixed **before any result was reported**.

## Compute
- Hardware: RTX 5060 Laptop GPU (8 GB), 31 GB RAM. CUDA 12.8 torch 2.11.
- **Phase 2:** GPU, 10 folds x all models in about 41 min of fold time; peak VRAM 7.1 GB.
- **Phase 3:** 120 graph runs in 1.73 GPU-h, plus 120 ablation runs in 1.87 GPU-h. Peak VRAM 2.0 GB; median early-stopped epochs 17 and 23.5. Checkpoints for every run are in `model_cache/checkpoints/` (manifest with SHA-256: `results/phase7/checkpoint_manifest.csv`, 245 files).
- **Phase 6:** about 90 s of pretraining per fold.
- **Interruptions.** The host stopped CPU runs three times for low system memory, because other applications held about 40 of the 50 GB commit limit. From then on, all training ran on the GPU in resumable foreground chunks with per-run checkpoints.

## Limitations
1. **Small, dependent sample.** 13 development tumours give wide CIs. The external units are not verified independent patients.
2. **Annotation shift.** Models were trained on author labels and tested on marker labels. The tumour own-covariate baseline breaks under that shift.
3. **Cross-type contribution is model-internal.** It is measured within a graph model that does not beat the aggregated baselines on program targets. The same-type-only ablation also removes the *positions* of other-type cells, so the contribution mixes neighbour expression with neighbour identity and layout.
4. **Segmentation spillover is not modelled.** It can make a macrophage look like its CAF neighbours. Neighbour-to-receiver co-variation is therefore not evidence of signalling.
5. **Measurement limits.** Programs are short CosMx gene lists (5-10 genes) above a type-aware background. T cells were not testable as receivers.
6. **Recorded deviations** (all logged before the affected results were seen):
   - fold rule and tumour cap (Phase 1 amendments);
   - XGBoost-CUDA in place of sklearn HGB and GPU execution (Phase 2 amendment);
   - XGBoost and the MLP were not refitted on all training folds per the plan's scope, and XGBoost was fitted to program targets only.

## Model card and data card
`reports/MODEL_CARD.md` and `reports/DATA_CARD.md`.

## Next steps (at most 5)
1. Add a segmentation-spillover control (neighbour-transcript contamination model, or a protein-based identity check) before reading any cross-type co-variation biologically.
2. Repeat the cross-type ablation with aggregated baselines: XGBoost with versus without other-type neighbour blocks, so the cross-type share is measured on the best model rather than on the graph model.
3. Obtain patient-level metadata for GSE310352, or a third cohort, so the external test has verified independent patients and more than 6 units.
4. Test the CAF -> C1Q/SELENOP macrophage hypothesis experimentally, for example macrophage-CAF co-culture with CXCL12 blocking and C1Q/SELENOP/SPP1 readout. It is a hypothesis from a model, not a finding.
5. Scale data before scaling models. The learning curve suggests about 0.06 R2 per tenfold more cells, so more tumours, not bigger networks, are the bottleneck.

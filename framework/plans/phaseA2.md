# Phase A2 Plan: Unit Tests and Thesis Reproduction Suite

## 1. Objectives
Validate the `ccc_audit` package through two distinct testing suites:
1. **Unit Tests (Synthetic Only):** Purely algorithmic and statistical tests on generated synthetic data (no real patient data inside the test suite).
2. **Reproduction Tests (Thesis Benchmarks):** Run against the real dataset outputs from `F:\Thesis` specified via a local configuration to confirm numerical and methodological fidelity to prior workstreams.

## 2. Unit Test Suite Design (`framework/tests/test_unit.py`)

### Test 1: Null Model Calibration Under the Null (Type I Error)
- **Protocol:** Generate synthetic tissue with random coordinates and independent cell labels with identical expression distributions ($H_0$).
- **Nulls tested:** Global label permutation and block label permutation.
- **Tolerance:** Under 100 simulations with $\alpha = 0.05$, the empirical rejection rate must lie near 0.05 (binomial 95% acceptance range: $[0.01, 0.11]$).

### Test 2: Planted Effect Sensitivity (Power)
- **Protocol:** Plant a strong contact effect ($f = 2.0$) on receiver cells in contact with active senders.
- **Expectation:** Empirical rejection rate (power) for proximity ratio $Q$ under global and block nulls $\ge 0.80$.

### Test 3: Patient Leakage Guard & Grouped CV
- **Protocol:** Synthetic dataset with 5 patients. Run `evaluate_neighbour_prediction`.
- **Validation:**
  - Verify that train and test folds have strictly disjoint `patient_id` sets for every split.
  - Leakage guard: When a patient-level target is evaluated, ensure no patient appears across fold boundaries.

### Test 4: Spillover Detection & Distance Binning
- **Case A (Contaminated):** Introduce synthetic transcript mixing where touching receivers ($\le 12\,\mu\text{m}$) receive sender marker transcripts.
  - **Verdict:** `is_explained_by_spillover == True`, touching excess $\ge 10\%$.
- **Case B (Clean):** Receivers and senders express distinct marker sets without mixing.
  - **Verdict:** `is_explained_by_spillover == False`, touching excess $< 5\%$.

### Test 5: Exhaustive Verdict Decision Branch Coverage
- Synthetic test verifying that each of the 5 branches of `AuditVerdict` is triggered by its declared criteria:
  - `NOT_TESTABLE`: cell count < 25 or edges < 10.
  - `EXPLAINED_BY_ABUNDANCE`: $C$ or $M$ significant, $Q$ non-significant.
  - `EXPLAINED_BY_ZONATION`: Global $Q$ significant, Block $Q$ non-significant, or zonation FPR $> 10\%$.
  - `EXPLAINED_BY_SPILLOVER`: Excess contamination $\ge 10\%$ and signal confined to touching distance.
  - `SUPPORTED`: Significant $Q$ under block null, low spillover, power $\ge 0.80$.

### Test 6: Input Contract Validation
- `validate_spatial_data` must raise `DataContractError` when required columns (`cell_type`, `patient_id`, `section_id`, `x_um`, `y_um`) are missing, coordinates are non-numeric or contain NaNs.

---

## 3. Reproduction Test Suite Design (`framework/tests/test_reproduction.py`)

### Pre-declared Benchmark Values and Tolerances
Reproduce benchmarks recorded across thesis workstreams:

1. **v3 Tumour-Level CAF->T Proximity Result (`results/pdac/v3_phase1/tumour_level_summary_by_null.csv`):**
   - Interaction: `consensus_T`, headline `CXCL12->CXCR4`, `CAF_fibroblast -> T_cell`, Score $Q$.
   - Evaluated across 13 tumours.
   - Benchmark median tumour mean $z_Q$:
     - `global`: $0.433$ (Tolerance: $\pm 0.02$)
     - `block200`: $0.863$ (Tolerance: $\pm 0.02$)
     - `block100`: $0.811$ (Tolerance: $\pm 0.02$)
     - `block50`: $1.298$ (Tolerance: $\pm 0.02$)
     - `dist`: $0.579$ (Tolerance: $\pm 0.02$)
   - Number of tumours with $z_Q > 0$: 10/13 (global), 11/13 (block200), 12/13 (block100), 11/13 (block50), 10/13 (dist).

2. **model-v2 Spillover Reference $Q$ for C1Q_TAM (`model_dev/results/v2_phase1/phase1_summary.csv`):**
   - Target: `myeloid:C1Q_TAM` score under XGBoost.
   - Benchmark CAF-expression contribution $Q = +0.004$ (exact summary value: $0.00392$, tolerance $\pm 0.0005$).

3. **model-v2 XGBoost Cross-Type Share (`model_dev/results/v2_phase2/phase2_quantities.csv`):**
   - Target: `myeloid:C1Q_TAM` cross-type contribution ($R^2_{\text{all}} - R^2_{\text{same\_type}}$).
   - Benchmark cross-type contribution = $-0.005$ (exact summary value: $-0.00495$, tolerance $\pm 0.0005$).

4. **model-v1 Ridge Neighbour Gains (`model_dev/results/v2_phase2/reproduction_check_vs_v1.csv`):**
   - Target: `myeloid,ridge (programs set),myeloid:C1Q_TAM`.
   - Benchmark held-out $R^2$: $0.332645$ (tolerance $\pm 0.0005$).
   - Own covariates baseline $R^2$: $0.3000 \pm 0.001$.
   - Neighbour gain = $+0.0326$ ($+0.033$).

### Local Config File
Stored at `framework/tests/local_test_config.yaml` (git-ignored to ensure machine portability).

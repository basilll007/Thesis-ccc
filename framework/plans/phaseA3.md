# Phase A3 Plan: End-to-End Demonstration on Mendeley PDAC-SMI Development Set

## 1. Objectives
Execute the full `ccc-audit` command-line interface on the Mendeley PDAC-SMI development set (representative section U5-a) for headline cell-cell interactions and receiver programs to verify end-to-end operation, generate all deliverables (tables, figures, markdown report), compare verdicts with thesis conclusions, and record runtime and peak memory.

## 2. Dataset & Configuration
- **Dataset:** Mendeley PDAC-SMI (section U5-a, standard contract formatted as `data/processed/mendeley_u5a_contract.h5ad`).
- **Contract columns:**
  - `cell_type`: `annotation_majortypes`
  - `patient_id`: `donor`
  - `section_id`: `sample_id`
  - `x_um`, `y_um`: physical coordinates in micrometres (converted from pixels via $0.18\,\mu\text{m/px}$).
- **Graph:** 6-NN Euclidean contact graph on physical coordinates.
- **Null models:**
  - Global label permutation: 500 permutations.
  - Spatially constrained block permutation: $100\,\mu\text{m}$ square tiles with random offset, 500 permutations.
  - Random seed: 42.

## 3. Headline Interactions Audited
1. **`CAF -> CD8+ T` (`CXCL12` -> `CXCR4`):**
   - Headline thesis interaction.
   - Thesis conclusion: In U5-a, $Q$ does not shrink under block nulls but simulation calibration shows it is confounded by regional zonation (42–88% FPR) and underpowered in small sections.
   - Expected verdict: `EXPLAINED_BY_ZONATION` or `NOT_TESTABLE`.
2. **`CAF -> Macrophage` (`COL1A1` -> `CD44` / Matrix-signaling):**
   - High-confidence contact candidate in initial analyses.
   - Thesis conclusion (model-v2): Explained by segmentation spillover (23% marker excess in touching macrophages, signal drops to zero in separated cells).
   - Expected verdict: `EXPLAINED_BY_SPILLOVER`.
3. **`CAF -> Macrophage` (`FN1` -> `CD44`):**
   - Thesis conclusion (model-v2): Explained by segmentation spillover / abundance.
   - Expected verdict: `EXPLAINED_BY_SPILLOVER` or `EXPLAINED_BY_ABUNDANCE`.
4. **`Macrophage -> Malignant` (`SPP1` -> `CD44`):**
   - High abundance marker interaction.
   - Thesis conclusion: Saturated by cell abundance, $Q \le 1.0$.
   - Expected verdict: `EXPLAINED_BY_ABUNDANCE`.

## 4. Resource & Telemetry Constraints
- CPU threads capped at 4.
- Memory telemetry: Measure initial RAM, peak RSS memory via `psutil`, and total elapsed wall time.
- Stop rule: If a phase fails twice with the same error, log, mark BLOCKED, continue.

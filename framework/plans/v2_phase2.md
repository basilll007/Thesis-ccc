# v2 Phase 2 plan: fix the framework

Committed before any Phase 2 result is computed. The package is still called `ccc_audit` here; it is
renamed in Phase 3.

## 2.1 Reproduction tests recompute from inputs

The v1 tests (`framework/tests/test_reproduction.py`) only re-read saved result CSVs, so they could not fail.
New rule: every reproduction test loads original inputs and recomputes with the package.

Original inputs:
- raw Mendeley h5ad `F:/Thesis/data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad`;
- the thesis cell annotations `results/pdac/phaseB/sections/<section>_consensus_labels.csv.gz` (labels are an
  input to the audit, not a result of it).

Loader: new `ccc_audit.datasets.load_mendeley_section`. It is a port of the thesis loader and QC:
- Base probes only;
- 0.18 um/px;
- cells with >= 10 counts, genes in >= 3 cells;
- log1p(CP10K).

| ID | Thesis value (file) | Recomputation with the package | Tolerance (declared now) |
|---|---|---|---|
| R1a-e | v3 tumour-level CAF_fibroblast->T_cell CXCL12-CXCR4, annotation consensus_T, score Q: `median_tumour_mean_z` and `n_mean_z_pos` for nulls global, block200, block100, block50, dist (`results/pdac/v3_phase1/tumour_level_summary_by_null.csv`) | All 16 sections. Package 6-NN graph and `compute_ccc_scores`, 500 perms, seed 42. Nulls: global (`global_null`), block (`block_null`, `random_offset=False`, same tiling as thesis), dist (new port of the thesis distance-band null, `nulls.distance_band_null`). Tumour = mean of section z_Q over the tumour's sections. Median over 13 tumours. | median z +- 0.02; n_mean_z_pos exact; n_tumours = 13 |
| R2 | per-section observed `ratio_obs`, `matched_obs`, `edge_obs`, `n_qualifying_edges` of the same headline under the global null (`results/pdac/v3_phase1/headline_by_null.csv`, annotation consensus_T) | same run as R1 | relative 1e-6; edges exact |
| R3 | model-v2 spillover reference Q (C1Q_TAM, +0.00398) | **Cannot be recomputed by the package.** The package does not implement the model_dev XGBoost feature pipeline. | removed from the reproduction suite; reported as NOT REPRODUCED |
| R4 | model-v2 XGBoost cross-type share (-0.00515) | same reason as R3 | NOT REPRODUCED |
| R5 | model-v1 ridge R2 (0.3326) | The package ridge (`neighbour_models`) uses different features from model_dev | NOT REPRODUCED |

If R1/R2 fail, the test fails and the failure is reported. The tolerance is not changed after seeing the result.

## 2.2 Claim audit with every check reported (`ccc_audit.contact`, `ccc_audit.verdicts.audit_claim`)

### Claim

"Sender type S signals to receiver type R through ligand L and receptor Rg." Operationally: Rg expression in R cells is
higher when the cell touches an L-expressing S cell.

### Contact graph

- Symmetric 6-NN graph (union of directions, duplicates removed).
- S->R edges are the undirected pairs with one end in S and the other in R.
- d = median nearest-neighbour distance (cell-diameter proxy).
- Expression: `log` = log1p(CP10K); `lin` = CP10K. Counts come from `layers['counts']` if present, otherwise
  expm1(X).

### Statistics and nulls

All nulls are permutation nulls. Default 1000 permutations; add-one p = (k+1)/(n+1); minimum attainable
p = 1/1001 = 0.000999.

- **Contact test, global** (`p_contact_global`):
  - Q = E/M with E = mean over S->R edges of L_log[i]*Rg_log[j] and M = mean(L_log over S) * mean(Rg_log over R).
  - Null: permute L among S cells and Rg among R cells (cell types and positions fixed). This keeps type-level
    expression (abundance) and breaks only who-touches-whom.
- **Contact test, block** (`p_contact_block`): same, but permutation is within 100 um square tiles with a random
  offset.
- **Ring-control test** (`p_ring`):
  - Exposed receivers X: R cells with >= 1 contact neighbour that is an S cell with L count > 0.
  - Controls C: R cells not in X whose nearest L+ S cell lies within 4d.
  - Statistic: mean Rg_log(X) - mean Rg_log(C).
  - Null: permute Rg_log among X ∪ C within 100 um tiles (random offset). Not run if |X| < 5 or |C| < 5.
- **Naive abundance score** (`p_M_label`, thesis style): M under global cell-type label permutation. Used only to
  separate EXPLAINED_BY_ABUNDANCE from NOT_SUPPORTED.
- **Power** (plant-and-shift):
  - Multiply lin-scale Rg of exposed receivers by f = 1.5 (deterministic). Recompute each of the three statistics.
    shift = planted - observed.
  - Power of each test = fraction of its null draws with null + shift >= the null's (1 - alpha) quantile.
  - Chain power = min of the three (an optimistic bound, stated as such).
- **Spillover** (claim-specific: computed for this S->R pair on the full gene panel):
  - T = R cells with >= 1 S contact neighbour. I = R cells with no S contact neighbour and nearest S > 2d (if < 10 such
    cells: R cells with no S contact neighbour; if still < 10: NOT_RUN).
  - Markers: genes other than L and Rg with lin mean(S) / mean(I) >= 10 and mean(S) >= 1 CP10K; top 10 by ratio; fewer
    than 3 -> NOT_RUN.
  - phi = median over markers of (mean_T - mean_I)/(mean_S - mean_I), lin scale.
  - Contamination detected if phi >= 0.005 and the one-sided permutation p of phi (T/I labels permuted, 1000 draws)
    is < 0.05.
  - dR_obs = mean_T(Rg_lin) - mean_I(Rg_lin); dR_exp = phi * (mean_S(Rg_lin) - mean_I(Rg_lin)).

### Checks

Each check is PASS / FAIL / NOT_RUN, and every check is computed and reported.

| Check | PASS | FAIL | NOT_RUN |
|---|---|---|---|
| power | n_R >= 25, S->R edges >= 10, chain power >= 0.5 | otherwise | – |
| abundance | p_contact_global < 0.05 | otherwise | no S->R edges |
| zonation | p_contact_block < 0.05 **and** p_ring < 0.05 | otherwise | ring not runnable |
| spillover | contamination not detected, or dR_obs <= 0, or dR_exp < 0.5 * dR_obs | contamination detected and dR_obs > 0 and dR_exp >= 0.5 * dR_obs | markers or isolated receivers unavailable |

### Primary verdict, pre-declared priority rule (first match wins)

1. power FAIL -> NOT_TESTABLE
2. spillover FAIL -> EXPLAINED_BY_SPILLOVER
3. abundance FAIL -> EXPLAINED_BY_ABUNDANCE if p_M_label < 0.05, else NOT_SUPPORTED
4. zonation FAIL -> EXPLAINED_BY_ZONATION
5. any of zonation / spillover NOT_RUN -> NOT_TESTABLE (the claim cannot be cleared)
6. otherwise -> SUPPORTED

Rationale for the order:
- Testability comes first, because a non-result means nothing without power.
- Spillover comes before abundance and zonation because it is a data-quality defect that can create any of the
  other signals.
- Abundance comes before zonation because, if there is no contact enrichment at all, zonation has nothing to
  explain.

New verdict value NOT_SUPPORTED: no contact enrichment and no naive signal either. v1 mapped this case to
EXPLAINED_BY_ABUNDANCE. The v1 `decide_verdict` stays for backward compatibility, marked deprecated.

## 2.3 Permutations

- Default `n_perms = 1000` everywhere: CLI, calibration, and the `global_null` / `block_null` defaults. Minimum
  attainable p = 1/1001.
- The full-section U5-a demo is re-run at 1000 permutations.

## 2.4 Null-calibration unit test

Synthetic data:
- 200 datasets, 600 cells, uniform positions in 600 x 600 um.
- 3 types, clustered by a smooth spatial field so that type composition is spatially structured.
- L and Rg drawn as independent NB counts with type-specific means. No contact effect.

Run the contact tests with 1000 permutations. Seeds 0..199.

Test rule: the number of rejections at alpha = 0.05 must lie inside the central binomial 95 % interval for
Binomial(200, 0.05), computed in the test as `scipy.stats.binom.ppf([0.025, 0.975], 200, 0.05)`. This is done
separately for:
- global contact test;
- block contact test;
- ring test.

Diagnostic only, not a test: the same datasets under the v1 label-permutation Q test, to show whether it is
calibrated.

## 2.5 Demos

Re-run U5-a crop (8,000 cells) and full section with the v2 engine and 1000 permutations.
- Config claims are unchanged: CAF->CD8+ T CXCL12-CXCR4, CAF->Macrophage COL1A1-CD44, CAF->Macrophage FN1-CD44,
  Macrophage->Malignant SPP1-CD44. Labels come from the contract h5ad.
- Output: the multi-label table per claim and the primary verdict.
- Explain the v1 crop-vs-section disagreement for CAF->Macrophage. The hypothesis to check against the saved v1
  outputs, written down before re-running: the v1 CLI ran the spillover audit only for the **first** interaction's
  pair (CAF->CD8+ T: 11 receivers in the crop), used only the 6 LR genes as the marker universe, and applied the
  result to every claim. The verdict tree also stopped at the first failing check.

## 2.6 Tests

All tests must pass. Report counts (unit, calibration, reproduction). Pytest markers:
- `repro`: needs the local data;
- `slow`: calibration.

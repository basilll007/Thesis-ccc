# v2 Phase 2 summary: framework fixes

Plan: `framework/plans/v2_phase2.md` + `v2_phase2_amendment1.md` (both committed before these results).
Status: **DONE**.

## 1. Reproduction tests now recompute from inputs

`framework/tests/test_reproduction.py` loads the raw Mendeley h5ad and the thesis `consensus_T` labels, then
recomputes with the package: 16 sections, 5 nulls, 500 permutations, seed 42.

- **Run 1: all 6 tests FAILED** (`pytest_reproduction_run1_knn_bug.txt`, `reproduction_recomputed_sections_run1_knn_bug.csv`).
  - The matched baseline M agreed exactly in every section. The qualifying-edge counts did not (T1-b: 574 vs thesis 542).
  - Running the thesis loader and the thesis `knn_edges` in the thesis environment gave 542 again, with identical
    coordinates and labels.
- **Cause: a bug in the v1 package's `graphs.knn_graph`.**
  - It called `NearestNeighbors.kneighbors()` without passing the query points. scikit-learn then already excludes
    each point from its own neighbour list.
  - The function also dropped column 0, so it **discarded every cell's true nearest neighbour** and used neighbours
    2 to 7.
  - Fixed, with a regression test (`test_knn_graph_keeps_nearest_neighbour`).
  - Everything in v1 that used `knn_graph` was affected: the v1 demo C/E/M/Q tables and v1 verdicts. The new v2 claim
    audit (`contact.SpatialContext`) and the simulators build their own correct graph.
- **Run 2, after the fix: 6 / 6 PASS** (`pytest_reproduction_after_knn_fix.txt`, `reproduction_recomputed_sections.csv`):
  - R1: median tumour mean z_Q within ±0.02 and the number of tumours with z > 0 exact, for global, block200, block100,
    block50 and dist;
  - R2: per-section Q, M, E within relative 1e-6, and edge counts exact.
- **Not reproduced (removed from the suite):** model-v2 spillover reference Q (+0.004), XGBoost cross-type share
  (−0.005) and the model-v1 ridge R² (0.333).
  - The package does not implement the model_dev feature pipeline. The v1 "passes" for these three were reloads of
    saved CSVs, not reproductions.

## 2. Every check is reported, with one primary verdict by a declared rule

`verdicts.audit_claim` returns PASS / FAIL / NOT_RUN for power, abundance, zonation and spillover, plus one primary
verdict chosen by `PRIORITY_RULE`:

1. power FAIL → NOT_TESTABLE
2. spillover FAIL → EXPLAINED_BY_SPILLOVER
3. abundance FAIL → EXPLAINED_BY_ABUNDANCE (or NOT_SUPPORTED if the naive score is also null)
4. zonation FAIL → EXPLAINED_BY_ZONATION
5. a required control NOT_RUN → NOT_TESTABLE
6. otherwise SUPPORTED

The CLI writes a multi-label `claims_register.csv` and `claim_metrics.csv`. Spillover is computed separately for
each claim's own sender/receiver pair on all genes.

### U5-a demo (v2 engine, 1000 permutations)

| Claim | Crop (8,000 cells) checks [power / abundance / zonation / spillover] → verdict | Full section (74,793 cells) |
|---|---|---|
| CAF→CD8+ T CXCL12-CXCR4 | FAIL / FAIL / NOT_RUN / NOT_RUN → NOT_TESTABLE | FAIL / PASS / FAIL / PASS → NOT_TESTABLE |
| CAF→Macrophage COL1A1-CD44 | FAIL / FAIL / FAIL / PASS → NOT_TESTABLE | PASS / FAIL / FAIL / PASS → EXPLAINED_BY_ABUNDANCE |
| CAF→Macrophage FN1-CD44 | FAIL / FAIL / FAIL / PASS → NOT_TESTABLE | PASS / FAIL / FAIL / PASS → EXPLAINED_BY_ABUNDANCE |
| Macrophage→Malignant SPP1-CD44 | PASS / PASS / FAIL / PASS → EXPLAINED_BY_ZONATION | PASS / PASS / FAIL / PASS → EXPLAINED_BY_ZONATION |

### Why v1 disagreed between crop and section for CAF→Macrophage (v1: SPILLOVER in the crop, ABUNDANCE in the full section)

The saved v1 outputs (`framework/results/demo_u5a*/spillover_summary.csv`) show three defects:

1. **Wrong cell pair.** The v1 CLI ran the spillover audit once, for the **first** interaction's pair (CAF→CD8+ T),
   and applied the result to every claim. In the crop that pair has 11 CD8+ T cells.
2. **Wrong marker universe.** Markers were chosen only among the 6 ligand/receptor genes. In the crop, 2 "markers"
   gave a 99.8 % "excess". In the full section, 0 markers were found, so the excess was 0 %.
3. **Only the first failing check was reported.** The v1 tree checked spillover before abundance and stopped there.

Under v2, CAF→Macrophage gets the same direction in both: no contact enrichment.
- Q is 0.91 in the crop and 0.93 in the full section; the global contact p is 0.94 and 1.0.
- It is NOT_TESTABLE in the crop (chain power 0.18 < 0.5) and EXPLAINED_BY_ABUNDANCE in the full section (power 0.91).
- Spillover is PASS in both, although CAF-marker contamination of touching macrophages is real (phi = 0.15 crop,
  0.13 full). The reason: CD44 is lower in CAFs than in macrophages, so spillover cannot create a CD44 excess, and
  touching macrophages have *less* CD44 than isolated ones.

v1 also reported SUPPORTED for FN1-CD44 and SPP1-CD44 in the full section (`demo_u5a_full/claims_register.csv`).
FRAMEWORK_REPORT_v1 did not mention this. Under v2 these are EXPLAINED_BY_ABUNDANCE and EXPLAINED_BY_ZONATION.

## 3. Permutations

- Default is 1000 everywhere. The minimum attainable p is 1/1001 = 0.000999.
- Full-section demo at 1000 permutations: 46.3 s, peak RSS 1,664 MB (`demo_u5a_full_v2/runtime.txt`).
- Crop: 5.2 s, 397 MB.

## 4. Null calibration (200 datasets x 1000 permutations; `null_calibration_summary.csv`)

The binomial 95 % interval for 200 tests at alpha 0.05 is [4, 16] rejections.

| Test | Rejections / 200 | Inside interval |
|---|---|---|
| v2 contact, global within-type permutation | 8 | yes |
| v2 within-tile receiver contrast | 8 | yes |
| v2 ring control | 10 | yes |
| v2 block-null Q (info only) | 9 | yes |
| v1 label-permutation Q, global (diagnostic) | 0 | **no (conservative)** |
| v1 label-permutation Q, block 100 (diagnostic) | 11 | yes |

## 5. Test counts

29 / 29 pass:
- 6 v1 unit tests;
- 13 v2 audit tests (8 priority-rule branches, planted contact, null, claim-specific spillover, minimum p, kNN
  regression);
- 4 calibration tests;
- 6 reproduction tests (after the kNN fix).

# v2 Phase 4b summary: post-hoc engine v0.2.1 on fresh replicates

Plan: `framework/plans/v2_phase4b.md` (commit 2720a9e).
- The plan was written **after** the frozen v0.2.0 results on T2-b/T6-b were seen and before the U5-a results
  were seen.
- Fresh seeds, 50 replicates per setting x 3 sections x 14 settings x 2 claims = 4,200 audits, 0 errors.
- CPU, 4 workers.
- Tests passed before the run (`pytest_v021_before_run.txt`: 29 passed).

This is a **development-stage** evaluation, not an independent validation: v0.2.1 was designed using the outcomes of
the same scenario design.

## Pass criterion (as pre-declared in v2_phase4.md; n = 150 per setting per claim)

| | False-SUPPORTED <= 0.10 in every artefact setting | S1 2x sensitivity >= 0.80 | Result |
|---|---|---|---|
| **C1** CAF_fibroblast -> myeloid FN1-ITGB1 | met (max 0.040, S4 0.1) | **0.913** (T2-b 0.84, T6-b 0.90, U5-a 1.00) | **PASS** |
| **C2** myeloid -> tumour SPP1-CD44 | met (max 0.020) | **0.673** (T2-b 1.00, T6-b 0.02, U5-a 1.00) | **FAIL** |
| Both claims pooled (n = 300) | met (max 0.020 pooled) | 0.793 | **FAIL** (just below 0.80) |

Sensitivity above the power floor, i.e. the user's wording P(SUPPORTED | power PASS) at 2x:

| Claim | P(SUPPORTED \| power PASS) | Replicates testable |
|---|---|---|
| C1 | 0.913 | 100 % |
| C2 | 0.962 | 70 % |
| pooled | 0.933 | 85 % |

The C2 failure comes from T6-b, which has only 222 tumour receivers. There the engine correctly reports NOT_TESTABLE
in most replicates (S1 2x: 45 of 50 NOT_TESTABLE).

**No setting fools v0.2.1.** The highest false-SUPPORTED rate is 0.040 (C1, S4 spillover 0.1; Wilson upper bound
0.085).

## Confusion matrices (pooled over 3 sections, 150 per row)

C1 (CAF -> myeloid FN1-ITGB1):

| Setting | SUPP | ABUND | ZONE | SPILL | NOT_SUPP | NOT_TEST | Accuracy |
|---|---|---|---|---|---|---|---|
| S0 null | 2 | 25 | 45 | 4 | 74 | 0 | 0.99 |
| S1 1.25x | 49 | 17 | 40 | 0 | 44 | 0 | 0.33 |
| S1 1.5x | 83 | 6 | 47 | 0 | 14 | 0 | 0.55 |
| S1 2x | 137 | 0 | 13 | 0 | 0 | 0 | 0.91 |
| S2 abundance | 0 | 109 | 41 | 0 | 0 | 0 | 0.73 |
| S3 50 um | 0 | 25 | 102 | 0 | 23 | 0 | 0.68 |
| S3 200 um | 2 | 4 | 141 | 3 | 0 | 0 | 0.94 |
| S3 400 um | 1 | 4 | 145 | 0 | 0 | 0 | 0.97 |
| S4 0.1 | 6 | 1 | 32 | 65 | 46 | 0 | 0.43 |
| S4 0.2 | 3 | 1 | 17 | 94 | 35 | 0 | 0.63 |
| S4 0.3 | 3 | 3 | 15 | 109 | 20 | 0 | 0.73 |
| S5 n = 25 / 50 / 100 | 0 | 0 / 0 / 1 | 0 | 0 | 0 | 150 / 150 / 149 | 1.00 / 1.00 / 0.99 |

C2 (myeloid -> tumour SPP1-CD44):

| Setting | SUPP | ABUND | ZONE | SPILL | NOT_SUPP | NOT_TEST | Accuracy |
|---|---|---|---|---|---|---|---|
| S0 null | 0 | 81 | 55 | 0 | 0 | 14 | 1.00 |
| S1 1.25x | 99 | 18 | 5 | 0 | 0 | 28 | 0.66 |
| S1 1.5x | 100 | 7 | 8 | 0 | 0 | 35 | 0.67 |
| S1 2x | 101 | 1 | 3 | 0 | 0 | 45 | 0.67 |
| S2 abundance | 1 | 91 | 56 | 0 | 0 | 2 | 0.61 |
| S3 50 um | 1 | 84 | 63 | 0 | 0 | 2 | 0.42 |
| S3 200 um | 0 | 43 | 103 | 0 | 0 | 4 | 0.69 |
| S3 400 um | 0 | 30 | 115 | 0 | 0 | 5 | 0.77 |
| S4 0.1 / 0.2 / 0.3 | 0 / 3 / 2 | 75 / 68 / 71 | 56 / 61 / 62 | 1 / 3 / 2 | 0 | 18 / 15 / 13 | 0.01 / 0.02 / 0.01 |
| S5 n = 25 / 50 / 100 | 0 | 0 / 0 / 4 | 0 / 0 / 9 | 0 | 0 | 150 / 150 / 137 | 1.00 / 1.00 / 0.91 |

## Readings

1. **Not fooled.** False-SUPPORTED is <= 0.04 everywhere: 6 of 150 at C1 S4 0.1, and 2-3 of 150 in a few other
   settings.
2. **Sensitive where contact is identifiable.** Above the power floor, the 2x detection rate is 0.91-0.96. At 1.25x
   and 1.5x, sensitivity for C1 is lower (0.33, 0.55).
3. **Low-data honesty.** S5: 0 SUPPORTED and 99-100 % NOT_TESTABLE at 25, 50 and 100 receivers for C1. For C2 at
   100 receivers, 91 % are NOT_TESTABLE and the remaining 13 are ABUNDANCE/ZONATION. Those 13 are false refutations:
   9 % of that setting.
4. **Cause attribution is partly wrong.**
   - S2 (abundance) is labelled ZONATION in 27-37 % of replicates.
   - S0 is labelled ZONATION in 30-37 %.
   - Likely cause: the regional co-variation statistic responds to spatially structured library size in the real
     geometry (log1p of sparse counts depends on detection depth). This is a cause not planted in the simulation.
     **Not verified.**
   - C2 S4 is almost never called SPILLOVER. This was expected: myeloid CD44 < tumour CD44, so spillover cannot create
     the signal and the claim correctly ends as not supported.
5. **Power estimate is optimistic.** For C1, estimated chain power at 1.5x has a median of about 0.9, while the
   actual SUPPORTED rate at 1.5x is 0.55. The power floor therefore admits some underpowered tests. The S1 1.25x and
   1.5x misses show up as ZONATION / NOT_SUPPORTED rather than NOT_TESTABLE.

Files: `chunks/`, `benchmark_all_replicates.csv.gz`, `confusion_matrix_*.csv`, `metrics_by_setting.csv`,
`pass_criterion.json`.

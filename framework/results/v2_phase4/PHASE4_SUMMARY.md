# v2 Phase 4 summary: frozen engine v0.2.0 on the ground-truth benchmark

Plan: `framework/plans/v2_phase4.md` (commit 39b5cbc, written before any replicate was run).
Engine: cellalibi v0.2.0 (fc9bab8). 4,200 audits: 3 sections (T6-b, T2-b, U5-a) x 14 settings x 100 replicates.
Claim: CAF_fibroblast -> myeloid via FN1 -> ITGB1. All replicates completed (0 errors).
Files: `chunks/`, `benchmark_all_replicates.csv.gz`, `confusion_matrix_*.csv`, `metrics_by_setting.csv`,
`pass_criterion.json`.

## Verdict: **FAIL** (pre-declared criterion)

1. False-SUPPORTED rate <= 0.10 in every artefact setting: **met**. 0/300 in each of S0, S2, S3 x 3 and S4 x 3;
   Wilson upper bound 0.013.
2. S1 sensitivity at 2x >= 0.80: **not met**. 0/300 (0.00) in every section.

The engine is not fooled by any artefact, but only because it **almost never decides**:
- NOT_TESTABLE in 67-100 % of audits per setting;
- SUPPORTED in 0 of 4,200 audits;
- EXPLAINED_BY_SPILLOVER in 0 of 900 spillover replicates.

## Pooled confusion matrix (rows = true setting, n = 300 each)

| Setting | SUPPORTED | ABUNDANCE | ZONATION | SPILLOVER | NOT_SUPPORTED | NOT_TESTABLE |
|---|---|---|---|---|---|---|
| S0 null | 0 | 97 | 3 | 0 | 0 | 200 |
| S1 contact 1.25x | 0 | 1 | 12 | 0 | 0 | 287 |
| S1 contact 1.5x | 0 | 0 | 0 | 0 | 0 | 300 |
| S1 contact 2x | 0 | 0 | 0 | 0 | 0 | 300 |
| S2 abundance 2x | 0 | 95 | 5 | 0 | 0 | 200 |
| S3 zonation 50 um | 0 | 82 | 18 | 0 | 0 | 200 |
| S3 zonation 200 um | 0 | 76 | 24 | 0 | 0 | 200 |
| S3 zonation 400 um | 0 | 69 | 31 | 0 | 0 | 200 |
| S4 spillover 0.1 | 0 | 14 | 72 | 0 | 0 | 214 |
| S4 spillover 0.2 | 0 | 0 | 50 | 0 | 0 | 250 |
| S4 spillover 0.3 | 0 | 0 | 22 | 0 | 0 | 278 |
| S5 n = 25 | 0 | 0 | 0 | 0 | 0 | 300 |
| S5 n = 50 | 0 | 0 | 0 | 0 | 0 | 300 |
| S5 n = 100 | 0 | 0 | 0 | 0 | 0 | 300 |

Per-setting accuracy (exact expected verdict):

| Setting | Accuracy |
|---|---|
| S0 | 1.00 |
| S1 (all effect sizes) | 0.00 |
| S2 | 0.32 |
| S3 50 / 200 / 400 um | 0.06 / 0.08 / 0.10 |
| S4 (all levels) | 0.00 |
| S5 (NOT_TESTABLE is correct) | 1.00 |

## Why: three diagnosed failures (evidence in the chunk files)

1. **Biased power estimate.** The plant multiplied the *observed* receptor values by 1.5, so the many zero-count
   cells of a sparse gene stayed zero.
   - T2-b and T6-b: estimated chain power was about 0.07-0.08 in every replicate, so power FAIL in 100 % of them.
   - The actual joint detection rate of the three decision tests at 1.5x was 0.28 (T6-b) and 0.06 (T2-b). At 2x it
     was 0.68 and 0.36.
2. **Spillover check could not run.** The marker rule needs a sender / isolated-receiver ratio >= 10. The simulated
   receiver means are taken from real CosMx data, which already contain contamination, so no CAF marker reached
   10x.
   - Spillover was NOT_RUN in 100 % of T6-b and U5-a replicates ("fewer than 3 sender markers").
   - A required control that cannot run forces NOT_TESTABLE. In U5-a the power check passed (median 0.62-0.71), but
     every S1 replicate was still NOT_TESTABLE for this reason.
   - The same NOT_RUN status explains why no S4 replicate was called EXPLAINED_BY_SPILLOVER.
3. **Weak global Q test when exposure is saturated.** 70 % (U5-a), 87 % (T6-b) and 95 % (T2-b) of myeloid cells
   touch an FN1+ CAF. Q = E/M then barely moves.
   - T2-b, S1 2x: Q-test p < 0.05 in 42 % of replicates, versus 92 % for the within-tile receiver contrast.

## Which scenarios fool it

None. No artefact setting produced a single SUPPORTED. The weakness is the opposite one: **zero sensitivity**,
caused by the three issues above. Cause attribution was also poor:
- S3 was mostly labelled EXPLAINED_BY_ABUNDANCE, not ZONATION;
- S4 was mostly EXPLAINED_BY_ZONATION or NOT_TESTABLE, never SPILLOVER.

## Follow-up

These failures motivated the **post-hoc** engine v0.2.1. Its plan, `v2_phase4b.md`, was committed in 2720a9e after
the T2-b/T6-b results and before the U5-a results were inspected. v0.2.1 is benchmarked on fresh seeds with a second
claim. The v0.2.0 results above are kept unchanged as the pre-registered outcome.

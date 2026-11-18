# v2 Phase 4b plan: post-hoc engine revision v0.2.1 and a fresh-seed benchmark

**Status: post-hoc.** This revision was designed after seeing the frozen v0.2.0 benchmark (Phase 4) on T2-b and
T6-b: all 2,800 audits there were NOT_TESTABLE. It was committed before any U5-a v0.2.0 result was inspected.
The v0.2.0 results stay the official pre-registered Phase 4 result and are reported unchanged.

v0.2.1 is benchmarked on **fresh replicates** (new seeds) with the same scenario design, plus a second claim. Because
the design of v0.2.1 was informed by Phase 4, this is a development-stage test, **not an independent validation**.

## Diagnosed v0.2.0 failures and the v0.2.1 changes

| # | v0.2.0 failure (evidence from T2-b / T6-b chunks) | v0.2.1 change |
|---|---|---|
| F1 | The power estimate multiplied *observed* receptor values by 1.5, so zero-count cells stayed zero. ITGB1 is sparse, so estimated chain power was about 0.08 while the actual joint detection rate at 1.5x was 0.28 (T6-b). Result: power FAIL in 100 % of replicates, including S1 2x. | Plant on the **mean**: add Poisson((f-1) x lambda_j) receptor counts to exposed receivers. lambda_j = mean receiver receptor count x the cell's library-size ratio. Use 20 stochastic plants, same shift method. Without counts: additive (f-1) x mean. |
| F2 | The spillover marker rule (sender/isolated-receiver ratio >= 10) found fewer than 3 markers in every T6-b replicate. Real contamination already sits in the real-derived receiver means, so spillover was NOT_RUN and the claim NOT_TESTABLE. | Marker ratio threshold 10 -> **3** (other marker rules unchanged). |
| F3 | The global Q test was the weakest component: in T2-b at S1 2x, p < 0.05 in 42 % of replicates versus 92 % for the within-tile receiver contrast. Q is diluted when 87-95 % of receivers are exposed. | The abundance check uses the **global receiver contrast**: exposed vs unexposed receivers, Rg permuted among all receivers. It uses the same exchangeability as the other v0.2 tests. Q is still computed and reported. |
| F4 | With the receiver contrast, zonation that is independent of exposure no longer shows up as contact enrichment, so it would be mislabelled as abundance. | New **regional co-variation** statistic: the Pearson r across 100 um tiles between the sender mean L and the receiver mean Rg, with receiver tile means permuted (`p_zone`). In priority step 3, when the abundance check fails: p_zone < 0.05 -> EXPLAINED_BY_ZONATION, else p_M_label < 0.05 -> EXPLAINED_BY_ABUNDANCE, else NOT_SUPPORTED. |

Unchanged:
- power floor (chain power >= 0.5 at 1.5x, >= 25 receivers, >= 10 edges);
- chain power = min of the per-test powers (now rc_global, block contrast, ring);
- the zonation check (within-tile contrast AND ring);
- the spillover logic apart from the marker threshold;
- 1000 permutations, alpha 0.05;
- the rest of the priority rule.

An abundance NOT_RUN (< 5 exposed or < 5 unexposed receivers) now leads to NOT_TESTABLE through step 5.

Engine:
- `AuditConfig(engine="0.2.1")`, now the package default (version 0.2.1). `engine="0.2.0"` reproduces the frozen
  engine.
- Code commit: the commit containing this plan.
- Before running, the unit and calibration tests must pass. The calibration now also covers `p_rc_global` and
  `p_zone`. A failing test is fixed before the run, and the fix is logged with its commit.

## Benchmark design (fresh seeds)

Claims:
- **C1**: CAF_fibroblast -> myeloid, FN1 -> ITGB1 (the Phase 4 claim; high exposure, 87-95 %).
- **C2**: myeloid -> epithelial_ductal_tumor, SPP1 -> CD44 (thesis catalogue pair; most tumour cells do not touch a
  myeloid cell, so exposure is low). C2 tests the engine where contact is identifiable. Note: CD44 is lower in myeloid
  than in tumour cells in U5-a, so S4 spillover is not expected to create a false signal for C2.

Sections: the Phase 4 rule applied per claim (eligible = >= 100 sender and >= 100 receiver cells; smallest /
lower-median / largest).

Settings: the same 14 (S0; S1 1.25 / 1.5 / 2; S2; S3 50 / 200 / 400 um; S4 0.1 / 0.2 / 0.3; S5 25 / 50 / 100).
- **50 replicates** per setting, section and claim (declared reduction for runtime): 2 x 3 x 14 x 50 = 4,200 audits.
- Seeds: CRC32 of ("v021", claim, section, setting, level, replicate). There is no overlap with Phase 4 seeds.

Expected verdicts and pass criterion: as declared in `v2_phase4.md`. The criterion is evaluated per claim, pooled over
its 3 sections (n = 150 per setting), and also over both claims.

Additionally reported (the user's wording, "sensitivity above the power floor"):
- P(SUPPORTED | power PASS) for S1 at 2x;
- the fraction of S1 replicates that are testable.

## Secondary: v0.2.1 on the Phase 5 tool replicates

v0.2.1 is also run on the exported Phase 4 replicates used for the tool comparison: C1, 440 datasets, the same
counts the tools see. This gives a side-by-side with the tools.

These replicates come from the same generator as the data that informed v0.2.1. They are labelled
**development set**.

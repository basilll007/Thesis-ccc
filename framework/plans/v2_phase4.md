# v2 Phase 4 plan: ground-truth benchmark on real tissue geometry

Committed before any benchmark replicate is simulated or audited. Engine under test: `cellalibi` v2 claim audit
(the former `ccc_audit`, renamed in Phase 3), exactly as committed at the end of Phase 3. No engine parameter may be
changed after benchmark results are seen. Any later change is reported as a separate, post-hoc version, and the
original results are kept.

## Claim under test (fixed)

**CAF_fibroblast -> myeloid via FN1 -> ITGB1.**
- The pair comes from the thesis ligand-receptor catalogue (`src/pdac/lr_catalog.py`, PDAC_ADDITIONS).
- The sender and receiver are the thesis `consensus_T` labels.
- Reason for this pair: CAF->macrophage was the thesis spillover case, and ITGB1 is also expressed by CAFs (real
  mean CP10K about 2x myeloid in T6-b and U5-a, checked before this plan), so simulated spillover can create a false
  contact signal.

## Sections (fixed rule)

- Eligible: Mendeley sections whose `consensus_T` labels have >= 100 myeloid and >= 100 CAF_fibroblast cells after
  QC. That leaves room for the S5 subsample of 100 receivers and is well above the 25-receiver floor.
- Chosen: the eligible sections with the smallest, the median (lower median when the count is even) and the largest
  total cell count.
- Written to `framework/results/v2_phase4/section_selection.json` before the run.

## Semi-synthetic data (`cellalibi.simulate.SemiSyntheticSection`)

What is kept real:
- coordinates;
- every cell's real `consensus_T` label;
- every cell's real library-size factor.

Simulated panel:
- FN1, ITGB1, the top 5 markers of every type with >= 20 cells, and 10 high-mean background genes.
- Per-type means are the real mean raw counts in that section.
- Counts ~ NB(mean = size factor x type mean, theta = 3), independent across cells.

Contact graph: symmetric 6-NN.

| Setting | Planted truth |
|---|---|
| S0 | none |
| S1 f = 1.25 / 1.5 / 2 | ITGB1 mean x f only in myeloid cells touching a CAF with FN1 count > 0 |
| S2 | FN1 mean x 2 in all CAFs and ITGB1 mean x 2 in all myeloid cells (no spatial dependence) |
| S3 s = 50 / 200 / 400 um | Square tiles of side s (random offset); each tile is "high" with probability 0.5. In high tiles, FN1 in CAFs x 2 and ITGB1 in myeloid x 2. No adjacency dependence. |
| S4 p = 0.1 / 0.2 / 0.3 | S0 counts; each CAF loses a binomial fraction p of every gene's transcripts, and each moved transcript goes to a uniformly chosen contact neighbour (simulated segmentation error). Total transcripts are conserved. No biology. |
| S5 n = 25 / 50 / 100 | Myeloid cells subsampled (without replacement) to n; all other cells kept; graph rebuilt; S1 at 2x |

Replication:
- 100 replicates per setting per section, so 14 settings x 3 sections x 100 = 4,200 audits.
- Each replicate has fresh simulated counts (and fresh tiles or subsamples) from a seed derived from
  (section, setting, replicate).
- If memory or time forces fewer replicates, this is logged.

Audit configuration: package defaults. 1000 permutations (minimum p = 1/1001), alpha = 0.05, power floor 0.5 at
1.5x, 25 receivers, 10 edges, tile 100 um, ring 4 d.

## Expected verdicts (declared now)

| Scenario | Correct verdict | Counted as "fooled" |
|---|---|---|
| S0 | anything except SUPPORTED | SUPPORTED |
| S1 | SUPPORTED; NOT_TESTABLE is acceptable (not an error, but not a detection) | any EXPLAINED_BY_* or NOT_SUPPORTED is a miss |
| S2 | EXPLAINED_BY_ABUNDANCE | SUPPORTED |
| S3 | EXPLAINED_BY_ZONATION | SUPPORTED |
| S4 | EXPLAINED_BY_SPILLOVER | SUPPORTED |
| S5 | NOT_TESTABLE, or SUPPORTED (the effect is real) | EXPLAINED_BY_* / NOT_SUPPORTED (a false refutation) |

## Reported

- Confusion matrix (true setting x verdict), pooled over sections and per section.
- Per-setting accuracy, i.e. the fraction with the exact correct verdict as tabulated above.
- S1 sensitivity, P(SUPPORTED), by effect size, with Wilson 95 % CI. Also P(SUPPORTED | power PASS).
- False-SUPPORTED rate with Wilson 95 % CI for S0, S2, each S3 level and each S4 level.
- How often each check fails in each setting (the multi-label view).
- Runtime per audit.

## Pass criterion (declared now)

The benchmark PASSES only if both hold, using replicates pooled over the 3 sections (n = 300 per setting):
1. The false-SUPPORTED rate is <= 0.10 in **every** artefact setting: S0, S2, S3 at 50 / 200 / 400 um, S4 at
   0.1 / 0.2 / 0.3.
2. S1 sensitivity at 2x is >= 0.80.

Otherwise it FAILS, and every setting that fools it (false-SUPPORTED > 0.10), or misses (S1 2x < 0.80), is listed.
Per-section results are also reported. A setting that passes pooled but fails in one section is listed as a
section-specific weakness.

## Known risks, stated before running

- S3 at 50 um is close to contact scale. The within-tile and ring controls may not separate it from contact.
- S4 only creates a false signal if CAF ITGB1 exceeds myeloid ITGB1 in the section.
- S1 power depends on how many myeloid cells touch FN1+ CAFs.

## Frozen engine and execution details

- Engine: `framework/src/cellalibi` at commit fc9bab8 (end of Phase 3).
- Runner: `framework/scripts/v2/phase4_benchmark.py`.
  - 4 worker processes, 1 BLAS thread each.
  - One resumable CSV chunk per (section, setting) in `framework/results/v2_phase4/chunks/`.
  - Seeds = CRC32 of (section, setting, level, replicate).
- Section selection: `framework/scripts/v2/phase4_select_sections.py` writes `section_selection.json`.
- Before the full run, a smoke test runs 1 replicate of S0 and S4 0.3 on the smallest section, writing to a scratch
  directory outside the results tree. It only checks that the code runs and how long it takes. It is not a result,
  and the engine is not changed because of it.

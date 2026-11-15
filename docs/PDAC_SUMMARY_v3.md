# PDAC Extension v3 — Summary (2026-11-15)

## Does the proximity-specific signal survive spatially constrained nulls, and could this audit detect a real effect of the size we would care about?

**It survives the constrained nulls, but the audit could not have told us much either way.** Under block permutation within 200, 100 and 50 um tiles and under a distance-stratified null, the CXCL12→CXCR4 CAF→T ratio Q = E/M does not shrink toward zero; its z grows because the null variance falls while the null mean of Q stays near 1. However, the simulation calibration shows that the Q audit has **zero power** in sections with about 25 T cells (two of the three calibration sections), and that in the one large section every null except block200 calls a pure 200 um expression zonation "positive" in 42–88% of simulations. The only null that is not fooled detects a planted 2x contact effect only half the time. So the positive trend is not a label-zonation artefact at ≥ 50 um, but the same scores cannot separate contact from expression zonation, and most sections are below the detection floor. This is an audit and diagnostic statement, not a biological or clinical claim.

Branch `thesis-v3` (from `pdac-v2` @ 1173e52, with `main` @ ed86917 merged). Every number below is in `results/pdac/v3_phase6/v3_key_numbers.csv` (key, value, source file, column). Plans were committed before computation in `docs/v3_plans/`. The full record, including the one crash, is in `docs/PDAC_PROGRESS_LOG.md` (section "PDAC v3"). v1: `docs/PDAC_SUMMARY.md`; v2: `docs/PDAC_SUMMARY_v2.md`. p-values are add-one, (k+1)/(n+1), throughout.

## Results

### Phase 1 — contact vs zonation (16 sections, 3 annotations, 5 nulls, 500 permutations)

Nulls: `global` (v2 reference; reproduces Phase B3 exactly), `block200` / `block100` / `block50` (labels permuted within square tiles), `dist` (labels permuted within signed-distance bands from the tumour compartment of the same annotation). Source: `results/pdac/v3_phase1/headline_summary_by_null.csv`, `tumour_level_summary_by_null.csv`, `nonbg_summary_by_null.csv`, `phase1_summary.json`. Figure: `figures/pdac/v3_phase1/phase1_Qz_vs_null`.

Headline CXCL12→CXCR4 CAF→T, consensus_T labels (observed median Q 1.320, fixed by construction):

| null | median z_Q (16 sections) | sections Q p<0.05 (BH) | median null mean Q | median null sd Q | tumours z_Q>0 (of 13) | tumours BH-sig | median tumour z_Q |
|---|---|---|---|---|---|---|---|
| global | 0.36 | 2 (0) | 1.067 | 0.433 | 10 | 0 | 0.43 |
| block200 | 0.84 | 3 (1) | 1.009 | 0.293 | 11 | 2 | 0.86 |
| block100 | 0.88 | 6 (0) | 0.954 | 0.257 | 12 | 0 | 0.81 |
| block50 | 1.35 | 6 (2) | 1.056 | 0.236 | 11 | 3 | 1.30 |
| dist | 0.38 | 3 (0) | 1.136 | 0.331 | 10 | 0 | 0.58 |

- **Shrink verdict** (rule declared in the plan): `shrinks_toward_zero = False` for consensus_T, marker and author_mapped. Marker: median tumour z_Q 0.90 (global) → 1.23 (block50); author_mapped: 0.78 → 1.49.
- The null mean of Q is close to 1 under every null, so label zonation at the 50–200 um scale does not reproduce the observed Q ≈ 1.3. Median z_E falls from 14.7 to 8.3 (global → block50), showing that the constrained nulls do remove part of the compositional component of E.
- Full non-background catalogue (7 axes, all sender→receiver rows, consensus_T, 17,500 rows): fraction of rows with Q p<0.05 is 0.138 (global) → 0.131 → 0.126 → 0.122 (block50) → 0.110 (dist); BH-significant Q rows 1,234 → 565 (block50) → 671 (dist). Spearman of z_Q with the global-null z_Q: 0.79 / 0.71 / 0.57 (block200/100/50), 0.86 (dist). The constrained nulls re-rank the catalogue but do not bring the small-p excess near 5%.

### Phase 2 — calibration by simulation (planted effects on real tissue)

Three real sections by a fixed rule: T6-b (smallest; 26 T cells, 83 CAF→T edges), T2-b (lower median; 23 T, 74 edges), U5-a (largest; 261 T, 448 edges). 100 simulations x 500 permutations per setting. PROXIMITY: CXCR4 counts of T cells adjacent to CXCL12+ CAFs redrawn as Poisson(f·c). ZONATION-ONLY: CXCL12 and CXCR4 counts of all cells in random 200 um "high" tiles redrawn as Poisson(f·c). The planted perturbation is a computational stand-in for a stimulation experiment; nothing experimental is claimed. Source: `results/pdac/v3_phase2/calibration_table.csv`, `calibration_deliverable.csv`. Figure: `figures/pdac/v3_phase2/phase2_calibration`.

**Calibration table, score Q, section U5-a** (fraction of 100 simulations with p < 0.05):

| null | control f=1 | power f=1.1 | f=1.25 | f=1.5 | f=2 | zonation FPR f=1.1 | f=1.25 | f=1.5 | f=2 |
|---|---|---|---|---|---|---|---|---|---|
| global | 0.02 | 0.12 | 0.34 | 0.78 | 1.00 | 0.42 | 0.44 | 0.56 | 0.54 |
| block200 | 0.00 | 0.00 | 0.00 | 0.02 | 0.49 | 0.10 | 0.06 | 0.06 | 0.08 |
| block100 | 0.06 | 0.18 | 0.47 | 0.85 | 1.00 | 0.56 | 0.60 | 0.73 | 0.74 |
| block50 | 0.00 | 0.01 | 0.11 | 0.52 | 1.00 | 0.35 | 0.32 | 0.46 | 0.43 |
| dist | 0.18 | 0.39 | 0.73 | 0.92 | 1.00 | 0.74 | 0.76 | 0.88 | 0.86 |

- **T6-b and T2-b**: Q detection rate 0.00 for every null, every f up to 2, both scenarios. The headline row is not testable at ~25 T cells.
- **C, M, E**: detection rate 1.00 (E: 0.75–0.96) already at f = 1 and under zonation. They are saturated by the real composition signal and cannot distinguish contact from zonation or from nothing.
- **Pooled deliverable** (`calibration_deliverable.csv`): no score/null reaches power 0.8 at any f (max pooled Q power 0.33). Fooled by zonation (pooled FPR > 0.10 at f = 1.5): Q under global (0.19), block100 (0.24), block50 (0.15), dist (0.29); not fooled: block200 (0.02), which is the null matching the planted zonation scale.
- **Which combination can detect what**: Q with global, block100 or dist detects a 1.5x contact effect in the largest section (power 0.78–0.92) but is fooled by 200 um zonation; Q with block200 is not fooled but detects only a 2x effect, at 0.49. No combination is both sensitive at ≤ 1.5x and robust to zonation.

### Phase 3 — effect sizes and tumour-level model

Source: `results/pdac/v3_phase3/Q_block_bootstrap_summary.csv`, `Q_block_bootstrap_by_section.csv`, `mixed_model_logQ.csv`, `z_vs_edge_count_summary.csv`. Figure: `figures/pdac/v3_phase3/phase3_effect_sizes`.

- **Block-bootstrap 95% CI for Q** (100 um tiles, B = 1000), headline: consensus_T — CI above 1 in 4/16 sections, below 1 in 1/16, median Q 1.320, median CI width 0.767, median 448 edges. marker 8/16 above (1 below); author_mapped 7/16 above (0 below). Among the 7 non-background axes CAF→T (consensus_T), CXCL12-CXCR4 is the only one with CIs above 1; MIF-CD74 is below 1 in 6/16.
- **Mixed model log Q ~ 1 + (1 | tumour)** (exploratory; 16 sections, 13 tumours, 10 singleton tumours; statsmodels REML, converged): geometric-mean Q 1.242 [1.013, 1.522] consensus_T (tumour variance 0.007, boundary warning), 1.315 [1.176, 1.471] marker, 1.348 [1.219, 1.491] author_mapped. The section-level t-interval is nearly identical, because most tumours contribute one section.
- **z vs edge count** (7 non-bg axes, 17,500 rows, global null, consensus_T): Spearman(|z_E|, n_edges) = 0.56; after dividing by sqrt(n_edges) the signed correlation moves from −0.12 to 0.25. For the headline across 16 sections, Spearman(z_E, n_edges) = 0.71 while z_Q is edge-count independent (−0.08). Normalised headline z_E: median 14.68 raw → 0.69 per sqrt(edge).

### Phases 4 and 5 — BLOCKED

- **Phase 4 (GSE310352 replication)**: `LICENSE_CONFIRMED_GSE310352.txt` absent → BLOCKED, licence unconfirmed (`results/pdac/v3_phase4/phase4_gate_status.json`). Nothing was extracted or analysed. Also outstanding: um/px unverified, no author labels, patient overlap with the Mendeley set unverified.
- **Phase 5 (GSE277783 programs)**: no GEO SOFT/series-matrix file under `data/`; organ site unknown for 49/57 slides → BLOCKED (`results/pdac/v3_phase5/phase5_gate_status.json`). Program gene lists are pre-declared in `docs/v3_plans/phase5_visium_programs_gse277783.md`.

## Calibration table (compact)

Pooled over the 3 calibration sections, Q only (`calibration_deliverable.csv`); per-section values are above and in `calibration_table.csv`.

| null | power f=1.1 / 1.25 / 1.5 / 2 | zonation FPR f=1.1 / 1.25 / 1.5 / 2 | min f with power ≥ 0.8 | fooled (FPR > 0.10 at f=1.5) |
|---|---|---|---|---|
| global | 0.04 / 0.11 / 0.26 / 0.33 | 0.14 / 0.15 / 0.19 / 0.18 | none | yes |
| block200 | 0.00 / 0.00 / 0.01 / 0.16 | 0.03 / 0.02 / 0.02 / 0.03 | none | no |
| block100 | 0.06 / 0.16 / 0.28 / 0.33 | 0.19 / 0.20 / 0.24 / 0.25 | none | yes |
| block50 | 0.00 / 0.04 / 0.17 / 0.33 | 0.12 / 0.11 / 0.15 / 0.14 | none | yes |
| dist | 0.13 / 0.24 / 0.31 / 0.33 | 0.25 / 0.25 / 0.29 / 0.29 | none | yes |

Pooled power is capped at 0.33 because two of the three sections contribute zero detections at every f.

## Claims register

| # | claim | evidence file | status |
|---|---|---|---|
| 1 | The headline Q signal does not shrink toward zero under block (200/100/50 um) or distance-stratified nulls; median tumour z_Q rises from 0.43 (global) to 1.30 (block50), consensus_T | `results/pdac/v3_phase1/tumour_level_summary_by_null.csv`, `phase1_summary.json` | supported (audit level) |
| 2 | Label zonation at ≥ 50 um does not reproduce the observed Q ≈ 1.3 (null mean Q 0.95–1.14 under every null) | `results/pdac/v3_phase1/headline_summary_by_null.csv` | supported (audit level) |
| 3 | Tumours with BH-significant Q under constrained nulls: 0–3 of 13 depending on null and annotation | `results/pdac/v3_phase1/tumour_level_summary_by_null.csv` | exploratory (depends on null choice; no pre-registered null) |
| 4 | The constrained nulls do not remove the catalogue-wide excess of small Q p-values (0.11–0.14 of rows at p<0.05) | `results/pdac/v3_phase1/nonbg_summary_by_null.csv` | supported (audit level) |
| 5 | The Q audit has zero power in sections with ~25 T cells, even for a planted 2x contact effect | `results/pdac/v3_phase2/calibration_table.csv` | supported (simulation) |
| 6 | In the largest section, Q under global/block100/dist detects a 1.5x contact effect (power 0.78–0.92) but calls 200 um expression zonation positive in 42–88% of simulations | `results/pdac/v3_phase2/calibration_table.csv` | supported (simulation; zonation at one scale only) |
| 7 | block200 is the only null not fooled by the planted zonation, at the cost of detecting only a 2x effect at power 0.49 | `results/pdac/v3_phase2/calibration_deliverable.csv`, `calibration_table.csv` | supported (simulation; scale matches the planted zonation by design) |
| 8 | C, M and E are saturated (detection rate ~1.0 at f = 1) and cannot distinguish contact from zonation | `results/pdac/v3_phase2/calibration_deliverable.csv` | supported (simulation) |
| 9 | Geometric-mean Q across tumours is 1.24 [1.01, 1.52] (consensus_T) with a random intercept per tumour | `results/pdac/v3_phase3/mixed_model_logQ.csv` | exploratory (16 sections, 13 tumours, boundary warning) |
| 10 | Per section, the block-bootstrap CI for the headline Q lies above 1 in 4/16 (consensus_T), 8/16 (marker), 7/16 (author_mapped) sections | `results/pdac/v3_phase3/Q_block_bootstrap_summary.csv` | supported (audit level) |
| 11 | \|z_E\| grows with edge count (Spearman 0.52–0.56); z_Q does not (headline −0.08) | `results/pdac/v3_phase3/z_vs_edge_count_summary.csv` | supported (audit level) |
| 12 | A proximity-specific CXCL12→CXCR4 CAF→T effect exists in PDAC tissue | — | not supported (not testable with this audit: see claims 5–7) |
| 13 | The result replicates in an independent cohort (GSE310352) | `results/pdac/v3_phase4/phase4_gate_status.json` | not supported (BLOCKED, not run) |
| 14 | Primary vs metastasis differences in gene programs (GSE277783) | `results/pdac/v3_phase5/phase5_gate_status.json` | not supported (BLOCKED, not run) |

The same table is in `results/pdac/v3_phase6/claims_register.csv`.

## Limitations

- The nulls are constrained on labels, not on expression. The simulation shows that expression zonation at 200 um fools every null except the one tiled at that scale; zonation at other scales was not simulated.
- Two of the three calibration sections have 23–26 T cells; the headline row has zero power there. Fifteen of the 16 sections have fewer than 1,000 CAF→T edges under consensus_T (6 have fewer than 300, and 5 have fewer than 100 T cells; `results/pdac/v3_phase3/Q_block_bootstrap_by_section.csv`), so the calibration result for the small sections applies to most of the dataset.
- `consensus_T` is post-hoc (v2) and weakly protein-supported; marker and author labels give larger Q and more sections with CIs above 1.
- Tumour-level inference rests on a donor-id mapping (16 sections → 13 tumours, 10 singletons); the mixed model is exploratory.
- Block tiles are anchored at the section origin; no random offset was used for the nulls (the simulation's zonation used random offsets).
- The 50 um null has a median of 5 cells per tile and many singleton tiles, which limits how much it can shuffle.
- Phases 4 and 5 did not run; nothing in this summary is replicated outside the Mendeley PDAC-SMI dataset.

## BLOCKED items

- Phase 4: licence for GSE310352 unconfirmed (gate file absent); um/px and patient overlap unverified; no author labels.
- Phase 5: no GEO SOFT/series-matrix file with site labels for GSE277783 in `data/`.
- Push of `thesis-v3`: not blocked (pushed on the first attempt at 818d39d; see the final log entry).

## Next steps (at most 5)

1. Add an expression-level null (e.g. permute expression within tiles while keeping labels) and re-run the calibration with zonation planted at 50, 100 and 400 um, to find a null that is both sensitive at 1.5x and not fooled.
2. Pre-register one null (block200 is the current best candidate) and one annotation before any further headline test.
3. Report the headline only for sections with ≥ 100 T cells and ≥ 300 CAF→T edges (the calibration floor), and state the floor in every table.
4. Obtain the GSE310352 licence confirmation and sample metadata, then run Phase 4 with the pre-registered null.
5. Obtain the GSE277783 series matrix with site labels, then run Phase 5 with the pre-declared gene lists.

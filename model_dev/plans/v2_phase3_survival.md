# v2 Phase 3 plan: survival analysis (amendment to docs/translational_plans/phase3.md)

Committed before any survival model is fitted, i.e. before any Cox model, PH diagnostic or Kaplan-Meier curve on patient outcomes. The committed translational Phase 3 plan (`docs/translational_plans/phase3.md`, including its "Computational resolution" section) is executed **unchanged**, apart from the amendments below.

The code reuses, without editing:
- `src/translational/survival.py`: scoring, fixed Cox models with BH, the validated PH diagnostic, Kaplan-Meier steps and correlations;
- `src/translational/phase3.py`: gene sets, expression loading with the SHA-256 identity check, and cohort configs.

Cohorts (prepared in translational Phase 2, `results/translational/phase2/20261116T043554395029Z/`):
- primary: strict-histology TCGA, 141 patients / 81 events;
- validation: CPTAC explicit ductal, 129 / 72;
- sensitivity: all-primary cohorts, each its own BH family, as in the plan.

## Amendments and why

**A1. Three additional pre-declared signatures from the model work.** Gene lists are copied verbatim from `results/chain/phase1b/programs_used.json` (the program definitions used as v1 model targets):
- `Macrophage_C1Q_TAM`: C1QA C1QB C1QC CD163 MRC1 SELENOP (key `myeloid:C1Q_TAM`);
- `Macrophage_SPP1_TAM`: SPP1 GPNMB CD9 FABP5 LGALS3 (key `myeloid:SPP1_TAM`);
- `Macrophage_C1Q_minus_SPP1` = C1Q mean gene z minus SPP1 mean gene z, then standardised. This is the same construction as the plan's Immune_exclusion_proxy. Each component must meet the plan's coverage rule (>= 5 genes and >= 50 % of the declared genes).
- No other gene list is added, and no list is changed after seeing any outcome.
- The v2 Phase 1 spillover result is not used to edit these lists: bulk RNA has no cell segmentation. Its finding that SELENOP is expressed in CAFs at about 0.8x the isolated-macrophage level is carried as an interpretation caveat: in bulk, the C1Q signature also partly reflects stroma.

**A2. Multiplicity.** As requested for this phase, the primary BH family is **all 9 signatures per cohort and model family**: the plan's six plus the three above, i.e. 9 tests per family.
- In the strict TCGA cohort this gives 3 model families x 9 = 27 primary models.
- The plan's original six-signature BH within each family is also reported, as a secondary column, for continuity.
- All-primary sensitivity cohorts form separate BH families. PH tests use BH across all PH tests per cohort and model family, as in the plan.

**A3. CPTAC validation.** The same models are applied, unchanged.
- For every signature and model family: CPTAC HR per SD with 95 % CI, and same-direction yes/no (plan).
- A signature's adjusted association "holds in CPTAC" only if the direction is the same AND the CPTAC purity-adjusted CI excludes 1. The direction-only result is also reported.
- Purity measures differ between cohorts: TCGA ABSOLUTE DNA purity versus CPTAC pathology cellularity.

**A4. Power statement** (no outcome data used beyond event counts):
- Schoenfeld approximation for a standardised continuous covariate: minimum detectable HR per SD at 80 % power = exp((z_{1-a/2} + z_{0.8}) / sqrt(D x (1 - R2))).
  - D = events in the fitted model;
  - R2 = R-squared of the signature on the model's adjustment covariates (0 for unadjusted).
- Reported for two-sided a = 0.05, and for a = 0.05/9 (Bonferroni bound for the 9-test family).

**A5. Outputs.**
- Results go to `results/translational/phase3/<UTC stamp>/`, with a run-local `MANIFEST_v2.csv` (path, bytes, SHA-256).
- The global `results/translational/MANIFEST.csv` and `docs/TRANSLATIONAL_PROGRESS_LOG.md` are not modified; progress is logged in `model_dev/reports/PROGRESS_LOG_v2.md`.
- Figure: KM panels at the cohort median (descriptive only, no inferential log-rank tests, no optimised cut-points), a forest plot of adjusted HRs in both cohorts, and the purity/inter-signature correlation matrix.

## Decision rule (from the plan, with A2)
- **Relates to survival after adjustment:** a signature qualifies only if the age + sex + stage + purity model is estimable (>= 10 events per coefficient, full rank, converged) AND its BH q < 0.05 in the 9-test family of the strict TCGA cohort, with PH status disclosed. A PH failure makes the constant-HR reading exploratory.
- **No signature passes:** the negative result is reported.
- **No final model can be fitted:** the question is reported as BLOCKED or inconclusive.
- **Purity proxy:** a signature is called mostly a purity proxy if |Spearman rho| with purity >= 0.7.
- Results are prognostic expression associations: not causal effects, not clinical recommendations, and not validation of a spatial mechanism.

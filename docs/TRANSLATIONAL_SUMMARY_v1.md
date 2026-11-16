# Translational v1 — draft PR hand-off

## Is the hypoxia association robust?

**Not robust under the committed decision rule.** The saved v3 result was reproduced exactly, including all 54 raw-slide hypoxia scores. The effect remains positive under every declared sensitivity, but tumour-rich/host-filtered spots give paired p=0.068115234375 (13 pairs), and common 2,000-UMI depth gives p=0.193359375 (10 pairs). Both fail the required p<0.05 criterion. This does not establish absence of an association: the depth restriction also reduces the paired sample size.

Excluding CytAssist retains 13 pairs (p=0.0478515625). All 12 leave-one-gene-out effects are positive, and the separately curated MSigDB Hallmark set is positive (p=0.00048828125). The original effect exceeds most size/expression-matched random sets (one-sided add-one p=0.00999000999000999). These checks do not override the prespecified failed b/c gates. Host/tumour marker rules are screening proxies, not validated annotations.

Sources: `../results/translational/phase1/20261116T040752758116Z/robustness_tests.csv`, `../results/translational/phase1/20261116T040752758116Z/raw_reproduction.csv`, `../results/translational/phase1/20261116T040752758116Z/random_set_effects.csv`. All values are generated results; patients are the inference unit.

## Do spatial program signatures relate to survival after adjusting for purity and stage?

**Not yet tested.** The user redirected this task to PR preparation before Phase 3 execution. Fixed signature/model code and plans are included for review, but there are no completed Cox, PH-diagnostic or Kaplan-Meier outcome analyses and no survival-association claim.

## Cohorts prepared

| Strict clinical-histology subset | Patients with expression and valid OS | Events | Available purity/cellularity |
|---|---|---|---|
| TCGA, ICD-O-3 8500/3 | 141 | 81 | 122 (ABSOLUTE) |
| CPTAC, explicit ductal diagnosis | 129 | 72 | 127 (pathology cellularity) |

All-primary sensitivity cohorts are also prepared. Source: `../results/translational/phase2/20261116T043554395029Z/cohort_summary.csv`; missingness and every exclusion are adjacent. The TCGA published curated-list supplement hit CAPTCHA; the strict clinical-code fallback was used rather than treating the broad universal PAAD label as histology verification. CPTAC semicolon-separated pathology review values are averaged; their percentage scale is inferred and differs from TCGA ABSOLUTE DNA purity. Treatment is not modelled. Bulk and spatial expression are different measurements and cannot establish a spatial mechanism or clinical prediction.

Source access, data types, exact URLs, counts, citation requirements and restrictions are in `../results/translational/cohort_vetting.csv` and the acquisition records. Raw downloads remain outside Git under `F:/Thesis/data/tcga/`.

## Claims register

| Claim | Status | Evidence |
|---|---|---|
| The v3 hypoxia calculation reproduces | Supported computational result | Baseline and raw reproduction tables |
| The association is robust to every declared a-d check | Not supported under the fixed rule | Robustness tests: b/c fail |
| All sensitivity directions remain positive | Supported descriptive result | Robustness tests |
| The observed effect is unusual relative to matched random sets | Exploratory expression association | Random-set effects and summary |
| Spatial program signatures predict survival after purity/stage adjustment | Not assessed | Phase 3 not run |
| Expression data establish physiological hypoxia, causation or clinical utility | Not supported | Outside this design |

## Checklist at user-requested PR hand-off

- Phase 1: DONE — NOT ROBUST by the declared rule.
- Phase 2: DONE — two cohorts prepared; curated-list download blocked, clinical-field fallback documented.
- Phase 3: NOT RUN — stopped for PR preparation; code is unvalidated on outcome data.
- Phase 4: NOT RUN — no chain integration performed. The requested ranked-hypotheses path was absent at the inspected remote snapshot.
- Phase 5: PARTIAL — this draft hand-off, source register and progress log only; full post-survival synthesis remains pending.

## Next steps

1. Review tumour/host proxies and the sensitivity estimands before interpreting composition effects.
2. Review the TCGA clinical-code cohort and CPTAC cellularity/endpoint assumptions.
3. Validate and run the committed fixed survival models and PH diagnostics.
4. Consume a chain hand-off only when explicit ranked hypotheses and gene lists are available.
5. Complete the final synthesis after outcome analyses; retain all negative and blocked results.

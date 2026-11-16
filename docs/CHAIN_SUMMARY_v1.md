# Chain-v1 summary: tracing A -> B -> C signalling chains in PDAC CosMx

Branch `chain-v1`. Every number below is read from a file in `results/chain/`; the claims register is `results/chain/phase7/claims_register.csv` and the key numbers with sources are in `results/chain/phase7/key_numbers.csv`. The data measure **expression, not mutation**. Results are associations and ranked hypotheses, not causal or mechanistic claims, and the tumour (patient) is the unit of inference. Plans were committed before each phase was computed (`docs/chain_plans/`); deviations are logged in `docs/CHAIN_PROGRESS_LOG.md`.

## Lead question
**Do neighbours predict cell state beyond the cell's own type, and is there support for any A -> B -> C chain within the calibrated range?**

- **Neighbours predict state: yes, but almost entirely through same-type neighbours.**
  - In Mendeley (13 tumours, leave-one-tumour-out), neighbourhood features raise held-out R2 beyond own covariates and a 200 um zonation control for 6 of 7 receiver programs. For tumour classical the ridge delta R2 is 0.141 [0.105, 0.164].
  - The gain comes from same-type neighbour programs (drop-group delta R2 0.003-0.117) and composition. Programs of *other* cell types add at most 0.0005.
  - T-cell cytotoxic state is not predicted (ridge 0.0005 [-0.0027, 0.0046]).
- **Chains: no support, and no calibrated range in which a chain could have been supported.**
  - All 44 pre-declared chain x program combinations are "not supported" under both annotations.
  - The simulation calibration shows the chain test has power <= 0.04 to detect planted chains of 1.1x-2x per hop, on 3 real sections and pooled. The calibrated range is therefore empty, and no chain call could fall inside it.

## Results by phase

### Phase 0: related work (`docs/chain_plans/related_work.md`)
- **Methods covered.** MISTy, NCEM, CellNEST and HoloNet were fetched; HoloNet as abstract only.
- **What they lack.** MISTy and NCEM model one-hop neighbourhood dependence. CellNEST counts two-hop ligand-receptor co-occurrences but does not test mediation. None calibrates a chain test against zonation.
- **Fetches that were blocked** were logged and not worked around: a PMC CAPTCHA, a bioRxiv 429 and a PubMed cookie wall.

### Phase 1: cell states (`results/chain/phase1`, `results/chain/phase1b`)
- **Original rule.** The pre-declared section-wide background rule flags 630/960 panel genes and removes every T-cell program gene, so T-containing chains are **not testable** under that rule.
- **Amendment 1.** The flaw is structural: rare-type genes cannot reach 2x negative-probe detection across all cells. Amendment 1 (type-aware rule) was committed before it was computed. It keeps 7 programs: tumour classical (8 genes) and basal_like (10); CAF myCAF (8) and a CXCL12 readout; T cytotoxic (9); myeloid C1Q_TAM (6) and SPP1_TAM (5).
- **Dropped under both rules:** tumour EMT, iCAF, T exhaustion and myeloid inflammatory.
- **Validity.** Classical minus basal_like separates the authors' CLS and BSL malignant cells with a median AUC of 0.858 across 16 sections.
- **Chain coverage.** The primary chain tumour -> CAF -> T passes the floor in 14/16 sections (11/13 tumours). The 5 secondary chains, chosen by counts only, contain no T cells.

### Phase 2: does the neighbourhood predict state? (`results/chain/phase2`)
| Receiver program (author labels) | Ridge delta R2 [95 % tumour CI] | HGB delta R2 |
|---|---|---|
| tumour:classical | 0.141 [0.105, 0.164] | 0.151 |
| tumour:basal_like | 0.057 [0.034, 0.074] | 0.060 |
| CAF:myCAF | 0.034 [0.021, 0.045] | 0.039 |
| CAF:CXCL12 (single gene) | 0.0085 [0.0027, 0.0147] | 0.015 |
| myeloid:C1Q_TAM | 0.030 [0.016, 0.045] | 0.038 |
| myeloid:SPP1_TAM | 0.056 [0.042, 0.065] | 0.068 |
| T:cytotoxic | 0.0005 [-0.0027, 0.0046] | 0.0057 [-0.0016, 0.0087] |

- **Spillover covariates** change the gains little (classical 0.141 -> 0.123).
- **Graph model.** The pre-declared GNN rule triggered, so an SGC-style graph model was fitted. It adds nothing: SPP1_TAM 0.048 and classical 0.135, against ridge 0.056 and 0.141. A first version had a self-leakage bug (implausible delta R2 = 0.67). It was found, fixed, and the invalid files are kept and labelled.

### Phase 3: chain test (`results/chain/phase3`)
- **Design.**
  - 6 chains (primary plus 5 secondary) and 44 program combinations.
  - Nulls: 200 expression and 200 label block permutations per section (200 um tiles, random offsets).
  - Paired permutations aggregate the result at tumour level.
  - BH within each family: hop C<-B 22 tests, hop B<-A 20, indirect 44.
- **Calls:** **44/44 not supported** (author labels) and **44/44 not supported** (marker labels). No forward chain, reversed chain, random-B or shuffled-B control meets criteria 1-3 under author labels.
- **Primary chain tumour -> CAF -> T:**
  - hop T <- CAF: beta 0.003 to 0.010;
  - indirect: -0.008 to 0.008;
  - every CI spans 0.
- **Amplification.** Per-hop effects are all near zero (|beta| <= 0.03); |beta_CB| / |beta_BA| has a median of 0.80. There is no amplification pattern.
- **Strongest single hop (exploratory).** CAFs adjacent to basal-like tumour cells have lower CXCL12: beta -0.028 [-0.039, -0.018], with the same sign in 11/11 tumours. Under the expression null p = 0.005 and q = 0.050; under the label null p = 1.0. The spillover covariates do not change it (-0.027). In GSE310352 the same hop is -0.030 [-0.047, -0.017], with p = 0.0099 under both nulls (Phase 6). This is a single-hop adjacency association, not a chain.
- **Label null.** Its means sit far from 0 (-0.24 to 0.02 for hop betas), because relabelled cells carry another type's expression. The pre-declared two-sided |statistic| test is therefore ill-posed for this null. Reported here as a limitation; the rule was not changed.

### Phase 4: calibration by simulation (`results/chain/phase4`)
Planted triple: tumour:basal_like -> CAF:myCAF -> T:cytotoxic. 100 simulations per setting, 100 permutations each.

| Section (analysed C cells, median) | Scenario | f | Chain-call rate | Hop C<-B | Hop B<-A |
|---|---|---|---|---|---|
| T4-a (230) | chain | 1.1 / 1.25 / 1.5 / 2.0 | 0.01 / 0.00 / 0.01 / 0.02 | 0.11 / 0.11 / 0.11 / 0.25 | 0.83 / 0.99 / 1.00 / 1.00 |
| T4-a | zonation | 1.5 / 2.0 | 0.00 / 0.05 | 0.71 / 1.00 | 0.70 / 0.97 |
| T4-a | direct A->C | 1.1-2.0 | 0.00-0.01 | 0.06-0.16 | 0.21-0.25 |
| T4-a | null | 1 | 0.00 | 0.17 | 0.23 |
| T1-b (53) | chain | 1.1 / 1.25 / 1.5 / 2.0 | 0.02 / 0.02 / 0.04 / 0.01 | 0.14 / 0.12 / 0.12 / 0.20 | 0.60 / 1.00 / 1.00 / 1.00 |
| T1-b | zonation | 1.5 / 2.0 | 0.00 / 0.00 | 0.02 / 0.00 | 0.22 / 0.79 |
| U2-b (65) | chain | 1.1 / 1.25 / 1.5 / 2.0 | 0.00 for all | 0.00-0.01 | 0.16 / 0.94 / 1.00 / 1.00 |
| pooled 3 sections | chain | 1.1-2.0 | 0.00 for all | 0.00-0.01 | 0.99-1.00 |

- **Why the chain is undetectable.** The upstream hop is detected, but the downstream hop is not. Only a minority of qualifying T cells are perturbed: 129 of 3,088 in T4-a, 27 of 610 in T1-b and 46 of 106 in U2-b. The neighbour-mean CAF exposure tracks those perturbed neighbours only weakly.
- **Zonation fools the single hops** in the largest section (T4-a at 2x: C<-B 1.00, B<-A 0.97). The full chain call is protected only by the indirect criterion.
- **Calibrated range: none.** No (f, section) reaches power >= 0.8 with false calls <= 0.10. Every Phase 3 call is therefore outside the calibrated range.

### Phase 5: candidate signals (`results/chain/ranked_hypotheses.csv`, de-duplicated `ranked_hypotheses_dedup.csv`)
**0 of 88 hops meet the support rule**, so every candidate is "hop not supported" and unranked. The top 10 rows by |hop effect| are listed only to show what was examined. They are **candidates from prior ligand-receptor knowledge, not evidence**, and all are outside the calibrated range.

| Chain | Hop | Programs | Candidate LR | Hop effect [95 % CI] | Tumours same sign | Hop support |
|---|---|---|---|---|---|---|
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | TGFB1-TGFBR1 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | TGFB1-TGFBR2 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | JAG1-NOTCH1 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | MIF-CD74 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | LGALS9-HAVCR2 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | SPP1-CD44 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| CAF->tumour->myeloid | C<-B | tumour:classical -> myeloid:C1Q_TAM | FN1-ITGB1 | -0.0304 [-0.0816, 0.0065] | 7/12 | no |
| tumour->CAF->T | B<-A | tumour:basal_like -> CAF:CXCL12 | TGFB1-TGFBR1 | -0.0283 [-0.0385, -0.0181] | 11/11 | no |
| tumour->CAF->T | B<-A | tumour:basal_like -> CAF:CXCL12 | TGFB1-TGFBR2 | -0.0283 [-0.0385, -0.0181] | 11/11 | no |
| tumour->CAF->T | B<-A | tumour:basal_like -> CAF:CXCL12 | JAG1-NOTCH1 | -0.0283 [-0.0385, -0.0181] | 11/11 | no |

### Phase 6: consistency in GSE310352 (`results/chain/phase6`; not replication, because patient independence is UNVERIFIED)
- **Setup.** 8 slides form 6 units (slides 5, 5b2 and 5b3 grouped, an UNVERIFIED assumption). Labels are marker labels. Tiles are scaled by nearest-neighbour distance, since um/px is UNVERIFIED. All 7 programs keep >= 5 usable genes (or the CXCL12 readout).
- **Phase 2 repeat.** All 14 model x program delta R2 values have the same sign as Mendeley, and 12/14 have CIs excluding 0. As in Mendeley, the ridge gains for T:cytotoxic and CAF:CXCL12 have CIs spanning 0.
- **Phase 3 repeat.** All 6 chains pass the floor in >= 3 units. Same-direction rates against Mendeley are close to chance: hop B<-A 0.61, hop C<-B 0.59, indirect 0.45. This is what essentially null effects predict. The one single-hop association with consistent direction is CXCL12-CAF <- basal-like tumour (above).

## Claims register (`results/chain/phase7/claims_register.csv`)
| ID | Claim | Status |
|---|---|---|
| C1 | No T-cell program is measurable under the pre-declared section-wide background rule, so T-containing chains are not testable under it | not testable |
| C2 | With the type-aware rule (Amendment 1), 7 programs are measurable; classical vs basal-like agrees with the authors' subtypes (median AUC 0.858) | supported |
| C3 | Neighbourhood improves held-out prediction beyond own covariates for 6/7 receiver programs | supported |
| C4 | The gain comes from same-type neighbours and composition; other types' programs add at most 0.0005 | supported |
| C5 | T-cell cytotoxic state is predicted by neighbours | not supported (marker-label HGB only: exploratory) |
| C6 | A graph model adds value over ridge or HGB | not supported |
| C7 | Any A -> B -> C chain meets the pre-declared chain-evidence rule | not supported |
| C8 | Primary tumour -> CAF -> T chain (hop or indirect) | not supported |
| C9 | CAFs adjacent to basal-like tumour cells have lower CXCL12 (single hop; direction-consistent in GSE310352) | exploratory |
| C10 | CAF <-> myeloid C1Q_TAM adjacency association (marker labels only, both directions) | exploratory |
| C11 | The chain test detects planted 1.1x-2x chains | not supported |
| C12 | Single-hop tests are robust to zonation | not supported |
| C13 | Any ligand-receptor candidate is attached to a supported hop | not supported |
| C14 | Phase 2/3 directions are consistent in GSE310352 | exploratory (Phase 2 consistent; Phase 3 near chance) |

## Limitations
1. **The chain test has essentially no power** at realistic sparse effects (Phase 4), so the Phase 3 "not supported" results are uninformative about whether chains exist. They show only that this design cannot find them. Undirected adjacency also cannot identify direction: reversed chains are symmetric by construction.
2. **The measurement is weak.** 630/960 genes are at background section-wide. T-cell programs exist only under a logged post-hoc amendment, and exhaustion, iCAF, EMT and inflammatory-myeloid programs are not measurable. The mediator (B program score) is a noisy proxy for a B "state".
3. **The same-type neighbourhood gain may partly reflect segmentation spillover** between same-type neighbours and sub-200 um zonation. The spillover covariates cover only other-type lineage genes.
4. **The label block null is ill-posed** with a two-sided |statistic| test (null means far from 0). This was reported, not corrected after the fact.
5. **GSE310352 is not independent:** patient independence, um/px and the slide-5 grouping are UNVERIFIED. The program gene lists have UNVERIFIED literature sources (recalled, not fetched).

## Next steps (at most 5)
1. Redesign the exposure so it can detect sparse transmission, for example the maximum or the fraction of "high" B neighbours rather than the mean. Recalibrate with the same Phase 4 harness before touching real data again.
2. Add a centred label null, comparing |obs - null mean|, as a pre-declared revision, and recalibrate it in simulation.
3. Model segmentation spillover explicitly, for example a same-type contamination covariate or transcript-level reassignment, before interpreting the same-type neighbourhood gain.
4. Test the exploratory CXCL12-CAF / basal-like adjacency association with a pre-registered, single-hop analysis in an independent cohort with verified patient IDs. A wet-lab version would be co-culture of basal-like vs classical tumour cells with CAFs, measuring CXCL12. That is a hypothesis to test, not a finding.
5. Use a higher-plex or protein-validated T-cell panel, so that T-cell states (including exhaustion) are measurable above background without an amendment.

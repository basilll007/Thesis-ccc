# Chain-v1 progress log

Branch `chain-v1`, created 2026-11-15 from origin/thesis-v3 (8021098). Rules: never push to main, never force-push; seed 42; relative paths; add-one permutation p-values; the tumour is the unit of inference; expression, not mutation; associations and ranked hypotheses only. Each phase's plan is committed under docs/chain_plans/ before its results are computed. Files are registered in results/chain/MANIFEST.csv, figures go in figures/chain/<phase>/, and run logs go in logs/chain_*.log.

Environment: conda env `thesis-spatial-ccc` (D:/Anaconda/envs/thesis-spatial-ccc/python.exe; Python 3.11, numpy 1.26.4, scikit-learn 1.5.2, pandas 2.2.3). The base Anaconda Python is not used: its NumPy 2 is incompatible with the stack (see the v3 log). torch is not installed.

Interpretation logged at the start: the brief's "use author labels as primary annotation" is read as `author_mapped` (the v1 AUTHOR_MAP of the authors' majortypes onto our 10 classes plus NK, schwann_nerve and unassigned). "Marker labels" are the v2 `marker` labels.

## 2026-11-15 Phase 0: literature check (DONE)
- Wrote docs/chain_plans/related_work.md. Four methods were fetched through the Europe PMC REST API: MISTy (full text), NCEM (full text), CellNEST (full text and GitHub README) and HoloNet (abstract only).
- BLOCKED and not worked around:
  - pmc.ncbi.nlm.nih.gov returned a reCAPTCHA page for PMC9011939 and PMC10017508.
  - biorxiv.org returned HTTP 429 for the HoloNet and NEST preprints.
  - pubmed.ncbi.nlm.nih.gov returned a cookie wall.
  - nature.com and link.springer.com redirected to authentication.
- Methods in the brief's "and similar" that were not fetched (NicheCompass, stMLnet, COMMOT) are UNVERIFIED and not cited.
- System note: free RAM was 3.4 GB of 31 GB at the start, so the memory-safety rules (one section at a time, at most 6 workers) are applied conservatively.

## 2026-11-15 Phase 1: cell states (DONE, with logged Amendment 1)
- The plan (docs/chain_plans/phase1_cell_states.md) was committed in c93fb80 before computing. Script: src/chain/chain_phase1_cell_states.py. Outputs: results/chain/phase1/ (pre-declared rule) and results/chain/phase1b/ (Amendment 1). All numbers below are from those files.
- **Pre-declared section-wide background rule:**
  - 630 of 960 panel genes are background-level (phase1/phase1_summary.json).
  - Programs kept: tumour classical (5 genes), basal_like (5), EMT (6); CAF myCAF (7) and the CXCL12 readout; myeloid SPP1_TAM (5).
  - Dropped: iCAF (2 usable), T cytotoxic (0), T exhaustion (0), myeloid inflammatory (1), C1Q_TAM (3).
  - Under this rule no T-cell program is measurable, so every T-containing chain, including the primary, is "not testable".
- **Amendment 1** (committed in cd9af18 before it was computed; a logged deviation):
  - The section-wide rule removes rare-type genes by construction. The maximum negative-probe detection is 0.039-0.101 per section, and even CD3D is near background in 16/16 sections.
  - The type-aware rule compares detection within own-type cells against negative-probe detection in the same cells.
  - Kept: tumour classical (8), basal_like (10); CAF myCAF (8), CXCL12 (1); T cytotoxic (9); myeloid C1Q_TAM (6), SPP1_TAM (5).
  - Dropped: tumour EMT (4), iCAF (2), T exhaustion (4), myeloid inflammatory (3).
  - Phases 2-4 use Amendment 1. Chain calls carry the flag.
- **Validity check:** classical minus basal_like separates the authors' CLS and BSL malignant cells with a median AUC of 0.792 (range 0.626-0.857) under the original rule and 0.858 (0.645-0.936) under the amended rule. Both are above 0.5 in all 16 sections (phase1*/validity_auc.csv). The myCAF/iCAF check is not computable because iCAF was dropped.
- **Floors** (>= 100 cells per type, >= 300 A-B and B-C edges; author_mapped):
  - The primary chain tumour -> CAF -> T passes in 14/16 sections and 11/13 tumours. T2-b (94 T cells) and T6-b (67 T cells) fail.
  - Secondary chains chosen by the outcome-free rule: myeloid->CAF->tumour, tumour->CAF->myeloid, CAF->myeloid->tumour, CAF->tumour->myeloid, myeloid->tumour->CAF (phase1b/chain_selection.csv).
  - No secondary chain contains T cells.
- **Median nearest-neighbour distance:** 9.12 um across sections (range 7.93-10.39), from phase1b/phase1_summary.json. This is used to scale GSE310352 tiles in Phase 6.
- **Phase 4 section choice** (plan rule; fixed now, before outcomes). Among the 14 passing sections:
  - largest: T4-a (3,333 T cells);
  - median, taking the lower middle of an even count: T1-b (733);
  - smallest: U2-b (120).
- **Phase 2-4 plans** were updated before commit to the realised programs and chains. The Phase 4 plant uses CAF:myCAF instead of the draft's CAF:iCAF, which was dropped. They are committed together with the Phase 1 results, before any Phase 2-4 computation.

## 2026-11-15 Phase 2: neighbourhood prediction (DONE)
- The plan was committed in ab6ef88 before computing. Scripts: src/chain/chain_phase2_neighbourhood.py, and src/chain/chain_phase2b_gnn.py (SGC-style model; the pre-declared GNN rule triggered, see below). Results: results/chain/phase2/delta_r2_by_program.csv and feature_group_attribution.csv.
- Neighbours predict cell state beyond own covariates (tumour-bootstrap 95 % CI > 0) for 6 of 7 receiver programs. Author_mapped, no spillover covariates, ridge delta R2:
  - tumour classical 0.141 [0.105, 0.164]
  - tumour basal_like 0.057 [0.034, 0.074]
  - CAF myCAF 0.034 [0.021, 0.045]
  - CAF CXCL12 0.0085 [0.0027, 0.0147]
  - myeloid C1Q_TAM 0.030 [0.016, 0.045]
  - myeloid SPP1_TAM 0.056 [0.042, 0.065]
- **T cytotoxic: no gain.** Ridge 0.0005 [-0.0027, 0.0046]; HGB 0.0057 [-0.0016, 0.0087]. Under marker labels the HGB gain is 0.0129 [0.0046, 0.0268], so this result is annotation-sensitive.
- **Which features carry the gain** (ridge drop/add-group): almost entirely same-type neighbour programs (e.g. classical drop-group 0.117) and composition (classical 0.012; SPP1_TAM 0.015). Other-type neighbour programs contribute <= 0.0006 delta R2 for every receiver program. Their CIs include 0 or sit at about 0.0005, e.g. myeloid programs -> CAF CXCL12 0.0005 [0.0002, 0.0008].
- Adding the spillover covariates changes delta R2 little (e.g. classical 0.141 -> 0.123).
- **GNN rule triggered** (phase2/gnn_decision.json) for SPP1_TAM and classical: HGB exceeded ridge by >= 0.01 and the lower CI was > 0.01. The SGC-style model was run as declared.
- Operational note: Phase 4 tasks for T4-a failed within seconds on the first pass. A direct retry succeeded, and the likely cause is transient memory pressure from 5 concurrent T4-a loads (about 2.6 GB peak each) while Phase 2 was still finishing. This is logged as failure 1 for that task; the tasks are rerun through the per-task resume.
- **Phase 2b SGC model.** The first run (gnn_sgc_delta_r2.csv) gave an implausible delta R2 of 0.67 for SPP1_TAM. Cause: a leakage bug, where the second propagation step included i -> j -> i paths, so each cell's own score entered its features. Those files are kept unmodified and marked INVALID in the manifest. The fixed v2 excludes return-to-self paths (gnn_sgc_v2_delta_r2.csv):
  - SPP1_TAM 0.048 [0.031, 0.059], vs ridge 0.056 and HGB 0.068;
  - classical 0.135 [0.097, 0.164], vs 0.141 and 0.151;
  - classical with spillover covariates 0.116 [0.085, 0.150].
  So the graph model adds nothing over the simpler models. The main Phase 2 features were re-checked: the neighbour means and the leave-one-out tile mean exclude the cell itself.

## 2026-11-15/16 Phase 4: calibration by simulation (DONE; result: no trustworthy setting)
- The plan was committed in ab6ef88. Script: src/chain/chain_phase4_simulation.py. Sections: T4-a (3,333 T cells), T1-b (733) and U2-b (120). Planted triple: tumour:basal_like -> CAF:myCAF -> T:cytotoxic. Settings: 100 simulations each with 100 expression-null permutations (logged reduction from 200); label null on 20 simulations for 4 settings.
- **Failures:** the first pass failed for 11 T4-a tasks with `_ArrayMemoryError` (failure 1, memory pressure from 5 concurrent loads). They were rerun with 2+3 workers and all 33 tasks completed. No second failure.
- **Results** (results/chain/phase4/calibration_table.csv, calibration_pooled_3sections.csv, calibrated_range.json):
  - **Power of the full chain call is <= 0.04 in every section and fold, including 2x per hop.** It is 0.00 for every fold when the 3 sections are pooled.
  - The upstream hop B<-A is detected well: rate 0.83-1.00 at f >= 1.1 in T4-a, and >= 0.94 at f >= 1.25 in all sections.
  - The downstream hop C<-B is detected at about the null rate (T4-a 0.11-0.25; T1-b 0.12-0.20; U2-b 0.00-0.01).
  - Reason: only a minority of qualifying C cells are perturbed (median 129 of 3,088 T cells adjacent to CAFs in T4-a; 27 of 610 in T1-b; 46 of 106 in U2-b), and the neighbour-mean B exposure tracks the perturbed B neighbours only weakly.
  - **Zonation fools the single hops in the large section.** T4-a zonation at 2x: hop C<-B rate 1.00 and hop B<-A 0.97. The full chain call stays low (0.05) only because the indirect criterion fails. The label null also flags C<-B under zonation (0.95, 20 simulations).
  - Direct A -> C is never called a chain (0.00-0.01).
  - Under the null (f = 1) the indirect-only rejection rate is 0.29 in T1-b. The indirect p-value is therefore not well calibrated at about 50 analysed C cells.
- **Calibrated range: none.** Every Phase 3 chain call is outside the calibrated range (criterion c5 = False).

## 2026-11-16 Phase 3: chain test (DONE)
- The plan was committed in ab6ef88. Script: src/chain/chain_phase3_chain_test.py. Settings: 200 expression and 200 label block-null permutations per section (200 um tiles, random offsets); 16 sections, 13 tumours. Final files: results/chain/phase3/*_author_mapped.csv (pre-calibration versions kept with the suffix _precalibration).
- Number of tests after deduplication (forward): hop C<-B 22, hop B<-A 20, indirect 44. The negative-control families are listed in n_tests_author_mapped.json.
- **All 44 forward chain x program combinations: "not supported".** No forward or negative-control variant meets decision criteria 1-3.
- The smallest forward BH q under the expression null is 0.050 (hop B<-A) and 0.055 (hop C<-B), and 0.109 for indirect. Every hop with p_expr < 0.05 fails the label null.
- **Primary chain tumour -> CAF -> T:**
  - hop T <- CAF: beta 0.003 to 0.010, CIs span 0;
  - indirect: -0.008 to 0.008, CIs span 0.
  - The strongest primary-chain hop is CAF:CXCL12 <- tumour:basal_like, beta -0.028 [-0.039, -0.018], negative in 11/11 tumours, p_expr = 0.005 (the minimum), q = 0.050, p_label = 1.0. It is unchanged with spillover covariates (-0.027).
- **Label-null caveat:** the label block null has means far from 0 (for example -0.05 to -0.24 for hop betas), because relabelled cells carry another type's expression. The pre-declared two-sided test on |statistic| is therefore conservative or ill-posed for this null. This is reported as a limitation and the rule was not changed after the fact.
- **Amplification** (|beta_CB| / |beta_BA|): median 0.80, IQR 0.52-1.48, across 44 combinations of effects that are all near zero. No amplification pattern.

## 2026-11-16 Phase 5: candidate signals (DONE)
- The plan was committed in 7c3a6fc. Scripts: src/chain/chain_phase5_candidates.py and chain_phase5b_dedup.py.
- 0 of 88 hops meet the support rule (q_expr < 0.05 and p_label < 0.05).
- results/chain/ranked_hypotheses.csv has 516 rows. Each hop is repeated per program context, so results/chain/ranked_hypotheses_dedup.csv is provided with 270 unique hop x ligand-receptor rows. All rows are "hop not supported" and unranked (rank = NA), with calibration status "outside calibrated range".

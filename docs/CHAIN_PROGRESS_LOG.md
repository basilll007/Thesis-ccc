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

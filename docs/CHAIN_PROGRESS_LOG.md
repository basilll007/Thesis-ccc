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

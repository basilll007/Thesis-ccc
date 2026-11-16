# model_dev progress log

Branch `model-v1`, created 2026-11-16 from origin/main (cfb111c). Rules: never push to main, never force-push; push only `model-v1`; plans committed before each phase runs; every number from a generated file; expression not mutation; associations not causal claims.

## 2026-11-16 Phase 0: setup (in progress)
- Working tree F:/Thesis was clean and no process was using it (another agent works in the separate worktree F:/Thesis-ccc-translational), so the branch was created in place rather than in a new worktree.
- No conda env on the machine has a CUDA torch (`scgpt-zeroshot` has torch 2.3 CPU-only). Created env `thesis-model` with a CUDA 12.8 torch build (RTX 5060 Laptop, Blackwell). pip temp and cache redirected to F:\Thesis\model_cache because C: has ~0.1 GB free; nothing on C: was deleted.
- Panel check: Mendeley 960 Base genes are all present in GSE310352 (970); 10 GSE-only genes are dropped. No Mendeley gene lost.
- Type-aware gene flags (dev set): 539 genes usable for >= 1 typed class; 160 for all four. These 539 form the model gene set G.

## 2026-11-16 Phase 0: data build (DONE)
- Script `model_dev/src/data.py`, run in env `thesis-spatial-ccc` (no torch needed; numpy 1.26). Outputs in `model_dev/results/phase0/`; datasets in `model_cache/datasets/` (24 files, checksums in `artefacts.csv`).
- First attempt failed on a module-name clash (`model_dev/src/common.py` shadowed `src/pdac/common.py`): renamed to `mdcommon.py`. Second attempt was slow on the external slides because the panel set was rebuilt inside a loop: fixed and the whole build re-run from scratch (nothing kept from the earlier attempts).
- Program scores reproduce the chain cache exactly in all 16 Mendeley sections (max abs diff 0.0; `score_reproduction_check.csv`).
- Panel: 960/960 Mendeley genes in GSE310352; 10 GSE-only genes dropped; model gene set G = 539 genes (`gene_set.csv`, `panel_intersection.json`).
- External set first loaded at commit 58ce3c9 for the data card only (`external_first_load.json`); no model touches it until Phase 4.

## 2026-11-16 Phase 1: targets and splits
- First run of the declared greedy rule produced folds with a single tumour (2626 and 2762 alone), failing the >= 2 tumours floor. Plan Amendment 1 (snake deal of the first 10 tumours, then greedy) was written and committed before any model was fitted; the rejected check is kept as `fold_floor_check_v1_rejected.csv`.

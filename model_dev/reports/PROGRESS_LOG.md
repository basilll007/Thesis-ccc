# model_dev progress log

Branch `model-v1`, created 2026-11-16 from origin/main (cfb111c). Rules: never push to main, never force-push; push only `model-v1`; plans committed before each phase runs; every number from a generated file; expression not mutation; associations not causal claims.

## 2026-11-16 Phase 0: setup (in progress)
- Working tree F:/Thesis was clean and no process was using it (another agent works in the separate worktree F:/Thesis-ccc-translational), so the branch was created in place rather than in a new worktree.
- No conda env on the machine has a CUDA torch (`scgpt-zeroshot` has torch 2.3 CPU-only). Created env `thesis-model` with a CUDA 12.8 torch build (RTX 5060 Laptop, Blackwell). pip temp and cache redirected to F:\Thesis\model_cache because C: has ~0.1 GB free; nothing on C: was deleted.
- Panel check: Mendeley 960 Base genes are all present in GSE310352 (970); 10 GSE-only genes are dropped. No Mendeley gene lost.
- Type-aware gene flags (dev set): 539 genes usable for >= 1 typed class; 160 for all four. These 539 form the model gene set G.

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
- Amended rule run: folds 0-4 hold 2/2/2/3/4 tumours with 15,642-25,193 myeloid cells each; all pass the floor (`fold_floor_check.csv`). 11 unit tests pass (`tests/`: leakage guard incl. planted target, tile LOO, train-only scaler, 2-hop seed exclusion, fold rules).

## 2026-11-16 Phase 2: baselines (running)
- Reproduction check PASSED exactly: delta R2 myeloid C1Q_TAM 0.0299 vs chain 0.0299; SPP1_TAM 0.0561 vs 0.0561; tumour classical 0.1411 vs 0.1411; basal_like 0.0568 vs 0.0568 (|diff| 0.0000, tolerance 0.01; `results/phase2/reproduction_check.csv`).

## 2026-11-16 Phase 3: graph models (first attempt partly failed; resumed)
- Attempt 1 completed all 60 myeloid runs (2 architectures x 2 learning rates x 3 seeds x 5 folds; per-tumour SSE, curves and checkpoints saved), then failed with `_ArrayMemoryError` while assembling the tumour receiver matrix (246,249 x 2,695 float32) with Phase 2 running in parallel. Failure 1 for this step.
- Fix: Phase 1 plan Amendment 2 (tumour receiver capped at 8,000 cells per section, seeded; declared before any tumour result; myeloid uncapped). The capped tumour cache (102,546 cells) was built before Phase 2 reached the tumour receiver, so Phase 2 uses the same capped rows.
- `train.py` now saves results per receiver and resumes per receiver. The myeloid results/compute tables lost in the crash were reconstructed by `src/reconstruct_phase3_myeloid.py` from the saved curves (best val MSE and epochs, cross-checked against the log), the per-tumour SSE file and the run log (wall time, peak VRAM). No model was retrained or re-evaluated for this.

## 2026-11-16 Runs stopped by the host for low system memory (not a code failure)
- Claude Code stopped the Phase 2 baselines run and the Phase 3 relaunch at the same time because the system was critically low on memory while both ran in parallel (plus the earlier tumour cache build). Per the host's instruction they are NOT restarted automatically; waiting for the user.
- State at the stop: Phase 2 finished all 5 myeloid folds (per-tumour SSE saved, MLP curves saved) but the results table with inner-validation scores is written only at the end, so it was lost; the tumour receiver had started. Phase 3: myeloid complete (resumable); tumour seed 0 complete for all 20 runs (curves and checkpoints saved), seeds 1-2 not run; the tumour per-tumour file was not yet written.
- Code changes (no computation): Phase 2 now saves results per receiver and resumes per receiver; Phase 3 saves per-seed partial files and resumes per seed. Recommendation: run phases sequentially, not in parallel.

## 2026-11-16 Restart (user approved: sequential)
- Phase 2 (`baselines.py --skip-reproduction`; the reproduction check had already passed and is saved) then Phase 3 (`train.py`, resuming myeloid) are run one after the other in a single chain, never in parallel.

## 2026-11-16 Third stop and move to GPU (user instruction)
- The sequential Phase 2 run was also stopped by the host for low system memory during the tumour receiver (myeloid baselines had completed). Diagnosis: with none of our jobs running, Windows commit charge was 40.6 of 50.1 GB (other applications); our CPU run added several GB of host-RAM feature copies.
- User instruction: train in batches with checkpoints, on the GPU only (no CPU training), complete the work, update the HTML dashboard, commit to main.
- Phase 2 plan Amendment 1 (GPU execution; XGBoost-CUDA replaces sklearn HGB; ridge closed-form on GPU) written before any Phase 2 result was inspected. CPU myeloid results archived in `results/phase2/cpu_v1_myeloid/`. 13 unit tests pass (2 new GPU-vs-CPU equivalence tests).

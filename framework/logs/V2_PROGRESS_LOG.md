# framework-v2 progress log

Append-only. Each phase records free RAM before start, decisions taken without asking, and
DONE / BLOCKED status.

## Setup (2026-11-17)
- `git fetch`; branch `framework-v2` created from origin/framework-v1 (8a9d8ea) in a separate
  worktree `F:/Thesis-framework-v2`, because an OpenAI Codex app-server/exec-server process was
  running on this machine and could be using the `F:/Thesis` working tree.
- `rtk` is not on PATH in this shell; plain `git` used.
- Data is read from `F:/Thesis/data` (git-ignored, shared). New downloads go to
  `F:/Thesis/data/v2_cache/` (outside git).
- Python environments found: `F:/Thesis/model_cache/envs/thesis-model-v2` (lifelines 0.30.3,
  pytest, psutil, ccc_audit 0.1.0) for package, survival and tests; conda `pdac-bench`
  (squidpy 1.6.5, liana 1.5.1, spatialdm 0.3.1, commot) for Phase 5. R / Rscript not on PATH.
- Free RAM at first check: 7.35 GB (Win32 FreePhysicalMemory) / 7.5 GB available. Below the
  8 GB floor, so no phase was started; re-checked a few minutes later at 8.5 GB available.

## Phase 1
- Free RAM before start: 8.5 GB available (`\Memory\Available MBytes` = 8507).
- 1a overlap: all 90 GSE28735 arrays (45 tumours) are inside GSE62452 (same running sample number, mutual best
  match). The pre-declared raw r > 0.99 rule passes only 21/90 arrays (13 tumours) because of Partek re-normalisation;
  reported as-is. A post-hoc probe-centred correlation was added: matched r >= 0.983 vs non-matched <= 0.833.
  41 patients are in both Phase B survival analyses.
- 1c GSE57495: GPL15048 coverage 12/12. Primary stage-adjusted HR 1.347 [0.994, 1.823], one-sided p = 0.027
  -> CONFIRMED by the pre-declared rule (two-sided CI includes 1).
- 1d overlap: no GSM/title/CEL overlap. Cross-platform screen vs TCGA/CPTAC flagged nothing; max r 0.778 (PDAC19)
  is followed up post-hoc (excluding it: HR 1.31, one-sided p = 0.044).
- Literature: PMC returned reCAPTCHA; DOAJ 403; HAL bot-check; Springer auth redirect. None were worked around.
  Verified via the Frontiers publisher page (Abou Khouzam 2021) and the Europe PMC REST API (Ren 2023, Connor 2019).
- Phase 1 status: **DONE**.

## Phase 2
- Free RAM before start: 10182 MB available.
- Amendment 1 committed before Phase 2 results (synthetic unit-test development showed low power of block-null Q
  when exposure is widespread; zonation check now uses the within-tile receiver contrast plus the ring control).
- Reproduction run 1: 6/6 FAILED -> found the v1 `knn_graph` bug (dropped each cell's nearest neighbour). Fixed;
  run 2: 6/6 PASS. Run-1 outputs kept with suffix `_run1_knn_bug`.
- v2 demo outputs were regenerated after the kNN fix (the claim verdicts did not change; only the thesis-style
  C/E/M/Q tables use `knn_graph`).
- Tests: 29/29 pass. Phase 2 status: **DONE**.

## Phase 3
- Renamed `ccc_audit` -> `cellalibi` (import, distribution, CLI), v0.2.0. Installed editable into the
  thesis-model-v2 venv with `uv pip install -e --no-deps` (the venv has no pip).
- PyPI 404 for `cellalibi` and `cell-alibi`; GitHub repo/user search 0 hits; GitHub code search not run (no auth).
  See framework/docs/NAMING_LOG.md. Nothing published.
- After rename: 23 non-repro tests pass; 6/6 reproduction tests pass (pytest_reproduction_after_rename.txt);
  `cellalibi run` works on the crop config.
- Phase 3 status: **DONE**.

## Phase 4
- Free RAM before start: 9908 MB available.
- Section selection: T6-b / T2-b / U5-a. Smoke test (scratch dir) passed.
- Frozen v0.2.0 benchmark complete (4,200 audits, 0 errors): **FAIL**. False-SUPPORTED 0/300 in every artefact
  setting, but S1 2x sensitivity 0/300. Causes diagnosed: zero-trap power plant, 10x marker rule (spillover
  NOT_RUN), weak Q under saturated exposure. See framework/results/v2_phase4/PHASE4_SUMMARY.md.
- Post-hoc engine v0.2.1 + plan v2_phase4b.md committed (2720a9e) after T2-b/T6-b results, before U5-a results.
- Phase 4 (pre-registered) status: **DONE (result: FAIL)**.

## Phase 4b (post-hoc v0.2.1)
- Free RAM before start: 7587 MB. Tests before run: 29 passed, 6 deselected in 28.79s
- RAM reached 8111 MB at 00:05:56 after waiting (was 7589 MB); benchmark launched.
- Phase 4b complete (4,200 audits, 0 errors). C1 PASS (false-SUPPORTED max 0.04; S1 2x 0.913). C2 FAIL on pooled
  S1 2x sensitivity (0.673; 0.962 given power PASS; T6-b has 222 tumour cells). Both pooled 0.793 (FAIL by the
  pre-declared pooled rule). No setting fools v0.2.1. Phase 4b status: **DONE**.
- User request (2026-11-18 ~00:30): use the GPU wherever possible. Added cellalibi/gpu.py (torch CUDA permutation
  nulls, opt-in device="cuda"); CPU/GPU agreement tests 3/3 pass. Phase 4b left on CPU (already running).

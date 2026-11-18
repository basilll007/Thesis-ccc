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

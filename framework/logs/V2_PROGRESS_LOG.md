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

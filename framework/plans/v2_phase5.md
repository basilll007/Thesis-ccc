# v2 Phase 5 plan: existing tools on the same benchmark replicates

Committed before any tool is run on benchmark data.

## Question

**Would a user of an existing tool be misled?** For each replicate, record whether the tool reports the planted claim,
CAF_fibroblast -> myeloid via FN1 -> ITGB1, as significant. Then report:
- the detection rate under S1 (true contact signalling);
- the false-positive rate under the artefact settings S2 (abundance), S3 (zonation) and S4 (spillover), plus S0.

These tools answer different questions from cellalibi. They score or rank interactions; cellalibi audits a specific
claim against named artefacts. This comparison therefore measures whether a user would be misled, not overall tool
quality.

## Replicates

- These are the Phase 4 replicates, regenerated with the identical seeds and code (`cellalibi.simulate`, numpy in
  the thesis-model-v2 env).
- They are exported once to `F:/Thesis/data/v2_cache/phase5_replicates/` (outside git), so the tools in the
  `pdac-bench` env see exactly the counts cellalibi saw.
- Export check: re-auditing 3 exported replicates with cellalibi must reproduce the Phase 4 verdicts.
- Settings: S0, S1 (1.25 / 1.5 / 2), S2, S3 (50 / 200 / 400), S4 (0.1 / 0.2 / 0.3). S5 is not run for the tools.
- Subsample, declared because of runtime: replicates 0-19 of each setting for T6-b (smallest section) and
  replicates 0-9 for T2-b and U5-a. That is 11 x (20 + 10 + 10) = 440 datasets per tool, except COMMOT (below).
- Tool input: AnnData with X = raw = log1p(CP10K) of the simulated panel, `obs['cell_type']` = consensus_T labels,
  `obsm['spatial']` in um.

## Tools and settings

Statistical settings are the documented defaults. Where a parameter has no default, the value is declared here.
Execution-only settings (seed, number of jobs, threads) are fixed for reproducibility.

| Tool (version) | Call | Non-default values declared | "Reports the claim as significant" |
|---|---|---|---|
| Squidpy `gr.ligrec` (1.6.5) | `ligrec(adata, cluster_key='cell_type', interactions=[('FN1','ITGB1')])`; defaults: n_perms 1000, threshold 0.01, corr_method None | seed 0, n_jobs 1, numba_parallel False (the thesis found NaN p-values with the parallel kernel on > 32k cells) | p(FN1,ITGB1 \| CAF_fibroblast -> myeloid) < 0.05 |
| SpatialDM (0.3.1) | `weight_matrix(l=25, single_cell=True)` (cutoff default 0.1); `extract_lr('human', datahost='package')` subset to the FN1 -> ITGB1 entry; `spatialdm_global` defaults (z-score, n_perm 1000); `sig_pairs` defaults (z-score, FDR, threshold 0.1) | l = 25 um (no default; the thesis value) | `selected` = True by sig_pairs defaults. Global, per pair, not cell-type specific. |
| COMMOT (0.0.3) | `spatial_communication(df_ligrec=[FN1, ITGB1], dis_thr=50, heteromeric=False)` then `cluster_communication(clustering='cell_type')`; defaults: n_permutations 100, random_seed 1 | dis_thr = 50 um (default None means all pairs: a dense n x n problem) | cluster p(CAF_fibroblast -> myeloid) < 0.05 |
| LIANA+ `mt.bivariate` (liana 1.5.1) | `ut.spatial_neighbors(bandwidth=25)` (cutoff default 0.1); `mt.bivariate(resource=[FN1, ITGB1], local_name='cosine', global_name='morans', n_perms=100)` | bandwidth 25 um; global_name 'morans' and n_perms 100 as in the LIANA+ bivariate tutorial (both are None by default) | global Moran's p < 0.05. Global, not cell-type specific. |
| MISTy, Statial/SpatioMark | – | – | **BLOCKED: R is not installed** (no `R`/`Rscript` on PATH, no `C:/Program Files/R`) |

COMMOT memory:
- COMMOT builds dense n x n cost matrices, which the thesis also noted.
- It runs only on T6-b (12,473 cells), 1 worker at a time. T2-b and U5-a are logged as BLOCKED (memory).
- If T6-b exceeds the 8 GB headroom, COMMOT is BLOCKED entirely.

cellalibi: SUPPORTED counts as "reports the claim". Its values come from the Phase 4 chunks for the same
(section, setting, replicate).

## Output

`framework/results/v2_phase5/`:
- one resumable chunk per (section, setting, tool);
- `tool_detection_rates.csv` (rate with Wilson 95 % CI per tool x setting);
- `TOOL_COMPARISON.md`;
- tool versions and runtime per dataset.

No pass/fail criterion is set for the other tools. Rates are reported as observed.

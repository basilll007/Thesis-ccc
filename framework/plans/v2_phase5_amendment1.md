# v2 Phase 5 plan, amendment 1 (committed before the full tool run)

Date 2026-11-18. Based on a one-replicate smoke test (T6-b, S1 2x, replicate 0), which only checked that each tool
runs.

1. **Squidpy threshold workaround.**
   - Observed: with the documented default `threshold=0.01`, squidpy 1.6.5 (pandas 2.2.3, pdac-bench env) returned
     NaN p-values for 167 of 169 cluster pairs, including CAF_fibroblast -> myeloid, whose mean was 2.30.
   - Yet FN1 is detected in 53 % of CAFs and ITGB1 in 34 % of myeloid cells, far above 1 %. The same NaN appears with
     dense/sparse X, with raw/X, and with the parallel/serial kernel. The cause inside squidpy was not identified.
   - Workaround: run `threshold=0` and apply the default 1 % detection rule ourselves to the claim pair. If either gene
     is detected in < 1 % of its cluster, p = NaN (not significant). For this pair, this is what a working default
     would compute.
   - Logged per replicate in the `note` column.
2. **Replicate subsample (runtime).** Smoke timings on T6-b: squidpy 12 s, SpatialDM 7 s, LIANA+ 2 s, COMMOT 57 s per
   dataset.
   - T6-b: 20 replicates per setting for squidpy, SpatialDM and LIANA+.
   - T2-b and U5-a: **5** replicates per setting (was 10), because those sections are 3.5x and 6x larger.
   - COMMOT: T6-b only, **5** replicates per setting (55 runs).
   - Total datasets per tool: 11 x (20 + 5 + 5) = 330 for squidpy / SpatialDM / LIANA+, and 55 for COMMOT.
   - Exported replicates beyond these are not used by the tools.
3. **SpatialDM pair entry.**
   - After subunit filtering, 5 CellChatDB entries reduce to FN1 -> ITGB1. The first is used:
     `FN1_ITGAV_ITGB1`.
   - The receptor complex is then represented only by ITGB1 (the other subunits are not on the simulated panel).
     All 5 entries are identical on this panel.

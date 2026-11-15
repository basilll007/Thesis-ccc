# v3 Phase 1 plan — contact vs zonation (spatially constrained nulls)

Written 2026-11-15 BEFORE any Phase 1 result was computed. Branch `thesis-v3`.

## Problem
The v1/v2 null shuffles cell labels across the whole section. Under that null a positive Q = E/M
(edge product over matched cell-type product) can arise from tissue-scale co-zonation of the sender and
receiver types (CAFs and T cells occupying the same region) and does not require cell–cell contact.

## Data and scoring (unchanged from v2)
- Mendeley PDAC-SMI (CosMx), 16 sections; loader `src/pdac/pdac_loaders.py`; QC min_counts=10, min_cells=3;
  normalize_total(1e4) + log1p; expression from `adata.raw`.
- Graph: 6-NN directed (`phase3_pdac_audit.knn_edges`), identical to v1/v2.
- Scores C, E, M, Q from `src/pdac/ccc_audit.py`; 500 permutations; seed 42; add-one p = (k+1)/(n+1) everywhere.
- Labels: `consensus_T` (primary; post-hoc, best-validated T set), `marker` and `author_mapped` (secondary).
  All are read from `results/pdac/phaseB/sections/<section>_labels.csv.gz` / `_consensus_labels.csv.gz` (v2 files, not recomputed).
- Catalogue: all 31 testable pairs are computed for every sender x receiver; reporting is restricted to
  (a) the headline CXCL12->CXCR4, CAF_fibroblast -> T_cell and (b) the 7 axes whose
  `background_level_gene_diagnostic` is False in `results/pdac/phase3/background_flag_per_axis.csv`
  (TGFB1-TGFBR1, TGFB1-TGFBR2, CXCL12-CXCR4, JAG1-NOTCH1, MIF-CD74, SPP1-CD44, FN1-ITGB1).

## Nulls (all applied to the same observed labels, graph and expression)
1. `global` — reference: the v2 null (`spatial_null.generate_label_permutations`, seed 42). Reproduces Phase B3.
2. `block50`, `block100`, `block200` — block permutation: a square grid of side s um (s = 50, 100, 200), anchored at
   the section's minimum x/y; labels are permuted independently within each tile (numpy default_rng(42)).
   Tiles with one cell keep their label. Global type counts are preserved exactly.
3. `dist` — distance-stratified permutation. Boundary = the tumour compartment of the *same annotation*
   (class `epithelial_ductal_tumor`). Signed distance d_i: for a non-tumour cell, +distance to the nearest tumour cell;
   for a tumour cell, -distance to the nearest non-tumour cell. Bands (um): [0,25), [25,50), [50,100), [100,200),
   [200,400), [400,inf) on each side (12 bands). Labels are permuted within bands. Consequence: tumour labels are
   never exchanged with non-tumour labels, so the tumour/stroma partition and each stromal type's radial profile
   relative to tumour are preserved; only the arrangement of non-tumour types at the same distance is randomised.
   Skip rule: if a section has < 50 tumour or < 50 non-tumour cells under that annotation, `dist` is skipped and logged.

Rationale of ordering: global -> block200 -> block100 -> block50 is increasingly constrained
(smaller tiles preserve finer zonation). `dist` constrains a different axis (radial position). A contact-specific
signal should survive all of them; a zonation artefact should shrink as tiles shrink.

## Per-section outputs (results/pdac/v3_phase1/)
- `sections/<section>__<annotation>__<null>.csv.gz`: all interactions, C/E/M/Q obs, null mean/sd, z, add-one p,
  n_qualifying_edges, n_perms_zero_edges, n_null_valid. Written immediately after each run (resume on restart).
- `headline_by_null.csv`: headline row per section x annotation x null.
- `nonbg_catalogue_by_null.csv`: the 7 non-background axes x all sender/receiver per section x annotation x null,
  with BH within section x annotation x null.
- `nonbg_summary_by_null.csv`: per annotation x null: rows evaluated, Q rows with p<0.05, with BH q<0.05;
  Spearman(z_Q global, z_Q null) over shared rows; same for E.

## Tumour-level outputs
- Section -> tumour mapping: donor id with `_2nd`/`_B` removed (identical to Phase D2; 16 sections -> 13 tumours).
- `tumour_level_headline_by_null.csv`: per annotation x null x score (E, M, Q) x tumour: Fisher-combined section add-one p,
  BH across tumours, tumour-mean z, tumour-mean observed value, n sections.
- `tumour_level_summary_by_null.csv`: per annotation x null x score: n tumours, n significant (uncorrected, BH),
  n tumours with mean z>0, median tumour-mean z, one-sided Wilcoxon p of tumour-mean z>0 (exploratory),
  n tumours with mean observed Q>1.
- Key question: does the median tumour-mean z_Q (and the count of z>0 tumours) move toward zero from global to block50?

## Figure (one per phase): figures/pdac/v3_phase1/phase1_Qz_vs_null.{png,pdf}
Panel A: headline z_Q per section (thin lines) and tumour-mean (points) against null (global, block200, block100, block50, dist),
consensus_T; Panel B: same for marker and author_mapped tumour medians; Panel C: non-background catalogue:
fraction of rows with Q add-one p<0.05 per null.

## Interpretation rules (declared now)
- Report the trend either way. "Shrinks toward zero" = median tumour-mean z_Q under block50 is closer to 0 than under
  global AND the count of tumours with z>0 does not increase. Otherwise "does not shrink".
- Audit/diagnostic framing only; no biological or clinical claims.
- Every reported number must come from a file registered in `results/pdac/MANIFEST.csv`.

## Memory
One section x annotation x null at a time; results written after each; `mem_status()` logged per section.

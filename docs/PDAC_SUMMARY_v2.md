# PDAC Extension v2 — Summary (2026-10-25)

## Does a proximity-specific signal exist at tumour level, or is the edge score mainly mirroring cell-type products?

**The edge score mainly mirrors cell-type products.** A proximity-specific excess for CXCL12→CXCR4 CAF→T exists only as a weak, directionally consistent trend across tumours. **No single tumour or section shows it after multiple-testing correction.** This is an audit and diagnostic statement, not a biological or clinical claim.

- **Mirroring.**
  - Across the full catalogue (16 sections), the edge z and the z of the matched cell-type product M = L̄s·R̄r correlate at Spearman 0.66–0.68. On the 7 axes without background flags it is 0.80–0.81.
  - 77% of edge-significant section-rows are also M-significant.
  - Only about 10% have a ratio Q = E/M that is significant after BH (marker 0.098, author 0.104, consensus 0.096).
  - For the CAF→T headline, the edge z is below the M z in 16/16 sections (consensus T labels).
- **Proximity-specific Q, CAF→T** (with the best-validated T labels, `consensus_T`):
  - Sections: 2/16 significant uncorrected, **0/16 after BH within section**.
  - Tumours: 2/13 significant (Fisher-combined p), **0/13 after BH across tumours** (smallest q 0.285).
- **Consistency across tumours** (exploratory, not pre-registered):
  - Consensus labels: tumour-mean z_Q > 0 in 10/13 tumours (one-sided Wilcoxon p = 0.016); observed Q > 1 in 11/13.
  - Marker labels: 12/13 (p = 0.0085). Author labels: 12/13 (p = 0.0004).
  - The label-permutation null breaks all spatial structure. A positive Q may therefore reflect tissue-scale zonation rather than contact (see Limitations).
- **Independent dataset: not tested.** GSE310352 is BLOCKED.

Branch `pdac-v2` (from `main` @ df3867c). Every number is in `results/pdac/phaseD/v2_key_numbers.csv` or `results/pdac/phaseD/tumour_level_headline_summary.csv` (value plus source), or in the file cited in brackets. The full record, including failures and corrections, is in `docs/PDAC_PROGRESS_LOG.md` (section "PDAC v2"). v1: `docs/PDAC_SUMMARY.md`.

## Results

### A. Housekeeping

- **Tag**: `aaai27-submitted` → 3cde438. It was not on the remote, so only the tag was pushed.
- **Spot-check** [`results/pdac/phaseA/spotcheck.csv`]: the five v1 headline numbers were recomputed from row-level files. All 5 match:
  - breast gate 3/9
  - CAF→T edge 16/16 and 16/16
  - Q 4/16 and 2/16; 0/16 after BH
  - 24/31 axes background-flagged
  - Spearman(edge z, ligrec) on the crop 0.022
- **p-value conventions** [`results/pdac/phaseA/`; note in `docs/RESULTS_SETTINGS_NOTE.md`]: legacy is mean(null ≥ obs); add-one is (k+1)/(n+1).
  - **No uncorrected verdict flips when all 500 nulls are valid.** Breast: 0 flips (292 rows). PDAC: 9 of 532,952 rows flip, all Q rows with fewer than 500 valid nulls.
  - Legacy p = 0 makes BH more liberal: 3 breast and 4,530 PDAC BH verdicts flip.
  - CAF→T edge after BH: 16/16 → 13/16 (marker) and 16/16 → 14/16 (author) under add-one.
  - The breast gate is 3/9 under both conventions.

### B. T-cell label purity

**Pre-registered protein rule.** [`docs/PDAC_V2_PROTEIN_RULE.md`, committed before any labels or CCC results]
- Rule: per-section GMM or MAD thresholds on CD3, CD45 and PanCK; T = CD3+ CD45+ PanCK−.
- **It failed validation** [`results/pdac/phaseB/validation_precision_recall.csv`, `protein_labels_summary.json`]:
  - T precision / recall vs author T: 0.022 / 0.036 in U5-a; 0.077 / 0.049 pooled.
  - v1 marker scoring scores 0.276 / 0.297 and 0.353 / 0.318.
- Cause: CD3 is bimodal (Ashman D ≥ 2) in only 3/16 sections and CD45 in 2/16. PanCK works (15/16).

**Post-hoc `consensus_T`.** Defined before any v2 CCC run.
- T = author T **and** marker T; 3,676 cells pooled, 261 in U5-a.
- Protein support is weak for every T set. Median AUC vs non-T cells, CD3 / CD45:
  - consensus T: 0.598 / 0.769
  - author-only T: 0.580 / 0.739
  - marker-only T: 0.588 / 0.684
- So "validated" here means two RNA annotations agree; it does not mean protein-confirmed.

**CCC re-run.** C/E/M/Q; 500 permutations, seed 42; 6-NN; 31 pairs; 16 sections. Sources: `results/pdac/phaseB/headline_CAF_T_summary_by_annotation.csv`, `edge_vs_product_mirroring_by_annotation.csv`, `annotation_pairwise_comparison.csv`.

| CXCL12→CXCR4 CAF→T, sections | evaluable | E sig (BH) | Q sig (BH) | median z_E | median z_Q | median Q |
|---|---|---|---|---|---|---|
| marker (v1) | 16 | 16 (13) | 4 (0) | 20.87 | 0.75 | 1.347 |
| author (v1) | 16 | 16 (14) | 2 (0) | 21.63 | 0.74 | 1.346 |
| protein_gmm (failed validation) | 14 | 10 (6) | 4 (0) | 5.00 | 1.17 | 1.623 |
| protein_mad (failed validation) | 13 | 9 (4) | 5 (0) | 4.69 | 1.56 | 1.799 |
| **consensus_T (post-hoc)** | 16 | 15 (11) | 2 (**0**) | 14.68 | 0.36 | 1.320 |

- **U5-a, consensus_T**: z_C / z_M / z_E / z_Q = 21.6 / 24.6 / 17.5 / 1.94; Q = 1.96, p = 0.038, BH q = 0.314 (448 edges).
- **Annotation agreement** (Spearman of z on shared rows, pooled):
  - author vs consensus: E 0.970, Q 0.951
  - marker vs author: E 0.655, Q 0.394
  - author vs protein_gmm: E 0.348, Q 0.188
  - Q is the most annotation-sensitive score.

### Tumour level (Phase D2)

Source: `results/pdac/phaseD/tumour_level_headline_CAF_T.csv`, `tumour_level_headline_summary.csv`; figure `figures/pdac/phaseD/phaseD2_tumour_level_headline`.

- The 16 sections map to 13 tumours (donor id without the suffixes `_2nd` / `_B`).
- Tumour p = Fisher combination of the section add-one p; BH across tumours.
- Tumour mean z; one-sided Wilcoxon of tumour-mean z > 0.

| CAF→T, 13 tumours | E: tumours sig BH | M: tumours sig BH | Q: sig uncorr. (BH) | Q: z > 0 | Q: Wilcoxon p | Q: mean obs > 1 |
|---|---|---|---|---|---|---|
| marker | 13 | 13 | 3 (0) | 12 | 0.0085 | 12 |
| author | 13 | 13 | 2 (0) | 12 | 0.0004 | 13 |
| **consensus_T** | 13 | 13 | 2 (**0**) | 10 | 0.016 | 11 |
| protein_gmm (11 tumours evaluable) | 9 | 8 | 4 (0) | 10 | 0.0068 | 10 |
| protein_mad (10 tumours evaluable) | 7 | 5 | 6 (0) | 8 | 0.080 | 8 |

Synthesis figure: `figures/pdac/phaseD/phaseD_v2_synthesis.{png,pdf}`.

## Limitations

- **No ground-truth T cells.** "Validated" = agreement of two RNA-based annotations, with weak protein support. `consensus_T` is post-hoc: it was defined after the pre-registered rule failed, though before any CCC was run on it.
- **The tumour-level trend is exploratory.**
  - The Wilcoxon test was added after the section-level results were known, and it is run under five annotations.
  - The tumour grouping (`_2nd`/`_B` = same tumour) is an assumption consistent with the stated 13 tumours.
  - Fisher combination of discrete add-one p is conservative.
- **Label-permutation null.** It keeps cell positions and expression but shuffles labels, so it destroys spatial structure. A positive Q can come from tissue-scale co-zonation (CAF CXCL12 and T-cell CXCR4 both higher in the same regions) as well as from contact. The audit cannot separate the two.
- **Q is unstable when M is small.** The minimum add-one p is 0.002 with 500 permutations. Non-significance after BH is not evidence of absence.
- **One dataset.** All of it is a single CosMx dataset (16 sections, 13 tumours). No independent replication (BLOCKED).
- **Background-level genes.** 24/31 PDAC axes (all checkpoint axes) involve near-background genes (v1). CXCL12–CXCR4 is not flagged.

## BLOCKED items

1. **GSE310352 replication (Phase C).**
   - Data licence not confirmed: the GEO terms page sat behind reCAPTCHA and the publisher page returned 403. Neither was worked around.
   - The article (Lyubetskaya et al., Cell Rep 2026, doi 10.1016/j.celrep.2025.116827) is CC BY-NC-ND per Europe PMC; that is the article's licence, not a data licence.
   - NCBI's general policy says NCBI places no restrictions but submitters may claim rights.
   - Author cell labels exist only in `seurat.RDS`, and there is no R on this host.
   - Ready: `src/pdac/phaseC_gse310352_audit.py` (gated by `--licence-confirmed`; self-tested on synthetic data only).
   - Needed from you: paste the paper's "Data and code availability" section or the GEO disclaimer, or confirm NCBI's policy is sufficient. Also an exported per-slide `meta.data` CSV, or R + SeuratObject. Details and an R snippet are in the progress log.
2. **Push of `pdac-v2` after 161418a.** A TLS read error during the large upload (`SEC_E_MESSAGE_ALTERED`) happened twice; an OpenSSL-backend attempt also failed. TLS was not disabled. Run `git push origin pdac-v2` from a stable network.
3. **Carried over from v1:**
   - full-catalogue radius sweep (stopped by the host for low memory)
   - full-section Squidpy ligrec p-values (all NaN in 2 attempts)

## Next steps (at most 5)

1. **Unblock and run GSE310352.** Before running, pre-register the tumour-level directional test (Wilcoxon of slide-level z_Q > 0, consensus-style labels) as the replication criterion.
2. **Add a spatially structured null**, for example permutations restricted within FOV or local windows, or toroidal shifts of one cell type. This separates contact-scale Q from tissue-scale co-zonation, the main alternative explanation for the tumour-level trend.
3. **Raise Q's resolution on the headline only.** Use 5,000 permutations for CAF→T per section, so the add-one floor (0.002) and the BH burden are not what limits detection.
4. **Improve T-cell purity beyond RNA agreement.** For example, flag cells whose transcripts look like neighbour contamination (CAF-marker counts in T calls), or use transcript-level re-segmentation where raw files allow. The CosMx protein channels here cannot do it.
5. **Make the background flag part of the audit** (v1 next step 1), and re-run the two v1 BLOCKED items when there is memory headroom (chunked radius sweep; minimal ligrec reproducer).

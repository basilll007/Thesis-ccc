# PDAC Extension v2 — Summary (2026-10-25)

Branch `pdac-v2` (from `main` @ df3867c). This is an audit and diagnostic summary. Nothing in it supports a biological or clinical claim.

Every number below is in `results/pdac/phaseD/v2_key_numbers.csv` (value plus source file) or in the file cited in brackets. The full record, including failures and corrections, is in `docs/PDAC_PROGRESS_LOG.md` (section "PDAC v2"). v1 is summarised in `docs/PDAC_SUMMARY.md`.

## Key question

> Does the CAF→T proximity-specific signal (Q = E/M for CXCL12→CXCR4) survive validated T-cell labels and an independent dataset, or is the edge score mostly mirroring cell-type products?

**Answer:**

- **Validated labels: it does not survive.**
  - With the best-validated T-cell class available, Q is significant in 2/16 sections uncorrected and in **0/16 after BH within section**. The best class is `consensus_T`: author T and marker T agree. It is post-hoc, see Phase B.
  - The marker labels give 4/16 and 0/16; the author labels give 2/16 and 0/16.
  - The edge score E stays significant in 15/16 sections (11/16 after BH). In every section its z is below the z of the matched cell-type product M.
- **Independent dataset: not tested.** Phase C is BLOCKED; see below.
- **Full catalogue: the edge score mostly mirrors the cell-type product.**
  - Spearman(z_E, z_M) is 0.66–0.68; 0.80–0.81 on the 7 axes without a background flag.
  - 77% of edge-significant section-rows are also M-significant.
  - Only about 10% are Q-significant after BH (marker 0.098, author 0.104, consensus 0.096).
- **Caveat: non-significance of Q is not evidence of absence.** The observed Q is above 1 in 14/16 consensus sections (median 1.32), but it does not exceed the label-permutation null at 500 permutations after correction.

## Checklist

| Phase | Status | Result |
|---|---|---|
| A. Housekeeping | **DONE** | Tag `aaai27-submitted` → 3cde438 (pushed; tag only). 5/5 headline numbers re-derived from row-level files; no summary correction needed. p-value conventions compared side by side. |
| B. T-cell label purity | **DONE** (negative result) | Pre-registered protein rule **failed validation**. Post-hoc `consensus_T` used for the key question. CCC re-run on 16 sections × 3 new annotations. |
| C. Replication on GSE310352 | **BLOCKED** | Licence not confirmed (GEO terms page behind reCAPTCHA; publisher page 403). Author labels exist only in `seurat.RDS` (no R on this host). A gated, self-tested script is ready. |
| D. Synthesis | **DONE** | This file, `figures/pdac/phaseD/phaseD_v2_synthesis.{png,pdf}`, README links. |

## Phase A — housekeeping

- **Tag**: `aaai27-submitted` at 3cde438. It did not exist on the remote (`git ls-remote --tags` was empty), so only the tag was pushed.
- **Spot-check** [`results/pdac/phaseA/spotcheck.csv`]. Each number was recomputed from row-level tables, not copied from summary files. All 5 match:
  - breast gate 3/9
  - CAF→T edge-significant sections: 16/16 (marker) and 16/16 (author)
  - Q-significant sections: 4/16 (marker) and 2/16 (author); 0/16 and 0/16 after BH
  - background-flagged axes: 24/31
  - Spearman(edge z, ligrec mean) on the crop: 0.022
- **p-value conventions** [`results/pdac/phaseA/pvalue_conventions_*.csv`; note in `docs/RESULTS_SETTINGS_NOTE.md`]:
  - With all 500 nulls valid, legacy p = mean(null ≥ obs) and add-one p = (k+1)/(n+1) give **identical uncorrected verdicts**. The threshold #{null ≥ obs} ≤ 24 is the same for both.
  - Breast: 0 uncorrected flips (292 rows, 289 with a finite p). The Phase 1 gate is 3/9 under both conventions.
  - PDAC: 9 uncorrected flips out of 532,952 rows. All are ratio rows with fewer than 500 valid nulls.
  - BH verdicts differ. Legacy p can be 0, which makes BH more liberal: 3 breast flips and 4,530 PDAC flips, all significant only under legacy.
  - Headline CAF→T after BH within section: edge 16/16 (legacy) vs 13/16 (add-one) for marker labels, and 16/16 vs 14/16 for author labels. Q 0/16 vs 0/16 (marker) and 1/16 vs 0/16 (author).
  - Breast `extra` family: CD86→CTLA4 fibroblast→immune (edge) is BH-significant only under legacy (q 0.032 vs 0.064).
  - The recomputation from saved nulls matches the stored columns to ≤ 1.1e-16.

## Phase B — T-cell label purity

**Rule.** `docs/PDAC_V2_PROTEIN_RULE.md` was committed (5b6089e) before any label or CCC result existed.
- Thresholds per section on log1p(`Mean.CD3`/`Mean.CD45`/`Mean.PanCK`): a two-component GMM, accepted when Ashman D ≥ 2; otherwise median + 3·MAD.
- Gates: T = CD3+ CD45+ PanCK−. Other compartments are subtyped by RNA marker argmax.
- Sensitivity variant `protein_mad`: MAD thresholds for every channel.

**Validation against author labels** [`results/pdac/phaseB/validation_precision_recall.csv`, `protein_labels_summary.json`, `confusion_*`; figure `figures/pdac/phaseB/phaseB1_protein_gating_validation`]:

| T_cell vs author T (CD8/CD4/Treg) | U5-a precision / recall | 16 sections pooled |
|---|---|---|
| protein_gmm (pre-registered) | **0.022 / 0.036** (32 of 1,467) | 0.077 / 0.049 (563 of 7,273) |
| protein_mad | 0 T cells called | 0.031 / 0.005 |
| v1 marker | 0.276 / 0.297 (261 of 945) | 0.353 / 0.318 (3,676 of 10,427) |

- **Cause.** In this dataset the CD3 and CD45 channels barely separate the author T cells:
  - The GMM was accepted for CD3 in 3/16 sections, for CD45 in 2/16 and for PanCK in 15/16.
  - U5-a CD3 Ashman D = 0.43.
  - Per-section AUC for author T vs other cells: CD3 0.45–0.92, CD45 0.50–0.89 (computed inline; recomputable from the per-cell tables).
- **PanCK works.** It separates author malignant cells (AUC 0.74–0.99).
- **Other classes are also worse.** Pooled per-class precision under `protein_gmm` is below marker for every class with predictions, because the CD45 gate is weak.
- **Verdict: the pre-registered protein rule does not give validated T-cell labels.**

**Post-hoc `consensus_T`.** Documented and committed (161418a) before any v2 CCC run.
- `T_cell` = author T **and** v1 marker T. Author T without marker support is kept as a separate type. All other cells keep their author labels.
- Size: 3,676 cells pooled; 261 in U5-a.
- Independent protein check [`consensus_T_protein_triangulation_median.csv`; figure `phaseB2_consensus_T_protein_triangulation`]. Median AUC vs non-T cells, CD3 / CD45:
  - consensus T: 0.598 / 0.769
  - author-only T: 0.580 / 0.739
  - marker-only T: 0.588 / 0.684
- Protein support is highest for the consensus set but weak for every set. "Validated" therefore means two RNA-based annotations agree, **not** protein-confirmed.

**CCC re-run.** C/E/M/Q; 500 permutations, seed 42; 6-NN; 31 pairs × all sender-receiver combinations; 16 sections; add-one p; BH within section × annotation. Sources: `results/pdac/phaseB/v2_ccc_6NN_all_rows.csv.gz`, `headline_CAF_T_summary_by_annotation.csv`, `headline_CAF_T_by_annotation_section.csv`; figure `phaseB4_headline_and_annotation_agreement`.

| CXCL12→CXCR4 CAF→T_cell | evaluable | E sig (BH) | **Q sig (BH)** | median z_E | median z_Q | median Q |
|---|---|---|---|---|---|---|
| marker (v1) | 16 | 16 (13) | 4 (**0**) | 20.87 | 0.75 | 1.347 |
| author (v1) | 16 | 16 (14) | 2 (**0**) | 21.63 | 0.74 | 1.346 |
| protein_gmm (failed validation) | 14 | 10 (6) | 4 (**0**) | 5.00 | 1.17 | 1.623 |
| protein_mad (failed validation) | 13 | 9 (4) | 5 (**0**) | 4.69 | 1.56 | 1.799 |
| **consensus_T (post-hoc)** | 16 | 15 (11) | 2 (**0**) | 14.68 | 0.36 | 1.320 |

- **U5-a under consensus_T**:
  - z_C / z_M / z_E / z_Q = 21.6 / 24.6 / 17.5 / 1.94
  - Q = 1.96, p = 0.038, BH q = 0.314; 448 edges
  - The smallest consensus BH q over all sections is 0.314.
  - z_M > z_E in 16/16 consensus sections. This is a descriptive comparison: z values of different scores have different null spreads.
- **Annotation agreement** (Spearman of z on shared rows, 16 sections) [`annotation_pairwise_comparison.csv`]:
  - author vs consensus: E 0.970, Q 0.951
  - marker vs author: E 0.655, Q 0.394
  - author vs protein_gmm: E 0.348, Q 0.188
  - Q is the most annotation-sensitive score, as in v1.

## Phase C — independent replication (BLOCKED)

**Confirmed:**
- File format: per-slide counts CSV (`fov`, `cell_ID`, 970 genes + 19 NegPrb).
- Metadata CSV: centroids in px, plus `Mean./Max.` for PanCK/CD45/CD3/MembraneStain/DAPI.
- Sizes: 16 small CSVs, 124,527,604 B in total.
- No cell-type labels in the CSVs.
- Article: Lyubetskaya et al., Cell Reports 2026, doi 10.1016/j.celrep.2025.116827. Its Europe PMC record lists `license: "cc by-nc-nd"` (the article's licence).
- NCBI's general Molecular Data Usage policy (fetched) says NCBI places no restrictions but submitters may claim rights. That is not a dataset licence.

**Not confirmed:**
- The data licence. The GEO terms page sat behind reCAPTCHA, and the publisher page returned 403. Neither was worked around.
- The authors' labels are only in `seurat.RDS`, and this host has no R. Phase B shows the protein channels cannot substitute for them.

**Ready:**
- `src/pdac/phaseC_gse310352_audit.py` (gated by `--licence-confirmed`; exits with BLOCKED otherwise).
- `pdac_loaders.load_geo_gse310352`.
- Self-tested end to end on synthetic data only. All synthetic outputs were discarded; no GSE310352 numbers exist in the repository.

**What is needed from you**: listed in the Phase C entry of `docs/PDAC_PROGRESS_LOG.md`. In short:
1. Licence: paste the paper's "Data and code availability" text or the GEO disclaimer text, or confirm that NCBI's policy is sufficient.
2. Labels: R + SeuratObject, or an exported `meta.data` CSV per slide (an R snippet is in the log).

## Other failures, corrections and process notes

- **Correction (Phase A numbers).** The settings note and the log first gave the PDAC BH-flip count as 4,533 and the PDAC row count as 533,244. Both are breast + PDAC totals. Corrected to 4,530 and 532,952, with the corrections marked in place.
- **Manifest race.** Two processes registered files concurrently and 2 rows were lost. `src/pdac/phaseD_manifest_check.py` found and re-registered them; one stale row pointing to a deleted log was removed [`results/pdac/phaseD/manifest_check.csv`].
- **Push of `pdac-v2`: BLOCKED for commits after 161418a.**
  - Pushing the Phase B commit (136 files, about 86 MB of results) failed twice with `curl 56 schannel: SEC_E_MESSAGE_ALTERED`.
  - A third attempt with Git's OpenSSL backend (certificate verification on) failed with "remote end hung up".
  - TLS was not disabled.
  - Phase A and the two rule commits are on the remote. Later commits are local only. To push from another network: `git push origin pdac-v2`.
- **v1 BLOCKED items not revisited:** the full-catalogue radius sweep and the full-section Squidpy ligrec p-values.

## Limitations

- **"Validated" T cells** = agreement of two RNA-based annotations, weakly supported by protein. There is no ground truth.
- **consensus_T is post-hoc.** It was defined after the pre-registered rule failed, though before any CCC was run on it.
- **No independent replication.** All results come from one CosMx dataset: 16 sections from 13 tumours, which are not fully independent.
- **Q is unstable** when M is small, and with 500 permutations the minimum add-one p is 0.002. Non-significance after BH limits what can be detected; it does not show that no effect exists.
- **Background flags.** Checkpoint axes remain flagged as background-level (v1, 24/31 axes). The CXCL12–CXCR4 axis is not flagged.

## Five most important findings (including negative ones)

1. **The CAF→T proximity-specific (Q) signal does not survive.** It is 0/16 after BH within section under every annotation, including the post-hoc `consensus_T`. Uncorrected it is 2/16 (consensus).
2. **The edge score mostly mirrors the cell-type product.** Spearman(z_E, z_M) is 0.66–0.68 (0.80–0.81 on unflagged axes); 77% of edge-significant rows are M-significant; about 10% are Q-significant after BH.
3. **The CosMx protein channels cannot fix T-cell label purity in this dataset.** The pre-registered rule's T precision is 0.022 in U5-a and 0.077 pooled, worse than marker scoring (0.276 / 0.353). CD3 is bimodal in only 3/16 sections.
4. **The independent replication is BLOCKED.** The GSE310352 licence is unconfirmed and its labels need R. The pipeline is ready behind an explicit gate.
5. **The p-value convention changes BH calls but not uncorrected calls.** Legacy p = 0 inflates BH significance (4,530 PDAC and 3 breast verdicts). Headline CAF→T edge calls after BH drop from 16/16 to 13/16 (marker) and 14/16 (author) under add-one. No uncorrected breast or headline verdict flips.

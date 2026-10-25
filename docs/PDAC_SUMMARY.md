# PDAC Extension — Summary (2026-10-24)

This is an audit/diagnostic summary. Nothing in it supports a biological or clinical claim. Every number comes from a file listed in `results/pdac/MANIFEST.csv`; the source file is given in brackets. The step-by-step record is in `docs/PDAC_PROGRESS_LOG.md`.

## What was done

1. **Phase 0 — hygiene.** Fixed the PDCD1 caption in `ablation.html` (12 → 5 of 7,163 cells). Replaced the hard-coded `F:\` paths in 5 scripts. Wrote `docs/RESULTS_SETTINGS_NOTE.md`, which separates 15 µm contact-restricted results from radius-sweep peaks. Renamed the margin/core split to "internal regional consistency check". The paper's 15 µm numbers reproduce exactly: CD86→CTLA4 z = 3.56, p = 0.035 [`results/pdac/phase0/multi_axis_15um_recheck.csv`].
2. **Phase 1 — matched dyadic baseline (breast).** Four scores share one label-permutation null (500 perms, seed 42, 6-NN):
   - compositional C = (L̄s + R̄r)/2
   - edge E (mean L·R over s→r edges)
   - matched M = L̄s·R̄r (product form, no graph)
   - ratio Q = E/M

   C vs M isolates score form; M vs E isolates spatial restriction.
3. **Phase 2 — dataset vetting.** 9 candidates vetted [`results/pdac/dataset_vetting.csv`]. Selected **Mendeley PDAC-SMI** (Shiau, Cao et al.; CosMx; CC BY 4.0; doi 10.17632/kx6b69n3cb.1): 717,493 cells, 16 sections from 13 tumours, 960-gene base panel, 37/39 catalogue genes present.
4. **Phase 3 — audit on PDAC.**
   - Same pipeline and QC thresholds as breast.
   - Three annotations: marker scoring, Leiden consensus, and the authors' labels as a reference.
   - 6-NN with C/E/M/Q and 500 permutations.
   - k ∈ {4, 6, 8, 12}; radius 15–300 µm (sensitivity only).
   - Primary section U5-a, plus all 16 sections as a cross-section consistency check.
   - Added two diagnostics: label purity against the author annotation, and gene detection against negative-control probes.
5. **Phase 4 — benchmark.**
   - Squidpy `gr.ligrec`, COMMOT 0.0.3 and SpatialDM 0.3.1 on the same section.
   - COMMOT ran on a deterministic 8,000-cell crop because it needs a dense n × n distance matrix.
   - Both tools were installed in a cloned env, so the paper env is unchanged.

## Results table: breast vs PDAC, per axis

Breast: 6-NN, per-cell labels [`results/pdac/phase1/`]. PDAC: primary section U5-a, marker labels, 6-NN, plus counts over 16 sections [`results/pdac/phase5/breast_vs_pdac_axes.csv`]. Significance = add-one permutation p < 0.05, uncorrected unless stated.

| Axis | Breast (s→r): z C / M / E / Q | Breast verdict | PDAC U5-a (s→r): z C / M / E / Q | PDAC E sig, sections (marker / author) | PDAC Q sig, sections (marker / author; BH) | Flags |
|---|---|---|---|---|---|---|
| CXCL12→CXCR4, stroma→immune | fibroblast→immune: 31.0 / 67.4 / 29.5 / 0.45 | C, M, E sig; **Q not** | CAF→T_cell: 33.2 / 36.1 / 23.4 / 2.81 | 16/16 / 16/16 | 4/16 / 2/16; **0/16 after BH** | marker T-cell label impure |
| CXCL12→CXCR4, stroma→tumour | fibroblast→epithelial: 42.9 / 42.7 / 4.6 / −14.6 | C, M, E sig; Q not | CAF→tumour: −20.2 / −19.7 / −15.2 / −8.9 | 0/16 / 0/16 | 1/16 / 1/16 | — |
| CD274→PDCD1, tumour→immune | epithelial→immune: not evaluable (PDCD1 in 0/361 immune cells) | — | tumour→T_cell: 16.8 / 17.2 / 17.3 / 3.02 | 11/16 / 8/16 | 4/16 / 5/16 | **background-level genes** (PDCD1 fold 0.96 over neg-probe max) |
| CD86→CTLA4, APC/immune→T | immune→immune: 17.9 / 25.9 / 1.85 / −0.07 | C, M sig; E, Q not | myeloid→T_cell: 42.5 / 81.2 / 23.1 / 0.40 | 12/16 / 16/16 | 0/16 / 0/16 | **background-level genes** (CTLA4 1.00, CD86 1.08) |
| PTN→SDC4, stroma→tumour | fibroblast→epithelial: −2.6 / −5.3 / −0.54 / 7.47 | Q sig only (negative M z) | not on panel | — | — | — |

Breast CD86→CTLA4 at 15 µm (paper setting, radius graph, 200 perms) is z = 3.56, p = 0.035 [`results/pdac/phase0/multi_axis_15um_recheck.csv`]. That is a different setting from the 6-NN row above.

## Which results held under the matched baseline

- **Breast (Phase 1):**
  - The matched baseline explains **3/9** compositional-vs-edge verdict disagreements (score form). The other 6 are proximity-type, so the pre-registered gate ("explains most") is **not met** [`results/pdac/phase1/phase1_gate.json`].
  - 4/9 disagreements are CD274–PDCD1, a sparse-gene case (PDCD1 in 5 cells), and 2 more rest on fewer than 30 edges [`results/pdac/phase1/phase1_gate_details.json`].
  - **Every edge-significant breast interaction (8/8 paper, 4/4 extra) is also matched-significant, and none has a significant ratio Q.** In breast, the edge score adds no detectable proximity-specific enrichment beyond the cell-type mean product. Its value is in *removing* compositional calls.
- **PDAC:**
  - In U5-a (marker labels), 141 of 385 edge-significant rows also have a significant Q [`results/pdac/phase3/mendeley_pdac_smi__U5-a_phase3_summary.json`]. Across 16 sections, 2,127 of 6,058 edge-significant section-rows are Q-significant.
  - Most of those rows involve background-level genes (24/31 axes flagged) [`results/pdac/phase3/background_flag_per_axis.csv`].
  - For the headline CXCL12→CXCR4 CAF→T axis, Q is significant in 4/16 sections (marker) and 2/16 (author), with 0/16 after BH within section. The edge signal is largely carried by M.

## What was sensitive to annotation, radius and k

- **Annotation** [`results/pdac/phase3/mendeley_pdac_smi__U5-a_annotation_sensitivity_summary.csv`]: marker vs Leiden, 496 shared rows in U5-a.

  | Score | Spearman z | verdict changes |
  |---|---|---|
  | C | 0.92 | 49 |
  | M | 0.89 | 50 |
  | E | 0.83 | 77 |
  | Q | 0.65 | 106 |

  Q is the most annotation-sensitive score. Leiden consensus produced only 4 types in U5-a, with no T-cell class, so the CAF→T headline is untestable under Leiden. Our marker "T_cell" class (945 cells) contains only 290 author-called T/NK cells, and 467 of the authors' 879 T cells are labelled CAF by marker scoring [`results/pdac/phase3/mendeley_pdac_smi__U5-a_author_vs_marker_crosstab.csv`]. Even so, the headline E result held under author labels (16/16 sections) [`results/pdac/phase3/author_annotation_headline_summary.csv`].
- **k** [`results/pdac/phase3/mendeley_pdac_smi__U5-a_k_sensitivity_summary.csv`]: edge z Spearman vs k = 6 is 0.945–0.974. Edge verdict changes are 80–134 of 3,100; ratio verdict changes are 161–239. The headline CAF→T z_Q rises 2.54 → 4.11 from k = 4 to 12 as edges grow from 1,429 to 4,274.
- **Radius (sensitivity only)** [`results/pdac/phase3/mendeley_pdac_smi__U5-a_radius_sensitivity_headline_pair.csv`]: for CAF→T, Q falls from 1.74 (15 µm) to about 1.5 (150–300 µm), while z_Q rises from 2.09 to 6.84 as edges grow about 40×. In breast and PDAC alike, larger neighbourhoods inflate z through edge count. Radius-sweep z must not be read as effect size.

## Benchmark agreement (Phase 4)

On the 8,000-cell crop, at cluster level (3,100 rows) [`results/pdac/phase4/comparison_*_crop.csv`]:
- Spearman between our edge z and ligrec mean is 0.02; between our edge z and COMMOT is 0.12; between ligrec and COMMOT is 0.56.
- Top-20 Jaccard between method families is ≤ 0.053.
- Significance kappa: ligrec–compositional 0.45, ligrec–edge 0.00, COMMOT–edge 0.44.

Agreement follows score design, and the top-ranked interactions are method-dependent. At L–R level (31 pairs, full section), SpatialDM global z correlates 0.26 with our max edge z and 0.54 with ligrec max mean [`results/pdac/phase4/comparison_lr_level_spearman_full.csv`].

## What failed or was blocked

- **Full-catalogue radius sweep (PDAC):** the host stopped the run for low system memory at the 200 µm step. Saved: 6-NN, annotation and k analyses, and radius nulls for 15–150 µm. A headline-only sweep (15–300 µm) was completed, and its 15–150 µm nulls reproduce the stopped run exactly. The full sweep is **BLOCKED pending a user go-ahead** to re-run.
- **Squidpy ligrec, full section:** all 3,100 p-values were NaN in two attempts (parallel and serial numba kernels). Means are valid. **BLOCKED**, cause not identified. Crop-level ligrec worked.
- **Web vetting:**
  - 10x Genomics pages returned HTTP 429, so 10x licences, sizes and URLs are UNVERIFIED.
  - The GEO terms page and PMC sat behind reCAPTCHA, so GEO data licences are not stated.
  - A Cancer Letters 2026 page returned 403.

  None of these was worked around.
- **COMMOT install:** the first attempt upgraded numpy to 2.2.6 and import failed (`np.Inf` removed). The second attempt pinned numpy < 2 and works.

## Limitations

- **Breast data:** 10x's format-validation subset (one section). PDAC: one public CosMx dataset, analysed section by section. The 16 sections come from 13 tumours, so the cross-section check is not fully independent.
- **Background detection:** 24/31 PDAC axes involve a gene detected at background level relative to the negative probes. That includes **all checkpoint axes**, so the PDAC checkpoint results are a "background-level gene diagnostic", not signal.
- **Label purity:** marker-scoring labels are coarse and, for T cells, impure. Conclusions about T-cell receivers depend on annotation.
- **Pixel scale:** 0.18 µm/px is inferred from the source paper's Methods ("60 pixels (10.8 μm)"). It affects only the radius sweep.
- **QC:** the breast QC thresholds removed 0 cells in PDAC, because the authors had already filtered.
- **Ratio Q:** unstable when M is small. All Q-significant breast rows have negative matched z.
- **Benchmark scope:** COMMOT parameters (`dis_thr` = 50 µm) and SpatialDM parameters (l = 25 µm) were pre-declared but not tuned. Method agreement may depend on them.
- **p-value resolution:** the minimum add-one p with 500 permutations is 0.002.

## Next steps (ordered by value)

1. **Make background-awareness part of the audit.** Add the negative-probe detection flag (and, optionally, a background-subtracted expression) to the pipeline. Re-score the PDAC checkpoint axes only where genes clear background.
2. **Fix T-cell label purity.** Use the CosMx protein channels already in the file (`Mean.CD3`, `Mean.CD45`, `Mean.PanCK`) or the author labels to define receivers. Re-test CAF→T under a validated T-cell class.
3. **Independent PDAC replication** on GSE310352 (CosMx, 8 untreated sections, same 37/39 gene coverage) once its licence is confirmed.
4. **Re-run the full-catalogue radius sweep** with memory headroom (or chunked per L–R pair) — needs user go-ahead.
5. **Build a minimal reproducer for the Squidpy ligrec all-NaN p-values** on the full section and report it upstream, then re-run the full-section significance comparison.

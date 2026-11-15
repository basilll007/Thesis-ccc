# v3 Phase 5 plan — Visium gene programs on GSE277783 (GATED)

Written 2026-11-15 BEFORE any Phase 5 computation. Gate: a GEO SOFT or series-matrix file for GSE277783 must exist in `data/`
with verifiable anatomic-site labels per GSM. Checked at run time by `src/pdac/v3_phase5_gate.py`.
Status at plan time: `eda_gse277783/EDA_REPORT.md` (main branch) states that no SOFT/series-matrix file exists in `data/` and
that the organ site is UNKNOWN for 49/57 slides -> expected BLOCKED.

## Pre-declared program gene lists (fixed now; only genes present in the Visium matrix are used; a program needs >= 5 present genes)
- CAF: COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRB, ACTA2, POSTN, THY1, SFRP2, VCAN
- EMT: VIM, FN1, CDH2, ZEB1, ZEB2, SNAI1, SNAI2, TWIST1, MMP2, SPARC, TGFBI, S100A4
- Hypoxia: VEGFA, SLC2A1, CA9, LDHA, PGK1, ADM, NDRG1, BNIP3, ANKRD37, P4HA1, ENO1, HK2
- Interferon: ISG15, IFI6, IFI27, IFIT1, IFIT2, IFIT3, MX1, OAS1, STAT1, IRF7, IFI44L, BST2
- Immune infiltration: PTPRC, CD3E, CD3D, CD2, CD8A, CD4, MS4A1, CD79A, CD68, CD163, LYZ, NKG7
These are commonly used marker sets assembled from general knowledge; they are not taken from one specific publication
(no citation is claimed).

## Method if the gate opens
- Per slide: spots with >= 500 UMIs; slide excluded if > 50% of spots are under 500 UMIs (pre-set rule).
- Score = mean log1p(normalize_total 1e4) of the program genes per spot, minus the mean of 50 expression-matched control genes
  (scanpy `score_genes`, seed 42); aggregate = median spot score per slide, then mean per patient x site.
- Compare primary vs metastasis with patient as the unit (paired Wilcoxon where a patient has both; otherwise Mann-Whitney),
  add-one not applicable (exact tests); report n patients per arm. Covariates: median UMI per slide (depth) and chemistry
  (Visium v1 vs CytAssist) in an OLS of slide score on site + depth + chemistry (exploratory).
- Outputs under `results/pdac/v3_phase5/`; figure `figures/pdac/v3_phase5/phase5_programs.{png,pdf}`.

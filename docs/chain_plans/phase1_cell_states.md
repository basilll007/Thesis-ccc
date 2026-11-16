# Chain Phase 1 plan: cell states (what "changed" means)

Committed before any Phase 1 result is computed. Seed 42. Data: Mendeley PDAC-SMI CosMx, 16 sections and 13 tumours. The section-to-tumour map is `v3_nulls.section_to_tumour()`, which strips `_2nd`/`_B` from the donor: U4-a/U4-b -> 011, U6-a_1/U6-a_2 -> 2626, T4-a/T4-b -> 2762. QC and normalisation reuse v3 `load_section` exactly: phase3 `qc` (min_counts 10, min_cells 3), then `normalize_total(1e4)` and `log1p`. The primary annotation is `author_mapped`; the sensitivity annotation is `marker` (both from results/pdac/phaseB/sections). Wording: these are expression programs, not mutations; everything is association.

## 1. Background-level genes (applies to all 960 panel genes)
- Same rule as v3 phase3f, extended from catalogue genes to every Base gene. For each section, a gene is *near background* if its detection rate (fraction of cells with count > 0) is < 2 x the maximum detection rate of the 19 negative probes in the same section.
- A gene is **background-level** if it is near background in >= 8 of 16 sections. Background-level genes are removed from every program before scoring.
- Output: `results/chain/phase1/gene_background_flags.csv`.

## 2. Pre-declared programs
Source note: these lists were written from recall of the PDAC subtype, CAF-subtype and T-cell/myeloid state literature (e.g. the Moffitt-style classical/basal-like split, the Elyada/Öhlund myCAF/iCAF split, standard cytotoxic/exhaustion markers, and complement/SPP1 tumour-associated macrophage markers). **None of these source papers was fetched in this session, so every list is UNVERIFIED as to its source.** The lists were restricted to genes on the 960-gene panel before any expression was looked at; only panel membership was checked. A program is dropped if it has < 5 usable genes after the background filter.

| Cell type (author_mapped / marker class) | Program | Genes (panel members) |
|---|---|---|
| tumour (`epithelial_ductal_tumor`) | classical | CEACAM6, AGR2, S100P, KRT20, LCN2, PIGR, OLFM4, SPINK1, CEACAM1 |
| tumour | basal_like | KRT5, KRT6A, KRT6B, KRT14, KRT17, S100A2, SERPINB5, AREG, SLC2A1, KRT15, LY6D, ITGA3, MET |
| tumour | EMT | VIM, SNAI1, SNAI2, TWIST2, FN1, MMP2, VCAN, CAV1 |
| CAF (`CAF_fibroblast`) | myCAF | ACTA2, TAGLN, MYL9, TPM2, TPM1, COL11A1, TWIST1, SOX4 |
| CAF | iCAF | IL6, CXCL1, CXCL2, PDGFRA, CFD, LMNA, CCL2, LIF, CXCL8 |
| CAF | CXCL12 (single-gene readout) | CXCL12; named explicitly in the brief, so kept as a gene-level readout exempt from the >= 5 rule and labelled as such everywhere |
| T cell (`T_cell`) | cytotoxic | GZMA, GZMB, GZMK, GZMH, PRF1, NKG7, GNLY, IFNG, CCL5, CST7, CD69 |
| T cell | exhaustion | PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, TOX, ENTPD1 |
| myeloid (`myeloid`) | inflammatory | IL1B, TNF, CCL3, CCL4, CXCL9, CXCL10, NFKBIA, PTGS2 |
| myeloid | C1Q_TAM | C1QA, C1QB, C1QC, CD163, MRC1, SELENOP |
| myeloid | SPP1_TAM | SPP1, GPNMB, CD9, FABP5, LGALS3 |

**Lineage marker sets.** These are used only as the spillover covariate in Phases 2 and 3, not as outcomes: tumour {KRT19, KRT8, KRT18, EPCAM, KRT7}; CAF {COL1A1, COL1A2, COL3A1, DCN, LUM}; T {CD3D, CD3E, CD3G, CD2, CD8A}; myeloid {CD68, CD14, LYZ, CSF1R, FCGR3A}. The same background filter applies.

## 3. Scoring (documented normalisation)
- Per section, on the log-normalised matrix, use a scanpy-`score_genes`-equivalent re-implementation so that control genes are recorded and can be reused in the Phase 4 simulations. Genes are binned into 25 bins by rank of mean expression over all cells of the section. For each program gene, 50 control genes are drawn from its bin (excluding the program's genes) with `numpy.random.default_rng(42)`. The control set is the union of the draws.
- score = mean(program genes) - mean(control genes). Every program is scored for every cell, whatever its type; analyses use only the cells of the program's own type.
- The CXCL12 readout is its log-normalised expression minus the mean of its 50 controls.
- Outputs:
  - `program_genes_used.csv`: declared genes, panel / background / used status and n_used per program.
  - `program_controls.csv`: the controls for each section and program.
  - Per-section and per-tumour distribution tables (n, mean, sd, quantiles 5/25/50/75/95) for each type x program.

## 4. Validity check (sanity only; cannot change the lists)
- Among author `Malignant` cells: AUC of (classical - basal_like) for author subtype CLS vs BSL.
- Among author `CAF` cells: AUC of (myCAF - iCAF) for myCAF vs iCAF.
- Reported per section. Expected direction: AUC > 0.5. A failure is reported as a validity concern and does not lead to a change of lists.

## 5. Cell counts and floors
- Graph: 6-NN in um coordinates (phase3 `knn_edges`), symmetrised: i~j if j in kNN(i) or i in kNN(j). This graph is used in Phases 1 and 3. Phase 2 additionally uses a 30 um radius graph.
- Cell counts per type per section, for both annotations.
- **Floor, per section and candidate chain A -> B -> C** (pre-declared default): >= 100 cells of each of A, B and C, AND >= 300 B-C edges AND >= 300 A-B edges in the symmetrised 6-NN graph. These edges are the "qualifying neighbour relations".
- Candidate chains: all 24 ordered triples of distinct types from {tumour, CAF, T cell, myeloid}. Pass/fail is recorded per section, chain and annotation.

## 6. Chain selection rule for Phase 3 (outcome-free; applied by the Phase 1 script)
- Primary chain (fixed): tumour (A) -> CAF (B) -> T cell (C).
- Excluded from secondary selection: the primary chain, and its reversal T -> CAF -> tumour, which is a Phase 3 negative control.
- Secondary chains: up to 5 of the remaining 22 triples. Rank by (1) the number of sections passing the floor under author_mapped (descending), (2) the number of distinct tumours passing (descending), and (3) the summed minimum type count over passing sections (descending). Take the top 5. A chain passing in < 3 tumours is not selected (minimum for any tumour-level inference). Output: `results/chain/phase1/chain_selection.csv`.
- This uses only cell counts and edge counts. No program score enters the selection.

## 7. Outputs and figure
- results/chain/phase1/: gene_background_flags.csv, program_genes_used.csv, program_controls.csv, cell_counts_by_section.csv, program_distribution_by_section.csv, program_distribution_by_tumour.csv, floor_by_section_chain.csv, chain_selection.csv, validity_auc.csv, phase1_summary.json.
- Cell-level cache (scores, coordinates, labels, graphs): data/chain_cache/ (git-ignored, deterministic, re-derivable; not results).
- Figure figures/chain/phase1/: per-tumour program score distributions for each type, plus the floor pass matrix.

## Decision rules
- Programs with < 5 usable genes: dropped (CXCL12 readout excepted, as stated).
- Sections failing the floor for a chain: excluded from that chain's Phase 3 test. This is reported, not hidden.

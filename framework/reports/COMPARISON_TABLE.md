# Positioning: cellalibi versus existing spatial CCC and segmentation tools

Plan: `framework/plans/v2_phase7.md`. All sources were fetched on 2026-11-17/18, by me or by a research subagent in
this session. I re-fetched three key items myself:
- the Europe PMC record of bioRxiv 10.1101/2025.01.02.631135;
- the Statial vignette `contamination = TRUE` passage;
- the Statial `minCells = 100` setting.

Legend: **Y** yes · **P** partial (the limited form is stated) · **N** no, meaning absent from the documentation or
paper sections that were read (it does not mean every section was checked) · **U** UNVERIFIED (the fetched source
does not settle it) · **n/a** the tool performs no interaction test, so the property does not apply.

Columns are defined in the plan:
1. Abundance control
2. Spatial/zonation null
3. Spillover check
4. Calibration / power floor
5. Per-claim verdict naming the artefact
6. Patient-level inference

| Method | 1 Abundance control | 2 Zonation null | 3 Spillover check | 4 Calibration / power floor | 5 Per-claim verdict | 6 Patient-level | Source(s) |
|---|---|---|---|---|---|---|---|
| **cellalibi v0.2 (this work)** | Y: contact tests permute expression within cell type, so type abundance and type-level expression are preserved (`contact.py`) | Y: within-tile receiver contrast + ring control against nearby non-touching receivers (`contact.py`) | Y: claim-specific sender-marker contamination and expected receptor excess (`contact.spillover_check`) | Y: null calibration test (200 datasets x 1000 permutations); plant-and-shift power floor, NOT_TESTABLE below it | Y: power / abundance / zonation / spillover PASS-FAIL-NOT_RUN plus a primary verdict naming the artefact (`verdicts.audit_claim`) | P: per-tumour aggregation of section verdicts with a majority rule (Phase 6); no mixed model | this repository, `framework/src/cellalibi/`, `framework/results/v2_phase2/` |
| SpatioMark / Statial | P: fits within one cell type; abundance is an alternative predictor ("function of abundance or proximity"), not a control | N: fitted "individually for each image, pairwise cell type, and marker"; no region term | Y: cell-type probability covariates (`contamination = TRUE`), described as "not a perfect solution" | P: relationships with < 100 cells are excluded (`minCells = 100`); no calibration reported | N: outputs coefficients / t / FDR | P: coefficients averaged per patient, then Cox | [4][5] |
| MISTy (mistyR) | P: intraview baseline "independent of the spatial localization"; cell types not used | N: only an unconstrained coordinate permutation | N: not mentioned | P: cross-validation p-value of the R² gain; no data floor | N: importances / R² gain | P: results aggregated; same-patient samples treated as independent | [6][7] |
| NCEM | P: compared with a non-spatial baseline conditioned on cell type; no composition null | P: image/patient domain covariates only | P: robustness to simulated segmentation errors; no correction | P: Wald test + FDR; null simulations; no explicit floor | N: counts of significant coefficients per pair | P: patient as a domain covariate | [8] |
| COMMOT | P: cluster-level p-value from permutations of cluster labels | P: option to permute spots within each cluster | N: not mentioned | N: only "a large number of cells or spots is preferred" | N: scores; edges kept at P < 0.05 | N: not described | [9] |
| SpatialDM | N: random shuffling of binding partners; no cell-type conditioning | N: no constrained null described | N | Y: p-values "well calibrated to a uniform distribution" on simulated null data | N: global / local p-values per L-R pair | Y: regression on per-sample z-scores with and without condition | [10] |
| Squidpy `gr.ligrec` | N: CellPhoneDB-style test across cell-type combinations; no spatial term | N: no spatial input | N | P: test skipped below the `threshold` expression fraction (NaN p-values) | N: returns means / p-values | N: no sample argument | [11][12] |
| LIANA+ (spatial + multi-sample) | N: the bivariate null shuffles spot labels | P: multi-view models can mask views by region labels | N | P: expression filters (nz_prop); systematic benchmarks "still pending" | N: per-interaction tables | Y: per-sample inference, then MOFA / tensor factorisation and pseudobulk DE | [13][14][15][16] |
| Proseg | n/a | n/a | Y: segmentation that reduces "spurious coexpression and transcript misassignment" | N: per-transcript assignment probabilities; no contamination QC | N | N | [17][18] |
| segger | n/a | n/a | Y: transcript-to-cell link prediction to fix misassignment | U: README silent; full text not fetched | N | U | [19][20] |
| ResolVI (scvi-tools) | n/a (differential niche abundance tests composition rather than controlling for it) | n/a | Y: reassigns expression to neighbouring cells and background | N: posterior samples, no calibration or power criterion | N | P: optional batch/donor covariates | [21] |
| SPLIT (Bilous et al.) | n/a | n/a | Y: RCTD-based decomposition plus per-cell singlet/doublet class; Xenium only in the fetched paper | N: no reliability cutoff found | N, but it urges users to "critically evaluate the outputs of cell–cell interaction inference tools" | N: samples annotated independently | [22][23] |
| Mitchel et al. 2025 preprint (bioRxiv 10.1101/2025.01.02.631135) / cellAdmix | n/a | n/a | Y: neighbourhood factorisation isolates "molecular admixtures". The abstract says segmentation errors "confound most downstream analysis … including … inference of neighboring cell influence, and ligand-receptor interactions" | P: a Bayesian metric of how well the correction worked | N: removes molecules; no per-result test | U | [1][2][3] |

## Sources (fetched 2026-11-17/18)

1. Europe PMC record, DOI 10.1101/2025.01.02.631135: Mitchel J, Gao T, Cole E, Petukhov V, Kharchenko PV.
   *Impact of Segmentation Errors in Analysis of Spatial Transcriptomics Data.* bioRxiv 2025.
   https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1101/2025.01.02.631135&resultType=core&format=json
   (re-fetched by me).
2. Europe PMC record: Mitchel J, Gao T, Petukhov V, Cole E, Kharchenko PV. *Impact and correction of segmentation
   errors in spatial transcriptomics.* Nat Genet 2026;58(2):434-444. doi:10.1038/s41588-025-02497-4.
3. https://github.com/kharchenkolab/cellAdmix
4. Iyengar SS, Qin AR, Robertson N, Harman AN, Patrick E. *SpatioMark: quantifying the impact of spatial proximity on
   cell phenotype.* Bioinformatics 2025;41(8):btaf409. doi:10.1093/bioinformatics/btaf409. Full text via Europe PMC:
   https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12343047/fullTextXML
5. Statial vignette source: https://raw.githubusercontent.com/SydneyBioX/Statial/devel/vignettes/Statial.Rmd
   (re-fetched by me; lines on `minCells = 100` and `contamination = TRUE`).
6. mistyR vignette: https://saezlab.github.io/mistyR/articles/mistyR.html
7. Tanevski J et al. *Explainable multiview framework for dissecting spatial relationships from highly multiplexed
   data.* Genome Biol 2022;23:97. doi:10.1186/s13059-022-02663-5.
   https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9011939/fullTextXML
8. Fischer DS, Schaar AC, Theis FJ. *Modeling intercellular communication in tissues using spatial graphs of cells.*
   Nat Biotechnol 2023;41(3):332-336. doi:10.1038/s41587-022-01467-z.
   https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10017508/fullTextXML
9. Cang Z et al. *Screening cell-cell communication in spatial transcriptomics via collective optimal transport.* Nat
   Methods 2023;20(2):218-228. doi:10.1038/s41592-022-01728-4.
   https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9911355/fullTextXML
10. Li Z et al. *SpatialDM for rapid identification of spatially co-expressed ligand-receptor and revealing cell-cell
    communication patterns.* Nat Commun 2023. doi:10.1038/s41467-023-39608-w.
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10325966/fullTextXML
11. squidpy `gr.ligrec` documentation: https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.ligrec.html
12. Palla G et al. *Squidpy: a scalable framework for spatial omics analysis.* Nat Methods 2022;19(2):171-178.
    doi:10.1038/s41592-021-01358-2.
13. LIANA+ bivariate tutorial: https://liana.readthedocs.io/en/stable/tutorials/notebooks/bivariate.html
14. https://liana.readthedocs.io/en/stable/tutorials/notebooks/inflow_score.html
15. https://liana.readthedocs.io/en/stable/tutorials/notebooks/targeted.html
16. Dimitrov D et al. *LIANA+ provides an all-in-one framework for cell-cell communication inference.* Nat Cell Biol
    2024;26(9):1613-1622. doi:10.1038/s41556-024-01469-w.
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11392821/fullTextXML
17. Jones DC et al. *Cell simulation as cell segmentation.* Nat Methods 2025;22(6):1331-1342.
    doi:10.1038/s41592-025-02697-0 (Europe PMC record).
18. https://github.com/dcjones/proseg
19. Heidari E et al. *Segger: Fast and accurate cell segmentation of imaging-based spatial transcriptomics data.*
    bioRxiv 2025. doi:10.1101/2025.03.14.643160 (Europe PMC record).
20. https://github.com/EliHei2/segger_dev ; https://github.com/dpeerlab/segger
21. ResolVI user guide: https://docs.scvi-tools.org/en/stable/user_guide/models/resolvi.html (cites Ergen & Yosef
    2025, bioRxiv doi:10.1101/2025.01.20.634005).
22. Bilous M et al. *Resolving sensitivity, specificity and signal contamination in Xenium spatial transcriptomics.*
    Nat Methods 2026;23(6):1152-1162. doi:10.1038/s41592-026-03089-8.
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13259927/fullTextXML
23. Bilous M et al. preprint, bioRxiv 2025, doi:10.1101/2025.04.23.649965 (Europe PMC record).

## Caveats

- Long full texts were read only in part (the first ~100k characters). An **N** therefore means "not found in what was
  read".
- NCEM was published in Nat Biotechnol, not Nat Methods.
- The SpatialDM paper describes its null in two ways that do not obviously agree ("shuffling binding partners" vs
  "shuffling w_ij").
- The cellalibi row describes features, not proven superiority. Benchmark performance is reported in
  FRAMEWORK_REPORT_v2.md.

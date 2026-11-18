# Search log v2 (Phases 1, 3 and 7)

All dates are 2026-11-17 (local) / 2026-11-18 UTC. Blocked fetches were not retried and not worked around.

## Phase 1: hypoxia literature and GEO records

| Query / URL | Outcome |
|---|---|
| GEO brief text records GSE62452, GSE28735, GSE57495, GPL15048 (`acc.cgi?...&form=text&view=brief`) | found. The GSE62452 summary contains the "sample number 1-90 ... GSE28735" statement. |
| GEO series matrices GSE62452 / GSE28735 / GSE57495; GPL15048 `view=data` table | found (cached in `F:/Thesis/data/v2_cache/`) |
| WebSearch "hypoxia gene signature prognostic overall survival pancreatic ductal adenocarcinoma" | found Abou Khouzam et al. 2021 and Ren et al. 2023 |
| PMC copies (PMC8173254, PMC10464224) | **blocked**: reCAPTCHA page |
| Frontiers publisher page, doi 10.3389/fimmu.2021.680435 | found (Abou Khouzam et al. 2021) |
| HAL hal-03253996 | **blocked**: bot check |
| DOAJ article page | **blocked**: HTTP 403 |
| wjso.biomedcentral.com -> link.springer.com | redirected to a login (idp.springer.com); not followed |
| Europe PMC REST, DOI:10.1186/s12957-023-03142-2 | found (Ren et al. 2023 abstract) |
| Europe PMC REST, title search "Integration of Genomic and Transcriptional Features in Pancreatic Cancer ..." | found (Connor et al. 2019, Cancer Cell, PMID 30686769). A first attempt with a guessed DOI returned an unrelated paper and was discarded. |

## Phase 3: name availability

| Query | Outcome |
|---|---|
| https://pypi.org/pypi/cellalibi/json | 404 (free) |
| https://pypi.org/pypi/cell-alibi/json | 404 (free) |
| https://api.github.com/search/repositories?q=cellalibi | total_count 0 |
| https://api.github.com/search/users?q=cellalibi | total_count 0 |
| GitHub code search | not run: `gh` not authenticated |

## Phase 7: tool documentation and novelty searches

The research subagent's queries, summarised from its report:

| Query / URL | Outcome |
|---|---|
| Europe PMC record and full-text lookups for SpatioMark, MISTy, NCEM, COMMOT, SpatialDM, Squidpy, LIANA+, Proseg, segger, SPLIT (both versions), Mitchel et al. (preprint and Nat Genet version), CellWHISPER, Bhadury et al., TRACER, InSituCor, Wu et al. review | found |
| Bioconductor Statial HTML | truncated; the vignette .Rmd on GitHub was used instead (re-fetched by me) |
| liana-py.readthedocs.io | 301 redirect to liana.readthedocs.io; followed; found |
| mistyR vignette; squidpy ligrec docs; Proseg, segger, cellAdmix GitHub; scvi-tools ResolVI guide | found |
| bioRxiv full text of SpatioMark; bioRxiv page for 10.1101/2025.01.02.631135v1; TRACER bioRxiv page | **blocked**: HTTP 429 |
| WebSearch "SPLIT spillover correction Xenium CosMx segmentation RCTD" | found Bilous et al. |
| WebSearch for the Mitchel paper | found cellAdmix; a proxy mirror also appeared and was **not used** |
| WebSearch "cell-cell communication spatial transcriptomics segmentation spillover false positive" | reviews and Proseg; no audit tool |
| WebSearch "ligand receptor spatial confounding abundance null calibration" | not relevant (SpatialDM, an MR-based CCC paper, an LR-evaluation audit) |
| WebSearch "audit cell-cell communication claims spatial transcriptomics artefacts verdict" | not relevant (reviews, a proximity-mismatch benchmark) |
| WebSearch "spatial colocalization confounding cell type composition null model" | InSituCor, KAMP, C-SIDE |
| Europe PMC query combining CCC, spatial, segmentation, confounding and abundance | CellWHISPER, Bhadury et al. |
| WebSearch (extended) for a diagnostic framework | TRACER, CONCISE, a 2026 segmentation benchmark |
| WebSearch combining zonation, segmentation, abundance and a per-interaction test | no match |

Seen only in search results, not fetched: arXiv 2606.09675; https://www.nature.com/articles/s41592-026-03272-x;
CONCISE (PMC13320749).

## Closest related work found (abstracts fetched by the research subagent; not re-fetched by me)

None of these audits a specific claim against all three artefacts.

| Work | What it covers |
|---|---|
| CellWHISPER (Kumar et al., bioRxiv 2026, doi:10.64898/2026.01.07.697982) | cell-type-specific expression and spatial structure; outputs significance scores |
| Bhadury, Frankel, Rao (bioRxiv 2026, doi:10.64898/2026.10.08.757632) | composition, density gradients and patient level; no segmentation |
| InSituCor (Danaher et al., Genome Biol 2025, doi:10.1186/s13059-025-03554-1) | conditions on the cell-type landscape |
| TRACER (Yuan et al., bioRxiv 2026, doi:10.64898/2026.03.08.710395) | segmentation purity/conflict metrics only |
| Wu, Su, Cui (Brief Bioinform 2026, doi:10.1093/bib/bbag477) | review noting that heterogeneity and spatial architecture bias permutation inference |

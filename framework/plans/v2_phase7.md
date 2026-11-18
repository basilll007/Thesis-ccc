# v2 Phase 7 plan: positioning (comparison table, novelty statement, search log)

Committed before the comparison table is filled.

## Rows

cellalibi, SpatioMark/Statial, MISTy, NCEM, COMMOT, SpatialDM, Squidpy `gr.ligrec`, LIANA+, the
segmentation-correction tools Proseg, segger, ResolVI and SPLIT, and the 2025 segmentation-error preprint
(bioRxiv doi 10.1101/2025.01.02.631135).

## Columns, each with a strict definition

1. **Abundance control.** Does the method separate a spatial/contact signal from cell-type abundance or type-level
   expression? Examples: a null that preserves type composition, or explicit conditioning on abundance.
2. **Spatial/zonation null.** Does it use a spatially constrained null, or otherwise control for regional
   co-variation (zonation) when testing an interaction?
3. **Spillover check.** Does it diagnose or correct segmentation spillover / transcript misassignment between
   neighbouring cells?
4. **Calibration / power floor.** Does it report a null calibration or a minimum-data/power criterion telling the user
   when a test is uninformative?
5. **Per-claim verdict.** Does it return a verdict for a specific claim, naming which artefact explains it, rather
   than a score or p-value per interaction?
6. **Patient-level inference.** Does it treat patient/sample as the unit of inference (multi-sample models,
   patient-level aggregation)?

## Cell values

- **Yes**, **Partial** or **No**, each with a short quote or paraphrase and the URL of the documentation page or paper
  that supports it.
- If the fetched source does not settle the question, the cell is **UNVERIFIED**, never guessed.
- The cellalibi row is filled from this repository (code and reports), and says so.

## Novelty statement and search log

- Five sentences, using "to our knowledge". Only claims supported by the table.
- A search log with date, query, source fetched and what was found or not found, including blocked fetches
  (CAPTCHA / 403). These are not worked around.

## Outputs

- `framework/reports/COMPARISON_TABLE.md`
- `framework/reports/SEARCH_LOG_v2.md`

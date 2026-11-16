# v3 resume: execution details (2026-11-15)

This addendum is written before resumed Phase 4/5 outcome computation. The
existing phase plans and program gene lists remain the scientific specification.

## Inputs and authorisation

- The user authorised continuing the unfinished work and explicitly confirmed
  having reviewed the GSE310352 terms. The confirmation is recorded in
  `LICENSE_CONFIRMED_GSE310352.txt`.
- GEO family SOFT metadata for both series were downloaded from NCBI HTTPS;
  URLs, sizes and SHA-256 digests are in `data/geo_metadata_download_manifest.json`.
- GSE310352: all 16 expression/metadata CSVs were extracted from the existing
  archive; their sizes match the GEO file list, gzip checksums pass, and SHA-256
  digests are recorded in `data/GSE310352_install_manifest.json`.
- GSE277783: 57 existing Visium matrices match GEO sample accessions and tissue
  labels. Seven additional CosMx samples in the SuperSeries metadata are not
  Visium inputs. Two PDX slides are excluded from patient tissue inference.

## Phase 5 implementation choices resolving unspecified details

- Run the repository's existing pinned Python environment and Scanpy scoring
  method, without substituting a new workflow engine or analysis method.
- Match samples to SOFT by GSM accession. Parse patient number from the GEO
  sample title (`Pt-N`); preserve both title and parsing evidence. Use GEO tissue
  as anatomical site, never infer organ from filename suffix letters.
- Chemistry is the existing EDA's inferred Visium v1/CytAssist batch, supported
  by the spatial archive inventory (`spatial_enrichment.csv` in CytAssist
  batches); it is not independently verified manufacturer metadata.
- Apply only the planned UMI thresholds: retain spots with >=500 UMIs and
  exclude a slide if >50% of its original spots fail this criterion. PDX slides
  are recorded in QC and excluded from program comparisons.
- Use Scanpy 1.10.3 `score_genes`, seed 42, `ctrl_size=50`, `n_bins=25`,
  `ctrl_as_ref=True`, `use_raw=False`, after normalize_total(1e4) and log1p.
  Scanpy samples up to 50 controls **per expression bin**, not 50 in total;
  record the actual controls. Do not silently replace this declared method
  with a different 50-total-controls algorithm. A program needs >=5 present
  genes. Retain the reader's standard unique gene-name suffixes.
- Use slide medians, then equal-weight slide means per patient x anatomical
  site. For the combined metastasis arm, average those site means within each
  patient so patients and sites cannot gain weight from extra slides.
- Primary vs combined metastasis is the main comparison. Use paired complete
  patients when available; report unmatched patients separately. If no pairing
  exists use a two-sided Mann-Whitney test. Also report organ-specific paired
  comparisons as exploratory. P-values are exact exhaustive permutations for
  Wilcoxon (including ties); BH correction across the five main programs, and
  separately across the exploratory organ-by-program family.
- Fit the planned exploratory slide-level OLS with anatomical site, raw median
  UMI depth and chemistry, using patient-clustered standard errors. Report design
  rank and chemistry/site imbalance. This does not establish causation or
  cell-intrinsic programs in mixed Visium spots.
- Save spot scores, controls, QC, slide and patient aggregates, tests, covariate
  results, input digests and a PNG/PDF figure; register all outputs. Resume only
  against matching input identities and parameters.

## Phase 4 implementation choices

- Use the saved Phase 3 QC, preprocess and annotate functions: marker and
  Leiden-consensus annotations, with their limitations made explicit.
- Use 6-NN, 500 permutations, seed 42, the existing C/E/M/Q computation and
  add-one empirical p-values; apply BH across all catalogue rows per slide,
  annotation and null, as in Phase 1.
- Coordinates remain pixels. No verified um/px exists: label blocks by 5x,
  10x and 20x the median positive nearest-neighbour distance, not micrometres.
  Distance strata likewise use nearest-neighbour units (original Phase 1 band
  ratios scaled to that distance); record the actual pixel boundaries.
- Preserve original cell metadata and negative-probe summaries. Evaluate each
  axis against negative-probe background diagnostically; do not convert missing
  or unvalidated labels into proof of absence.
- GEO metadata do not establish patient non-overlap with Mendeley. Call the
  outcome a cross-dataset consistency audit, with independence UNVERIFIED;
  do not label it independent replication.

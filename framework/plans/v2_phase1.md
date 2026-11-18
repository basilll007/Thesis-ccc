# v2 Phase 1 plan: correct the hypoxia validation and run a fresh confirmatory test

Branch `framework-v2` (from origin/framework-v1 at 8a9d8ea). Written and committed **before** any
GSE57495 expression value is scored or related to survival, and before the GSE28735/GSE62452
mapping is computed.

What has been looked at before this plan (so it is not "blind"):
- GEO text records (brief view) of GSE62452, GSE28735, GSE57495 and GPL15048, and the sample
  metadata rows (titles, characteristics, CEL file names) of the three series matrices.
  Cached under `F:/Thesis/data/v2_cache/` (outside git).
- GSE57495 covariates visible in the series matrix: `overall survival (month)`, `vital.status`,
  `Stage`. No age, sex, grade or purity fields were seen.
- GSE57495 survival has never been related to any signature in this project (cohort_vetting.csv
  only counted rows: 63 tumours, 42 events).
- GSE62452 and GSE28735 hypoxia results were seen in Phase B (framework-v1).

## 1a. Overlap between GSE28735 and GSE62452

Fact to verify: GSE62452 `!Series_summary` (fetched 2026-11-17) says "Affymetrix data from
sample number 1-90 were earlier submited by us as GEO accession#: GSE28735" (same contributor,
S.P. Hussain).

Method (all three are computed and reported; a tumour is called "the same patient sample" if
rule E holds, and rules T and S are reported as corroboration):
- E (expression identity): both series are GPL6244. Over all probes present in both matrices,
  compute the Pearson correlation of every GSE28735 sample with every GSE62452 sample. A pair is
  an identity match if r > 0.99 **and** it is the best match for both samples (mutual best).
  GSE62452 was batch-corrected in Partek, so if no pair passes r > 0.99 the Spearman correlation
  is also reported and the threshold result is reported as-is (no threshold change).
- T (title / order): GSE62452 CEL file names carry a running sample number (`<n>_<id>_T/N`).
  Report whether GSE28735 sample k maps to GSE62452 sample number k (1..90).
- S (survival): report whether matched pairs have equal survival months (GSE28735 is rounded to
  integer months) and equal event status.

Output: `results/translational/phaseB/v2_correction/gse28735_in_gse62452_map.csv`
(one row per GSE28735 sample) and a count of GSE28735 tumours (and tumours with OS) inside
GSE62452.

Decision rule: if >= 1 GSE28735 tumour is inside GSE62452, the Phase B two-cohort
meta-analysis, its I2, its pooled-events power calculation and the 4-cohort secondary pooled
result are declared invalid (patients counted twice). GSE62452 alone then becomes the only
independent evidence, and this is labelled a **post-hoc correction** (GSE62452 results were
already seen).

## 1b. Correction note

`results/translational/phaseB/HYPOXIA_VALIDATION_v1_CORRECTION.md` (new file; the v1 report is
not edited). Contents: the overlap count, which v1 numbers are invalid, GSE62452 alone
(stage-adjusted HR 1.48 [1.08, 2.03], as recorded in `hypoxia_models.csv`), prior literature
showing hypoxia signatures are prognostic in PDAC (only pages actually fetched; otherwise
UNVERIFIED), and the GSE57495 result.

## 1c. Confirmatory test in GSE57495 (pre-registered)

Hypothesis: higher tumour expression of the 12-gene hypoxia signature is associated with worse
overall survival. H0: beta <= 0 vs H1: beta > 0, Wald test, **one-sided alpha = 0.05**.
This is an association, not a causal or clinical-utility claim.

Gene list (unchanged from Phase B): VEGFA, SLC2A1, CA9, LDHA, PGK1, ADM, NDRG1, BNIP3, ANKRD37,
P4HA1, ENO1, HK2.

Probe annotation: the GPL15048 platform data table from GEO (`acc.cgi?acc=GPL15048&view=data`,
column `GeneSymbol`, NCBI-generated). A probe maps to a gene with confidence only if its
`GeneSymbol` field is exactly that symbol (no `///` multi-gene entries). Coverage = genes with
>= 1 confident probe present in the series matrix. If coverage < 10/12 (80 %) -> log
**BLOCKED** and do not score.

Probe collapsing (outcome-blind): if a gene has several confident probes, use the probe with the
highest mean expression across all GSE57495 samples (MaxMean). Sensitivity (reported, not used
for the decision): mean of the within-cohort z-scores of all confident probes of the gene.

Scale: if the matrix maximum is > 100 the values are treated as linear and log2(x + 1)
transformed; otherwise used as provided.

Scoring (unchanged from Phase B): z-score each gene across tumours within the cohort, average the
available genes, re-standardise the average to mean 0, SD 1. HR is per SD.

Samples: all GSE57495 samples with OS months > 0 and a parseable vital status
(DEAD = 1, ALIVE = 0). The series describes all 63 as pancreatic cancer tissue; no tissue filter.

Covariates: stage is the only clinical covariate in the series. Stage is binarised as in Phase B:
early = stage I/IA/IB/IIA (written `1*`, `2A`), advanced = IIB/III/IV (`2B`, `3*`, `4*`).
Primary model: Cox(OS ~ hypoxia_z + stage_advanced), complete cases, lifelines CoxPHFitter,
Breslow baseline, Efron ties (lifelines default). If stage is missing for > 20 % of samples or a
stage group has < 5 events, the primary model becomes the unadjusted model (logged). Unadjusted
model always reported as secondary.

Assumption check (reported, not a gate): Schoenfeld-residual proportional-hazards test for the
hypoxia term.

Power, declared now (from the vetting count of 42 events; recomputed after loading):
with D events and a per-SD covariate, power = Phi(log(HR) * sqrt(D) - 1.645).
For D = 42: power ~ 0.56 at HR 1.32 (TCGA estimate) and ~ 0.81 at HR 1.48 (GSE62452 estimate).
So a non-significant result is **inconclusive**, not a refutation, unless the upper 95 % CI bound
is below 1.32.

Decision table (primary model):
- one-sided p < 0.05 and HR > 1 -> CONFIRMED in an independent cohort.
- one-sided p >= 0.05, HR > 1 -> NOT CONFIRMED (direction consistent, underpowered).
- HR <= 1 -> NOT CONFIRMED (direction inconsistent).

Secondary (descriptive, labelled post-hoc because GSE62452 was already seen): inverse-variance
fixed-effect and DerSimonian-Laird pooling of GSE62452 (stage-adjusted) with GSE57495 (primary
model), the two cohorts with no shared patients. No decision is taken from it.

## 1d. Overlap check for GSE57495

Cohorts used anywhere in the project: TCGA PAAD (cBioPortal pancan), CPTAC PDA (cBioPortal 2021),
GSE62452, GSE28735.
- Identifier check: GSM accessions, sample titles and CEL file names of GSE57495 against the other
  GEO series; institution (Moffitt, contributors Chen and Malafa) against GSE62452/GSE28735
  (Hussain, NCI).
- Cross-platform expression screen: on genes measured in both cohorts (gene-level, each gene
  centred and scaled within cohort), Spearman correlation of each GSE57495 tumour with each
  tumour of the other cohort. A pair is flagged if its correlation exceeds the 99.9th percentile
  of all cross-pairs **and** is > 0.8. Cross-platform identity detection is weak; this is
  reported as a screen, not proof.
- TCGA barcodes carry a tissue-source-site code; if a source-site table is available locally it
  is used, otherwise logged UNVERIFIED.
If any overlap is found, the overlapping patients are removed before the primary model and this
is logged.

## Outputs
`results/translational/phaseB/v2_correction/` : mapping CSV, GSE57495 probe map, coverage,
cohort table, Cox results, overlap screen, run log. Code: `framework/scripts/v2/phase1_*.py`.

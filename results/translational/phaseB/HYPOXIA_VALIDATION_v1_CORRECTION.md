# Correction to HYPOXIA_VALIDATION_v1.md (post-hoc correction)

**Date:** 2026-11-17 · **Branch:** `framework-v2` · **Corrects:** `results/translational/phaseB/HYPOXIA_VALIDATION_v1.md`,
`meta_analysis_primary.csv`, `meta_analysis_secondary.csv`, `power_and_detectable_hr.csv`,
`leave_one_out.csv` and the "Phase B" rows of `FINDINGS.md` on `framework-v1`. Those files are left unchanged; this
note supersedes the parts listed below.
**Plan (committed before the new analyses):** `framework/plans/v2_phase1.md` (commit d897e19).
**Evidence files:** `results/translational/phaseB/v2_correction/`.

This is a **post-hoc correction**: the GSE62452 and GSE28735 survival results were already seen when the overlap was
found. Only the GSE57495 test below was pre-registered before its data were scored.

---

## 1. The two "independent" cohorts share patients

The GEO record for GSE62452 (`!Series_summary`, fetched 2026-11-17) says: "Affymetrix data from sample number 1-90
were earlier submited by us as GEO accession#: GSE28735." Both series have the same contributor (S.P. Hussain) and
platform (GPL6244). I checked this with the expression data (`v2_correction/gse28735_in_gse62452_map*.csv`,
`overlap_summary*.json`):

| Check | Result |
|---|---|
| Pre-declared rule E: raw Pearson r > 0.99 and mutual best match (28,869 shared probes) | 21 / 90 arrays pass, of which 13 / 45 tumours. All 90 best matches are mutual. Best-match r ranges from 0.984 to 0.991 (median 0.990), because GSE62452 was re-normalised and batch-corrected in Partek. |
| Rule T: the best match carries the same running sample number (CEL file name `<n>_...`) | 90 / 90 |
| Post-hoc diagnostic: each probe centred within its own series before correlating (this removes the shared tissue profile) | All 90 arrays are mutual best matches with the same sample number. Matched pairs have r 0.983–0.997 (median 0.992); every non-matched pair has r ≤ 0.833 (99.9th percentile 0.778). |
| Rule S: same survival months (±1) and event status | 38 / 45 tumours. Of the other 7: 3 have no survival data in either series, and 4 have updated or conflicting follow-up (see below). |

**Conclusion: all 45 GSE28735 tumours are inside GSE62452.** Of these, 42 have survival data in GSE28735 and 41 have
it in GSE62452, so 41 patients appear in both Phase B survival analyses. They make up 41 of the 65 GSE62452 patients
(63 %) and 41 of the 42 GSE28735 patients (98 %). The pre-declared r > 0.99 threshold alone confirms only 13 tumours.
I report that result as it came out; the threshold was set too high for arrays that had been re-normalised. The order
and survival checks, and the post-hoc centred correlation, are unambiguous.

Survival follow-up differs for 4 shared tumours (GSE28735 → GSE62452): sample 33 had 38 → 40 months (still
censored); sample 35 had 1 month censored → 1.2 months with an event; sample 49 had 14 months with an event → missing;
sample 55 had 24 months alive → 9.8 months dead. The last of these is a genuine inconsistency between the two
records. It cannot be resolved from GEO.

## 2. Which v1 results are invalid

Because the same patients were counted twice, these results from HYPOXIA_VALIDATION_v1.md are **invalid** and must
not be cited:

1. **The two-cohort "independent" random-effects meta-analysis**: HR 1.383 [1.067, 1.794], one-sided p = 0.0072
   (`meta_analysis_primary.csv`). Its inputs are not independent.
2. **Its heterogeneity**: Q = 0.532, I² = 0 %. The two estimates agree largely because they come from the same
   patients.
3. **The power calculation**: "78 pooled independent events" and minimum detectable HR 1.325
   (`power_and_detectable_hr.csv`). At least 28 of the 29 GSE28735 events are the same patients as GSE62452 events.
   The unique count is 49 + 1 = 50 events, not 78.
4. **The 4-cohort secondary pooled result**: HR 1.367 [1.191, 1.570], p = 4.4×10⁻⁶, I² = 0 %
   (`meta_analysis_secondary.csv`). It double-counts the same patients and also includes the two discovery cohorts.
5. **The leave-one-out table** (`leave_one_out.csv`) and the statement "107 primary PDAC tumours and 78 overall
   survival events" in two "fully independent" cohorts. The true figure is 66 unique patients with survival data.
6. **The exploratory pooled HRs** for the 8 secondary signatures (`exploratory_signatures.csv`). They were pooled in
   the same way, so they are also not two independent replications. Their null conclusion is unchanged, but their
   intervals are too narrow.

**What remains valid:** **GSE62452 alone**, stage-adjusted Cox: **HR 1.48 [1.08, 2.03] per SD**, β = 0.3921
(SE 0.1619), one-sided p = 0.0077, n = 65, 49 events (`hypoxia_models.csv`). GSE62452 contains GSE28735, so GSE28735
adds no independent evidence. Until Section 3, this was the **only independent evidence**. It is a post-hoc
re-labelling of a result that was already seen.

## 3. New pre-registered confirmatory test: GSE57495

The plan was committed in d897e19 before any GSE57495 expression value was scored.

- **Cohort:** GSE57495 (Moffitt Cancer Center; contributors Chen and Malafa; PMID 26247463). 63 microdissected frozen
  PDAC tissues; platform GPL15048 (Rosetta/Merck Human RSTA Custom Affymetrix 2.0). All 63 have overall survival;
  42 events. The only clinical covariate in the series is stage (2B: 33, 2A: 17, 1B: 12, "1": 1).
- **Probe annotation:** taken from the GPL15048 platform table on GEO (column `GeneSymbol`, NCBI-generated; 60,607
  probes, 41,024 with a symbol). A probe counts only if it maps to exactly one gene symbol. Coverage is **12 / 12**
  (100 %), so the ≥ 80 % gate passed. Where a gene had several probes, the one with the highest mean expression was
  used (`gse57495_gpl15048_probe_map.csv`). Values are already on log2 scale (range 1.46–15.06), so no transform was
  applied.
- **Overlap with other cohorts used** (`gse57495_overlap_check.json`):
  - No shared GSM accessions, titles or CEL file names with GSE62452 or GSE28735, and the series come from a
    different institution.
  - Cross-platform Spearman screen against TCGA PAAD (177 tumours, 18,060 shared genes) and CPTAC PDA (140 tumours,
    17,752 genes): no pair exceeded the pre-declared flag (r > 0.8 and above the 99.9th percentile). Maximum r was
    0.778 for TCGA and 0.522 for CPTAC.
  - The TCGA maximum belongs to one GSE57495 sample (PDAC19). That sample also correlates 0.63, 0.60 and 0.57 with
    three other TCGA tumours, which fits a shared expression profile better than a duplicated patient. TCGA
    tissue-source-site institutions were not checked: **UNVERIFIED**.
  - A post-hoc sensitivity analysis excluding PDAC19 is reported below.

| Model (`gse57495_cox_results.csv`) | n | events | HR per SD [95 % CI] | one-sided p | PH test p |
|---|---|---|---|---|---|
| **Stage-adjusted (PRIMARY)** | 63 | 42 | **1.347 [0.994, 1.823]** | **0.027** | 0.80 |
| Unadjusted | 63 | 42 | 1.280 [0.941, 1.743] | 0.058 | 0.51 |
| Stage-adjusted, mean of all probes per gene (sensitivity) | 63 | 42 | 1.365 [1.002, 1.859] | 0.024 | 0.98 |
| Stage-adjusted, excluding PDAC19 (post-hoc sensitivity) | 62 | 42 | 1.311 [0.961, 1.788] | 0.044 | – |

**Pre-declared decision: CONFIRMED in an independent cohort** (one-sided p = 0.027 < 0.05, HR > 1).

Honest caveats:
- The two-sided 95 % CI just includes 1, and the two-sided p is 0.054. The test passes only because it was declared
  one-sided in advance.
- With 42 events, pre-declared power was 0.56 at HR 1.32 and 0.81 at HR 1.48.
- The unadjusted model does not reach significance (one-sided p = 0.058).

**Secondary, descriptive and post-hoc** (GSE62452 was already seen): fixed-effect pooling of the two cohorts that
share no patients, GSE62452 (stage-adjusted) and GSE57495 (primary). Pooled HR is **1.41 [1.13, 1.75]**, one-sided
p = 0.0011, Q = 0.18, I² = 0 %. DerSimonian-Laird gives the same values (τ² = 0). No decision is taken from this.
With two cohorts, I² carries almost no information.

## 4. This is a confirmation, not a discovery

Hypoxia expression signatures have been reported as prognostic in PDAC before. Each source below was opened:

- Abou Khouzam R, …, Chouaib S. *An Eight-Gene Hypoxia Signature Predicts Survival in Pancreatic Cancer and Is
  Associated With an Immunosuppressed Tumor Microenvironment.* Front Immunol 2021;12:680435,
  doi:10.3389/fimmu.2021.680435. Source: publisher page
  https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.680435/full, fetched 2026-11-17.
  - The abstract says the signature "was significantly associated with overall survival in two distinct PDA
    datasets" (TCGA PAAD and E-MTAB-6134) and "showed independent prognostic value in multivariate analysis".
  - 6 of its 8 genes (VEGFA, NDRG1, SLC2A1, P4HA1, LDHA, PGK1) are in our 12-gene panel.
  - It already used TCGA PAAD, which is one of our discovery cohorts.
- Ren M, Feng L, Zong R, Sun H. *Novel prognostic gene signature for pancreatic ductal adenocarcinoma based on
  hypoxia.* World J Surg Oncol 2023, doi:10.1186/s12957-023-03142-2, PMID 37605192. Abstract fetched from the Europe
  PMC REST API on 2026-11-17. It describes a TCGA-derived hypoxia risk model validated in ICGC PDAC, which it reports
  as "an independent factor for prognosis prediction".
- Connor AA, et al. *Integration of Genomic and Transcriptional Features in Pancreatic Cancer Reveals Increased Cell
  Cycle Progression in Metastases.* Cancer Cell 2019, doi:10.1016/j.ccell.2018.12.010, PMID 30686769. Abstract
  fetched from Europe PMC. It reports "Half of the cases were hypoxic by expression markers" and "Established clinical
  co-variates dominated survival analyses, although CCP and hypoxia may inform clinical practice".

Some sources could not be opened:
- The PMC copies of these articles returned a reCAPTCHA page, and DOAJ and HAL returned 403 / bot-check pages. I did
  not retry or work around them.
- Other papers cited in Abou Khouzam et al. (e.g. Miller et al. 2015, EMBO Mol Med) were not opened: **UNVERIFIED**.

The thesis therefore **confirms**, in one independent cohort with pre-registration, a known association between
tumour hypoxia expression and worse overall survival in PDAC. It is not a new prognostic discovery. It is an
expression association only, with no causal or clinical-utility claim, and the patient is the unit of inference.

## 5. Corrected status (one line)

**Hypoxia → worse OS:**
- Discovery cohorts: TCGA and CPTAC (post-hoc selected).
- Independent evidence: GSE62452 alone, HR 1.48 [1.08, 2.03] (post-hoc re-labelled).
- Pre-registered confirmation: GSE57495, HR 1.35 [0.99, 1.82], one-sided p = 0.027. This passes the pre-declared
  one-sided test, but the two-sided CI touches 1.
- GSE28735 is a subset of GSE62452 and adds nothing.

# Phase B Plan: Independent Validation of the Hypoxia Survival Lead

## 1. Context and Hypothesis
In discovery analyses across TCGA (strict PDAC, 141 patients, 81 events) and CPTAC (129 patients, 72 events), a 12-gene hypoxia signature emerged as the strongest prognostic lead:
- TCGA: adjusted HR 1.32 [1.07, 1.62] per SD ($p = 0.009$, BH $q = 0.078$).
- CPTAC: adjusted HR 1.43 [1.10, 1.86] per SD ($p = 0.007$).

Because hypoxia was selected post-hoc following inspection of these cohorts, TCGA and CPTAC cannot serve as their own validation. Phase B conducts an independent validation on public PDAC cohorts that were **not** used in model development.

### Single Primary Hypothesis
- Higher tumour expression of the hypoxia signature (exact v2 12-gene panel) associates with worse overall survival in independent PDAC cohorts:
  $$H_0: \beta_{\text{hypoxia}} \le 0 \quad \text{vs.} \quad H_1: \beta_{\text{hypoxia}} > 0 \quad (\alpha = 0.05, \text{one-sided})$$

---

## 2. Phase B1: Cohort Candidate Vetting
Candidates to evaluate from their official public portals/repositories without bypassing access controls:
1. `ICGC PACA-AU` (Australia, Bailey et al. Nature 2016)
2. `ICGC PACA-CA` (Canada)
3. `GSE71729` (Moffitt et al. Nat Genet 2015)
4. `GSE62452` (Yang et al. Cancer Cell 2016)
5. `GSE57495` (Chen et al. Clin Cancer Res 2015)
6. `GSE28735` (Zhang et al. PLoS ONE 2012)
7. `E-MTAB-6134` (ArrayExpress PDAC cohort)

### Deliverable: `results/translational/phaseB/cohort_vetting.csv`
Columns: `cohort_id`, `platform`, `n_total`, `n_with_os`, `n_events`, `access_terms`, `citation`, `status` (`OBTAINED`, `BLOCKED`, `INSUFFICIENT_CLINICAL`), `reason_if_blocked`, `source_url`.

---

## 3. Phase B2: Statistical Protocol & Primary Meta-Analysis

### 3.1 Gene Signature Scoring
- **Gene list (12 genes):**
  `VEGFA`, `SLC2A1`, `CA9`, `LDHA`, `PGK1`, `ADM`, `NDRG1`, `BNIP3`, `ANKRD37`, `P4HA1`, `ENO1`, `HK2`.
- **Coverage threshold:** Cohort must measure $\ge 80\%$ of genes (at least 10 of 12 genes).
- **Standardisation:** Genes are z-scored across patients within cohort, then averaged to form the raw signature score, and finally standard-normalised within the cohort ($z \sim \mathcal{N}(0, 1)$).

### 3.2 Survival Modeling per Cohort
- Cox proportional hazards regression:
  $$\lambda(t) = \lambda_0(t) \exp(\beta_{\text{hypoxia}} \cdot Z_{\text{hypoxia}} + \boldsymbol{\gamma}^T \mathbf{X})$$
- Adjustment covariates: Adjust for whatever of `age`, `stage`, and `tumour_purity` each cohort provides (declared per cohort).
- Extract log hazard ratio $\hat{\beta}_i$, standard error $\text{SE}(\hat{\beta}_i)$, and hazard ratio $\text{HR}_i = \exp(\hat{\beta}_i)$ per SD.

### 3.3 Primary Meta-Analysis
- DerSimonian-Laird random-effects meta-analysis across NEW independent cohorts only:
  $$w_i = \frac{1}{\text{SE}_i^2 + \tau^2}, \quad \hat{\beta}_{\text{pooled}} = \frac{\sum w_i \hat{\beta}_i}{\sum w_i}$$
- One-sided p-value against $H_0: \beta \le 0$:
  $$p_{\text{one-sided}} = 1 - \Phi\left(\frac{\hat{\beta}_{\text{pooled}}}{\text{SE}_{\text{pooled}}}\right)$$
- Heterogeneity:
  - Cochran's $Q = \sum \frac{1}{\text{SE}_i^2} (\hat{\beta}_i - \hat{\beta}_{\text{fixed}})^2$
  - Higgins & Thompson $I^2 = \max\left(0, \frac{Q - (k - 1)}{Q}\right) \times 100\%$
- Leave-one-cohort-out sensitivity analysis.

### 3.4 Secondary Meta-Analysis
- Re-run including discovery cohorts TCGA ($\hat{\beta} = 0.278$, $\text{SE} = 0.106$) and CPTAC ($\hat{\beta} = 0.358$, $\text{SE} = 0.134$), explicitly flagged as `DISCOVERY (NOT INDEPENDENT)`.

---

## 4. Phase B3: Statistical Power & Exploratory Signatures
- **Pooled Events:** Sum of overall survival events $E = \sum E_i$.
- **Minimum Detectable HR:** At 80% power ($\beta_{\text{type2}} = 0.20$) and one-sided $\alpha = 0.05$:
  $$\text{MDHR} = \exp\left( \frac{1.645 + 0.842}{\sqrt{E}} \right)$$
- **Exploratory Signatures:**
  - 8 additional v2 signatures: `CAF`, `EMT`, `Interferon`, `Immune_infiltration`, `T_cell_proxy`, `Macrophage_C1Q_TAM`, `Macrophage_SPP1_TAM`, `C1Q_minus_SPP1`.
  - FDR corrected with Benjamini-Hochberg across the 8 exploratory tests.
  - Reported separately in exploratory tables; cannot alter the primary hypothesis verdict.

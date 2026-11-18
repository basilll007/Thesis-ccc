# Spatial Communication Audit Report: Mendeley PDAC-SMI Full U5-a Headline Demonstration

**Date:** 2026-11-17 18:42:42
**Execution Runtime:** 3.09 seconds
**Dataset:** `F:/Thesis/data/processed/mendeley_u5a_full_contract.h5ad` (74793 cells, 960 genes)

## Executive Summary & Verdicts

| ID | Target Interaction | Verdict | Primary Reason |
|---|---|---|---|
| CLAIM_001 | `CAF->CD8+ T (CXCL12-CXCR4)` | **EXPLAINED_BY_ABUNDANCE** | Matched dyadic baseline M is significant (p=0.0099), but proximity ratio Q is not (p=0.0594). Score is driven by cell abundance. |
| CLAIM_002 | `CAF->Macrophage (COL1A1-CD44)` | **EXPLAINED_BY_ABUNDANCE** | Matched dyadic baseline M is significant (p=0.0099), but proximity ratio Q is not (p=0.1881). Score is driven by cell abundance. |
| CLAIM_003 | `CAF->Macrophage (FN1-CD44)` | **SUPPORTED** | Passes all confounding controls: contact enrichment persists under spatially constrained nulls and outside spillover range. |
| CLAIM_004 | `Macrophage->Malignant (SPP1-CD44)` | **SUPPORTED** | Passes all confounding controls: contact enrichment persists under spatially constrained nulls and outside spillover range. |

## Diagnostic Details
- **Graph connectivity:** 6-NN with 448758 directed edges.
- **Null model calibrations:** 100 permutations; block tile size = 100.0 um.
- **Segmentation spillover excess:** 0.0% (explained by spillover: False).

## Generated Deliverables
- Detailed global scores: `framework\results\demo_u5a_full\scores_global_null.csv`
- Detailed block scores: `framework\results\demo_u5a_full\scores_block_null.csv`
- Claims register table: `framework\results\demo_u5a_full\claims_register.csv`
- Diagnostic figure: `framework\results\demo_u5a_full\figures\audit_summary_figure.png`
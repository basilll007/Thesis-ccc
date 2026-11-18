# Spatial Communication Audit Report: Mendeley PDAC-SMI U5-a Headline Demonstration

**Date:** 2026-11-17 18:42:16
**Execution Runtime:** 1.53 seconds
**Dataset:** `F:/Thesis/data/processed/mendeley_u5a_contract.h5ad` (8000 cells, 960 genes)

## Executive Summary & Verdicts

| ID | Target Interaction | Verdict | Primary Reason |
|---|---|---|---|
| CLAIM_001 | `CAF->CD8+ T (CXCL12-CXCR4)` | **NOT_TESTABLE** | Receiver cell count (11) is below power floor (25). |
| CLAIM_002 | `CAF->Macrophage (COL1A1-CD44)` | **EXPLAINED_BY_SPILLOVER** | Explained by segmentation spillover (marker excess 99.8% in touching cells; signal confined to touching distance). |
| CLAIM_003 | `CAF->Macrophage (FN1-CD44)` | **EXPLAINED_BY_SPILLOVER** | Explained by segmentation spillover (marker excess 99.8% in touching cells; signal confined to touching distance). |
| CLAIM_004 | `Macrophage->Malignant (SPP1-CD44)` | **EXPLAINED_BY_SPILLOVER** | Explained by segmentation spillover (marker excess 99.8% in touching cells; signal confined to touching distance). |

## Diagnostic Details
- **Graph connectivity:** 6-NN with 48000 directed edges.
- **Null model calibrations:** 500 permutations; block tile size = 100.0 um.
- **Segmentation spillover excess:** 99.8% (explained by spillover: True).

## Generated Deliverables
- Detailed global scores: `framework\results\demo_u5a\scores_global_null.csv`
- Detailed block scores: `framework\results\demo_u5a\scores_block_null.csv`
- Claims register table: `framework\results\demo_u5a\claims_register.csv`
- Diagnostic figure: `framework\results\demo_u5a\figures\audit_summary_figure.png`
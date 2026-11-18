# Spatial Communication Audit Report: Mendeley PDAC-SMI Full U5-a Headline Demonstration (v2 engine)

**Date:** 2026-11-17 22:38:52
**Execution Runtime:** 46.27 seconds
**Dataset:** `F:/Thesis/data/processed/mendeley_u5a_full_contract.h5ad` (74793 cells, 960 genes)

## Executive Summary & Verdicts

Every check is computed for every claim (PASS / FAIL / NOT_RUN). The primary verdict follows the
pre-declared priority rule: 1 power FAIL -> NOT_TESTABLE; 2 spillover FAIL -> EXPLAINED_BY_SPILLOVER; 3 abundance FAIL -> EXPLAINED_BY_ABUNDANCE if naive M (label null) significant else NOT_SUPPORTED; 4 zonation FAIL -> EXPLAINED_BY_ZONATION; 5 zonation or spillover NOT_RUN -> NOT_TESTABLE; 6 SUPPORTED.

| ID | Claim | power | abundance | zonation | spillover | Primary verdict | Failing checks |
|---|---|---|---|---|---|---|---|
| CLAIM_001 | `CAF->CD8+ T (CXCL12-CXCR4)` | FAIL | PASS | FAIL | PASS | **NOT_TESTABLE** | power;zonation |
| CLAIM_002 | `CAF->Macrophage (COL1A1-CD44)` | PASS | FAIL | FAIL | PASS | **EXPLAINED_BY_ABUNDANCE** | abundance;zonation |
| CLAIM_003 | `CAF->Macrophage (FN1-CD44)` | PASS | FAIL | FAIL | PASS | **EXPLAINED_BY_ABUNDANCE** | abundance;zonation |
| CLAIM_004 | `Macrophage->Malignant (SPP1-CD44)` | PASS | PASS | FAIL | PASS | **EXPLAINED_BY_ZONATION** | zonation |

### Check reasons

- `CAF->CD8+ T (CXCL12-CXCR4)`: power: n_receivers=476 (floor 25), edges=986 (floor 10), chain power at 1.5x=0.19 (floor 0.5); abundance: contact test vs within-type permutation p=0.0010; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.2138; ring-control p=0.3616; [info: block-null Q p=0.2657]; spillover: no receptor excess in touching receivers (phi=0.3220, p=0.0010, dR_obs=-4.237, dR_exp=-10.886)
- `CAF->Macrophage (COL1A1-CD44)`: power: n_receivers=7014 (floor 25), edges=12713 (floor 10), chain power at 1.5x=0.91 (floor 0.5); abundance: contact test vs within-type permutation p=1.0000; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.8402; ring-control p=0.8551; [info: block-null Q p=0.9101]; spillover: no receptor excess in touching receivers (phi=0.1275, p=0.0010, dR_obs=-3.293, dR_exp=-1.552)
- `CAF->Macrophage (FN1-CD44)`: power: n_receivers=7014 (floor 25), edges=12713 (floor 10), chain power at 1.5x=0.90 (floor 0.5); abundance: contact test vs within-type permutation p=1.0000; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.8591; ring-control p=0.8571; [info: block-null Q p=0.5914]; spillover: no receptor excess in touching receivers (phi=0.1275, p=0.0010, dR_obs=-3.293, dR_exp=-1.552)
- `Macrophage->Malignant (SPP1-CD44)`: power: n_receivers=41216 (floor 25), edges=14726 (floor 10), chain power at 1.5x=1.00 (floor 0.5); abundance: contact test vs within-type permutation p=0.0010; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.9870; ring-control p=0.9890; [info: block-null Q p=0.7453]; spillover: no receptor excess in touching receivers (phi=0.0830, p=0.0010, dR_obs=-0.321, dR_exp=-1.461)

## Diagnostic Details
- **Graph connectivity:** 6-NN with 448758 directed edges.
- **Permutations:** 1000 per null; minimum attainable p = 1/1001 = 0.000999; block tile = 100.0 um.
- **Spillover:** computed separately for each claim's sender/receiver pair on all 960 genes (see claims_register.csv).

## Generated Deliverables
- Detailed global scores: `framework\results\demo_u5a_full_v2\scores_global_null.csv`
- Detailed block scores: `framework\results\demo_u5a_full_v2\scores_block_null.csv`
- Claims register table (multi-label): `framework\results\demo_u5a_full_v2\claims_register.csv`
- Per-claim metrics: `framework\results\demo_u5a_full_v2\claim_metrics.csv`
- Diagnostic figure: `framework\results\demo_u5a_full_v2\figures\audit_summary_figure.png`
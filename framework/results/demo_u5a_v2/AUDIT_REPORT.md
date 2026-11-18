# Spatial Communication Audit Report: Mendeley PDAC-SMI U5-a Headline Demonstration (v2 engine)

**Date:** 2026-11-17 22:32:49
**Execution Runtime:** 5.20 seconds
**Dataset:** `F:/Thesis/data/processed/mendeley_u5a_contract.h5ad` (8000 cells, 960 genes)

## Executive Summary & Verdicts

Every check is computed for every claim (PASS / FAIL / NOT_RUN). The primary verdict follows the
pre-declared priority rule: 1 power FAIL -> NOT_TESTABLE; 2 spillover FAIL -> EXPLAINED_BY_SPILLOVER; 3 abundance FAIL -> EXPLAINED_BY_ABUNDANCE if naive M (label null) significant else NOT_SUPPORTED; 4 zonation FAIL -> EXPLAINED_BY_ZONATION; 5 zonation or spillover NOT_RUN -> NOT_TESTABLE; 6 SUPPORTED.

| ID | Claim | power | abundance | zonation | spillover | Primary verdict | Failing checks |
|---|---|---|---|---|---|---|---|
| CLAIM_001 | `CAF->CD8+ T (CXCL12-CXCR4)` | FAIL | FAIL | NOT_RUN | NOT_RUN | **NOT_TESTABLE** | power;abundance |
| CLAIM_002 | `CAF->Macrophage (COL1A1-CD44)` | FAIL | FAIL | FAIL | PASS | **NOT_TESTABLE** | power;abundance;zonation |
| CLAIM_003 | `CAF->Macrophage (FN1-CD44)` | FAIL | FAIL | FAIL | PASS | **NOT_TESTABLE** | power;abundance;zonation |
| CLAIM_004 | `Macrophage->Malignant (SPP1-CD44)` | PASS | PASS | FAIL | PASS | **EXPLAINED_BY_ZONATION** | zonation |

### Check reasons

- `CAF->CD8+ T (CXCL12-CXCR4)`: power: n_receivers=11 (floor 25), edges=19 (floor 10), chain power at 1.5x=0.05 (floor 0.5); abundance: contact test vs within-type permutation p=1.0000; naive M label-null p=0.0010; zonation: local controls not runnable (exposed=4, ring controls=3); spillover: too few isolated receivers (3)
- `CAF->Macrophage (COL1A1-CD44)`: power: n_receivers=318 (floor 25), edges=450 (floor 10), chain power at 1.5x=0.18 (floor 0.5); abundance: contact test vs within-type permutation p=0.9401; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.5185; ring-control p=0.2657; [info: block-null Q p=0.7882]; spillover: contamination explains < 50% of the receptor excess (phi=0.1496, p=0.0010, dR_obs=0.760, dR_exp=-0.904)
- `CAF->Macrophage (FN1-CD44)`: power: n_receivers=318 (floor 25), edges=450 (floor 10), chain power at 1.5x=0.19 (floor 0.5); abundance: contact test vs within-type permutation p=0.9680; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.4316; ring-control p=0.2288; [info: block-null Q p=0.7542]; spillover: contamination explains < 50% of the receptor excess (phi=0.1496, p=0.0010, dR_obs=0.760, dR_exp=-0.904)
- `Macrophage->Malignant (SPP1-CD44)`: power: n_receivers=5190 (floor 25), edges=1074 (floor 10), chain power at 1.5x=0.95 (floor 0.5); abundance: contact test vs within-type permutation p=0.0010; naive M label-null p=0.0010; zonation: within-tile (100 um) receiver contrast p=0.9690; ring-control p=0.9910; [info: block-null Q p=0.8102]; spillover: contamination explains < 50% of the receptor excess (phi=0.0689, p=0.0010, dR_obs=0.121, dR_exp=-1.737)

## Diagnostic Details
- **Graph connectivity:** 6-NN with 48000 directed edges.
- **Permutations:** 1000 per null; minimum attainable p = 1/1001 = 0.000999; block tile = 100.0 um.
- **Spillover:** computed separately for each claim's sender/receiver pair on all 960 genes (see claims_register.csv).

## Generated Deliverables
- Detailed global scores: `framework\results\demo_u5a_v2\scores_global_null.csv`
- Detailed block scores: `framework\results\demo_u5a_v2\scores_block_null.csv`
- Claims register table (multi-label): `framework\results\demo_u5a_v2\claims_register.csv`
- Per-claim metrics: `framework\results\demo_u5a_v2\claim_metrics.csv`
- Diagnostic figure: `framework\results\demo_u5a_v2\figures\audit_summary_figure.png`
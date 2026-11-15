# v3 Phase 2 plan — calibration by simulation (is the null result low power?)

Written 2026-11-15 BEFORE any Phase 2 result was computed. Branch `thesis-v3`.

## Question
Can the audit scores (C, M, E, Q) under the global and Phase 1 constrained nulls detect a *planted* proximity-specific
ligand–receptor effect of a given size, and which score/null combination is fooled by zonation alone?
The planted perturbation is a computational stand-in for a stimulation experiment; nothing experimental is claimed.

## Sections (fixed rule, from `results/pdac/phase3/cross_section_section_info.csv`, n_cells_after_qc)
smallest = T6-b (12,473 cells); lower median (rank 8 of 16) = T2-b (43,658); largest = U5-a (74,793).
Real positions, real 6-NN graph, real `consensus_T` labels, real counts are kept; only CXCL12 / CXCR4 counts are perturbed.

## Axis
CXCL12 -> CXCR4, sender CAF_fibroblast -> receiver T_cell (the headline). Scores computed exactly as in `ccc_audit.py`
(C = (L̄s + R̄r)/2, M = L̄s·R̄r, E = mean over CAF->T edges of L_i·R_j, Q = E/M) on log1p(normalize_total 1e4) expression
recomputed after each perturbation of the counts.

## Ground truths (perturbations on raw counts, then renormalise + log1p)
(i) PROXIMITY: target receivers = T cells that are the `col` of a 6-NN edge whose `row` is a CAF with CXCL12 count > 0.
    CXCR4 count of each target j is redrawn as Poisson(f · c_j). Expected target-set total is f × real total (fold change f).
    Non-target cells are untouched. f ∈ {1.0 (control), 1.1, 1.25, 1.5, 2.0}.
(ii) ZONATION-ONLY: square 200 um tiles with a random offset per simulation (uniform in [0,200) per axis);
    each tile is "high" with probability 0.5. In high tiles, CXCL12 and CXCR4 counts of ALL cells (any type) are redrawn as
    Poisson(f · c). No dependence on adjacency or on cell type. f ∈ {1.0 (control, identical procedure to (i) control), 1.1, 1.25, 1.5, 2.0}.
Simulation s uses numpy default_rng(42 + s) for the Poisson draws and tile offsets; s = 0..99 (100 simulations per setting).
If memory-bound, the number of simulations is reduced and logged.

## Nulls
global, block200, block100, block50, dist — exactly `v3_nulls.make_perms`, seed 42, 500 permutations, built once per
section x null (labels and positions do not change between simulations) and re-used for every simulation.
Implementation: per permutation we store the indices of CAF->T edges and of CAF / T cells; E, M, C, Q are then evaluated
by bincount for all 500 permutations at once. This is algebraically identical to `ccc_audit.run_null` restricted to one row.

## Statistics
z = (obs − null mean)/null sd (ddof 1); add-one p = (k+1)/(n+1), one-sided (null ≥ obs). Detection = add-one p < 0.05.
- Power (scenario i): fraction of the 100 simulations with p < 0.05, per score x null x f x section.
- False-positive rate (scenario ii): same fraction, per score x null x f x section.
- The f = 1.0 control gives the baseline rate of each score/null on the real (unperturbed-in-expectation) tissue;
  the real section may itself carry signal, so rates are reported both raw and alongside the control.
- Also reported: mean z, mean observed Q, number of target cells, number of CAF->T edges.

## Outputs (results/pdac/v3_phase2/)
- `sim_rows.csv.gz`: one row per section x scenario x f x simulation x null x score (obs, null mean/sd, z, p).
- `calibration_table.csv`: power / FPR per section x scenario x f x null x score (+ pooled over the 3 sections).
- `calibration_summary.json`: smallest f detected with power >= 0.8 per score x null (pooled), FPR at f = 2 zonation.
- Figure `figures/pdac/v3_phase2/phase2_calibration.{png,pdf}`: power vs f (proximity) and FPR vs f (zonation), one line per null, one panel per score.

## Interpretation rules (declared now)
- "Can detect" = pooled power >= 0.8 at that f. "Fooled by zonation" = pooled FPR > 0.10 at f = 1.5 (zonation).
- Audit/diagnostic framing only. Every number from a registered file.

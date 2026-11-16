# Chain Phase 4 plan: calibration by simulation (required before interpreting Phase 3)

Committed before any Phase 2, 3 or 4 result was computed. Seed 42.

## Sections
Three real sections that pass the primary-chain floor (tumour -> CAF -> T, author_mapped) in Phase 1: those with the largest, the median and the smallest number of T cells among the passing sections. If two sections tie, the one with the lower alphabetical name is taken. The point is to span the cell-count range, so that a minimum trustworthy cell count can be stated.

## Planted programs (fixed before Phase 3)
A = tumour:basal_like, B = CAF:myCAF, C = T:cytotoxic. (The first draft of this plan, before commit, named CAF:iCAF; iCAF was dropped in Phase 1/1b with only 2 usable genes, so myCAF, the only multi-gene CAF program kept, is used. This was fixed before any Phase 2-4 outcome.) Programs and controls come from Phase 1 Amendment 1 (results/chain/phase1b, data/chain_cache/mendeley_v1b). Planting acts on the raw counts of each program's own genes, in cells of the program's own type.

## Count model
- For every cell of type A, B or C, the counts of that type's planted-program genes are redrawn as Poisson(f_i * c), where f_i = 1 for unperturbed cells. The extra Poisson noise is therefore the same in every scenario.
- Then the per-cell normalize_total(1e4) and log1p are recomputed, and the cell's planted-program score is recomputed with the Phase 1 recorded controls (control genes are untouched).
- "A-high" means an A cell whose basal_like score is >= the 75th percentile of A cells in the section (real data).

## Scenarios (100 simulations each, per section)
- **(0) Null:** f = 1 everywhere, i.e. Poisson resampling only.
- **(i) True chain**, per-hop effect f in {1.1, 1.25, 1.5, 2.0}:
  - B cells with >= 1 A-high neighbour get f on their B-program genes.
  - Then C cells with >= 1 *perturbed* B neighbour get f on their C-program genes.
- **(ii) Zonation only**, f in {1.5, 2.0}:
  - The field is z(x, y) = sin(2 pi x / 500 um + phi1) + sin(2 pi y / 500 um + phi2), with phases redrawn in each simulation.
  - Every A, B and C cell in the region z > 0 gets f on its program genes.
  - There is no adjacency dependence.
- **(iii) Direct A -> C, no B involvement**, f in {1.1, 1.25, 1.5, 2.0}: C cells with >= 1 A-high neighbour get f; B is untouched.

## Test applied
This is the Phase 3 section-level test for the planted triple only, with the same exposures, controls and expression block null. Permutations are reduced to 100 per simulation for compute (logged). The section-level chain call mirrors decision-rule criteria 1-3, without BH, because there is one test per simulation:
- hop C<-B p < 0.05;
- hop B<-A p < 0.05;
- indirect p < 0.05;
- |beta_direct| < |beta_total|.

Two further checks:
- **Label block null:** calibrated on a reduced set of 20 simulations per setting, for the null, chain f = 1.5, zonation f = 2.0 and direct f = 1.5 (logged reduction). Criterion: hop p < 0.05 under the label null.
- **Pooled 3-section call:** simulations with the same index are combined across the 3 sections in the Phase 3 tumour-level way (mean of section statistics, with a paired null replicate). This gives a multi-section power estimate.

## Reported
- Power to call the chain: the fraction of scenario (i) runs called, per f and section.
- False-positive rate: the fraction of (0) and (ii) runs called.
- Wrongly called chain under (iii): the fraction of (iii) runs called a chain.
- Per-hop detection rates: hop C<-B and hop B<-A separately.
- Mean estimates of beta_CB, beta_BA and indirect, per scenario.

## Calibrated range (decision rule)
The test is "trustworthy" at (f, section cell counts) when all three hold:
1. Power >= 0.8 in (i).
2. Call rate <= 0.10 in (0) and in zonation (ii) at f = 2.0.
3. Call rate <= 0.10 in (iii) at the same f.

The smallest f and the smallest section (by number of analysed C cells) meeting these define the calibrated range. Phase 3 is interpreted only within that range, and the summary says so. If no setting qualifies, every Phase 3 chain call is "outside calibrated range".

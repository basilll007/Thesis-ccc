# Chain Phase 3 plan: the chain test (core)

Committed after Phase 1 (counts and graphs only) and before any Phase 2 or Phase 3 outcome was computed. Seed 42. Primary annotation `author_mapped`; `marker` is a sensitivity analysis, computed after the primary and labelled as such. Wording: associations and ranked hypotheses; expression, not mutation; no causal or mechanistic claims.

## Chains (pre-declared, outcome-free)
- **Primary:** tumour (A) -> CAF (B) -> T cell (C).
- **Secondary:** up to 5 chains, chosen by the Phase 1 rule (cell and edge counts only): results/chain/phase1/chain_selection.csv, rows with role == "secondary". The realised list, from results/chain/phase1b/chain_selection.csv (identical to phase1/, since selection uses counts only), in rank order:
  1. myeloid -> CAF -> tumour (16 sections, 13 tumours pass)
  2. tumour -> CAF -> myeloid (16, 13)
  3. CAF -> myeloid -> tumour (15, 12)
  4. CAF -> tumour -> myeloid (15, 12)
  5. myeloid -> tumour -> CAF (15, 12)
  The primary chain passes in 14 sections (11 tumours). No secondary chain contains T cells: T-cell counts and edges are lower, so T-containing triples rank lower. Pairs 1/2 and 3/5 are mutual reversals; each also serves as the other's reversed control.
- Programs (Phase 1 Amendment 1, results/chain/phase1b/programs_used.json): tumour {classical, basal_like}; CAF {myCAF, CXCL12 readout}; T {cytotoxic}; myeloid {C1Q_TAM, SPP1_TAM}. That is 4 combinations for the primary chain and 8 for each secondary chain, 44 in total. Under the original section-wide rule, the primary chain is not testable (no usable T-cell program).
- A section enters a chain's test only if it passes that chain's Phase 1 floor (>= 100 cells of each type, >= 300 A-B and >= 300 B-C edges, symmetrised 6-NN).
- Programs: every combination (A program x B program x C program) of the Phase 1 kept programs. The CAF CXCL12 readout counts as a CAF program.

## Exposures (graph W = symmetrised 6-NN, um coordinates)
- For each C cell c, the **B-exposure** e_B(c) is the mean B-program score over B cells adjacent to c. It is defined when c has >= 1 B neighbour.
- The **2-hop A-exposure** e_A2(c) is the path-weighted mean A-program score over paths c - b - a, where b is a B neighbour of c, a is an A neighbour of b, a is not adjacent to c, and a != c. An A cell reached through two B cells counts twice.
- For each B cell b, the **1-hop A-exposure** e_A1(b) is the mean A-program score over A neighbours of b.

## Tests (OLS within section)
Variables are z-scored within section over the analysis rows, so all coefficients are standardised.
- **Controls** (the Phase 2 zonation controls, extended to the exposure side):
  - log1p total counts of the outcome cell.
  - Local density: cells within 30 um.
  - Leave-one-out 200 um tile mean of the outcome program over same-type cells.
  - The 200 um tile mean of each exposure program over its own type's cells in the same tile.
- **(1) Hop C <- B:** among C cells with e_B defined, y_C ~ e_B + controls. The coefficient is beta_CB.
- **(2) Hop B <- A:** among B cells with e_A1 defined, y_B ~ e_A1 + controls. The coefficient is beta_BA.
- **(3) Chain:** among C cells with both e_B and e_A2 defined:
  - total: y_C ~ e_A2 + controls, giving beta_total;
  - direct: y_C ~ e_A2 + e_B + controls, giving beta_direct and b;
  - mediator: e_B ~ e_A2 + controls, giving a;
  - indirect = a * b, which equals beta_total - beta_direct in OLS on the same rows and controls;
  - the proportion mediated, indirect / beta_total, is reported descriptively only.
- **Spillover sensitivity:** the same tests plus the outcome cell's own lineage scores for the other two chain types. Observed values and tumour-bootstrap CIs are reported; no permutations are run.

## Nulls (per section; 200 permutations each; add-one two-sided p = (k+1)/(n+1) on |statistic|)
- **Expression block null:** within each (200 um tile, cell type) group, permute program-score rows (all of a cell's programs move together). The tile grid gets a random offset, drawn uniformly in [0, 200) um for x and y, at each permutation. Labels and geometry are fixed. All exposures and controls are recomputed from the permuted scores.
- **Label block null:** within 200 um tiles with a random offset, permute the annotation labels of all cells. Each cell keeps its expression, since every program is scored on every cell. Chain membership, exposures and controls are recomputed.
- Both nulls extend the v3 `_within_group_perms` and `tile_codes` with random offsets. The same permutation index b is used for every chain and program in a section, which lets the tumour-level null pair across sections.

## Tumour level
- Tumour estimate: the mean of the section statistics over that tumour's passing sections.
- Tumour-level statistic: the mean of the tumour estimates over the tumours with >= 1 passing section. A chain needs >= 3 such tumours; otherwise it is "not testable".
- Tumour-level null: replicate b is the same mean computed from the b-th null replicate in each section. Add-one p-value.
- CI: tumour bootstrap, 2,000 resamples of tumours, 95 % percentile.
- **BH** is applied separately per statistic family (hop C<-B, hop B<-A, indirect) and per null, across all chains and program combinations. Duplicate tests are removed first: a hop C<-B test is identified by (B, C, B-program, C-program) and a hop B<-A test by (A, B, A-program, B-program). The number of tests in each family is reported.

## Negative controls (run for every selected chain)
- **Reversed chain:** C -> B -> A, with the same programs and roles swapped. Because spatial adjacency is undirected, a reversed chain passing as often as the forward one means that direction is not identifiable from these data. This will be reported, not hidden.
- **Random B:** B is replaced by one other cell type, drawn once per chain with seed 42 from the annotation classes that are not in the chain and have >= 100 cells in >= 3 tumours. The replacement cells carry the B-program score computed on them, since every program is scored on every cell.
- **Shuffled B scores:** B-program scores are permuted once (seed 42) among B cells within each section, and then the full test is run.

## Decision rule (tumour level), applied to each chain x program combination
**"Chain evidence"** requires all of the following:
1. Hop C<-B has BH q < 0.05 under the expression null and p < 0.05 under the label null.
2. Hop B<-A has BH q < 0.05 under the expression null and p < 0.05 under the label null.
3. Indirect has BH q < 0.05 under the expression null, the tumour-bootstrap 95 % CI of the indirect effect excludes 0, and |beta_direct| < |beta_total|, i.e. the A->C association shrinks when B is added.
4. Under the spillover sensitivity, the indirect effect keeps the same sign with a CI excluding 0.
5. The chain lies inside the Phase 4 calibrated range: effect size and cell counts at which the rule had power >= 0.8 and false-positive rates <= 0.10 under zonation and under direct A->C.

**Labels:**
- "supported (within calibrated range)" when 1-5 all hold.
- "exploratory" when 1-3 hold but 4 or 5 fails.
- "not supported" otherwise.
- "not testable" when fewer than 3 tumours pass the floor.

Negative-control outcomes are reported next to each label. If a negative control meets 1-3, the chain's label is annotated "not specific (negative control also passes)".

**Per-hop effect sizes:** beta_BA (upstream hop) and beta_CB (downstream hop) are reported side by side with |beta_CB| / |beta_BA|, so amplification (small upstream, larger downstream) is visible or absent.

## Outputs
results/chain/phase3/:
- per-section statistics with null p-values (sections/);
- tumour_level_tests.csv with q-values;
- negative_controls.csv;
- per_hop_effects.csv;
- chain_calls.csv;
- n_tests.json.

Figure figures/chain/phase3/: per-hop effects and the indirect effect with CIs for each chain.

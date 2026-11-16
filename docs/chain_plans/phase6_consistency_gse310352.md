# Chain Phase 6 plan: consistency in GSE310352

Committed before any GSE310352 chain computation. Patient independence from the Mendeley set is **UNVERIFIED**, so this is a **consistency check, not replication**. Licence: LICENSE_CONFIRMED_GSE310352.txt.

## Data and labels
- 8 slides, loaded with `pdac_loaders` (`geo_gse310352:<slide>`) and the phase3 `qc`. Labels are the v3 Phase 4 **marker** labels (results/pdac/v3_phase4/slides/<slide>_labels.csv.gz), aligned on cell_id. If an id is missing, the slide stops with an error; there is no silent realignment.
- **Units of inference:** slides 5, 5b2 and 5b3 are treated as one specimen unit ("slide5_group") because their names suggest serial sections. This grouping is an UNVERIFIED assumption. The other slides are one unit each, giving 6 units.
- Coordinates are in pixels and the um/px factor is UNVERIFIED.
  - **Tile scaling by nearest neighbour.** The Mendeley median nearest-neighbour distance across sections is 9.12 um (results/chain/phase1b/phase1_summary.json). So 200 um corresponds to 200 / 9.12 = 21.9 median-NN distances; the tile side in each slide is 21.9 x its own median positive NN distance in px.
  - The 30 um radius graph is likewise 30 / 9.12 = 3.29 x the median NN px.

## Programs
- The same gene lists as Mendeley Amendment 1 (results/chain/phase1b/programs_used.json), unchanged.
- Each gene's GSE310352 type-aware background status is computed with the same rule (own-type detection vs 2 x the maximum negative-probe detection in the same cells; flagged if near background in >= 50 % of slides with >= 50 cells of the type). It is reported as a diagnostic. A program whose usable genes fall below 5 in GSE310352 is reported as "measurement not consistent". The program is still scored with the Mendeley list, so the comparison stays like for like.
- Scoring is the same scanpy-style procedure (25 bins, 50 controls per gene, seed 42), with controls recorded per slide.

## Analyses
- Floors (>= 100 cells per type, >= 300 A-B and B-C edges) are computed per slide for the 6 Phase 3 chains. A chain is analysed if it passes in >= 3 units.
- **Phase 2 repeat:** ridge and gradient boosting, leave-one-unit-out cross-validation, the same features (graphs scaled as above), the 20,000-cells-per-slide cap, and a unit-bootstrap 95 % CI of delta R2.
- **Phase 3 repeat:** the same engine. 200 um-equivalent tiles with random offsets; expression block null and label block null with **100** permutations each (a reduction from 200 for compute, logged). Unit-level aggregation as in Phase 3 (mean over slides within a unit, then over units). The spillover sensitivity is included. Negative controls are not repeated (logged), because this is a consistency check of the forward chains.

## Reported
For each chain x program combination, and each of hop C<-B, hop B<-A and the indirect effect:
- the GSE310352 estimate, its unit-bootstrap CI and its add-one p values;
- **same direction as Mendeley: yes / no**, judged on the sign of the point estimates;
- whether the GSE310352 CI excludes 0.

The Phase 2 delta R2 is reported per receiver program, with "same direction" judged against the Mendeley delta R2.

Outputs: results/chain/phase6/. Figure: figures/chain/phase6/.

# Data card (model_dev v1)

## Sources and licences (as recorded in this repository)
- Development set: Mendeley PDAC-SMI CosMx, doi 10.17632/kx6b69n3cb.1, CC BY 4.0 (as recorded in `src/pdac/pdac_loaders.py`). 16 sections, 13 tumours; authors' major-type labels mapped to 13 classes (`src/pdac/phase3d_author_annotation.AUTHOR_MAP`).
- External set: GEO GSE310352 CosMx (Lyubetskaya et al., Cell Rep 2026, doi 10.1016/j.celrep.2025.116827, as recorded in the loader); terms confirmed in `LICENSE_CONFIRMED_GSE310352.txt`. 8 slides; marker labels from `results/pdac/v3_phase4`. Patient identity, independence from the development set and the um/px scale are UNVERIFIED; slides 5/5b2/5b3 are grouped as one unit by name (UNVERIFIED).

## Processing
- QC: min 10 counts per cell, min 3 cells per gene (`src/pdac/phase3_pdac_audit.qc`). Normalisation: CP10k + log1p. Gene set G: 539 of the 960 panel-intersection genes passing the type-aware background rule for >= 1 typed class (`results/chain/phase1b`). 10 GSE310352-only genes dropped: CEACAM5, FERMT1, LGALS4, LOX, NT5E, SLC5A3, TFF1, TFF3, TNS4, TSPAN8.
- Program scores: scanpy-style (25 bins, 50 controls/gene, fixed seeds) per section; identical to the chain Phase 1b cache (max abs diff 0.0).
- Graphs: symmetrised 6-NN; 30 um radius (external: 30 um-equivalent via median nearest-neighbour distance); 200 um tiles.

## Cells per group and type

| split       | group        |   B_plasma |   CAF_fibroblast |   NK |   T_cell |   acinar |   endocrine |   endothelial |   epithelial_ductal_tumor |   mast |   myeloid |   pericyte_smc |   schwann_nerve |   unassigned |
|:------------|:-------------|-----------:|-----------------:|-----:|---------:|---------:|------------:|--------------:|--------------------------:|-------:|----------:|---------------:|----------------:|-------------:|
| development | 004          |       1431 |            13799 |  320 |      879 |        0 |          53 |          1838 |                     41216 |    189 |      7470 |            953 |              62 |         6583 |
| development | 011          |       3021 |            30616 |  439 |     1623 |        0 |          79 |          3013 |                     35470 |   1582 |     11896 |           2645 |             370 |         3476 |
| development | 2462         |        342 |             7862 |   19 |       67 |        0 |         137 |           457 |                       222 |    146 |      1652 |            384 |             150 |         1035 |
| development | 2540         |       1312 |            14895 |  146 |      733 |        0 |           1 |           952 |                      5507 |    309 |      4321 |           1518 |              76 |         1314 |
| development | 2626         |       1895 |            33604 |  153 |      525 |        0 |         261 |          3034 |                     72197 |    634 |     22040 |           2552 |             348 |         3485 |
| development | 2675         |       1164 |            25996 |   68 |       94 |        0 |           1 |           548 |                      7133 |    664 |      5416 |            534 |              40 |         2000 |
| development | 2692         |        354 |            14855 |  324 |      822 |        0 |         597 |          1105 |                     16099 |    540 |      5847 |           1014 |             147 |         2090 |
| development | 2762         |       4250 |            46799 |  988 |     4548 |        0 |         136 |          2449 |                     11218 |   1098 |     18742 |           1723 |             290 |         5560 |
| development | 2924         |       6862 |            22865 |  357 |     1045 |        0 |         358 |          3609 |                     22530 |    789 |      5775 |           2318 |             238 |         4139 |
| development | 2997         |        230 |            18635 |   96 |      230 |        0 |         274 |          1076 |                      5900 |     75 |      3153 |           1372 |            1498 |         1261 |
| development | B10          |        165 |             5696 |   60 |      243 |        0 |           0 |           679 |                      3646 |    115 |      3144 |            441 |              27 |          803 |
| development | C4           |        951 |             8661 |  180 |      628 |        0 |           0 |           981 |                      4920 |    244 |      3729 |           2357 |              86 |         2157 |
| development | D10          |        221 |             8734 |   61 |      120 |        0 |           1 |           658 |                     20191 |    274 |      2756 |            326 |              22 |          970 |
| external    | slide1       |       6158 |            29878 |    0 |     2184 |     3140 |        1816 |          2636 |                     48431 |   3392 |      7863 |           3632 |               0 |            0 |
| external    | slide2       |       6586 |            23919 |    0 |     3041 |    12165 |        2020 |          3073 |                     40769 |   3454 |     15902 |           3154 |               0 |            0 |
| external    | slide3       |       2023 |            14692 |    0 |     2078 |     2583 |         302 |          2612 |                     30573 |   1594 |      8978 |           1164 |               0 |            0 |
| external    | slide4       |      14451 |            19380 |    0 |     6528 |     3514 |        4255 |          2835 |                     18815 |   2228 |     12652 |           1801 |               0 |            0 |
| external    | slide5_group |       3596 |            32948 |    0 |     2966 |    22591 |        1889 |          6536 |                     70069 |   3827 |     12261 |          10763 |               0 |            0 |
| external    | slide6       |       1348 |            35088 |    0 |     1479 |     2671 |         655 |          2443 |                     28615 |   1472 |     10438 |           4334 |               0 |            0 |

## Exclusions
- Cells failing QC (see `results/phase0/section_info.csv` for post-QC counts).
- T cells are not a modelling target (declared not testable: too few cells per fold).
- No cells were excluded on the basis of any model output.

# Which setting produced which number (breast Xenium run)

Every number below comes from the **same pipeline** (same processed AnnData `data/processed/xenium_breast_baseline.h5ad`, same per-cell marker labels, same label-permutation null idea). The **graph construction and permutation settings differ**, so the numbers are **not interchangeable**. Do not quote a radius-sweep peak as if it were the contact-restricted result, and do not mix them in one table without the setting column.

## (a) Contact-restricted, 15 µm radius graph — the setting used in the paper

| Axis (sender → receiver) | z | p (legacy, `mean(null >= obs)`) | Source file |
|---|---|---|---|
| CD86 → CTLA4 (immune → immune) | +3.56 | 0.035 | `results/multi_axis_expansion/multi_axis_summary.csv` (row radius_um = 15) |
| CXCL12 → CXCR4 (fibroblast → immune) | +36.20 | 0.000 | same file |
| PTN → SDC4 (fibroblast → epithelial) | −1.36 | 0.920 | same file |
| CD274 → PDCD1 (epithelial → immune) | 0.00 | 1.000 | same file |

Settings: radius graph r = 15 µm (`sklearn` `radius_neighbors`, self-loops removed, 40,618 directed edges); 200 label permutations from `numpy.random.RandomState(42)`; edge score = mean of L(i)·R(j) over qualifying edges.
Reproduced exactly on 2026-10-24 by `src/pdac/phase0_repro_check.py` → `results/pdac/phase0/multi_axis_15um_recheck.csv`. With the add-one p-value convention, (k+1)/(n+1), CD86 → CTLA4 gives p = 0.0398 (same file). Both conventions are below 0.05, but 200 permutations give coarse p resolution.

## (b) Radius-sweep peaks — sensitivity analysis only, NOT the paper result

| Axis | Peak z | Radius at peak | Source file |
|---|---|---|---|
| CXCL12 → CXCR4 (fibroblast → immune) | +79.42 | 200 µm | `results/multi_axis_expansion/multi_axis_summary.csv` (200 permutations) |
| CXCL12 → CXCR4 (fibroblast → immune) | +79.79 | 200 µm | `results/distance_decay_horizon/distance_decay_summary.csv` |
| CXCL12 → CXCR4 (fibroblast → immune) | +83.70 | 300 µm | `results/distance_decay_horizon/distance_decay_summary.csv` (largest radius tested; z still rising) |
| CD86 → CTLA4 (immune → immune) | +19.96 | 200 µm | `results/multi_axis_expansion/multi_axis_summary.csv` (largest radius tested in that sweep) |

The sweep z-scores grow with radius mainly because the number of qualifying edges grows by orders of magnitude (560 → 83,580 edges for CXCL12 fibroblast → immune between 15 and 200 µm). That shrinks the null standard deviation. A larger z at a larger radius therefore does **not** mean stronger or more local communication. The z peaks at the largest radius tested, so these are not interior optima.

## (c) 6-nearest-neighbour graph (k = 6) — the main null-model comparison

`results/spatial_null_corrected.csv`, `results/spatial_edge_null.csv` and `results/null_model_comparison.csv` use the Squidpy 6-NN graph (42,978 directed edges), 500 permutations, `numpy.random.default_rng(42)`. Example: CXCL12 → CXCR4 fibroblast → immune, edge z = +29.51, compositional z = +31.03.

## Terminology

The margin-vs-core split (`src/run_cross_fov_replication.py`, `results/cross_fov_replication/`) cuts **one tissue section** in half at the median x coordinate. It is an **internal regional consistency check**, not a replication. File and directory names still say `replication` only so that existing links keep working.

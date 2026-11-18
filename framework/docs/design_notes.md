# cellalibi Design Notes: Assumptions and Defaults

This document records the empirical rationale for all default parameter values, spatial thresholds, and statistical assumptions used throughout the `cellalibi` framework.

---

## 1. Spatial Contact Graph: $k = 6$ Nearest Neighbours
- **Default:** $k = 6$ nearest neighbours on physical coordinates $(x, y)$ in micrometres.
- **Why:**
  - In close-packed 2D cellular monolayers and tissue sections, a central cell is adjacent to approximately 6 immediate topological neighbours (Euler's formula for planar hexagonal Voronoi packings).
  - $k=6$ directly matches standard literature baselines (Squidpy default, CellPhoneDB spatial extensions).
  - Ensures graph connectivity across sparse stromal regions without creating spurious long-range connections across empty tissue voids.
- **Alternative:** Radius graphs with threshold $r = 15-20\,\mu\text{m}$ (approximately 1.5 cell diameters) are supported when cell density varies dramatically.

---

## 2. Permutation Count: $n = 500$
- **Default:** 500 label permutations per section / run.
- **Why:**
  - With 500 permutations, the add-one empirical p-value resolution is $1 / (500 + 1) \approx 0.00199$.
  - This resolution allows standard significance testing at $\alpha = 0.05$ and multiple-testing corrections across typical targeted panel sizes ($10-30$ ligand-receptor pairs) without Monte Carlo noise destabilising rankings.
  - Across large tissue sections ($50,000-80,000$ cells), 500 vectorised permutations execute in $< 5$ seconds on 4 CPU threads while keeping RAM usage $< 250$ MB.
  - For quick interactive screening, $n=100$ is supported (resolution $0.01$).

---

## 3. Spatial Block Tile Sizes: $50\,\mu\text{m}$, $100\,\mu\text{m}$, and $200\,\mu\text{m}$
- **Default:** $100\,\mu\text{m}$ square tiles with random coordinate offset $(u_x, u_y) \sim \text{Uniform}(0, \text{tile\_size})$.
- **Why:**
  - **$200\,\mu\text{m}$:** Matches macroscopic tumour nest compartmentalisation and hypoxia diffusion limits (Krogh cylinder oxygen diffusion radius in solid tumours is approximately $100-150\,\mu\text{m}$).
  - **$100\,\mu\text{m}$:** Intermediate microenvironment scale (contains approximately $30-80$ cells). Controls for regional composition while providing enough cells per tile to permute without excessive singleton tiles.
  - **$50\,\mu\text{m}$:** Tight juxtacrine niche scale (contains approximately $10-25$ cells). Permuting within $50\,\mu\text{m}$ tiles isolates whether cells that touch have higher signal than cells separated by just 2–3 cell bodies in the same niche.
  - **Random coordinate offset:** Anchoring grid tiles to arbitrary fixed coordinates causes cells near tile boundaries to be artificially separated. Adding a random offset ensures that boundary effects average out across replicates.

---

## 4. Spillover Distance Threshold: 1 Median Cell Diameter ($\approx 12\,\mu\text{m}$)
- **Default:**
  - `touching`: Distance to nearest sender $\le d$ (1 median diameter, typically $10-14\,\mu\text{m}$).
  - `near`: $d < \text{distance} \le 2d$ ($12-24\,\mu\text{m}$).
  - `separated`: $\text{distance} > 2d$ ($> 24\,\mu\text{m}$).
- **Why:**
  - Optical point spread functions and automated nuclear expansion algorithms (e.g. 2–3 $\mu\text{m}$ expansion around DAPI) cause cell boundaries to overlap.
  - Empirical measurements on CosMx PDAC show that cells within 1 diameter of a CAF carry 23% excess collagen/matrix RNA.
  - At $> 2$ cell diameters, physical membrane overlap drops to near-zero, providing a pure, uncorrupted reference population of isolated cells.

---

## 5. Statistical Power Detection Floor: 25 Receiver Cells
- **Default:** Interactions with $< 25$ receiver cells or $< 10$ qualifying contact edges are marked `NOT_TESTABLE`.
- **Why:**
  - Extensive ground-truth simulation calibration across 3 real PDAC tissue sections (T6-b with 26 T cells, T2-b with 23 T cells, and U5-a with 261 T cells) showed that statistical power is precisely **0.00** across all nulls for sections with $\approx 25$ T cells, even when a 2-fold ($200\%$) contact effect is planted.
  - Calling such interactions "non-significant" or "negative" is a false negative caused by statistical blindness. They are unpowered and untestable.

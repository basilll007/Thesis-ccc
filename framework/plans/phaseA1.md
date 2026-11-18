# Phase A1 Plan: Package Design (`ccc_audit`)

## 1. Objective and Architecture
Turn the thesis spatial analysis, null models, calibration, spillover audits, and neighbourhood models into an installable, tested, reusable Python package under `framework/` with the working name `ccc_audit`.

### Package Layout
```
framework/
├── pyproject.toml
├── requirements-dev.txt
├── plans/
│   └── phaseA1.md
├── src/
│   └── ccc_audit/
│       ├── __init__.py
│       ├── io.py
│       ├── graphs.py
│       ├── scores.py
│       ├── nulls.py
│       ├── calibration.py
│       ├── spillover.py
│       ├── neighbour_models.py
│       ├── verdicts.py
│       └── cli.py
```

## 2. Module Specifications & Contracts

### 2.1 `io.py`
- **Input Contract:** `AnnData` object with required `.obs` columns:
  - `cell_type` (str / categorical): cell type label.
  - `patient_id` (str / categorical): patient (or tumour) identifier, unit of independent grouping.
  - `section_id` (str / categorical): tissue section / slide identifier.
  - `x_um` (float): physical x coordinate in micrometres.
  - `y_um` (float): physical y coordinate in micrometres.
  - Optional: negative control probe channels (`neg_probes` or matching prefix `NegPrb`), protein channels, raw count layer `layers['counts']`.
- **Validation:** Clear error messages if columns are missing, coordinates non-numeric, or cell counts insufficient.

### 2.2 `graphs.py`
- Construct spatial graphs:
  - kNN graph (default k=6, Euclidean on physical `(x_um, y_um)`).
  - Radius graph (distance threshold in um).
- Tile assignment with random coordinate offset $(u_x, u_y) \sim \text{Uniform}(0, \text{tile\_size})$ to prevent artificial boundary edge artefacts.

### 2.3 `scores.py`
- Vectorised audit scores (citing `src/pdac/ccc_audit.py` and `src/pdac/phase1_matched_baseline.py`):
  - Compositional score $C_{sr} = \frac{\bar{L}_s + \bar{R}_r}{2}$ (Squidpy / CellPhoneDB style).
  - Edge score $E_{sr} = \frac{1}{|E_{sr}|} \sum_{(i \to j) \in G, t(i)=s, t(j)=r} L_i R_j$.
  - Matched dyadic baseline $M_{sr} = \bar{L}_s \cdot \bar{R}_r$ (product form matching $E$, but without spatial restriction).
  - Proximity ratio $Q_{sr} = \frac{E_{sr}}{M_{sr}}$ (edge enrichment over matched non-spatial expectation).
- Add-one empirical p-values: $p = \frac{k + 1}{n + 1}$ where $k$ is the number of null draws $\ge$ observed.
- Standardised z-scores: $z = \frac{\text{obs} - \mu_{\text{null}}}{\sigma_{\text{null}}}$.

### 2.4 `nulls.py`
- Permutation generators preserving global or local cell-type counts:
  - Global label permutation: random permutation of cell-type labels across all cells.
  - Block label permutation: random permutation of labels within square tiles (e.g. 50, 100, 200 $\mu$m) with random anchor offset.
  - Within-tile expression null: permute expression values within spatial tiles keeping labels fixed.

### 2.5 `calibration.py`
- Ground-truth simulation calibrations on real spatial geometries (citing `src/pdac/v3_phase2_calibration.py`):
  - Planted proximity: scale receptor counts by factor $f \in [1.0, 2.0]$ in receiver cells adjacent to ligand-positive senders.
  - Planted zonation-only: scale ligand and receptor counts by factor $f$ for all cells in random $200\,\mu$m regional tiles, without spatial contact dependence.
  - Outputs: empirical power and false-positive rate (FPR) tables.

### 2.6 `spillover.py`
- Segmentation spillover diagnostics (citing `model_dev/src/v2_phase1.py` and `model_dev/reports/MODEL_REPORT_v2.md`):
  - Automatic detection of sender contamination markers (e.g., top CAF collagen / matrix genes vs isolated receivers).
  - Physical distance binning based on median cell diameter $d$ (touching $\le d$, near $d < \text{dist} \le 2d$, separated $> 2d$).
  - Per-cell contamination covariate score $s_c$.
  - Target adjustment: removal of shared sender/receiver genes.
  - Negative-control spillover-free genes $S$ (high expression in receiver, zero/negligible in sender).
  - Mixing expectation model comparing observed vs expected transcript sharing.

### 2.7 `neighbour_models.py`
- Machine learning neighbourhood prediction models (citing `model_dev/src/v2engine.py` and `baselines.py`):
  - Baseline own covariates vs neighbourhood features (same-type vs cross-type neighbours).
  - Grouped cross-validation grouped strictly by `patient_id`.
  - Paired patient/tumour-level bootstrap confidence intervals (2,000 resamples).
  - Ridge regression and XGBoost / tree-based baselines.

### 2.8 `verdicts.py`
- Unified decision logic mapping an interaction or receiver program into an evidence classification:
  - `NOT_TESTABLE`: receiver count < power floor (default 25 cells), qualifying contact edges < 10, or sparse gene detection < 20 cells.
  - `EXPLAINED_BY_ABUNDANCE`: compositional score $C$ or baseline $M$ is significant, but $Q$ fails significance ($p \ge 0.05$ or $z_Q \approx 0$).
  - `EXPLAINED_BY_ZONATION`: $Q$ is significant under global null but drops to non-significance under block null, or simulation indicates false-positive rate $> 10\%$ under regional zonation.
  - `EXPLAINED_BY_SPILLOVER`: predictive gain or enrichment is confined to touching cells ($\le 12\,\mu\text{m}$) and vanishes in separated cells; or fails when adjusting for contamination covariate $s_c$; or is matched on negative-control genes.
  - `SUPPORTED`: passed power floor, $Q > 1$ with $p < 0.05$ under spatially constrained nulls, persists in separated cells, and survives spillover controls.
- Output: structured claims-register table.

### 2.9 `cli.py`
- Single command-line interface: `ccc-audit run config.yaml` producing summary tables, diagnostic figures, and markdown reports.

## 3. Resource & Safety Constraints
- Execution limited to $\le 4$ CPU worker threads.
- Free RAM logged before runs, maintaining $\ge 8$ GB.
- Code reuse via explicit import or clean copying with comments naming original source files.

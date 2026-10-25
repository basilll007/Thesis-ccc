# PDAC v2 — protein-guided annotation rule (pre-registered)

Written on 2026-10-25, **before** protein labels or any CCC result under them were computed. It was committed on its own, ahead of the code that applies it (see `git log -- docs/PDAC_V2_PROTEIN_RULE.md`). The rule is audit tooling, not a biological claim.

## Inputs

- Dataset: Mendeley PDAC-SMI (CosMx), file `data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad`.
- Per-cell protein channels in `obs`: `Mean.CD3`, `Mean.CD45`, `Mean.PanCK`, `Mean.MembraneStain`, `Mean.DAPI`, plus the matching `Max.*` columns.
- No other protein channels are present (checked 2026-10-25).
- **Used for gating**: `Mean.CD3`, `Mean.CD45`, `Mean.PanCK`.
- **Not used for gating**: `Mean.MembraneStain` and `Mean.DAPI`. They are segmentation stains, not lineage markers. They are reported as diagnostics only.
- RNA: the same log-normalised expression and the same `PDAC_MARKERS` sets as the v1 marker scoring (`src/pdac/phase3_pdac_audit.py`), used only *within* a protein-defined compartment (see the labelling table).
- The author labels (`annotation_majortypes`) are **not** used to build the rule or to choose thresholds. They are used only for validation afterwards.

## Thresholds (per section, per channel)

Thresholds are set per section because staining intensity differs between slides.

For each section and each channel c ∈ {CD3, CD45, PanCK}, let x = log1p(`Mean.c`).

1. **Primary: two-component Gaussian mixture.**
   - Fit `sklearn.mixture.GaussianMixture(n_components=2, covariance_type="full", n_init=5, random_state=42)` on x.
   - If the section has more than 50,000 cells, fit on 50,000 cells sampled with `numpy.random.default_rng(42)`; otherwise fit on all cells.
   - The positive component is the one with the higher mean.
   - Threshold t_c = the smallest x on a 2,000-point grid between the two component means where the posterior of the positive component is ≥ 0.5.
   - **Acceptance check**: Ashman's D = √2·|μ₊ − μ₋| / √(σ₊² + σ₋²) ≥ 2, and the weight of the smaller component ≥ 0.005.
2. **Fallback, used only when the acceptance check fails**: t_c = median(x) + 3 × 1.4826 × MAD(x), a robust outlier rule that does not assume a prevalence. The diagnostics table records which method set each threshold.
3. **Positive**: a cell is c-positive when x > t_c.

**Pre-declared sensitivity variant `protein_mad`**: the same gates, with the robust rule (2) used for every channel in every section.

## Labelling (applied in this order)

Notation: P = PanCK+, H = CD45+, T = CD3+.

| Protein state | Label |
|---|---|
| P and (H or T) | `protein_ambiguous` (kept in the graph as its own type) |
| T and H and not P | **`T_cell`** |
| T and not H and not P | `T_cell_CD45neg` (separate type; not the headline receiver) |
| H and not T and not P | RNA argmax among {`myeloid`, `B_plasma`, `mast`} |
| P and not H and not T | RNA argmax among {`epithelial_ductal_tumor`, `acinar`, `endocrine`} |
| not P, not H, not T | RNA argmax among {`epithelial_ductal_tumor`, `acinar`, `endocrine`, `CAF_fibroblast`, `endothelial`, `pericyte_smc`} |

The RNA argmax uses the same per-cell mean-marker scores as v1. A type needs at least 2 of its markers on the panel. NK cells have no protein marker here, so they fall into the CD45+CD3− immune group. That is a known limitation.

Validation also reports a **loose T** definition, T and not P (CD45 not required), which merges `T_cell` and `T_cell_CD45neg`. It is not used for the CCC run.

## Validation (after labelling, before CCC interpretation)

- Author labels are mapped to coarse types with the v1 `AUTHOR_MAP` (`src/pdac/phase3d_author_annotation.py`). Author T = {CD8+ T, CD4+ T, Treg}; author T/NK adds Natural killer.
- Report a confusion matrix plus per-class precision and recall, for U5-a and for all 16 sections, for the protein labels and for the v1 marker labels.

## CCC re-run

- Same audit as v1: C, E, M, Q; 500 label permutations (`spatial_null.generate_label_permutations`, seed 42); 6-NN graph; within-section permutations; the full 31-pair testable catalogue; all sender × receiver combinations.
- Annotations: `protein_gmm` (primary) and `protein_mad` (sensitivity). Compare with v1 marker and author labels.
- Significance: add-one p < 0.05; BH within section (and within annotation).
- **Headline**: CXCL12→CXCR4, `CAF_fibroblast` → `T_cell`. The question is whether Q (proximity-specific enrichment over the matched product baseline) is significant, uncorrected and after BH within section. The answer is reported whichever way it goes.

---

## Addendum (2026-10-25 18:35): POST-HOC consensus annotation, added after the protein rule failed validation

Written **after** the protein labels were validated against the author labels and **before** any CCC was run on any v2 annotation. It is post-hoc and is labelled that way in every output.

**Validation result for the pre-registered rule** (`results/pdac/phaseB/protein_labels_summary.json`, `validation_precision_recall.csv`):
- `T_cell` precision and recall against author T: **U5-a 0.022 / 0.036**; pooled 16 sections 0.077 / 0.049.
- The v1 marker labels score 0.276 / 0.297 (U5-a) and 0.353 / 0.318 (pooled).
- `protein_mad` gives 0 T cells in U5-a.
- Cause: in this dataset, `Mean.CD3` and `Mean.CD45` separate the author T cells only weakly.
  - Per-section AUC for author T: CD3 0.45–0.92; CD45 0.50–0.89 (U5-a: 0.70 and 0.55).
  - GMM acceptance: CD3 3/16 sections, CD45 2/16, PanCK 15/16.
  - PanCK separates author malignant cells well (AUC 0.74–0.99).
- **The pre-registered protein rule does not yield validated T-cell labels.** Its CCC re-run is still carried out as pre-registered, but it is reported as "labels failed validation".

**Post-hoc annotation `consensus_T`**:
- Start from the author labels mapped with the v1 `AUTHOR_MAP`.
- A cell is **`T_cell`** only if the author label is T (CD8+ T, CD4+ T or Treg) **and** the v1 marker scoring also calls it `T_cell`. These are two independently derived RNA-based annotations that agree.
- An author T cell not confirmed by marker scoring becomes `T_cell_unconfirmed`, a separate type that is not the headline receiver.
- All other cells keep their author-mapped label.
- Precision against the author labels is 1 by construction, so it is not used as evidence. Instead, a **third, independent source** is checked: CD3 and CD45 protein intensity in consensus-T cells vs author-only T cells vs marker-only T cells vs all other cells.
- CCC: same audit as above (C/E/M/Q, 500 perms, seed 42, 6-NN, full catalogue, 16 sections).

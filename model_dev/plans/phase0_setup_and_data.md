# Phase 0 plan: setup and data build

Committed before the data build is run. Branch `model-v1` (from origin/main at cfb111c). All code and small outputs live in `model_dev/`; large artefacts go to `F:\Thesis\model_cache\` (git-ignored) with paths and SHA-256 checksums recorded in `model_dev/results/phase0/artefacts.csv`.

## Hardware and environment (recorded at setup)
- GPU: NVIDIA GeForce RTX 5060 Laptop GPU, 8151 MiB VRAM (driver 592.82, CUDA toolkit 12.8 on the machine). RAM: 31.1 GB. 32 logical CPUs.
- The existing env `thesis-spatial-ccc` has no torch. A new conda env `thesis-model` (Python 3.11) is created with a CUDA 12.8 torch build (the RTX 5060 is a Blackwell part and needs it), PyTorch Geometric, scanpy, scikit-learn, pytest. Pinned versions are exported to `model_dev/environment.yml` after installation.
- The system drive C: has about 0.1 GB free. Temp and pip caches are redirected to `F:\Thesis\model_cache\`. Nothing is deleted from C:.

## Data
- **Development set:** Mendeley PDAC-SMI CosMx, 16 sections, 13 tumours (donor with `_2nd`/`_B` stripped, as in `src/pdac/v3_nulls.section_to_tumour`). Labels: `author_mapped` (authors' major types mapped to 10 classes plus NK, schwann_nerve, unassigned; from `results/pdac/phaseB/sections`). QC and normalisation reuse `v3_nulls.load_section` (min_counts 10, min_cells 3; normalize_total 1e4; log1p).
- **External test set:** GSE310352 CosMx, 8 slides, marker labels from `results/pdac/v3_phase4/slides/<slide>_labels.csv.gz`. Loaded with `pdac_loaders.load("geo_gse310352:<slide>")` and the same QC. It is built in Phase 0 only to record its data card (cells per type per slide and the gene intersection). **No model, scaler, gene filter or threshold touches it until Phase 4.** The commit hash at first load is written to `model_dev/results/phase0/external_first_load.json`.

## Gene set (declared here)
- Panel intersection: Mendeley Base genes (960) and GSE310352 genes (970). Checked before this plan was written: the intersection is 960, so no Mendeley gene is lost; the 10 GSE-only genes (CEACAM5, FERMT1, LGALS4, LOX, NT5E, SLC5A3, TFF1, TFF3, TNS4, TSPAN8) are dropped. Negative probes are never model inputs.
- **Model gene set G:** genes that pass the type-aware background rule for at least one of the four typed classes (tumour, CAF, T, myeloid) in the development set (`results/chain/phase1b/gene_background_typeaware_flags.csv`, `background_level_for_type == False` for >= 1 type). Written to `model_dev/results/phase0/gene_set.csv` with the per-type flags. This rule was fitted on the development set only (all 16 sections). It is a gene *filter* fixed before any split; the inner-fold feature selection in later phases happens inside training folds.
- Program definitions are reused unchanged from `results/chain/phase1b/programs_used.json` (tumour classical and basal_like, CAF myCAF and CXCL12, T cytotoxic, myeloid C1Q_TAM and SPP1_TAM) with their lineage sets.

## Per-section graph dataset (one `.npz` per section/slide in `model_cache/datasets/<dataset>/<section>.npz`)
- `cell_id`, `xy` (um for Mendeley; px for GSE310352 with um/px UNVERIFIED), `label` (string class), `tumour` (Mendeley) or `unit` (GSE: slides 5, 5b2 and 5b3 grouped as one UNVERIFIED unit, as in chain Phase 6), `total_counts`, `n_genes`.
- `X_log` : log1p CP10k expression over G as float16 (cells x |G|), and `X_counts` as uint16 for the same genes.
- Program scores: scanpy-style score (25 expression bins, 50 controls per gene, seed by section index) via `src/chain/chain_phase1_cell_states.score_sets`, with controls recorded. For Mendeley the scores are asserted equal to `data/chain_cache/mendeley_v1b/<section>.npz` (a build-time reproduction check).
- Graphs: symmetrised 6-NN (`chain_phase1_cell_states.sym_knn`) and 30 um radius (GSE: 30 um-equivalent = 30 / 9.12 x median NN distance in px, as in chain Phase 6). Stored as CSR arrays plus edge lengths.
- 200 um tile codes (GSE: 200 um-equivalent) for the zonation covariate and the shuffled-neighbour ablation.

## Reported in Phase 0
- Cells per type per section and per tumour (both datasets), `model_dev/results/phase0/cell_counts.csv`.
- Genes in G, and per-type usable counts.
- Artefact registry with checksums.

## Success criterion
All 16 + 8 datasets build; Mendeley program scores reproduce the chain cache exactly (max abs diff < 1e-4); the gene intersection is as stated.

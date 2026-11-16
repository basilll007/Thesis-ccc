# Thesis Prototype v0 — Spatial Cell–Cell Communication

> **Key findings across all workstreams: [`FINDINGS.md`](FINDINGS.md).** Dashboard: `index.html` (incl. `neighbourhood_models.html`).

This repository provides a reproducible baseline for cell–cell communication analysis in publicly available, single-cell-resolution spatial transcriptomics data from breast cancer. It performs basic quality control and normalization, Leiden clustering, coarse marker-based cell-type annotation, Squidpy spatial-neighbor graph construction, and a baseline Squidpy ligand–receptor analysis. Functional-response evidence, reproducibility calibration, and Spatial Prize-Collecting Steiner Tree optimization are future work and are not part of v0.

## Public dataset

The data are intentionally **not stored in this repository**. Download the official 10x Genomics Xenium V1 Human Breast 2-FOV public output bundle:

- Dataset page and examples: [10x Genomics Xenium example data](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/2.0/resources/xenium-example-data)
- Direct public bundle: [Xenium V1 Human Breast 2-FOV outputs](https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip)
- License: Creative Commons Attribution 4.0 International (CC BY 4.0), as stated by 10x Genomics.

The archive is approximately 379 MB. Its URL, downloaded size, checksum, date, and access terms are recorded in `docs/PROJECT_CONTEXT.md`.

## Environment

The project targets Python 3.11. Create the environment with conda/mamba from `environment.yml`:

```powershell
conda env create -f environment.yml
conda activate thesis-spatial-ccc
```

Alternatively, use a Python 3.11 virtual environment and uv:

```powershell
uv venv .venv --python 3.11
.\.venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
```

## Download the data

From the repository root, run:

```powershell
python src\download_data.py
```

The downloader retrieves the official public bundle, unpacks its nested cell-feature matrix, and records provenance. Raw files are written under `data\raw\`; the complete `data\` directory is excluded from Git because it contains downloaded or derived data.

## Run the prototype

```powershell
python src\run_prototype.py
```

The fixed random seed is 42. The pipeline performs:

1. Xenium loading into AnnData and basic QC.
2. Normalization, Leiden clustering, and coarse marker-based labels.
3. Squidpy spatial-neighbor graph construction using six nearest neighbors.
4. Baseline `squidpy.gr.ligrec` analysis over a curated ligand–receptor panel.
5. Export of the ranked table, processed AnnData, run log, and required PNG/PDF figures.

Generated paths are:

- `data\processed\xenium_breast_baseline.h5ad` — excluded from Git.
- `results\baseline_lr_ranked.csv` — included as a lightweight baseline result.
- `figures\` — included PNG/PDF baseline figures.
- `logs\prototype.log` — included execution log.

The v0 ligand–receptor run uses one deterministic permutation as an end-to-end smoke validation. Its p-values are not calibrated significance estimates.

## Verification

```powershell
python -c "import scanpy, squidpy, anndata, liana, numpy, scipy, pandas, networkx, matplotlib, seaborn, sklearn, leidenalg, igraph; print('imports OK')"
uv pip check
```

See `docs/PROJECT_CONTEXT.md` for scientific context and provenance and `docs/PROGRESS_LOG.md` for the append-only implementation record.

## PDAC extension

This extends the spatial CCC audit to pancreatic ductal adenocarcinoma, using the public Mendeley PDAC-SMI CosMx dataset (Shiau, Cao et al.; CC BY 4.0; doi 10.17632/kx6b69n3cb.1). It also adds a matched dyadic baseline (L̄s·R̄r and E/(L̄s·R̄r)) to the breast analysis. Results are audit and diagnostic only.

- Summary of results, failures and limitations: `docs/PDAC_SUMMARY.md`
- **v2 (branch `pdac-v2`)**: [`docs/PDAC_SUMMARY_v2.md`](docs/PDAC_SUMMARY_v2.md). It covers:
  - spot-check of the v1 headline numbers
  - p-value conventions side by side (see also [`docs/RESULTS_SETTINGS_NOTE.md`](docs/RESULTS_SETTINGS_NOTE.md))
  - a pre-registered protein-guided T-cell rule ([`docs/PDAC_V2_PROTEIN_RULE.md`](docs/PDAC_V2_PROTEIN_RULE.md)), which failed validation, and a post-hoc consensus T class
  - the CCC re-run on 16 sections
  - the GSE310352 replication (BLOCKED)

  Lead question: *does a proximity-specific signal exist at tumour level, or is the edge score mainly mirroring cell-type products?*

  Answer: mainly mirroring. The CXCL12→CXCR4 CAF→T ratio Q is significant after BH in 0/16 sections and 0/13 tumours under every annotation. Across tumours it shows only a weak, exploratory directional trend (consensus T labels: z_Q > 0 in 10/13 tumours). Tumour-level tables: [`results/pdac/phaseD/tumour_level_headline_summary.csv`](results/pdac/phaseD/tumour_level_headline_summary.csv). Key numbers with source files: [`results/pdac/phaseD/v2_key_numbers.csv`](results/pdac/phaseD/v2_key_numbers.csv).
- **v3 (branch `thesis-v3`)**: [`docs/PDAC_SUMMARY_v3.md`](docs/PDAC_SUMMARY_v3.md). Lead question: does the proximity-specific signal survive spatially constrained nulls, and could this audit detect a real effect of the size we would care about? It covers:
  - Phase 1: block-permutation nulls (200 / 100 / 50 um tiles) and a distance-stratified null on the CXCL12→CXCR4 CAF→T headline and the 7 non-background axes, with tumour-level Fisher + BH. The Q signal does not shrink (median tumour z_Q 0.43 → 1.30 from global to block50, consensus_T); plan `docs/v3_plans/phase1_contact_vs_zonation.md`.
  - Phase 2: calibration by planting contact effects (1.1x–2x) and zonation-only effects on three real sections. Q has zero power at ~25 T cells; in the largest section it detects 1.5x but every null except block200 is fooled by 200 um zonation. Table: [`results/pdac/v3_phase2/calibration_deliverable.csv`](results/pdac/v3_phase2/calibration_deliverable.csv).
  - Phase 3: block-bootstrap CIs for Q, an exploratory mixed model (geometric-mean Q 1.24 [1.01, 1.52], consensus_T), and z vs edge-count inflation.
  - Phase 4 (GSE310352): DONE after the user confirmed the data terms, with 8 slides and 632,260 post-QC cells. Marker-label headline Q>1 in all 8 (median 1.408), with 0 global-null BH-significant slides and 0–1 under each constrained null. Patient independence and physical scale remain unverified; 3 slides have headline background-gene flags and Leiden recovers T cells on only 2 slides. Tables: [`results/pdac/v3_phase4/headline_summary_by_null.csv`](results/pdac/v3_phase4/headline_summary_by_null.csv); script: `src/pdac/v3_phase4_consistency.py`.
  - Phase 5 (GSE277783 programs): DONE after downloading GEO SOFT metadata. Fifty-four included human slides yield 13 paired patients; the hypoxia-associated program is higher in metastases (exact paired p=0.00244, BH q=0.01221). The other four programs do not pass BH. These are mixed-spot expression associations. Tables: [`results/pdac/v3_phase5/program_comparisons.csv`](results/pdac/v3_phase5/program_comparisons.csv); script: `src/pdac/v3_phase5_programs.py`.
  - Claims register and all quoted numbers with sources: [`results/pdac/v3_phase6/claims_register.csv`](results/pdac/v3_phase6/claims_register.csv), [`results/pdac/v3_phase6/v3_key_numbers.csv`](results/pdac/v3_phase6/v3_key_numbers.csv). Plans written before computation: `docs/v3_plans/`.
- Step-by-step log with commands, outputs and decisions: `docs/PDAC_PROGRESS_LOG.md`
- Manifest of every produced file (path, phase, description, script): `results/pdac/MANIFEST.csv`
- Which breast number comes from which setting (15 µm vs radius sweep vs 6-NN): `docs/RESULTS_SETTINGS_NOTE.md`
- Code: `src/pdac/` (run logs in `logs/<script>_<timestamp>.log`; figures in `figures/pdac/<phase>/`)

Data are not committed. Download `raw_meta_data_final.h5ad` from the Mendeley dataset page into `data/pdac/mendeley_pdac_smi/`. COMMOT and SpatialDM were run in a separate conda env, `pdac-bench`: a clone of `thesis-spatial-ccc` plus `commot` 0.0.3 and `SpatialDM` 0.3.1, with numpy pinned below 2.

## Chain extension (branch `chain-v1`)
Summary: [`docs/CHAIN_SUMMARY_v1.md`](docs/CHAIN_SUMMARY_v1.md). Plans (committed before computing): `docs/chain_plans/`. Log: `docs/CHAIN_PROGRESS_LOG.md`. Code: `src/chain/`. Results: `results/chain/` (manifest `results/chain/MANIFEST.csv`). Figures: `figures/chain/`.

Question: do neighbours predict a cell's expression-program state beyond its own type, and is there support for A -> B -> C chains (e.g. tumour -> CAF -> T cell) within a simulation-calibrated range? Answer from the generated files:
- Neighbourhood features improve held-out prediction for 6/7 programs (e.g. tumour classical ridge delta R2 0.141 [0.105, 0.164]). The gain comes from same-type neighbours and composition; other cell types' programs add <= 0.0005.
- All 44 chain x program tests are "not supported" (author and marker labels). The calibration shows the chain test cannot detect planted 1.1x-2x chains (power <= 0.04), so no calibrated range exists.
- Ranked candidate list: [`results/chain/ranked_hypotheses.csv`](results/chain/ranked_hypotheses.csv) (de-duplicated: `ranked_hypotheses_dedup.csv`). 0/88 hops supported, so all rows are unranked candidates, not evidence.
- Claims register: [`results/chain/phase7/claims_register.csv`](results/chain/phase7/claims_register.csv).
- Run order (env `thesis-spatial-ccc`): `chain_phase1_cell_states.py` (then `--amendment1`), `chain_phase2_neighbourhood.py`, `chain_phase2b_gnn.py`, `chain_phase4_simulation.py`, `chain_phase3_chain_test.py` (`--annotation marker` for the sensitivity analysis), `chain_phase5_candidates.py`, `chain_phase5b_dedup.py`, `chain_phase6_gse310352.py`, `chain_phase7_synthesis.py`.

## Model development (branch `model-v1`)
Folder [`model_dev/`](model_dev/README.md): baseline-stage models that predict a cell's expression state from its spatial neighbourhood (PDAC CosMx; Mendeley development set, GSE310352 external test). Report: [`model_dev/reports/MODEL_REPORT_v1.md`](model_dev/reports/MODEL_REPORT_v1.md); dashboard page: [`neighbourhood_models.html`](neighbourhood_models.html).
- Neighbours add +0.034 [0.018, 0.048] held-out R2 for the macrophage C1Q_TAM program (13 tumours, grouped CV); a GPU graph model (SAGE) is not better than simple aggregated baselines (+0.0025 [-0.013, 0.015]).
- Cross-type contribution inside the graph model: +0.019 [0.004, 0.028] for macrophages (CAF neighbours lead), +0.007 [-0.004, 0.024] for tumour cells.
- External test (single touch, 6 units of unverified independence): aggregated-neighbourhood gain transfers (XGBoost +0.072 [0.031, 0.107]); the graph model's does not (+0.014 [-0.014, 0.046]).
- Pretraining (masked expression, linear probe) does not help like for like.
- Environment: `model_dev/environment.yml` (CUDA torch, GPU training only); one command per phase in `model_dev/README.md`; 15 unit tests (`python -m pytest model_dev/tests -q`).

## Translational link (draft, separate track)

The prespecified Visium hypoxia robustness audit reproduced the original result but failed the robustness rule after tumour/host filtering and common-depth restriction. Effects remain positive; significance and paired sample size depend on the sensitivity. TCGA and CPTAC bulk PDAC cohorts are prepared with explicit histology restrictions and purity/cellularity provenance. Survival models are implemented and predeclared but **not run or validated on outcome data** at the user-requested PR hand-off. See [the draft translational summary](docs/TRANSLATIONAL_SUMMARY_v1.md), plans in `docs/translational_plans/`, and `results/translational/MANIFEST.csv`.

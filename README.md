# Thesis Prototype v0 — Spatial Cell–Cell Communication

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
- Step-by-step log with commands, outputs and decisions: `docs/PDAC_PROGRESS_LOG.md`
- Manifest of every produced file (path, phase, description, script): `results/pdac/MANIFEST.csv`
- Which breast number comes from which setting (15 µm vs radius sweep vs 6-NN): `docs/RESULTS_SETTINGS_NOTE.md`
- Code: `src/pdac/` (run logs in `logs/<script>_<timestamp>.log`; figures in `figures/pdac/<phase>/`)

Data are not committed. Download `raw_meta_data_final.h5ad` from the Mendeley dataset page into `data/pdac/mendeley_pdac_smi/`. COMMOT and SpatialDM were run in a separate conda env, `pdac-bench`: a clone of `thesis-spatial-ccc` plus `commot` 0.0.3 and `SpatialDM` 0.3.1, with numpy pinned below 2.

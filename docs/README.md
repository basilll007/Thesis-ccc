# Thesis Prototype v0 — Spatial Cell–Cell Communication

This prototype provides a reproducible baseline for cell–cell communication analysis in publicly available spatial transcriptomics data from breast cancer. It loads single-cell-resolution spatial expression data when available, performs basic quality control and normalization, clusters cells, builds a spatial neighbor graph, computes a baseline ligand–receptor analysis, and writes processed data, ranked tables, and figures. The prototype is infrastructure plus a baseline; functional-response evidence, reproducibility calibration, and Spatial Prize-Collecting Steiner Tree optimization are future work.

## Environment

The project targets Python 3.11. Create the environment with either conda/mamba from `environment.yml` or a Python virtual environment using the pinned `requirements.txt` with `uv pip`:

```powershell
conda env create -f environment.yml
conda activate thesis-spatial-ccc
# or, in a Python 3.11 venv:
uv venv .venv --python 3.11
.\.venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
```

## Data

Run the downloader in the requested order. It attempts the 10x Genomics Xenium breast-cancer public outputs, then GEO GSE243168, and only then uses a documented smoke-test fallback if both public sources cannot be downloaded:

```powershell
python src\download_data.py
```

The downloader records the exact source, date, size, and access terms in `docs\PROJECT_CONTEXT.md`. Files under `data\raw\` are never edited by the pipeline.

## Run the prototype

After environment setup and data download, run the complete pipeline with one command:

```powershell
python src\run_prototype.py
```

Outputs are written to `data\processed\`, `results\`, `figures\`, and `logs\`. The fixed random seed is 42. The pipeline logs cell/gene counts before and after QC, graph parameters, package versions, and the selected dataset.

## Verification

```powershell
python -c "import scanpy, squidpy, anndata, liana, numpy, scipy, pandas, networkx, matplotlib, seaborn, sklearn, leidenalg, igraph; print('imports OK')"
```

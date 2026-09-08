# Thesis Prototype v0 — Project Context

## Scientific Context

In spatial transcriptomics of tumors, most cell–cell communication (CCC) methods infer interactions from spatial proximity + ligand–receptor co-expression. These are plausible but not proof of a functional downstream response, and spatial autocorrelation inflates false positives, producing large, unstable interaction networks. The thesis will prioritize a **minimal, high-confidence set of communication signals** that are (a) supported by downstream functional-response evidence in receiver cells, (b) consistent with prior pathway knowledge, and (c) reproducible across patients, with calibrated uncertainty — eventually via a Spatial Prize-Collecting Steiner Tree (S-PCST) graph optimization. Perturbation data is only an optional external cross-check, not the core.

## Design Decisions

- Prefer single-cell-resolution spatial platforms (Xenium/CosMx) over spot-level Visium, because CCC needs cell-level resolution.
- Backbone data must be fully public (no controlled access) so nothing blocks us.
- v0 prototype = infrastructure + real data + a baseline model. The novel method (reproducibility calibration, S-PCST) is future work, **not** in this prototype.

## Glossary

- **AnnData**: Annotated data matrix container used for expression matrices and cell/gene metadata.
- **CCC**: Cell–cell communication; inferred signaling between sender and receiver cells.
- **Ligand–receptor (L–R) pair**: A signaling molecule and its cognate receptor used as a candidate interaction.
- **Spatial transcriptomics**: Transcriptome measurement retaining spatial coordinates.
- **Xenium**: 10x Genomics single-cell-resolution in situ transcriptomics platform.
- **CosMx**: NanoString single-cell/spatial molecular imaging platform.
- **Visium**: Spot-level spatial transcriptomics platform; each spot can contain multiple cells.
- **Leiden clustering**: Graph-based community detection used to group cells by expression profiles.
- **Spatial neighbor graph**: Graph connecting cells based on spatial proximity.
- **Baseline**: A simple reference analysis, not the proposed thesis method.
- **S-PCST**: Spatial Prize-Collecting Steiner Tree, planned future graph optimization.

## Data Provenance

- **Dataset**: 10x Genomics Xenium V1 Human Breast, 2 FOV public output bundle.
- **Source URL/accession**: https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip (no accession; official 10x public dataset).
- **Download date**: 2026-09-08 (UTC).
- **Size**: 379,018,013 bytes compressed; SHA-256 `cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b`.
- **License/access terms**: Creative Commons Attribution 4.0 International (CC BY 4.0), as stated on the official 10x dataset/support page; no access application required.

## Reproducibility

- Fixed random seed: 42.
- Raw data are stored under `data/raw/` and must not be edited.
- Processed AnnData objects are stored under `data/processed/`.
- Tables and metrics are stored under `results/`; figures under `figures/`; execution logs under `logs/`.

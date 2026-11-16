# model_dev: neighbourhood models of cell state (baseline stage, v1)

Predict a cell's expression state (program score, or expression over a declared gene set) from its spatial neighbourhood in PDAC CosMx data, and measure how much of the prediction comes from *other* cell types. These data measure expression; results are predictive associations, not causal or clinical claims.

- Plans (committed before each phase runs): `plans/`
- Report: `reports/MODEL_REPORT_v1.md`; data card and model card: `reports/`
- Progress log: `reports/PROGRESS_LOG.md`
- Small results and figures: `results/<phase>/`
- Large artefacts (datasets, checkpoints): `F:\Thesis\model_cache\` (git-ignored; paths and checksums in `results/phase0/artefacts.csv`)

## Environment
Conda env `thesis-model` (Python 3.11, CUDA torch). Pinned: `environment.yml`. Hardware: RTX 5060 Laptop (8 GB VRAM), 31 GB RAM.

```
conda env create -f model_dev/environment.yml
```

## One command per phase (run from the repository root)
| Phase | Command |
|---|---|
| 0 setup/data | `python model_dev/src/data.py --config model_dev/configs/data.yaml` |
| 1 targets/splits | `python model_dev/src/splits.py --config model_dev/configs/splits.yaml` |
| 2 baselines | `python model_dev/src/baselines.py --config model_dev/configs/baselines.yaml` |
| 3 graph models | `python model_dev/src/train.py --config model_dev/configs/graph.yaml` then `python model_dev/src/ablations.py --config model_dev/configs/ablations.yaml` |
| 4 external test | `python model_dev/src/evaluate.py --receiver myeloid` then `--receiver tumour` (each runs once) |
| 5 interpretation | `python model_dev/src/interpret.py --config model_dev/configs/interpret.yaml` |
| 6 pretraining | `python model_dev/src/pretrain.py --config model_dev/configs/pretrain.yaml` |
| 7 report | `python model_dev/src/report.py` |
| tests | `python -m pytest model_dev/tests -q` |

Every run writes `model_dev/results/<phase>/run_<timestamp>.log` with package versions, GPU model, seed and commit hash.

## Rules enforced in code
- Patient-level splits: all sections of one tumour share a fold (`results/phase1/folds.csv`, reused everywhere).
- No leakage: a cell's own expression (any gene in its target set or target program) is never an input; the feature registry carries provenance and `features.assert_no_target_leak` refuses offending matrices (see `tests/test_leakage.py`). Scalers, gene selection and hyperparameters are fitted on training folds only.
- External test set (GSE310352) is evaluated once, in Phase 4, with frozen preprocessing and a frozen model.

## Results (v1)
See `reports/MODEL_REPORT_v1.md` (lead question, all tables with CIs, ablations, external test, leakage tests, compute, limitations), `reports/MODEL_CARD.md`, `reports/DATA_CARD.md`, figures in `results/figures/`, and the dashboard page `neighbourhood_models.html` at the repository root.

## Running under memory pressure
All training is GPU-only. Long phases are resumable and can be run in foreground chunks:
`baselines.py --fold-budget N`, `ablations.py --time-budget S`, `pretrain.py --fold-budget N`, `learning_curve.py --time-budget S`; `train.py` resumes per run from `results/phase3/partial/` and checkpoints in `model_cache/checkpoints/`.

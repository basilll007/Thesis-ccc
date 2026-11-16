# Key findings across the thesis workstreams

One page per question, with the answer as stated in each workstream's own summary. Every number below is copied from that summary, which in turn cites the generated result files. All workstreams are expression-data analyses: results are audits, associations or predictive models, **not causal, mechanistic or clinical claims**. The patient (tumour) is the unit of inference throughout.

| Workstream | Question | Answer (short) | Full summary |
|---|---|---|---|
| Breast Xenium prototype | Do cluster-mean CCC scores reflect spatial contact? | The dashboard's calibrated edge-based score with a 500-permutation spatial null flags 9 of 17 previously significant interactions as not supported by spatial contacts. | `index.html`, `docs/PROJECT_CONTEXT.md` |
| PDAC v2 | Does a proximity-specific signal exist at tumour level, or does the edge score mainly mirror cell-type products? | Mainly mirroring. CXCL12->CXCR4 CAF->T ratio Q is BH-significant in 0/16 sections and 0/13 tumours; only a weak exploratory trend (z_Q > 0 in 10/13 tumours). | `docs/PDAC_SUMMARY_v2.md` |
| PDAC v3 | Does that signal survive spatially constrained nulls, and could the audit detect a real effect? | It does not shrink under block and distance nulls, but calibration shows zero power at about 25 T cells, and most nulls are fooled by 200 um zonation (42-88 % false positives). GSE310352: Q > 1 in 8/8 slides, BH-significant in at most 1. GSE277783 Visium: a hypoxia program is higher in metastases (paired p = 0.00244, q = 0.0122, 13 patients). | `docs/PDAC_SUMMARY_v3.md` |
| Translational v1 (draft) | Is the v3 hypoxia association robust? Do spatial signatures relate to survival? | **Not robust** under the prespecified rule: the tumour-rich filter gives p = 0.068 and common depth p = 0.193. All directions remain positive. TCGA (141 patients) and CPTAC (129) cohorts are prepared; **survival analysis was not run**. | `docs/TRANSLATIONAL_SUMMARY_v1.md` |
| Chain v1 | Do neighbours predict cell state, and is any A -> B -> C signalling chain supported? | Neighbours predict state almost only through **same-type** neighbours (tumour classical ridge delta R2 0.141; other types add <= 0.0005). **44/44 chains not supported**, and the chain test has power <= 0.04 even for 2x planted chains, so it could not have found them. | `docs/CHAIN_SUMMARY_v1.md` |
| Model development v1 | Does a graph model predict macrophage and tumour state from neighbours better than simple baselines, and how much comes from other cell types? | Neighbours add +0.034 [0.018, 0.048] R2 for macrophage C1Q_TAM. A GPU graph model is **not better** than simple aggregated baselines (+0.0025 [-0.013, 0.015]). The other-cell-type share is +0.019 [0.004, 0.028] for macrophages (CAFs lead) and +0.007 [-0.004, 0.024] for tumour cells. The aggregated-baseline gain transfers to GSE310352 (+0.072 [0.031, 0.107]); the graph model's does not. Pretraining did not help. | `model_dev/reports/MODEL_REPORT_v1.md`, dashboard `neighbourhood_models.html` |

## What holds across the workstreams
1. **Composition and zonation explain most "communication" signal.** Edge-based CCC scores mostly mirror cell-type products (v2). Constrained nulls are fooled by zonation (v3). Neighbourhood prediction is driven by same-type neighbours (chain, model v1).
2. **Power is the binding constraint.** The Q audit (v3) and the chain test (chain v1) were both calibrated by simulation and found unable to detect effects of plausible size at the available cell counts. The model learning curve suggests about +0.06 R2 per tenfold more cells.
3. **The only consistent cross-type lead is CAF <-> macrophage.** It appears as an exploratory adjacency association (chain claim C10), as the largest other-type contribution to macrophage-state prediction (model v1), and as a model-based hypothesis (CXCL12-high CAFs <-> C1Q/SELENOP-type macrophages). All three come from the same kind of data. It is a candidate for experimental testing, not a finding.
4. **The external cohorts support consistency, not replication.** GSE310352 patient identity, independence and um/px scale are unverified.

## Status
- **Branches:** `pdac-v2`, `thesis-v3`, `chain-v1`, `model-v1` and `translational-v1` are all merged into `main`.
- **Translational work:** a draft. Survival models are implemented and pre-declared but not run.
- **Remaining limits:** segmentation spillover is not modelled anywhere, and T-cell states are largely below the CosMx background.

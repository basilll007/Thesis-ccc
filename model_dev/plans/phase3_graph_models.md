# Phase 3 plan: graph models and ablations

Committed before `train.py` / `ablations.py` are run. Folds, targets, metrics and success criteria are those of the Phase 1 plan; baselines are from Phase 2.

## Implementation choice (logged reading of the brief)
PyTorch Geometric 2.8.1 installed cleanly, but the models are implemented as plain padded-tree message passing in PyTorch. Reason: the chain Phase 2b bug was a return path (i -> j -> i) that put a cell's own expression into its inputs. In a padded 2-hop tree the seed is removed from its own second hop by construction (`graphs.two_hop_batch`, guarded by `tests/test_leakage.py`), and ablation masks are explicit index masks rather than graph rewiring. PyG stays available for a later stage.

## Inputs
- **Seed (receiver) cell:** cell-type embedding (all annotation classes; 16 dims), the own covariates (log1p total counts, log1p detected genes, density within 30 um, mean 6-NN edge length) and the tile block (leave-one-out 200 um tile means of the targets over same-type cells, empty indicator). **Never its own expression.** This is the same information baseline (2) sees.
- **Neighbour cells (hop 1 and hop 2):** cell-type embedding, log1p total counts, log1p detected genes, and their 539-gene log expression (G), standardised per gene with training-fold statistics. Neighbour counts/genes are other cells' values.
- **Edges:** log1p edge length (um; um-equivalent for GSE310352) as an edge feature at both hops.
- Neighbourhood: symmetrised 6-NN graph; k1 = k2 = 12 nearest neighbours kept per node (the fraction of nodes with degree > 12 is reported). Seed excluded at hop 2.

## Architectures
- **SAGE-style:** hop-2 node state h_j = ReLU(W_self x_j + W_nb * mean_k phi(x_k, e_jk)); seed state h_i = ReLU(W_s c_i + W_n * mean_j psi(h_j, x_j, e_ij)); linear head to all target columns of the receiver (multi-output). Hidden 128, dropout 0.1.
- **Attention-style:** identical, but each mean over neighbours is replaced by 4-head dot-product attention with the receiving node's state as query and the neighbour messages (with edge feature) as keys/values, masked over padding. Attention weights are *not* reported as explanations (Phase 5).
- Loss: mean squared error over target columns. Mixed precision (torch autocast, bf16) on the GPU; batches of 512 seeds; if memory fails the batch is halved and logged; CPU fallback is logged if the GPU is unavailable.

## Training and selection
- Train on the inner training folds, early stopping on the inner validation fold (patience 8, max 60 epochs), evaluate on the outer test fold; no refit (stated with the results, as for the MLP).
- Grid (inner validation only): learning rate in {1e-3, 3e-4}; weight decay 1e-4; everything else fixed. Seeds 0, 1, 2 for every configuration. Every configuration x fold x seed is written to `graph_results.csv`; training curves to `curves/`.
- Architecture selection for ablations: the architecture with the lower mean inner-validation MSE on the primary receiver (myeloid) at its selected learning rate.

## Ablations (same folds, same seeds, selected architecture and learning rate)
(i) same-type neighbours only: all tree nodes not of the receiver's class are masked (hop 1 and hop 2);
(ii) other-type neighbours only: all tree nodes of the receiver's class are masked;
(iii) composition only: neighbour gene expression zeroed (type embedding, counts covariates and distances kept);
(iv) zonation control: each neighbour's expression (and counts) replaced by those of a random cell of the same class in the same 200 um tile (one fixed permutation per seed and section);
(v) no graph: all neighbours masked (seed covariates only; the graph model reduces to an MLP on the baseline (2) inputs).
Scope (compute bound, declared here): all five ablations for the primary receiver (myeloid); (i), (ii) and (v) for the positive control (tumour). 3 seeds each.

## Quantities and decision rules
- **Beats the strongest baseline:** R2(graph, seed-mean) - R2(strongest Phase 2 baseline for that target) with a paired tumour-bootstrap 95 % CI excluding 0. Otherwise the graph model is reported as "not better than baselines".
- **Cross-type contribution** = R2(full) - R2(same-type only), paired tumour bootstrap CI. "Other cell types carry signal" only if the CI excludes 0. Also reported: R2(full) - R2(other-type only), R2(full) - R2(composition only), R2(full) - R2(shuffled), and R2(no graph) versus baseline (2) as an implementation check (they should be close).
- Primary endpoint remains myeloid:C1Q_TAM.

## Outputs
`model_dev/results/phase3/`: `graph_results.csv`, `graph_summary.csv` (selected configurations with CIs and comparison to the strongest baseline), `ablation_results.csv`, `ablation_summary.csv` (cross-type contribution with CI), `curves/`, `compute.csv` (wall time, peak VRAM, device), figure. Checkpoints go to `model_cache/checkpoints/` (registered with checksums).

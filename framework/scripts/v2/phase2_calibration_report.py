"""v2 Phase 2: rejection counts of the v2 contact tests and the v1 label-permutation Q test on the
200 null datasets of tests/test_calibration.py (same seeds, 1000 permutations)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binom
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext, contact_statistics
from cellalibi.simulate import synthetic_tissue
from cellalibi.graphs import knn_graph
from cellalibi.nulls import global_null, block_null
from cellalibi.scores import prepare_lr_setup, compute_scores_for_labels
import anndata as ad

rows = []
for s in range(200):
    xy, labels, counts, genes = synthetic_tissue(seed=s)
    ctx = SpatialContext(xy, labels, AuditConfig(n_perms=1000, tile_um=100.0, seed=s))
    ex = ExpressionData(counts, genes)
    m = contact_statistics(ctx, ex, "Sender", "Receiver", "L", "R", cfg=ctx.cfg, seed=1000 + s, run_spillover=False)
    # v1 label-permutation Q test (thesis/v1 engine) on the same data, log1p CP10K
    tot = counts.sum(1, keepdims=True); X = np.log1p(counts / tot * 1e4)
    a = ad.AnnData(X=X); a.var_names = genes
    r, c = knn_graph(xy, k=6)
    types = sorted(set(labels)); lab = np.array([types.index(x) for x in labels])
    setup = prepare_lr_setup(a, [("L", "R")], r, c, use_raw=False)
    si, ri = types.index("Sender"), types.index("Receiver")
    obs = compute_scores_for_labels(lab, len(types), setup)["ratio"][0, si, ri]
    res = {}
    for name, perms in (("v1_label_global", global_null(lab, 1000, seed=s)),
                        ("v1_label_block100", block_null(lab, xy, 100.0, 1000, seed=s)[0])):
        nullq = np.array([compute_scores_for_labels(p, len(types), setup, zero_fill_edge=True)["ratio"][0, si, ri] for p in perms])
        res[name] = (np.sum(np.nan_to_num(nullq, nan=-np.inf) >= obs) + 1) / 1001
    rows.append({"seed": s, "v2_global": m["p_contact_global"], "v2_block_contrast": m["p_block_contrast"],
                 "v2_ring": m["p_ring"], "v2_block_q_info": m["p_contact_block"], **res})
df = pd.DataFrame(rows)
od = Path(__file__).resolve().parents[2] / "results/v2_phase2"
df.to_csv(od / "null_calibration_pvalues.csv", index=False)
lo, hi = binom.ppf([0.025, 0.975], 200, 0.05)
summ = pd.DataFrame([{"test": c, "n_datasets": 200, "n_rejections_p_lt_0.05": int((df[c] < 0.05).sum()),
                      "rejection_rate": float((df[c] < 0.05).mean()), "binom95_lo": int(lo), "binom95_hi": int(hi),
                      "inside_interval": bool(lo <= (df[c] < 0.05).sum() <= hi)} for c in df.columns if c != "seed"])
summ.to_csv(od / "null_calibration_summary.csv", index=False)
print(summ.to_string(index=False))

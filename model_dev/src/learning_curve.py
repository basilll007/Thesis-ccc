"""Phase 6 (data-scale estimate): learning curves over the number of training tumours.

What (plan: model_dev/plans/phase6_pretraining.md): for each outer fold, train the selected supervised graph
model on nested subsets of the inner-training tumours (3, 6, all; 9 where available) with one seed, evaluate on
the outer test fold, and fit a descriptive log-linear curve of pooled R2 against the number of training cells.
The probe curve re-fits the ridge probe on the same tumour subsets using the fold's pretrained encoder (the
encoder itself saw all inner-training cells; this is logged as a simplification). The extrapolation is rough.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import assemble as A  # noqa: E402
import features as F  # noqa: E402
import train as T  # noqa: E402
from baselines import pooled_r2, sse_sst_by_group  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

PHASE = "phase6"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "graph.yaml"))
    args = ap.parse_args()
    gcfg = C.load_config(args.config)
    log = C.start_run(PHASE, "learning_curve", gcfg)
    sel = json.loads((C.RESULTS / "phase3" / "architecture_selection.json").read_text())
    arch, lr = sel["selected_arch"], float(sel["selected_lr"])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    fold_table = A.folds("grouped"); n_folds = int(fold_table.fold.max()) + 1
    spec = A.targets_spec()
    store0, offsets, classes = T.build_store("development", max(gcfg["k1"], gcfg["k2"]), device, log)
    rows = []
    for receiver in gcfg["receivers"]:
        s = spec[receiver]
        data = A.build_receiver("development", receiver, log)
        rc = classes.index(s["class"])
        seed_global = np.array([offsets[x] + r for x, r in zip(data.section, data.row)])
        raw = data.matrix(["own", "tile"])
        j = data.target_names.index(s["programs"][0])   # first program of the receiver (C1Q_TAM / classical)
        for k in range(n_folds):
            tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
            train_secs = set(fold_table[fold_table.fold != k].section)
            gmask = torch.as_tensor(np.isin(store0.section_of, list(train_secs))).to(device)
            mu, sd, cmu, csd = T.fold_standardisation(store0, gmask)
            store = T.StandardisedStore(store0, mu, sd, cmu, csd)
            tumours = sorted(set(data.group[tr_in]))
            rng = np.random.default_rng([C.SEED, 6, k])
            order = list(rng.permutation(tumours))
            sizes = sorted({n for n in (3, 6, 9, len(order)) if n <= len(order)})
            ytr_mean = data.Y[tr_all].mean(axis=0)
            for n in sizes:
                sub = tr_in & np.isin(data.group, order[:n])
                smu, ssd = F.fit_scaler(raw[sub]); cov = F.apply_scaler(raw, smu, ssd)
                pred, hist, ep, bv, peak, _ = T.train_one(store, seed_global, cov, data.Y, (sub, val, test), rc, "full", arch, lr, 0,
                                                          gcfg, len(classes), device, log)
                p = sse_sst_by_group(data.Y[test][:, j:j + 1], pred[:, j:j + 1], data.group[test], ytr_mean[j:j + 1])
                rows.append({"receiver": receiver, "target": data.target_names[j], "fold": k, "n_train_tumours": n,
                             "n_train_cells": int(sub.sum()), "model": f"graph/{arch}", "test_r2_fold": pooled_r2(p), "epochs": ep})
                log.info("%s fold %d n_tumours=%d (%d cells): R2 %.4f", receiver, k, n, int(sub.sum()), pooled_r2(p))
    df = pd.DataFrame(rows)
    C.save_table(df, PHASE, "learning_curve_runs.csv", log)
    # descriptive log-linear fit per receiver: R2 ~ a + b * log10(n_cells), averaged over folds at each size index
    out = []
    for (receiver, t), d in df.groupby(["receiver", "target"]):
        agg = d.groupby("n_train_tumours").agg(mean_r2=("test_r2_fold", "mean"), mean_cells=("n_train_cells", "mean"), n_folds=("fold", "nunique")).reset_index()
        b, a = np.polyfit(np.log10(agg.mean_cells), agg.mean_r2, 1)
        for r in agg.itertuples():
            out.append({"receiver": receiver, "target": t, "n_train_tumours": r.n_train_tumours, "mean_train_cells": r.mean_cells,
                        "mean_r2": r.mean_r2, "n_folds": r.n_folds, "loglinear_slope_per_decade": b,
                        "r2_at_10x_cells_extrapolated": a + b * np.log10(agg.mean_cells.max() * 10),
                        "note": "descriptive log-linear extrapolation; rough"})
    C.save_table(pd.DataFrame(out), PHASE, "learning_curve.csv", log)
    log.info("done")


if __name__ == "__main__":
    main()

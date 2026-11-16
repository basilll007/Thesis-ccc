"""Phase 4: the single external test on GSE310352 with frozen preprocessing and frozen models.

What (plan: model_dev/plans/phase4_external_test.md): retrain the declared model list once on the full
development set with hyperparameters fixed from the internal cross-validation, then evaluate every model once on
the external receivers (marker labels). Reports pooled R2 with unit-bootstrap CIs (6 units, UNVERIFIED grouping and
independence), per-slide breakdown, the internal-to-external drop and the gain over frozen own covariates.

Why: the external set is touched exactly once, by this script, and nothing in it is tuned here.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import assemble as A  # noqa: E402
import features as F  # noqa: E402
import models as M  # noqa: E402
import train as T  # noqa: E402
from baselines import FEATURE_SETS, merge_pairs, pooled_r2, sse_sst_by_group  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

PHASE = "phase4"


def modal_config(results: pd.DataFrame, receiver: str, model: str, fset: str, key: str = "selected") -> str:
    d = results[(results.receiver == receiver) & (results.model == model) & (results.feature_set == fset)]
    d = d.drop_duplicates(["fold", "config"])
    if "selected" in d and d.selected.astype(str).eq("True").any():
        d = d[d.selected.astype(str) == "True"]
    else:  # choose by mean inner-val MSE over folds
        best = d.groupby("config").inner_val_mse.mean().idxmin()
        return best
    counts = d.config.value_counts()
    top = counts[counts == counts.max()].index
    return sorted(top, key=lambda c: float(c.split("=")[-1]))[-1]


def ridge_alpha_from(config: str) -> float:
    return float(config.split("alpha=")[1])


def evaluate_models(receiver: str, cfg: dict, log, device: str) -> tuple[list[dict], dict]:
    spec = A.targets_spec()[receiver]
    dev = A.build_receiver("development", receiver, log)
    ext = A.build_receiver("external", receiver, log)
    F.assert_no_target_leak(ext.feature_set(list(ext.X)), None, ext.Y, set(spec["genes"]), set(spec["programs"]))
    # external classes may differ from development classes: align composition columns by name
    for b in ("composition",):
        if ext.names[b] != dev.names[b]:
            idx = {n: i for i, n in enumerate(ext.names[b])}
            cols = np.column_stack([ext.X[b][:, idx[n]] if n in idx else np.zeros(len(ext.Y), np.float32) for n in dev.names[b]])
            ext.X[b], ext.names[b], ext.prov[b] = cols, dev.names[b], dev.prov[b]
            log.info("%s: external composition columns aligned to development classes", receiver)
    tnames, n_prog = dev.target_names, len(spec["programs"])
    ymean = dev.Y.mean(axis=0)                       # frozen mean predictor
    p2 = pd.read_csv(C.RESULTS / "phase2" / "baseline_results.csv")
    best = pd.read_csv(C.RESULTS / "phase2" / "strongest_baseline_per_target.csv")
    rows, pairs = [], {}

    def record(model_name, pred, internal_r2_lookup):
        for j, t in enumerate(tnames + (["genes_pooled"] if len(tnames) > n_prog else [])):
            if t == "genes_pooled":
                yj, pj, mj = ext.Y[:, n_prog:], pred[:, n_prog:], ymean[n_prog:]
            else:
                yj, pj, mj = ext.Y[:, j:j + 1], pred[:, j:j + 1], ymean[j:j + 1]
            pu = sse_sst_by_group(yj, pj, ext.group, mj)
            ps = sse_sst_by_group(yj, pj, ext.section, mj)
            pairs[(model_name, t)] = pu
            r2, lo, hi = C.tumour_bootstrap(pu, cfg["n_boot"])
            rows.append({"receiver": receiver, "model": model_name, "target": t, "external_r2": r2, "unit_ci_low": lo,
                         "unit_ci_high": hi, "n_units": len(pu), "internal_r2": internal_r2_lookup(t),
                         "per_slide_r2": json.dumps({k: round(1 - v[0] / v[1], 4) for k, v in ps.items()}),
                         "n_external_cells": len(ext.Y), "note": "units = slides with 5/5b2/5b3 grouped (UNVERIFIED); not independent patients"})

    # 1. mean predictor
    record("mean", np.tile(ymean, (len(ext.Y), 1)), lambda t: 0.0)
    summ2 = pd.read_csv(C.RESULTS / "phase2" / "baseline_summary.csv")

    def internal(model, fset):
        def f(t):
            d = summ2[(summ2.receiver == receiver) & (summ2.model == model) & (summ2.feature_set == fset) & (summ2.target == t)]
            return float(d.r2.iloc[0]) if len(d) else np.nan
        return f

    # 2./3. ridge (own covariates and the strongest baselines' feature sets), XGBoost-CUDA and MLP - all on the GPU
    import gpu as GP
    for fset in ["own"] + sorted({f for f in best[best.receiver == receiver].feature_set.unique() if f != "none"}):
        blocks = FEATURE_SETS[fset]
        Xd, Xe = GP.to_gpu(dev.matrix(blocks)), GP.to_gpu(ext.matrix(blocks))
        mu, sd = GP.fit_scaler(Xd)                                   # development statistics only (frozen)
        Xd, Xe = GP.apply_scaler(Xd, mu, sd), GP.apply_scaler(Xe, mu, sd)
        Ydg = GP.to_gpu(dev.Y)
        alpha = ridge_alpha_from(modal_config(p2, receiver, "ridge", fset))
        rp = GP.RidgePath(Xd, Ydg)
        record(f"ridge/{fset}", rp.predict(Xe, alpha).cpu().numpy(), internal("ridge", fset))
        del rp
        torch.cuda.empty_cache()
        log.info("%s ridge/%s alpha=%g evaluated", receiver, fset, alpha)
        for b in best[(best.receiver == receiver) & (best.feature_set == fset)].itertuples():
            if b.model == "xgb" and (f"xgb/{fset}", tnames[0]) not in pairs:
                cfgname = modal_config(p2, receiver, "xgb", fset)
                lr = float(cfgname.split("lr=")[1].split(",")[0]); leaves = int(cfgname.split("leaves=")[1])
                d = p2[(p2.receiver == receiver) & (p2.model == "xgb") & (p2.feature_set == fset) & (p2.config == cfgname)]
                rows_tr = np.random.default_rng(C.SEED).choice(len(dev.Y), min(len(dev.Y), cfg["hgb_max_rows"]), replace=False)
                Xd_np = Xd[torch.as_tensor(rows_tr, device="cuda")].cpu().numpy(); Xe_np = Xe.cpu().numpy()
                pred = np.tile(ymean, (len(ext.Y), 1)).astype(np.float32)
                for j in range(n_prog):
                    iters = int(np.median([int(x.split(";")[j]) for x in d.drop_duplicates("fold").best_iters]))
                    pred[:, j] = GP.xgb_fit_predict(Xd_np, dev.Y[rows_tr, j], Xe_np, lr, leaves, iters, C.SEED)
                record(f"xgb/{fset}", pred, internal("xgb", fset))
                del Xd_np, Xe_np
            if b.model == "mlp" and (f"mlp/{fset}", tnames[0]) not in pairs:
                d = p2[(p2.receiver == receiver) & (p2.model == "mlp") & (p2.feature_set == fset)]
                lr = float(modal_config(p2, receiver, "mlp", fset).split("lr=")[1].split(",")[0])
                epochs = int(np.median(d.drop_duplicates(["fold", "config"]).epochs))
                preds = [GP.GPUMLP(Xd.shape[1], len(tnames), (256, 128), 0.1, lr, 1e-4, seed, 1024).fit_fixed(Xd, Ydg, epochs)
                         .predict(Xe).cpu().numpy() for seed in cfg["seeds"]]
                record(f"mlp/{fset}", np.mean(preds, axis=0), internal("mlp", fset))
        del Xd, Xe, Ydg
        torch.cuda.empty_cache()
    # 4. selected graph model
    sel = json.loads((C.RESULTS / "phase3" / "architecture_selection.json").read_text())
    gcfg = C.load_config(C.MODEL_DEV / "configs" / "graph.yaml")
    comp = pd.read_csv(C.RESULTS / "phase3" / "compute_graph.csv")
    comp = comp[(comp.receiver == receiver) & (comp.arch == sel["selected_arch"]) & (comp.lr == float(sel["selected_lr"])) & (comp["mode"] == "full")]
    epochs = int(np.median(comp.epochs))
    gcfg = {**gcfg, "max_epochs": epochs, "patience": 10**9}
    summ3 = pd.read_csv(C.RESULTS / "phase3" / "graph_summary.csv")
    preds = []
    for seed in cfg["seeds"]:
        store_d, off_d, classes = T.build_store("development", max(gcfg["k1"], gcfg["k2"]), device, log)
        gmask = torch.ones(store_d.n, dtype=torch.bool, device=device)
        mu, sd, cmu, csd = T.fold_standardisation(store_d, gmask)
        std_d = T.StandardisedStore(store_d, mu, sd, cmu, csd)
        raw_d = dev.matrix(["own", "tile"]); smu, ssd = F.fit_scaler(raw_d)
        g_d = np.array([off_d[s] + r for s, r in zip(dev.section, dev.row)])
        rc = classes.index(A.targets_spec()[receiver]["class"])
        tr = np.ones(len(dev.Y), bool); val = np.zeros(len(dev.Y), bool); val[:2048] = True
        _, hist, ep, bv, peak, state = T.train_one(std_d, g_d, F.apply_scaler(raw_d, smu, ssd), dev.Y, (tr, val, val), rc, "full",
                                                   sel["selected_arch"], float(sel["selected_lr"]), seed, gcfg, len(classes), device, log)
        del store_d, std_d; torch.cuda.empty_cache()
        # external store: classes mapped by NAME onto the development class list (unknown names -> 'unassigned')
        store_e, off_e, classes_e = T.build_store("external", max(gcfg["k1"], gcfg["k2"]), device, log)
        remap = torch.as_tensor([classes.index(c) if c in classes else classes.index("unassigned") for c in classes_e]).to(device)
        store_e.cls = remap[store_e.cls]
        std_e = T.StandardisedStore(store_e, mu, sd, cmu, csd)       # development statistics (frozen)
        model = M.TreeModel(std_e.X.shape[1], len(classes), raw_d.shape[1], len(tnames), gcfg["hidden"], gcfg["emb"],
                            "attention" if sel["selected_arch"] == "attention" else "mean", gcfg["dropout"], gcfg["heads"]).to(device)
        model.load_state_dict(state); model.eval()
        g_e = torch.as_tensor(np.array([off_e[s] + r for s, r in zip(ext.section, ext.row)])).long().to(device)
        cov_e = torch.as_tensor(F.apply_scaler(ext.matrix(["own", "tile"]), smu, ssd)).to(device)
        out = []
        with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device == "cuda"):
            for i in range(0, len(g_e), 1024):
                tb = M.TreeBatch(std_e, g_e[i:i + 1024], cov_e[i:i + 1024], rc, "full", gcfg["k1"], gcfg["k2"])
                out.append(model(tb).float().cpu().numpy())
        preds.append(np.vstack(out))
        del store_e, std_e; torch.cuda.empty_cache()
        log.info("%s graph %s seed %d: trained %d epochs on all development cells, evaluated externally", receiver, sel["selected_arch"], seed, ep)
    def internal_graph(t):
        d = summ3[(summ3.receiver == receiver) & (summ3.arch == sel["selected_arch"]) & (summ3.target == t)]
        return float(d.r2.iloc[0]) if len(d) else np.nan
    record(f"graph/{sel['selected_arch']}", np.mean(preds, axis=0), internal_graph)
    # gains over frozen own covariates (paired unit bootstrap)
    for r in rows:
        if r["model"] != "ridge/own" and (("ridge/own", r["target"]) in pairs):
            g, lo, hi = C.tumour_bootstrap_diff(pairs[(r["model"], r["target"])], pairs[("ridge/own", r["target"])], cfg["n_boot"])  # R2(model) - R2(own)
            r.update({"gain_vs_own_external": g, "gain_unit_ci_low": lo, "gain_unit_ci_high": hi, "transfers": bool(lo > 0)})
        r["internal_to_external_drop"] = r["internal_r2"] - r["external_r2"] if r["internal_r2"] == r["internal_r2"] else np.nan
    return rows, pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "external.yaml"))
    ap.add_argument("--receiver", required=True, help="one receiver per call (keeps each call short); each is evaluated once")
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    out = C.RESULTS / PHASE / f"external_results_{args.receiver}.csv"
    if out.exists():
        raise FileExistsError(f"{C.rel(out)} exists: the external set is evaluated once per receiver; delete deliberately and log it to re-run")
    log = C.start_run(PHASE, f"external_{args.receiver}", cfg)
    log.info("EXTERNAL TEST: single touch for %s at commit %s", args.receiver, C.git_commit())
    if not torch.cuda.is_available():
        raise RuntimeError("GPU-only policy: CUDA not available")
    rows, _ = evaluate_models(args.receiver, cfg, log, "cuda")
    df = pd.DataFrame(rows)
    C.save_table(df, PHASE, out.name, log)
    C.save_json({"receiver": args.receiver, "commit": C.git_commit(), "evaluated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                 "models": sorted(df.model.unique()), "note": "single touch; see plan phase4_external_test.md"},
                PHASE, f"external_touch_{args.receiver}.json", log)
    parts = [C.RESULTS / PHASE / f"external_results_{r}.csv" for r in cfg["receivers"]]
    if all(p.exists() for p in parts):
        C.save_table(pd.concat([pd.read_csv(p) for p in parts], ignore_index=True), PHASE, "external_results.csv", log)
    cols = [c for c in ["receiver", "model", "target", "internal_r2", "external_r2", "unit_ci_low", "unit_ci_high",
                        "gain_vs_own_external", "transfers"] if c in df]
    log.info("external results:
%s", df[cols].round(4).to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()

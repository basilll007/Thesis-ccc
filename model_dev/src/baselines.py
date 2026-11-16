"""Phase 2: baselines on identical patient-level folds, plus the chain Phase 2 reproduction check.

Models (plan: model_dev/plans/phase2_baselines.md): mean predictor; own-covariates ridge; ridge and gradient
boosting on aggregated neighbourhood features (`programs` and `genes` sets); torch MLP on `genes`.
Hyperparameters come from a small pre-declared grid evaluated on the inner validation fold only; every
configuration tried is written to baseline_results.csv.
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

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase2"
FEATURE_SETS = {"own": ["own", "tile"],
                "programs": ["own", "tile", "composition", "nb_programs"],
                "genes": ["own", "tile", "composition", "nb_programs", "nb_genes"]}


# ----------------------------------------------------------------------------- metrics
def sse_sst_by_group(y: np.ndarray, pred: np.ndarray, group: np.ndarray, y_train_mean: np.ndarray) -> dict:
    """Per-group (SSE, SST) pooled over target columns; SST uses the training-fold mean."""
    out = {}
    for g in np.unique(group):
        m = group == g
        out[str(g)] = (float(((y[m] - pred[m]) ** 2).sum()), float(((y[m] - y_train_mean) ** 2).sum()))
    return out


def pooled_r2(pairs: dict) -> float:
    sse = sum(v[0] for v in pairs.values())
    sst = sum(v[1] for v in pairs.values())
    return 1 - sse / sst


def merge_pairs(list_of_pairs: list[dict]) -> dict:
    out = {}
    for p in list_of_pairs:
        for k, v in p.items():
            a, b = out.get(k, (0.0, 0.0))
            out[k] = (a + v[0], b + v[1])
    return out


# ----------------------------------------------------------------------------- models
def fit_ridge(Xtr, ytr, alpha: float):
    from sklearn.linear_model import Ridge
    return Ridge(alpha=alpha).fit(Xtr, ytr)


def fit_hgb_staged(Xtr, ytr, Xval, yval, lr: float, leaves: int, max_iter: int, seed: int, max_rows: int):
    """One HGB fit per grid point; the iteration count is chosen on the inner validation loss (staged predictions)."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    rng = np.random.default_rng(seed)
    if len(ytr) > max_rows:
        k = rng.choice(len(ytr), max_rows, replace=False)
        Xtr, ytr = Xtr[k], ytr[k]
    m = HistGradientBoostingRegressor(learning_rate=lr, max_leaf_nodes=leaves, max_iter=max_iter, min_samples_leaf=50,
                                      early_stopping=False, random_state=seed).fit(Xtr, ytr)
    losses = [float(((yval - p) ** 2).mean()) for p in m.staged_predict(Xval)]
    best_iter = int(np.argmin(losses)) + 1
    return m, best_iter, losses[best_iter - 1]


class TorchMLP:
    """Small multi-output MLP with early stopping on the inner validation fold (torch, GPU if available)."""

    def __init__(self, n_in: int, n_out: int, hidden=(256, 128), dropout=0.1, lr=1e-3, wd=1e-4, seed=0,
                 batch=1024, max_epochs=100, patience=10):
        import torch
        import torch.nn as nn
        C.set_seed(seed)
        layers, d = [], n_in
        for h in hidden:
            layers += [nn.Linear(d, h), nn.ReLU(), nn.Dropout(dropout)]
            d = h
        layers.append(nn.Linear(d, n_out))
        self.net = nn.Sequential(*layers)
        self.dev = "cuda" if torch.cuda.is_available() else "cpu"
        self.net.to(self.dev)
        self.lr, self.wd, self.batch, self.max_epochs, self.patience = lr, wd, batch, max_epochs, patience
        self.history: list[dict] = []
        self.seed = seed

    def fit(self, Xtr, ytr, Xval, yval):
        import torch
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=self.wd)
        Xt, yt = torch.as_tensor(Xtr), torch.as_tensor(ytr)
        Xv, yv = torch.as_tensor(Xval).to(self.dev), torch.as_tensor(yval).to(self.dev)
        g = torch.Generator().manual_seed(self.seed)
        best, best_state, bad = np.inf, None, 0
        for epoch in range(self.max_epochs):
            self.net.train()
            perm = torch.randperm(len(Xt), generator=g)
            tot = 0.0
            for i in range(0, len(perm), self.batch):
                idx = perm[i:i + self.batch]
                xb, yb = Xt[idx].to(self.dev), yt[idx].to(self.dev)
                opt.zero_grad()
                loss = ((self.net(xb) - yb) ** 2).mean()
                loss.backward()
                opt.step()
                tot += float(loss) * len(idx)
            self.net.eval()
            with torch.no_grad():
                vl = float(((self.net(Xv) - yv) ** 2).mean())
            self.history.append({"epoch": epoch, "train_mse": tot / len(Xt), "val_mse": vl})
            if vl < best - 1e-6:
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in self.net.state_dict().items()}
            else:
                bad += 1
                if bad >= self.patience:
                    break
        self.net.load_state_dict(best_state)
        self.best_val, self.epochs = best, len(self.history)
        return self

    def predict(self, X):
        import torch
        self.net.eval()
        out = []
        with torch.no_grad():
            for i in range(0, len(X), 8192):
                out.append(self.net(torch.as_tensor(X[i:i + 8192]).to(self.dev)).cpu().numpy())
        return np.vstack(out)


# ----------------------------------------------------------------------------- evaluation loop
def run_receiver(receiver: str, data: A.ReceiverData, fold_table: pd.DataFrame, cfg: dict, log) -> tuple[list, dict]:
    n_folds = int(fold_table.fold.max()) + 1
    tnames = data.target_names
    spec = A.targets_spec()[receiver]
    n_prog = len(spec["programs"])
    results, pairs_store = [], {}  # pairs_store[(model, fset, config, target)] -> merged per-tumour pairs

    def record(model, fset, config, fold, pred_test, test, ytr_mean, inner_val_mse, extra=None):
        yte = data.Y[test]
        for j, t in enumerate(tnames):
            pairs = sse_sst_by_group(yte[:, j:j + 1], pred_test[:, j:j + 1], data.group[test], ytr_mean[j:j + 1])
            key = (model, fset, config, t)
            pairs_store[key] = merge_pairs([pairs_store.get(key, {}), pairs])
            results.append({"receiver": receiver, "model": model, "feature_set": fset, "config": config, "fold": fold,
                            "target": t, "target_kind": "program" if j < n_prog else "gene",
                            "n_test": int(test.sum()), "inner_val_mse": inner_val_mse, "test_r2_fold": pooled_r2(pairs),
                            **(extra or {})})
        # pooled over gene columns (the gene-vector target)
        if len(tnames) > n_prog:
            gj = slice(n_prog, None)
            pairs = sse_sst_by_group(yte[:, gj], pred_test[:, gj], data.group[test], ytr_mean[gj])
            key = (model, fset, config, "genes_pooled")
            pairs_store[key] = merge_pairs([pairs_store.get(key, {}), pairs])
            results.append({"receiver": receiver, "model": model, "feature_set": fset, "config": config, "fold": fold,
                            "target": "genes_pooled", "target_kind": "gene_vector", "n_test": int(test.sum()),
                            "inner_val_mse": inner_val_mse, "test_r2_fold": pooled_r2(pairs), **(extra or {})})

    for k in range(n_folds):
        tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
        ytr_mean = data.Y[tr_all].mean(axis=0)
        yin_mean = data.Y[tr_in].mean(axis=0)
        t0 = time.time()
        # (1) mean predictor
        record("mean", "none", "train_mean", k, np.tile(ytr_mean, (test.sum(), 1)), test, ytr_mean,
               float(((data.Y[val] - yin_mean) ** 2).mean()))
        for fset, blocks in FEATURE_SETS.items():
            Xall = data.matrix(blocks)
            F.assert_no_target_leak(data.feature_set(blocks), None, data.Y, set(spec["genes"]), set(spec["programs"]))
            mu, sd = F.fit_scaler(Xall[tr_in])             # inner-train statistics for selection
            Xin, Xval = F.apply_scaler(Xall[tr_in], mu, sd), F.apply_scaler(Xall[val], mu, sd)
            mu2, sd2 = F.fit_scaler(Xall[tr_all])           # refit statistics on all training folds
            Xtr2, Xte2 = F.apply_scaler(Xall[tr_all], mu2, sd2), F.apply_scaler(Xall[test], mu2, sd2)
            # (2)/(3) ridge: grid on inner val, refit on train_all
            val_scores = {}
            for alpha in cfg["ridge_alphas"]:
                m = fit_ridge(Xin, data.Y[tr_in], alpha)
                val_scores[alpha] = float(((m.predict(Xval) - data.Y[val]) ** 2).mean())
            best_a = min(val_scores, key=val_scores.get)
            for alpha, vmse in val_scores.items():
                m = fit_ridge(Xtr2, data.Y[tr_all], alpha)
                record("ridge", fset, f"alpha={alpha}", k, m.predict(Xte2), test, ytr_mean, vmse,
                       {"selected": alpha == best_a, "n_features": Xall.shape[1]})
            # (3) HGB on program targets only, fsets programs/genes
            if fset in ("programs", "genes") and cfg.get("run_hgb", True):
                for lr in cfg["hgb"]["learning_rate"]:
                    for leaves in cfg["hgb"]["max_leaf_nodes"]:
                        preds, vm, iters = np.zeros((test.sum(), len(tnames)), np.float32), [], []
                        for j in range(n_prog):
                            m, it, vmse = fit_hgb_staged(Xin, data.Y[tr_in, j], Xval, data.Y[val, j], lr, leaves,
                                                         cfg["hgb"]["max_iter"], C.SEED, cfg["hgb"]["max_rows"])
                            from sklearn.ensemble import HistGradientBoostingRegressor
                            rng = np.random.default_rng(C.SEED)
                            rows = np.where(tr_all)[0]
                            if len(rows) > cfg["hgb"]["max_rows"]:
                                rows = rng.choice(rows, cfg["hgb"]["max_rows"], replace=False)
                            m2 = HistGradientBoostingRegressor(learning_rate=lr, max_leaf_nodes=leaves, max_iter=it,
                                                               min_samples_leaf=50, early_stopping=False,
                                                               random_state=C.SEED).fit(F.apply_scaler(Xall[rows], mu2, sd2), data.Y[rows, j])
                            preds[:, j] = m2.predict(Xte2)
                            vm.append(vmse); iters.append(it)
                        preds[:, n_prog:] = ytr_mean[n_prog:]  # gene columns not modelled by HGB (scope)
                        record("hgb", fset, f"lr={lr},leaves={leaves}", k, preds, test, ytr_mean, float(np.mean(vm)),
                               {"best_iters": ";".join(map(str, iters)), "n_features": Xall.shape[1]})
            # (4) MLP on genes
            if fset == "genes" and cfg.get("run_mlp", True):
                for lr in cfg["mlp"]["learning_rate"]:
                    for seed in cfg["mlp"]["seeds"]:
                        mlp = TorchMLP(Xin.shape[1], len(tnames), tuple(cfg["mlp"]["hidden"]), cfg["mlp"]["dropout"], lr,
                                       cfg["mlp"]["weight_decay"], seed, cfg["mlp"]["batch"], cfg["mlp"]["max_epochs"],
                                       cfg["mlp"]["patience"]).fit(Xin, data.Y[tr_in], Xval, data.Y[val])
                        record("mlp", fset, f"lr={lr},seed={seed}", k, mlp.predict(F.apply_scaler(Xall[test], mu, sd)),
                               test, ytr_mean, mlp.best_val, {"epochs": mlp.epochs, "n_features": Xall.shape[1]})
                        pd.DataFrame(mlp.history).assign(receiver=receiver, fold=k, lr=lr, seed=seed).to_csv(
                            C.RESULTS / PHASE / f"mlp_curve_{receiver}_fold{k}_lr{lr}_seed{seed}.csv", index=False)
        log.info("%s fold %d done in %.0f s", receiver, k, time.time() - t0)
    return results, pairs_store


def run_receiver_gpu(receiver: str, data: A.ReceiverData, fold_table: pd.DataFrame, cfg: dict, log) -> tuple[list, dict]:
    """GPU version of `run_receiver` (Phase 2 plan, Amendment 1): same folds, grid, metrics and records.

    Feature blocks are moved to the GPU once (host copies freed); ridge uses `gpu.RidgePath` (all alphas from one
    eigendecomposition), gradient boosting uses XGBoost-CUDA ('xgb'), the MLP trains on GPU-resident tensors.
    """
    import torch
    import gpu as GP
    n_folds = int(fold_table.fold.max()) + 1
    tnames = data.target_names
    spec = A.targets_spec()[receiver]
    n_prog = len(spec["programs"])
    results, pairs_store = [], {}
    for fset, blocks in FEATURE_SETS.items():
        F.assert_no_target_leak(data.feature_set(blocks), None, data.Y, set(spec["genes"]), set(spec["programs"]))
    Xblk = {}
    for b in list(data.X):                                           # each block once on the GPU; host copy freed
        Xblk[b] = GP.to_gpu(data.X[b])
        data.X[b] = None
    Yg = GP.to_gpu(data.Y)

    def record(model, fset, config, fold, pred_test, test, ytr_mean, inner_val_mse, extra=None):
        yte = data.Y[test]
        cols = [(t, slice(j, j + 1), "program" if j < n_prog else "gene") for j, t in enumerate(tnames)]
        if len(tnames) > n_prog:
            cols.append(("genes_pooled", slice(n_prog, None), "gene_vector"))
        for t, sl, kind in cols:
            pairs = sse_sst_by_group(yte[:, sl], pred_test[:, sl], data.group[test], ytr_mean[sl])
            key = (model, fset, config, t)
            pairs_store[key] = merge_pairs([pairs_store.get(key, {}), pairs])
            results.append({"receiver": receiver, "model": model, "feature_set": fset, "config": config, "fold": fold,
                            "target": t, "target_kind": kind, "n_test": int(test.sum()), "inner_val_mse": inner_val_mse,
                            "test_r2_fold": pooled_r2(pairs), **(extra or {})})

    for k in range(n_folds):
        tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
        it, iv, ia, ie = (torch.as_tensor(np.where(m)[0], device="cuda") for m in (tr_in, val, tr_all, test))
        ytr_mean = data.Y[tr_all].mean(axis=0)
        yin_mean = data.Y[tr_in].mean(axis=0)
        t0 = time.time()
        record("mean", "none", "train_mean", k, np.tile(ytr_mean, (test.sum(), 1)), test, ytr_mean,
               float(((data.Y[val] - yin_mean) ** 2).mean()))
        for fset, blocks in FEATURE_SETS.items():
            Xall = torch.cat([Xblk[b] for b in blocks], 1)
            mu, sd = GP.fit_scaler(Xall[it])
            Xin, Xval = GP.apply_scaler(Xall[it], mu, sd), GP.apply_scaler(Xall[iv], mu, sd)
            mu2, sd2 = GP.fit_scaler(Xall[ia])
            Xtr2, Xte2 = GP.apply_scaler(Xall[ia], mu2, sd2), GP.apply_scaler(Xall[ie], mu2, sd2)
            del Xall
            # ridge: alpha chosen on inner val; every alpha refitted on all training folds and recorded
            rp = GP.RidgePath(Xin, Yg[it])
            val_scores = {a: float(((rp.predict(Xval, a) - Yg[iv]) ** 2).mean()) for a in cfg["ridge_alphas"]}
            best_a = min(val_scores, key=val_scores.get)
            rp2 = GP.RidgePath(Xtr2, Yg[ia])
            for alpha, vmse in val_scores.items():
                record("ridge", fset, f"alpha={alpha}", k, rp2.predict(Xte2, alpha).cpu().numpy(), test, ytr_mean, vmse,
                       {"selected": alpha == best_a, "n_features": Xin.shape[1]})
            del rp, rp2
            torch.cuda.empty_cache()   # return torch's cached blocks so XGBoost (own allocator) has GPU memory
            # gradient boosting (XGBoost CUDA) on program targets; <= max_rows training rows
            if fset in ("programs", "genes") and cfg.get("run_hgb", True):
                rng = np.random.default_rng(C.SEED)
                sub_in = rng.choice(len(it), min(len(it), cfg["hgb"]["max_rows"]), replace=False)
                sub_all = np.random.default_rng(C.SEED).choice(len(ia), min(len(ia), cfg["hgb"]["max_rows"]), replace=False)
                Xin_np, Xval_np = Xin[torch.as_tensor(sub_in, device="cuda")].cpu().numpy(), Xval.cpu().numpy()
                Xtr_np, Xte_np = Xtr2[torch.as_tensor(sub_all, device="cuda")].cpu().numpy(), Xte2.cpu().numpy()
                yin_np, yval_np = data.Y[tr_in][sub_in], data.Y[val]
                ytr_np = data.Y[tr_all][sub_all]
                for lr in cfg["hgb"]["learning_rate"]:
                    for leaves in cfg["hgb"]["max_leaf_nodes"]:
                        preds = np.tile(ytr_mean, (test.sum(), 1)).astype(np.float32)   # gene columns: not modelled (scope)
                        vm, iters = [], []
                        for j in range(n_prog):
                            best_it, vmse = GP.xgb_fit_select(Xin_np, yin_np[:, j], Xval_np, yval_np[:, j], lr, leaves,
                                                              cfg["hgb"]["max_iter"], C.SEED)
                            preds[:, j] = GP.xgb_fit_predict(Xtr_np, ytr_np[:, j], Xte_np, lr, leaves, best_it, C.SEED)
                            vm.append(vmse); iters.append(best_it)
                        record("xgb", fset, f"lr={lr},leaves={leaves}", k, preds, test, ytr_mean, float(np.mean(vm)),
                               {"best_iters": ";".join(map(str, iters)), "n_features": Xin.shape[1]})
                del Xin_np, Xval_np, Xtr_np, Xte_np
                import gc; gc.collect()
            # MLP (GPU) on the genes set; early stopping on inner val; no refit
            if fset == "genes" and cfg.get("run_mlp", True):
                Xte_inner = GP.apply_scaler(torch.cat([Xblk[b] for b in blocks], 1)[ie], mu, sd)
                for lr in cfg["mlp"]["learning_rate"]:
                    for seed in cfg["mlp"]["seeds"]:
                        mlp = GP.GPUMLP(Xin.shape[1], len(tnames), tuple(cfg["mlp"]["hidden"]), cfg["mlp"]["dropout"], lr,
                                        cfg["mlp"]["weight_decay"], seed, cfg["mlp"]["batch"], cfg["mlp"]["max_epochs"],
                                        cfg["mlp"]["patience"]).fit(Xin, Yg[it], Xval, Yg[iv])
                        record("mlp", fset, f"lr={lr},seed={seed}", k, mlp.predict(Xte_inner).cpu().numpy(), test, ytr_mean,
                               mlp.best_val, {"epochs": mlp.epochs, "n_features": Xin.shape[1]})
                        pd.DataFrame(mlp.history).assign(receiver=receiver, fold=k, lr=lr, seed=seed).to_csv(
                            C.RESULTS / PHASE / f"mlp_curve_{receiver}_fold{k}_lr{lr}_seed{seed}.csv", index=False)
                        del mlp
                del Xte_inner
            del Xin, Xval, Xtr2, Xte2
            torch.cuda.empty_cache()
        log.info("%s fold %d done in %.0f s (peak VRAM %.0f MiB)", receiver, k, time.time() - t0,
                 torch.cuda.max_memory_allocated() / 2**20)
    del Xblk, Yg
    torch.cuda.empty_cache()
    return results, pairs_store


def summarise(results: pd.DataFrame, pairs_store: dict, receiver: str, cfg: dict) -> pd.DataFrame:
    """Select per model/feature set the configuration with the lowest mean inner-validation MSE (averaged over
    folds); for the MLP, report the mean over seeds of the selected learning rate. Then pooled R2 + CI and gain."""
    rows = []
    r = results[results.receiver == receiver]
    sel = {}
    for (model, fset), d in r.groupby(["model", "feature_set"]):
        if model == "mlp":
            # average inner val over seeds per lr; choose lr; report each seed and the seed-mean
            d2 = d.assign(lr=d.config.str.extract(r"lr=([^,]+)")[0])
            best_lr = d2.groupby("lr").inner_val_mse.mean().idxmin()
            sel[(model, fset)] = [c for c in d2[d2.lr == best_lr].config.unique()]
        else:
            best = d.groupby("config").inner_val_mse.mean().idxmin()
            sel[(model, fset)] = [best]
    own_pairs = {t: pairs_store[("ridge", "own", sel[("ridge", "own")][0], t)] for t in r.target.unique()}
    for (model, fset), configs in sel.items():
        for t in sorted(r.target.unique()):
            per_cfg = [pairs_store[(model, fset, c, t)] for c in configs if (model, fset, c, t) in pairs_store]
            if not per_cfg:
                continue
            # seed-mean for MLP: average SSE across seeds per tumour (same SST)
            pairs = {g: (float(np.mean([p[g][0] for p in per_cfg])), per_cfg[0][g][1]) for g in per_cfg[0]}
            r2, lo, hi = C.tumour_bootstrap(pairs, cfg["n_boot"])
            gain, glo, ghi = C.tumour_bootstrap_diff(own_pairs[t], pairs, cfg["n_boot"])
            rows.append({"receiver": receiver, "model": model, "feature_set": fset, "selected_config": ";".join(configs),
                         "target": t, "r2": r2, "ci_low": lo, "ci_high": hi, "gain_vs_own": gain, "gain_ci_low": glo,
                         "gain_ci_high": ghi, "gain_ci_excludes_0": bool(glo > 0 or ghi < 0),
                         "n_tumours": len(pairs), "n_seeds": len(configs)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- reproduction check
def reproduction_check(cfg: dict, log) -> pd.DataFrame:
    """Chain Phase 2 ridge, re-run from the Phase 0 datasets with the chain definitions (see plan)."""
    from sklearn.linear_model import RidgeCV
    chain = pd.read_csv(C.REPO_ROOT / "results/chain/phase2/delta_r2_by_program.csv")
    chain = chain[(chain.annotation == "author_mapped") & (~chain.spillover_covariates.astype(bool)) & (chain.model == "ridge")]
    fold_table = A.folds("loto")
    spec = A.targets_spec()
    classes = A.all_classes("development")
    rows = []
    type_index = {"tumour": 0, "CAF": 1, "T": 2, "myeloid": 3}
    for receiver in ("myeloid", "tumour"):
        s = spec[receiver]
        # build chain-style features: cap 20k per section with chain seeds, demean within section
        parts = []
        for s_i, path in enumerate(A.section_files("development")):
            d = A.load_section(path)
            fs, y = F.build_section_features(d, s["class"], s["programs"], [], classes, include_gene_block=False)
            if len(fs.rows) == 0:
                continue
            rng = np.random.default_rng([C.SEED, s_i, type_index[receiver]])
            idx = np.arange(len(fs.rows))
            if len(idx) > 20000:
                idx = np.sort(rng.choice(idx, 20000, replace=False))
            # chain Phase 2 used, per typed class, only that class's OWN programs as neighbour features
            nb_cols = [j for j, n in enumerate(fs.names["nb_programs"]) if n.split("__")[0].split("_")[-1] == n.split("__")[1].split(":")[0]]
            for p_i, prog in enumerate(s["programs"]):
                own = fs.blocks["own"][idx][:, [0, 2]]                                # log1p total counts, density
                tile = fs.blocks["tile"][idx][:, [p_i, -1]]                            # tile LOO mean of this program, empty flag
                tile = np.where(np.isnan(tile), y[:, p_i].mean(), tile)                # chain: section mean of the class
                base = np.column_stack([own, tile])
                nbr = np.column_stack([fs.blocks["composition"][idx], fs.blocks["nb_programs"][idx][:, nb_cols]])
                yy = y[idx, p_i:p_i + 1]
                M = np.column_stack([yy, base, nbr]).astype(np.float64)
                M -= M.mean(axis=0, keepdims=True)
                parts.append((prog, str(d["group"]), M, base.shape[1]))
        for prog in s["programs"]:
            pp = [p for p in parts if p[0] == prog]
            M = np.vstack([p[2] for p in pp]); tum = np.concatenate([[p[1]] * len(p[2]) for p in pp]); nb = pp[0][3]
            y, Xb, Xn = M[:, 0], M[:, 1:1 + nb], M[:, 1:]
            sse_b = sse_n = sst = 0.0
            for held in sorted(set(tum)):
                te, tr = tum == held, tum != held
                for X, acc in ((Xb, "b"), (Xn, "n")):
                    mu, sd = X[tr].mean(0), X[tr].std(0); sd[sd == 0] = 1
                    m = RidgeCV(alphas=np.logspace(-2, 4, 13)).fit((X[tr] - mu) / sd, y[tr])
                    e = float(((y[te] - m.predict((X[te] - mu) / sd)) ** 2).sum())
                    if acc == "b":
                        sse_b += e
                    else:
                        sse_n += e
                sst += float((y[te] ** 2).sum())
            delta = (sse_b - sse_n) / sst
            ref = float(chain[chain.program == prog].delta_r2.iloc[0])
            rows.append({"receiver": receiver, "program": prog, "delta_r2_reproduced": delta, "delta_r2_chain": ref,
                         "abs_diff": abs(delta - ref), "tolerance": cfg["reproduction_tolerance"],
                         "pass": abs(delta - ref) <= cfg["reproduction_tolerance"], "n_cells": len(y)})
            log.info("reproduction %s: delta R2 %.4f vs chain %.4f (|diff| %.4f)", prog, delta, ref, abs(delta - ref))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "baselines.yaml"))
    ap.add_argument("--skip-reproduction", action="store_true")
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "baselines", cfg)
    if not args.skip_reproduction:
        rep = reproduction_check(cfg, log)
        C.save_table(rep, PHASE, "reproduction_check.csv", log)
        if not rep["pass"].all():
            log.error("reproduction check FAILED; stopping per plan:\n%s", rep.to_string(index=False))
            sys.exit(2)
    fold_table = A.folds("grouped")
    all_results, summaries = [], []
    for receiver in cfg["receivers"]:
        res_path = C.RESULTS / PHASE / f"baseline_results_{receiver}.csv"
        pairs_path = C.RESULTS / PHASE / f"per_tumour_sse_{receiver}.csv"
        if res_path.exists() and pairs_path.exists():      # resume unit: one receiver
            df = pd.read_csv(res_path)
            pt = pd.read_csv(pairs_path)
            pairs = {}
            for r in pt.itertuples(index=False):
                pairs.setdefault((r.model, r.feature_set, r.config, r.target), {})[str(r.tumour)] = (r.sse, r.sst)
            log.info("%s: resumed from %s", receiver, C.rel(res_path))
        else:
            data = A.build_receiver("development", receiver, log)
            log.info("%s: %d cells, targets %s", receiver, len(data.Y), data.target_names)
            res, pairs = run_receiver_gpu(receiver, data, fold_table, cfg, log)
            del data
            df = pd.DataFrame(res)
            df.to_csv(res_path, index=False)
            pd.DataFrame([{"receiver": receiver, "model": k[0], "feature_set": k[1], "config": k[2], "target": k[3],
                           "tumour": g, "sse": v[0], "sst": v[1]} for k, p in pairs.items() for g, v in p.items()]
                         ).to_csv(pairs_path, index=False)
        all_results.append(df)
        summaries.append(summarise(df, pairs, receiver, cfg))
    res = pd.concat(all_results, ignore_index=True)
    summ = pd.concat(summaries, ignore_index=True)
    C.save_table(res, PHASE, "baseline_results.csv", log)
    C.save_table(summ, PHASE, "baseline_summary.csv", log)
    prog = summ[summ.target.isin([t for r in cfg["receivers"] for t in A.targets_spec()[r]["programs"]] + ["genes_pooled"])]
    log.info("summary (selected configs):\n%s", prog.round(4).to_string(index=False))
    best = (summ[summ.model != "mean"].sort_values("r2", ascending=False).groupby(["receiver", "target"]).head(1))
    C.save_table(best, PHASE, "strongest_baseline_per_target.csv", log)
    figure(summ, log)
    log.info("done")


def figure(summ: pd.DataFrame, log) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    targets = [t for t in summ.target.unique() if ":" in t or t == "genes_pooled"]
    fig, axes = plt.subplots(1, len(targets), figsize=(3.2 * len(targets), 4.5), sharey=False)
    for ax, t in zip(np.atleast_1d(axes), targets):
        d = summ[summ.target == t].copy()
        d["label"] = d.model + "/" + d.feature_set
        d = d.sort_values("r2")
        ax.errorbar(d.r2, range(len(d)), xerr=[d.r2 - d.ci_low, d.ci_high - d.r2], fmt="o", capsize=2)
        ax.set_yticks(range(len(d)), d.label, fontsize=7)
        ax.axvline(0, color="k", lw=.6)
        ax.set_title(t, fontsize=9)
        ax.set_xlabel("pooled held-out R2 (95% tumour CI)")
    fig.tight_layout()
    d = C.RESULTS / PHASE
    fig.savefig(d / "baselines.png", dpi=200); fig.savefig(d / "baselines.pdf")
    plt.close(fig)
    log.info("saved figure %s", C.rel(d / "baselines.png"))


if __name__ == "__main__":
    main()

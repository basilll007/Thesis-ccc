"""Phase 3: train and evaluate the padded-tree graph models (and their ablations) on the Phase 1 folds.

What: for each receiver, outer fold, architecture, learning rate, seed and ablation mode, train on the inner
training folds with early stopping on the inner validation fold, evaluate on the outer test fold, and record
per-tumour (SSE, SST) pairs for every target column. `ablations.py` reuses `run_experiments` with the modes of the
Phase 3 plan. Summaries compare the seed-mean graph model with the strongest Phase 2 baseline using a paired
tumour bootstrap.
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
import graphs as G  # noqa: E402
import models as M  # noqa: E402
from baselines import merge_pairs, pooled_r2, sse_sst_by_group  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

PHASE = "phase3"


# ----------------------------------------------------------------------------- store construction
def build_store(split: str, kmax: int, device: str, log, shuffle_seed: int | None = None) -> tuple[M.GraphStore, dict, list[str]]:
    """Concatenate all sections of a split into one GraphStore (raw log expression; standardised later per fold)."""
    classes = A.all_classes(split)
    cls_id = {c: i for i, c in enumerate(classes)}
    Xs, covs, cls, nbrs, dists, secs, offsets, tiles = [], [], [], [], [], [], {}, []
    n0 = 0
    trunc_num = trunc_den = 0
    for path in A.section_files(split):
        d = A.load_section(path)
        n = len(d["label"])
        offsets[str(d["section"])] = n0
        nbr, dist = G.padded_neighbours(d["knn_indptr"], d["knn_indices"], d["knn_length"], kmax)
        deg = np.diff(d["knn_indptr"])
        trunc_num += int((deg > kmax).sum()); trunc_den += n
        nbr = np.where(nbr >= 0, nbr + n0, -1)
        scale = float(d["scale_units_per_um"])
        Xs.append(d["X_log"]); covs.append(np.column_stack([np.log1p(d["total_counts"]), np.log1p(d["n_genes_detected"])]))
        cls.append(np.array([cls_id[x] for x in d["label"].astype(str)])); nbrs.append(nbr)
        dists.append(np.log1p(dist / scale).astype(np.float32)); secs.append(np.array([str(d["section"])] * n))
        tiles.append(d["tile"] + 10_000_000 * len(tiles))
        n0 += n
    X = torch.as_tensor(np.vstack(Xs)).to(device)                       # float16
    cov = torch.as_tensor(np.vstack(covs).astype(np.float32)).to(device)
    cls_t = torch.as_tensor(np.concatenate(cls)).long().to(device)
    nbr_t = torch.as_tensor(np.vstack(nbrs)).long().to(device)
    dist_t = torch.as_tensor(np.vstack(dists)).to(device)
    store = M.GraphStore(X, cov, cls_t, nbr_t, dist_t, np.concatenate(secs))
    if shuffle_seed is not None:
        tile_all, cls_all = np.concatenate(tiles), np.concatenate(cls)
        rng = np.random.default_rng([C.SEED, 7, shuffle_seed])
        perm = np.arange(n0)
        key = tile_all.astype(np.int64) * 64 + cls_all
        order = np.argsort(key, kind="stable")
        sorted_key = key[order]
        starts = np.r_[0, np.where(np.diff(sorted_key) != 0)[0] + 1, len(order)]
        for a, b in zip(starts[:-1], starts[1:]):
            grp = order[a:b]
            perm[grp] = grp[rng.permutation(len(grp))]
        store.tile_perm = torch.as_tensor(perm).long().to(device)
    log.info("store %s: %d cells, %d classes, truncated-degree fraction %.4f", split, n0, len(classes), trunc_num / trunc_den)
    return store, offsets, classes


def fold_standardisation(store: M.GraphStore, train_mask_global: torch.Tensor):
    """Per-gene mean/sd of log expression and of the counts covariates over TRAINING-fold cells only."""
    Xt = store.X[train_mask_global].float()
    mu, sd = Xt.mean(0), Xt.std(0)
    sd[sd == 0] = 1
    cmu, csd = store.cov[train_mask_global].mean(0), store.cov[train_mask_global].std(0)
    csd[csd == 0] = 1
    return mu, sd, cmu, csd


class StandardisedStore(M.GraphStore):
    """View of a GraphStore whose X/cov are standardised on the fly with training-fold statistics."""

    def __init__(self, base: M.GraphStore, mu, sd, cmu, csd):
        self.__dict__.update(base.__dict__)
        self._mu, self._sd, self._cmu, self._csd = mu, sd, cmu, csd
        self._X_raw, self._cov_raw = base.X, base.cov

    class _Std:
        def __init__(self, raw, mu, sd):
            self.raw, self.mu, self.sd = raw, mu, sd

        def __getitem__(self, idx):
            return (self.raw[idx].float() - self.mu) / self.sd

        @property
        def shape(self):
            return self.raw.shape

    @property
    def X(self):
        return self._Std(self._X_raw, self._mu, self._sd)

    @X.setter
    def X(self, v):
        self._X_raw = v

    @property
    def cov(self):
        return self._Std(self._cov_raw, self._cmu, self._csd)

    @cov.setter
    def cov(self, v):
        self._cov_raw = v


# ----------------------------------------------------------------------------- training
def train_one(store, seed_global: np.ndarray, seed_cov: np.ndarray, Y: np.ndarray, masks, receiver_cls: int, mode: str,
              arch: str, lr: float, seed: int, cfg: dict, n_classes: int, device: str, log):
    """Train one model; returns (predictions on test, history, epochs, best val mse, peak VRAM, model state)."""
    C.set_seed(seed)
    tr, val, te = masks
    k1, k2 = cfg["k1"], cfg["k2"]
    n_out = Y.shape[1]
    model = M.TreeModel(store.X.shape[1], n_classes, seed_cov.shape[1], n_out, cfg["hidden"], cfg["emb"],
                        "attention" if arch == "attention" else "mean", cfg["dropout"], cfg["heads"]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=cfg["weight_decay"])
    cov_t = torch.as_tensor(seed_cov).to(device)
    Y_t = torch.as_tensor(Y).to(device)
    g_t = torch.as_tensor(seed_global).long().to(device)
    tr_idx, val_idx, te_idx = (np.where(m)[0] for m in (tr, val, te))
    gen = torch.Generator().manual_seed(seed)
    use_amp = device == "cuda"
    batch = cfg["batch"]

    def predict(idx: np.ndarray) -> np.ndarray:
        model.eval()
        out = []
        with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
            for i in range(0, len(idx), batch * 2):
                b = torch.as_tensor(idx[i:i + batch * 2]).to(device)
                tb = M.TreeBatch(store, g_t[b], cov_t[b], receiver_cls, mode, k1, k2)
                out.append(model(tb).float().cpu().numpy())
        return np.vstack(out)

    history, best, bad, best_state = [], np.inf, 0, None
    for epoch in range(cfg["max_epochs"]):
        model.train()
        perm = tr_idx[torch.randperm(len(tr_idx), generator=gen).numpy()]
        tot = 0.0
        for i in range(0, len(perm), batch):
            b = torch.as_tensor(perm[i:i + batch]).to(device)
            tb = M.TreeBatch(store, g_t[b], cov_t[b], receiver_cls, mode, k1, k2)
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                pred = model(tb)
            loss = ((pred.float() - Y_t[b]) ** 2).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            tot += float(loss) * len(b)
        vl = float(((predict(val_idx) - Y[val_idx]) ** 2).mean())
        history.append({"epoch": epoch, "train_mse": tot / len(tr_idx), "val_mse": vl})
        if vl < best - 1e-6:
            best, bad = vl, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= cfg["patience"]:
                break
    model.load_state_dict(best_state)
    pred_te = predict(te_idx)
    peak = torch.cuda.max_memory_allocated() / 2**20 if device == "cuda" else 0.0
    return pred_te, history, len(history), best, peak, best_state


class TimeBudgetReached(RuntimeError):
    """A foreground chunk used its time budget; completed runs are saved per run (resume on the next call)."""


def run_experiments(cfg: dict, modes: list[str], archs: list[str], lrs: list[float], phase: str, tag: str, log,
                    receivers_modes: dict | None = None) -> tuple[pd.DataFrame, dict]:
    t_start = time.time()
    if not torch.cuda.is_available():
        raise RuntimeError("GPU-only policy: CUDA not available, refusing to train on CPU")
    device = "cuda"
    fold_table = A.folds("grouped")
    n_folds = int(fold_table.fold.max()) + 1
    spec = A.targets_spec()
    results, pairs_store, compute = [], {}, []
    need_shuffle = "shuffled" in modes
    for receiver in cfg["receivers"]:
        s = spec[receiver]
        done = [C.RESULTS / phase / f"{n}_{tag}_{receiver}.csv" for n in ("results", "compute", "per_tumour_sse")]
        if all(p.exists() for p in done):   # resume unit: one receiver
            results += pd.read_csv(done[0]).to_dict("records")
            compute += pd.read_csv(done[1]).to_dict("records")
            pairs_store.update(pairs_from_csv(phase, tag, [receiver]))
            log.info("%s: resumed from %s", receiver, ", ".join(C.rel(p) for p in done))
            continue
        n_res0, n_comp0 = len(results), len(compute)
        if not A.section_files("development"):
            raise RuntimeError("no datasets")
        if not (C.CACHE / "features" / f"development__{receiver}.npz").exists():
            A.build_receiver("development", receiver, log)      # build the full cache once if missing
        data = A.build_receiver("development", receiver, log, blocks=["own", "tile"])
        F.assert_no_target_leak(data.feature_set(["own", "tile"]), None, data.Y, set(s["genes"]), set(s["programs"]))
        seed_cov_raw = data.matrix(["own", "tile"])
        tnames, n_prog = data.target_names, len(s["programs"])
        r_modes = (receivers_modes or {}).get(receiver, modes)
        for seed in cfg["seeds"]:
            store0, offsets, classes = build_store("development", max(cfg["k1"], cfg["k2"]), device, log,
                                                   shuffle_seed=seed if need_shuffle else None)
            receiver_cls = classes.index(s["class"])
            seed_global = np.array([offsets[sec] + r for sec, r in zip(data.section, data.row)])
            for k in range(n_folds):
                tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
                ytr_mean = data.Y[tr_all].mean(axis=0)
                # standardisation: training-fold cells of all types (section membership, no labels)
                train_secs = set(fold_table[fold_table.fold != k].section)
                gmask = torch.as_tensor(np.isin(store0.section_of, list(train_secs))).to(device)
                mu, sd, cmu, csd = fold_standardisation(store0, gmask)
                store = StandardisedStore(store0, mu, sd, cmu, csd)
                smu, ssd = F.fit_scaler(seed_cov_raw[tr_in])
                seed_cov = F.apply_scaler(seed_cov_raw, smu, ssd)
                for mode in r_modes:
                    for arch in archs:
                        for lr in lrs:
                            run_file = C.RESULTS / phase / "partial" / f"run_{tag}_{receiver}_{mode}_{arch}_lr{lr}_seed{seed}_fold{k}.json"
                            if run_file.exists():          # resume unit: one training run (checkpoint + result)
                                rec = json.loads(run_file.read_text())
                                results += rec["results"]; compute.append(rec["compute"])
                                for t, pr in rec["pairs"].items():
                                    key = (receiver, mode, arch, lr, seed, t)
                                    pairs_store[key] = merge_pairs([pairs_store.get(key, {}), {g: tuple(v) for g, v in pr.items()}])
                                continue
                            if cfg.get("_time_budget") and time.time() - t_start > cfg["_time_budget"]:
                                raise TimeBudgetReached(f"time budget {cfg['_time_budget']} s used; rerun to continue")
                            n_r0 = len(results)
                            run_pairs: dict = {}
                            t0 = time.time()
                            torch.cuda.reset_peak_memory_stats() if device == "cuda" else None
                            pred, hist, ep, bv, peak, state = train_one(store, seed_global, seed_cov, data.Y, (tr_in, val, test),
                                                                        receiver_cls, mode, arch, lr, seed, cfg, len(classes), device, log)
                            yte = data.Y[test]
                            cfg_name = f"{arch}|lr={lr}|seed={seed}|mode={mode}"
                            for j, t in enumerate(tnames):
                                pairs = sse_sst_by_group(yte[:, j:j + 1], pred[:, j:j + 1], data.group[test], ytr_mean[j:j + 1])
                                key = (receiver, mode, arch, lr, seed, t)
                                pairs_store[key] = merge_pairs([pairs_store.get(key, {}), pairs])
                                run_pairs[t] = pairs
                                results.append({"receiver": receiver, "mode": mode, "arch": arch, "lr": lr, "seed": seed, "fold": k,
                                                "target": t, "target_kind": "program" if j < n_prog else "gene",
                                                "inner_val_mse": bv, "epochs": ep, "test_r2_fold": pooled_r2(pairs)})
                            if len(tnames) > n_prog:
                                pairs = sse_sst_by_group(yte[:, n_prog:], pred[:, n_prog:], data.group[test], ytr_mean[n_prog:])
                                key = (receiver, mode, arch, lr, seed, "genes_pooled")
                                pairs_store[key] = merge_pairs([pairs_store.get(key, {}), pairs])
                                run_pairs["genes_pooled"] = pairs
                                results.append({"receiver": receiver, "mode": mode, "arch": arch, "lr": lr, "seed": seed, "fold": k,
                                                "target": "genes_pooled", "target_kind": "gene_vector", "inner_val_mse": bv,
                                                "epochs": ep, "test_r2_fold": pooled_r2(pairs)})
                            cdir = C.RESULTS / phase / "curves"
                            cdir.mkdir(parents=True, exist_ok=True)
                            pd.DataFrame(hist).to_csv(cdir / f"{tag}_{receiver}_{mode}_{arch}_lr{lr}_seed{seed}_fold{k}.csv", index=False)
                            ck = C.CACHE / "checkpoints" / phase / f"{tag}_{receiver}_{mode}_{arch}_lr{lr}_seed{seed}_fold{k}.pt"
                            ck.parent.mkdir(parents=True, exist_ok=True)
                            torch.save({"state": state, "mu": mu.cpu(), "sd": sd.cpu(), "cmu": cmu.cpu(), "csd": csd.cpu(),
                                        "smu": smu, "ssd": ssd, "classes": classes, "cfg": cfg, "arch": arch, "lr": lr,
                                        "seed": seed, "mode": mode, "epochs": ep}, ck)
                            compute.append({"receiver": receiver, "mode": mode, "arch": arch, "lr": lr, "seed": seed, "fold": k,
                                            "wall_s": time.time() - t0, "epochs": ep, "peak_vram_mib": peak, "device": device,
                                            "n_train": int(tr_in.sum()), "checkpoint": C.rel(ck)})
                            run_file.parent.mkdir(parents=True, exist_ok=True)
                            run_file.write_text(json.dumps({"results": results[n_r0:], "compute": compute[-1], "pairs": run_pairs}, default=float))
                            log.info("%s %s %s lr=%g seed=%d fold %d: %d epochs, val mse %.4f, %.0fs, peak %.0f MiB",
                                     receiver, mode, arch, lr, seed, k, ep, bv, time.time() - t0, peak)
            del store0
            torch.cuda.empty_cache() if device == "cuda" else None
        pd.DataFrame([{"receiver": k[0], "mode": k[1], "arch": k[2], "lr": k[3], "seed": k[4], "target": k[5], "tumour": g,
                       "sse": v[0], "sst": v[1]} for k, p in pairs_store.items() if k[0] == receiver for g, v in p.items()]
                     ).to_csv(C.RESULTS / phase / f"per_tumour_sse_{tag}_{receiver}.csv", index=False)
        pd.DataFrame(results[n_res0:]).to_csv(C.RESULTS / phase / f"results_{tag}_{receiver}.csv", index=False)
        pd.DataFrame(compute[n_comp0:]).to_csv(C.RESULTS / phase / f"compute_{tag}_{receiver}.csv", index=False)
        del data
    pd.DataFrame(compute).to_csv(C.RESULTS / phase / f"compute_{tag}.csv", index=False)
    return pd.DataFrame(results), pairs_store


def baseline_pairs(receiver: str, target: str) -> tuple[str, dict]:
    """Per-tumour pairs of the strongest Phase 2 baseline for this receiver/target."""
    best = pd.read_csv(C.RESULTS / "phase2" / "strongest_baseline_per_target.csv")
    b = best[(best.receiver == receiver) & (best.target == target)].iloc[0]
    per = pd.read_csv(C.RESULTS / "phase2" / f"per_tumour_sse_{receiver}.csv")
    cfgs = b.selected_config.split(";")
    d = per[(per.model == b.model) & (per.feature_set == b.feature_set) & per.config.isin(cfgs) & (per.target == target)]
    pairs = {g: (float(x.sse.mean()), float(x.sst.iloc[0])) for g, x in d.groupby("tumour")}
    return f"{b.model}/{b.feature_set}", pairs


def seed_mean_pairs(pairs_store: dict, receiver: str, mode: str, arch: str, lr: float, seeds: list[int], target: str) -> dict:
    per = [pairs_store[(receiver, mode, arch, lr, sd, target)] for sd in seeds if (receiver, mode, arch, lr, sd, target) in pairs_store]
    return {g: (float(np.mean([p[g][0] for p in per])), per[0][g][1]) for g in per[0]} if per else {}


def summarise_graph(res: pd.DataFrame, pairs_store: dict, cfg: dict, log) -> pd.DataFrame:
    rows = []
    for receiver in cfg["receivers"]:
        r = res[(res.receiver == receiver) & (res["mode"] == "full")]
        for arch in r.arch.unique():
            d = r[r.arch == arch]
            best_lr = d.groupby("lr").inner_val_mse.mean().idxmin()
            for t in sorted(d.target.unique()):
                pairs = seed_mean_pairs(pairs_store, receiver, "full", arch, best_lr, cfg["seeds"], t)
                r2, lo, hi = C.tumour_bootstrap(pairs, cfg["n_boot"])
                bname, bpairs = baseline_pairs(receiver, t)
                gain, glo, ghi = C.tumour_bootstrap_diff(pairs, bpairs, cfg["n_boot"])   # R2(graph) - R2(baseline)
                per_seed = [C.tumour_bootstrap(pairs_store[(receiver, "full", arch, best_lr, sd, t)], 200)[0] for sd in cfg["seeds"]]
                rows.append({"receiver": receiver, "arch": arch, "selected_lr": best_lr, "target": t, "r2": r2, "ci_low": lo,
                             "ci_high": hi, "r2_per_seed": ";".join(f"{x:.4f}" for x in per_seed), "strongest_baseline": bname,
                             "r2_baseline": C.tumour_bootstrap(bpairs, 200)[0], "gain_vs_baseline": gain,
                             "gain_ci_low": glo, "gain_ci_high": ghi, "beats_baseline": bool(glo > 0),
                             "mean_inner_val_mse": float(d[d.lr == best_lr].inner_val_mse.mean())})
    return pd.DataFrame(rows)


def pairs_from_csv(phase: str, tag: str, receivers: list[str]) -> dict:
    """Rebuild the pairs store from the per-tumour CSVs written by run_experiments (for re-summarising)."""
    out = {}
    for r in receivers:
        df = pd.read_csv(C.RESULTS / phase / f"per_tumour_sse_{tag}_{r}.csv")
        for row in df.itertuples(index=False):
            out.setdefault((row.receiver, row.mode, row.arch, float(row.lr), int(row.seed), row.target), {})[str(row.tumour)] = (row.sse, row.sst)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "graph.yaml"))
    ap.add_argument("--summarise-only", action="store_true")
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "train" if not args.summarise_only else "train_summary", cfg)
    if args.summarise_only:
        res = pd.read_csv(C.RESULTS / PHASE / "graph_results.csv")
        pairs = pairs_from_csv(PHASE, "graph", cfg["receivers"])
    else:
        res, pairs = run_experiments(cfg, ["full"], cfg["archs"], cfg["learning_rates"], PHASE, "graph", log)
        C.save_table(res, PHASE, "graph_results.csv", log)
    summ = summarise_graph(res, pairs, cfg, log)
    C.save_table(summ, PHASE, "graph_summary.csv", log)
    # architecture selection for ablations: lower mean inner-val MSE on the primary receiver
    prim = summ[(summ.receiver == cfg["primary_receiver"])].drop_duplicates("arch")
    sel = prim.sort_values("mean_inner_val_mse").iloc[0]
    C.save_json({"selected_arch": sel.arch, "selected_lr": float(sel.selected_lr),
                 "rule": "lower mean inner-validation MSE on the primary receiver at its selected learning rate",
                 "candidates": prim[["arch", "selected_lr", "mean_inner_val_mse"]].to_dict("records")},
                PHASE, "architecture_selection.json", log)
    log.info("graph summary:\n%s", summ.round(4).to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()

"""v2 fitting engine: ridge and XGBoost-CUDA on column subsets of one design matrix, same folds and protocol as v1.

What: a `Group` is a set of receiver rows, named column subsets ("variants") of a GPU-resident design matrix, ridge
targets (multi-output) and XGBoost targets (one model per column). For every outer fold the engine standardises with
training statistics, fits ridge for every alpha (alpha is selected later by the mean inner-validation MSE over folds,
as v1 `baselines.summarise`), and fits XGBoost per grid configuration with the round count chosen on the inner
validation fold and a refit on all training folds (<= 60,000 seeded rows, as v1). Held-out predictions per
(group, fold) are saved to disk, so runs resume and every metric (overall, per bin, per subset) is computed later
from the same predictions.

Why: the v2 questions (spillover control, cross-type share) are all "same model, different columns / rows /
targets" comparisons. Ridge for many variants is solved from ONE Gram matrix per fold: centring and scaling are
per column, so the Gram matrix of a column subset is the sub-block of the full one (exact; tests/test_v2.py).
"""
from __future__ import annotations

import gc
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import features as F  # noqa: E402
import gpu as GP  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

RIDGE_ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0]
MAX_ROWS = 60000
MAX_ROUNDS = 500


@dataclass
class Group:
    name: str
    variants: dict[str, np.ndarray]                 # variant -> column indices into the design
    Y: np.ndarray                                    # ridge targets (n x k), all rows of the design
    y_names: list[str]
    xgb_targets: list[str] = field(default_factory=list)    # subset of y_names fitted with XGBoost
    xgb_variants: list[str] = field(default_factory=list)
    xgb_grid: list[tuple[float, int]] = field(default_factory=lambda: [(0.05, 15)])
    rows: np.ndarray | None = None                   # boolean mask over design rows (None = all)
    target_genes: set = field(default_factory=set)   # for the leakage guard
    target_programs: set = field(default_factory=set)


class GramRidge:
    """Centred Gram statistics of a (scaled) matrix; ridge for any column subset and alpha from sub-blocks."""

    def __init__(self, X: torch.Tensor, Y: torch.Tensor, chunk: int = 8192):
        n = len(X)
        self.xm = sum(X[i:i + chunk].double().sum(0) for i in range(0, n, chunk)) / n
        self.ym = Y.double().mean(0)
        p, k = X.shape[1], Y.shape[1]
        self.G = torch.zeros((p, p), dtype=torch.float64, device=X.device)
        self.XtY = torch.zeros((p, k), dtype=torch.float64, device=X.device)
        for i in range(0, n, chunk):
            xc = X[i:i + chunk].double() - self.xm
            self.G += xc.T @ xc
            self.XtY += xc.T @ (Y[i:i + chunk].double() - self.ym)
            del xc

    def coefs(self, cols: torch.Tensor, alphas: list[float]) -> list[torch.Tensor]:
        evals, V = torch.linalg.eigh(self.G[cols][:, cols])
        b = V.T @ self.XtY[cols]
        return [V @ (b / (evals.clamp(min=0) + a)[:, None]) for a in alphas]

    def predict(self, X: torch.Tensor, cols: torch.Tensor, W: torch.Tensor, chunk: int = 16384) -> torch.Tensor:
        xm = self.xm[cols]
        return torch.cat([((X[i:i + chunk][:, cols].double() - xm) @ W + self.ym).float() for i in range(0, len(X), chunk)])


def _xgb_params(lr: float, leaves: int, seed: int) -> dict:
    # identical to gpu.xgb_fit_select / xgb_fit_predict (v1 Phase 2 Amendment 1)
    return {"tree_method": "hist", "device": "cuda", "learning_rate": lr, "max_leaves": leaves, "max_depth": 0,
            "grow_policy": "lossguide", "min_child_weight": 50, "objective": "reg:squarederror", "seed": seed,
            "nthread": 4}   # nthread: CPU-worker cap only; does not change the fitted model


def _train_retry(params, dtrain, rounds, **kw):
    """xgb.train with one retry after releasing GPU caches if CUDA runs out of memory (same inputs, same result)."""
    import xgboost as xgb
    try:
        return xgb.train(params, dtrain, rounds, **kw)
    except xgb.core.XGBoostError as e:
        if "out of memory" not in str(e) and "cudaErrorMemoryAllocation" not in str(e):
            raise
        gc.collect()
        torch.cuda.empty_cache()
        return xgb.train(params, dtrain, rounds, **kw)


def fold_masks(section: np.ndarray, sec2fold: dict, k: int, n_folds: int):
    """(train_inner, val, train_all, test) for outer fold k; inner val = fold (k+1) mod n (v1 rule)."""
    fold_of = np.array([sec2fold[s] for s in section])
    test = fold_of == k
    val = fold_of == ((k + 1) % n_folds)
    return ~test & ~val, val, ~test, test


def leak_check(names, prov, Xhost_cols: dict[int, np.ndarray], Y, group: Group) -> None:
    """Run the v1 leakage guard on the self-sourced and tile columns of every variant (neighbour columns are
    always allowed by the guard's rules, so they are not materialised)."""
    for v, cols in group.variants.items():
        keep = [int(j) for j in cols if prov[int(j)].source in ("self_covariate", "tile")]
        if not keep:
            continue
        fs = F.FeatureSet(rows=np.arange(len(Y)))
        fs.add("checked", np.column_stack([Xhost_cols[j] for j in keep]), [names[j] for j in keep], [prov[j] for j in keep])
        F.assert_no_target_leak(fs, None, Y, set(group.target_genes), set(group.target_programs))


def run_group(X: torch.Tensor, names: list[str], prov: list, section: np.ndarray, sec2fold: dict, group: Group,
              out_dir: Path, log, fold_budget: int | None = None) -> int:
    """Fit every variant of a group on every outer fold; saves out_dir/<group>_fold<k>.npz. Returns new folds run."""
    import xgboost as xgb
    out_dir.mkdir(parents=True, exist_ok=True)
    n_folds = max(sec2fold.values()) + 1
    rows = np.ones(len(section), bool) if group.rows is None else group.rows
    # leakage guard (self/tile columns only; host copies of those few columns)
    selfcols = sorted({int(j) for cols in group.variants.values() for j in cols if prov[int(j)].source in ("self_covariate", "tile")})
    host = {j: X[:, j].cpu().numpy() for j in selfcols}
    leak_check(names, prov, host, group.Y, group)
    del host
    new = 0
    xgb_idx = [group.y_names.index(t) for t in group.xgb_targets]
    for k in range(n_folds):
        part = out_dir / f"{group.name}_fold{k}.npz"
        if part.exists():
            continue
        if fold_budget is not None and new >= fold_budget:
            log.info("%s: fold budget %d reached", group.name, fold_budget)
            return new
        t0 = time.time()
        tr_in, val, tr_all, test = (m & rows for m in fold_masks(section, sec2fold, k, n_folds))
        it, iv, ia, ie = (torch.as_tensor(np.where(m)[0], device="cuda") for m in (tr_in, val, tr_all, test))
        Yg = GP.to_gpu(group.Y)
        mu, sd = GP.fit_scaler(X[it])
        Xin, Xval = GP.apply_scaler(X[it], mu, sd), GP.apply_scaler(X[iv], mu, sd)
        mu2, sd2 = GP.fit_scaler(X[ia])
        Xtr2, Xte2 = GP.apply_scaler(X[ia], mu2, sd2), GP.apply_scaler(X[ie], mu2, sd2)
        out = {"test_idx": np.where(test)[0].astype(np.int64), "ytr_mean": group.Y[tr_all].mean(0).astype(np.float64),
               "y_names": np.array(group.y_names), "alphas": np.array(RIDGE_ALPHAS)}
        # ridge: one Gram per fold, sub-blocks per variant
        g_in, g_all = GramRidge(Xin, Yg[it]), GramRidge(Xtr2, Yg[ia])
        for v, cols in group.variants.items():
            c = torch.as_tensor(np.asarray(cols), device="cuda")
            Ws = g_in.coefs(c, RIDGE_ALPHAS)
            out[f"ridge_valmse__{v}"] = np.stack([((g_in.predict(Xval, c, W) - Yg[iv]) ** 2).mean(0).cpu().numpy() for W in Ws])
            Ws = g_all.coefs(c, RIDGE_ALPHAS)
            out[f"ridge_pred__{v}"] = np.stack([g_all.predict(Xte2, c, W).cpu().numpy() for W in Ws]).astype(np.float32)
        del g_in, g_all
        # XGBoost (v1 protocol: <= 60,000 seeded rows; rounds chosen on inner validation; refit on all training
        # folds). Low-memory order (v2 Phase 2 attempt 3): scaled matrices stay on the GPU; each DMatrix is built
        # from a host copy of one column slice that is freed at once; at most two DMatrices exist at a time
        # (selection pass: inner-train + validation; refit pass: all-train + test). Every fit has the same
        # inputs and seed as before, so the fitted models are unchanged.
        if xgb_idx and group.xgb_variants:
            rng = np.random.default_rng(C.SEED)
            sub_in = rng.choice(len(it), min(len(it), MAX_ROWS), replace=False)
            sub_all = np.random.default_rng(C.SEED).choice(len(ia), min(len(ia), MAX_ROWS), replace=False)
            Xin_s = Xin[torch.as_tensor(sub_in, device="cuda")]
            Xtr_s = Xtr2[torch.as_tensor(sub_all, device="cuda")]
            del Xin, Xtr2
            torch.cuda.empty_cache()
            yin = group.Y[np.where(tr_in)[0][sub_in]][:, xgb_idx]
            yval = group.Y[val][:, xgb_idx]
            ytr = group.Y[np.where(tr_all)[0][sub_all]][:, xgb_idx]

            def dm(T, cols):
                a = T[:, torch.as_tensor(cols, device="cuda")].cpu().numpy()
                d = xgb.DMatrix(a)
                del a
                return d
            for v in group.xgb_variants:
                cols = np.asarray(group.variants[v])
                G, J = len(group.xgb_grid), len(xgb_idx)
                preds = np.zeros((G, int(test.sum()), J), np.float32)
                vm = np.zeros((G, J))
                rounds = np.zeros((G, J), np.int32)
                d_in, d_val = dm(Xin_s, cols), dm(Xval, cols)
                for j in range(J):
                    d_in.set_label(yin[:, j]); d_val.set_label(yval[:, j])
                    for g, (lr, leaves) in enumerate(group.xgb_grid):
                        ev: dict = {}
                        p = {**_xgb_params(lr, leaves, C.SEED), "eval_metric": "rmse"}
                        _train_retry(p, d_in, MAX_ROUNDS, evals=[(d_val, "val")], evals_result=ev, verbose_eval=False)
                        rmse = np.array(ev["val"]["rmse"])
                        best = int(np.argmin(rmse)) + 1
                        vm[g, j], rounds[g, j] = float(rmse[best - 1] ** 2), best
                del d_in, d_val
                gc.collect()
                d_tr, d_te = dm(Xtr_s, cols), dm(Xte2, cols)
                for j in range(J):
                    d_tr.set_label(ytr[:, j])
                    for g, (lr, leaves) in enumerate(group.xgb_grid):
                        booster = _train_retry(_xgb_params(lr, leaves, C.SEED), d_tr, int(rounds[g, j]))
                        preds[g, :, j] = booster.predict(d_te)
                        del booster
                del d_tr, d_te
                gc.collect()
                out[f"xgb_pred__{v}"], out[f"xgb_valmse__{v}"], out[f"xgb_rounds__{v}"] = preds, vm, rounds
            out["xgb_targets"] = np.array(group.xgb_targets)
            out["xgb_grid"] = np.array(group.xgb_grid, dtype=float)
            del Xin_s, Xtr_s, Xval, Xte2
        else:
            del Xin, Xval, Xtr2, Xte2
        del Yg
        torch.cuda.empty_cache()
        gc.collect()
        tmp = part.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **out)
        tmp.replace(part)
        new += 1
        log.info("%s fold %d done in %.0f s (n_test %d; peak VRAM %.0f MiB)", group.name, k, time.time() - t0,
                 int(test.sum()), torch.cuda.max_memory_allocated() / 2**20)
    return new


# ----------------------------------------------------------------------------- evaluation from saved predictions
def load_parts(out_dir: Path, group: str, n_folds: int = 5) -> list[dict]:
    parts = []
    for k in range(n_folds):
        z = np.load(out_dir / f"{group}_fold{k}.npz", allow_pickle=False)
        parts.append({key: z[key] for key in z.files})
    return parts


def select_config(parts: list[dict], model: str, variant: str, cols: list[int]) -> int:
    """Index of the configuration (ridge alpha or XGBoost grid point) with the lowest inner-validation MSE averaged
    over folds and over the given target columns (v1 `summarise` rule)."""
    vm = np.mean([p[f"{model}_valmse__{variant}"][:, cols].mean(1) for p in parts], axis=0)
    return int(np.argmin(vm))


def pairs_for(parts: list[dict], group_arr: np.ndarray, model: str, variant: str, target: str | list[str],
              cell_mask: np.ndarray | None = None, Y: np.ndarray | None = None, config: int | None = None) -> dict:
    """Per-tumour (SSE, SST) over held-out cells (optionally restricted by a boolean mask over design rows), pooled
    over the target columns. SST uses the training-fold mean (v1 metric)."""
    names = list(parts[0]["y_names"]) if model == "ridge" else list(parts[0]["xgb_targets"])
    tl = [target] if isinstance(target, str) else list(target)
    cols = [names.index(t) for t in tl]
    ycols = [list(parts[0]["y_names"]).index(t) for t in tl]
    cfg = select_config(parts, model, variant, cols) if config is None else config
    out: dict = {}
    for p in parts:
        idx = p["test_idx"]
        pred = p[f"{model}_pred__{variant}"][cfg][:, cols]
        y = Y[idx][:, ycols]
        keep = np.ones(len(idx), bool) if cell_mask is None else cell_mask[idx]
        r = (y - pred) ** 2
        s = (y - p["ytr_mean"][ycols]) ** 2
        for g in np.unique(group_arr[idx][keep]):
            m = keep & (group_arr[idx] == g)
            a, b = out.get(str(g), (0.0, 0.0))
            out[str(g)] = (a + float(r[m].sum()), b + float(s[m].sum()))
    return out


def ratio_bootstrap(num_a: dict, num_b: dict, den_a: dict, den_b: dict, n_boot: int = 2000, seed: int = C.SEED):
    """Share = (SSE_num_b - SSE_num_a) / (SSE_den_b - SSE_den_a) summed over resampled tumours (paired).

    With num = (all, same) and den = (all, own): share = (R2(all) - R2(same)) / (R2(all) - R2(own)), because all
    pairs share the per-tumour SST, which cancels. Inputs are per-tumour (SSE, SST) dicts."""
    keys = sorted(set(num_a) & set(num_b) & set(den_a) & set(den_b))
    na = np.array([num_a[k][0] for k in keys]); nb = np.array([num_b[k][0] for k in keys])
    da = np.array([den_a[k][0] for k in keys]); db = np.array([den_b[k][0] for k in keys])
    point = (nb.sum() - na.sum()) / (db.sum() - da.sum())
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), size=(n_boot, len(keys)))
    with np.errstate(divide="ignore", invalid="ignore"):
        boots = (nb[idx].sum(1) - na[idx].sum(1)) / (db[idx].sum(1) - da[idx].sum(1))
    lo, hi = np.nanquantile(boots, [0.025, 0.975])
    return float(point), float(lo), float(hi)


def summary_row(parts, group_arr, model, variant, target, Y, ref=None, cell_mask=None, n_boot=2000) -> dict:
    """R2 with CI, and optionally the paired difference to a reference (model, variant) on the same cells."""
    pr = pairs_for(parts, group_arr, model, variant, target, cell_mask, Y)
    r2, lo, hi = C.tumour_bootstrap(pr, n_boot)
    row = {"model": model, "variant": variant, "target": target if isinstance(target, str) else "+".join(target),
           "r2": r2, "ci_low": lo, "ci_high": hi, "n_tumours": len(pr)}
    if ref is not None:
        rp = pairs_for(parts if ref[2] is None else ref[2], group_arr, ref[0], ref[1], target, cell_mask, Y)
        d, dlo, dhi = C.tumour_bootstrap_diff(pr, rp, n_boot)
        row.update({"ref": f"{ref[0]}:{ref[1]}", "diff": d, "diff_ci_low": dlo, "diff_ci_high": dhi,
                    "diff_ci_excludes_0": bool(dlo > 0 or dhi < 0)})
    return row

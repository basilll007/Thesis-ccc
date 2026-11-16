"""Chain Phase 2b: SGC-style graph model, run because the pre-declared GNN rule triggered (phase2/gnn_decision.json).

Plan (phase2 plan, GNN decision): neighbour features propagated over 1 and 2 hops of the symmetrised 6-NN graph
(row-normalised adjacency), fitted with sklearn MLPRegressor; LOTO; tumour-bootstrap delta R2 vs an MLP on the
baseline features only. Only the triggering (annotation, program, spillover) combinations are fitted.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402
import chain_features as cf  # noqa: E402
import chain_phase2_neighbourhood as p2  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

PHASE = "phase2"
ANN = "author_mapped"
MAX_TRAIN = 60000


def two_hop_no_return(W, lab, classes, S, progs, d):
    """2-hop path-weighted means over i - j - k with k != i (v1 bug: k == i paths leaked the cell's own score).

    For symmetric W with zero diagonal, the k == i contribution of sum_j W_ij W_jk x_k is deg_i * x_i.
    """
    deg = np.asarray(W.sum(axis=1)).ravel()
    n = len(lab)
    cols = []
    onehot = np.zeros((n, len(classes)))
    onehot[np.arange(n), [classes.index(x) for x in lab]] = 1
    num = W @ (W @ onehot) - deg[:, None] * onehot
    den = (W @ (W @ np.ones(n))) - deg
    cols.append(num / np.maximum(den, 1)[:, None])
    for t2, cls2 in cf.TYPES.items():
        pl = [p for p in progs if p.startswith(t2 + ":")]
        mem = (lab == cls2).astype(float)
        X = S[:, [d["score_names"].index(p) for p in pl]] * mem[:, None]
        num = W @ (W @ X) - deg[:, None] * X
        cnt = W @ (W @ mem) - deg * mem
        fill = X[mem > 0].mean(axis=0) if mem.any() else np.zeros(len(pl))
        with np.errstate(invalid="ignore", divide="ignore"):
            m = num / cnt[:, None]
        cols.append(np.where(cnt[:, None] > 0, m, fill[None, :]))
    return np.column_stack(cols)


def build_hop2(t, log):
    """Per section: 1-hop features (nb6 means, comp6) for all cells, then propagate one more hop (A_norm @ X)."""
    out = cc.CACHE / "phase2_features" / f"{ANN}__{t}__hop2_v2.npz"
    if out.exists():
        return out
    progs = cf.programs_used()
    sections, s2t = p2.sections_and_tumours()
    classes = sorted(pd.read_csv(cc.RESULTS / "phase1b" / "cell_counts_by_section.csv")
                     .query("annotation == @ANN").cell_type.unique())
    z = np.load(cc.CACHE / "phase2_features" / f"{ANN}__{t}.npz")
    sec_rows = z["section"]
    parts = []
    for s_i, s in enumerate(sections):
        d = cf.load_cache(s)
        lab = d[ANN].astype(str)
        n = len(lab)
        W = d["W"]
        deg = np.asarray(W.sum(axis=1)).ravel()
        An = sp.diags(1 / np.maximum(deg, 1)) @ W
        S = d["scores"].astype(np.float64)
        onehot = np.zeros((n, len(classes)))
        onehot[np.arange(n), [classes.index(x) for x in lab]] = 1
        feats = [An @ onehot]
        names = [f"hop2comp__{c}" for c in classes]
        for t2, cls2 in cf.TYPES.items():
            pl = [p for p in progs if p.startswith(t2 + ":")]
            mem = lab == cls2
            cols = [d["score_names"].index(p) for p in pl]
            m, _ = cf.neighbour_mean(W, S[:, cols], mem)
            fill = S[mem][:, cols].mean(axis=0) if mem.any() else np.zeros(len(cols))
            m = np.where(np.isnan(m), fill[None, :], m)
            feats.append(m)
            names += [f"hop1nb__{p}" for p in pl]
        H1 = np.column_stack(feats)
        H2 = two_hop_no_return(W, lab, classes, S, progs, d)  # v2: excludes paths returning to the cell itself
        # select the same rows as the Phase 2 feature matrix (same seeds -> same sampled cells)
        mem_t = np.where(lab == cf.TYPES[t])[0]
        rng = np.random.default_rng([cc.SEED, s_i, list(cf.TYPES).index(t)])
        idx = mem_t if len(mem_t) <= p2.CAP_PER_SECTION else np.sort(rng.choice(mem_t, p2.CAP_PER_SECTION, replace=False))
        assert len(idx) == int((sec_rows == s).sum()), s
        M = H2[idx] - H2[idx].mean(axis=0, keepdims=True)
        parts.append(M)
        log.info("hop2 %s %s done", t, s)
    np.savez_compressed(out, X=np.vstack(parts).astype(np.float32), names=np.array([f"hop2__{x}" for x in names]))
    return out


def fit_mlp(Xtr, ytr, Xte, seed):
    from sklearn.neural_network import MLPRegressor
    rng = np.random.default_rng(seed)
    if len(ytr) > MAX_TRAIN:
        k = rng.choice(len(ytr), MAX_TRAIN, replace=False)
        Xtr, ytr = Xtr[k], ytr[k]
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0)
    sd[sd == 0] = 1
    m = MLPRegressor(hidden_layer_sizes=(64, 32), alpha=1e-3, learning_rate_init=1e-3, max_iter=200,
                     early_stopping=True, validation_fraction=0.1, n_iter_no_change=10, random_state=cc.SEED)
    m.fit((Xtr - mu) / sd, ytr)
    return m.predict((Xte - mu) / sd)


def main():
    log = cc.start_run("chain_phase2b_gnn", {"max_train": MAX_TRAIN, "mlp": "(64,32), alpha 1e-3, early stopping"})
    dec = json.loads((cc.RESULTS / PHASE / "gnn_decision.json").read_text())
    rows = []
    for trig in dec["triggering_programs"]:
        prog, spill = trig["program"], trig["spill"]
        t = prog.split(":")[0]
        build_hop2(t, log)
        z = np.load(cc.CACHE / "phase2_features" / f"{ANN}__{t}.npz")
        h = np.load(cc.CACHE / "phase2_features" / f"{ANN}__{t}__hop2_v2.npz")
        names = [str(x) for x in z["names"]]
        ni = {n: i for i, n in enumerate(names)}
        y = z["X"][:, ni[f"y__{prog}"]].astype(np.float64)
        base, nbr, _ = p2.feature_sets(names, prog, spill)
        Xb = z["X"][:, [ni[c] for c in base]].astype(np.float64)
        Xn = np.column_stack([z["X"][:, [ni[c] for c in nbr]], h["X"]]).astype(np.float64)
        tum = z["tumour"]
        for f_i, held in enumerate(sorted(set(tum))):
            te, tr = tum == held, tum != held
            sst = float(np.sum(y[te] ** 2))
            for cname, X in (("baseline", Xb), ("sgc_neighbour", Xn)):
                p = fit_mlp(X[tr], y[tr], X[te], [cc.SEED, f_i])
                rows.append({"annotation": ANN, "program": prog, "spillover_covariates": spill, "config": cname,
                             "model": "sgc_mlp", "tumour": held, "n": int(te.sum()),
                             "sse": float(np.sum((y[te] - p) ** 2)), "sst": sst, "n_features": X.shape[1]})
            log.info("%s spill=%s fold %s done", prog, spill, held)
    df = pd.DataFrame(rows)
    cc.save_csv(df, PHASE, "gnn_sgc_v2_loto_per_tumour_sse.csv", "SGC-style MLP LOTO per-tumour SSE/SST (v2: no return-to-self 2-hop paths)", log)
    rng = np.random.default_rng(cc.SEED)
    res, prev = [], pd.read_csv(cc.RESULTS / PHASE / "delta_r2_by_program.csv")
    for (prog, spill), d in df.groupby(["program", "spillover_covariates"]):
        pt, ci, per_t, r2b, r2n = p2.boot_delta(d[d.config == "baseline"], d[d.config == "sgc_neighbour"], rng)
        pv = prev[(prev.annotation == ANN) & (prev.program == prog) & (prev.spillover_covariates == spill)].set_index("model")
        res.append({"annotation": ANN, "program": prog, "spillover_covariates": spill, "model": "sgc_mlp",
                    "r2_baseline": r2b, "r2_neighbour": r2n, "delta_r2": pt, "ci_low": ci[0], "ci_high": ci[1],
                    "n_tumours_delta_pos": int((per_t > 0).sum()), "ridge_delta_r2": pv.loc["ridge", "delta_r2"],
                    "hgb_delta_r2": pv.loc["hgb", "delta_r2"], "r2_neighbour_ridge": pv.loc["ridge", "r2_neighbour"],
                    "r2_neighbour_hgb": pv.loc["hgb", "r2_neighbour"]})
    res = pd.DataFrame(res)
    cc.save_csv(res, PHASE, "gnn_sgc_v2_delta_r2.csv", "SGC-style graph model delta R2 vs ridge/HGB (v2, leakage fixed; tumour bootstrap CI)", log)
    log.info("\n%s", res.to_string(index=False))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

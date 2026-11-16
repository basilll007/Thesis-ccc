"""Chain Phase 2: does the neighbourhood predict a cell's program state beyond own covariates?

Plan: docs/chain_plans/phase2_neighbourhood_prediction.md (committed before running).
Steps: (1) build per-(annotation, receiver type) feature matrices one section at a time (cached);
(2) leave-one-tumour-out ridge + gradient boosting per receiver program (parallel tasks, max 6 workers,
per-task outputs saved for resume); (3) tumour-bootstrap delta R2, feature-group attribution, GNN decision, figure.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "5")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "5")
os.environ.setdefault("MKL_NUM_THREADS", "5")

import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402
import chain_features as cf  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase2"
CAP_PER_SECTION = 20000
HGB_MAX_TRAIN = 60000
ALPHAS = np.logspace(-2, 4, 13)
N_BOOT = 2000
ANNOTS = ["author_mapped", "marker"]
FEAT_DIR = cc.CACHE / "phase2_features"


def sections_and_tumours():
    import v3_nulls as vn
    s2t = vn.section_to_tumour()
    return sorted(s2t), s2t


# ----------------------------------------------------------------------------- features
def build_features(ann: str, log):
    progs = cf.programs_used()
    sections, s2t = sections_and_tumours()
    classes = sorted(pd.read_csv(cc.RESULTS / "phase1b" / "cell_counts_by_section.csv")
                     .query("annotation == @ann").cell_type.unique())
    by_type = {t: [p for p in progs if p.startswith(t + ":")] for t in cf.TYPES}
    todo = [t for t in cf.TYPES if by_type[t] and not (FEAT_DIR / f"{ann}__{t}.npz").exists()]
    if not todo:
        return
    parts = {t: [] for t in todo}
    for s_i, s in enumerate(sections):
        d = cf.load_cache(s)
        lab = d[ann].astype(str)
        n = len(lab)
        W, R = d["W"], d["R"]
        S = d["scores"].astype(np.float64)
        sidx = {p: d["score_names"].index(p) for p in progs}
        deg6 = np.asarray(W.sum(axis=1)).ravel()
        onehot = np.zeros((n, len(classes)), np.float32)
        onehot[np.arange(n), [classes.index(x) for x in lab]] = 1
        comp6 = (W @ onehot) / np.maximum(deg6, 1)[:, None]
        compr = R @ onehot
        dens = np.asarray(R.sum(axis=1)).ravel()
        tiles = cf.tile_codes(d["xy"])
        nb = {}
        for t2, cls2 in cf.TYPES.items():
            if not by_type[t2]:
                continue
            mem = lab == cls2
            cols = [sidx[p] for p in by_type[t2]]
            fill = S[mem][:, cols].mean(axis=0) if mem.any() else np.zeros(len(cols))
            for gname, A in (("nb6", W), ("nbr30", R)):
                m, _ = cf.neighbour_mean(A, S[:, cols], mem)
                m = np.where(np.isnan(m), fill[None, :], m)
                for j, p in enumerate(by_type[t2]):
                    nb[f"{gname}__{t2}__{p.split(':')[1]}"] = m[:, j]
        for t in todo:
            mem = lab == cf.TYPES[t]
            idx = np.where(mem)[0]
            if len(idx) == 0:
                continue
            rng = np.random.default_rng([cc.SEED, s_i, list(cf.TYPES).index(t)])
            if len(idx) > CAP_PER_SECTION:
                idx = np.sort(rng.choice(idx, CAP_PER_SECTION, replace=False))
            cols, names = [], []
            for p in by_type[t]:
                cols.append(S[idx, sidx[p]]); names.append(f"y__{p}")
                tm, has = cf.group_mean_loo(S[:, sidx[p]], tiles, mem)
                secmean = S[mem, sidx[p]].mean()
                cols.append(np.where(has, tm, secmean)[idx]); names.append(f"base__tilemean__{p}")
                cols.append((~has).astype(float)[idx]); names.append(f"base__tileempty__{p}")
            cols += [np.log1p(d["total_counts"][idx].astype(float)), dens[idx]]
            names += ["base__log_total_counts", "base__density_r30"]
            for j, c in enumerate(classes):
                cols.append(comp6[idx, j]); names.append(f"comp6__{c}")
                cols.append(compr[idx, j]); names.append(f"compr30__{c}")
            for k, v in nb.items():
                cols.append(v[idx]); names.append(k)
            for j, ln in enumerate(d["lineage_names"]):
                cols.append(d["lineage"][idx, j].astype(float)); names.append(f"spill__{ln}")
            M = np.column_stack(cols).astype(np.float64)
            M -= M.mean(axis=0, keepdims=True)  # section fixed effect (demean within section)
            parts[t].append((s, M, names))
        log.info("features %s %s built: %s", ann, s, cc.mem_status())
        del d, W, R, S
    FEAT_DIR.mkdir(parents=True, exist_ok=True)
    for t in todo:
        names = parts[t][0][2]
        assert all(p[2] == names for p in parts[t])
        X = np.vstack([p[1] for p in parts[t]]).astype(np.float32)
        sec = np.concatenate([[p[0]] * len(p[1]) for p in parts[t]])
        np.savez_compressed(FEAT_DIR / f"{ann}__{t}.npz", X=X, names=np.array(names), section=sec,
                            tumour=np.array([s2t[x] for x in sec]))
        log.info("saved features %s__%s: %s rows x %d cols", ann, t, X.shape[0], X.shape[1])


# ----------------------------------------------------------------------------- models
def _ridge_fit_predict(Xtr, ytr, Xte):
    from sklearn.linear_model import RidgeCV
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0)
    sd[sd == 0] = 1
    m = RidgeCV(alphas=ALPHAS).fit((Xtr - mu) / sd, ytr)
    return m.predict((Xte - mu) / sd), float(m.alpha_)


def _hgb_fit_predict(Xtr, ytr, Xte, seed):
    from sklearn.ensemble import HistGradientBoostingRegressor
    rng = np.random.default_rng(seed)
    if len(ytr) > HGB_MAX_TRAIN:
        k = rng.choice(len(ytr), HGB_MAX_TRAIN, replace=False)
        Xtr, ytr = Xtr[k], ytr[k]
    m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, early_stopping=True,
                                      validation_fraction=0.1, random_state=cc.SEED).fit(Xtr, ytr)
    return m.predict(Xte)


def feature_sets(names, prog, spill: bool):
    t = prog.split(":")[0]
    base = [f"base__tilemean__{prog}", f"base__tileempty__{prog}", "base__log_total_counts", "base__density_r30"]
    if spill:
        base += [n for n in names if n.startswith("spill__")]
    groups = {"composition": [n for n in names if n.startswith(("comp6__", "compr30__"))],
              "same_type_programs": [n for n in names if n.startswith(("nb6__", "nbr30__")) and n.split("__")[1] == t]}
    for t2 in cf.TYPES:
        if t2 != t:
            g = [n for n in names if n.startswith(("nb6__", "nbr30__")) and n.split("__")[1] == t2]
            if g:
                groups[f"other_type_programs__{t2}"] = g
    nbr = base + [n for g in groups.values() for n in g]
    return base, nbr, groups


def run_task(ann, prog, spill):
    """LOTO for one (annotation, receiver program, spillover flag). Writes a per-task CSV (resume marker)."""
    out = cc.out_dir(f"{PHASE}/tasks") / f"{ann}__{prog.replace(':', '-')}__spill{int(spill)}.csv"
    if out.exists():
        return str(out)
    t = prog.split(":")[0]
    z = np.load(FEAT_DIR / f"{ann}__{t}.npz")
    X, names, tum = z["X"], [str(x) for x in z["names"]], z["tumour"]
    ni = {n: i for i, n in enumerate(names)}
    y = X[:, ni[f"y__{prog}"]].astype(np.float64)
    base, nbr, groups = feature_sets(names, prog, spill)
    configs = {"baseline": base, "neighbour": nbr}
    if ann == "author_mapped":
        for g, cols in groups.items():
            configs[f"ridge_only__drop__{g}"] = [c for c in nbr if c not in cols]
            configs[f"ridge_only__add__{g}"] = base + cols
    rows = []
    for f_i, held in enumerate(sorted(set(tum))):
        te, tr = tum == held, tum != held
        sst = float(np.sum(y[te] ** 2))  # within-section demeaned outcome
        for cname, cols in configs.items():
            ci = [ni[c] for c in cols]
            Xtr, Xte = X[tr][:, ci].astype(np.float64), X[te][:, ci].astype(np.float64)
            p, alpha = _ridge_fit_predict(Xtr, y[tr], Xte)
            rows.append({"annotation": ann, "program": prog, "spillover_covariates": spill, "config": cname,
                         "model": "ridge", "tumour": held, "n": int(te.sum()), "sse": float(np.sum((y[te] - p) ** 2)),
                         "sst": sst, "alpha": alpha, "n_features": len(ci)})
            if cname in ("baseline", "neighbour"):
                p = _hgb_fit_predict(Xtr, y[tr], Xte, [cc.SEED, f_i])
                rows.append({"annotation": ann, "program": prog, "spillover_covariates": spill, "config": cname,
                             "model": "hgb", "tumour": held, "n": int(te.sum()),
                             "sse": float(np.sum((y[te] - p) ** 2)), "sst": sst, "alpha": np.nan,
                             "n_features": len(ci)})
    df = pd.DataFrame(rows)
    tmp = out.with_suffix(".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(out)
    return str(out)


# ----------------------------------------------------------------------------- aggregation
def boot_delta(d_base: pd.DataFrame, d_nbr: pd.DataFrame, rng):
    a = d_base.set_index("tumour")[["sse", "sst"]]
    b = d_nbr.set_index("tumour")[["sse"]].loc[a.index]
    sb, sn, st = a.sse.to_numpy(), b.sse.to_numpy(), a.sst.to_numpy()
    point = (sb.sum() - sn.sum()) / st.sum()
    k = rng.integers(0, len(st), size=(N_BOOT, len(st)))
    boots = (sb[k].sum(1) - sn[k].sum(1)) / st[k].sum(1)
    per_t = (sb - sn) / st
    r2b, r2n = 1 - sb.sum() / st.sum(), 1 - sn.sum() / st.sum()
    return point, np.quantile(boots, [.025, .975]), per_t, r2b, r2n


def aggregate(log):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tasks = pd.concat([pd.read_csv(f) for f in sorted((cc.RESULTS / PHASE / "tasks").glob("*.csv"))], ignore_index=True)
    rng = np.random.default_rng(cc.SEED)
    rows, arows = [], []
    for (ann, prog, spill), d in tasks.groupby(["annotation", "program", "spillover_covariates"]):
        for model in ("ridge", "hgb"):
            b = d[(d.config == "baseline") & (d.model == model)]
            n = d[(d.config == "neighbour") & (d.model == model)]
            pt, ci, per_t, r2b, r2n = boot_delta(b, n, rng)
            rows.append({"annotation": ann, "program": prog, "spillover_covariates": spill, "model": model,
                         "n_cells": int(b.n.sum()), "n_tumours": len(b), "r2_baseline": r2b, "r2_neighbour": r2n,
                         "delta_r2": pt, "ci_low": ci[0], "ci_high": ci[1],
                         "n_tumours_delta_pos": int((per_t > 0).sum()),
                         "median_tumour_delta_r2": float(np.median(per_t)),
                         "neighbours_predict_beyond_own": bool(ci[0] > 0)})
        if ann == "author_mapped":
            nb = d[(d.config == "neighbour") & (d.model == "ridge")]
            base = d[(d.config == "baseline") & (d.model == "ridge")]
            for cfg in sorted(c for c in d.config.unique() if c.startswith("ridge_only__")):
                _, kind, grp = cfg.split("__", 2)
                x = d[(d.config == cfg) & (d.model == "ridge")]
                if kind == "drop":  # R2 lost when dropping the group from the full neighbour model
                    pt, ci, _, _, _ = boot_delta(x, nb, rng)
                else:  # R2 gained when adding only this group to baseline
                    pt, ci, _, _, _ = boot_delta(base, x, rng)
                arows.append({"annotation": ann, "program": prog, "spillover_covariates": spill, "kind": kind,
                              "group": grp, "delta_r2": pt, "ci_low": ci[0], "ci_high": ci[1]})
    res = pd.DataFrame(rows)
    att = pd.DataFrame(arows)
    cc.save_csv(tasks, PHASE, "loto_per_tumour_sse.csv", "Per held-out tumour SSE/SST for every config/model", log)
    cc.save_csv(res, PHASE, "delta_r2_by_program.csv",
                "Held-out (LOTO) delta R2 neighbour vs baseline, tumour-bootstrap 95% CI (2000)", log)
    cc.save_csv(att, PHASE, "feature_group_attribution.csv",
                "Ridge drop-group (R2 lost) and add-group (R2 gained) per neighbour feature group, tumour-bootstrap CI", log)
    # GNN decision
    prim = res[(res.annotation == "author_mapped")]
    trig = []
    for (prog, spill), d in prim.groupby(["program", "spillover_covariates"]):
        r, h = d.set_index("model").loc["ridge"], d.set_index("model").loc["hgb"]
        if max(r.ci_low, h.ci_low) > 0.01 and h.delta_r2 - r.delta_r2 >= 0.01:
            trig.append({"program": prog, "spill": bool(spill), "ridge": r.delta_r2, "hgb": h.delta_r2})
    gnn = {"rule": "GNN only if some program has delta R2 lower CI > 0.01 AND hgb delta R2 >= ridge + 0.01 (author_mapped)",
           "triggered": bool(trig), "triggering_programs": trig,
           "decision": "GNN warranted (SGC-style per plan)" if trig else "GNN not warranted; not fitted"}
    cc.save_json(gnn, PHASE, "gnn_decision.json", "GNN decision per pre-declared rule", log)
    log.info("GNN decision: %s", gnn)
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, ann in zip(axes, ANNOTS):
        d = res[res.annotation == ann].sort_values(["program", "spillover_covariates", "model"])
        progs = sorted(d.program.unique())
        for k, (model, spill, mk, col) in enumerate([("ridge", False, "o", "#1f77b4"), ("hgb", False, "s", "#ff7f0e"),
                                                     ("ridge", True, "o", "#9ecae1"), ("hgb", True, "s", "#fdd0a2")]):
            e = d[(d.model == model) & (d.spillover_covariates == spill)].set_index("program").loc[progs]
            x = np.arange(len(progs)) + (k - 1.5) * 0.18
            ax.errorbar(x, e.delta_r2, yerr=[e.delta_r2 - e.ci_low, e.ci_high - e.delta_r2], fmt=mk, color=col,
                        label=f"{model}{' + spillover cov.' if spill else ''}", ms=5, capsize=2)
        ax.axhline(0, color="k", lw=.7)
        ax.set_xticks(range(len(progs)), progs, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"{ann}: held-out delta R2 (LOTO, 95% tumour-bootstrap CI)", fontsize=9)
    axes[0].set_ylabel("delta R2 (neighbour - baseline), within-section")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, "phase2_delta_r2", "Phase 2 held-out delta R2 per receiver program", log)
    summ = {"n_receiver_programs": int(res.program.nunique()), "n_combinations": len(res),
            "combinations": "programs x annotation(2) x spillover(2) x model(2)",
            "author_mapped_no_spill": prim[~prim.spillover_covariates][["program", "model", "delta_r2", "ci_low",
                                                                        "ci_high"]].to_dict("records"),
            "gnn": gnn["decision"]}
    cc.save_json(summ, PHASE, "phase2_summary.json", "Phase 2 key numbers", log)
    return res


def main():
    log = cc.start_run("chain_phase2_neighbourhood", {"cap_per_section": CAP_PER_SECTION, "hgb_max_train": HGB_MAX_TRAIN,
                                                      "alphas": ALPHAS.tolist(), "n_boot": N_BOOT})
    for ann in ANNOTS:
        build_features(ann, log)
    progs = cf.programs_used()
    tasks = [(a, p, s) for a in ANNOTS for p in progs for s in (False, True)]
    n_workers = int(os.environ.get("CHAIN_WORKERS", "6"))
    with ProcessPoolExecutor(max_workers=min(6, n_workers)) as ex:
        futs = {ex.submit(run_task, *t): t for t in tasks}
        for f in as_completed(futs):
            log.info("task %s done -> %s; %s", futs[f], f.result(), cc.mem_status())
    for f in sorted((cc.RESULTS / PHASE / "tasks").glob("*.csv")):
        cc.register(f, PHASE, "Per-task LOTO SSE/SST (resume unit)")
    res = aggregate(log)
    log.info("delta R2 (author_mapped):\n%s", res[res.annotation == "author_mapped"].to_string(index=False))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

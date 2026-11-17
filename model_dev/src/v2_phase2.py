"""v2 Phase 2: cross-type share of the neighbourhood gain on the aggregated baselines (plan: plans/v2_phase2_crosstype.md).

Subcommands:
  fit --run {myeloid,tumour,adjusted}   one run at a time on the GPU (resumable per group x fold)
  summarise                             quantities with paired tumour-bootstrap CIs, agreement with the v1 graph
                                        estimate, reproduction check against v1 Phase 2, figure data
Variants split the v1 feature columns by the neighbour cell type they describe (v2common.column_type). Ridge is
evaluated on `programs`-set columns (r_* variants), XGBoost on `genes`-set columns (x_* variants) with the full v1
grid. Development set only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import v2common as V  # noqa: E402
import features as F  # noqa: E402
import assemble as A  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "v2_phase2"
OUT = C.RESULTS / PHASE
CACHE = V.V2_CACHE / "phase2"
P1 = C.RESULTS / "v2_phase1"
XGB_GRID = [(0.05, 15), (0.05, 31), (0.1, 15), (0.1, 31)]   # v1 Phase 2 grid (lr x max_leaves)
SAME = {"myeloid": "myeloid", "tumour": "tumour"}
OTHERS = {"myeloid": ["CAF", "tumour", "T", "other"], "tumour": ["CAF", "myeloid", "T", "other"]}
PROGRAMS_BLOCKS = {"own", "tile", "composition", "nb_programs"}


def variant_columns(names: list[str], block_of: list[str], receiver: str, extra_own: list[int]) -> dict[str, np.ndarray]:
    """r_* (programs set) and x_* (genes set) variants. `extra_own`: columns always kept (e.g. s_c, new tile means)."""
    same = SAME[receiver]
    ctype = [V.column_type(n) for n in names]
    kind = [V.column_kind(n) for n in names]
    base = [j for j, b in enumerate(block_of) if b in ("own", "tile")] + list(extra_own)
    out = {}
    for pfx, allowed in (("r", PROGRAMS_BLOCKS), ("x", PROGRAMS_BLOCKS | {"nb_genes"})):
        nb = [j for j, b in enumerate(block_of) if b in allowed and b not in ("own", "tile")]
        sel = {
            "all": nb,
            "same": [j for j in nb if ctype[j] == same],
            "other_only": [j for j in nb if ctype[j] != same],
            "same_plus_other_comp": [j for j in nb if ctype[j] == same or kind[j] == "composition"],
        }
        for t in OTHERS[receiver]:
            sel[f"minus_{t}"] = [j for j in nb if ctype[j] != t]
        out[f"{pfx}_own"] = np.array(base)
        for k, v in sel.items():
            out[f"{pfx}_{k}"] = np.array(base + v)
    return out


def load(run: str, log):
    """Design and groups for one run."""
    import torch
    import gpu as GP
    import v2engine as E
    receiver = "tumour" if run == "tumour" else "myeloid"
    data = A.build_receiver("development", receiver)
    spec = A.targets_spec()[receiver]
    blocks, names, prov, block_of = [], [], [], []
    for b in ["own", "tile", "composition", "nb_programs", "nb_genes"]:
        blocks.append(data.X[b]); names += data.names[b]; prov += data.prov[b]; block_of += [b] * len(data.names[b])
        data.X[b] = None
    groups = []
    if run in ("myeloid", "tumour"):
        Y, y_names = data.Y, list(data.target_names)
        vv = variant_columns(names, block_of, receiver, [])
        groups.append(E.Group(f"{run}_orig", vv, Y, y_names, list(spec["programs"]), [k for k in vv if k.startswith("x_")],
                              XGB_GRID, target_genes=set(spec["genes"]), target_programs=set(spec["programs"])))
    else:
        r = np.load(V.V2_CACHE / "phase1_receivers.npz", allow_pickle=False)
        sets = json.loads((P1 / "gene_sets_and_adjusted_targets.json").read_text())
        build = json.loads((P1 / "receiver_build.json").read_text())
        nn = [str(x) for x in r["new_target_names"]]
        want = [t for t in ["C1Q_TAM_adj", "SPP1_TAM_adj", "C1Q_TAM_mean", "SPP1_TAM_mean"] if t in nn]
        ti = [nn.index(t) for t in want]
        j0 = len(names)
        blocks.append(r["new_tile"][:, ti].astype(np.float32)); names += [f"tile_loo_mean__{t}" for t in want]
        prov += [F.Provenance("tile", genes=frozenset(build["new_targets"][t])) for t in want]; block_of += ["tile_new"] * len(want)
        blocks.append(r["s_c"].astype(np.float32)[:, None]); names.append("s_c_contamination")
        prov.append(F.Provenance("self_covariate", genes=frozenset(sets["M"]["genes"]))); block_of.append("s_c")
        tile_new, sc = list(range(j0, j0 + len(want))), [j0 + len(want)]
        gi = {g: i for i, g in enumerate(data.target_names)}
        adj_t = [t for t in ["C1Q_TAM_adj", "SPP1_TAM_adj"] if t in want]
        vec_adj = [g for g in sets["gene_vector_adj"] if g in gi]
        Yadj = np.column_stack([r["new_targets"][:, nn.index(t)] for t in adj_t] + [data.Y[:, gi[g]] for g in vec_adj]).astype(np.float32)
        vv = variant_columns(names, block_of, receiver, tile_new + sc)
        genes_adj = {g for t in adj_t for g in build["new_targets"][t]} | set(vec_adj)
        groups.append(E.Group("myeloid_adj_sc", vv, Yadj, adj_t + vec_adj, adj_t, [k for k in vv if k.startswith("x_")],
                              XGB_GRID, target_genes=genes_adj, target_programs=set()))
        ref_t = [t for t in ["C1Q_TAM_mean", "SPP1_TAM_mean"] if t in want]
        Yref = np.column_stack([r["new_targets"][:, nn.index(t)] for t in ref_t]).astype(np.float32)
        vr = variant_columns(names, block_of, receiver, tile_new)
        groups.append(E.Group("myeloid_mean_ref", vr, Yref, ref_t, ref_t, [k for k in vr if k.startswith("x_")],
                              XGB_GRID, target_genes={g for t in ref_t for g in build["new_targets"][t]}, target_programs=set()))
        log.info("adjusted run: targets %s, gene_vector_adj %s; reference %s", adj_t, vec_adj, ref_t)
    Xg = torch.cat([GP.to_gpu(b) for b in blocks], 1)
    del blocks
    log.info("%s design: %d x %d on GPU", run, Xg.shape[0], Xg.shape[1])
    return Xg, names, prov, data, groups


def fit(args) -> None:
    import v2engine as E
    log = C.start_run(PHASE, f"v2_phase2_fit_{args.run}")
    V.memory_gate(log)
    Xg, names, prov, data, groups = load(args.run, log)
    sec2fold = dict(zip(A.folds().section, A.folds().fold))
    for g in groups:
        log.info("group %s: %s; XGBoost targets %s", g.name, {k: len(v) for k, v in g.variants.items()}, g.xgb_targets)
        E.run_group(Xg, names, prov, data.section, sec2fold, g, CACHE, log, args.fold_budget)
    log.info("fit %s finished", args.run)


# ----------------------------------------------------------------------------- summarise
def group_targets(gname: str, ynames: list[str]) -> np.ndarray:
    receiver = "tumour" if gname.startswith("tumour") else "myeloid"
    z = np.load(C.CACHE / "features" / f"development__{receiver}.npz", allow_pickle=False)
    spec = A.targets_spec()[receiver]
    names = spec["programs"] + spec["genes"]
    cols = {t: z["Y"][:, names.index(t)] for t in names}
    if receiver == "myeloid":
        r = np.load(V.V2_CACHE / "phase1_receivers.npz", allow_pickle=False)
        for i, t in enumerate(r["new_target_names"]):
            cols[str(t)] = r["new_targets"][:, i]
    return np.column_stack([cols[t] for t in ynames]).astype(np.float32)


def quantities(parts, grp, Y, model: str, pfx: str, target, receiver: str, own_variant: str, n_boot: int = 2000) -> list[dict]:
    import v2engine as E
    # configuration selection exactly as v1 `summarise`: XGBoost jointly over the group's program targets, ridge
    # jointly over all target columns of the group (lowest inner-validation MSE averaged over folds)
    nx = len(parts[0]["xgb_targets"]) if "xgb_targets" in parts[0] else 0
    ncols = list(range(nx if model == "xgb" else len(parts[0]["y_names"])))

    def pp(m, v):
        cfg = E.select_config(parts, m, v, ncols if m == model else list(range(len(parts[0]["y_names"]))))
        return E.pairs_for(parts, grp, m, v, target, Y=Y, config=cfg)
    P = {v: pp(model, f"{pfx}_{v}") for v in
         ["all", "same", "other_only", "same_plus_other_comp"] + [f"minus_{t}" for t in OTHERS[receiver]]}
    own = pp("ridge", own_variant)
    rows = []

    def add(q, a, b):
        d, lo, hi = C.tumour_bootstrap_diff(a, b, n_boot)
        rows.append({"quantity": q, "value": d, "ci_low": lo, "ci_high": hi, "ci_excludes_0": bool(lo > 0 or hi < 0)})
    for v, pr in P.items():
        r2, lo, hi = C.tumour_bootstrap(pr, n_boot)
        rows.append({"quantity": f"R2 {v}", "value": r2, "ci_low": lo, "ci_high": hi, "ci_excludes_0": None})
    r2, lo, hi = C.tumour_bootstrap(own, n_boot)
    rows.append({"quantity": "R2 own (ridge own+tile)", "value": r2, "ci_low": lo, "ci_high": hi, "ci_excludes_0": None})
    add("neighbour gain (all - own)", P["all"], own)
    add("cross-type contribution (all - same)", P["all"], P["same"])
    add("same-type contribution (all - other_only)", P["all"], P["other_only"])
    add("same-type only gain (same - own)", P["same"], own)
    add("who is nearby: other-type composition (same_plus_other_comp - same)", P["same_plus_other_comp"], P["same"])
    add("what they express: other-type expression (all - same_plus_other_comp)", P["all"], P["same_plus_other_comp"])
    for t in OTHERS[receiver]:
        add(f"unique contribution of {t} (all - minus_{t})", P["all"], P[f"minus_{t}"])
    gain_row = [r for r in rows if r["quantity"].startswith("neighbour gain")][0]
    if gain_row["ci_low"] > 0 or gain_row["ci_high"] < 0:
        s, lo, hi = E.ratio_bootstrap(P["all"], P["same"], P["all"], own, n_boot)
        rows.append({"quantity": "cross-type share ((all - same) / (all - own))", "value": s, "ci_low": lo, "ci_high": hi,
                     "ci_excludes_0": bool(lo > 0 or hi < 0)})
    else:
        rows.append({"quantity": "cross-type share ((all - same) / (all - own))", "value": np.nan, "ci_low": np.nan,
                     "ci_high": np.nan, "ci_excludes_0": None, "note": "not defined: neighbour-gain CI includes 0"})
    return rows


def summarise(args) -> None:
    import v2engine as E
    log = C.start_run(PHASE, "v2_phase2_summarise")
    all_rows = []
    specs = [("myeloid_orig", "myeloid", "", "r_own"), ("tumour_orig", "tumour", "", "r_own"),
             ("myeloid_adj_sc", "myeloid", "adjusted (s_c covariate + target-gene removal)", "r_own"),
             ("myeloid_mean_ref", "myeloid", "unadjusted gene-mean reference", "r_own")]
    for gname, receiver, label, own_v in specs:
        if not (CACHE / f"{gname}_fold4.npz").exists():
            log.warning("%s: not complete; skipped", gname)
            continue
        parts = E.load_parts(CACHE, gname)
        ynames = list(parts[0]["y_names"])
        Y = group_targets(gname, ynames)
        z = np.load(C.CACHE / "features" / f"development__{receiver}.npz", allow_pickle=False)
        grp = z["group"].astype(str)
        spec = A.targets_spec()[receiver]
        progs = [t for t in ynames if t in spec["programs"] or t.endswith(("_adj", "_mean"))]
        genes = [t for t in ynames if t not in progs]
        tsets = [(t, t) for t in progs] + ([("gene_vector" if gname.endswith("orig") else "gene_vector_adj", genes)] if genes else [])
        xt = list(parts[0]["xgb_targets"])
        for tname, tcols in tsets:
            for model, pfx in (("xgb", "x"), ("ridge", "r"), ("ridge", "x")):
                if model == "xgb" and (not isinstance(tcols, str) or tcols not in xt):
                    continue
                for row in quantities(parts, grp, Y, model, pfx, tcols, receiver, own_v):
                    all_rows.append({"group": gname, "receiver": receiver, "label": label, "target": tname,
                                     "model": f"{model} ({'genes' if pfx == 'x' else 'programs'} set)", **row})
        log.info("summarised %s", gname)
    res = pd.DataFrame(all_rows)
    C.save_table(res, PHASE, "phase2_quantities.csv", log)
    agreement(res, log)
    reproduction(res, log)


def agreement(res: pd.DataFrame, log) -> None:
    g = pd.read_csv(C.RESULTS / "phase3" / "ablation_summary.csv")
    g = g[(g["mode"] == "same_type") & (g.arch == "sage")]
    pairs = [("myeloid", "myeloid:C1Q_TAM", "myeloid:C1Q_TAM", "xgb (genes set)"),
             ("myeloid", "myeloid:SPP1_TAM", "myeloid:SPP1_TAM", "xgb (genes set)"),
             ("myeloid", "gene_vector", "genes_pooled", "ridge (programs set)"),
             ("tumour", "tumour:classical", "tumour:classical", "xgb (genes set)")]
    rows = []
    for rec, t, gt, model in pairs:
        gr = g[(g.receiver == rec) & (g.target == gt)]
        a = res[(res.group == f"{rec}_orig") & (res.target == t) & (res.model == model) &
                (res.quantity == "cross-type contribution (all - same)")]
        if gr.empty or a.empty:
            continue
        gr, a = gr.iloc[0], a.iloc[0]
        gv, glo, ghi = gr.full_minus_mode, gr.diff_ci_low, gr.diff_ci_high
        agree = (a.ci_low <= gv <= a.ci_high) and (glo <= a.value <= ghi)
        rows.append({"receiver": rec, "target": t, "aggregated_model": model, "aggregated_cross_type": a.value,
                     "aggregated_ci_low": a.ci_low, "aggregated_ci_high": a.ci_high, "graph_cross_type": gv,
                     "graph_ci_low": glo, "graph_ci_high": ghi, "agrees": bool(agree),
                     "direction_if_not": "" if agree else ("larger" if a.value > gv else "smaller")})
    C.save_table(pd.DataFrame(rows), PHASE, "agreement_with_graph_estimate.csv", log)


def reproduction(res: pd.DataFrame, log) -> None:
    v1 = pd.read_csv(C.RESULTS / "phase2" / "baseline_summary.csv")
    rows = []
    for rec in ("myeloid", "tumour"):
        for model, fset, tol, mlabel in (("xgb", "genes", 0.005, "xgb (genes set)"), ("ridge", "programs", 0.001, "ridge (programs set)")):
            for t in v1[(v1.receiver == rec) & (v1.model == model) & (v1.feature_set == fset)].target.unique():
                tt = "gene_vector" if t == "genes_pooled" else t
                a = res[(res.group == f"{rec}_orig") & (res.target == tt) & (res.model == mlabel) & (res.quantity == "R2 all")]
                if a.empty:
                    continue
                ref = v1[(v1.receiver == rec) & (v1.model == model) & (v1.feature_set == fset) & (v1.target == t)].r2.iloc[0]
                rows.append({"receiver": rec, "model": mlabel, "target": tt, "v1_r2": ref, "v2_r2": a.value.iloc[0],
                             "abs_diff": abs(ref - a.value.iloc[0]), "tolerance": tol, "pass": abs(ref - a.value.iloc[0]) <= tol})
    df = pd.DataFrame(rows)
    C.save_table(df, PHASE, "reproduction_check_vs_v1.csv", log)
    if len(df):
        log.info("reproduction: %d/%d within tolerance", int(df["pass"].sum()), len(df))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["fit", "summarise"])
    ap.add_argument("--run", choices=["myeloid", "tumour", "adjusted"])
    ap.add_argument("--fold-budget", type=int)
    args = ap.parse_args()
    {"fit": fit, "summarise": summarise}[args.step](args)


if __name__ == "__main__":
    main()

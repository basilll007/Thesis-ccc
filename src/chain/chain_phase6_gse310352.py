"""Chain Phase 6: consistency check of Phases 2-3 in GSE310352 (marker labels, NN-scaled tiles).

Plan: docs/chain_plans/phase6_consistency_gse310352.md (committed before running).
Patient independence from Mendeley is UNVERIFIED -> consistency check, not replication.
Steps: (1) per slide: load, QC, align v3 marker labels, score Mendeley Amendment-1 programs, type-aware background,
graphs scaled by median NN px, cache; (2) Phase 2 repeat (leave-one-unit-out); (3) Phase 3 repeat for forward chains
passing the floor in >= 3 units (100 perms per null); (4) same-direction comparison with Mendeley.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402
import chain_engine as ce  # noqa: E402
import chain_features as cf  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

PHASE = "phase6"
SLIDES = ["slide1", "slide2", "slide3", "slide4", "slide5", "slide5b2", "slide5b3", "slide6"]
UNIT = {s: ("slide5_group" if s.startswith("slide5") else s) for s in SLIDES}  # UNVERIFIED grouping (plan)
CACHE_SUB = "gse310352"
N_PERM = 100
N_BOOT = 2000
ANN = "marker"


def mendeley_nn_um():
    return float(np.median(list(json.loads((cc.RESULTS / "phase1b" / "phase1_summary.json").read_text())
                                ["median_nn_um_by_section"].values())))


def prep_slide(slide, log):
    import scanpy as sc
    from sklearn.neighbors import NearestNeighbors
    import chain_phase1_cell_states as p1
    import pdac_loaders
    from phase3_pdac_audit import qc
    cache = cc.CACHE / CACHE_SUB / f"{slide}.npz"
    sdir = cc.out_dir(f"{PHASE}/slides")
    if cache.exists() and (sdir / f"{slide}_info.json").exists():
        return json.loads((sdir / f"{slide}_info.json").read_text())
    adata, meta = pdac_loaders.load(f"geo_gse310352:{slide}", cc.DATA_DIR, log)
    adata, _, _ = qc(adata, log)
    labs = pd.read_csv(cc.REPO_ROOT / "results/pdac/v3_phase4/slides" / f"{slide}_labels.csv.gz",
                       usecols=["cell_id", "marker"]).set_index("cell_id")
    ids = adata.obs_names.astype(str)
    missing = ~ids.isin(labs.index)
    if missing.any():
        raise ValueError(f"{slide}: {int(missing.sum())} cells without v3 marker labels; no silent realignment")
    lab = labs.loc[ids, "marker"].astype(str).to_numpy()
    counts = adata.X.tocsr().copy()
    neg = adata.obsm["negprobes"]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    progs = cf.programs_used()
    lineage = json.loads((cc.RESULTS / "phase1b" / "programs_used.json").read_text())["lineage"]
    var = list(adata.var_names)
    pres = {k: [g for g in v if g in var] for k, v in progs.items()}
    lpres = {k: [g for g in v if g in var] for k, v in lineage.items()}
    s_idx = SLIDES.index(slide)
    scores, controls = p1.score_sets(adata.X.tocsr(), var, pres, 500 + s_idx)
    lin, _ = p1.score_sets(adata.X.tocsr(), var, lpres, 1500 + s_idx)
    xy = np.asarray(adata.obsm["spatial"], float)
    nn = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]
    med_nn = float(np.median(nn[nn > 0]))
    scale = med_nn / mendeley_nn_um()  # px per Mendeley-um-equivalent
    W = p1.sym_knn(xy)
    R = p1.radius_graph(xy, 30.0 * scale)
    # type-aware background (diagnostic)
    bg_rows = []
    det_all = (counts > 0)
    negdet = neg > 0
    for t, cls in cf.TYPES.items():
        m = lab == cls
        if m.sum() < 50:
            continue
        det = np.asarray(det_all[m].mean(axis=0)).ravel()
        nmax = float(negdet[m].mean(axis=0).max())
        for g in sorted({x for k, v in progs.items() if k.startswith(t + ":") for x in v} |
                        set(lineage.get(f"lineage:{t}", []))):
            dr = float(det[var.index(g)]) if g in var else np.nan
            bg_rows.append({"slide": slide, "type": t, "gene": g, "on_panel": g in var, "detection_in_type": dr,
                            "negprobe_max_in_type": nmax, "near_background": bool(dr < 2 * nmax) if g in var else True})
    pd.DataFrame(bg_rows).to_csv(sdir / f"{slide}_background_typeaware.csv", index=False)
    pd.DataFrame([{"slide": slide, "set": k, "n_genes_present": len(pres[k]), "controls": ";".join(v)}
                  for k, v in controls.items()]).to_csv(sdir / f"{slide}_controls.csv", index=False)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, cell_id=ids.to_numpy().astype("U"), xy=xy, marker=lab.astype("U"),
                        total_counts=np.asarray(counts.sum(axis=1)).ravel().astype(np.float32),
                        score_names=np.array(list(scores), dtype="U"), scores=np.column_stack(list(scores.values())),
                        lineage_names=np.array(list(lin), dtype="U"), lineage=np.column_stack(list(lin.values())),
                        knn_indptr=W.indptr, knn_indices=W.indices, r30_indptr=R.indptr, r30_indices=R.indices)
    cl, cn = np.unique(lab, return_counts=True)
    info = {"slide": slide, "unit": UNIT[slide], "n_cells": int(len(lab)), "median_nn_px": med_nn,
            "px_per_mendeley_um_equivalent": scale, "tile_px": 200.0 * scale, "radius_px": 30.0 * scale,
            "um_per_px": "UNVERIFIED", "label_counts": dict(zip(cl.tolist(), cn.tolist())),
            "program_genes_present": {k: len(v) for k, v in pres.items()},
            "program_genes_missing_from_panel": {k: sorted(set(progs[k]) - set(v)) for k, v in pres.items()}}
    (sdir / f"{slide}_info.json").write_text(json.dumps(info, indent=2))
    for f in (f"{slide}_background_typeaware.csv", f"{slide}_controls.csv", f"{slide}_info.json"):
        cc.register(sdir / f, PHASE, f"GSE310352 per-slide {f.split('_', 1)[1]}")
    log.info("prepared %s: %d cells, median NN %.2f px, tile %.1f px; %s", slide, len(lab), med_nn, 200 * scale,
             cc.mem_status())
    return info


def load(slide):
    d = cf.load_cache(slide, CACHE_SUB)
    info = json.loads((cc.RESULTS / PHASE / "slides" / f"{slide}_info.json").read_text())
    return d, info


# ----------------------------------------------------------------------------- Phase 2 repeat
def phase2_repeat(log):
    import chain_phase2_neighbourhood as p2
    progs = cf.programs_used()
    by_type = {t: [p for p in progs if p.startswith(t + ":")] for t in cf.TYPES}
    classes = sorted({c for s in SLIDES for c in load(s)[1]["label_counts"]})
    parts = {t: [] for t in cf.TYPES if by_type[t]}
    for s_i, s in enumerate(SLIDES):
        d, info = load(s)
        lab = d[ANN].astype(str)
        n = len(lab)
        W, R = d["W"], d["R"]
        S = d["scores"].astype(np.float64)
        sidx = {p: d["score_names"].index(p) for p in progs}
        deg6 = np.asarray(W.sum(axis=1)).ravel()
        onehot = np.zeros((n, len(classes)), np.float32)
        onehot[np.arange(n), [classes.index(x) for x in lab]] = 1
        comp6, compr = (W @ onehot) / np.maximum(deg6, 1)[:, None], R @ onehot
        dens = np.asarray(R.sum(axis=1)).ravel()
        tiles = cf.tile_codes(d["xy"], info["tile_px"])
        nb = {}
        for t2, cls2 in cf.TYPES.items():
            mem = lab == cls2
            cols = [sidx[p] for p in by_type[t2]]
            fill = S[mem][:, cols].mean(axis=0) if mem.any() else np.zeros(len(cols))
            for gname, A in (("nb6", W), ("nbr30", R)):
                m, _ = cf.neighbour_mean(A, S[:, cols], mem)
                m = np.where(np.isnan(m), fill[None, :], m)
                for j, p in enumerate(by_type[t2]):
                    nb[f"{gname}__{t2}__{p.split(':')[1]}"] = m[:, j]
        for t in parts:
            mem = lab == cf.TYPES[t]
            idx = np.where(mem)[0]
            if len(idx) == 0:
                continue
            rng = np.random.default_rng([cc.SEED, 600 + s_i, list(cf.TYPES).index(t)])
            if len(idx) > p2.CAP_PER_SECTION:
                idx = np.sort(rng.choice(idx, p2.CAP_PER_SECTION, replace=False))
            cols, names = [], []
            for p in by_type[t]:
                cols.append(S[idx, sidx[p]]); names.append(f"y__{p}")
                tm, has = cf.group_mean_loo(S[:, sidx[p]], tiles, mem)
                cols.append(np.where(has, tm, S[mem, sidx[p]].mean())[idx]); names.append(f"base__tilemean__{p}")
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
            parts[t].append((s, M - M.mean(axis=0, keepdims=True), names))
    rows = []
    for t, pl in parts.items():
        names = pl[0][2]
        X = np.vstack([p[1] for p in pl])
        unit = np.concatenate([[UNIT[p[0]]] * len(p[1]) for p in pl])
        ni = {n: i for i, n in enumerate(names)}
        for prog in by_type[t]:
            y = X[:, ni[f"y__{prog}"]]
            for spill in (False, True):
                base, nbr, _ = p2.feature_sets(names, prog, spill)
                for f_i, held in enumerate(sorted(set(unit))):
                    te, tr = unit == held, unit != held
                    sst = float(np.sum(y[te] ** 2))
                    for cname, cs in (("baseline", base), ("neighbour", nbr)):
                        ci = [ni[c] for c in cs]
                        pr, _ = p2._ridge_fit_predict(X[tr][:, ci], y[tr], X[te][:, ci])
                        rows.append({"program": prog, "spillover_covariates": spill, "config": cname, "model": "ridge",
                                     "tumour": held, "n": int(te.sum()), "sse": float(np.sum((y[te] - pr) ** 2)), "sst": sst})
                        ph = p2._hgb_fit_predict(X[tr][:, ci], y[tr], X[te][:, ci], [cc.SEED, 60 + f_i])
                        rows.append({"program": prog, "spillover_covariates": spill, "config": cname, "model": "hgb",
                                     "tumour": held, "n": int(te.sum()), "sse": float(np.sum((y[te] - ph) ** 2)), "sst": sst})
            log.info("phase2 repeat %s done", prog)
    df = pd.DataFrame(rows)
    cc.save_csv(df, PHASE, "phase2_repeat_per_unit_sse.csv", "GSE310352 leave-one-unit-out SSE/SST", log)
    rng = np.random.default_rng(cc.SEED)
    men = pd.read_csv(cc.RESULTS / "phase2" / "delta_r2_by_program.csv").query("annotation == 'marker'")
    out = []
    for (prog, spill, model), d in df.groupby(["program", "spillover_covariates", "model"]):
        pt, ci, per_t, r2b, r2n = p2.boot_delta(d[d.config == "baseline"], d[d.config == "neighbour"], rng)
        m = men[(men.program == prog) & (men.spillover_covariates == spill) & (men.model == model)].delta_r2.item()
        out.append({"program": prog, "spillover_covariates": spill, "model": model, "n_units": len(per_t),
                    "r2_baseline": r2b, "r2_neighbour": r2n, "delta_r2": pt, "ci_low": ci[0], "ci_high": ci[1],
                    "mendeley_marker_delta_r2": m, "same_direction_as_mendeley": bool(np.sign(pt) == np.sign(m)),
                    "ci_excludes_0": bool(ci[0] > 0 or ci[1] < 0)})
    out = pd.DataFrame(out)
    cc.save_csv(out, PHASE, "phase2_repeat_delta_r2.csv", "GSE310352 delta R2 (unit bootstrap) vs Mendeley marker-label delta R2", log)
    return out


# ----------------------------------------------------------------------------- Phase 3 repeat
def chains_and_floor(log):
    import chain_phase3_chain_test as p3
    rows = []
    for s in SLIDES:
        d, info = load(s)
        lab = d[ANN].astype(str)
        W = d["W"]
        for ch in p3.selected_chains():
            A, B, C = (cf.TYPES[x] for x in ch)
            n = [int((lab == x).sum()) for x in (A, B, C)]
            e = [ce.edge_count(W, lab, A, B), ce.edge_count(W, lab, B, C)]
            rows.append({"chain": "->".join(ch), "slide": s, "unit": UNIT[s], "n_A": n[0], "n_B": n[1], "n_C": n[2],
                         "edges_AB": e[0], "edges_BC": e[1], "passes_floor": min(n) >= 100 and min(e) >= 300})
    fl = pd.DataFrame(rows)
    cc.save_csv(fl, PHASE, "floor_by_slide_chain.csv", "GSE310352 floor per slide x chain (marker labels)", log)
    nu = fl[fl.passes_floor].groupby("chain").unit.nunique()
    return fl, [c for c in nu.index if nu[c] >= 3]


def run_slide_chains(slide, chains):
    import chain_phase3_chain_test as p3
    out = cc.out_dir(f"{PHASE}/sections") / f"{slide}.npz"
    if out.exists():
        return str(out)
    log = cc.start_run(f"chain_phase6_slide_{slide}", {"n_perm": N_PERM})
    d, info = load(slide)
    lab = d[ANN].astype(str)
    S = d["scores"].astype(np.float64)
    col = {n: j for j, n in enumerate(d["score_names"])}
    progs = list(cf.programs_used())
    ctx = ce.SectionCtx(d["xy"], d["W"], d["R"], d["total_counts"], d["lineage"], d["lineage_names"], tile=info["tile_px"])
    fl = pd.read_csv(cc.RESULTS / PHASE / "floor_by_slide_chain.csv")
    s_idx = SLIDES.index(slide)
    var = []
    for cs in chains:
        r = fl[(fl.slide == slide) & (fl.chain == cs)]
        if not bool(r.passes_floor.item()):
            continue
        ch = cs.split("->")
        var.append((cs, *(cf.TYPES[x] for x in ch), p3.triples_for(ch, progs)))
    obs, obs_sp, geoms = {}, {}, {}
    for cs, A, B, C, tr in var:
        r, g = ce.chain_stats(ctx, lab, S, col, A, B, C, tr, return_geom=True)
        geoms[cs] = g
        obs.update(p3.flatten("forward", cs, r))
        obs_sp.update(p3.flatten("forward", cs, ce.chain_stats(ctx, lab, S, col, A, B, C, tr, spill=True, geom=g)))
    keys = sorted(obs)
    ki = {k: i for i, k in enumerate(keys)}
    ne, nl = np.full((N_PERM, len(keys)), np.nan), np.full((N_PERM, len(keys)), np.nan)
    for b in range(N_PERM):
        pe = ctx.expr_perm(lab, np.random.default_rng([cc.SEED, 60, s_idx, 1, b]))
        pl = ctx.label_perm(np.random.default_rng([cc.SEED, 60, s_idx, 2, b]))
        for cs, A, B, C, tr in var:
            for k, v in p3.flatten("forward", cs, ce.chain_stats(ctx, lab, S[pe], col, A, B, C, tr, geom=geoms[cs])).items():
                ne[b, ki[k]] = v
            for k, v in p3.flatten("forward", cs, ce.chain_stats(ctx, lab[pl], S, col, A, B, C, tr)).items():
                nl[b, ki[k]] = v
        if b % 25 == 0:
            log.info("%s perm %d %s", slide, b, cc.mem_status())
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, keys=np.array(keys), obs=np.array([obs[k] for k in keys]),
                        obs_spill=np.array([obs_sp[k] for k in keys]), null_expr=ne, null_label=nl)
    tmp.replace(out)
    cc.register(out, PHASE, f"GSE310352 {slide}: observed chain stats + {N_PERM} expression/label null replicates")
    cc.finish_run(log, PHASE)
    return str(out)


def unit_aggregate(log):
    per = {}
    for s in SLIDES:
        f = cc.RESULTS / PHASE / "sections" / f"{s}.npz"
        if f.exists():
            z = np.load(f)
            per[s] = {"keys": [str(k) for k in z["keys"]], "obs": z["obs"], "sp": z["obs_spill"], "ne": z["null_expr"],
                      "nl": z["null_label"]}
    allkeys = sorted({k for p in per.values() for k in p["keys"]})
    rng = np.random.default_rng(cc.SEED)
    rows = []
    for key in allkeys:
        variant, cs, pa, pb, pc, st = key.split("|")
        if st.startswith("n_"):
            continue
        U = {}
        for s, p in per.items():
            if key in p["keys"]:
                j = p["keys"].index(key)
                if np.isfinite(p["obs"][j]):
                    U.setdefault(UNIT[s], []).append((p["obs"][j], p["sp"][j], p["ne"][:, j], p["nl"][:, j]))
        if not U:
            continue
        units = sorted(U)
        est = np.array([np.mean([x[0] for x in U[u]]) for u in units])
        esp = np.array([np.nanmean([x[1] for x in U[u]]) for u in units])
        ne = np.nanmean(np.array([np.nanmean(np.vstack([x[2] for x in U[u]]), axis=0) for u in units]), axis=0)
        nl = np.nanmean(np.array([np.nanmean(np.vstack([x[3] for x in U[u]]), axis=0) for u in units]), axis=0)
        k = rng.integers(0, len(units), size=(N_BOOT, len(units)))
        T = est.mean()
        rows.append({"chain": cs, "prog_A": pa, "prog_B": pb, "prog_C": pc, "stat": st, "n_units": len(units),
                     "units": ";".join(units), "estimate": T, "ci_low": np.quantile(est[k].mean(1), .025),
                     "ci_high": np.quantile(est[k].mean(1), .975), "p_expr": cc.add_one_p(ne, T),
                     "p_label": cc.add_one_p(nl, T), "estimate_spill": float(np.nanmean(esp))})
    tl = pd.DataFrame(rows)
    men = pd.read_csv(cc.RESULTS / "phase3" / "tumour_level_tests_author_mapped.csv").query("variant == 'forward'")
    menm_path = cc.RESULTS / "phase3" / "tumour_level_tests_marker.csv"
    menm = pd.read_csv(menm_path).query("variant == 'forward'") if menm_path.exists() else None
    key = ["chain", "prog_A", "prog_B", "prog_C", "stat"]
    tl = tl.merge(men[key + ["estimate", "ci_low", "ci_high"]].rename(columns={
        "estimate": "mendeley_author_estimate", "ci_low": "mendeley_author_ci_low", "ci_high": "mendeley_author_ci_high"}),
        on=key, how="left")
    if menm is not None:
        tl = tl.merge(menm[key + ["estimate"]].rename(columns={"estimate": "mendeley_marker_estimate"}), on=key, how="left")
        tl["same_direction_vs_mendeley_marker"] = np.sign(tl.estimate) == np.sign(tl.mendeley_marker_estimate)
    tl["same_direction_vs_mendeley_author"] = np.sign(tl.estimate) == np.sign(tl.mendeley_author_estimate)
    tl["gse_ci_excludes_0"] = (tl.ci_low > 0) | (tl.ci_high < 0)
    cc.save_csv(tl, PHASE, "phase3_repeat_unit_level.csv",
                "GSE310352 unit-level chain stats (100-perm add-one p), same-direction vs Mendeley", log)
    return tl


def figure(p2r, tl, log):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    d = p2r[~p2r.spillover_covariates]
    for model, mk in (("ridge", "o"), ("hgb", "s")):
        e = d[d.model == model]
        axes[0].errorbar(e.mendeley_marker_delta_r2, e.delta_r2, yerr=[e.delta_r2 - e.ci_low, e.ci_high - e.delta_r2],
                         fmt=mk, label=model, capsize=2)
        for r in e.itertuples():
            axes[0].annotate(r.program, (r.mendeley_marker_delta_r2, r.delta_r2), fontsize=6)
    lim = [min(-0.01, d.delta_r2.min()), max(d.delta_r2.max(), d.mendeley_marker_delta_r2.max()) + 0.01]
    axes[0].plot(lim, lim, "k:", lw=.7)
    axes[0].axhline(0, color="grey", lw=.6)
    axes[0].set_xlabel("Mendeley delta R2 (marker labels)")
    axes[0].set_ylabel("GSE310352 delta R2 (unit bootstrap 95% CI)")
    axes[0].set_title("Phase 2 repeat: neighbourhood gain", fontsize=10)
    axes[0].legend(fontsize=7)
    e = tl[tl.stat.isin(["beta_BA", "beta_CB", "indirect"])].dropna(subset=["mendeley_author_estimate"])
    for st, c in (("beta_BA", "#2ca02c"), ("beta_CB", "#1f77b4"), ("indirect", "#d62728")):
        x = e[e.stat == st]
        axes[1].errorbar(x.mendeley_author_estimate, x.estimate, yerr=[x.estimate - x.ci_low, x.ci_high - x.estimate],
                         fmt="o", ms=3, color=c, label=st, capsize=1.5, alpha=.8)
    axes[1].axhline(0, color="grey", lw=.6)
    axes[1].axvline(0, color="grey", lw=.6)
    axes[1].set_xlabel("Mendeley tumour-level estimate (author_mapped)")
    axes[1].set_ylabel("GSE310352 unit-level estimate (95% CI)")
    axes[1].set_title("Phase 3 repeat: per-hop and indirect effects (consistency check)", fontsize=10)
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, "phase6_consistency", "GSE310352 consistency of Phase 2 gains and Phase 3 chain effects", log)


def main():
    log = cc.start_run("chain_phase6_gse310352", {"slides": SLIDES, "units": UNIT, "n_perm": N_PERM,
                                                  "mendeley_median_nn_um": mendeley_nn_um()})
    infos = [prep_slide(s, log) for s in SLIDES]
    cc.save_csv(pd.DataFrame([{k: v for k, v in i.items() if not isinstance(v, dict)} for i in infos]), PHASE,
                "slide_info.csv", "GSE310352 slide sizes, NN scaling, tile sizes (um/px UNVERIFIED)", log,
                overwrite_ok=False) if not (cc.RESULTS / PHASE / "slide_info.csv").exists() else None
    bg = pd.concat([pd.read_csv(cc.RESULTS / PHASE / "slides" / f"{s}_background_typeaware.csv") for s in SLIDES])
    fr = bg.groupby(["type", "gene"]).near_background.mean().reset_index(name="frac_slides_near_background")
    fr["background_level_gse"] = fr.frac_slides_near_background >= .5
    progs = cf.programs_used()
    cons = []
    for p, genes in progs.items():
        t = p.split(":")[0]
        g = fr[(fr.type == t) & fr.gene.isin(genes)]
        nu = int((~g.background_level_gse).sum())
        cons.append({"program": p, "n_genes_mendeley": len(genes), "n_usable_gse": nu,
                     "measurement_consistent": nu >= 5 or (p == "CAF:CXCL12" and nu == 1)})
    if not (cc.RESULTS / PHASE / "program_measurement_consistency.csv").exists():
        cc.save_csv(fr, PHASE, "gene_background_typeaware_flags.csv", "GSE310352 type-aware background per gene", log)
        cc.save_csv(pd.DataFrame(cons), PHASE, "program_measurement_consistency.csv",
                    "Programs with >= 5 usable genes in GSE310352 (Mendeley lists unchanged)", log)
    p2r = phase2_repeat(log)
    fl, chains = chains_and_floor(log)
    log.info("chains analysed (floor in >= 3 units): %s", chains)
    with ProcessPoolExecutor(max_workers=min(6, int(os.environ.get("CHAIN_WORKERS", "6")))) as ex:
        futs = {ex.submit(run_slide_chains, s, chains): s for s in SLIDES}
        for f in as_completed(futs):
            log.info("slide %s done", futs[f])
    tl = unit_aggregate(log)
    figure(p2r, tl, log)
    summ = {"chains_analysed": chains, "chains_not_analysed": sorted(set(fl.chain) - set(chains)),
            "phase2_same_direction": p2r[~p2r.spillover_covariates][["program", "model", "delta_r2", "ci_low", "ci_high",
                                                                     "same_direction_as_mendeley"]].to_dict("records"),
            "phase3_same_direction_rate": {st: float(tl[tl.stat == st].same_direction_vs_mendeley_author.mean())
                                           for st in ("beta_BA", "beta_CB", "indirect")},
            "independence": "UNVERIFIED (consistency check, not replication)"}
    cc.save_json(summ, PHASE, "phase6_summary.json", "Phase 6 key numbers", log)
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

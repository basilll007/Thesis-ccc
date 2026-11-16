"""Chain Phase 4: calibration of the chain test by simulation on 3 real sections.

Plan: docs/chain_plans/phase4_calibration_simulation.md (committed before running).
Planted triple: A = tumour:basal_like, B = CAF:myCAF, C = T:cytotoxic (author_mapped labels, real geometry).
Counts of each type's planted-program genes are redrawn Poisson(f_i * c) (f_i = 1 unperturbed), per-cell
normalisation and the program score (Phase 1b recorded controls) are recomputed, then the Phase 3 section-level test
(same exposures/controls; expression block null, 100 permutations) is applied. Label block null on 20 sims for 4
settings. One task per (section, setting) -> per-task npz (resume unit); max 6 workers.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

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

PHASE = "phase4"
SECTIONS = ["T4-a", "T1-b", "U2-b"]  # largest / median (lower-middle) / smallest T-cell count among passing sections
PLANT = ("tumour:basal_like", "CAF:myCAF", "T:cytotoxic")
A_CLS, B_CLS, C_CLS = cf.TYPES["tumour"], cf.TYPES["CAF"], cf.TYPES["T"]
N_SIM, N_PERM, N_SIM_LABEL = 100, 100, 20
SETTINGS = ([("null", 1.0)] + [("chain", f) for f in (1.1, 1.25, 1.5, 2.0)] + [("zonation", f) for f in (1.5, 2.0)]
            + [("direct", f) for f in (1.1, 1.25, 1.5, 2.0)])
LABEL_SETTINGS = {("null", 1.0), ("chain", 1.5), ("zonation", 2.0), ("direct", 1.5)}
WAVELENGTH = 500.0
STATS = ["beta_CB", "beta_BA", "beta_total", "beta_direct", "indirect"]


class Planter:
    def __init__(self, section, log=None):
        import v3_nulls as vn
        d = cf.load_cache(section)
        self.d = d
        self.lab = d["author_mapped"].astype(str)
        self.col = {n: j for j, n in enumerate(d["score_names"])}
        self.ctx = ce.SectionCtx(d["xy"], d["W"], d["R"], d["total_counts"], d["lineage"], d["lineage_names"])
        progs = cf.programs_used()
        ctrl = pd.read_csv(cc.RESULTS / "phase1b" / "program_controls.csv").query("section == @section and kind == 'program'")
        adata, _, _, _, _ = vn.load_section(section, log or _NullLog())
        if not (adata.obs_names.astype(str).to_numpy() == d["cell_id"]).all():
            raise ValueError("cache/adata cell order mismatch")
        var = list(adata.var_names)
        counts = adata.layers["counts"].tocsr()
        self.L = np.asarray(counts.sum(axis=1)).ravel().astype(float)
        self.sets = {}
        for p, cls in zip(PLANT, (A_CLS, B_CLS, C_CLS)):
            genes = progs[p]
            cg = ctrl.loc[ctrl.set == p, "controls"].iloc[0].split(";")
            cells = np.where(self.lab == cls)[0]
            gi, ci = [var.index(g) for g in genes], [var.index(g) for g in cg]
            sub = counts[cells]
            self.sets[p] = {"cells": cells, "G": sub[:, gi].toarray().astype(float), "Cc": sub[:, ci].toarray().astype(float)}
            # reproduce the Phase 1b score exactly from counts (check)
            s0 = self.score(p, self.sets[p]["G"])
            ref = d["scores"][cells, self.col[p]]
            if not np.allclose(s0, ref, atol=1e-4):
                raise ValueError(f"rescoring mismatch for {p}: max diff {np.abs(s0 - ref).max()}")
        del adata, counts
        W = self.ctx.W
        a = self.sets[PLANT[0]]["cells"]
        sa = d["scores"][a, self.col[PLANT[0]]]
        self.a_high = np.zeros(len(self.lab), bool)
        self.a_high[a[sa >= np.quantile(sa, 0.75)]] = True
        self.nb_ahigh = np.asarray(W @ self.a_high.astype(np.float32)).ravel() > 0
        self.base_S = d["scores"].astype(np.float64)
        _, self.geom = ce.chain_stats(self.ctx, self.lab, self.base_S, self.col, A_CLS, B_CLS, C_CLS, [PLANT],
                                      return_geom=True)

    def score(self, p, Gnew):
        st = self.sets[p]
        Lnew = self.L[st["cells"]] - st["G"].sum(axis=1) + Gnew.sum(axis=1)
        f = 1e4 / np.maximum(Lnew, 1)[:, None]
        return np.log1p(Gnew * f).mean(axis=1) - np.log1p(st["Cc"] * f).mean(axis=1)

    def simulate(self, scenario, fold, rng):
        """Return a perturbed score matrix S and the per-type perturbed-cell masks."""
        n = len(self.lab)
        pert = {p: np.zeros(n, bool) for p in PLANT}
        if scenario == "chain":
            pert[PLANT[1]] = (self.lab == B_CLS) & self.nb_ahigh
            pert[PLANT[2]] = (self.lab == C_CLS) & (np.asarray(self.ctx.W @ pert[PLANT[1]].astype(np.float32)).ravel() > 0)
        elif scenario == "zonation":
            ph = rng.uniform(0, 2 * np.pi, 2)
            xy = self.d["xy"]
            z = np.sin(2 * np.pi * xy[:, 0] / WAVELENGTH + ph[0]) + np.sin(2 * np.pi * xy[:, 1] / WAVELENGTH + ph[1])
            for p, cls in zip(PLANT, (A_CLS, B_CLS, C_CLS)):
                pert[p] = (self.lab == cls) & (z > 0)
        elif scenario == "direct":
            pert[PLANT[2]] = (self.lab == C_CLS) & self.nb_ahigh
        S = self.base_S.copy()
        for p in PLANT:
            st = self.sets[p]
            fvec = np.where(pert[p][st["cells"]], fold, 1.0)
            Gnew = rng.poisson(st["G"] * fvec[:, None]).astype(float)
            S[st["cells"], self.col[p]] = self.score(p, Gnew)
        return S, {p: int(pert[p].sum()) for p in PLANT}


class _NullLog:
    def info(self, *a, **k):
        pass

    warning = info


def stats_vec(r):
    return np.array([r[PLANT][s] for s in STATS] + [r[PLANT]["n_chain"], r[PLANT]["n_hop1"], r[PLANT]["n_hop2"]])


def run_task(section, s_idx, scenario, fold):
    out = cc.out_dir(f"{PHASE}/tasks") / f"{section}__{scenario}__{fold}.npz"
    if out.exists():
        return str(out)
    log = cc.start_run(f"chain_phase4_{section}_{scenario}_{fold}", {"n_sim": N_SIM, "n_perm": N_PERM})
    pl = Planter(section, log)
    ctx, col, lab = pl.ctx, pl.col, pl.lab
    obs = np.full((N_SIM, len(STATS) + 3), np.nan)
    null = np.full((N_SIM, N_PERM, len(STATS)), np.nan)
    do_label = (scenario, fold) in LABEL_SETTINGS
    lnull = np.full((N_SIM_LABEL, N_PERM, len(STATS)), np.nan)
    npert = []
    sc_i = ["null", "chain", "zonation", "direct"].index(scenario)
    for sim in range(N_SIM):
        rng = np.random.default_rng([cc.SEED, 4, s_idx, sc_i, int(round(fold * 100)), sim])
        S, npt = pl.simulate(scenario, fold, rng)
        npert.append(npt)
        obs[sim] = stats_vec(ce.chain_stats(ctx, lab, S, col, A_CLS, B_CLS, C_CLS, [PLANT], geom=pl.geom))
        for b in range(N_PERM):
            pe = ctx.expr_perm(lab, np.random.default_rng([cc.SEED, 41, s_idx, sc_i, sim, b]))
            null[sim, b] = stats_vec(ce.chain_stats(ctx, lab, S[pe], col, A_CLS, B_CLS, C_CLS, [PLANT],
                                                    geom=pl.geom))[:len(STATS)]
        if do_label and sim < N_SIM_LABEL:
            for b in range(N_PERM):
                p_l = ctx.label_perm(np.random.default_rng([cc.SEED, 42, s_idx, sc_i, sim, b]))
                lnull[sim, b] = stats_vec(ce.chain_stats(ctx, lab[p_l], S, col, A_CLS, B_CLS, C_CLS,
                                                         [PLANT]))[:len(STATS)]
        if sim % 20 == 0:
            log.info("%s %s %.2f sim %d %s", section, scenario, fold, sim, cc.mem_status())
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, obs=obs, null=null, lnull=lnull, do_label=do_label,
                        n_perturbed=np.array([[d_[p] for p in PLANT] for d_ in npert]))
    tmp.replace(out)
    cc.register(out, PHASE, f"Simulation task {section} {scenario} f={fold}: observed stats + expression-null replicates")
    cc.finish_run(log, PHASE)
    return str(out)


def p2(null, obs):
    """Add-one two-sided p on |stat| per simulation: null (sims x perms), obs (sims,)."""
    k = (np.abs(null) >= np.abs(obs)[:, None]).sum(axis=1)
    return (k + 1) / (np.isfinite(null).sum(axis=1) + 1)


def aggregate(log):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows, sims = [], {}
    for s_idx, sec in enumerate(SECTIONS):
        for scen, f in SETTINGS:
            z = np.load(cc.RESULTS / PHASE / "tasks" / f"{sec}__{scen}__{f}.npz")
            obs, null = z["obs"], z["null"]
            sims[(sec, scen, f)] = (obs, null)
            pv = {s: p2(null[:, :, j], obs[:, j]) for j, s in enumerate(STATS)}
            shrink = np.abs(obs[:, 3]) < np.abs(obs[:, 2])
            call = (pv["beta_CB"] < .05) & (pv["beta_BA"] < .05) & (pv["indirect"] < .05) & shrink
            r = {"section": sec, "scenario": scen, "fold": f, "n_sim": len(obs),
                 "median_n_chain": float(np.nanmedian(obs[:, 5])), "median_n_hop1": float(np.nanmedian(obs[:, 6])),
                 "median_n_hop2": float(np.nanmedian(obs[:, 7])),
                 "median_n_perturbed_B": float(np.median(z["n_perturbed"][:, 1])),
                 "median_n_perturbed_C": float(np.median(z["n_perturbed"][:, 2])),
                 "chain_call_rate": float(call.mean()), "hop_CB_rate": float((pv["beta_CB"] < .05).mean()),
                 "hop_BA_rate": float((pv["beta_BA"] < .05).mean()), "indirect_rate": float((pv["indirect"] < .05).mean()),
                 **{f"mean_{s}": float(np.nanmean(obs[:, j])) for j, s in enumerate(STATS)},
                 "mean_abs_indirect": float(np.nanmean(np.abs(obs[:, 4])))}
            if bool(z["do_label"]):
                ln = z["lnull"]
                k = ln.shape[0]
                pl = {s: p2(ln[:, :, j], obs[:k, j]) for j, s in enumerate(STATS)}
                r.update({"label_null_n_sim": k, "label_hop_CB_rate": float((pl["beta_CB"] < .05).mean()),
                          "label_hop_BA_rate": float((pl["beta_BA"] < .05).mean()),
                          "label_indirect_rate": float((pl["indirect"] < .05).mean()),
                          "chain_call_rate_with_label_criterion": float(
                              (call[:k] & (pl["beta_CB"] < .05) & (pl["beta_BA"] < .05)).mean())})
            rows.append(r)
        log.info("aggregated %s", sec)
    tab = pd.DataFrame(rows)
    # pooled 3-section call (one section per tumour; same sim index)
    prow = []
    for scen, f in SETTINGS:
        o = np.stack([sims[(s, scen, f)][0][:, :len(STATS)] for s in SECTIONS])   # 3 x sims x stats
        nl = np.stack([sims[(s, scen, f)][1] for s in SECTIONS])                  # 3 x sims x perms x stats
        T, TN = np.nanmean(o, axis=0), np.nanmean(nl, axis=0)
        pv = {s: p2(TN[:, :, j], T[:, j]) for j, s in enumerate(STATS)}
        call = (pv["beta_CB"] < .05) & (pv["beta_BA"] < .05) & (pv["indirect"] < .05) & (np.abs(T[:, 3]) < np.abs(T[:, 2]))
        prow.append({"scenario": scen, "fold": f, "pooled_sections": ";".join(SECTIONS), "chain_call_rate": float(call.mean()),
                     "hop_CB_rate": float((pv["beta_CB"] < .05).mean()), "hop_BA_rate": float((pv["beta_BA"] < .05).mean()),
                     "indirect_rate": float((pv["indirect"] < .05).mean()), "mean_indirect": float(np.nanmean(T[:, 4]))})
    pooled = pd.DataFrame(prow)
    cc.save_csv(tab, PHASE, "calibration_table.csv",
                "Per section x scenario x fold: chain-call rate (power / FPR), per-hop rates, mean effects (add-one p, 100 perms)", log)
    cc.save_csv(pooled, PHASE, "calibration_pooled_3sections.csv", "Pooled 3-section (tumour-level style) call rates", log)
    # calibrated range
    null_ok = {s: tab[(tab.section == s) & (tab.scenario == "null")].chain_call_rate.item() <= .10 for s in SECTIONS}
    zon_ok = {s: tab[(tab.section == s) & (tab.scenario == "zonation") & (tab.fold == 2.0)].chain_call_rate.item() <= .10
              for s in SECTIONS}
    trust = []
    for r in tab[tab.scenario == "chain"].itertuples():
        dr = tab[(tab.section == r.section) & (tab.scenario == "direct") & (tab.fold == r.fold)].chain_call_rate.item()
        ok = r.chain_call_rate >= .8 and null_ok[r.section] and zon_ok[r.section] and dr <= .10
        trust.append({"section": r.section, "fold": r.fold, "power": r.chain_call_rate, "direct_call_rate": dr,
                      "null_ok": null_ok[r.section], "zonation_ok": zon_ok[r.section], "trustworthy": ok,
                      "median_n_chain": r.median_n_chain, "mean_abs_indirect": r.mean_abs_indirect})
    trust = pd.DataFrame(trust)
    cc.save_csv(trust, PHASE, "trustworthy_settings.csv", "Decision-rule check per section x fold (plan: calibrated range)", log)
    tt = trust[trust.trustworthy]
    if len(tt):
        best = tt.sort_values(["fold", "median_n_chain"]).iloc[0]
        cal = {"trustworthy_setting_found": True, "min_fold_trustworthy": float(tt.fold.min()),
               "min_n_chain_trustworthy": float(tt.median_n_chain.min()),
               "min_abs_indirect_trustworthy": float(tt[tt.fold == tt.fold.min()].mean_abs_indirect.min()),
               "smallest_trustworthy_setting": best.to_dict(),
               "note": "Phase 3 calls are interpreted only within this range; calibration used the primary chain triple"}
    else:
        cal = {"trustworthy_setting_found": False,
               "note": "No (fold, section) met power >= 0.8 with null/zonation/direct call rates <= 0.10; every Phase 3 "
                       "chain call is outside the calibrated range"}
    cal.update({"null_call_rate": {s: tab[(tab.section == s) & (tab.scenario == "null")].chain_call_rate.item() for s in SECTIONS},
                "zonation2_call_rate": {s: tab[(tab.section == s) & (tab.scenario == "zonation") & (tab.fold == 2.0)]
                                        .chain_call_rate.item() for s in SECTIONS}})
    cc.save_json(cal, PHASE, "calibrated_range.json", "Calibrated range for interpreting Phase 3 (decision rule in plan)", log)
    # figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, scen in zip(axes, ["chain", "zonation", "direct"]):
        for sec, mk in zip(SECTIONS, "osd"):
            d = tab[(tab.section == sec) & (tab.scenario == scen)].sort_values("fold")
            ax.plot(d.fold, d.chain_call_rate, marker=mk, label=f"{sec} (n_chain~{d.median_n_chain.median():.0f})")
            if scen == "chain":
                ax.plot(d.fold, d.hop_CB_rate, marker=mk, ls=":", alpha=.6, color=ax.lines[-1].get_color())
        pdd = pooled[pooled.scenario == scen].sort_values("fold")
        ax.plot(pdd.fold, pdd.chain_call_rate, "k--", marker="x", label="pooled 3 sections")
        ax.axhline(.8 if scen == "chain" else .10, color="grey", lw=.7)
        ax.set_title({"chain": "(i) planted chain: power", "zonation": "(ii) zonation only: false calls",
                      "direct": "(iii) direct A->C: wrongly called chain"}[scen], fontsize=10)
        ax.set_xlabel("per-hop fold change f")
    axes[0].set_ylabel("chain-call rate (100 sims; dotted = hop C<-B rate)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, "phase4_calibration", "Chain-test calibration: power, zonation FPR, direct-effect misclassification", log)
    return tab, pooled, cal


def main():
    if "--aggregate-only" not in sys.argv:
        tasks = [(s, i, scen, f) for i, s in enumerate(SECTIONS) for scen, f in SETTINGS]
        if os.environ.get("CHAIN_REVERSE"):  # second scheduler instance takes tasks from the end (6-worker cap)
            tasks = [t for t in tasks[::-1] if not (cc.RESULTS / PHASE / "tasks" / f"{t[0]}__{t[2]}__{t[3]}.npz").exists()]
        with ProcessPoolExecutor(max_workers=min(6, int(os.environ.get("CHAIN_WORKERS", "6")))) as ex:
            futs = {ex.submit(run_task, *t): t for t in tasks}
            for f in as_completed(futs):
                print("done", futs[f], f.result(), flush=True)
        if os.environ.get("CHAIN_REVERSE"):
            return
    log = cc.start_run("chain_phase4_aggregate", {"sections": SECTIONS, "plant": PLANT, "settings": SETTINGS,
                                                  "n_sim": N_SIM, "n_perm": N_PERM, "n_sim_label": N_SIM_LABEL})
    log.info("Logged reductions: 100 expression-null permutations per simulation (Phase 3 uses 200); label null "
             "on 20 sims for 4 settings")
    tab, pooled, cal = aggregate(log)
    log.info("calibration table:\n%s", tab[["section", "scenario", "fold", "median_n_chain", "chain_call_rate",
                                             "hop_CB_rate", "hop_BA_rate", "indirect_rate"]].to_string(index=False))
    log.info("pooled:\n%s", pooled.to_string(index=False))
    log.info("calibrated range: %s", json.dumps(cal, default=str))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

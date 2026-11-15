"""v3 Phase 2: calibration by simulation on real tissue (plan: docs/v3_plans/phase2_calibration_simulation.md).

Three real sections (smallest / lower-median / largest cell count), real positions, real 6-NN graph, real consensus_T labels.
Two planted ground truths on the CXCL12->CXCR4 CAF->T_cell axis (perturbations of raw counts, Poisson(f*c), then renormalise):
  PROXIMITY : CXCR4 of T cells adjacent (6-NN col) to CXCL12+ CAF senders scaled by f.
  ZONATION  : CXCL12 and CXCR4 of ALL cells in random "high" 200 um tiles scaled by f (no adjacency, no cell-type dependence).
Scores C, M, E, Q under nulls global / block200 / block100 / block50 / dist (500 perms, seed 42; the permutation sets are
built once per section x null and re-used; evaluation by bincount is algebraically identical to ccc_audit.run_null for one row).
100 simulations per setting; add-one p everywhere. The planted perturbation is a computational stand-in for a stimulation
experiment; no experimental claim is made.
Usage: python src/pdac/v3_phase2_calibration.py [--sections ...] [--n-sims 100] [--aggregate-only]
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import SEED, finish_run, mem_status, out_dir, register, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import ccc_audit as ca  # noqa: E402
import v3_nulls as vn  # noqa: E402

PHASE = "v3_phase2"
N_PERMS = 500
SECTIONS = {"T6-b": "smallest", "T2-b": "lower_median", "U5-a": "largest"}
LIG, REC, SENDER, RECEIVER = "CXCL12", "CXCR4", "CAF_fibroblast", "T_cell"
EFFECTS = [1.0, 1.1, 1.25, 1.5, 2.0]
ZONE_TILE = 200.0
SCORES = ["comp", "matched", "edge", "ratio"]
SCORE_LAB = {"comp": "C", "matched": "M", "edge": "E", "ratio": "Q"}
P_NOTE = "add-one (k+1)/(n+1), one-sided"


class PermStore:
    """Per-permutation index sets for one (sender, receiver) row, concatenated with permutation ids."""

    def __init__(self, perms, rows, cols, T, iS, iR):
        n = len(perms)
        e_idx, e_pid, s_idx, s_pid, r_idx, r_pid = [], [], [], [], [], []
        self.e_count = np.zeros(n, np.int64)
        self.s_count = np.zeros(n, np.int64)
        self.r_count = np.zeros(n, np.int64)
        for k, pl in enumerate(perms):
            e = np.flatnonzero((pl[rows] == iS) & (pl[cols] == iR)).astype(np.int32)
            s = np.flatnonzero(pl == iS).astype(np.int32)
            r = np.flatnonzero(pl == iR).astype(np.int32)
            self.e_count[k], self.s_count[k], self.r_count[k] = len(e), len(s), len(r)
            e_idx.append(e); e_pid.append(np.full(len(e), k, np.int32))
            s_idx.append(s); s_pid.append(np.full(len(s), k, np.int32))
            r_idx.append(r); r_pid.append(np.full(len(r), k, np.int32))
        self.n = n
        self.e_idx, self.e_pid = np.concatenate(e_idx), np.concatenate(e_pid)
        self.s_idx, self.s_pid = np.concatenate(s_idx), np.concatenate(s_pid)
        self.r_idx, self.r_pid = np.concatenate(r_idx), np.concatenate(r_pid)

    def null_scores(self, L, R, prods):
        """(n_perms, 4) array of C, M, E, Q under the stored permutations (same conventions as ccc_audit: E=0 if no edges)."""
        with np.errstate(invalid="ignore", divide="ignore"):
            Ls = np.bincount(self.s_pid, weights=L[self.s_idx], minlength=self.n) / self.s_count
            Rr = np.bincount(self.r_pid, weights=R[self.r_idx], minlength=self.n) / self.r_count
            Ls[self.s_count == 0] = 0.0
            Rr[self.r_count == 0] = 0.0
            E = np.bincount(self.e_pid, weights=prods[self.e_idx], minlength=self.n) / self.e_count
            E[self.e_count == 0] = 0.0
            M = Ls * Rr
            Q = np.where(M > 0, E / M, np.nan)
        return np.column_stack([(Ls + Rr) / 2.0, M, E, Q])

    def nbytes(self):
        return sum(a.nbytes for a in (self.e_idx, self.e_pid, self.s_idx, self.s_pid, self.r_idx, self.r_pid))


def observed_scores(L, R, prods, e_obs, s_obs, r_obs):
    Ls, Rr = L[s_obs].mean(), R[r_obs].mean()
    E = prods[e_obs].mean() if len(e_obs) else np.nan
    M = Ls * Rr
    return np.array([(Ls + Rr) / 2.0, M, E, E / M if M > 0 else np.nan])


def lognorm(cL, cR, total):
    with np.errstate(invalid="ignore", divide="ignore"):
        L = np.log1p(cL / total * 1e4)
        R = np.log1p(cR / total * 1e4)
    return np.nan_to_num(L), np.nan_to_num(R)


def run_section(section, n_sims, log, od):
    outp = od / f"sim_rows_{section}.csv.gz"
    if outp.exists():
        log.info("section %s: %s exists (resume)", section, outp.name)
        return
    t0 = time.time()
    adata, meta, labs, pairs, det = vn.load_section(section, log)
    xy = np.asarray(adata.obsm["spatial"], float)
    rows, cols = vn.knn_edges(xy, 6)
    labels_str = labs["consensus_T"].astype(str).to_numpy()
    types, lab = vn.encode(labels_str)
    T, iS, iR = len(types), types.index(SENDER), types.index(RECEIVER)
    vnames = list(adata.var_names)
    cL = adata.layers["counts"][:, vnames.index(LIG)].toarray().ravel().astype(float)
    cR = adata.layers["counts"][:, vnames.index(REC)].toarray().ravel().astype(float)
    total = np.asarray(adata.layers["counts"].sum(axis=1)).ravel().astype(float)
    # check our renormalisation reproduces adata.raw (normalize_total 1e4 + log1p)
    L0, R0 = lognorm(cL, cR, total)
    rawL = ca.gene_matrix(adata, [LIG, REC])
    err = max(np.abs(L0 - rawL[:, 0]).max(), np.abs(R0 - rawL[:, 1]).max())
    log.info("section %s: renormalisation check max|diff| vs adata.raw = %.2e", section, err)
    if err > 1e-5:
        raise RuntimeError("renormalisation mismatch")
    e_obs = np.flatnonzero((lab[rows] == iS) & (lab[cols] == iR))
    s_obs, r_obs = np.flatnonzero(lab == iS), np.flatnonzero(lab == iR)
    targets = np.unique(cols[e_obs[cL[rows[e_obs]] > 0]])  # T cells adjacent to CXCL12+ CAF senders
    log.info("section %s: %d cells, %d edges, CAF %d, T %d, CAF->T edges %d, proximity targets %d (CXCR4 counts in targets: mean %.3f, frac>0 %.3f)",
             section, adata.n_obs, len(rows), len(s_obs), len(r_obs), len(e_obs), len(targets), cR[targets].mean(), (cR[targets] > 0).mean())
    del adata
    gc.collect()
    stores, null_info = {}, {}
    for n in vn.NULL_ORDER:
        t1 = time.time()
        perms, info = vn.make_perms(n, labels_str, types, lab, xy, N_PERMS)
        if perms is None:
            log.warning("section %s null %s skipped: %s", section, n, info)
            null_info[n] = info
            continue
        stores[n] = PermStore(perms, rows, cols, T, iS, iR)
        del perms
        gc.collect()
        null_info[n] = {**info, "store_MB": round(stores[n].nbytes() / 1e6, 1), "median_null_edges": float(np.median(stores[n].e_count))}
        log.info("section %s null %s: store built in %.1fs (%s); %s", section, n, time.time() - t1, null_info[n], mem_status())
    save_json({"section": section, "role": SECTIONS[section], "n_cells": int(len(lab)), "n_edges": int(len(rows)), "n_CAF": int(len(s_obs)),
               "n_T": int(len(r_obs)), "n_CAF_T_edges": int(len(e_obs)), "n_proximity_targets": int(len(targets)),
               "mean_CXCR4_count_targets": float(cR[targets].mean()), "null_info": null_info, "n_perms": N_PERMS, "n_sims": n_sims},
              f"{PHASE}/sections", f"{section}_setup.json", f"Phase 2 setup for section {section}: counts, targets, null stores", log)

    recs = []
    settings = [("proximity", f) for f in EFFECTS] + [("zonation", f) for f in EFFECTS if f != 1.0]
    for scen, f in settings:
        t1 = time.time()
        for s in range(n_sims):
            rng = np.random.default_rng(SEED + s)
            nL, nR = cL.copy(), cR.copy()
            if scen == "proximity":
                nR[targets] = rng.poisson(f * cR[targets])
                n_pert = len(targets)
            else:
                off = rng.uniform(0, ZONE_TILE, size=2)
                t = np.floor((xy - xy.min(axis=0) + off) / ZONE_TILE).astype(np.int64)
                code = t[:, 0] * (t[:, 1].max() + 1) + t[:, 1]
                uniq, inv = np.unique(code, return_inverse=True)
                high = rng.random(len(uniq)) < 0.5
                cells = np.flatnonzero(high[inv])
                nL[cells] = rng.poisson(f * cL[cells])
                nR[cells] = rng.poisson(f * cR[cells])
                n_pert = len(cells)
            tot = total + (nL - cL) + (nR - cR)
            L, R = lognorm(nL, nR, tot)
            prods = L[rows] * R[cols]
            obs = observed_scores(L, R, prods, e_obs, s_obs, r_obs)
            for n, st in stores.items():
                null = st.null_scores(L, R, prods)
                stt = ca.stats(obs, null)
                for k, sc in enumerate(SCORES):
                    recs.append({"section": section, "scenario": scen, "effect": f, "sim": s, "null": n, "score": SCORE_LAB[sc],
                                 "obs": obs[k], "null_mean": stt["null_mean"][k], "null_sd": stt["null_sd"][k], "z": stt["z"][k],
                                 "p_plus1": stt["p_plus1"][k], "n_null_valid": int(stt["n_null_valid"][k]),
                                 "n_perturbed_cells": n_pert, "n_edges_obs": int(len(e_obs))})
        sub = pd.DataFrame(recs[-n_sims * len(stores) * len(SCORES):])
        q = sub[(sub.score == "Q")].groupby("null").agg(power=("p_plus1", lambda p: float((p < 0.05).mean())), mean_z=("z", "mean"))
        log.info("section %s %s f=%.2f: %d sims in %.1fs; Q p<.05 frac by null: %s; mean zQ: %s; %s", section, scen, f, n_sims,
                 time.time() - t1, q.power.round(2).to_dict(), q.mean_z.round(2).to_dict(), mem_status())
        pd.DataFrame(recs).to_csv(od / f"sim_rows_{section}.partial.csv.gz", index=False)  # partial save after every setting
    df = pd.DataFrame(recs)
    df["p_value_convention"] = P_NOTE
    df.to_csv(outp, index=False)
    (od / f"sim_rows_{section}.partial.csv.gz").unlink(missing_ok=True)
    register(outp, PHASE, f"Phase 2 simulation rows for section {section}: scenario x effect x sim x null x score")
    del stores
    gc.collect()
    log.info("section %s done in %.1fs; %s", section, time.time() - t0, mem_status())


def aggregate(log, od):
    df = pd.concat([pd.read_csv(f) for f in sorted(od.glob("sim_rows_*.csv.gz")) if ".partial" not in f.name], ignore_index=True)
    p = out_dir(PHASE) / "sim_rows.csv.gz"
    df.to_csv(p, index=False)
    register(p, PHASE, "All Phase 2 simulation rows (3 sections): obs, null mean/sd, z, add-one p per score x null x scenario x effect x sim")
    df["detect"] = df.p_plus1 < 0.05
    keys = ["scenario", "effect", "null", "score"]
    per = df.groupby(["section"] + keys).agg(n_sims=("sim", "count"), rate_p05=("detect", "mean"), mean_z=("z", "mean"),
                                             sd_z=("z", "std"), mean_obs=("obs", "mean"), mean_null_mean=("null_mean", "mean"),
                                             n_edges_obs=("n_edges_obs", "first"), n_perturbed_cells=("n_perturbed_cells", "mean")).reset_index()
    pooled = df.groupby(keys).agg(n_sims=("sim", "count"), rate_p05=("detect", "mean"), mean_z=("z", "mean"), sd_z=("z", "std"),
                                  mean_obs=("obs", "mean")).reset_index()
    pooled.insert(0, "section", "pooled_3_sections")
    tab = pd.concat([per, pooled], ignore_index=True)
    tab["rate_meaning"] = np.where(tab.scenario == "proximity", np.where(tab.effect == 1.0, "control_rate", "power"), "false_positive_rate")
    tab["null"] = pd.Categorical(tab["null"], vn.NULL_ORDER, ordered=True)
    tab["score"] = pd.Categorical(tab["score"], ["C", "M", "E", "Q"], ordered=True)
    tab = tab.sort_values(["section", "scenario", "effect", "score", "null"]).reset_index(drop=True)
    tab["p_value_convention"] = P_NOTE
    save_csv(tab, PHASE, "calibration_table.csv",
             "Power (proximity) / false-positive rate (zonation) per section and pooled, per scenario x effect x null x score; 100 sims each", log)
    # compact deliverable: which score/null detects which effect size; which is fooled by zonation
    pp = tab[tab.section == "pooled_3_sections"]
    wide = pp.pivot_table(index=["score", "null"], columns=["scenario", "effect"], values="rate_p05", observed=True)
    wide.columns = [f"{a}_f{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    prox_cols = [c for c in wide.columns if c.startswith("proximity") and not c.endswith("f1.0")]

    def min_detectable(r):
        for c in prox_cols:
            if r[c] >= 0.8:
                return float(c.split("_f")[1])
        return np.nan

    wide["min_effect_power_ge_0.8"] = wide.apply(min_detectable, axis=1)
    wide["fooled_by_zonation_f1.5"] = wide["zonation_f1.5"] > 0.10
    wide["control_rate"] = wide["proximity_f1.0"]
    save_csv(wide, PHASE, "calibration_deliverable.csv",
             "Pooled over 3 sections: detection rate per score x null for each planted proximity effect and zonation effect; min detectable effect (power>=0.8); fooled-by-zonation flag (FPR>0.10 at f=1.5)", log)
    W = wide.to_dict(orient="records")
    for r in W:
        log.info("CALIB %s %s: control %.2f; power f1.1 %.2f f1.25 %.2f f1.5 %.2f f2 %.2f; zonation FPR f1.1 %.2f f1.25 %.2f f1.5 %.2f f2 %.2f; min f(power>=.8) %s; fooled %s",
                 r["score"], r["null"], r["control_rate"], r["proximity_f1.1"], r["proximity_f1.25"], r["proximity_f1.5"], r["proximity_f2.0"],
                 r["zonation_f1.1"], r["zonation_f1.25"], r["zonation_f1.5"], r["zonation_f2.0"], r["min_effect_power_ge_0.8"], r["fooled_by_zonation_f1.5"])
    summary = {"n_sims_per_setting": int(df.groupby(["section", "scenario", "effect", "null", "score"]).sim.count().min()),
               "n_perms": N_PERMS, "sections": SECTIONS, "effects": EFFECTS, "nulls": vn.NULL_ORDER,
               "min_effect_power_ge_0.8": {f"{r['score']}/{r['null']}": (None if np.isnan(r["min_effect_power_ge_0.8"]) else r["min_effect_power_ge_0.8"]) for r in W},
               "fooled_by_zonation_f1.5": {f"{r['score']}/{r['null']}": bool(r["fooled_by_zonation_f1.5"]) for r in W},
               "note": "planted perturbation is a computational stand-in for a stimulation experiment; no experimental claim"}
    save_json(summary, PHASE, "calibration_summary.json", "Phase 2 calibration summary (min detectable effect per score/null; zonation false positives)", log)

    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharey=True)
    colors = {"global": "#000000", "block200": "#1f77b4", "block100": "#2ca02c", "block50": "#d62728", "dist": "#9467bd"}
    for j, sc in enumerate(["C", "M", "E", "Q"]):
        for i, scen in enumerate(["proximity", "zonation"]):
            ax = axes[i, j]
            for n in vn.NULL_ORDER:
                g = pp[(pp.scenario == scen) & (pp.score == sc) & (pp["null"] == n)].sort_values("effect")
                if scen == "zonation":
                    ctrl = pp[(pp.scenario == "proximity") & (pp.effect == 1.0) & (pp.score == sc) & (pp["null"] == n)]
                    g = pd.concat([ctrl.assign(scenario="zonation"), g])
                ax.plot(g.effect, g.rate_p05, "o-", color=colors[n], label=n)
            ax.axhline(0.05, color="grey", lw=0.6, ls="--")
            ax.axhline(0.8 if scen == "proximity" else 0.10, color="grey", lw=0.6, ls=":")
            ax.set_xticks(EFFECTS)
            ax.set_title(f"{sc}: {'power, planted proximity' if scen == 'proximity' else 'false-positive rate, zonation only'}", fontsize=9)
            ax.set_xlabel("planted fold change f (f=1 is the control)")
            if j == 0:
                ax.set_ylabel("fraction of sims with add-one p < 0.05")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("v3 Phase 2 calibration on real tissue (T6-b, T2-b, U5-a pooled; 100 sims x 500 perms; CXCL12→CXCR4 CAF→T; audit stand-in, not an experiment)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase2_calibration", "Power vs planted proximity effect and false-positive rate vs zonation effect, per score and null", log)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--n-sims", type=int, default=100)
    ap.add_argument("--aggregate-only", action="store_true")
    ap.add_argument("--no-aggregate", action="store_true")
    args = ap.parse_args()
    log = start_run("v3_phase2_calibration")
    log.info("sections=%s effects=%s nulls=%s n_perms=%d n_sims=%d; %s", SECTIONS, EFFECTS, vn.NULL_ORDER, N_PERMS, args.n_sims, mem_status())
    od = out_dir(f"{PHASE}/sections")
    if not args.aggregate_only:
        for s in SECTIONS:
            if args.sections and s not in args.sections:
                continue
            run_section(s, args.n_sims, log, od)
    if not args.no_aggregate:
        aggregate(log, od)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

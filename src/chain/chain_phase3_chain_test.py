"""Chain Phase 3: A -> B -> C chain test with block nulls, negative controls and tumour-level inference.

Plan: docs/chain_plans/phase3_chain_test.md (committed before running).
Per section (one worker per section, max 6 workers): observed statistics (with and without spillover covariates),
200 expression-block and 200 label-block null replicates (200 um tiles, random offsets), saved per section
(resume unit). Then tumour-level aggregation (paired null replicates across sections), BH, decision rule.
Usage: python chain_phase3_chain_test.py [--annotation author_mapped|marker] [--aggregate-only]
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

PHASE = "phase3"
N_PERM = 200
N_BOOT = 2000
FLOOR_CELLS, FLOOR_EDGES = 100, 300
CHAIN_SEL = cc.RESULTS / "phase1b" / "chain_selection.csv"
VARIANTS = ["forward", "reversed", "randomB", "shuffledB"]
STAT_TESTED = ["beta_CB", "beta_BA", "indirect", "beta_total", "beta_direct"]


def selected_chains():
    sel = pd.read_csv(CHAIN_SEL)
    sel = sel[sel.role.isin(["primary", "secondary"])]
    order = ["primary"] + ["secondary"] * 5
    sel = pd.concat([sel[sel.role == "primary"], sel[sel.role == "secondary"]])
    return [tuple(c.split("->")) for c in sel.chain]


def random_b_choice(ann):
    counts = pd.read_csv(cc.RESULTS / "phase1b" / "cell_counts_by_section.csv").query("annotation == @ann")
    ok = counts[counts.n_cells >= 100].groupby("cell_type").tumour.nunique()
    cands_all = sorted(ok.index[ok >= 3])
    rng = np.random.default_rng(cc.SEED)
    out = {}
    for ch in selected_chains():
        excl = {cf.TYPES[x] for x in ch}
        cands = [c for c in cands_all if c not in excl and c not in ("unassigned",)]
        out["->".join(ch)] = str(rng.choice(cands))
    return out, cands_all


def triples_for(ch, progs):
    by = {t: [p for p in progs if p.startswith(t + ":")] for t in cf.TYPES}
    A, B, C = ch
    return [(a, b, c) for a in by[A] for b in by[B] for c in by[C]]


def floor_pass(ann):
    fl = pd.read_csv(cc.RESULTS / "phase1b" / "floor_by_section_chain.csv").query("annotation == @ann")
    return {(r.chain, r.section): bool(r.passes_floor) for r in fl.itertuples()}


def build_variants(section, s_idx, ann, d, lab, S, col, progs, rb, fp, ctx):
    """List of (variant, chain_str, A_cls, B_cls, C_cls, triples, S_base) for chains passing the floor here."""
    out = []
    for c_i, ch in enumerate(selected_chains()):
        cs = "->".join(ch)
        if not fp.get((cs, section), False):
            continue
        A, B, C = (cf.TYPES[x] for x in ch)
        tr = triples_for(ch, progs)
        out.append(("forward", cs, A, B, C, tr, S))
        out.append(("reversed", cs, C, B, A, [(c, b, a) for a, b, c in tr], S))
        Bp = rb[cs]
        nB = int((lab == Bp).sum())
        if nB >= FLOOR_CELLS and min(ce.edge_count(ctx.W, lab, A, Bp), ce.edge_count(ctx.W, lab, Bp, C)) >= FLOOR_EDGES:
            out.append(("randomB", cs, A, Bp, C, tr, S))
        Ssh = S.copy()
        bi = np.where(lab == B)[0]
        rng = np.random.default_rng([cc.SEED, s_idx, c_i, 7])
        bcols = [col[p] for p in progs if p.startswith(ch[1] + ":")]
        Ssh[np.ix_(bi, bcols)] = S[np.ix_(rng.permutation(bi), bcols)]
        out.append(("shuffledB", cs, A, B, C, tr, Ssh))
    return out


def flatten(variant, cs, res):
    return {f"{variant}|{cs}|{t[0]}|{t[1]}|{t[2]}|{st}": v for t, r in res.items() for st, v in r.items()}


def run_section(section, s_idx, ann, n_perm=N_PERM):
    sdir = cc.out_dir(f"{PHASE}/sections_{ann}")
    out = sdir / f"{section}.npz"
    if out.exists():
        return str(out)
    log = cc.start_run(f"chain_phase3_section_{section}_{ann}", {"n_perm": n_perm})
    d = cf.load_cache(section)
    lab = d[ann].astype(str)
    S = d["scores"].astype(np.float64)
    col = {n: j for j, n in enumerate(d["score_names"])}
    progs = list(cf.programs_used())
    ctx = ce.SectionCtx(d["xy"], d["W"], d["R"], d["total_counts"], d["lineage"], d["lineage_names"])
    rb, _ = random_b_choice(ann)
    variants = build_variants(section, s_idx, ann, d, lab, S, col, progs, rb, floor_pass(ann), ctx)
    if not variants:
        np.savez_compressed(out, keys=np.array([], dtype="U"), obs=np.array([]), obs_spill=np.array([]),
                            null_expr=np.zeros((n_perm, 0)), null_label=np.zeros((n_perm, 0)))
        cc.finish_run(log, PHASE)
        return str(out)
    obs, obs_sp, geoms = {}, {}, {}
    for v in variants:
        r, g = ce.chain_stats(ctx, lab, v[6], col, v[2], v[3], v[4], v[5], spill=False, return_geom=True)
        geoms[(v[0], v[1])] = g
        obs.update(flatten(v[0], v[1], r))
        obs_sp.update(flatten(v[0], v[1], ce.chain_stats(ctx, lab, v[6], col, v[2], v[3], v[4], v[5], spill=True,
                                                          geom=g)))
    keys = sorted(obs)
    ne = np.full((n_perm, len(keys)), np.nan)
    nl = np.full((n_perm, len(keys)), np.nan)
    ki = {k: i for i, k in enumerate(keys)}
    for b in range(n_perm):
        pe = ctx.expr_perm(lab, np.random.default_rng([cc.SEED, s_idx, 1, b]))
        pl = ctx.label_perm(np.random.default_rng([cc.SEED, s_idx, 2, b]))
        lab_l = lab[pl]
        for v in variants:
            g = geoms[(v[0], v[1])]
            if g is not None:
                for k, val in flatten(v[0], v[1], ce.chain_stats(ctx, lab, v[6][pe], col, v[2], v[3], v[4], v[5],
                                                                  geom=g)).items():
                    ne[b, ki[k]] = val
            for k, val in flatten(v[0], v[1], ce.chain_stats(ctx, lab_l, v[6], col, v[2], v[3], v[4], v[5])).items():
                nl[b, ki[k]] = val
        if b % 50 == 0:
            log.info("%s %s perm %d/%d %s", section, ann, b, n_perm, cc.mem_status())
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, keys=np.array(keys), obs=np.array([obs[k] for k in keys]),
                        obs_spill=np.array([obs_sp.get(k, np.nan) for k in keys]), null_expr=ne, null_label=nl)
    tmp.replace(out)
    cc.register(out, PHASE, f"Per-section observed stats + {n_perm} expression/label block-null replicates ({ann})")
    cc.finish_run(log, PHASE)
    return str(out)


# ----------------------------------------------------------------------------- tumour level
def aggregate(ann, log):
    import v3_nulls as vn
    s2t = vn.section_to_tumour()
    sdir = cc.RESULTS / PHASE / f"sections_{ann}"
    per = {}
    for s in sorted(s2t):
        z = np.load(sdir / f"{s}.npz")
        per[s] = {"keys": [str(k) for k in z["keys"]], "obs": z["obs"], "obs_spill": z["obs_spill"],
                  "ne": z["null_expr"], "nl": z["null_label"]}
    allkeys = sorted({k for p in per.values() for k in p["keys"]})
    rng = np.random.default_rng(cc.SEED)
    sec_rows, rows = [], []
    for key in allkeys:
        variant, cs, pa, pb, pc, st = key.split("|")
        if st.startswith("n_"):
            for s, p in per.items():
                if key in p["keys"]:
                    j = p["keys"].index(key)
                    sec_rows.append({"annotation": ann, "variant": variant, "chain": cs, "prog_A": pa, "prog_B": pb,
                                     "prog_C": pc, "stat": st, "section": s, "tumour": s2t[s], "obs": p["obs"][j]})
            continue
        tum_obs, tum_sp, tum_ne, tum_nl = {}, {}, {}, {}
        for s, p in per.items():
            if key not in p["keys"]:
                continue
            j = p["keys"].index(key)
            o = p["obs"][j]
            if not np.isfinite(o):
                continue
            t = s2t[s]
            sec_rows.append({"annotation": ann, "variant": variant, "chain": cs, "prog_A": pa, "prog_B": pb,
                             "prog_C": pc, "stat": st, "section": s, "tumour": t, "obs": o,
                             "obs_spill": p["obs_spill"][j],
                             "p_expr_section": cc.add_one_p(p["ne"][:, j], o), "p_label_section": cc.add_one_p(p["nl"][:, j], o)})
            for dct, val in ((tum_obs, o), (tum_sp, p["obs_spill"][j]), (tum_ne, p["ne"][:, j]), (tum_nl, p["nl"][:, j])):
                dct.setdefault(t, []).append(val)
        tums = sorted(tum_obs)
        if not tums:
            continue
        est = np.array([np.mean(tum_obs[t]) for t in tums])
        est_sp = np.array([np.nanmean(tum_sp[t]) for t in tums])
        T = est.mean()
        # paired null: replicate b = mean over tumours of mean over sections of replicate b (NaN-robust)
        ne = np.array([np.nanmean(np.vstack(tum_ne[t]), axis=0) for t in tums])
        nl = np.array([np.nanmean(np.vstack(tum_nl[t]), axis=0) for t in tums])
        null_e, null_l = np.nanmean(ne, axis=0), np.nanmean(nl, axis=0)
        k = rng.integers(0, len(tums), size=(N_BOOT, len(tums)))
        bo, bs = est[k].mean(axis=1), est_sp[k].mean(axis=1)
        rows.append({"annotation": ann, "variant": variant, "chain": cs, "prog_A": pa, "prog_B": pb, "prog_C": pc,
                     "stat": st, "n_tumours": len(tums), "n_sections": int(sum(len(tum_obs[t]) for t in tums)),
                     "tumours": ";".join(tums), "estimate": T, "ci_low": np.quantile(bo, .025),
                     "ci_high": np.quantile(bo, .975), "n_tumours_pos": int((est > 0).sum()),
                     "p_expr": cc.add_one_p(null_e, T), "p_label": cc.add_one_p(null_l, T),
                     "n_null_expr": int(np.isfinite(null_e).sum()), "n_null_label": int(np.isfinite(null_l).sum()),
                     "null_mean_expr": float(np.nanmean(null_e)), "null_mean_label": float(np.nanmean(null_l)),
                     "estimate_spill": float(np.nanmean(est_sp)), "ci_low_spill": np.nanquantile(bs, .025),
                     "ci_high_spill": np.nanquantile(bs, .975)})
    sec = pd.DataFrame(sec_rows)
    tl = pd.DataFrame(rows)
    # BH per statistic family x null across chains and programs (forward variant = the tested family);
    # negative-control variants get their own BH within variant. Duplicated hop tests are removed first.
    tl["q_expr"] = np.nan
    tl["q_label"] = np.nan
    ntests = {}
    for variant in tl.variant.unique():
        for st in ["beta_CB", "beta_BA", "indirect"]:
            m = (tl.variant == variant) & (tl.stat == st)
            d = tl[m]
            if st == "beta_CB":
                ident = d.chain.str.split("->").str[1:].str.join("->") + "|" + d.prog_B + "|" + d.prog_C
            elif st == "beta_BA":
                ident = d.chain.str.split("->").str[:2].str.join("->") + "|" + d.prog_A + "|" + d.prog_B
            else:
                ident = d.chain + "|" + d.prog_A + "|" + d.prog_B + "|" + d.prog_C
            first = ~ident.duplicated()
            for null in ("expr", "label"):
                q = cc.bh(d.loc[first, f"p_{null}"].to_numpy())
                qmap = dict(zip(ident[first], q))
                tl.loc[m, f"q_{null}"] = ident.map(qmap).to_numpy()
            ntests[f"{variant}|{st}"] = int(first.sum())
    return sec, tl, ntests


def decide(tl, calib, ann, sec):
    """Apply the plan decision rule per (chain, triple) for each variant.

    c5 (calibrated range, Phase 4): median analysed C cells per section >= calib min_n_chain AND
    |tumour-level indirect| >= mean simulated |indirect| at the smallest trustworthy fold change. Calibration was
    run on the primary chain's planted triple only; for other chains/triples c5 is an extrapolation (flagged).
    """
    rows = []
    nch = sec[sec.stat == "n_chain"].groupby(["variant", "chain", "prog_A", "prog_B", "prog_C"]).obs.median()
    key = ["variant", "chain", "prog_A", "prog_B", "prog_C"]
    w = tl.pivot_table(index=key, columns="stat", values=["estimate", "ci_low", "ci_high", "p_expr", "p_label",
                                                           "q_expr", "estimate_spill", "ci_low_spill", "ci_high_spill",
                                                           "n_tumours"], aggfunc="first")
    for idx, r in w.iterrows():
        g = lambda f, s: r[(f, s)] if (f, s) in r.index else np.nan  # noqa: E731
        c1 = g("q_expr", "beta_CB") < .05 and g("p_label", "beta_CB") < .05
        c2 = g("q_expr", "beta_BA") < .05 and g("p_label", "beta_BA") < .05
        ci_ex0 = g("ci_low", "indirect") > 0 or g("ci_high", "indirect") < 0
        shrink = abs(g("estimate", "beta_direct")) < abs(g("estimate", "beta_total"))
        c3 = g("q_expr", "indirect") < .05 and ci_ex0 and shrink
        sp_same = np.sign(g("estimate_spill", "indirect")) == np.sign(g("estimate", "indirect"))
        sp_ci = g("ci_low_spill", "indirect") > 0 or g("ci_high_spill", "indirect") < 0
        c4 = bool(sp_same and sp_ci)
        med_n = float(nch.get(idx, np.nan))
        if calib and calib.get("trustworthy_setting_found"):
            c5 = bool(med_n >= calib["min_n_chain_trustworthy"] and
                      abs(g("estimate", "indirect")) >= calib["min_abs_indirect_trustworthy"])
        else:
            c5 = False
        ntum = g("n_tumours", "indirect")
        if not (ntum >= 3):
            label = "not testable"
        elif c1 and c2 and c3 and c4 and c5:
            label = "supported (within calibrated range)"
        elif c1 and c2 and c3:
            label = "exploratory"
        else:
            label = "not supported"
        rows.append({"annotation": ann, **dict(zip(key, idx)), "n_tumours": ntum,
                     "beta_BA": g("estimate", "beta_BA"), "beta_CB": g("estimate", "beta_CB"),
                     "amplification_abs_CB_over_BA": abs(g("estimate", "beta_CB")) / abs(g("estimate", "beta_BA"))
                     if g("estimate", "beta_BA") not in (0, np.nan) else np.nan,
                     "beta_total": g("estimate", "beta_total"), "beta_direct": g("estimate", "beta_direct"),
                     "indirect": g("estimate", "indirect"), "indirect_ci_low": g("ci_low", "indirect"),
                     "indirect_ci_high": g("ci_high", "indirect"),
                     "q_expr_CB": g("q_expr", "beta_CB"), "p_label_CB": g("p_label", "beta_CB"),
                     "q_expr_BA": g("q_expr", "beta_BA"), "p_label_BA": g("p_label", "beta_BA"),
                     "q_expr_indirect": g("q_expr", "indirect"), "p_label_indirect": g("p_label", "indirect"),
                     "indirect_spill": g("estimate_spill", "indirect"),
                     "indirect_spill_ci": f"[{g('ci_low_spill', 'indirect'):.4f}, {g('ci_high_spill', 'indirect'):.4f}]",
                     "c1_hop_CB": c1, "c2_hop_BA": c2, "c3_indirect_shrink": c3, "c4_spillover_robust": c4,
                     "median_n_chain_per_section": med_n, "c5_calibrated": c5,
                     "c5_note": "" if (idx[1] == "tumour->CAF->T") else "calibration extrapolated from primary chain triple",
                     "label": label, "flag": "uses Amendment 1 (type-aware background rule)"})
    calls = pd.DataFrame(rows)
    fwd = calls[calls.variant == "forward"].copy()
    neg = calls[calls.variant != "forward"]
    notes = []
    for r in fwd.itertuples():
        m = neg[(neg.chain == r.chain) & (neg.c1_hop_CB & neg.c2_hop_BA & neg.c3_indirect_shrink)]
        if r.variant == "forward":
            # reversed variant uses swapped programs (pc, pb, pa); match on program set
            mm = m[((m.variant == "reversed") & (m.prog_A == r.prog_C) & (m.prog_B == r.prog_B) & (m.prog_C == r.prog_A))
                   | ((m.variant != "reversed") & (m.prog_A == r.prog_A) & (m.prog_B == r.prog_B) & (m.prog_C == r.prog_C))]
        notes.append(";".join(sorted(mm.variant)) if len(mm) else "")
    fwd["negative_controls_meeting_1to3"] = notes
    fwd.loc[(fwd.negative_controls_meeting_1to3 != "") & (fwd.label != "not supported"), "label"] += \
        " - not specific (negative control also passes)"
    return calls, fwd


def figure(tl, ann, log):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = tl[(tl.variant == "forward") & tl.stat.isin(["beta_BA", "beta_CB", "indirect"])].copy()
    d["triple"] = d.chain + "\n" + d.prog_A.str.split(":").str[1] + "/" + d.prog_B.str.split(":").str[1] + "/" + \
        d.prog_C.str.split(":").str[1]
    trip = sorted(d.triple.unique(), key=lambda x: (not x.startswith("tumour->CAF->T"), x))
    fig, ax = plt.subplots(figsize=(max(10, len(trip) * 0.35), 5.5))
    for k, (st, col) in enumerate([("beta_BA", "#2ca02c"), ("beta_CB", "#1f77b4"), ("indirect", "#d62728")]):
        e = d[d.stat == st].set_index("triple").reindex(trip)
        x = np.arange(len(trip)) + (k - 1) * 0.25
        ax.errorbar(x, e.estimate, yerr=[e.estimate - e.ci_low, e.ci_high - e.estimate], fmt="o", ms=3, color=col,
                    capsize=1.5, label={"beta_BA": "hop B<-A (upstream)", "beta_CB": "hop C<-B (downstream)",
                                        "indirect": "indirect a*b"}[st])
    ax.axhline(0, color="k", lw=.6)
    ax.set_xticks(range(len(trip)), trip, rotation=90, fontsize=5.5)
    ax.set_ylabel("tumour-level standardised coefficient (95% tumour bootstrap CI)")
    ax.set_title(f"Phase 3 chain test ({ann}): per-hop effects and indirect effect, forward chains", fontsize=10)
    ax.legend(fontsize=7)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, f"phase3_chain_effects_{ann}", f"Per-hop and indirect effects ({ann})", log)


def main():
    ann = sys.argv[sys.argv.index("--annotation") + 1] if "--annotation" in sys.argv else "author_mapped"
    import v3_nulls as vn
    sections = sorted(vn.section_to_tumour())
    if "--aggregate-only" not in sys.argv:
        n_workers = min(6, int(os.environ.get("CHAIN_WORKERS", "6")))
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(run_section, s, i, ann): s for i, s in enumerate(sections)}
            for f in as_completed(futs):
                print("done", futs[f], f.result(), flush=True)
    log = cc.start_run(f"chain_phase3_aggregate_{ann}", {"n_perm": N_PERM, "n_boot": N_BOOT})
    rb, cands = random_b_choice(ann)
    sec, tl, ntests = aggregate(ann, log)
    calib_path = cc.RESULTS / "phase4" / "calibrated_range.json"
    calib = json.loads(calib_path.read_text()) if calib_path.exists() else None
    in_range = calib
    calls, fwd = decide(tl, calib, ann, sec)
    suffix = "" if calib else "_precalibration"
    cc.save_csv(sec, PHASE, f"section_level_stats_{ann}{suffix}.csv", "Per-section statistics with section-level add-one p (expr/label nulls)", log)
    cc.save_csv(tl, PHASE, f"tumour_level_tests_{ann}{suffix}.csv",
                "Tumour-level estimates, bootstrap CI, add-one p (paired block nulls, 200 perms), BH q per family", log)
    cc.save_csv(calls, PHASE, f"all_variant_calls_{ann}{suffix}.csv", "Decision-rule criteria for forward and negative-control variants", log)
    cc.save_csv(fwd, PHASE, f"chain_calls_{ann}{suffix}.csv", "Forward-chain calls with negative-control annotation", log)
    hop = tl[(tl.variant == "forward") & tl.stat.isin(["beta_BA", "beta_CB"])].pivot_table(
        index=["chain", "prog_A", "prog_B", "prog_C"], columns="stat", values=["estimate", "ci_low", "ci_high"]).reset_index()
    hop.columns = ["_".join([c for c in col if c]) for col in hop.columns]
    hop["abs_CB_over_abs_BA"] = hop.estimate_beta_CB.abs() / hop.estimate_beta_BA.abs()
    cc.save_csv(hop, PHASE, f"per_hop_effects_{ann}{suffix}.csv", "Per-hop tumour-level effects (amplification visible/absent)", log)
    cc.save_json({"n_tests_after_dedup_per_variant_family": ntests, "n_perm": N_PERM, "random_B_choice": rb,
                  "random_B_candidates": cands, "calibration_file_present": bool(calib), "calibration": in_range},
                 PHASE, f"n_tests_{ann}{suffix}.json", "Number of tests per BH family; random-B choices", log)
    figure(tl, ann, log)
    log.info("forward calls:\n%s", fwd[["chain", "prog_A", "prog_B", "prog_C", "beta_BA", "beta_CB", "indirect",
                                         "q_expr_CB", "q_expr_BA", "q_expr_indirect", "label"]].to_string(index=False))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

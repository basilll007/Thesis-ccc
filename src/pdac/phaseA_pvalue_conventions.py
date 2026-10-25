"""Phase A (v2): put both empirical p-value conventions side by side for breast and PDAC.

  legacy  p = mean(null >= obs)            (breast paper files, src/spatial_null.py; can be exactly 0)
  add-one p = (#(null >= obs) + 1)/(n + 1) (PDAC extension, Phase 1+; minimum 1/(n+1))

Verdict = p < 0.05 (uncorrected) and BH q < 0.05 within the family/section the original analysis used.
Where null distributions were saved they are re-read and p is recomputed from scratch (and compared with the
stored columns); otherwise the stored p_legacy/p_plus1 columns (written by the same ccc_audit.stats call) are used.

  breast  : Phase 1 families paper/extra/leiden (6-NN, 500 perms; nulls in results/pdac/phase1/*.npz);
            15 um contact-restricted multi-axis recheck (200 perms; results/pdac/phase0/multi_axis_15um_recheck.csv)
  PDAC    : U5-a 6-NN marker + Leiden (500 perms; nulls in results/pdac/phase3/*_null_6NN_*.npz);
            16 sections marker/Leiden (cross_section_6NN_all_rows.csv) and author (author_annotation_6NN_all_sections.csv)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS_PDAC, finish_run, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import ccc_audit as ca  # noqa: E402

PHASE = "phaseA"
ALPHA = 0.05


def bh_by(df, col, groups):
    if not groups:
        return pd.Series(ca.bh(df[col].to_numpy()), index=df.index)
    return df.groupby(groups)[col].transform(lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))


def long_table(df, dataset, family, groups, id_cols, scores=ca.SCORES, recomputed=None):
    """One row per (interaction, score) with both conventions, verdicts, BH under both, flip flags."""
    out = []
    for key in scores:
        pl, pp = f"{key}_p_legacy", f"{key}_p_plus1"
        if pl not in df:
            continue
        d = df[id_cols].copy()
        d["dataset"], d["family"], d["score"] = dataset, family, key
        d["z"] = df.get(f"{key}_z")
        d["n_null"] = df.get(f"{key}_n_null_valid", np.nan)
        d["p_legacy"], d["p_plus1"] = df[pl].to_numpy(), df[pp].to_numpy()
        d["q_bh_legacy"] = bh_by(df.assign(_p=df[pl]), "_p", groups).to_numpy()
        d["q_bh_plus1"] = bh_by(df.assign(_p=df[pp]), "_p", groups).to_numpy()
        d["sig_legacy"], d["sig_plus1"] = d.p_legacy < ALPHA, d.p_plus1 < ALPHA
        d["sig_bh_legacy"], d["sig_bh_plus1"] = d.q_bh_legacy < ALPHA, d.q_bh_plus1 < ALPHA
        d["flip_uncorrected"] = d.sig_legacy != d.sig_plus1
        d["flip_bh"] = d.sig_bh_legacy != d.sig_bh_plus1
        d["p_source"] = "recomputed from saved nulls" if recomputed else "stored columns"
        out.append(d)
    return pd.concat(out, ignore_index=True)


def recompute_from_npz(df, npz_path, log, tag):
    z = np.load(npz_path)
    diffs = {}
    for key in ca.SCORES:
        nm = z[f"{key}_null"].astype(np.float64)
        obs = df[f"{key}_obs"].to_numpy(dtype=np.float64)
        st = ca.stats(obs, nm)
        for conv in ("p_legacy", "p_plus1"):
            stored = df[f"{key}_{conv}"].to_numpy()
            ok = np.isfinite(stored) & np.isfinite(st[conv])
            diffs[f"{key}_{conv}_max_abs_diff_vs_stored"] = float(np.max(np.abs(stored[ok] - st[conv][ok]))) if ok.any() else np.nan
            diffs[f"{key}_{conv}_n_verdict_diff_vs_stored"] = int(((stored[ok] < ALPHA) != (st[conv][ok] < ALPHA)).sum())
            df[f"{key}_{conv}"] = st[conv]
        df[f"{key}_n_null_valid"] = st["n_null_valid"]
    log.info("[%s] recomputed from %s: %s", tag, npz_path.name, diffs)
    return df, diffs


def main() -> None:
    log = start_run("phaseA_pvalue_conventions")
    tabs, checks = [], {}
    id4 = ["sender", "receiver", "ligand", "receptor"]

    # ---------------- breast, Phase 1 families (6-NN, 500 perms) ----------------
    for fam, csv in (("paper", "breast_matched_baseline_paper22.csv"),
                     ("extra", "breast_matched_baseline_extra_CD86CTLA4_PTNSDC4.csv"),
                     ("leiden", "breast_matched_baseline_leiden.csv")):
        df = pd.read_csv(RESULTS_PDAC / "phase1" / csv)
        npz = RESULTS_PDAC / "phase1" / f"null_distributions_{fam}.npz"
        lab = np.load(npz)["interaction_labels"]
        assert list(lab) == [f"{a}|{b}|{c}|{d}" for a, b, c, d in df[id4].itertuples(index=False)], "npz order mismatch"
        df, diffs = recompute_from_npz(df, npz, log, f"breast/{fam}")
        checks[f"breast_{fam}"] = diffs
        tabs.append(long_table(df, "breast_xenium_6NN", fam, [], id4, recomputed=True))

    # breast gate under both conventions
    p22 = pd.read_csv(RESULTS_PDAC / "phase1" / "breast_matched_baseline_paper22.csv")
    gate = {}
    for conv in ("p_legacy", "p_plus1"):
        ev = p22[p22.edge_status == "evaluated"]
        c, e, m = (ev[f"comp_{conv}"] < ALPHA), (ev[f"edge_{conv}"] < ALPHA), (ev[f"matched_{conv}"] < ALPHA)
        dis = c != e
        gate[conv] = {"n_C_E_disagreements": int(dis.sum()), "n_score_form": int((dis & (m == e)).sum())}
    log.info("breast Phase 1 gate under both conventions: %s", gate)

    # ---------------- breast, 15 um contact-restricted recheck (200 perms; paper setting) --------------
    r15 = pd.read_csv(RESULTS_PDAC / "phase0" / "multi_axis_15um_recheck.csv")
    r15 = r15.rename(columns={"p_value_legacy_mean_ge": "edge_p_legacy", "p_value_plus1": "edge_p_plus1", "z_score": "edge_z"})
    r15["edge_n_null_valid"] = 200
    tabs.append(long_table(r15, "breast_xenium_15um_radius", "multi_axis_15um", [], ["axis", "sender", "receiver"],
                           scores=("edge",)))

    # ---------------- PDAC U5-a 6-NN (nulls saved) ----------------
    u5 = pd.read_csv(RESULTS_PDAC / "phase3" / "mendeley_pdac_smi__U5-a_6NN_scores_both_annotations.csv")
    for ann in ("cell_type_marker", "cell_type_leiden"):
        d = u5[u5.annotation == ann].reset_index(drop=True).copy()
        d, diffs = recompute_from_npz(d, RESULTS_PDAC / "phase3" / f"mendeley_pdac_smi__U5-a_null_6NN_{ann}.npz", log,
                                      f"U5-a/{ann}")
        checks[f"pdac_U5-a_{ann}"] = diffs
        tabs.append(long_table(d, "pdac_U5-a_6NN", ann, [], id4, recomputed=True))

    # ---------------- PDAC 16 sections (stored columns; nulls not saved) ----------------
    xs = pd.read_csv(RESULTS_PDAC / "phase3" / "cross_section_6NN_all_rows.csv")
    for ann in ("cell_type_marker", "cell_type_leiden"):
        d = xs[xs.annotation == ann].reset_index(drop=True)
        tabs.append(long_table(d, "pdac_16_sections_6NN", ann, ["section"], ["section"] + id4))
    au = pd.read_csv(RESULTS_PDAC / "phase3" / "author_annotation_6NN_all_sections.csv")
    tabs.append(long_table(au, "pdac_16_sections_6NN", "author_mapped", ["section"], ["section"] + id4))

    allt = pd.concat(tabs, ignore_index=True)
    breast = allt[allt.dataset.str.startswith("breast")]
    pdac = allt[allt.dataset.str.startswith("pdac")]
    save_csv(breast, PHASE, "pvalue_conventions_breast.csv",
             "Breast: legacy vs add-one p, uncorrected and BH verdicts per interaction x score, flip flags", log)
    save_csv(pdac, PHASE, "pvalue_conventions_pdac.csv.gz",
             "PDAC: legacy vs add-one p, uncorrected and BH (within section) verdicts per interaction x score, flip flags", log)
    summ = allt.groupby(["dataset", "family", "score"]).agg(
        n_rows=("p_plus1", lambda x: int(np.isfinite(x).sum())),
        n_sig_legacy=("sig_legacy", "sum"), n_sig_plus1=("sig_plus1", "sum"),
        n_flip_uncorrected=("flip_uncorrected", "sum"),
        n_sig_bh_legacy=("sig_bh_legacy", "sum"), n_sig_bh_plus1=("sig_bh_plus1", "sum"),
        n_flip_bh=("flip_bh", "sum"), min_n_null=("n_null", "min")).reset_index()
    save_csv(summ, PHASE, "pvalue_conventions_flip_summary.csv",
             "Per dataset x family x score: significant counts under each convention and verdict flips", log)
    flips = allt[allt.flip_uncorrected | allt.flip_bh]
    save_csv(flips, PHASE, "pvalue_conventions_flipped_rows.csv", "Every row whose verdict differs between conventions", log)

    # headline rows, both conventions
    hl = pdac[(pdac.ligand == "CXCL12") & (pdac.receptor == "CXCR4") & (pdac.sender == "CAF_fibroblast")
              & (pdac.receiver == "T_cell") & (pdac.dataset == "pdac_16_sections_6NN")]
    hls = hl.groupby(["family", "score"]).agg(n_sections=("section", "nunique"), sig_legacy=("sig_legacy", "sum"),
                                              sig_plus1=("sig_plus1", "sum"), sig_bh_legacy=("sig_bh_legacy", "sum"),
                                              sig_bh_plus1=("sig_bh_plus1", "sum")).reset_index()
    save_csv(hls, PHASE, "pvalue_conventions_headline_CAF_T.csv",
             "PDAC CAF->T CXCL12-CXCR4: sections significant under each convention (uncorrected, BH within section)", log)

    res = {"alpha": ALPHA, "breast_gate_both_conventions": gate, "recompute_vs_stored_checks": checks,
           "n_rows_total": int(len(allt)), "n_flip_uncorrected": int(allt.flip_uncorrected.sum()),
           "n_flip_bh": int(allt.flip_bh.sum()),
           "note": "With n=500 valid nulls, legacy p<0.05 <=> #ge<=24 <=> add-one p<0.05, so uncorrected verdicts cannot "
                   "flip; with n=200 both thresholds are #ge<=9. Flips are possible only when n_null_valid < 500 (ratio "
                   "score, NaN nulls) or after BH, where legacy p=0 vs add-one p=1/(n+1) changes the ranks/q-values."}
    save_json(res, PHASE, "pvalue_conventions_summary.json", "Summary of p-value convention comparison", log)
    log.info("TOTAL rows %d; uncorrected flips %d; BH flips %d", len(allt), res["n_flip_uncorrected"], res["n_flip_bh"])
    for r in summ.itertuples():
        if r.n_flip_uncorrected or r.n_flip_bh:
            log.info("flips: %s", r._asdict())

    # figure: legacy vs add-one p, all rows, coloured by flip
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (name, d) in zip(axes, (("breast", breast), ("PDAC", pdac))):
        ok = np.isfinite(d.p_legacy) & np.isfinite(d.p_plus1)
        d = d[ok]
        col = np.where(d.flip_uncorrected, "#d62728", np.where(d.flip_bh, "#ff7f0e", "#9e9e9e"))
        ax.scatter(d.p_legacy + 1e-4, d.p_plus1, s=6, c=col, alpha=0.6, rasterized=True)
        ax.axvline(ALPHA, color="k", ls=":", lw=0.8)
        ax.axhline(ALPHA, color="k", ls=":", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("legacy p = mean(null >= obs)  (+1e-4 for log axis)")
        ax.set_ylabel("add-one p = (k+1)/(n+1)")
        ax.set_title(f"{name}: {len(d)} rows; flips uncorrected {int(d.flip_uncorrected.sum())} (red), "
                     f"BH only {int((d.flip_bh & ~d.flip_uncorrected).sum())} (orange)", fontsize=9)
    fig.suptitle("Phase A: p-value conventions side by side (verdict at 0.05)")
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseA_pvalue_conventions_scatter", "Legacy vs add-one p for breast and PDAC rows", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

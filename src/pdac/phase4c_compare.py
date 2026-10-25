"""Phase 4c: compare our scores with Squidpy ligrec, COMMOT and SpatialDM.

Cluster level (L, R, sender, receiver), crop (all methods) and full section (ours vs ligrec):
  Spearman correlation between method scores, top-20 overlap (Jaccard), significance agreement (p < 0.05).
L-R level: per pair, ours = max z over sender-receiver pairs; ligrec = max mean; COMMOT = total signal;
  SpatialDM = global z (CellChatDB entry whose ligand set is exactly {L} and whose receptor set contains R;
  exact single-gene receptor match preferred). Mapping table saved.
Missing method outputs are reported as BLOCKED/FAILED, never imputed.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import finish_run, out_dir, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

PHASE = "phase4"
K4 = ["ligand", "receptor", "sender", "receiver"]
TOPN = 20


def _read(name):
    p = out_dir(PHASE) / name
    return pd.read_csv(p) if p.exists() else None


def cluster_compare(tag, log):
    ours = _read(f"ours_CEMQ_6NN_{tag}.csv")
    lr = _read(f"squidpy_ligrec_{tag}.csv")
    cm = _read("commot_cluster_crop.csv") if tag == "crop" else None
    m = ours[K4 + ["comp_z", "matched_z", "edge_z", "ratio_z", "comp_sig", "edge_sig", "matched_sig", "ratio_sig",
                   "n_qualifying_edges"]].copy()
    if lr is not None:
        m = m.merge(lr, on=K4, how="left")
        m["ligrec_sig"] = m["ligrec_pvalue"] < 0.05
    if cm is not None:
        m = m.merge(cm, on=K4, how="left")
        m["commot_sig"] = m["commot_cluster_pvalue"] < 0.05
    scores = {"ours_comp_z": "comp_z", "ours_matched_z": "matched_z", "ours_edge_z": "edge_z", "ours_ratio_z": "ratio_z"}
    if lr is not None:
        scores["ligrec_mean"] = "ligrec_mean"
    if cm is not None:
        scores["commot_cluster_score"] = "commot_cluster_score"
    names = list(scores)
    rho = pd.DataFrame(index=names, columns=names, dtype=float)
    nobs = pd.DataFrame(index=names, columns=names, dtype=float)
    jac = pd.DataFrame(index=names, columns=names, dtype=float)
    tops = {}
    for a in names:
        v = m[scores[a]]
        tops[a] = set(m.loc[v.dropna().sort_values(ascending=False).index[:TOPN]].apply(lambda r: tuple(r[K4]), axis=1))
    for a in names:
        for b in names:
            x, y = m[scores[a]], m[scores[b]]
            ok = np.isfinite(x) & np.isfinite(y)
            rho.loc[a, b] = spearmanr(x[ok], y[ok]).statistic if ok.sum() > 2 else np.nan
            nobs.loc[a, b] = ok.sum()
            jac.loc[a, b] = len(tops[a] & tops[b]) / len(tops[a] | tops[b]) if tops[a] | tops[b] else np.nan
    sigcols = {"ours_comp": "comp_sig", "ours_matched": "matched_sig", "ours_edge": "edge_sig", "ours_ratio": "ratio_sig"}
    if lr is not None and lr["ligrec_pvalue"].notna().any():
        sigcols["ligrec"] = "ligrec_sig"
    elif lr is not None:
        log.warning("[%s] all squidpy ligrec p-values are NaN -> ligrec excluded from significance agreement (BLOCKED)", tag)
    if cm is not None:
        sigcols["commot"] = "commot_sig"
    agree = []
    for a, ca_ in sigcols.items():
        for b, cb in sigcols.items():
            if a >= b:
                continue
            x, y = m[ca_], m[cb]
            ok = x.notna() & y.notna()
            x, y = x[ok].astype(bool), y[ok].astype(bool)
            po = float((x == y).mean())
            pe = float(x.mean() * y.mean() + (1 - x.mean()) * (1 - y.mean()))
            agree.append({"method_a": a, "method_b": b, "n_rows": int(ok.sum()), "n_sig_a": int(x.sum()), "n_sig_b": int(y.sum()),
                          "n_both_sig": int((x & y).sum()), "agreement": po,
                          "cohen_kappa": (po - pe) / (1 - pe) if pe < 1 else np.nan})
    save_csv(m, PHASE, f"comparison_per_axis_{tag}.csv", f"Per-axis (L,R,sender,receiver) scores of all methods ({tag})", log)
    save_csv(rho, PHASE, f"comparison_spearman_{tag}.csv", f"Spearman between method scores across rows ({tag})", log, index=True)
    save_csv(nobs, PHASE, f"comparison_spearman_n_{tag}.csv", f"Rows used for each Spearman ({tag})", log, index=True)
    save_csv(jac, PHASE, f"comparison_top{TOPN}_jaccard_{tag}.csv", f"Jaccard overlap of top-{TOPN} rows per method ({tag})", log, index=True)
    save_csv(pd.DataFrame(agree), PHASE, f"comparison_significance_agreement_{tag}.csv",
             f"Pairwise agreement/kappa of p<0.05 calls ({tag})", log)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, mat, t in ((axes[0], rho, "Spearman ρ"), (axes[1], jac, f"top-{TOPN} Jaccard")):
        im = ax.imshow(mat.to_numpy(float), vmin=-1 if t.startswith("S") else 0, vmax=1, cmap="RdBu_r" if t.startswith("S") else "viridis")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=8)
        for i in range(len(names)):
            for j in range(len(names)):
                ax.text(j, i, f"{mat.iat[i, j]:.2f}", ha="center", va="center", fontsize=7, color="k")
        ax.set_title(f"{t} ({tag}, cluster-level rows)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    save_figure(fig, PHASE, f"phase4_method_agreement_cluster_level_{tag}", f"Method agreement heatmaps ({tag})", log)
    return {"tag": tag, "n_rows": int(len(m)), "spearman": rho.round(3).to_dict(), "top20_jaccard": jac.round(3).to_dict(),
            "significance_agreement": agree}


def lr_compare(log):
    out = {}
    for tag in ("crop", "full"):
        ours = _read(f"ours_CEMQ_6NN_{tag}.csv")
        lr = _read(f"squidpy_ligrec_{tag}.csv")
        g = ours.groupby(["ligand", "receptor"]).agg(ours_max_comp_z=("comp_z", "max"), ours_max_edge_z=("edge_z", "max"),
                                                      ours_max_matched_z=("matched_z", "max"), ours_max_ratio_z=("ratio_z", "max")).reset_index()
        if lr is not None:
            g = g.merge(lr.groupby(["ligand", "receptor"]).agg(ligrec_max_mean=("ligrec_mean", "max")).reset_index(), how="left")
        if tag == "crop":
            cg = _read("commot_global_crop.csv")
            if cg is not None:
                g = g.merge(cg[["ligand", "receptor", "commot_total_signal"]], how="left")
        sd = _read(f"spatialdm_global_{tag}.csv")
        mapping = []
        if sd is not None:
            zs = []
            for l, r in g[["ligand", "receptor"]].itertuples(index=False):
                c = sd[(sd["ligand_genes"] == l) & sd["receptor_genes"].str.split(";").apply(lambda s: r in s)]
                exact = c[c["receptor_genes"] == r]
                pick = exact.iloc[0] if len(exact) else (c.iloc[0] if len(c) else None)
                zs.append(np.nan if pick is None else pick["z"])
                mapping.append({"ligand": l, "receptor": r, "tag": tag,
                                "cellchat_id": None if pick is None else pick["cellchat_id"],
                                "cellchat_ligand_genes": None if pick is None else pick["ligand_genes"],
                                "cellchat_receptor_genes": None if pick is None else pick["receptor_genes"],
                                "match": "none" if pick is None else ("exact" if len(exact) else "receptor-in-complex")})
            g["spatialdm_global_z"] = zs
            save_csv(pd.DataFrame(mapping), PHASE, f"spatialdm_mapping_{tag}.csv",
                     f"Mapping of catalogue L-R pairs to SpatialDM/CellChatDB entries ({tag})", log)
        save_csv(g, PHASE, f"comparison_lr_level_{tag}.csv", f"L-R-level summary scores of all methods ({tag})", log)
        cols = [c for c in g.columns if c not in ("ligand", "receptor")]
        rho = pd.DataFrame(index=cols, columns=cols, dtype=float)
        for a in cols:
            for b in cols:
                ok = np.isfinite(g[a]) & np.isfinite(g[b])
                rho.loc[a, b] = spearmanr(g.loc[ok, a], g.loc[ok, b]).statistic if ok.sum() > 2 else np.nan
        save_csv(rho, PHASE, f"comparison_lr_level_spearman_{tag}.csv", f"Spearman across L-R pairs ({tag})", log, index=True)
        out[tag] = {"n_pairs": int(len(g)), "n_spatialdm_mapped": int(np.isfinite(g.get("spatialdm_global_z", pd.Series(dtype=float))).sum()),
                    "spearman": rho.round(3).to_dict()}
    return out


def main() -> None:
    log = start_run("phase4c_compare")
    status = {}
    p = out_dir(PHASE) / "phase4b_status.json"
    if p.exists():
        import json
        status = json.loads(p.read_text())
    summary = {"phase4b_status": status, "cluster_level": [cluster_compare("crop", log), cluster_compare("full", log)],
               "lr_level": lr_compare(log)}
    save_json(summary, PHASE, "phase4_agreement_summary.json", "Agreement summary (Spearman, top-20 Jaccard, kappa)", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

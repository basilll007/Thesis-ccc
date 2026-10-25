"""Phase B1 (v2): protein-guided annotation of the PDAC-SMI sections (rule: docs/PDAC_V2_PROTEIN_RULE.md).

Per section (processed one at a time; per-section outputs are written immediately and re-used on restart):
  * per-channel thresholds on log1p(Mean.CD3 / Mean.CD45 / Mean.PanCK): GMM (primary) or median+3*MAD (fallback,
    and for every channel in the 'protein_mad' sensitivity variant), with diagnostics
  * labels protein_gmm, protein_mad (gating + RNA argmax within compartment), v1 marker labels (recomputed with the
    identical v1 rule), author labels (raw + v1 AUTHOR_MAP)
After all sections: confusion matrices and per-class precision/recall vs author labels (U5-a and all 16 sections),
threshold diagnostics table, figures. No CCC is run here.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA_DIR, SEED, finish_run, mem_status, out_dir, register, save_csv, save_figure,  # noqa: E402
                    save_json, start_run)

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from sklearn.mixture import GaussianMixture  # noqa: E402

import pdac_loaders  # noqa: E402
from phase3_pdac_audit import PDAC_MARKERS, qc  # noqa: E402
from phase3d_author_annotation import AUTHOR_MAP  # noqa: E402

PHASE = "phaseB"
CHANNELS = ["CD3", "CD45", "PanCK"]
GMM_MAX_CELLS = 50_000
ASHMAN_MIN = 2.0
MIN_WEIGHT = 0.005
IMMUNE_NONT = ["myeloid", "B_plasma", "mast"]
EPI = ["epithelial_ductal_tumor", "acinar", "endocrine"]
NONIMMUNE = EPI + ["CAF_fibroblast", "endothelial", "pericyte_smc"]
AUTHOR_T = {"CD8+ T", "CD4+ T", "Treg"}
SEC_DIR = "phaseB/sections"


def mad_threshold(x):
    med = np.median(x)
    return float(med + 3 * 1.4826 * np.median(np.abs(x - med)))


def gmm_threshold(x, rng):
    xs = x if len(x) <= GMM_MAX_CELLS else x[rng.choice(len(x), GMM_MAX_CELLS, replace=False)]
    g = GaussianMixture(n_components=2, covariance_type="full", n_init=5, random_state=SEED).fit(xs.reshape(-1, 1))
    mu = g.means_.ravel()
    sd = np.sqrt(g.covariances_.ravel())
    w = g.weights_.ravel()
    hi, lo = int(np.argmax(mu)), int(np.argmin(mu))
    grid = np.linspace(mu[lo], mu[hi], 2000)
    post = g.predict_proba(grid.reshape(-1, 1))[:, hi]
    above = np.nonzero(post >= 0.5)[0]
    t = float(grid[above[0]]) if len(above) else np.nan
    ashman = float(np.sqrt(2) * abs(mu[hi] - mu[lo]) / np.sqrt(sd[hi] ** 2 + sd[lo] ** 2))
    return {"gmm_mu_neg": float(mu[lo]), "gmm_mu_pos": float(mu[hi]), "gmm_sd_neg": float(sd[lo]),
            "gmm_sd_pos": float(sd[hi]), "gmm_w_neg": float(w[lo]), "gmm_w_pos": float(w[hi]),
            "ashman_d": ashman, "gmm_threshold": t, "n_fit": int(len(xs))}


def rna_scores(adata):
    scores = {}
    for t, genes in PDAC_MARKERS.items():
        g = [x for x in genes if x in adata.var_names]
        if len(g) >= 2:
            scores[t] = np.asarray(adata[:, g].X.mean(axis=1)).ravel()
    return pd.DataFrame(scores, index=adata.obs_names)


def gate_labels(pos, sf):
    P, H, T = pos["PanCK"], pos["CD45"], pos["CD3"]
    n = len(P)
    lab = np.empty(n, dtype=object)

    def argmax_in(cands, mask):
        c = [x for x in cands if x in sf.columns]
        return sf.loc[mask, c].idxmax(axis=1).to_numpy()

    amb = P & (H | T)
    tcell = T & H & ~P
    t45neg = T & ~H & ~P
    imm = H & ~T & ~P
    epi = P & ~H & ~T
    neg = ~P & ~H & ~T
    lab[amb] = "protein_ambiguous"
    lab[tcell] = "T_cell"
    lab[t45neg] = "T_cell_CD45neg"
    if imm.any():
        lab[imm] = argmax_in(IMMUNE_NONT, imm)
    if epi.any():
        lab[epi] = argmax_in(EPI, epi)
    if neg.any():
        lab[neg] = argmax_in(NONIMMUNE, neg)
    assert all(isinstance(v, str) for v in lab)
    return lab


def process_section(s, log):
    t0 = time.time()
    adata, meta = pdac_loaders.load(f"mendeley_pdac_smi:{s}", DATA_DIR, log)
    adata, _, _ = qc(adata, log)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sf = rna_scores(adata)
    marker = sf.idxmax(axis=1).to_numpy()  # identical rule to v1 phase3_pdac_audit.annotate (cell_type_marker)
    rng = np.random.default_rng(SEED)
    diag, pos_gmm, pos_mad = [], {}, {}
    for c in CHANNELS:
        x = np.log1p(adata.obs[f"Mean.{c}"].to_numpy(dtype=float))
        g = gmm_threshold(x, rng)
        tm = mad_threshold(x)
        accept = (g["ashman_d"] >= ASHMAN_MIN) and (min(g["gmm_w_neg"], g["gmm_w_pos"]) >= MIN_WEIGHT) and np.isfinite(g["gmm_threshold"])
        t_primary = g["gmm_threshold"] if accept else tm
        pos_gmm[c] = x > t_primary
        pos_mad[c] = x > tm
        diag.append({"section": s, "channel": c, **g, "gmm_accepted": bool(accept), "mad_threshold": tm,
                     "threshold_primary": t_primary, "primary_method": "gmm" if accept else "mad_fallback",
                     "frac_pos_primary": float(pos_gmm[c].mean()), "frac_pos_mad": float(pos_mad[c].mean()),
                     "median_x": float(np.median(x)), "p99_x": float(np.percentile(x, 99))})
    lab_gmm = gate_labels(pos_gmm, sf)
    lab_mad = gate_labels(pos_mad, sf)
    author = adata.obs["author_majortype"].astype(str).to_numpy()
    cells = pd.DataFrame({
        "section": s, "cell_id": adata.obs_names.astype(str), "treatment_status": meta["treatment_status"],
        "log1p_CD3": np.log1p(adata.obs["Mean.CD3"].to_numpy(float)),
        "log1p_CD45": np.log1p(adata.obs["Mean.CD45"].to_numpy(float)),
        "log1p_PanCK": np.log1p(adata.obs["Mean.PanCK"].to_numpy(float)),
        "CD3_pos_gmm": pos_gmm["CD3"], "CD45_pos_gmm": pos_gmm["CD45"], "PanCK_pos_gmm": pos_gmm["PanCK"],
        "CD3_pos_mad": pos_mad["CD3"], "CD45_pos_mad": pos_mad["CD45"], "PanCK_pos_mad": pos_mad["PanCK"],
        "protein_gmm": lab_gmm, "protein_mad": lab_mad, "marker": marker,
        "author_majortype": author, "author_mapped": pd.Series(author).map(AUTHOR_MAP).to_numpy(),
    })
    d = out_dir(SEC_DIR)
    cells.to_csv(d / f"{s}_labels.csv.gz", index=False)
    pd.DataFrame(diag).to_csv(d / f"{s}_thresholds.csv", index=False)
    register(d / f"{s}_labels.csv.gz", PHASE, f"Per-cell protein-guided, marker and author labels, section {s}")
    register(d / f"{s}_thresholds.csv", PHASE, f"Per-channel threshold diagnostics, section {s}")
    log.info("section %s: %d cells; T_cell(gmm)=%d T_cell(mad)=%d; methods=%s; %.1fs; %s", s, len(cells),
             int((lab_gmm == "T_cell").sum()), int((lab_mad == "T_cell").sum()),
             {r["channel"]: r["primary_method"] for r in diag}, time.time() - t0, mem_status())
    del adata
    return cells, pd.DataFrame(diag)


def prf(pred, truth_set_mask, pred_label):
    pm = pred == pred_label
    tp = int((pm & truth_set_mask).sum())
    return {"n_pred": int(pm.sum()), "n_true": int(truth_set_mask.sum()), "tp": tp,
            "precision": tp / pm.sum() if pm.sum() else np.nan, "recall": tp / truth_set_mask.sum() if truth_set_mask.sum() else np.nan}


def validation(cells, scope):
    """Per-class precision/recall of each annotation vs author labels (mapped to coarse types)."""
    rows = []
    am = cells["author_mapped"].to_numpy()
    araw = cells["author_majortype"].to_numpy()
    classes = ["T_cell", "myeloid", "B_plasma", "mast", "CAF_fibroblast", "epithelial_ductal_tumor", "endothelial",
               "pericyte_smc", "endocrine", "acinar"]
    for ann in ("protein_gmm", "protein_mad", "marker"):
        pred = cells[ann].to_numpy()
        for c in classes:
            truth = am == c
            rows.append({"scope": scope, "annotation": ann, "class": c, "truth": f"author_mapped=={c}", **prf(pred, truth, c)})
        # T vs author T/NK
        rows.append({"scope": scope, "annotation": ann, "class": "T_cell", "truth": "author T or NK",
                     **prf(pred, np.isin(araw, list(AUTHOR_T | {"Natural killer"})), "T_cell")})
        # loose T (protein variants only): T_cell + T_cell_CD45neg
        if ann.startswith("protein"):
            loose = np.where(np.isin(pred, ["T_cell", "T_cell_CD45neg"]), "T_loose", pred)
            rows.append({"scope": scope, "annotation": ann + "_looseT", "class": "T_cell (CD3+ PanCK-, any CD45)",
                         "truth": "author_mapped==T_cell", **prf(loose, am == "T_cell", "T_loose")})
    return pd.DataFrame(rows)


def main() -> None:
    log = start_run("phaseB1_protein_labels", extra_packages=[])
    log.info("Rule: docs/PDAC_V2_PROTEIN_RULE.md; GMM max cells %d; Ashman D >= %.1f; min weight %.3f; %s",
             GMM_MAX_CELLS, ASHMAN_MIN, MIN_WEIGHT, mem_status())
    full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    sections = sorted(full.obs["sample_id"].astype(str).unique())
    del full
    allc, alld = [], []
    d = out_dir(SEC_DIR)
    for s in sections:
        lp, tp = d / f"{s}_labels.csv.gz", d / f"{s}_thresholds.csv"
        if lp.exists() and tp.exists():
            log.info("section %s: re-using saved per-section outputs (resume)", s)
            allc.append(pd.read_csv(lp))
            alld.append(pd.read_csv(tp))
            continue
        c, dg = process_section(s, log)
        allc.append(c)
        alld.append(dg)
    cells = pd.concat(allc, ignore_index=True)
    diag = pd.concat(alld, ignore_index=True)
    save_csv(diag, PHASE, "protein_threshold_diagnostics.csv",
             "Per section x channel: GMM fit, Ashman D, accepted?, thresholds (GMM / MAD), fraction positive", log)

    # counts per label
    cnt = pd.concat([cells.groupby("section")[a].value_counts().rename("n").reset_index().rename(columns={a: "label"})
                     .assign(annotation=a) for a in ("protein_gmm", "protein_mad", "marker", "author_mapped")])
    save_csv(cnt, PHASE, "label_counts_by_section.csv", "Cells per label per section for each annotation", log)

    # validation: U5-a and all sections pooled, and per section
    # per-section scopes (U5-a is one of them) + pooled
    val = pd.concat([validation(cells[cells.section == s], s) for s in sections]
                    + [validation(cells, "all_16_sections")], ignore_index=True)
    save_csv(val, PHASE, "validation_precision_recall.csv",
             "Per-class precision/recall vs author labels for protein_gmm, protein_mad and marker (U5-a, pooled, per section)", log)
    for scope in ("U5-a", "all_16_sections"):
        sub = cells if scope == "all_16_sections" else cells[cells.section == scope]
        for ann in ("protein_gmm", "protein_mad", "marker"):
            cm = pd.crosstab(sub["author_majortype"], sub[ann], margins=True)
            save_csv(cm, PHASE, f"confusion_{ann}_vs_author_{scope}.csv",
                     f"Confusion matrix: author annotation_majortypes (rows) vs {ann} (columns), {scope}", log, index=True)
    key = val[(val["class"].str.startswith("T_cell"))].copy()
    for r in key[key.scope.isin(["U5-a", "all_16_sections"])].itertuples():
        log.info("VALIDATION %s %s %s vs %s: n_pred=%d n_true=%d tp=%d precision=%.3f recall=%.3f", r.scope, r.annotation,
                 r._3, r.truth, r.n_pred, r.n_true, r.tp, r.precision, r.recall)

    # summary json
    def pr(scope, ann, cls="T_cell", truth="author_mapped==T_cell"):
        r = val[(val.scope == scope) & (val.annotation == ann) & (val["class"] == cls) & (val.truth == truth)].iloc[0]
        return {"n_pred": int(r.n_pred), "tp": int(r.tp), "n_true": int(r.n_true), "precision": float(r.precision),
                "recall": float(r.recall)}
    summ = {"n_cells": int(len(cells)), "sections": sections,
            "gmm_accepted_by_channel": diag.groupby("channel").gmm_accepted.sum().astype(int).to_dict(),
            "T_cell_vs_author_T": {sc_: {a: pr(sc_, a) for a in ("protein_gmm", "protein_mad", "marker")}
                                   for sc_ in ("U5-a", "all_16_sections")}}
    save_json(summ, PHASE, "protein_labels_summary.json", "Summary: GMM acceptance, T-cell precision/recall vs author", log)

    make_figures(cells, diag, val, sections, log)
    finish_run(log, PHASE)


def make_figures(cells, diag, val, sections, log):
    u5 = cells[cells.section == "U5-a"]
    dg = diag[diag.section == "U5-a"].set_index("channel")
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for ax, c in zip(axes[0], CHANNELS):
        x = u5[f"log1p_{c}"]
        ax.hist(x, bins=150, color="#bbbbbb", density=True, label="all cells")
        at = u5.author_majortype.isin(AUTHOR_T)
        ax.hist(x[at], bins=150, histtype="step", color="#d62728", density=True, label="author T (CD8/CD4/Treg)")
        ax.hist(x[u5.author_majortype == "Malignant"], bins=150, histtype="step", color="#1f77b4", density=True,
                label="author Malignant")
        ax.axvline(dg.loc[c, "threshold_primary"], color="k", lw=1.2, label=f"primary ({dg.loc[c, 'primary_method']})")
        ax.axvline(dg.loc[c, "mad_threshold"], color="k", ls="--", lw=1, label="median+3MAD")
        ax.set_xlabel(f"log1p(Mean.{c})")
        ax.set_title(f"U5-a {c}: Ashman D={dg.loc[c, 'ashman_d']:.2f}", fontsize=10)
        ax.legend(fontsize=7)
    ax = axes[1, 0]
    smp = u5.sample(n=min(len(u5), 40000), random_state=SEED)
    ax.scatter(smp.log1p_CD45, smp.log1p_CD3, s=1, c="#cccccc", rasterized=True, label="all (40k sample)")
    tt = smp[smp.author_majortype.isin(AUTHOR_T)]
    ax.scatter(tt.log1p_CD45, tt.log1p_CD3, s=3, c="#d62728", rasterized=True, label="author T")
    ax.axvline(dg.loc["CD45", "threshold_primary"], color="k", lw=1)
    ax.axhline(dg.loc["CD3", "threshold_primary"], color="k", lw=1)
    ax.set_xlabel("log1p(Mean.CD45)")
    ax.set_ylabel("log1p(Mean.CD3)")
    ax.legend(fontsize=7, markerscale=4)
    ax.set_title("U5-a: CD45 vs CD3 with primary thresholds")
    ax = axes[1, 1]
    v = val[(val["class"] == "T_cell") & (val.truth == "author_mapped==T_cell") & val.scope.isin(sections)]
    for ann, col in (("marker", "#7f7f7f"), ("protein_gmm", "#d62728"), ("protein_mad", "#ff7f0e")):
        vv = v[v.annotation == ann].set_index("scope").reindex(sections)
        ax.scatter(vv.recall, vv.precision, c=col, label=ann, s=30)
    ax.set_xlabel("recall vs author T")
    ax.set_ylabel("precision vs author T")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.set_title("T_cell precision/recall per section (16)")
    ax = axes[1, 2]
    vv = val[(val.scope == "all_16_sections") & (val.truth.str.startswith("author_mapped==")) &
             val.annotation.isin(["marker", "protein_gmm"])]
    piv = vv.pivot_table(index="class", columns="annotation", values="precision", dropna=False)
    pir = vv.pivot_table(index="class", columns="annotation", values="recall", dropna=False)
    idx = sorted(set(piv.index) | set(pir.index))
    piv, pir = piv.reindex(idx), pir.reindex(idx)
    y = np.arange(len(piv))
    ax.barh(y - 0.2, piv["marker"], height=0.2, color="#bbbbbb", label="marker precision")
    ax.barh(y, piv["protein_gmm"], height=0.2, color="#d62728", label="protein_gmm precision")
    ax.barh(y + 0.2, pir["protein_gmm"], height=0.2, color="#fcae91", label="protein_gmm recall")
    ax.set_yticks(y)
    ax.set_yticklabels(piv.index, fontsize=8)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=7)
    ax.set_title("Pooled 16 sections, per class vs author")
    fig.suptitle("Phase B1: protein-guided labels — threshold diagnostics and validation against author labels")
    fig.tight_layout()
    save_figure(fig, PHASE, "phaseB1_protein_gating_validation", "Protein channel separation, thresholds, and validation", log)


if __name__ == "__main__":
    main()

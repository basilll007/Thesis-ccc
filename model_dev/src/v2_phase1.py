"""v2 Phase 1: segmentation-spillover control for the CAF -> macrophage signal (plan: plans/v2_phase1_spillover.md).

Subcommands (run in order, one at a time):
  prepare    gene statistics (CAF vs isolated / touching macrophages), contamination markers M, spillover-free genes
             S, target-gene removal, per-receiver arrays (distance bins, 3-nearest-CAF block, contamination score,
             new targets and their tile means, protein flags), descriptive parts (b), (e), (f).
  fit        model groups on the GPU (resumable per group x fold; --fold-budget N for foreground chunks).
  summarise  R2 / neighbour gain / CAF-expression contribution Q with tumour-bootstrap CIs, per bin and per subset;
             decision JSON; figure.
Why: see the plan. Development set only; GSE310352 is not touched.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")   # CPU workers capped at 4 (set before numpy is imported)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import v2common as V  # noqa: E402
import features as F  # noqa: E402
import assemble as A  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "v2_phase1"
OUT = C.RESULTS / PHASE
CACHE = V.V2_CACHE / "phase1"
RECV = V.V2_CACHE / "phase1_receivers.npz"
BINS = ["touching", "near", "separated"]
LINEAGE_MYELOID = ["CD68", "CD14", "LYZ", "CSF1R", "FCGR3A"]
PROGRAM_GENES = {"C1Q_TAM": ["C1QA", "C1QB", "C1QC", "CD163", "MRC1", "SELENOP"],
                 "SPP1_TAM": ["SPP1", "GPNMB", "CD9", "FABP5", "LGALS3"]}
HYPOTHESIS_GENES = ["SELENOP", "CD163", "C1QA", "SPP1"]
N_BOOT = 2000


def diameter() -> float:
    return float(json.loads((OUT / "preplan" / "preplan_summary.json").read_text())["median_eq_diameter_um"])


def sections():
    return sorted((C.CACHE / "datasets" / "mendeley_pdac_smi").glob("*.npz"))


# ----------------------------------------------------------------------------- prepare
def gene_statistics(log) -> tuple[pd.DataFrame, list[str]]:
    """Per tumour and class (CAF, isolated macrophage, touching macrophage): mean linear CP10k and detection."""
    D = diameter()
    acc: dict = {}
    genes = None
    for f in sections():
        z = np.load(f, allow_pickle=False)
        lab, grp = z["label"].astype(str), str(z["group"])
        genes = [str(g) for g in z["genes"]]
        d, _ = V.nearest_type_distance(z["xy"].astype(float), lab, V.CAF, k=1)
        dmin = d[:, 0]
        masks = {"CAF": lab == V.CAF, "iso_macrophage": (lab == V.MYELOID) & (dmin > 2 * D),
                 "touching_macrophage": (lab == V.MYELOID) & (dmin <= D)}
        for cls, m in masks.items():
            if not m.any():
                continue
            lin = np.expm1(z["X_log"][m].astype(np.float32)).astype(np.float64)
            det = (z["X_counts"][m] > 0).sum(0).astype(np.float64)
            s, dd, n = acc.get((grp, cls), (0.0, 0.0, 0))
            acc[(grp, cls)] = (s + lin.sum(0), dd + det, n + int(m.sum()))
    rows = []
    for (grp, cls), (s, dd, n) in acc.items():
        rows.append(pd.DataFrame({"tumour": grp, "class": cls, "gene": genes, "mean_cp10k": s / n, "det_frac": dd / n,
                                  "n_cells": n}))
    long = pd.concat(rows, ignore_index=True)
    log.info("gene statistics: %d tumours x %d classes x %d genes", long.tumour.nunique(), long["class"].nunique(), len(genes))
    return long, genes


def summarise_genes(long: pd.DataFrame) -> pd.DataFrame:
    med = long.groupby(["gene", "class"])[["mean_cp10k", "det_frac"]].median().unstack("class")
    med.columns = [f"{a}__{b}" for a, b in med.columns]
    med = med.reset_index()
    with np.errstate(divide="ignore", invalid="ignore"):
        med["ratio_caf_over_iso"] = med["mean_cp10k__CAF"] / med["mean_cp10k__iso_macrophage"]
        med["ratio_iso_over_caf"] = med["mean_cp10k__iso_macrophage"] / med["mean_cp10k__CAF"]
    return med


def select_sets(med: pd.DataFrame, targets: list[str], log) -> dict:
    out = {}
    excl = set(targets) | set(LINEAGE_MYELOID)
    for thr in (10.0, 5.0):
        m = med[(~med.gene.isin(excl)) & (med["det_frac__CAF"] >= 0.25) & (med.ratio_caf_over_iso >= thr)]
        if len(m) >= 3:
            break
    m = m.sort_values("mean_cp10k__CAF", ascending=False).head(20)
    out["M"] = {"genes": m.gene.tolist(), "ratio_threshold": thr, "n_passing_before_cap": None}
    for thr_s in (10.0, 5.0):
        s = med[(med["det_frac__iso_macrophage"] >= 0.25) & (med.ratio_iso_over_caf >= thr_s)]
        if len(s) >= 3:
            break
    s = s.sort_values("mean_cp10k__iso_macrophage", ascending=False).head(20)
    out["S"] = {"genes": s.gene.tolist(), "ratio_threshold": thr_s}
    removed = med[med.gene.isin(targets) & (med.ratio_caf_over_iso >= 0.5)].gene.tolist()
    out["removed_target_genes"] = removed
    out["retained"] = {p: [g for g in gs if g not in removed] for p, gs in PROGRAM_GENES.items()}
    out["estimable"] = {p: len(v) >= 2 for p, v in out["retained"].items()}
    out["gene_vector_adj"] = [g for g in targets if g not in removed]
    out["blocked"] = {"b_c_i": len(out["M"]["genes"]) < 3, "d": len(out["S"]["genes"]) < 3}
    log.info("M (%d, ratio>=%g): %s", len(out["M"]["genes"]), thr, out["M"]["genes"])
    log.info("S (%d, ratio>=%g): %s", len(out["S"]["genes"]), thr_s, out["S"]["genes"])
    log.info("removed target genes (CAF/iso >= 0.5): %s; retained %s", removed, out["retained"])
    return out


def prepare(args) -> None:
    log = C.start_run(PHASE, "v2_phase1_prepare")
    V.memory_gate(log)
    D = diameter()
    spec = A.targets_spec()["myeloid"]
    targets = spec["genes"]
    long, genes = gene_statistics(log)
    gidx = {g: i for i, g in enumerate(genes)}
    CACHE.mkdir(parents=True, exist_ok=True)
    long.to_csv(V.V2_CACHE / "phase1_gene_stats_per_tumour.csv.gz", index=False)
    V.register_v2_artefact(V.V2_CACHE / "phase1_gene_stats_per_tumour.csv.gz", PHASE, "per-tumour gene statistics (CAF, isolated/touching macrophages)")
    med = summarise_genes(long)
    C.save_table(med, PHASE, "gene_statistics_median_over_tumours.csv", log)
    sets = select_sets(med, targets, log)
    C.save_json(sets, PHASE, "gene_sets_and_adjusted_targets.json", log)
    M, S = sets["M"]["genes"], sets["S"]["genes"]
    Mi, Si = [gidx[g] for g in M], [gidx[g] for g in S]
    meta = V.cell_metadata(log).set_index(["section", "row"])
    # per-receiver arrays, in the v1 feature-cache row order
    z0 = np.load(C.CACHE / "features" / "development__myeloid.npz", allow_pickle=False)
    cache_section, cache_row = z0["section"].astype(str), z0["row"]
    new_t = {"C1Q_TAM_mean": PROGRAM_GENES["C1Q_TAM"], "SPP1_TAM_mean": PROGRAM_GENES["SPP1_TAM"],
             "C1Q_TAM_adj": sets["retained"]["C1Q_TAM"], "SPP1_TAM_adj": sets["retained"]["SPP1_TAM"],
             "S_mean": S, "M_mean": M}
    new_t = {k: v for k, v in new_t.items() if len(v) >= (2 if k.endswith("_adj") else 1)}
    if sets["blocked"]["d"]:
        new_t.pop("S_mean", None); log.warning("S has < 3 genes: part (d) BLOCKED")
    if sets["blocked"]["b_c_i"]:
        raise RuntimeError("M has < 3 genes: parts (b) and (c-i) BLOCKED (plan); stopping before model fits")
    pergene = [g for g in dict.fromkeys(S + M) if g not in targets]
    cols = {k: [] for k in ["section", "group", "row", "d_caf", "n_caf6", "s_c", "marker_cp10k", "caf3_marker_cp10k",
                            "area_um2", "cd45", "cd45_gt_caf75", "cd45_gt_caf90", "ncaf_dist_min", "ncaf_dist_mean3"]}
    ncaf, newY, newTile, pg = [], [], [], []
    caf_marker = {}  # tumour -> (sum marker CP10k over CAF cells, n)
    from scipy.spatial import cKDTree
    for f in sections():
        z = np.load(f, allow_pickle=False)
        sec, grp = str(z["section"]), str(z["group"])
        lab = z["label"].astype(str)
        xy = z["xy"].astype(float)
        n = len(lab)
        Xlog = z["X_log"]
        tot = z["total_counts"].astype(np.float64)
        mk = 1e4 * z["X_counts"][:, Mi].astype(np.float64).sum(1) / np.maximum(tot, 1)
        cafm = lab == V.CAF
        s_, n_ = caf_marker.get(grp, (0.0, 0))
        caf_marker[grp] = (s_ + float(mk[cafm].sum()), n_ + int(cafm.sum()))
        rec = np.where(lab == V.MYELOID)[0]
        sel = cache_section == sec
        assert np.array_equal(cache_row[sel], rec), f"{sec}: receiver order differs from the v1 feature cache"
        d3, j3 = V.nearest_type_distance(xy, lab, V.CAF, k=3)
        _, nn6 = cKDTree(xy).query(xy[rec], k=7)
        n_caf6 = (lab[nn6[:, 1:]] == V.CAF).sum(1)
        ncaf.append(np.stack([Xlog[j3[rec, i]].astype(np.float32) for i in range(3)]).mean(0))
        # new targets for all cells (tile LOO needs every same-type cell), then receivers
        full_new = np.column_stack([Xlog[:, [gidx[g] for g in gs]].astype(np.float32).mean(1) for gs in new_t.values()])
        tm, _ = F.tile_loo_mean(full_new, z["tile"], lab == V.MYELOID)
        newY.append(full_new[rec]); newTile.append(tm[rec]); pg.append(Xlog[rec][:, [gidx[g] for g in pergene]].astype(np.float32))
        mrow = meta.loc[sec]
        caf_cd45 = mrow.loc[np.where(cafm)[0], "Mean.CD45"]
        q75, q90 = float(caf_cd45.quantile(0.75)), float(caf_cd45.quantile(0.90))
        cd45 = mrow.loc[rec, "Mean.CD45"].to_numpy(float)
        cols["section"] += [sec] * len(rec); cols["group"] += [grp] * len(rec); cols["row"] += rec.tolist()
        cols["d_caf"] += d3[rec, 0].tolist(); cols["n_caf6"] += n_caf6.tolist()
        cols["s_c"] += np.log1p(mk[rec]).tolist(); cols["marker_cp10k"] += mk[rec].tolist()
        cols["caf3_marker_cp10k"] += mk[j3[rec]].mean(1).tolist()
        cols["area_um2"] += mrow.loc[rec, "area_um2"].tolist(); cols["cd45"] += cd45.tolist()
        cols["cd45_gt_caf75"] += (cd45 > q75).tolist(); cols["cd45_gt_caf90"] += (cd45 > q90).tolist()
        cols["ncaf_dist_min"] += d3[rec, 0].tolist(); cols["ncaf_dist_mean3"] += d3[rec].mean(1).tolist()
        log.info("%s: %d receivers", sec, len(rec))
    df = pd.DataFrame(cols)
    assert np.array_equal(df.section.to_numpy(), cache_section) and np.array_equal(df.row.to_numpy(), cache_row)
    df["bin"] = pd.cut(df.d_caf, [0, D, 2 * D, np.inf], labels=BINS, right=True).astype(str)
    np.savez_compressed(RECV, **{k: df[k].to_numpy() for k in df.columns if k not in ("section", "group", "bin")},
                        section=df.section.to_numpy().astype("U"), group=df.group.to_numpy().astype("U"),
                        bin=df.bin.to_numpy().astype("U"), ncaf_expr=np.vstack(ncaf).astype(np.float32),
                        new_targets=np.vstack(newY), new_target_names=np.array(list(new_t)),
                        new_tile=np.vstack(newTile), pergene=np.vstack(pg), pergene_names=np.array(pergene),
                        genes=np.array(genes))
    V.register_v2_artefact(RECV, PHASE, "per-receiver arrays: distance bins, 3-nearest-CAF block, contamination score, new targets, protein flags")
    C.save_json({"new_targets": new_t, "pergene_targets": pergene,
                 "caf_marker_cp10k_by_tumour": {g: v[0] / v[1] for g, v in caf_marker.items()},
                 "caf_marker_sum_n": {g: list(v) for g, v in caf_marker.items()}}, PHASE, "receiver_build.json", log)
    descriptives(df, long, sets, caf_marker, log)


def boot_mean_by_tumour(df: pd.DataFrame, value: str, n_boot: int = N_BOOT) -> tuple[float, float, float]:
    """Pooled mean over cells with a tumour-bootstrap CI (resample tumours, sum / count)."""
    g = df.groupby("group")[value].agg(["sum", "count"])
    s, c = g["sum"].to_numpy(), g["count"].to_numpy()
    rng = np.random.default_rng(C.SEED)
    idx = rng.integers(0, len(g), size=(n_boot, len(g)))
    b = s[idx].sum(1) / c[idx].sum(1)
    return float(s.sum() / c.sum()), float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))


def descriptives(df: pd.DataFrame, long: pd.DataFrame, sets: dict, caf_marker: dict, log) -> None:
    rows = []
    # (b) marker CP10k by bin and by number of CAFs among the 6 nearest cells
    df = df.copy()
    df["n_caf6_cat"] = np.where(df.n_caf6 >= 3, ">=3", df.n_caf6.astype(str))
    for col, levels in (("bin", BINS), ("n_caf6_cat", ["0", "1", "2", ">=3"])):
        for lv in levels:
            sub = df[df[col] == lv]
            m, lo, hi = boot_mean_by_tumour(sub, "marker_cp10k")
            rows.append({"stratifier": col, "level": lv, "n_cells": len(sub), "marker_cp10k": m, "ci_low": lo, "ci_high": hi})
    caf = pd.DataFrame([{"group": g, "sum": v[0], "count": v[1]} for g, v in caf_marker.items()]).set_index("group")
    rows.append({"stratifier": "CAF cells", "level": "all CAFs", "n_cells": int(caf["count"].sum()),
                 "marker_cp10k": float(caf["sum"].sum() / caf["count"].sum()), "ci_low": np.nan, "ci_high": np.nan})
    C.save_table(pd.DataFrame(rows), PHASE, "b_contamination_by_adjacency.csv", log)
    # spillover rate per bin: (m_bin - m_iso) / (m_caf - m_iso), tumour bootstrap (tumours resampled jointly)
    tum = sorted(df.group.unique())
    agg = {b: df[df.bin == b].groupby("group").marker_cp10k.agg(["sum", "count"]).reindex(tum, fill_value=0) for b in BINS}
    cafa = caf.reindex(tum)

    def rate(idx, b):
        mb = agg[b]["sum"].to_numpy()[idx].sum(-1) / agg[b]["count"].to_numpy()[idx].sum(-1)
        mi = agg["separated"]["sum"].to_numpy()[idx].sum(-1) / agg["separated"]["count"].to_numpy()[idx].sum(-1)
        mc = cafa["sum"].to_numpy()[idx].sum(-1) / cafa["count"].to_numpy()[idx].sum(-1)
        return (mb - mi) / (mc - mi)
    rng = np.random.default_rng(C.SEED)
    bidx = rng.integers(0, len(tum), size=(N_BOOT, len(tum)))
    rr = []
    for b in ["touching", "near"]:
        pt = float(rate(np.arange(len(tum)), b))
        bs = rate(bidx, b)
        rr.append({"bin": b, "spillover_rate": pt, "ci_low": float(np.quantile(bs, 0.025)), "ci_high": float(np.quantile(bs, 0.975)),
                   "definition": "(marker CP10k in bin - isolated) / (CAF - isolated)"})
    C.save_table(pd.DataFrame(rr), PHASE, "b_spillover_rate.csv", log)
    # Spearman per tumour (touching receivers): s_c vs 3-nearest-CAF marker CP10k
    from scipy.stats import spearmanr
    sp = []
    for g, d in df[df.bin == "touching"].groupby("group"):
        r = spearmanr(d.marker_cp10k, d.caf3_marker_cp10k).statistic if len(d) > 10 else np.nan
        sp.append({"tumour": g, "n_touching": len(d), "spearman_marker_vs_caf3_marker": r})
    sp = pd.DataFrame(sp)
    C.save_table(sp, PHASE, "b_marker_correlation_with_nearest_cafs.csv", log)
    # (e) descriptive: marker CP10k by CD45 status and by area quartile
    erows = []
    for flag in ["cd45_gt_caf75", "cd45_gt_caf90"]:
        for v in (True, False):
            sub = df[df[flag] == v]
            m, lo, hi = boot_mean_by_tumour(sub, "marker_cp10k")
            erows.append({"stratifier": flag, "level": str(v), "n_cells": len(sub), "frac_cells": len(sub) / len(df),
                          "marker_cp10k": m, "ci_low": lo, "ci_high": hi})
    df["area_q"] = pd.qcut(df.area_um2, 4, labels=["Q1 smallest", "Q2", "Q3", "Q4 largest"]).astype(str)
    for lv in ["Q1 smallest", "Q2", "Q3", "Q4 largest"]:
        sub = df[df.area_q == lv]
        m, lo, hi = boot_mean_by_tumour(sub, "marker_cp10k")
        erows.append({"stratifier": "area quartile", "level": lv, "n_cells": len(sub), "frac_cells": len(sub) / len(df),
                      "marker_cp10k": m, "ci_low": lo, "ci_high": hi})
    C.save_table(pd.DataFrame(erows), PHASE, "e_protein_area_descriptives.csv", log)
    # (f) hypothesis and target genes: CAF vs isolated vs touching macrophages
    spec = A.targets_spec()["myeloid"]["genes"]
    frows = []
    for gname in spec:
        d = long[long.gene == gname].pivot(index="tumour", columns="class", values="mean_cp10k")
        det = long[long.gene == gname].pivot(index="tumour", columns="class", values="det_frac")
        rc = d["CAF"].to_numpy(); ri = d["iso_macrophage"].to_numpy()
        idx = rng.integers(0, len(d), size=(N_BOOT, len(d)))
        with np.errstate(divide="ignore", invalid="ignore"):
            bs = np.median(rc[idx], 1) / np.median(ri[idx], 1)
        frows.append({"gene": gname, "hypothesis_gene": gname in HYPOTHESIS_GENES,
                      "removed_by_rule": gname in sets["removed_target_genes"],
                      "caf_mean_cp10k_median": float(np.median(rc)), "caf_mean_cp10k_min": float(rc.min()), "caf_mean_cp10k_max": float(rc.max()),
                      "iso_mac_mean_cp10k_median": float(np.median(ri)), "iso_mac_mean_cp10k_min": float(ri.min()), "iso_mac_mean_cp10k_max": float(ri.max()),
                      "touch_mac_mean_cp10k_median": float(d["touching_macrophage"].median()),
                      "caf_det_median": float(det["CAF"].median()), "iso_mac_det_median": float(det["iso_macrophage"].median()),
                      "touch_mac_det_median": float(det["touching_macrophage"].median()),
                      "ratio_caf_over_iso": float(np.median(rc) / np.median(ri)),
                      "ratio_ci_low": float(np.nanquantile(bs, 0.025)), "ratio_ci_high": float(np.nanquantile(bs, 0.975)),
                      "touch_minus_iso_cp10k": float(d["touching_macrophage"].median() - np.median(ri))})
    C.save_table(pd.DataFrame(frows), PHASE, "f_target_genes_caf_vs_macrophage.csv", log)
    log.info("descriptives saved")


# ----------------------------------------------------------------------------- fit
def load_design(log):
    """GPU design: v1 myeloid cache blocks + s_c + tile means of new targets + 3-nearest-CAF block."""
    import torch
    import gpu as GP
    data = A.build_receiver("development", "myeloid")
    r = np.load(RECV, allow_pickle=False)
    sets = json.loads((OUT / "gene_sets_and_adjusted_targets.json").read_text())
    genes = [str(g) for g in r["genes"]]
    new_names = [str(x) for x in r["new_target_names"]]
    build = json.loads((OUT / "receiver_build.json").read_text())
    blocks, names, prov = [], [], []
    for b in ["own", "tile", "composition", "nb_programs", "nb_genes"]:
        blocks.append(data.X[b]); names += data.names[b]; prov += data.prov[b]
        data.X[b] = None
    blocks.append(r["s_c"].astype(np.float32)[:, None]); names.append("s_c_contamination")
    prov.append(F.Provenance("self_covariate", genes=frozenset(sets["M"]["genes"])))
    blocks.append(r["new_tile"].astype(np.float32)); names += [f"tile_loo_mean__{t}" for t in new_names]
    prov += [F.Provenance("tile", genes=frozenset(build["new_targets"][t])) for t in new_names]
    blocks.append(r["ncaf_expr"]); names += [f"ncaf_expr__{g}" for g in genes]
    prov += [F.Provenance("neighbour", genes=frozenset([g])) for g in genes]
    blocks.append(np.column_stack([r["ncaf_dist_min"], r["ncaf_dist_mean3"]]).astype(np.float32))
    names += ["ncaf_dist_min", "ncaf_dist_mean3"]; prov += [F.Provenance("neighbour")] * 2
    Xg = torch.cat([GP.to_gpu(b) for b in blocks], 1)
    del blocks
    # targets: v1 programs + 11 genes, new targets, per-gene S/M targets not among the 11
    Y = np.column_stack([data.Y, r["new_targets"], r["pergene"]]).astype(np.float32)
    y_names = list(data.target_names) + new_names + [f"g:{g}" for g in r["pergene_names"]]
    log.info("design: %d rows x %d columns on GPU; %d target columns", Xg.shape[0], Xg.shape[1], Y.shape[1])
    return Xg, names, prov, data, r, sets, Y, y_names


def variants(names: list[str]) -> dict[str, np.ndarray]:
    kind = [V.column_kind(n) for n in names]
    ctype = [V.column_type(n) for n in names]
    own = [j for j, n in enumerate(names) if kind[j] == "own" and n != "s_c_contamination"]
    sc = [names.index("s_c_contamination")]
    nb = [j for j in range(len(names)) if kind[j] != "own"]
    caf_expr = {j for j in nb if ctype[j] == "CAF" and kind[j] == "expression"}
    full = own + nb
    minus = own + [j for j in nb if j not in caf_expr]
    return {"own": np.array(own), "full": np.array(full), "minuscaf": np.array(minus),
            "own_sc": np.array(own + sc), "full_sc": np.array(full + sc), "minuscaf_sc": np.array(minus + sc)}


def groups(names, data, r, sets, Y, y_names) -> list:
    import v2engine as E
    vv = variants(names)
    spec = A.targets_spec()["myeloid"]
    M = set(sets["M"]["genes"])
    new_names = [str(x) for x in r["new_target_names"]]
    build = json.loads((OUT / "receiver_build.json").read_text())
    all_genes = set(spec["genes"]) | {g for t in new_names for g in build["new_targets"][t]} | set(map(str, r["pergene_names"]))
    xgb_ref = [t for t in ["myeloid:C1Q_TAM", "myeloid:SPP1_TAM", "C1Q_TAM_mean", "SPP1_TAM_mean", "C1Q_TAM_adj",
                           "SPP1_TAM_adj", "S_mean", "M_mean"] if t in y_names]
    gs = [E.Group("ref", {k: vv[k] for k in ("own", "full", "minuscaf")}, Y, y_names, xgb_ref, ["full", "minuscaf"],
                  target_genes=all_genes, target_programs=set(spec["programs"]))]
    # with s_c: targets must not involve M genes (the guard would refuse them)
    keep = [j for j, t in enumerate(y_names) if not (t == "M_mean" or (t.startswith("g:") and t[2:] in M))]
    Ysc, nsc = Y[:, keep], [y_names[j] for j in keep]
    xgb_sc = [t for t in ["myeloid:C1Q_TAM", "myeloid:SPP1_TAM", "C1Q_TAM_adj", "SPP1_TAM_adj"] if t in nsc]
    sc_genes = set(spec["genes"]) | {g for t in new_names if t != "M_mean" for g in build["new_targets"][t]} | \
        {t[2:] for t in nsc if t.startswith("g:")}
    gs.append(E.Group("sc", {k: vv[k] for k in ("own_sc", "full_sc", "minuscaf_sc")}, Ysc, nsc, xgb_sc,
                      ["full_sc", "minuscaf_sc"], target_genes=sc_genes, target_programs=set(spec["programs"])))
    ncore = len(data.target_names)
    for flag in ("cd45_gt_caf75", "cd45_gt_caf90"):
        gs.append(E.Group(flag, {k: vv[k] for k in ("own", "full", "minuscaf")}, Y[:, :ncore], y_names[:ncore],
                          ["myeloid:C1Q_TAM"], ["full", "minuscaf"], rows=r[flag].astype(bool),
                          target_genes=set(spec["genes"]), target_programs=set(spec["programs"])))
    return gs


def fit(args) -> None:
    import v2engine as E
    log = C.start_run(PHASE, "v2_phase1_fit")
    V.memory_gate(log)
    Xg, names, prov, data, r, sets, Y, y_names = load_design(log)
    sec2fold = dict(zip(A.folds().section, A.folds().fold))
    for g in groups(names, data, r, sets, Y, y_names):
        if args.group and g.name not in args.group:
            continue
        log.info("group %s: variants %s; %d ridge targets; XGBoost %s on %s; rows %s", g.name,
                 {k: len(v) for k, v in g.variants.items()}, len(g.y_names), g.xgb_targets, g.xgb_variants,
                 "all" if g.rows is None else int(g.rows.sum()))
        E.run_group(Xg, names, prov, data.section, sec2fold, g, CACHE, log, args.fold_budget)
    log.info("fit finished")


# ----------------------------------------------------------------------------- summarise
def summarise(args) -> None:
    import v2engine as E
    log = C.start_run(PHASE, "v2_phase1_summarise")
    r = np.load(RECV, allow_pickle=False)
    sets = json.loads((OUT / "gene_sets_and_adjusted_targets.json").read_text())
    grp = r["group"].astype(str)
    bins = r["bin"].astype(str)
    rows = []
    for gname, sc in (("ref", ""), ("sc", "_sc"), ("cd45_gt_caf75", ""), ("cd45_gt_caf90", "")):
        parts = E.load_parts(CACHE, gname)
        ynames = list(parts[0]["y_names"])
        xt = list(parts[0]["xgb_targets"]) if "xgb_targets" in parts[0] else []
        Yfull = None
        # rebuild Y for this group from the saved names (targets are deterministic functions of the receiver data)
        Yfull = group_targets(ynames)
        targets = [t for t in ynames if not t.startswith("g:")]
        vec = [t for t in A.targets_spec()["myeloid"]["genes"] if t in ynames]
        vec_adj = [g for g in sets["gene_vector_adj"] if g in ynames]
        tsets = [(t, t) for t in targets] + [("gene_vector", vec), ("gene_vector_adj", vec_adj)] + \
                [(t, t) for t in ynames if t.startswith("g:")]
        for strat, levels in (("all", [None]), ("bin", BINS)):
            if strat == "bin" and gname != "ref":
                continue
            for lv in levels:
                mask = None if lv is None else (bins == lv)
                for tname, tcols in tsets:
                    if isinstance(tcols, list) and len(tcols) < 1:
                        continue
                    for model in ("ridge", "xgb"):
                        if model == "xgb" and (not isinstance(tcols, str) or tcols not in xt):
                            continue
                        full, minus, own = f"full{sc}", f"minuscaf{sc}", f"own{sc}"
                        rf = E.summary_row(parts, grp, model, full, tcols, Yfull, ref=("ridge", own, None), cell_mask=mask)
                        rm = E.summary_row(parts, grp, model, full, tcols, Yfull, ref=(model, minus, None), cell_mask=mask)
                        ro = E.summary_row(parts, grp, "ridge", own, tcols, Yfull, cell_mask=mask)
                        rows.append({"group": gname, "stratum": strat, "level": lv or "all", "target": tname, "model": model,
                                     "r2_own": ro["r2"], "r2_full": rf["r2"], "r2_full_ci_low": rf["ci_low"], "r2_full_ci_high": rf["ci_high"],
                                     "r2_minuscaf": rf["r2"] - rm["diff"],
                                     "neighbour_gain": rf["diff"], "gain_ci_low": rf["diff_ci_low"], "gain_ci_high": rf["diff_ci_high"],
                                     "Q_caf_expr": rm["diff"], "Q_ci_low": rm["diff_ci_low"], "Q_ci_high": rm["diff_ci_high"],
                                     "n_tumours": rf["n_tumours"],
                                     "n_cells": int(sum(len(p["test_idx"]) if mask is None else int(mask[p["test_idx"]].sum()) for p in parts))})
        log.info("summarised group %s", gname)
    res = pd.DataFrame(rows)
    C.save_table(res, PHASE, "phase1_summary.csv", log)
    decide(res, sets, log)


_TARGET_CACHE: dict = {}


def group_targets(ynames: list[str]) -> np.ndarray:
    """Target matrix for the given names, from the v1 cache and the Phase 1 receiver file."""
    if "Y" not in _TARGET_CACHE:
        z = np.load(C.CACHE / "features" / "development__myeloid.npz", allow_pickle=False)
        r = np.load(RECV, allow_pickle=False)
        spec = A.targets_spec()["myeloid"]
        names = spec["programs"] + spec["genes"] + [str(x) for x in r["new_target_names"]] + [f"g:{g}" for g in r["pergene_names"]]
        _TARGET_CACHE["Y"] = np.column_stack([z["Y"], r["new_targets"], r["pergene"]]).astype(np.float32)
        _TARGET_CACHE["names"] = names
    idx = [_TARGET_CACHE["names"].index(t) for t in ynames]
    return _TARGET_CACHE["Y"][:, idx]


def decide(res: pd.DataFrame, sets: dict, log) -> None:
    def get(group, target, model="xgb", stratum="all", level="all"):
        d = res[(res.group == group) & (res.target == target) & (res.model == model) & (res.stratum == stratum) & (res.level == level)]
        return None if d.empty else d.iloc[0]

    def pos(row):
        return row is not None and row.Q_ci_low > 0

    ref = get("ref", "myeloid:C1Q_TAM")
    bin_rows = {b: get("ref", "myeloid:C1Q_TAM", stratum="bin", level=b) for b in BINS}
    a_pass = pos(bin_rows["near"]) or pos(bin_rows["separated"])
    if sets["estimable"]["C1Q_TAM"]:
        c_row, c_target = get("sc", "C1Q_TAM_adj"), "C1Q_TAM_adj (XGBoost, with s_c)"
    elif len(sets["gene_vector_adj"]) >= 2:
        c_row, c_target = get("sc", "gene_vector_adj", model="ridge"), "gene_vector_adj (ridge, with s_c)"
    else:
        c_row, c_target = None, "not estimable"
    c_pass = pos(c_row)
    d_row = get("ref", "S_mean")
    d_pass = pos(d_row)
    spill_pattern = pos(bin_rows["touching"]) and not pos(bin_rows["near"]) and not pos(bin_rows["separated"])
    failed = [n for n, ok in (("a", a_pass), ("c", c_pass), ("d", d_pass)) if not ok]
    if ref is None or not pos(ref):
        verdict = "inconclusive: no reference CAF-expression signal in the aggregated model"
    elif not failed:
        verdict = "survives spillover control"
    elif ("c" in failed or "d" in failed) and spill_pattern:
        verdict = "explained by spillover"
    else:
        verdict = "inconclusive"

    def q(row):
        return None if row is None else {"Q": row.Q_caf_expr, "ci": [row.Q_ci_low, row.Q_ci_high], "gain": row.neighbour_gain,
                                         "gain_ci": [row.gain_ci_low, row.gain_ci_high], "r2_full": row.r2_full}
    out = {"verdict": verdict, "failed_checks": failed, "check_a_pass": a_pass, "check_c_pass": c_pass, "check_d_pass": d_pass,
           "a_separated_bin_strict_pass": pos(bin_rows["separated"]), "spillover_distance_pattern": spill_pattern,
           "reference_Q_C1Q_TAM": q(ref), "Q_by_bin": {b: q(v) for b, v in bin_rows.items()},
           "c_target": c_target, "c_adjusted_Q": q(c_row), "d_S_mean_Q": q(d_row), "d_M_mean_Q": q(get("ref", "M_mean")),
           "rule": "plans/v2_phase1_spillover.md, Decision rule"}
    C.save_json(out, PHASE, "phase1_decision.json", log)
    log.info("DECISION: %s (failed: %s)", verdict, failed)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "fit", "summarise"])
    ap.add_argument("--group", nargs="*")
    ap.add_argument("--fold-budget", type=int)
    args = ap.parse_args()
    {"prepare": prepare, "fit": fit, "summarise": summarise}[args.step](args)


if __name__ == "__main__":
    main()

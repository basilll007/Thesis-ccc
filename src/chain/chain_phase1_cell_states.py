"""Chain Phase 1: cell-state programs, background filter, cell counts, floors and outcome-free chain selection.

Plan: docs/chain_plans/phase1_cell_states.md (committed before this script was run).
One section at a time; per-section outputs are saved as they finish (resume skips finished sections).
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402  (sets sys.path for src/pdac)

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

PHASE = "phase1"
TYPES = {"tumour": "epithelial_ductal_tumor", "CAF": "CAF_fibroblast", "T": "T_cell", "myeloid": "myeloid"}
PROGRAMS = {  # (type key, program) -> declared genes (panel members), see plan section 2
    ("tumour", "classical"): ["CEACAM6", "AGR2", "S100P", "KRT20", "LCN2", "PIGR", "OLFM4", "SPINK1", "CEACAM1"],
    ("tumour", "basal_like"): ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17", "S100A2", "SERPINB5", "AREG", "SLC2A1",
                               "KRT15", "LY6D", "ITGA3", "MET"],
    ("tumour", "EMT"): ["VIM", "SNAI1", "SNAI2", "TWIST2", "FN1", "MMP2", "VCAN", "CAV1"],
    ("CAF", "myCAF"): ["ACTA2", "TAGLN", "MYL9", "TPM2", "TPM1", "COL11A1", "TWIST1", "SOX4"],
    ("CAF", "iCAF"): ["IL6", "CXCL1", "CXCL2", "PDGFRA", "CFD", "LMNA", "CCL2", "LIF", "CXCL8"],
    ("CAF", "CXCL12"): ["CXCL12"],
    ("T", "cytotoxic"): ["GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "GNLY", "IFNG", "CCL5", "CST7", "CD69"],
    ("T", "exhaustion"): ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1"],
    ("myeloid", "inflammatory"): ["IL1B", "TNF", "CCL3", "CCL4", "CXCL9", "CXCL10", "NFKBIA", "PTGS2"],
    ("myeloid", "C1Q_TAM"): ["C1QA", "C1QB", "C1QC", "CD163", "MRC1", "SELENOP"],
    ("myeloid", "SPP1_TAM"): ["SPP1", "GPNMB", "CD9", "FABP5", "LGALS3"],
}
SINGLE_GENE_EXEMPT = {("CAF", "CXCL12")}
LINEAGE = {"tumour": ["KRT19", "KRT8", "KRT18", "EPCAM", "KRT7"], "CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
           "T": ["CD3D", "CD3E", "CD3G", "CD2", "CD8A"], "myeloid": ["CD68", "CD14", "LYZ", "CSF1R", "FCGR3A"]}
BG_FOLD, BG_MIN_SECTIONS, MIN_GENES = 2.0, 8, 5
N_BINS, CTRL_SIZE, K, RADIUS = 25, 50, 6, 30.0
FLOOR_CELLS, FLOOR_EDGES = 100, 300
ANNOTS = ["author_mapped", "marker"]
CACHE_SUB = "mendeley"  # Amendment 1 run uses "mendeley_v1b"
RULE = "section"  # "section" (pre-declared) or "typeaware" (Amendment 1)
TYPEAWARE_MIN_CELLS, TYPEAWARE_FRAC = 50, 0.5
H5 = cc.DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad"


def prog_name(key):
    return f"{key[0]}:{key[1]}"


# ----------------------------------------------------------------------------- background
def background_flags(log) -> pd.DataFrame:
    out = cc.out_dir("phase1") / "gene_background_by_section.csv"
    if out.exists():
        return pd.read_csv(out)
    a = ad.read_h5ad(H5, backed="r")
    pt = a.var["probe_type"].astype(str).to_numpy()
    base, neg = np.where(pt == "Base")[0], np.where(pt == "Negprob")[0]
    rows = []
    for s in sorted(a.obs["sample_id"].astype(str).unique()):
        idx = np.where((a.obs["sample_id"].astype(str) == s).to_numpy())[0]
        X = a[idx].to_memory().X
        X = X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
        det = np.asarray((X > 0).mean(axis=0)).ravel()
        nmax = float(det[neg].max())
        for j in base:
            rows.append({"section": s, "gene": a.var_names[j], "detection_rate": float(det[j]),
                         "negprobe_max_detection_rate": nmax, "near_background": bool(det[j] < BG_FOLD * nmax)})
        log.info("background %s: %d cells, neg-probe max detection %.4f", s, len(idx), nmax)
        del X
    df = pd.DataFrame(rows)
    cc.save_csv(df, PHASE, "gene_background_by_section.csv",
                "Per section x panel gene: detection rate vs 2x max negative-probe detection (near_background)", log)
    return df


def typeaware_flags(log) -> pd.DataFrame:
    """Amendment 1: per section x type: gene detection among own-type cells vs 2x max neg-probe detection in the same cells."""
    out = cc.out_dir(PHASE) / "gene_background_typeaware_by_section.csv"
    if out.exists():
        return pd.read_csv(out)
    a = ad.read_h5ad(H5, backed="r")
    pt = a.var["probe_type"].astype(str).to_numpy()
    base, neg = np.where(pt == "Base")[0], np.where(pt == "Negprob")[0]
    obs_names = a.obs_names.astype(str).to_numpy()
    rows = []
    for s in sorted(a.obs["sample_id"].astype(str).unique()):
        z = np.load(cc.CACHE / "mendeley" / f"{s}.npz")  # original-rule cache: QC'd cell ids + author_mapped labels
        lab = pd.Series(z["author_mapped"], index=z["cell_id"])
        idx = np.where((a.obs["sample_id"].astype(str) == s).to_numpy())[0]
        keep = idx[np.isin(obs_names[idx], lab.index)]
        X = a[keep].to_memory().X
        X = (X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)) > 0
        lk = lab.loc[obs_names[keep]].to_numpy()
        for tkey, tcls in TYPES.items():
            m = lk == tcls
            if m.sum() < TYPEAWARE_MIN_CELLS:
                continue
            det = np.asarray(X[m].mean(axis=0)).ravel()
            nmax = float(det[neg].max())
            for j in base:
                rows.append({"section": s, "type": tkey, "n_type_cells": int(m.sum()), "gene": a.var_names[j],
                             "detection_rate_in_type": float(det[j]), "negprobe_max_detection_in_type": nmax,
                             "near_background": bool(det[j] < BG_FOLD * nmax)})
        log.info("type-aware background %s done", s)
        del X
    df = pd.DataFrame(rows)
    cc.save_csv(df, PHASE, "gene_background_typeaware_by_section.csv",
                "Amendment 1: per section x type x gene: detection in own-type cells vs 2x max neg-probe detection in them", log)
    return df


# ----------------------------------------------------------------------------- scoring
def score_sets(X: sp.csr_matrix, var: list[str], sets: dict, section_seed: int):
    """scanpy-score_genes-equivalent with recorded controls. Returns (scores n x P, controls dict)."""
    means = pd.Series(np.asarray(X.mean(axis=0)).ravel(), index=var)
    n_items = int(np.round(len(means) / (N_BINS - 1)))
    cut = (means.rank(method="min") // n_items).astype(int)
    vi = {g: i for i, g in enumerate(var)}
    Xc = X.tocsc()
    scores, controls = {}, {}
    for p_i, (name, genes) in enumerate(sets.items()):
        genes = [g for g in genes if g in vi]
        if not genes:
            scores[name] = np.full(X.shape[0], np.nan, np.float32)
            controls[name] = []
            continue
        rng = np.random.default_rng([cc.SEED, section_seed, p_i])
        ctrl = set()
        for g in genes:
            pool = np.array(sorted(set(cut.index[cut == cut[g]]) - set(genes)))
            if len(pool):
                ctrl.update(rng.choice(pool, size=min(CTRL_SIZE, len(pool)), replace=False).tolist())
        ctrl = sorted(ctrl)
        gi, ci = [vi[g] for g in genes], [vi[g] for g in ctrl]
        s = np.asarray(Xc[:, gi].mean(axis=1)).ravel() - np.asarray(Xc[:, ci].mean(axis=1)).ravel()
        scores[name] = s.astype(np.float32)
        controls[name] = ctrl
    return scores, controls


def sym_knn(xy, k=K) -> sp.csr_matrix:
    from phase3_pdac_audit import knn_edges
    r, c = knn_edges(xy, k)
    n = xy.shape[0]
    A = sp.csr_matrix((np.ones(len(r), np.int8), (r, c)), shape=(n, n))
    A = ((A + A.T) > 0).astype(np.int8).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A


def radius_graph(xy, r=RADIUS) -> sp.csr_matrix:
    A = NearestNeighbors(radius=r).fit(xy).radius_neighbors_graph(xy, mode="connectivity").astype(np.int8).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A


def type_edge_counts(A: sp.csr_matrix, lab: np.ndarray, classes: list[str]) -> pd.DataFrame:
    T = sp.csr_matrix((np.ones(len(lab)), (np.arange(len(lab)), [classes.index(x) for x in lab])),
                      shape=(len(lab), len(classes)))
    M = (T.T @ A.astype(np.float64) @ T).toarray()  # M[b,c] = undirected b-c edges for b != c
    return pd.DataFrame(M, index=classes, columns=classes)


# ----------------------------------------------------------------------------- per section
def run_section(section, used_sets, lineage_sets, s_idx, log):
    sdir = cc.out_dir(f"{PHASE}/sections")
    done = sdir / f"{section}_counts.csv"
    cache = cc.CACHE / CACHE_SUB / f"{section}.npz"
    if done.exists() and cache.exists():
        log.info("section %s already done; skipping (resume)", section)
        return
    import v3_nulls as vn
    adata, meta, labs, _, _ = vn.load_section(section, log)
    X = adata.X.tocsr()
    var = list(adata.var_names)
    xy = adata.obsm["spatial"].astype(np.float64)
    scores, controls = score_sets(X, var, used_sets, s_idx)
    lin, lin_ctrl = score_sets(X, var, lineage_sets, 1000 + s_idx)
    A = sym_knn(xy)
    R = radius_graph(xy)
    total = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    author_major = adata.obs["annotation_majortypes"].astype(str).to_numpy()
    author_sub = adata.obs["annotation_subtypes"].astype(str).to_numpy()
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, cell_id=adata.obs_names.astype(str).to_numpy().astype("U"), xy=xy,
                        author_mapped=labs["author_mapped"].astype(str).to_numpy().astype("U"),
                        marker=labs["marker"].astype(str).to_numpy().astype("U"),
                        author_major=author_major.astype("U"), author_sub=author_sub.astype("U"),
                        total_counts=total.astype(np.float32),
                        score_names=np.array(list(scores), dtype="U"), scores=np.column_stack(list(scores.values())),
                        lineage_names=np.array(list(lin), dtype="U"), lineage=np.column_stack(list(lin.values())),
                        knn_indptr=A.indptr, knn_indices=A.indices, r30_indptr=R.indptr, r30_indices=R.indices)
    # controls
    crow = [{"section": section, "set": k, "kind": "program", "n_controls": len(v), "controls": ";".join(v)}
            for k, v in controls.items()]
    crow += [{"section": section, "set": k, "kind": "lineage", "n_controls": len(v), "controls": ";".join(v)}
             for k, v in lin_ctrl.items()]
    pd.DataFrame(crow).to_csv(sdir / f"{section}_controls.csv", index=False)
    # distributions (own-type cells), counts, edges, validity
    dist_rows, count_rows, edge_rows, auc_rows = [], [], [], []
    for ann in ANNOTS:
        lab = labs[ann].astype(str).to_numpy()
        classes = sorted(set(lab))
        for c, n in zip(*np.unique(lab, return_counts=True)):
            count_rows.append({"section": section, "annotation": ann, "cell_type": c, "n_cells": int(n)})
        E = type_edge_counts(A, lab, classes)
        for b, c in itertools.product(classes, classes):
            if b < c:
                edge_rows.append({"section": section, "annotation": ann, "type_1": b, "type_2": c,
                                  "n_edges_6nn_sym": int(E.loc[b, c])})
        for name, s in scores.items():
            tkey = name.split(":")[0]
            m = lab == TYPES[tkey]
            v = s[m]
            q = np.quantile(v, [.05, .25, .5, .75, .95]) if m.sum() else [np.nan] * 5
            dist_rows.append({"section": section, "annotation": ann, "program": name, "cell_type": TYPES[tkey],
                              "n": int(m.sum()), "mean": float(np.mean(v)) if m.sum() else np.nan,
                              "sd": float(np.std(v)) if m.sum() else np.nan,
                              **{f"q{int(p * 100):02d}": float(x) for p, x in zip([.05, .25, .5, .75, .95], q)}})
    for (maj, a_, b_, pa, pb) in [("Malignant", "CLS", "BSL", "tumour:classical", "tumour:basal_like"),
                                  ("CAF", "myCAF", "iCAF", "CAF:myCAF", "CAF:iCAF")]:
        m = (author_major == maj) & np.isin(author_sub, [a_, b_])
        if pa in scores and pb in scores and m.sum() and len(set(author_sub[m])) == 2:
            y = author_sub[m] == a_
            d = scores[pa][m] - scores[pb][m]
            auc_rows.append({"section": section, "comparison": f"{pa} - {pb} for {a_} vs {b_}",
                             "n_first": int(y.sum()), "n_second": int((~y).sum()), "auc": float(roc_auc_score(y, d))})
        else:
            auc_rows.append({"section": section, "comparison": f"{pa} - {pb} for {a_} vs {b_}", "n_first": 0,
                             "n_second": 0, "auc": np.nan})
    pd.DataFrame(dist_rows).to_csv(sdir / f"{section}_distribution.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(sdir / f"{section}_edges.csv", index=False)
    pd.DataFrame(auc_rows).to_csv(sdir / f"{section}_validity.csv", index=False)
    info = {"section": section, "n_cells": int(adata.n_obs), "n_genes_after_qc": int(adata.n_vars),
            "mean_degree_6nn_sym": float(A.nnz / adata.n_obs), "mean_degree_r30": float(R.nnz / adata.n_obs),
            "median_nn_um": float(np.median(NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]))}
    (sdir / f"{section}_info.json").write_text(json.dumps(info, indent=2))
    pd.DataFrame(count_rows).to_csv(done, index=False)  # written last = completion marker
    for f in ("controls", "distribution", "edges", "validity", "counts"):
        cc.register(sdir / f"{section}_{f}.csv", PHASE, f"Per-section {f} (Phase 1)")
    cc.register(sdir / f"{section}_info.json", PHASE, "Per-section graph/size info (Phase 1)")
    log.info("section %s done: %d cells; %s", section, adata.n_obs, cc.mem_status())
    del adata, X, A, R


def main():
    params = {"bg_fold": BG_FOLD, "bg_min_sections": BG_MIN_SECTIONS, "min_genes": MIN_GENES, "n_bins": N_BINS,
              "ctrl_size": CTRL_SIZE, "k": K, "radius_um": RADIUS, "floor_cells": FLOOR_CELLS,
              "floor_edges": FLOOR_EDGES}
    log = cc.start_run("chain_phase1_cell_states", {**params, "rule": RULE, "phase_dir": PHASE, "cache": CACHE_SUB})
    import v3_nulls as vn
    s2t = vn.section_to_tumour()
    sections = sorted(s2t)
    bg = background_flags(log)
    nb = bg.groupby("gene").near_background.sum()
    bg_genes = set(nb.index[nb >= BG_MIN_SECTIONS])
    gflag = pd.DataFrame({"gene": nb.index, "n_sections_near_background": nb.values,
                          "background_level": nb.values >= BG_MIN_SECTIONS})
    if not (cc.out_dir(PHASE) / "gene_background_flags.csv").exists():
        cc.save_csv(gflag, PHASE, "gene_background_flags.csv",
                    f"Per panel gene: sections near background; background_level if >= {BG_MIN_SECTIONS}/16", log)
    panel = set(nb.index)
    if RULE == "typeaware":
        ta = typeaware_flags(log)
        frac = ta.groupby(["type", "gene"]).near_background.mean()
        bad = {t: set(frac.loc[t].index[frac.loc[t] >= TYPEAWARE_FRAC]) for t in TYPES}
        summ_ta = ta.groupby(["type", "gene"]).agg(n_sections_eligible=("section", "nunique"),
                                                   frac_sections_near_background=("near_background", "mean"),
                                                   median_detection_in_type=("detection_rate_in_type", "median")).reset_index()
        summ_ta["background_level_for_type"] = summ_ta.frac_sections_near_background >= TYPEAWARE_FRAC
        if not (cc.out_dir(PHASE) / "gene_background_typeaware_flags.csv").exists():
            cc.save_csv(summ_ta, PHASE, "gene_background_typeaware_flags.csv",
                        "Amendment 1: per type x gene: fraction of eligible sections near background (>= 0.5 -> flagged)", log)
    else:
        bad = {t: bg_genes for t in TYPES}
    rows, used = [], {}
    for key, genes in PROGRAMS.items():
        bg_genes_t = bad[key[0]]
        u = [g for g in genes if g in panel and g not in bg_genes_t]
        for g in genes:
            rows.append({"program": prog_name(key), "gene": g, "on_panel": g in panel, "background_level": g in bg_genes_t,
                         "used": g in u})
        keep = len(u) >= MIN_GENES or key in SINGLE_GENE_EXEMPT and len(u) == 1
        if keep:
            used[prog_name(key)] = u
        log.info("program %s: %d declared, %d usable -> %s", prog_name(key), len(genes), len(u),
                 "kept" if keep else "DROPPED (<5 usable genes)")
    lin_used = {}
    for t, genes in LINEAGE.items():
        bg_genes_t = bad[t]
        u = [g for g in genes if g in panel and g not in bg_genes_t]
        lin_used[f"lineage:{t}"] = u
        for g in genes:
            rows.append({"program": f"lineage:{t}", "gene": g, "on_panel": g in panel, "background_level": g in bg_genes_t,
                         "used": g in u})
    pg = pd.DataFrame(rows)
    n_used = pg[pg.used].groupby("program").size()
    pg["n_used_in_program"] = pg.program.map(n_used).fillna(0).astype(int)
    pg["program_kept"] = pg.program.isin(list(used) + list(lin_used))
    if not (cc.out_dir(PHASE) / "programs_used.json").exists():
        cc.save_csv(pg, PHASE, "program_genes_used.csv", "Declared program genes: panel, background-level, used, kept", log)
        cc.save_json({"programs": used, "lineage": lin_used}, PHASE, "programs_used.json", "Programs actually scored", log)

    for i, s in enumerate(sections):
        run_section(s, used, lin_used, i, log)

    sdir = cc.out_dir(f"{PHASE}/sections")
    cat = lambda f: pd.concat([pd.read_csv(sdir / f"{s}_{f}.csv") for s in sections], ignore_index=True)  # noqa: E731
    counts, dist, edges, val, ctrl = cat("counts"), cat("distribution"), cat("edges"), cat("validity"), cat("controls")
    for d in (counts, dist, edges, val):
        d.insert(1, "tumour", d.section.map(s2t))
    cc.save_csv(counts, PHASE, "cell_counts_by_section.csv", "Cell counts per type per section (author_mapped, marker)", log)
    cc.save_csv(dist, PHASE, "program_distribution_by_section.csv",
                "Program score distribution among own-type cells per section (scanpy-style score, log1p CP10k)", log)
    cc.save_csv(edges, PHASE, "type_edge_counts_by_section.csv", "Undirected symmetrised-6NN edges per type pair", log)
    cc.save_csv(val, PHASE, "validity_auc.csv", "Sanity: program contrast vs author subtypes (AUC; expected > 0.5)", log)
    cc.save_csv(ctrl, PHASE, "program_controls.csv", "Recorded control genes per section and gene set", log)
    # per-tumour distributions from the cache (pooled cells)
    trows = []
    for tum in sorted(set(s2t.values())):
        secs = [s for s in sections if s2t[s] == tum]
        for ann in ANNOTS:
            pooled = {}
            for s in secs:
                z = np.load(cc.CACHE / CACHE_SUB / f"{s}.npz")
                lab = z[ann]
                for j, name in enumerate(z["score_names"]):
                    pooled.setdefault(name, []).append(z["scores"][lab == TYPES[name.split(":")[0]], j])
            for name, parts in pooled.items():
                v = np.concatenate(parts)
                q = np.quantile(v, [.05, .25, .5, .75, .95]) if v.size else [np.nan] * 5
                trows.append({"tumour": tum, "sections": ";".join(secs), "annotation": ann, "program": name, "n": v.size,
                              "mean": float(v.mean()) if v.size else np.nan, "sd": float(v.std()) if v.size else np.nan,
                              **{f"q{int(p * 100):02d}": float(x) for p, x in zip([.05, .25, .5, .75, .95], q)}})
    tdist = pd.DataFrame(trows)
    cc.save_csv(tdist, PHASE, "program_distribution_by_tumour.csv", "Program score distribution per tumour (pooled)", log)

    # floors
    def n_of(sec, ann, ct):
        r = counts[(counts.section == sec) & (counts.annotation == ann) & (counts.cell_type == ct)]
        return int(r.n_cells.sum())

    def e_of(sec, ann, t1, t2):
        a_, b_ = sorted([t1, t2])
        r = edges[(edges.section == sec) & (edges.annotation == ann) & (edges.type_1 == a_) & (edges.type_2 == b_)]
        return int(r.n_edges_6nn_sym.sum())

    frows = []
    for A_, B_, C_ in itertools.permutations(list(TYPES), 3):
        for ann in ANNOTS:
            for s in sections:
                na, nb_, nc = (n_of(s, ann, TYPES[x]) for x in (A_, B_, C_))
                ebc, eab = e_of(s, ann, TYPES[B_], TYPES[C_]), e_of(s, ann, TYPES[A_], TYPES[B_])
                frows.append({"chain": f"{A_}->{B_}->{C_}", "A": A_, "B": B_, "C": C_, "annotation": ann, "section": s,
                              "tumour": s2t[s], "n_A": na, "n_B": nb_, "n_C": nc, "edges_BC": ebc, "edges_AB": eab,
                              "passes_floor": bool(min(na, nb_, nc) >= FLOOR_CELLS and min(ebc, eab) >= FLOOR_EDGES)})
    fl = pd.DataFrame(frows)
    cc.save_csv(fl, PHASE, "floor_by_section_chain.csv",
                f"Floor (>= {FLOOR_CELLS} cells of A,B,C and >= {FLOOR_EDGES} A-B and B-C edges) per section x chain", log)
    au = fl[fl.annotation == "author_mapped"].copy()
    au["minc"] = au[["n_A", "n_B", "n_C"]].min(axis=1)
    rank = au.groupby("chain").apply(lambda d: pd.Series({
        "n_sections_pass": int(d.passes_floor.sum()), "n_tumours_pass": int(d.loc[d.passes_floor, "tumour"].nunique()),
        "sum_min_count_pass": int(d.loc[d.passes_floor, "minc"].sum()),
        "sections_pass": ";".join(d.loc[d.passes_floor, "section"])}), include_groups=False).reset_index()
    rank = rank.sort_values(["n_sections_pass", "n_tumours_pass", "sum_min_count_pass"], ascending=False)
    primary, reversed_primary = "tumour->CAF->T", "T->CAF->tumour"
    rank["role"] = "not selected"
    rank.loc[rank.chain == primary, "role"] = "primary"
    rank.loc[rank.chain == reversed_primary, "role"] = "negative control (reversed primary)"
    elig = rank[(rank.role == "not selected") & (rank.n_tumours_pass >= 3)]
    rank.loc[elig.index[:5], "role"] = "secondary"
    rank["selection_rank"] = np.arange(1, len(rank) + 1)
    cc.save_csv(rank, PHASE, "chain_selection.csv", "Outcome-free chain selection (counts/edges only; plan section 6)", log)
    log.info("chain selection:\n%s", rank.to_string(index=False))
    val_summary = val.groupby("comparison").auc.agg(["median", "min", "max", "count"]).reset_index()
    log.info("validity AUC:\n%s", val_summary.to_string(index=False))

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1.1, 1]})
    td = tdist[tdist.annotation == "author_mapped"].pivot(index="program", columns="tumour", values="q50")
    tz = td.sub(td.mean(axis=1), axis=0).div(td.std(axis=1).replace(0, np.nan), axis=0)
    im = axes[0].imshow(tz.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    axes[0].set_yticks(range(len(tz)), tz.index, fontsize=8)
    axes[0].set_xticks(range(tz.shape[1]), tz.columns, rotation=60, fontsize=8)
    axes[0].set_title("Median program score among own-type cells per tumour\n(author_mapped; z across tumours)",
                      fontsize=10)
    plt.colorbar(im, ax=axes[0], shrink=.7, label="z (per program)")
    sel = rank[rank.role.isin(["primary", "secondary", "negative control (reversed primary)"])].chain.tolist()
    pm = au[au.chain.isin(sel)].pivot(index="chain", columns="section", values="passes_floor").loc[sel].astype(float)
    axes[1].imshow(pm.to_numpy(), aspect="auto", cmap="Greens", vmin=0, vmax=1.3)
    axes[1].set_yticks(range(len(pm)), [f"{c} ({rank.set_index('chain').role[c][:9]})" for c in pm.index], fontsize=8)
    axes[1].set_xticks(range(pm.shape[1]), pm.columns, rotation=60, fontsize=8)
    axes[1].set_title(f"Floor pass (>= {FLOOR_CELLS} cells/type, >= {FLOOR_EDGES} A-B and B-C edges), author_mapped",
                      fontsize=10)
    for i in range(pm.shape[0]):
        for j in range(pm.shape[1]):
            axes[1].text(j, i, "pass" if pm.iat[i, j] else "", ha="center", va="center", fontsize=6)
    fig.tight_layout()
    cc.save_figure(fig, PHASE, "phase1_programs_and_floors", "Per-tumour program medians and chain floor pass matrix", log)
    summ = {"n_sections": len(sections), "n_tumours": len(set(s2t.values())), "n_background_level_genes": len(bg_genes),
            "programs_kept": {k: len(v) for k, v in used.items()},
            "programs_dropped": [prog_name(k) for k in PROGRAMS if prog_name(k) not in used],
            "chains": rank[rank.role != "not selected"][["chain", "role", "n_sections_pass", "n_tumours_pass"]].to_dict("records"),
            "validity_auc_median": val_summary.set_index("comparison")["median"].to_dict(),
            "median_nn_um_by_section": {s: json.loads((sdir / f"{s}_info.json").read_text())["median_nn_um"] for s in sections}}
    cc.save_json(summ, PHASE, "phase1_summary.json", "Phase 1 key numbers", log)
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    if "--amendment1" in sys.argv:
        PHASE, CACHE_SUB, RULE = "phase1b", "mendeley_v1b", "typeaware"
    main()

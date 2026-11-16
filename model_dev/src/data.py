"""Phase 0: build per-section graph datasets for the development set (Mendeley) and the external set (GSE310352).

What: for each section/slide, cells, positions, labels, group (tumour / unit), total counts, expression over the
model gene set G (log1p CP10k as float16 and raw counts as uint16), program scores with recorded controls, the
symmetrised 6-NN graph and the 30 um radius graph with edge lengths, and 200 um tile codes. One .npz per section
in model_cache/datasets/<dataset>/ (git-ignored), registered with a checksum.

Why: every later phase reads these files, so the gene filter, QC, normalisation and graphs are fixed once here and
cannot drift between baselines, graph models and the external test. Reuses the existing loaders, QC, scoring and
type-aware background flags (nothing outside model_dev/ is modified).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

PHASE = "phase0"


def model_gene_set(cfg: dict) -> pd.DataFrame:
    """Genes passing the type-aware background rule for >= 1 typed class, with per-type usability flags."""
    flags = pd.read_csv(C.REPO_ROOT / cfg["gene_flags"])
    usable = (~flags.background_level_for_type).groupby([flags.gene, flags.type]).any().unstack(fill_value=False)
    usable = usable.reindex(columns=list(cfg["typed_classes"]), fill_value=False)
    usable["n_types_usable"] = usable.sum(axis=1)
    usable = usable[usable.n_types_usable >= 1].sort_index()
    return usable.reset_index()


def edge_arrays(A: sp.csr_matrix, xy: np.ndarray) -> dict[str, np.ndarray]:
    """CSR indptr/indices plus Euclidean edge lengths, for a symmetric adjacency matrix with zero diagonal."""
    A = A.tocsr()
    rows = np.repeat(np.arange(A.shape[0]), np.diff(A.indptr))
    length = np.linalg.norm(xy[rows] - xy[A.indices], axis=1).astype(np.float32)
    return {"indptr": A.indptr.astype(np.int64), "indices": A.indices.astype(np.int64), "length": length}


def build_one(name: str, adata, lab: np.ndarray, group: str, genes: list[str], cfg: dict, s_idx: int,
              scale: float, dataset: str, log) -> dict:
    """Build and save one section dataset. `scale` converts the Mendeley um parameters to this section's units."""
    import chain_phase1_cell_states as p1  # reused: score_sets, sym_knn, radius_graph (same seeds/definitions)
    import chain_features as cf

    programs = json.loads((C.REPO_ROOT / cfg["programs"]).read_text())
    var = list(adata.var_names)
    vi = {g: i for i, g in enumerate(var)}
    gi = [vi[g] for g in genes]
    Xlog = adata.X.tocsc()[:, gi].toarray().astype(np.float16)
    Xcnt = adata.layers["counts"].tocsc()[:, gi].toarray()
    if Xcnt.max() >= 65535:
        raise ValueError("counts exceed uint16")
    Xcnt = Xcnt.astype(np.uint16)
    sets = {**programs["programs"], **programs["lineage"]}
    present = {k: [g for g in v if g in vi] for k, v in sets.items()}
    seed_base = 0 if dataset == "mendeley_pdac_smi" else 500  # same seeds as chain Phase 1b / Phase 6
    scores, controls = p1.score_sets(adata.X.tocsr(), var, present, seed_base + s_idx)
    xy = np.asarray(adata.obsm["spatial"], float)
    W = p1.sym_knn(xy, cfg["knn_k"])
    R = p1.radius_graph(xy, cfg["radius_um"] * scale)
    nn = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]
    out = {
        "cell_id": adata.obs_names.astype(str).to_numpy().astype("U"), "xy": xy.astype(np.float32),
        "label": lab.astype("U"), "group": np.array(group), "section": np.array(name),
        "total_counts": np.asarray(adata.layers["counts"].sum(axis=1)).ravel().astype(np.float32),
        "n_genes_detected": np.asarray((adata.layers["counts"] > 0).sum(axis=1)).ravel().astype(np.int32),
        "genes": np.array(genes, dtype="U"), "X_log": Xlog, "X_counts": Xcnt,
        "score_names": np.array(list(scores), dtype="U"), "scores": np.column_stack(list(scores.values())).astype(np.float32),
        "tile": cf.tile_codes(xy, cfg["tile_um"] * scale).astype(np.int64),
        "scale_units_per_um": np.array(scale, dtype=np.float32), "median_nn": np.array(float(np.median(nn[nn > 0]))),
    }
    for key, A in (("knn", W), ("r30", R)):
        for k, v in edge_arrays(A, xy).items():
            out[f"{key}_{k}"] = v
    d = C.CACHE / "datasets" / dataset
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.npz"
    np.savez_compressed(path, **out)
    C.register_artefact(path, PHASE, f"{dataset} {name} graph dataset (cells, X over G, scores, graphs)")
    ctrl = pd.DataFrame([{"dataset": dataset, "section": name, "set": k, "n_genes_present": len(present[k]),
                          "controls": ";".join(v)} for k, v in controls.items()])
    log.info("%s %s: %d cells, %d genes in G, mean knn degree %.2f, mean r30 degree %.2f, median NN %.2f",
             dataset, name, len(lab), len(genes), W.nnz / len(lab), R.nnz / len(lab), out["median_nn"])
    return {"path": path, "controls": ctrl, "scores": out["scores"], "score_names": out["score_names"],
            "cell_id": out["cell_id"], "label": lab, "n_cells": len(lab), "group": group,
            "mean_knn_degree": W.nnz / len(lab), "mean_r30_degree": R.nnz / len(lab), "median_nn": out["median_nn"]}


def build_development(cfg: dict, genes: list[str], log) -> tuple[list[dict], pd.DataFrame]:
    import v3_nulls as vn
    s2t = vn.section_to_tumour()
    sections = sorted(s2t)
    infos, checks = [], []
    for i, s in enumerate(sections):
        adata, _, labs, _, _ = vn.load_section(s, log)  # QC + normalize_total + log1p; counts layer kept
        lab = labs[cfg["development"]["label_column"]].astype(str).to_numpy()
        info = build_one(s, adata, lab, s2t[s], genes, cfg, i, 1.0, cfg["development"]["dataset"], log)
        info["section"] = s
        # reproduction check against the chain cache (same scoring function and seeds -> identical scores)
        z = np.load(C.REPO_ROOT / cfg["chain_cache_mendeley"] / f"{s}.npz")
        names_chain = [str(x) for x in z["score_names"]]
        assert (z["cell_id"] == info["cell_id"]).all(), f"cell order differs from chain cache for {s}"
        mine = {n: info["scores"][:, j] for j, n in enumerate(info["score_names"])}
        diff = max(float(np.abs(mine[n] - z["scores"][:, j]).max()) for j, n in enumerate(names_chain))
        checks.append({"section": s, "max_abs_score_diff_vs_chain_cache": diff, "n_cells": info["n_cells"]})
        log.info("%s score reproduction vs chain cache: max abs diff %.2e", s, diff)
        infos.append(info)
        del adata
    return infos, pd.DataFrame(checks)


def build_external(cfg: dict, genes: list[str], log) -> list[dict]:
    import pdac_loaders
    import scanpy as sc
    from phase3_pdac_audit import qc
    ext = cfg["external"]
    infos = []
    first = C.RESULTS / PHASE / "external_first_load.json"
    if not first.exists():
        C.save_json({"commit": C.git_commit(), "loaded_at": pd.Timestamp.now().isoformat(),
                     "purpose": "Phase 0 data build and data card only; no model touches it until Phase 4"},
                    PHASE, "external_first_load.json", log)
    for i, slide in enumerate(ext["slides"]):
        adata, _ = pdac_loaders.load(f"{ext['dataset']}:{slide}", C.DATA_DIR, log)
        adata, _, _ = qc(adata, log)
        # restrict to the panel intersection BEFORE normalisation/scoring so preprocessing matches Mendeley
        keep = [g for g in adata.var_names if g in set(genes_all_mendeley(cfg))]
        adata = adata[:, keep].copy()
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        labs = pd.read_csv(C.REPO_ROOT / ext["label_file_pattern"].format(slide=slide),
                           usecols=["cell_id", ext["label_column"]]).set_index("cell_id")
        ids = adata.obs_names.astype(str)
        if not ids.isin(labs.index).all():
            raise ValueError(f"{slide}: cells without v3 marker labels; no silent realignment")
        lab = labs.loc[ids, ext["label_column"]].astype(str).to_numpy()
        xy = np.asarray(adata.obsm["spatial"], float)
        nn = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]
        scale = float(np.median(nn[nn > 0])) / cfg["mendeley_median_nn_um"]  # px per um-equivalent (um/px UNVERIFIED)
        unit = ext["units"].get(slide, slide)
        info = build_one(slide, adata, lab, unit, genes, cfg, i, scale, ext["dataset"], log)
        info["section"] = slide
        infos.append(info)
        del adata
    return infos


def genes_all_mendeley(cfg: dict) -> list[str]:
    import anndata as ad
    a = ad.read_h5ad(C.DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    return list(a.var_names[(a.var["probe_type"].astype(str) == "Base").to_numpy()])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "data.yaml"))
    ap.add_argument("--skip-external", action="store_true")
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "data", cfg)
    gs = model_gene_set(cfg)
    genes = gs.gene.tolist()
    C.save_table(gs, PHASE, "gene_set.csv", log)
    men_all = genes_all_mendeley(cfg)
    gse_cols = pd.read_csv(C.DATA_DIR / "pdac" / "geo_gse310352" / "GSM9294399_slide1_cell_by_gene_counts.csv.gz",
                           nrows=1).columns
    gse_genes = [c for c in gse_cols if c not in ("fov", "cell_ID") and not c.startswith("NegPrb")]
    panel = {"n_mendeley_base": len(men_all), "n_gse310352": len(gse_genes),
             "n_intersection": len(set(men_all) & set(gse_genes)),
             "mendeley_only": sorted(set(men_all) - set(gse_genes)), "gse_only": sorted(set(gse_genes) - set(men_all)),
             "n_model_gene_set_G": len(genes), "G_subset_of_intersection": set(genes) <= (set(men_all) & set(gse_genes)),
             "usable_per_type": {t: int(gs[t].sum()) for t in cfg["typed_classes"]}}
    C.save_json(panel, PHASE, "panel_intersection.json", log)
    log.info("panel: %s", json.dumps(panel))
    dev, checks = build_development(cfg, genes, log)
    C.save_table(checks, PHASE, "score_reproduction_check.csv", log)
    ext = [] if args.skip_external else build_external(cfg, genes, log)
    rows, ctrls, secs = [], [], []
    for ds, infos in (("development", dev), ("external", ext)):
        for info in infos:
            cl, cn = np.unique(info["label"], return_counts=True)
            rows += [{"split": ds, "section": info["section"], "group": info["group"], "cell_type": c, "n_cells": int(n)}
                     for c, n in zip(cl, cn)]
            ctrls.append(info["controls"])
            secs.append({"split": ds, "section": info["section"], "group": info["group"], "n_cells": info["n_cells"],
                         "mean_knn_degree": info["mean_knn_degree"], "mean_r30_degree": info["mean_r30_degree"],
                         "median_nn": info["median_nn"]})
    counts = pd.DataFrame(rows)
    C.save_table(counts, PHASE, "cell_counts.csv", log)
    C.save_table(counts.groupby(["split", "group", "cell_type"]).n_cells.sum().reset_index(), PHASE,
                 "cell_counts_by_group.csv", log)
    C.save_table(pd.concat(ctrls, ignore_index=True), PHASE, "program_controls.csv", log)
    C.save_table(pd.DataFrame(secs), PHASE, "section_info.csv", log)
    log.info("max score diff vs chain cache over sections: %.2e", checks.max_abs_score_diff_vs_chain_cache.max())
    log.info("done")


if __name__ == "__main__":
    main()

"""GSE310352 consistency audit using the Phase 4 plan and resume addendum.

Independence from Mendeley is UNVERIFIED; this is not independent replication.
"""
from __future__ import annotations

import gc
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "data/pdac/numba_cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

import ccc_audit as ca
import pdac_loaders
import v3_nulls as vn
from common import (DATA_DIR, RESULTS_PDAC, finish_run, mem_status, out_dir, register,
                    rel, save_csv, save_figure, save_json, start_run)
from lr_catalog import catalogue
from phase3_pdac_audit import add_flags, annotate, knn_edges, preprocess, qc
from v3_phase5_programs import digest, read_soft

PHASE = "v3_phase4"
NULLS = ["global", "block20nn", "block10nn", "block5nn", "dist_nn"]
ANNOTATIONS = {"marker": "cell_type_marker", "leiden_consensus": "cell_type_leiden"}
PARAMS = {"n_perms": 500, "seed": 42, "k": 6, "blocks_in_nn": [20, 10, 5],
          "distance_bands_in_nn": [0, 2.5, 5, 10, 20, 40], "independence": "UNVERIFIED"}


def permutations(name, strings, types, lab, xy, scale):
    if name == "global":
        # Shuffling encoded labels gives exactly the same index permutation as
        # shuffling their strings with the original RandomState implementation.
        from spatial_null import generate_label_permutations
        return generate_label_permutations(lab, n_perms=500, seed=42), {"scale": "global"}
    if name.startswith("block"):
        multiple = int(name.removeprefix("block").removesuffix("nn"))
        codes = vn.tile_codes(xy, multiple * scale)
        _, count = np.unique(codes, return_counts=True)
        info = {"tile_px": multiple * scale, "tile_nn_multiple": multiple, "n_tiles": len(count),
                "median_cells_per_tile": float(np.median(count)), "n_singleton_tiles": int((count == 1).sum())}
    else:
        if vn.TUMOUR_CLASS not in types:
            return None, {"skipped": "no tumour class"}
        tumour = lab == types.index(vn.TUMOUR_CLASS)
        if min(tumour.sum(), (~tumour).sum()) < 50:
            return None, {"skipped": "fewer than 50 cells in a tumour/non-tumour compartment"}
        distance = vn.signed_boundary_distance(xy, tumour)
        bands = np.array([0, 2.5, 5, 10, 20, 40, np.inf]) * scale
        codes = np.searchsorted(bands, np.abs(distance), side="right") - 1 + 10 * tumour.astype(np.int64)
        _, count = np.unique(codes, return_counts=True)
        info = {"bands_px": bands[:-1].tolist(), "n_tumour": int(tumour.sum()),
                "n_bands_used": len(count), "min_band_cells": int(count.min())}
    return vn._within_group_perms(lab, codes, 500, 42), info


def run_slide(slide, log):
    od = out_dir(f"{PHASE}/slides")
    gsm = pdac_loaders.GSE310352_SLIDES[slide]
    data = DATA_DIR / "pdac/geo_gse310352"
    identity = {"inputs": {kind: digest(data / f"{gsm}_{slide}_{kind}.csv.gz")
                            for kind in ("cell_by_gene_counts", "cell_metadata")}, "params": PARAMS,
                "annotation_method": "phase3_qc_preprocess_annotate_v1"}
    state = od / f"{slide}_identity.json"
    if state.exists() and json.loads(state.read_text()) != identity:
        raise ValueError(f"Resume identity mismatch for {slide}")
    if state.exists() and all((od / f"{slide}__{a}__{n}.csv.gz").exists() for a in ANNOTATIONS for n in NULLS):
        log.info("%s: all 10 results present (verified resume)", slide)
        return
    t0 = time.time()
    adata, meta = pdac_loaders.load_geo_gse310352(data, log, slide)
    adata, before, _ = qc(adata, log)
    xy = np.asarray(adata.obsm["spatial"], float)
    if not np.isfinite(xy).all():
        raise ValueError("Nonfinite spatial coordinates")
    nearest = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]
    positive = nearest[nearest > 0]
    if not len(positive):
        raise ValueError("No positive nearest-neighbour distances")
    scale = float(np.median(positive))
    # Negative-probe diagnostics on the same post-QC cells and raw counts.
    pairs = [(l, r) for l, r, _ in catalogue() if l in adata.var_names and r in adata.var_names]
    genes = sorted({g for pair in pairs for g in pair})
    neg = np.asarray(adata.obsm["negprobes"])
    nmax = float((neg > 0).mean(axis=0).max())
    bg, det = [], {}
    for gene in genes:
        x = adata[:, gene].X.toarray().ravel()
        det[gene] = int((x > 0).sum())
        rate = float((x > 0).mean())
        bg.append({"slide": slide, "gsm": gsm, "gene": gene, "n_cells": adata.n_obs,
                   "detection_rate": rate, "mean_count": float(x.mean()), "negprobe_max_detection_rate": nmax,
                   "fold_over_negprobe_max": rate / nmax if nmax else np.inf, "near_background": rate < 2 * nmax})
    save_csv(pd.DataFrame(bg), f"{PHASE}/slides", f"{slide}_background.csv", "Gene detection vs negative probes (<2x max diagnostic)", log)
    save_csv(pd.DataFrame({"probe": adata.uns["negprobe_names"], "detection_rate": (neg > 0).mean(axis=0),
                           "mean_count": neg.mean(axis=0)}), f"{PHASE}/slides", f"{slide}_negative_probes.csv", "Actual negative-probe detection and mean counts", log)
    preprocess(adata, log)
    used, dropped, _ = annotate(adata, log)
    labs = {a: adata.obs[column].astype(str).to_numpy() for a, column in ANNOTATIONS.items()}
    labdf = pd.DataFrame({"cell_id": adata.obs_names, "fov": adata.obs.fov.to_numpy(),
                          "cell_ID": adata.obs.cell_ID.to_numpy(), "x_px": xy[:, 0], "y_px": xy[:, 1], **labs})
    save_csv(labdf, f"{PHASE}/slides", f"{slide}_labels.csv.gz", "Marker/Leiden labels; original FOV and cell identifiers", log)
    rows, cols = knn_edges(xy, 6)
    S = ca.setup(adata, pairs, rows, cols)
    background = {r["gene"]: r["near_background"] for r in bg}
    info = {"slide": slide, "gsm": gsm, "n_before_qc": before[0], "n_cells": adata.n_obs,
            "n_testable_axes": len(pairs), "n_edges": len(rows), "median_positive_nn_px": scale,
            "n_zero_nn": int((nearest == 0).sum()), "um_per_px": "UNVERIFIED",
            "patient_independence": "UNVERIFIED", "markers_used": used, "types_dropped": dropped,
            "label_counts": {a: pd.Series(v).value_counts().to_dict() for a, v in labs.items()}, "null_info": {}}
    for a, strings in labs.items():
        types, lab = vn.encode(strings)
        for n in NULLS:
            path = od / f"{slide}__{a}__{n}.csv.gz"
            if path.exists() and state.exists():
                continue
            perms, ni = permutations(n, strings, types, lab, xy, scale)
            info["null_info"][f"{a}/{n}"] = ni
            if perms is None:
                df = pd.DataFrame({"section": [slide], "annotation": [a], "null": [n], "skipped": [ni["skipped"]]})
            else:
                t = time.time()
                obs, null, nze = ca.run_null(lab, len(types), S, perms)
                interactions = [(s, r, l, rc) for l, rc in pairs for s in types for r in types]
                df, _ = ca.tidy(interactions, types, S, obs, null, nze)
                df = add_flags(df, det)
                df["background_level_gene_diagnostic"] = df.ligand.map(background) | df.receptor.map(background)
                df["n_receiver_cells"] = pd.Series(strings).value_counts().reindex(df.receiver).to_numpy()
                df["n_sender_cells"] = pd.Series(strings).value_counts().reindex(df.sender).to_numpy()
                df["meets_reference_calibration_floor"] = (df.n_receiver_cells >= 100) & (df.n_qualifying_edges >= 300)
                for key, value in [("null", n), ("annotation", a), ("section", slide)]:
                    df.insert(0, key, value)
                df["coordinate_units"] = "pixels; blocks and distance bands in nearest-neighbour units"
                log.info("%s/%s/%s %.1fs; %s", slide, a, n, time.time() - t, mem_status())
                del obs, null, perms
            df.to_csv(path, index=False)
            register(path, PHASE, f"GSE310352 full catalogue, {slide}/{a}/{n}; 500 add-one permutations")
            gc.collect()
    save_json(info, f"{PHASE}/slides", f"{slide}_info.json", "QC, scale, annotation, compartment and null-construction diagnostics", log)
    save_json(identity, f"{PHASE}/slides", state.name, "Input digests and immutable resume parameters", log)
    log.info("%s complete in %.1fs", slide, time.time() - t0)
    del adata, S
    gc.collect()


def aggregate(log):
    od = out_dir(f"{PHASE}/slides")
    frames, skipped = [], []
    for p in sorted(od.glob("*__*.csv.gz")):
        df = pd.read_csv(p)
        (skipped if "skipped" in df else frames).append(df)
    allrows = pd.concat(frames, ignore_index=True)
    caf_t = allrows[(allrows.sender == "CAF_fibroblast") & (allrows.receiver == "T_cell")].copy()
    headline = caf_t[(caf_t.ligand == "CXCL12") & (caf_t.receptor == "CXCR4")]
    summary = headline.groupby(["annotation", "null"], as_index=False).agg(
        n_slides=("section", "size"), n_evaluable=("ratio_obs", "count"), median_Q=("ratio_obs", "median"),
        median_zQ=("ratio_z", "median"), n_Q_gt1=("ratio_obs", lambda x: int((x > 1).sum())),
        n_p05=("ratio_p_plus1", lambda x: int((x < .05).sum())),
        n_bh05=("ratio_q_bh", lambda x: int((x < .05).sum())), median_T_cells=("n_receiver_cells", "median"),
        median_edges=("n_qualifying_edges", "median"), n_meet_reference_floor=("meets_reference_calibration_floor", "sum"),
        n_background_flagged=("background_level_gene_diagnostic", "sum"))
    axis = caf_t.groupby(["annotation", "null", "ligand", "receptor"], as_index=False).agg(
        n_evaluable=("ratio_obs", "count"), median_Q=("ratio_obs", "median"), median_zQ=("ratio_z", "median"),
        n_Q_gt1=("ratio_obs", lambda x: int((x > 1).sum())), n_bh05=("ratio_q_bh", lambda x: int((x < .05).sum())),
        n_background_flagged=("background_level_gene_diagnostic", "sum"))
    ref_frames = []
    for p in sorted((RESULTS_PDAC / "v3_phase1/sections").glob("*__marker__global.csv.gz")):
        ref = pd.read_csv(p)
        ref_frames.append(ref[(ref.sender == "CAF_fibroblast") & (ref.receiver == "T_cell")])
    ref = pd.concat(ref_frames).groupby(["ligand", "receptor"], as_index=False).ratio_obs.median().rename(columns={"ratio_obs": "mendeley_marker_median_Q"})
    axis = axis.merge(ref, on=["ligand", "receptor"], how="left", validate="many_to_one")
    axis["same_Q_direction_as_mendeley_marker"] = np.where(axis.mendeley_marker_median_Q.notna() & axis.median_Q.notna(),
        np.where(np.sign(axis.median_Q - 1) == np.sign(axis.mendeley_marker_median_Q - 1), "yes", "no"), "UNVERIFIED")
    axis["independence"] = "UNVERIFIED"
    bg = pd.concat([pd.read_csv(p) for p in od.glob("*_background.csv")], ignore_index=True)
    for frame, name, desc in [(headline, "headline_by_null.csv", "Headline CAF to T CXCL12-CXCR4 with calibration/background diagnostics"),
                              (summary, "headline_summary_by_null.csv", "Across-slide headline summaries (not patient-level inference)"),
                              (caf_t, "CAF_T_all_axes_by_null.csv", "Per-axis CAF to T results under every null"),
                              (axis, "axis_consistency_summary.csv", "Direction vs Mendeley marker Q; independence unverified"),
                              (bg, "background_diagnostics.csv", "Raw detection relative to measured negative probes")]:
        save_csv(frame, PHASE, name, desc, log)
    save_csv(pd.concat(skipped, ignore_index=True) if skipped else pd.DataFrame(columns=["section", "annotation", "null", "skipped"]),
             PHASE, "skipped_runs.csv", "Skipped nulls and reasons", log)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, (a, label) in zip(axes, [("marker", "Marker labels"), ("leiden_consensus", "Leiden consensus labels")]):
        h = headline[headline.annotation == a]
        for slide, data in h.groupby("section"):
            values = data.set_index("null").reindex(NULLS).ratio_z
            ax.plot(range(len(NULLS)), values, marker="o", ms=3, lw=.8, alpha=.6, label=slide)
        ax.axhline(0, color="black", ls="--", lw=.7)
        ax.set(xticks=range(len(NULLS)), xticklabels=NULLS, title=label, ylabel="Headline z(Q)")
        ax.tick_params(axis="x", rotation=35)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].legend(fontsize=7, ncol=2)
    fig.suptitle("GSE310352 consistency audit; independence and micrometre scale unverified")
    fig.tight_layout()
    save_figure(fig, PHASE, "phase4_consistency", "Per-slide headline zQ; blocks in multiples of median nearest-neighbour distance", log)
    info = [json.loads(p.read_text()) for p in sorted(od.glob("*_info.json"))]
    save_json({"status": "DONE", "framing": "cross-dataset consistency audit", "independence": "UNVERIFIED",
               "n_slides": len(info), "n_cells": sum(r["n_cells"] for r in info), "n_completed_null_runs": len(frames),
               "n_skipped_null_runs": len(skipped), "parameters": PARAMS, "um_per_px": "UNVERIFIED",
               "script_sha256": digest(__file__), "headline": summary.to_dict(orient="records")},
              PHASE, "phase4_summary.json", "Completion, dataset scope, limitations and headline findings", log)


def main():
    if not (ROOT / "LICENSE_CONFIRMED_GSE310352.txt").exists():
        raise PermissionError("Phase 4 licence confirmation absent")
    log = start_run("v3_phase4_consistency")
    register(log.log_path, PHASE, "Phase 4 run log (registered before computation)")
    log.info("Parameters: %s; %s", PARAMS, mem_status())
    geo = read_soft(DATA_DIR / "GSE310352_family.soft.gz")
    save_csv(geo, PHASE, "geo_sample_metadata.csv", "GEO sample titles and tissue metadata; no patient identities provided", log)
    for slide in pdac_loaders.GSE310352_SLIDES:
        run_slide(slide, log)
    aggregate(log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

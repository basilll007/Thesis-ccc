"""v2 Phase 1 pre-plan descriptives: cell size, distance from each myeloid receiver to its nearest CAF, and protein
channel distributions by cell type. No target (program/gene) values are read.

Why: the Phase 1 plan must choose distance bins from the distance distribution and the protein identity thresholds
from the protein distributions, before any spillover result is computed. Outputs go to
model_dev/results/v2_phase1/preplan/ and are committed with the plan.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import v2common as V  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = "v2_phase1/preplan"


def main() -> None:
    log = C.start_run(OUT, "v2_preplan_descriptives")
    V.memory_gate(log)
    meta = V.cell_metadata(log)
    q = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
    # cell size by type
    size = meta.groupby("label").eq_diam_um.describe(percentiles=q).reset_index()
    C.save_table(size, OUT, "cell_diameter_by_type.csv", log)
    D = float(meta.eq_diam_um.median())
    log.info("median equivalent cell diameter (all development cells) = %.2f um", D)
    # nearest-CAF distance for myeloid receivers (centroid to centroid)
    rows = []
    for f in sorted((C.CACHE / "datasets" / "mendeley_pdac_smi").glob("*.npz")):
        z = np.load(f, allow_pickle=False)
        lab = z["label"].astype(str)
        d, _ = V.nearest_type_distance(z["xy"].astype(float), lab, V.CAF, k=1)
        my = np.where(lab == V.MYELOID)[0]
        rows.append(pd.DataFrame({"section": str(z["section"]), "group": str(z["group"]), "row": my,
                                  "d_nearest_caf_um": d[my, 0]}))
    dist = pd.concat(rows, ignore_index=True)
    dq = dist.d_nearest_caf_um.quantile(q + [0.99]).rename("d_nearest_caf_um").reset_index().rename(columns={"index": "quantile"})
    C.save_table(dq, OUT, "myeloid_nearest_caf_distance_quantiles.csv", log)
    # candidate bins in units of the median cell diameter D (touching <= D < near <= 2D < separated)
    edges = [0, D, 2 * D, np.inf]
    dist["bin"] = pd.cut(dist.d_nearest_caf_um, edges, labels=["touching", "near", "separated"], right=True)
    bins = dist.groupby(["group", "bin"], observed=False).size().unstack(fill_value=0)
    bins.loc["ALL"] = bins.sum()
    bins["n"] = bins.sum(axis=1)
    for b in ["touching", "near", "separated"]:
        bins[f"frac_{b}"] = bins[b] / bins["n"]
    C.save_table(bins.reset_index(), OUT, "myeloid_distance_bins_by_tumour.csv", log)
    hist = np.histogram(np.clip(dist.d_nearest_caf_um, 0, 100), bins=np.arange(0, 101, 2))
    C.save_table(pd.DataFrame({"bin_left_um": hist[1][:-1], "n_myeloid": hist[0]}), OUT,
                 "myeloid_nearest_caf_distance_hist.csv", log)
    # protein channels by type (per section quantiles of CD45 and PanCK)
    prot = (meta.groupby(["section", "label"])[["Mean.CD45", "Mean.PanCK", "Mean.CD3"]]
            .quantile([0.1, 0.5, 0.9]).unstack().reset_index())
    prot.columns = ["_".join(str(c) for c in col if c != "") for col in prot.columns]
    C.save_table(prot, OUT, "protein_quantiles_by_section_type.csv", log)
    allp = meta.groupby("label")[["Mean.CD45", "Mean.PanCK", "Mean.CD3", "area_um2"]].median().reset_index()
    C.save_table(allp, OUT, "protein_median_by_type.csv", log)
    C.save_json({"median_eq_diameter_um": D, "n_cells": int(len(meta)), "n_myeloid": int(len(dist)),
                 "bin_edges_um": [0, D, 2 * D, "inf"],
                 "bin_counts_all": {b: int(bins.loc["ALL", b]) for b in ["touching", "near", "separated"]}},
                OUT, "preplan_summary.json", log)
    log.info("bins ALL: %s", bins.loc["ALL"].to_dict())


if __name__ == "__main__":
    main()

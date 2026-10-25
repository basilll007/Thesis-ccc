"""Phase 3e: low-memory radius sensitivity for the headline L-R pair (CXCL12-CXCR4) only. SENSITIVITY ONLY.

The full-catalogue radius sweep in phase3_pdac_audit.py was stopped by the host for low memory at 200 um.
This script repeats the sweep for one L-R pair (all sender/receiver type combinations), which needs a small
fraction of the memory. It reuses radius_edges() with the same RNG sequence (numpy default_rng(42), consumed
radius by radius), so for 15-150 um the edge nulls must equal the saved *_null_radius*_marker.npz columns;
this is checked and reported.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, SEED, finish_run, out_dir, save_csv, save_figure, save_json, save_npz, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import RADII, radius_edges, run_scores  # noqa: E402

PHASE = "phase3"
TAG = "mendeley_pdac_smi__U5-a"
PAIR = ("CXCL12", "CXCR4")


def main() -> None:
    log = start_run("phase3e_radius_headline")
    a = ad.read_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_processed.h5ad")  # written by phase4a (marker labels, log-norm raw)
    labels = np.asarray(a.obs["cell_type_marker"].astype(str).to_numpy())
    coords = np.asarray(a.obsm["spatial"], float)
    rng = np.random.default_rng(SEED)
    vn = list(a.var_names)
    all_pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
    pidx = all_pairs.index(PAIR)
    n_types = len(set(labels))
    tabs, checks = [], []
    for r in RADII:
        rows, cols, frac = radius_edges(coords, r, rng, log)
        df, nulls = run_scores(a, labels, rows, cols, [PAIR], log, f"r={r}", prod_dtype=np.float32)
        df.insert(0, "radius_um", r)
        df.insert(1, "query_cell_fraction", frac)
        df.insert(2, "n_directed_edges", len(rows))
        tabs.append(df)
        save_npz(PHASE, f"{TAG}_radius{int(r)}_headline_pair_nulls.npz",
                 f"SENSITIVITY ONLY: edge/ratio nulls, CXCL12-CXCR4, radius {r} um", log,
                 edge_null=nulls["edge"].astype(np.float32), ratio_null=nulls["ratio"].astype(np.float32))
        prev = out_dir(PHASE) / f"{TAG}_null_radius{int(r)}_marker.npz"
        if prev.exists():
            old = np.load(prev)["edge_null"]
            sl = old[:, pidx * n_types * n_types:(pidx + 1) * n_types * n_types]
            diff = float(np.nanmax(np.abs(sl - nulls["edge"].astype(np.float32))))
            checks.append({"radius_um": r, "max_abs_diff_vs_stopped_run_null": diff})
            log.info("r=%.0f: max |edge null - stopped-run null| = %.3g", r, diff)
        del rows, cols
    rs = pd.concat(tabs, ignore_index=True)
    save_csv(rs, PHASE, f"{TAG}_radius_sensitivity_headline_pair.csv",
             "SENSITIVITY ONLY: CXCL12-CXCR4 C/E/M/Q across radius 15-300 um (marker labels, 500 perms)", log)
    save_json(checks, PHASE, f"{TAG}_radius_headline_repro_check.json",
              "Edge-null agreement with the saved nulls of the stopped full-catalogue run (15-150 um)", log)
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))
    for rc, c in (("T_cell", "#d62728"), ("epithelial_ductal_tumor", "#1f77b4")):
        d = rs[(rs.sender == "CAF_fibroblast") & (rs.receiver == rc)]
        axes[0].plot(d.radius_um, d.edge_z, "s-", c=c, label=f"CAF→{rc}")
        axes[1].plot(d.radius_um, d.ratio_z, "^-", c=c, label=f"CAF→{rc}")
        axes[2].plot(d.radius_um, d.n_qualifying_edges, "o-", c=c, label=f"CAF→{rc}")
    for ax, t in zip(axes, ["edge z", "ratio z (E/M)", "qualifying edges"]):
        ax.set_xlabel("radius (µm) — SENSITIVITY ONLY")
        ax.set_title(t)
        ax.legend(fontsize=8)
    axes[0].axhline(1.645, color="#aaa", ls=":")
    axes[1].axhline(1.645, color="#aaa", ls=":")
    axes[2].set_yscale("log")
    fig.suptitle("U5-a CXCL12→CXCR4 radius sensitivity (receiver subsampling above 12M edges; see query_cell_fraction)")
    save_figure(fig, PHASE, f"{TAG}_3e_radius_sensitivity_headline", "Headline radius sensitivity (sensitivity only)", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

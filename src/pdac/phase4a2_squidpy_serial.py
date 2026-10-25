"""Phase 4a2: Squidpy ligrec re-run with numba_parallel=False (attempt 2 for the full section).

Attempt 1 (phase4a) on the full U5-a section (74,793 cells) returned NaN for all 3,100 p-values. Squidpy
switches to its numba-parallel kernel when n_cells >= 2**15; the 8,000-cell crop (serial kernel) returned
finite p-values for all rows not masked by the 1% expression threshold. This script re-runs ligrec with
numba_parallel=False on the full section and, as a control, on the crop (must reproduce attempt 1 exactly).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, SEED, finish_run, out_dir, save_csv, save_json, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import squidpy as sq  # noqa: E402

from lr_catalog import catalogue  # noqa: E402

PHASE = "phase4"
TAG = "mendeley_pdac_smi__U5-a"


def run(adata, pairs, log, tag):
    t0 = time.time()
    res = sq.gr.ligrec(adata, cluster_key="cell_type_marker", interactions=pd.DataFrame(pairs, columns=["source", "target"]),
                       n_perms=500, seed=SEED, copy=True, use_raw=True, threshold=0.01, n_jobs=1,
                       numba_parallel=False, show_progress_bar=False)
    means, pvals = res["means"], res["pvalues"]
    rows = [{"ligand": l, "receptor": r, "sender": s, "receiver": rc, "ligrec_mean": float(m[(s, rc)]),
             "ligrec_pvalue": float(pvals.loc[(l, r), (s, rc)])} for (l, r), m in means.iterrows() for (s, rc) in means.columns]
    df = pd.DataFrame(rows)
    log.info("[%s] ligrec serial: %d rows, %d NaN p, %d p<0.05, %.1fs", tag, len(df), int(df.ligrec_pvalue.isna().sum()),
             int((df.ligrec_pvalue < 0.05).sum()), time.time() - t0)
    return df


def main() -> None:
    log = start_run("phase4a2_squidpy_serial")
    full = ad.read_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_processed.h5ad")
    crop = ad.read_h5ad(DATA_DIR / "processed" / f"pdac_{TAG}_benchmark_crop.h5ad")
    vn = list(full.var_names)
    pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
    c2 = run(crop, pairs, log, "crop")
    c1 = pd.read_csv(out_dir(PHASE) / "squidpy_ligrec_crop.csv")
    k = ["ligand", "receptor", "sender", "receiver"]
    m = c1.merge(c2, on=k, suffixes=("_a1", "_a2"))
    same_nan = bool((m.ligrec_pvalue_a1.isna() == m.ligrec_pvalue_a2.isna()).all())
    maxdiff = float(np.nanmax(np.abs(m.ligrec_pvalue_a1 - m.ligrec_pvalue_a2)))
    log.info("crop control: identical NaN pattern=%s, max |p diff|=%.3g", same_nan, maxdiff)
    f2 = run(full, pairs, log, "full")
    save_csv(f2, PHASE, "squidpy_ligrec_full_serial.csv",
             "squidpy.gr.ligrec on full section, numba_parallel=False (attempt 2), 500 perms, seed 42", log)
    a1 = pd.read_csv(out_dir(PHASE) / "squidpy_ligrec_full.csv")
    save_json({"attempt1_full_nan_p": int(a1.ligrec_pvalue.isna().sum()), "attempt1_full_rows": int(len(a1)),
               "attempt2_full_nan_p": int(f2.ligrec_pvalue.isna().sum()), "attempt2_full_p_lt_0.05": int((f2.ligrec_pvalue < 0.05).sum()),
               "crop_control_same_nan_pattern": same_nan, "crop_control_max_abs_p_diff": maxdiff,
               "means_identical_attempt1_vs_2_full": bool(np.allclose(a1.merge(f2, on=k)["ligrec_mean_x"], a1.merge(f2, on=k)["ligrec_mean_y"]))},
              PHASE, "squidpy_ligrec_attempts.json", "Squidpy ligrec attempt 1 (parallel) vs attempt 2 (serial)", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

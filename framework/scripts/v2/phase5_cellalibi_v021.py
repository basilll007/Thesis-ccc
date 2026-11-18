"""v2 Phase 5 (secondary, development set): cellalibi engine v0.2.1 on the exported tool replicates.

Same counts the tools see (F:/Thesis/data/v2_cache/phase5_replicates). Plan: framework/plans/v2_phase4b.md
(secondary analysis). Output: framework/results/v2_phase5/chunks/<section>__<scen>__<level>__cellalibi_v021.csv
"""
from __future__ import annotations

import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import phase4_benchmark as pb  # noqa: E402

REP = Path("F:/Thesis/data/v2_cache/phase5_replicates")
OUT = pb.ROOT / "framework/results/v2_phase5/chunks"


def run(sec, scen, level, reps):
    from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
    from cellalibi.verdicts import audit_claim
    out = OUT / f"{sec}__{scen}__{level}__cellalibi_v021.csv"
    if out.exists():
        return f"skip {out.name}"
    g = np.load(REP / f"{sec}__geometry.npz")
    ctx = SpatialContext(g["xy"], g["labels"], AuditConfig(n_perms=1000, seed=pb.seed_of(sec, "geom"), engine="0.2.1"))
    rows = []
    t0 = time.time()
    for rep in reps:
        s = pb.seed_of(sec, scen, level, rep)
        counts = np.load(REP / f"{sec}__{scen}__{level}__{rep}.npz")["counts"]
        t1 = time.time()
        au = audit_claim(ctx, ExpressionData(counts, list(g["genes"])), cfg=AuditConfig(n_perms=1000, engine="0.2.1"),
                         seed=s + 3, **pb.CLAIM)
        rows.append({"section": sec, "scenario": scen, "level": level, "rep": rep, "tool": "cellalibi_v021",
                     "status": "ok", "seconds": time.time() - t1, "verdict": au.primary.value,
                     "significant": au.primary.value == "SUPPORTED",
                     **{f"check_{c}": v for c, v in au.checks.items()}})
    pd.DataFrame(rows).to_csv(out, index=False)
    return f"done {out.name} in {time.time() - t0:.0f}s"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(REP / "manifest.csv")
    man["level"] = man["level"].astype(str).replace({"nan": "None"})
    tasks = [(sec, scen, lv, sorted(g.rep.tolist())) for (sec, scen, lv), g in man.groupby(["section", "scenario", "level"], sort=False)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(run, *t) for t in tasks]
        for f in as_completed(futs):
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()

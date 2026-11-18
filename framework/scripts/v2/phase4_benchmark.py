"""v2 Phase 4: semi-synthetic ground-truth benchmark (plan: framework/plans/v2_phase4.md).

One task = one (section, scenario, level) setting with N_REPS replicates; output one CSV chunk per task
(resumable: existing chunks are skipped). 4 worker processes, 1 BLAS thread each.
Usage: python phase4_benchmark.py [--reps 100] [--workers 4] [--only SECTION]
"""
from __future__ import annotations

import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"

import argparse
import json
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "framework/src"))
H5AD = "F:/Thesis/data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad"
LABELS = ROOT / "results/pdac/phaseB/sections"
OUT = ROOT / "framework/results/v2_phase4"
CHUNKS = OUT / "chunks"
CLAIM = dict(sender="CAF_fibroblast", receiver="myeloid", ligand="FN1", receptor="ITGB1")
SETTINGS = ([("S0", None)] + [("S1", f) for f in (1.25, 1.5, 2.0)] + [("S2", 2.0)]
            + [("S3", s) for s in (50.0, 200.0, 400.0)] + [("S4", p) for p in (0.1, 0.2, 0.3)]
            + [("S5", n) for n in (25, 50, 100)])
EXPECTED = {"S0": "not SUPPORTED", "S1": "SUPPORTED (or NOT_TESTABLE)", "S2": "EXPLAINED_BY_ABUNDANCE",
            "S3": "EXPLAINED_BY_ZONATION", "S4": "EXPLAINED_BY_SPILLOVER", "S5": "NOT_TESTABLE (or SUPPORTED)"}

_CACHE: dict = {}


def seed_of(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode()) & 0x7FFFFFFF


def load_section(sec: str):
    if sec not in _CACHE:
        from cellalibi.datasets import load_mendeley_section
        from cellalibi.simulate import SemiSyntheticSection
        a = load_mendeley_section(H5AD, sec, LABELS, annotations=("consensus_T",))
        sim = SemiSyntheticSection(np.c_[a.obs.x_um, a.obs.y_um], a.obs["consensus_T"].to_numpy(),
                                   a.layers["counts"], a.var_names, **{k: CLAIM[k] for k in ("ligand", "receptor", "sender", "receiver")})
        _CACHE.clear()
        _CACHE[sec] = sim
    return _CACHE[sec]


def run_task(sec: str, scen: str, level, n_reps: int, chunk_dir: str | None = None) -> str:
    from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
    from cellalibi.verdicts import audit_claim
    tag = f"{sec}__{scen}__{level}"
    out = (Path(chunk_dir) if chunk_dir else CHUNKS) / f"{tag}.csv"
    if out.exists():
        return f"skip {tag}"
    t0 = time.time()
    sim = load_section(sec)
    cfg = AuditConfig(n_perms=1000, engine="0.2.0")
    base_ctx = None
    rows = []
    for rep in range(n_reps):
        s = seed_of(sec, scen, level, rep)
        if scen == "S5":
            rng = np.random.default_rng(s)
            rec = np.flatnonzero(sim.r_mask)
            keep_r = rng.choice(rec, size=int(level), replace=False)
            keep = ~sim.r_mask
            keep[keep_r] = True
            ss = sim.subset(keep)
            counts = ss.draw("S1", 2.0, s + 1)
            ctx = SpatialContext(ss.xy, ss.labels, AuditConfig(n_perms=1000, seed=s + 2, engine="0.2.0"))
        else:
            if base_ctx is None:
                base_ctx = SpatialContext(sim.xy, sim.labels, AuditConfig(n_perms=1000, seed=seed_of(sec, "geom"), engine="0.2.0"))
            counts = sim.draw(scen, level, s + 1)
            ctx = base_ctx
        ex = ExpressionData(counts, sim.genes)
        t1 = time.time()
        au = audit_claim(ctx, ex, cfg=cfg, seed=s + 3, **CLAIM)
        m = {k: v for k, v in au.metrics.items() if not isinstance(v, (list, dict))}
        rows.append({"section": sec, "scenario": scen, "level": level, "rep": rep, "seed": s,
                     **{f"check_{c}": v for c, v in au.checks.items()}, "verdict": au.primary.value,
                     "failing_checks": ";".join(au.failing_checks()), "audit_seconds": time.time() - t1, **m})
    pd.DataFrame(rows).to_csv(out, index=False)
    return f"done {tag} in {time.time() - t0:.0f}s"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--settings", nargs="*", default=None, help="e.g. S0:None S4:0.3")
    ap.add_argument("--outdir", default=None, help="chunk directory (default results/v2_phase4/chunks)")
    args = ap.parse_args()
    chunk_dir = Path(args.outdir) if args.outdir else CHUNKS
    chunk_dir.mkdir(parents=True, exist_ok=True)
    sel = json.loads((OUT / "section_selection.json").read_text())
    secs = args.sections or sel["selected"]
    settings = SETTINGS
    if args.settings:
        settings = [st for st in SETTINGS if f"{st[0]}:{st[1]}" in args.settings]
    tasks = [(sec, sc, lv) for sec in secs for sc, lv in settings]
    # group by section so each worker keeps one section cached
    tasks.sort(key=lambda t: (t[0], SETTINGS.index((t[1], t[2]))))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_task, *t, args.reps, str(chunk_dir)): t for t in tasks}
        for f in as_completed(futs):
            try:
                print(f.result(), flush=True)
            except Exception as e:  # log and continue
                print(f"ERROR {futs[f]}: {e!r}", flush=True)


if __name__ == "__main__":
    main()

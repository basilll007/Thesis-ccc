"""v2 Phase 4b: fresh-seed benchmark of engine v0.2.1 on two claims (plan: framework/plans/v2_phase4b.md).

One task = (claim, section, setting) with N replicates; one resumable chunk per task. 4 workers x 1 BLAS thread.
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
OUT = ROOT / "framework/results/v2_phase4b"
CLAIMS = {
    "C1": dict(sender="CAF_fibroblast", receiver="myeloid", ligand="FN1", receptor="ITGB1"),
    "C2": dict(sender="myeloid", receiver="epithelial_ductal_tumor", ligand="SPP1", receptor="CD44"),
}
SETTINGS = ([("S0", None)] + [("S1", f) for f in (1.25, 1.5, 2.0)] + [("S2", 2.0)]
            + [("S3", s) for s in (50.0, 200.0, 400.0)] + [("S4", p) for p in (0.1, 0.2, 0.3)]
            + [("S5", n) for n in (25, 50, 100)])
_CACHE: dict = {}


def seed_of(*parts) -> int:
    return zlib.crc32("|".join(map(str, ("v021",) + parts)).encode()) & 0x7FFFFFFF


def load_section(claim: str, sec: str):
    key = (claim, sec)
    if key not in _CACHE:
        from cellalibi.datasets import load_mendeley_section
        from cellalibi.simulate import SemiSyntheticSection
        a = load_mendeley_section(H5AD, sec, LABELS, annotations=("consensus_T",))
        sim = SemiSyntheticSection(np.c_[a.obs.x_um, a.obs.y_um], a.obs["consensus_T"].to_numpy(),
                                   a.layers["counts"], a.var_names, **CLAIMS[claim])
        _CACHE.clear()
        _CACHE[key] = sim
    return _CACHE[key]


def select_sections(claim: str) -> list[str]:
    c = CLAIMS[claim]
    rows = []
    for f in sorted(LABELS.glob("*_consensus_labels.csv.gz")):
        lab = pd.read_csv(f, usecols=["consensus_T"])["consensus_T"].astype(str)
        rows.append({"section": f.name.replace("_consensus_labels.csv.gz", ""), "n_cells": len(lab),
                     "n_sender": int((lab == c["sender"]).sum()), "n_receiver": int((lab == c["receiver"]).sum())})
    df = pd.DataFrame(rows).sort_values("n_cells").reset_index(drop=True)
    el = df[(df.n_sender >= 100) & (df.n_receiver >= 100)].reset_index(drop=True)
    k = len(el)
    return [el.section.iloc[0], el.section.iloc[(k - 1) // 2], el.section.iloc[k - 1]]


def run_task(claim: str, sec: str, scen: str, level, n_reps: int) -> str:
    from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
    from cellalibi.verdicts import audit_claim
    tag = f"{claim}__{sec}__{scen}__{level}"
    out = OUT / "chunks" / f"{tag}.csv"
    if out.exists():
        return f"skip {tag}"
    t0 = time.time()
    sim = load_section(claim, sec)
    cfg = AuditConfig(n_perms=1000, engine="0.2.1")
    base_ctx = None
    rows = []
    for rep in range(n_reps):
        s = seed_of(claim, sec, scen, level, rep)
        if scen == "S5":
            rng = np.random.default_rng(s)
            rec = np.flatnonzero(sim.r_mask)
            keep_r = rng.choice(rec, size=int(level), replace=False)
            keep = ~sim.r_mask
            keep[keep_r] = True
            ss = sim.subset(keep)
            counts = ss.draw("S1", 2.0, s + 1)
            ctx = SpatialContext(ss.xy, ss.labels, AuditConfig(n_perms=1000, seed=s + 2, engine="0.2.1"))
        else:
            if base_ctx is None:
                base_ctx = SpatialContext(sim.xy, sim.labels, AuditConfig(n_perms=1000, seed=seed_of(claim, sec, "geom"), engine="0.2.1"))
            counts = sim.draw(scen, level, s + 1)
            ctx = base_ctx
        ex = ExpressionData(counts, sim.genes)
        t1 = time.time()
        au = audit_claim(ctx, ex, cfg=cfg, seed=s + 3, **CLAIMS[claim])
        m = {k: v for k, v in au.metrics.items() if not isinstance(v, (list, dict))}
        rows.append({"claim": claim, "section": sec, "scenario": scen, "level": level, "rep": rep, "seed": s,
                     **{f"check_{c}": v for c, v in au.checks.items()}, "verdict": au.primary.value,
                     "failing_checks": ";".join(au.failing_checks()), "audit_seconds": time.time() - t1, **m})
    pd.DataFrame(rows).to_csv(out, index=False)
    return f"done {tag} in {time.time() - t0:.0f}s"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--claims", nargs="*", default=["C1", "C2"])
    args = ap.parse_args()
    (OUT / "chunks").mkdir(parents=True, exist_ok=True)
    sel = {c: select_sections(c) for c in args.claims}
    (OUT / "section_selection.json").write_text(json.dumps(sel, indent=2))
    print("sections", sel, flush=True)
    tasks = [(c, sec, sc, lv) for c in args.claims for sec in sel[c] for sc, lv in SETTINGS]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_task, *t, args.reps): t for t in tasks}
        for f in as_completed(futs):
            try:
                print(f.result(), flush=True)
            except Exception as e:
                print(f"ERROR {futs[f]}: {e!r}", flush=True)


if __name__ == "__main__":
    main()

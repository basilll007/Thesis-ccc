"""v2 Phase 6: cellalibi verdicts on all 16 Mendeley sections and 8 GSE310352 slides.

Plan: framework/plans/v2_phase6.md. One resumable chunk per section; 4 workers x 1 BLAS thread.
Usage: python phase6_real_verdicts.py [--workers 4] [--sections ...]
"""
from __future__ import annotations

import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"

import argparse
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "framework/src"))
H5AD = "F:/Thesis/data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad"
LABELS = ROOT / "results/pdac/phaseB/sections"
GEO = "F:/Thesis/data/pdac/geo_gse310352"
GEO_CACHE = "F:/Thesis/data/chain_cache/gse310352"
OUT = ROOT / "framework/results/v2_phase6"
DEVICE = os.environ.get("CELLALIBI_DEVICE", "cuda")  # amendment 2: GPU per user request
FOCUS = ["epithelial_ductal_tumor", "CAF_fibroblast", "T_cell", "myeloid"]

# thesis catalogue (src/pdac/lr_catalog.py + src/pipeline.py::LR_PAIRS), copied verbatim
CATALOGUE = [("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"), ("EGF", "EGFR"), ("FGF2", "FGFR1"), ("VEGFA", "KDR"),
             ("CXCL12", "CXCR4"), ("CCL2", "CCR2"), ("CCL5", "CCR5"), ("CXCL9", "CXCR3"), ("CXCL10", "CXCR3"),
             ("TNF", "TNFRSF1A"), ("IFNG", "IFNGR1"), ("IL6", "IL6R"), ("IL1B", "IL1R1"), ("HGF", "MET"),
             ("PDGFB", "PDGFRA"), ("PDGFB", "PDGFRB"), ("JAG1", "NOTCH1"), ("CD274", "PDCD1"), ("MIF", "CD74"),
             ("CD86", "CTLA4"), ("PTN", "SDC4"),
             ("CD80", "CTLA4"), ("CD80", "CD28"), ("CD86", "CD28"), ("PDCD1LG2", "PDCD1"), ("LGALS9", "HAVCR2"),
             ("SPP1", "CD44"), ("FN1", "ITGB1"), ("CXCL13", "CXCR5"), ("CCL19", "CCR7"), ("CCL21", "CCR7"),
             ("CXCL12", "ACKR3")]
# thesis chain programs (src/chain/chain_phase1_cell_states.py::PROGRAMS), copied verbatim for the ones used
PROGRAMS = {
    "classical": ["CEACAM6", "AGR2", "S100P", "KRT20", "LCN2", "PIGR", "OLFM4", "SPINK1", "CEACAM1"],
    "basal_like": ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17", "S100A2", "SERPINB5", "AREG", "SLC2A1", "KRT15", "LY6D",
                   "ITGA3", "MET"],
    "C1Q_TAM": ["C1QA", "C1QB", "C1QC", "CD163", "MRC1", "SELENOP"],
    "SPP1_TAM": ["SPP1", "GPNMB", "CD9", "FABP5", "LGALS3"],
}
PROGRAM_CLAIMS = [("CAF_fibroblast", "myeloid", "C1Q_TAM"), ("CAF_fibroblast", "myeloid", "SPP1_TAM"),
                  ("CAF_fibroblast", "epithelial_ductal_tumor", "classical"),
                  ("CAF_fibroblast", "epithelial_ductal_tumor", "basal_like"),
                  ("myeloid", "epithelial_ductal_tumor", "classical"), ("myeloid", "epithelial_ductal_tumor", "basal_like"),
                  ("T_cell", "epithelial_ductal_tumor", "classical"), ("T_cell", "epithelial_ductal_tumor", "basal_like")]


def run_section(dataset: str, sec: str) -> str:
    import psutil
    from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
    from cellalibi.datasets import load_gse310352_slide, load_mendeley_section, mendeley_section_to_tumour
    from cellalibi.verdicts import audit_claim
    out = OUT / "chunks" / f"{dataset}__{sec}.csv"
    if out.exists():
        return f"skip {out.name}"
    proc = psutil.Process(os.getpid())
    t0 = time.time()
    if dataset == "mendeley":
        a = load_mendeley_section(H5AD, sec, LABELS, annotations=("consensus_T",))
        labels = a.obs["consensus_T"].to_numpy()
        unit = mendeley_section_to_tumour(ROOT / "results/pdac/phase3/cross_section_section_info.csv")[sec]
    else:
        a = load_gse310352_slide(GEO, sec, GEO_CACHE)
        labels = a.obs["marker"].to_numpy()
        unit = sec
    xy = np.c_[a.obs.x_um.to_numpy(), a.obs.y_um.to_numpy()]
    cfg = AuditConfig(n_perms=1000, engine="0.2.1", device=DEVICE)
    ctx = SpatialContext(xy, labels, cfg)
    if dataset == "gse310352":  # px coordinates: tile = 10 x median NN distance (plan)
        cfg = AuditConfig(n_perms=1000, tile_um=10.0 * ctx.diameter, engine="0.2.1", device=DEVICE)
        ctx = SpatialContext(xy, labels, cfg)
    cfg0 = AuditConfig(n_perms=1000, tile_um=cfg.tile_um, engine="0.2.0", device=DEVICE)  # amendment 1: secondary engine
    ex = ExpressionData(a.layers["counts"], a.var_names)
    vn = set(a.var_names)
    # program pseudo-features (linear CP10K means) and the "any contact" pseudo-ligand
    for p, genes in PROGRAMS.items():
        g = [x for x in genes if x in vn]
        if g:
            lin = np.asarray(ex.lin[:, [ex.gi[x] for x in g]].mean(axis=1)).ravel()
            ex._cache[f"PROGRAM:{p}"] = (lin, np.log1p(lin), lin > 0)
    one = np.ones(a.n_obs)
    ex._cache["ANY_CONTACT"] = (one, np.log1p(one), one > 0)
    t_load = time.time() - t0
    rows = []
    claims = [(s, r, l, rc, "catalogue") for (l, rc) in CATALOGUE if l in vn and rc in vn
              for s in FOCUS for r in FOCUS if s != r]
    claims += [(s, r, "ANY_CONTACT", f"PROGRAM:{p}", "program") for (s, r, p) in PROGRAM_CLAIMS
               if f"PROGRAM:{p}" in ex._cache]
    headline = {("CAF_fibroblast", "T_cell", "CXCL12", "CXCR4"), ("CAF_fibroblast", "epithelial_ductal_tumor", "CXCL12", "CXCR4"),
                ("CAF_fibroblast", "myeloid", "FN1", "ITGB1"), ("myeloid", "epithelial_ductal_tumor", "SPP1", "CD44")}
    for (s, r, l, rc, kind) in claims:
        engines = [("0.2.1", cfg)] + ([("0.2.0", cfg0)] if (kind == "program" or (s, r, l, rc) in headline) else [])
        for eng, c in engines:
            t1 = time.time()
            try:
                au = audit_claim(ctx, ex, s, r, l, rc, cfg=c)
                m = {k: v for k, v in au.metrics.items() if not isinstance(v, (list, dict))}
                rows.append({"dataset": dataset, "section": sec, "unit": unit, "kind": kind, "engine_run": eng,
                             **{f"check_{cc}": v for cc, v in au.checks.items()}, "verdict": au.primary.value,
                             "failing_checks": ";".join(au.failing_checks()), "seconds": time.time() - t1, **m})
            except Exception as e:
                rows.append({"dataset": dataset, "section": sec, "unit": unit, "kind": kind, "engine_run": eng,
                             "sender": s, "receiver": r, "ligand": l, "receptor": rc, "verdict": "ERROR",
                             "error": f"{type(e).__name__}: {e}"[:300]})
    df = pd.DataFrame(rows)
    df["section_load_seconds"] = t_load
    df["section_total_seconds"] = time.time() - t0
    df["section_peak_rss_mb"] = proc.memory_info().peak_wset / 1e6 if hasattr(proc.memory_info(), "peak_wset") else proc.memory_info().rss / 1e6
    df["n_cells"] = a.n_obs
    df["tile_um_used"] = cfg.tile_um
    df["device"] = DEVICE
    df.to_csv(out, index=False)
    return f"done {out.name}: {len(df)} claims in {time.time() - t0:.0f}s"


def main():
    from cellalibi.datasets import GSE310352_SLIDES, MENDELEY_SECTIONS
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--sections", nargs="*", default=None)
    args = ap.parse_args()
    (OUT / "chunks").mkdir(parents=True, exist_ok=True)
    tasks = [("mendeley", s) for s in MENDELEY_SECTIONS] + [("gse310352", s) for s in GSE310352_SLIDES]
    if args.sections:
        tasks = [t for t in tasks if t[1] in args.sections]
    # max_tasks_per_child=1: a fresh process per section so peak RSS is per section
    with ProcessPoolExecutor(max_workers=args.workers, max_tasks_per_child=1) as pool:
        futs = {pool.submit(run_section, *t): t for t in tasks}
        for f in as_completed(futs):
            try:
                print(f.result(), flush=True)
            except Exception:
                print("ERROR", futs[f], traceback.format_exc()[-600:], flush=True)


if __name__ == "__main__":
    main()

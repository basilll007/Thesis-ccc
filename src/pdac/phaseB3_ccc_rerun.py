"""Phase B3 (v2): re-run the 6-NN C/E/M/Q audit (500 perms, seed 42, full testable catalogue, all sender x receiver)
on all 16 PDAC-SMI sections under the v2 annotations:
  protein_gmm (pre-registered primary; FAILED validation), protein_mad (pre-registered sensitivity; FAILED validation),
  consensus_T (POST-HOC: author T AND marker T).
Same pipeline as v1 (src/pdac/phase3d_author_annotation.py): loader -> QC -> normalise_total(1e4) -> log1p -> raw;
knn_edges(k=6); run_scores (spatial_null.generate_label_permutations, default_rng(42)).
Memory: one section x one annotation at a time; each result is written immediately to
results/pdac/phaseB/ccc_sections/<section>__<annotation>.csv.gz and re-used on restart.
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, RESULTS_PDAC, finish_run, mem_status, out_dir, register, save_csv, start_run  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

import ccc_audit as ca  # noqa: E402
import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import add_flags, knn_edges, qc, run_scores  # noqa: E402

PHASE = "phaseB"
ANNOTATIONS = ["protein_gmm", "protein_mad", "consensus_T"]
CCC_DIR = "phaseB/ccc_sections"


def main() -> None:
    only = sys.argv[1:]  # optional subset of sections
    log = start_run("phaseB3_ccc_rerun")
    log.info("annotations=%s; %s", ANNOTATIONS, mem_status())
    full = ad.read_h5ad(DATA_DIR / "pdac" / "mendeley_pdac_smi" / "raw_meta_data_final.h5ad", backed="r")
    sections = sorted(full.obs["sample_id"].astype(str).unique())
    del full
    if only:
        sections = [s for s in sections if s in only]
    lab_dir = RESULTS_PDAC / "phaseB" / "sections"
    od = out_dir(CCC_DIR)
    for s in sections:
        todo = [a for a in ANNOTATIONS if not (od / f"{s}__{a}.csv.gz").exists()]
        if not todo:
            log.info("section %s: all annotations already saved (resume)", s)
            continue
        t0 = time.time()
        adata, meta = pdac_loaders.load(f"mendeley_pdac_smi:{s}", DATA_DIR, log)
        adata, _, _ = qc(adata, log)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        labs = pd.read_csv(lab_dir / f"{s}_labels.csv.gz", usecols=["cell_id", "protein_gmm", "protein_mad"])
        cons = pd.read_csv(lab_dir / f"{s}_consensus_labels.csv.gz", usecols=["cell_id", "consensus_T"])
        if not (labs.cell_id.astype(str).to_numpy() == adata.obs_names.astype(str).to_numpy()).all() or \
                not (cons.cell_id.astype(str).to_numpy() == adata.obs_names.astype(str).to_numpy()).all():
            raise ValueError(f"cell order mismatch in section {s}")
        labs["consensus_T"] = cons["consensus_T"].to_numpy()
        vn = list(adata.var_names)
        pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
        det = {g: int((adata.layers["counts"][:, vn.index(g)].toarray().ravel() > 0).sum()) for g in {x for p in pairs for x in p}}
        rows, cols = knn_edges(np.asarray(adata.obsm["spatial"], float), 6)
        for a in todo:
            labels = labs[a].astype(str).to_numpy()
            df, _ = run_scores(adata, labels, rows, cols, pairs, log, f"{s}/{a}")
            df = add_flags(df, det)
            df.insert(0, "annotation", a)
            df.insert(0, "section", s)
            df.insert(1, "treatment_status", meta["treatment_status"])
            p = od / f"{s}__{a}.csv.gz"
            df.to_csv(p, index=False)
            register(p, PHASE, f"6-NN C/E/M/Q (500 perms) section {s}, annotation {a}")
            log.info("section %s / %s saved (%d rows; %d types); %s", s, a, len(df), len(set(labels)), mem_status())
            del df
            gc.collect()
        del adata
        gc.collect()
        log.info("section %s done in %.1fs; %s", s, time.time() - t0, mem_status())

    # aggregate (all sections present)
    files = sorted(od.glob("*.csv.gz"))
    allt = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    for key in ca.SCORES:
        allt[f"{key}_q_bh_within_section_annotation"] = allt.groupby(["section", "annotation"])[f"{key}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    allt = allt.drop(columns=[f"{k}_q_bh" for k in ca.SCORES])
    p = out_dir(PHASE) / "v2_ccc_6NN_all_rows.csv.gz"
    allt.to_csv(p, index=False)
    register(p, PHASE, "All sections x v2 annotations: C/E/M/Q obs, z, p (+1 and legacy), BH within section/annotation")
    log.info("aggregate: %d rows from %d files; sections=%d", len(allt), len(files), allt.section.nunique())
    hl = allt[(allt.ligand == "CXCL12") & (allt.receptor == "CXCR4") & (allt.sender == "CAF_fibroblast") & (allt.receiver == "T_cell")]
    save_csv(hl, PHASE, "v2_headline_CAF_T_rows.csv", "CXCL12-CXCR4 CAF->T_cell rows per section under v2 annotations", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

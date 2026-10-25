"""Phase C (v2): independent replication on GEO GSE310352 (CosMx, 8 slides, treatment-naive PDAC). GATED.

STATUS 2026-10-25: BLOCKED. The data licence could not be confirmed (GEO terms page behind reCAPTCHA; publisher page
403; see docs/PDAC_PROGRESS_LOG.md). The script refuses to download or analyse anything unless it is called with
--licence-confirmed. That flag is to be set only after the user confirms the terms of use.

What it does once unblocked (per slide, one at a time; every output written immediately; resumable):
  1. download the minimal file set (8 x cell_by_gene_counts.csv.gz + 8 x cell_metadata.csv.gz, 124,527,604 B per
     filelist.txt) from ftp.ncbi.nlm.nih.gov over HTTPS into data/pdac/geo_gse310352/ (git-ignored), with sha256
  2. load (pdac_loaders.load_geo_gse310352), QC (breast thresholds), normalise/log1p
  3. annotations: v1 marker scoring; pre-registered protein rule (docs/PDAC_V2_PROTEIN_RULE.md; failed validation on
     Mendeley, so reported as unvalidated); optional --labels CSV (columns slide, fov, cell_ID, label) e.g. the
     authors' cell types exported from the seurat.RDS files (needs R, not available on this host)
  4. 6-NN C/E/M/Q audit, 500 perms, seed 42, full testable catalogue; BH within slide x annotation
  5. headline CXCL12->CXCR4 CAF_fibroblast->T_cell summary
Usage: python src/pdac/phaseC_gse310352_audit.py --licence-confirmed [--labels path.csv] [--data-dir DIR] [--no-download]
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, SEED, finish_run, mem_status, out_dir, register, save_csv, save_json, start_run  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

import ccc_audit as ca  # noqa: E402
import pdac_loaders  # noqa: E402
from lr_catalog import catalogue  # noqa: E402
from phase3_pdac_audit import add_flags, knn_edges, qc, run_scores  # noqa: E402
from phaseB1_protein_labels import (ASHMAN_MIN, CHANNELS, MIN_WEIGHT, gate_labels, gmm_threshold,  # noqa: E402
                                    mad_threshold, rna_scores)

PHASE = "phaseC"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9294nnn/{gsm}/suppl/{gsm}_{slide}_{kind}.csv.gz"


def download(d: Path, log):
    d.mkdir(parents=True, exist_ok=True)
    rec = []
    for slide, gsm in pdac_loaders.GSE310352_SLIDES.items():
        for kind in ("cell_by_gene_counts", "cell_metadata"):
            url = BASE.format(gsm=gsm, slide=slide, kind=kind)
            f = d / url.rsplit("/", 1)[1]
            if not f.exists():
                log.info("downloading %s", url)
                tmp = f.with_suffix(".part")
                with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as o:  # TLS verified (default)
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        o.write(b)
                tmp.rename(f)
            h = hashlib.sha256(f.read_bytes()).hexdigest()
            rec.append({"slide": slide, "gsm": gsm, "file": f.name, "url": url, "bytes": f.stat().st_size, "sha256": h})
    return pd.DataFrame(rec)


def protein_labels(adata, sf):
    rng = np.random.default_rng(SEED)
    pos, diag = {}, []
    for c in CHANNELS:
        x = np.log1p(adata.obs[f"Mean.{c}"].to_numpy(dtype=float))
        g = gmm_threshold(x, rng)
        tm = mad_threshold(x)
        ok = g["ashman_d"] >= ASHMAN_MIN and min(g["gmm_w_neg"], g["gmm_w_pos"]) >= MIN_WEIGHT and np.isfinite(g["gmm_threshold"])
        t = g["gmm_threshold"] if ok else tm
        pos[c] = x > t
        diag.append({"channel": c, **g, "gmm_accepted": bool(ok), "mad_threshold": tm, "threshold_primary": t,
                     "frac_pos": float(pos[c].mean())})
    return gate_labels(pos, sf), pd.DataFrame(diag)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--licence-confirmed", action="store_true")
    ap.add_argument("--labels", default=None)
    ap.add_argument("--data-dir", default=str(DATA_DIR / "pdac" / "geo_gse310352"))
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out-phase", default=PHASE)
    args = ap.parse_args()
    if not args.licence_confirmed:
        print("BLOCKED: GSE310352 licence/terms not confirmed. Re-run with --licence-confirmed only after the user "
              "confirms. Nothing was downloaded or analysed.")
        sys.exit(2)
    phase = args.out_phase
    log = start_run("phaseC_gse310352_audit")
    d = Path(args.data_dir)
    if not args.no_download:
        files = download(d, log)
        save_csv(files, phase, "gse310352_files.csv", "Downloaded GSE310352 files: URL, bytes, sha256", log)
    user_labels = pd.read_csv(args.labels) if args.labels else None
    od = out_dir(f"{phase}/slides")
    for slide in pdac_loaders.GSE310352_SLIDES:
        anns = ["marker", "protein_gmm"] + (["user_labels"] if user_labels is not None else [])
        if all((od / f"{slide}__{a}.csv.gz").exists() for a in anns):
            log.info("slide %s: already done (resume)", slide)
            continue
        adata, meta = pdac_loaders.load_geo_gse310352(d, log, slide)
        adata, _, _ = qc(adata, log)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sf = rna_scores(adata)
        labs = {"marker": sf.idxmax(axis=1).to_numpy()}
        labs["protein_gmm"], diag = protein_labels(adata, sf)
        diag.insert(0, "slide", slide)
        diag.to_csv(od / f"{slide}_protein_thresholds.csv", index=False)
        register(od / f"{slide}_protein_thresholds.csv", phase, f"Protein threshold diagnostics, GSE310352 {slide}")
        if user_labels is not None:
            u = user_labels[user_labels.slide == slide].set_index(["fov", "cell_ID"])["label"]
            idx = pd.MultiIndex.from_arrays([adata.obs["fov"], adata.obs["cell_ID"]])
            labs["user_labels"] = u.reindex(idx).fillna("unlabelled").astype(str).to_numpy()
        lab_df = pd.DataFrame({"cell": adata.obs_names, **labs})
        lab_df.to_csv(od / f"{slide}_labels.csv.gz", index=False)
        register(od / f"{slide}_labels.csv.gz", phase, f"Per-cell labels (marker, protein_gmm[, user]) GSE310352 {slide}")
        vn = list(adata.var_names)
        pairs = [(l, r) for l, r, _ in catalogue() if l in vn and r in vn]
        det = {g: int((adata.layers["counts"][:, vn.index(g)].toarray().ravel() > 0).sum()) for g in {x for p in pairs for x in p}}
        rows, cols = knn_edges(np.asarray(adata.obsm["spatial"], float), 6)
        for a in anns:
            p = od / f"{slide}__{a}.csv.gz"
            if p.exists():
                continue
            df, _ = run_scores(adata, labs[a].astype(str), rows, cols, pairs, log, f"{slide}/{a}")
            df = add_flags(df, det)
            df.insert(0, "annotation", a)
            df.insert(0, "section", slide)
            df.to_csv(p, index=False)
            register(p, phase, f"6-NN C/E/M/Q (500 perms) GSE310352 {slide}, annotation {a}")
            log.info("slide %s / %s saved; %s", slide, a, mem_status())
        del adata
        gc.collect()
    allt = pd.concat([pd.read_csv(f) for f in sorted(od.glob("*__*.csv.gz"))], ignore_index=True)
    for key in ca.SCORES:
        allt[f"{key}_q_bh_within_section_annotation"] = allt.groupby(["section", "annotation"])[f"{key}_p_plus1"].transform(
            lambda p: pd.Series(ca.bh(p.to_numpy()), index=p.index))
    hl = allt[(allt.ligand == "CXCL12") & (allt.receptor == "CXCR4") & (allt.sender == "CAF_fibroblast") & (allt.receiver == "T_cell")]
    save_csv(hl, phase, "gse310352_headline_CAF_T_rows.csv", "Headline rows per slide x annotation, GSE310352", log)
    summ = hl.groupby("annotation").agg(n_slides=("section", "nunique"), n_edge_sig=("edge_sig", "sum"),
                                        n_ratio_sig=("ratio_sig", "sum"),
                                        n_ratio_sig_bh=("ratio_q_bh_within_section_annotation", lambda x: int((x < 0.05).sum())),
                                        median_edge_z=("edge_z", "median"), median_ratio_z=("ratio_z", "median")).reset_index()
    save_csv(summ, phase, "gse310352_headline_summary.csv", "Headline summary per annotation, GSE310352", log)
    save_json({"headline": summ.to_dict(orient="records")}, phase, "gse310352_summary.json", "GSE310352 summary", log)
    finish_run(log, phase)


if __name__ == "__main__":
    main()

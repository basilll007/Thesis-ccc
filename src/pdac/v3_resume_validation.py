"""Independent checks of resumed Phase 4/5 outputs and patient-level inference."""
from __future__ import annotations

import itertools
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "data/pdac/numba_cache"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import rankdata


def validate_phase5():
    p = ROOT / "results/pdac/v3_phase5"
    qc = pd.read_csv(p / "slide_qc.csv")
    scores = pd.read_csv(p / "slide_program_scores.csv")
    patient = pd.read_csv(p / "patient_site_program_scores.csv")
    arms = pd.read_csv(p / "patient_arm_program_scores.csv")
    tests = pd.read_csv(p / "program_comparisons.csv")
    metadata = pd.read_csv(p / "verified_sample_metadata.csv").set_index("gsm")
    for row in qc.itertuples(index=False):
        with h5py.File(ROOT / "eda_gse277783/extracted" / metadata.loc[row.gsm, "h5_file"], "r") as h:
            m = h["matrix"]
            x = csc_matrix((m["data"][:].astype(np.int64), m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
        counts = np.asarray(x.sum(axis=0)).ravel()
        assert len(counts) == row.n_spots
        assert np.sum(counts >= 500) == row.n_retained
        assert np.isclose(np.median(counts), row.median_umi)
        expected = "PDX" not in row.site and np.mean(counts < 500) <= .5
        assert bool(row.included) == expected
        if expected:
            spot = pd.read_csv(p / "slides" / f"{row.gsm}_spot_scores.csv.gz")
            assert len(spot) == row.n_retained and spot.barcode.is_unique
            assert (spot.total_umi >= 500).all()
            for program in scores.program.unique():
                value = scores[(scores.gsm == row.gsm) & (scores.program == program)].score.item()
                assert np.isclose(value, spot[program].median(), atol=1e-12)
    recalculated = scores.groupby(["patient", "site", "program"]).score.mean()
    stored = patient.set_index(["patient", "site", "program"]).score
    assert np.allclose(recalculated.sort_index(), stored.sort_index(), atol=1e-12)
    patient["arm"] = np.where(patient.site == "Primary PDAC", "primary", "metastasis")
    recalculated = patient.groupby(["patient", "program", "arm"]).score.mean()
    stored = arms.set_index(["patient", "program", "arm"]).score
    assert np.allclose(recalculated.sort_index(), stored.sort_index(), atol=1e-12)
    for row in tests.itertuples(index=False):
        d = arms[arms.program == row.program] if row.scope == "combined_metastasis" else patient[
            (patient.program == row.program) & patient.site.isin(["Primary PDAC", row.scope])]
        paired = d.pivot(index="patient", columns="arm", values="score").dropna()
        assert len(paired) == row.n_pairs
        # Enumerate all signs independently of scipy.stats.wilcoxon.
        delta = np.round((paired.metastasis - paired.primary).to_numpy(), 12)
        delta = delta[delta != 0]
        ranks = rankdata(np.abs(delta))
        plus_observed = ranks[delta > 0].sum()
        null = np.array([np.dot(signs, ranks) for signs in itertools.product([0, 1], repeat=len(ranks))])
        exact = min(1., 2 * min(np.mean(null <= plus_observed), np.mean(null >= plus_observed)))
        assert np.isclose(exact, row.p, atol=1e-12), (row.program, row.scope, exact, row.p)
    coverage = pd.read_csv(p / "program_gene_coverage_controls.csv")
    for row in coverage.itertuples(index=False):
        present, controls = row.present_genes.split(";"), row.control_genes.split(";")
        assert len(present) == row.n_present >= 5
        assert len(controls) == row.n_controls and not (set(present) & set(controls))
    # Recalculate BH without importing the analysis helper.
    for main in (True, False):
        d = tests[tests.scope.eq("combined_metastasis") == main]
        vals = d.p.to_numpy()
        order = np.argsort(vals)
        adj = np.minimum.accumulate((vals[order] * len(vals) / np.arange(1, len(vals) + 1))[::-1])[::-1]
        assert np.allclose(np.minimum(adj, 1), d.q_bh.to_numpy()[order])
    return {"phase": 5, "checks": "raw matrix QC; spot/slide/patient aggregation; control disjointness; exhaustive Wilcoxon; BH", "status": "PASS", "n_slides_checked": len(qc), "n_tests_checked": len(tests)}


def validate_phase4():
    p = ROOT / "results/pdac/v3_phase4"
    files = sorted((p / "slides").glob("*__*.csv.gz"))
    assert len(files) == 80, len(files)
    observations = {}
    n_evaluated = n_skipped = 0
    for file in files:
        d = pd.read_csv(file)
        if "skipped" in d:
            n_skipped += 1
            continue
        n_evaluated += 1
        for score in ["comp", "edge", "matched", "ratio"]:
            pv = d[f"{score}_p_plus1"].dropna()
            assert ((pv > 0) & (pv <= 1)).all()
            n = d[f"{score}_n_null_valid"].to_numpy()
            scaled = d[f"{score}_p_plus1"].to_numpy() * (n + 1)
            assert np.allclose(scaled[np.isfinite(scaled)], np.round(scaled[np.isfinite(scaled)]), atol=1e-8)
            assert (d[f"{score}_q_bh"].dropna() >= 0).all() and (d[f"{score}_q_bh"].dropna() <= 1).all()
        evaluated = (d.matched_obs > 0) & (d.n_qualifying_edges > 0)
        assert np.allclose(d.loc[evaluated, "edge_obs"] / d.loc[evaluated, "matched_obs"], d.loc[evaluated, "ratio_obs"])
        index = ["sender", "receiver", "ligand", "receptor"]
        obs = d.set_index(index)[[f"{score}_obs" for score in ["comp", "edge", "matched", "ratio"]]].sort_index()
        key = (d.section.iloc[0], d.annotation.iloc[0])
        if key in observations:
            assert np.allclose(obs, observations[key], equal_nan=True)
        else:
            observations[key] = obs
    # Explicitly tabulate absent receiver classes; absence is not a null result.
    availability = []
    for file in sorted((p / "slides").glob("*_info.json")):
        info = json.loads(file.read_text())
        for annotation, counts in info["label_counts"].items():
            nt, nc = counts.get("T_cell", 0), counts.get("CAF_fibroblast", 0)
            availability.append({"section": info["slide"], "annotation": annotation,
                                 "n_cells": info["n_cells"], "n_T_cells": nt, "n_CAF_cells": nc,
                                 "median_positive_nn_px": info["median_positive_nn_px"],
                                 "headline_types_present": bool(nt and nc),
                                 "reason_if_absent": "no T-cell class recovered" if not nt else (
                                     "no CAF class recovered" if not nc else "")})
    from common import save_csv
    save_csv(pd.DataFrame(availability), "v3_phase4", "annotation_availability.csv",
             "Explicit cell-type availability; missing receivers mean untestable headline")
    return {"phase": 4, "checks": "80 run markers; add-one p grid; Q=E/M; observed-score invariance across nulls", "status": "PASS", "n_evaluated": n_evaluated, "n_skipped": n_skipped}


if __name__ == "__main__":
    if "--phase5-only" not in sys.argv:
        from common import start_run, register, finish_run
        log = start_run("v3_resume_validation")
        register(log.log_path, "v3_phase6", "Independent resumed-results validation log")
    phase5 = validate_phase5()
    print(phase5, flush=True)
    if "--phase5-only" not in sys.argv:
        phase4 = validate_phase4()
        print(phase4, flush=True)
        from common import save_json
        save_json({"phase4": phase4, "phase5": phase5}, "v3_phase6", "resumed_validation.json",
                  "Independent checks of raw QC, aggregation, exact tests and CCC invariants")
        finish_run(log, "v3_phase6")

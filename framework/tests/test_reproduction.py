"""Reproduction tests that RECOMPUTE thesis numbers from the original inputs with the package.

Plan: framework/plans/v2_phase2.md section 2.1. The v1 version of this file only re-read saved result CSVs.

R1: v3 tumour-level CAF_fibroblast->T_cell CXCL12-CXCR4 (consensus_T, score Q), median tumour mean z_Q and
    number of tumours with mean z_Q > 0, for nulls global / block200 / block100 / block50 / dist.
R2: per-section observed Q, M, E and qualifying edges (global null rows).
R3-R5 (model-v2 XGBoost / ridge quantities) are NOT recomputable by the package and are not tested here
(see framework/reports/FRAMEWORK_REPORT_v2.md).

Inputs: raw Mendeley h5ad + thesis consensus_T label CSVs (an input annotation). Thesis values are read only
for comparison, after recomputation.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellalibi.datasets import MENDELEY_SECTIONS, load_mendeley_section, mendeley_section_to_tumour
from cellalibi.graphs import knn_graph
from cellalibi.nulls import block_null, distance_band_null, global_null
from cellalibi.scores import compute_ccc_scores

REPO = Path(__file__).resolve().parents[2]
H5AD = Path(os.environ.get("CCC_MENDELEY_H5AD", "F:/Thesis/data/pdac/mendeley_pdac_smi/raw_meta_data_final.h5ad"))
LABELS = REPO / "results/pdac/phaseB/sections"
SECTION_INFO = REPO / "results/pdac/phase3/cross_section_section_info.csv"
THESIS_TUMOUR = REPO / "results/pdac/v3_phase1/tumour_level_summary_by_null.csv"
THESIS_SECTION = REPO / "results/pdac/v3_phase1/headline_by_null.csv"
HEAD = ("CAF_fibroblast", "T_cell", "CXCL12", "CXCR4")
NULLS = ["global", "block200", "block100", "block50", "dist"]
N_PERMS, SEED = 500, 42
TOL_Z = 0.02  # declared in v1 and kept

pytestmark = [pytest.mark.repro, pytest.mark.skipif(not H5AD.exists(), reason="local Mendeley h5ad not available")]


def _perms(null, lab, xy, types):
    if null == "global":
        return global_null(lab, n_perms=N_PERMS, seed=SEED)
    if null.startswith("block"):
        return block_null(lab, xy, tile_size_um=float(null[5:]), n_perms=N_PERMS, random_offset=False, seed=SEED)[0]
    if null == "dist":
        if "epithelial_ductal_tumor" not in types:
            return None
        return distance_band_null(lab, xy, types.index("epithelial_ductal_tumor"), n_perms=N_PERMS, seed=SEED)[0]
    raise ValueError(null)


@pytest.fixture(scope="module")
def recomputed():
    tmap = mendeley_section_to_tumour(SECTION_INFO)
    rows = []
    for sec in MENDELEY_SECTIONS:
        a = load_mendeley_section(H5AD, sec, LABELS, annotations=("consensus_T",))
        a.obs["cell_type"] = a.obs["consensus_T"].astype(str)
        xy = np.c_[a.obs.x_um.to_numpy(), a.obs.y_um.to_numpy()]
        r, c = knn_graph(xy, k=6)
        types = sorted(set(a.obs.cell_type))
        lab = np.array([types.index(x) for x in a.obs.cell_type], dtype=np.int64)
        for null in NULLS:
            perms = _perms(null, lab, xy, types)
            if perms is None:
                continue
            df = compute_ccc_scores(a, [HEAD[2:]], [HEAD], r, c, perms).dataframe.iloc[0]
            rows.append({"section": sec, "tumour": tmap[sec], "null": null, "ratio_obs": df.ratio_obs,
                         "matched_obs": df.matched_obs, "edge_obs": df.edge_obs, "ratio_z": df.ratio_z,
                         "n_qualifying_edges": int(df.n_qualifying_edges), "edge_status": df.edge_status})
            del perms
        del a
    out = pd.DataFrame(rows)
    od = REPO / "framework/results/v2_phase2"
    od.mkdir(parents=True, exist_ok=True)
    out.to_csv(od / "reproduction_recomputed_sections.csv", index=False)
    return out


@pytest.mark.parametrize("null", NULLS)
def test_R1_tumour_level_caf_t(recomputed, null):
    th = pd.read_csv(THESIS_TUMOUR)
    th = th[(th.annotation == "consensus_T") & (th.score == "Q") & (th.null == null)].iloc[0]
    d = recomputed[(recomputed.null == null) & (recomputed.edge_status == "evaluated")]
    tz = d.groupby("tumour").ratio_z.mean()
    assert len(tz) == int(th.n_tumours) == 13
    med = float(np.median(tz))
    assert abs(med - float(th.median_tumour_mean_z)) <= TOL_Z, f"{null}: median z {med:.4f} vs thesis {th.median_tumour_mean_z:.4f}"
    assert int((tz > 0).sum()) == int(th.n_mean_z_pos), f"{null}: n z>0 {(tz > 0).sum()} vs thesis {th.n_mean_z_pos}"


def test_R2_section_observed_scores(recomputed):
    th = pd.read_csv(THESIS_SECTION)
    th = th[(th.annotation == "consensus_T") & (th.null == "global")].set_index("section")
    d = recomputed[recomputed.null == "global"].set_index("section")
    assert set(d.index) == set(th.index)
    for col in ("ratio_obs", "matched_obs", "edge_obs"):
        np.testing.assert_allclose(d.loc[th.index, col].to_numpy(), th[col].to_numpy(), rtol=1e-6, err_msg=col)
    assert (d.loc[th.index, "n_qualifying_edges"].to_numpy() == th.n_qualifying_edges.to_numpy()).all()

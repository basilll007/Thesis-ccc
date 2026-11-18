"""Null-calibration test for the v2 contact tests (plan: framework/plans/v2_phase2.md section 2.4).

200 synthetic tissues without any contact effect, 1000 permutations each. The number of rejections at
alpha = 0.05 must fall inside the central binomial 95 % interval of Binomial(200, 0.05).
"""
import numpy as np
import pytest
from scipy.stats import binom

from ccc_audit.contact import AuditConfig, ExpressionData, SpatialContext, contact_statistics
from ccc_audit.simulate import synthetic_tissue

N_DATASETS = 200
N_PERMS = 1000
ALPHA = 0.05


@pytest.fixture(scope="module")
def null_pvalues():
    cfg = AuditConfig(n_perms=N_PERMS, tile_um=100.0)
    rows = []
    for s in range(N_DATASETS):
        xy, labels, counts, genes = synthetic_tissue(seed=s)
        ctx = SpatialContext(xy, labels, AuditConfig(n_perms=N_PERMS, tile_um=100.0, seed=s))
        ex = ExpressionData(counts, genes)
        m = contact_statistics(ctx, ex, "Sender", "Receiver", "L", "R", cfg=ctx.cfg, seed=1000 + s,
                               run_spillover=False)
        rows.append((m["p_contact_global"], m.get("p_block_contrast", np.nan), m.get("p_ring", np.nan), m["p_contact_block"]))
    return np.array(rows), cfg


def _interval():
    lo, hi = binom.ppf([0.025, 0.975], N_DATASETS, ALPHA)
    return int(lo), int(hi)


@pytest.mark.slow
@pytest.mark.parametrize("col,name", [(0, "global"), (1, "block_contrast"), (2, "ring"), (3, "block_q_info")])
def test_contact_test_calibrated_under_null(null_pvalues, col, name):
    P, cfg = null_pvalues
    p = P[:, col]
    assert np.isfinite(p).sum() == N_DATASETS, f"{name}: test not runnable in every dataset"
    assert p.min() >= 1.0 / (N_PERMS + 1) - 1e-12
    k = int((p < ALPHA).sum())
    lo, hi = _interval()
    assert lo <= k <= hi, f"{name}: {k}/{N_DATASETS} rejections outside binomial 95% interval [{lo}, {hi}]"

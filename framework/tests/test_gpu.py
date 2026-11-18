"""GPU backend tests (skipped when CUDA is unavailable).

1. Null calibration on the GPU: same 200 synthetic null tissues and 1000 permutations as test_calibration.py;
   rejections must fall inside the binomial 95 % interval.
2. CPU/GPU agreement: identical statistics; p-values agree within Monte-Carlo error and verdicts match.
"""
import numpy as np
import pytest
from scipy.stats import binom

from cellalibi import gpu
from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext, contact_statistics
from cellalibi.simulate import synthetic_tissue
from cellalibi.verdicts import audit_claim

pytestmark = pytest.mark.skipif(not gpu.available(), reason="CUDA not available")
N, PERMS = 200, 1000


@pytest.fixture(scope="module")
def gpu_null_p():
    rows = []
    for s in range(N):
        xy, labels, counts, genes = synthetic_tissue(seed=s)
        cfg = AuditConfig(n_perms=PERMS, tile_um=100.0, seed=s, device="cuda")
        ctx = SpatialContext(xy, labels, cfg)
        m = contact_statistics(ctx, ExpressionData(counts, genes), "Sender", "Receiver", "L", "R", cfg=cfg,
                               seed=1000 + s, run_spillover=False)
        rows.append((m["p_contact_global"], m["p_block_contrast"], m["p_ring"], m["p_rc_global"]))
    return np.array(rows)


@pytest.mark.slow
@pytest.mark.parametrize("col,name", [(0, "global"), (1, "block_contrast"), (2, "ring"), (3, "rc_global")])
def test_gpu_null_calibration(gpu_null_p, col, name):
    lo, hi = binom.ppf([0.025, 0.975], N, 0.05)
    k = int((gpu_null_p[:, col] < 0.05).sum())
    assert lo <= k <= hi, f"{name}: {k}/{N} outside [{lo}, {hi}]"


@pytest.mark.parametrize("kw", [dict(contact_effect=3.0), dict(spill_frac=0.3), {}])
def test_cpu_gpu_agree(kw):
    xy, labels, counts, genes = synthetic_tissue(seed=21, n_cells=8000, size_um=1800.0, **kw)
    res = {}
    for dev in ("cpu", "cuda"):
        cfg = AuditConfig(n_perms=PERMS, seed=21, device=dev)
        ctx = SpatialContext(xy, labels, cfg)
        res[dev] = audit_claim(ctx, ExpressionData(counts, genes), "Sender", "Receiver", "L", "R", cfg=cfg, seed=5)
    a, b = res["cpu"].metrics, res["cuda"].metrics
    for k in ("E", "M", "Q", "ring_diff", "rc_global_diff", "spill_phi", "spill_dR_obs"):
        assert np.isclose(a[k], b[k], rtol=1e-9, equal_nan=True), k
    for k in ("p_contact_global", "p_rc_global", "p_block_contrast", "p_ring", "spill_p_phi"):
        assert abs(a[k] - b[k]) <= 0.06, (k, a[k], b[k])
    assert res["cpu"].primary == res["cuda"].primary

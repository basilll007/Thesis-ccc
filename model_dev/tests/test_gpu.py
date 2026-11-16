"""GPU baseline implementations must match the CPU references (sklearn Ridge; features.fit_scaler)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
torch = pytest.importorskip("torch")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")


def test_gpu_ridge_matches_sklearn():
    import gpu as GP
    from sklearn.linear_model import Ridge
    rng = np.random.default_rng(0)
    X = rng.normal(size=(500, 30)).astype(np.float32)
    Y = (X @ rng.normal(size=(30, 3)) + rng.normal(size=(500, 3))).astype(np.float32)
    rp = GP.RidgePath(GP.to_gpu(X), GP.to_gpu(Y))
    for a in (0.1, 10.0, 1000.0):
        ref = Ridge(alpha=a).fit(X, Y).predict(X)
        got = rp.predict(GP.to_gpu(X), a).cpu().numpy()
        assert np.allclose(got, ref, atol=1e-3), a


def test_gpu_scaler_matches_cpu_scaler_with_nan():
    import features as F
    import gpu as GP
    rng = np.random.default_rng(1)
    X = rng.normal(size=(300, 5)).astype(np.float32)
    X[::7, 2] = np.nan
    X[:, 4] = 3.0
    mu, sd = F.fit_scaler(X)
    gmu, gsd = GP.fit_scaler(GP.to_gpu(X))
    assert np.allclose(gmu.cpu().numpy(), mu, atol=1e-5) and np.allclose(gsd.cpu().numpy(), sd, atol=1e-5)
    assert np.allclose(GP.apply_scaler(GP.to_gpu(X), gmu, gsd).cpu().numpy(), F.apply_scaler(X, mu, sd), atol=1e-5)

"""Metric tests: pooled R2, its tumour bootstrap, and the SIGN of the paired difference (a bug once inverted it)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import mdcommon as C  # noqa: E402


def test_pooled_r2_and_bootstrap_point():
    pairs = {"t1": (2.0, 10.0), "t2": (3.0, 10.0)}
    r2, lo, hi = C.tumour_bootstrap(pairs, 200)
    assert np.isclose(r2, 0.75) and lo <= r2 <= hi


def test_paired_difference_sign_model_first():
    better = {"t1": (1.0, 10.0), "t2": (2.0, 10.0)}   # R2 = 0.85
    worse = {"t1": (4.0, 10.0), "t2": (4.0, 10.0)}    # R2 = 0.60
    d, lo, hi = C.tumour_bootstrap_diff(better, worse, 500)
    assert np.isclose(d, 0.25) and lo > 0, "gain of the better model must be positive"
    d2, _, _ = C.tumour_bootstrap_diff(worse, better, 500)
    assert np.isclose(d2, -0.25)

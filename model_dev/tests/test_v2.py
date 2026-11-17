"""Tests for the v2 code: column typing, nearest-type distances, Gram-ridge sub-blocks, metric plumbing, the share
bootstrap and the leakage guard on the contamination covariate."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import v2common as V  # noqa: E402
import features as F  # noqa: E402


def test_column_type_and_kind():
    assert V.column_type("knn_frac__CAF_fibroblast") == "CAF"
    assert V.column_type("r30_count__epithelial_ductal_tumor") == "tumour"
    assert V.column_type("r30_count__B_plasma") == "other"
    assert V.column_type("nb_knn_CAF__CAF:myCAF") == "CAF"
    assert V.column_type("nb_r30_myeloid__myeloid:C1Q_TAM") == "myeloid"
    assert V.column_type("nb_knn_T__CD3E") == "T"
    assert V.column_type("nb_knn_other__COL1A1") == "other"
    assert V.column_type("ncaf_expr__COL1A1") == "CAF"
    assert V.column_type("tile_loo_mean__C1QA") == "none"
    assert V.column_kind("ncaf_dist_min") == "composition"
    assert V.column_kind("ncaf_expr__DCN") == "expression"
    assert V.column_kind("knn_frac__myeloid") == "composition"
    assert V.column_kind("s_c_contamination") == "own"


def test_nearest_type_distance_excludes_self():
    xy = np.array([[0, 0], [1, 0], [5, 0], [6, 0]], float)
    lab = np.array(["CAF", "CAF", "m", "CAF"])
    d, j = V.nearest_type_distance(xy, lab, "CAF", k=1)
    assert j[0, 0] == 1 and d[0, 0] == pytest.approx(1.0)      # a CAF's nearest CAF is not itself
    assert j[2, 0] == 3 and d[2, 0] == pytest.approx(1.0)
    d3, j3 = V.nearest_type_distance(xy, lab, "CAF", k=2)
    assert set(j3[0]) == {1, 3}


def test_gram_ridge_subblock_matches_sklearn():
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available():
        pytest.skip("CUDA required")
    from sklearn.linear_model import Ridge
    import v2engine as E
    rng = np.random.default_rng(0)
    X = rng.standard_normal((500, 12)).astype(np.float32)
    Y = (X[:, :3] @ rng.standard_normal((3, 2)) + rng.standard_normal((500, 2))).astype(np.float32)
    g = E.GramRidge(torch.as_tensor(X, device="cuda"), torch.as_tensor(Y, device="cuda"), chunk=128)
    cols = np.array([0, 2, 5, 7, 11])
    c = torch.as_tensor(cols, device="cuda")
    for a, W in zip([0.1, 10.0], g.coefs(c, [0.1, 10.0])):
        pred = g.predict(torch.as_tensor(X, device="cuda"), c, W).cpu().numpy()
        ref = Ridge(alpha=a).fit(X[:, cols].astype(np.float64), Y).predict(X[:, cols].astype(np.float64))
        assert np.abs(pred - ref).max() < 1e-3


def test_pairs_for_and_ratio_bootstrap():
    import v2engine as E
    import mdcommon as C
    rng = np.random.default_rng(1)
    n = 200
    group = np.repeat(["a", "b", "c", "d"], 50)
    Y = rng.standard_normal((n, 1)).astype(np.float32)
    good = Y[:, 0] + 0.1 * rng.standard_normal(n)
    bad = Y[:, 0] + 0.5 * rng.standard_normal(n)
    parts = [{"test_idx": np.arange(n), "ytr_mean": np.zeros(1), "y_names": np.array(["t"]), "alphas": np.array([1.0]),
              "ridge_pred__good": good[None, :, None].astype(np.float32), "ridge_valmse__good": np.zeros((1, 1)),
              "ridge_pred__bad": bad[None, :, None].astype(np.float32), "ridge_valmse__bad": np.zeros((1, 1)),
              "ridge_pred__zero": np.zeros((1, n, 1), np.float32), "ridge_valmse__zero": np.zeros((1, 1))}]
    pg = E.pairs_for(parts, group, "ridge", "good", "t", Y=Y)
    r2 = C.tumour_bootstrap(pg, 200)[0]
    manual = 1 - ((Y[:, 0] - good) ** 2).sum() / (Y[:, 0] ** 2).sum()
    assert r2 == pytest.approx(manual, rel=1e-5)
    pb = E.pairs_for(parts, group, "ridge", "bad", "t", Y=Y)
    pz = E.pairs_for(parts, group, "ridge", "zero", "t", Y=Y)
    # share of (zero -> good) explained by (bad -> good): (R2good - R2bad) / (R2good - R2zero)
    share, lo, hi = E.ratio_bootstrap(pg, pb, pg, pz, 500)
    r = lambda p: C.tumour_bootstrap(p, 10)[0]
    assert share == pytest.approx((r(pg) - r(pb)) / (r(pg) - r(pz)), rel=1e-6)
    assert 0 < share < 1
    # mask restricts to the selected cells
    mask = group == "a"
    pa = E.pairs_for(parts, group, "ridge", "good", "t", cell_mask=mask, Y=Y)
    assert list(pa) == ["a"]


def test_leak_guard_refuses_contamination_covariate_on_marker_target():
    rng = np.random.default_rng(2)
    n = 100
    fs = F.FeatureSet(rows=np.arange(n))
    fs.add("own", rng.standard_normal((n, 1)), ["s_c_contamination"], [F.Provenance("self_covariate", genes=frozenset({"COL1A1", "DCN"}))])
    y = rng.standard_normal((n, 1))
    F.assert_no_target_leak(fs, None, y, {"C1QA"}, set())          # disjoint target: allowed
    with pytest.raises(F.LeakError):
        F.assert_no_target_leak(fs, None, y, {"DCN"}, set())       # marker gene as target: refused

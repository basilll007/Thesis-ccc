"""Leakage guard tests: these must FAIL if a cell's own target can reach its inputs, or if scaling peeks at test rows."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import features as F  # noqa: E402


def _toy_fs(n: int = 200, seed: int = 0):
    rng = np.random.default_rng(seed)
    fs = F.FeatureSet(rows=np.arange(n))
    fs.add("own", rng.normal(size=(n, 2)).astype(np.float32), ["log1p_total_counts", "density_r30"],
           [F.Provenance("self_covariate")] * 2)
    fs.add("nb", rng.normal(size=(n, 3)).astype(np.float32), ["nb_knn_CAF__C1QA", "nb_knn_tumour__KRT19", "nb_knn_myeloid__myeloid:C1Q_TAM"],
           [F.Provenance("neighbour", genes=frozenset(["C1QA"])), F.Provenance("neighbour", genes=frozenset(["KRT19"])),
            F.Provenance("neighbour", programs=frozenset(["myeloid:C1Q_TAM"]))])
    target = rng.normal(size=(n, 2)).astype(np.float32)
    return fs, target


def test_clean_features_pass():
    fs, y = _toy_fs()
    F.assert_no_target_leak(fs, None, y, {"C1QA", "C1QB"}, {"myeloid:C1Q_TAM"})


def test_planted_target_copy_is_detected_even_with_innocent_tag():
    fs, y = _toy_fs()
    fs.add("planted", y[:, :1] * 3.0 + 1.0, ["some_covariate"], [F.Provenance("self_covariate")])
    with pytest.raises(F.LeakError):
        F.assert_no_target_leak(fs, None, y, {"C1QA"}, {"myeloid:C1Q_TAM"})


def test_self_sourced_target_gene_is_detected_by_provenance():
    fs, y = _toy_fs()
    fs.add("own_expr", np.random.default_rng(1).normal(size=(len(y), 1)).astype(np.float32), ["own__C1QB"],
           [F.Provenance("self_covariate", genes=frozenset(["C1QB"]))])
    with pytest.raises(F.LeakError):
        F.assert_no_target_leak(fs, None, y, {"C1QA", "C1QB"}, set())


def test_self_sourced_target_program_is_detected_by_provenance():
    fs, y = _toy_fs()
    fs.add("own_prog", np.zeros((len(y), 1), np.float32) + 0.5, ["own_score"],
           [F.Provenance("self_covariate", programs=frozenset(["myeloid:C1Q_TAM"]))])
    with pytest.raises(F.LeakError):
        F.assert_no_target_leak(fs, None, y, set(), {"myeloid:C1Q_TAM"})


def test_neighbour_copy_of_target_gene_is_allowed():
    """Other cells' values of a target gene are legitimate inputs (that is the point of a neighbourhood model)."""
    fs, y = _toy_fs()
    F.assert_no_target_leak(fs, None, y, {"C1QA"}, {"myeloid:C1Q_TAM"})


def test_tile_loo_mean_excludes_the_cell_itself():
    tile = np.array([0, 0, 0, 1, 1])
    member = np.array([True, True, True, True, False])
    v = np.array([[1.0], [2.0], [6.0], [4.0], [100.0]])
    m, has = F.tile_loo_mean(v, tile, member)
    assert np.allclose(m[:3, 0], [4.0, 3.5, 1.5])      # (2+6)/2, (1+6)/2, (1+2)/2
    assert not has[3]                                   # alone in its tile -> no other member
    assert np.isclose(m[4, 0], 4.0) and has[4]          # non-member sees the member mean of its tile


def test_scaler_fit_on_training_rows_only():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(100, 3)).astype(np.float32)
    train, test = np.arange(60), np.arange(60, 100)
    X[test] += 50.0  # a shift that would show up in a scaler fitted on all rows
    mu, sd = F.fit_scaler(X[train])
    assert np.allclose(mu, X[train].mean(axis=0), atol=1e-5)
    assert np.all(np.abs(mu) < 5), "scaler mean contaminated by test rows"
    Z = F.apply_scaler(X[test], mu, sd)
    assert np.all(Z.mean(axis=0) > 20), "test rows must keep their shift after train-fitted scaling"


def test_two_hop_tree_excludes_seed():
    """Padded 2-hop neighbourhoods must never contain the seed cell at hop 2 (the chain SGC bug)."""
    import graphs as G
    # triangle 0-1-2 plus 3 attached to 2
    indptr = np.array([0, 2, 4, 7, 8]); indices = np.array([1, 2, 0, 2, 0, 1, 3, 2])
    nb1, nb2 = G.padded_two_hop(indptr, indices, seeds=np.array([0]), k1=3, k2=3)
    assert 0 not in nb2[0][nb2[0] >= 0], "seed leaked into its own 2-hop neighbourhood"

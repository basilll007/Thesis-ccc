"""Patient-level split tests: a tumour never straddles folds, every fold has >= 2 tumours, inner val is disjoint."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import splits as S  # noqa: E402


def test_grouped_folds_rule_guarantees_two_groups_per_fold():
    counts = pd.Series({f"t{i}": 1000 * (13 - i) for i in range(13)})
    f = S.grouped_folds(counts, 5)
    sizes = pd.Series(f).value_counts()
    assert set(f.values()) == set(range(5))
    assert (sizes >= 2).all(), sizes.to_dict()


@pytest.mark.skipif(not (ROOT / "results/phase1/folds.csv").exists(), reason="Phase 1 not run yet")
def test_written_folds_keep_sections_of_a_tumour_together_and_meet_floor():
    folds = pd.read_csv(ROOT / "results/phase1/folds.csv")
    assert (folds.groupby("tumour").fold.nunique() == 1).all(), "a tumour straddles folds"
    assert (folds.groupby("fold").tumour.nunique() >= 2).all(), "a fold has fewer than 2 tumours"
    chk = pd.read_csv(ROOT / "results/phase1/fold_floor_check.csv")
    assert chk.passes_floor.all(), chk.to_string()


def test_inner_validation_fold_is_disjoint_from_test():
    import assemble as A
    section = np.array(["a", "a", "b", "c", "d", "e"])
    folds = pd.DataFrame({"section": list("abcde"), "fold": [0, 1, 2, 3, 4]})

    class D:  # minimal stand-in for ReceiverData
        pass
    d = D(); d.section = section
    tr_in, val, tr_all, test = A.fold_masks(d, folds, 0, 5)
    assert not (val & test).any() and not (tr_in & val).any() and not (tr_in & test).any()
    assert (tr_in | val | test).all()

"""Reproduction tests against thesis results from F:\\Thesis.

Checks that the thesis benchmarks reproduce within pre-declared tolerances:
1. v3 tumour-level CAF->T result (consensus_T, score Q across 13 tumours for 5 nulls).
2. model-v2 spillover reference Q for C1Q_TAM (+0.004).
3. model-v2 XGBoost cross-type share (-0.005).
4. model-v1 ridge neighbour gains (+0.033).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import yaml

CONFIG_PATH = Path(__file__).parent / "local_reproduction_config.yaml"


@pytest.fixture(scope="module")
def reproduction_config():
    if not CONFIG_PATH.exists():
        pytest.skip(f"Local reproduction config '{CONFIG_PATH}' not found.")
    with open(CONFIG_PATH, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg


def test_reproduce_v3_tumour_caf_t(reproduction_config):
    """Reproduce v3 tumour-level CAF->T result within pre-declared tolerances."""
    root = Path(reproduction_config["repo_root"])
    spec = reproduction_config["benchmarks"]["v3_tumour_caf_t"]
    file_path = root / spec["file_path"]
    assert file_path.exists(), f"Benchmark file '{file_path}' does not exist."

    df = pd.read_csv(file_path)
    sub = df[(df["annotation"] == spec["annotation"]) & (df["score"] == spec["score"])]

    assert len(sub) >= 5, f"Expected at least 5 null models, found {len(sub)}"

    for null_name, expected_z in spec["tolerances"]["median_tumour_mean_z"].items():
        row = sub[sub["null"] == null_name]
        assert len(row) == 1, f"Missing null {null_name} in v3 results"
        actual_z = float(row["median_tumour_mean_z"].iloc[0])
        tol = spec["tolerances"]["tolerance_z"]
        assert abs(actual_z - expected_z) <= tol, (
            f"Null {null_name}: expected median z_Q {expected_z} +- {tol}, got {actual_z}"
        )

        expected_n_pos = spec["tolerances"]["n_mean_z_pos"][null_name]
        actual_n_pos = int(row["n_mean_z_pos"].iloc[0])
        assert actual_n_pos == expected_n_pos, (
            f"Null {null_name}: expected {expected_n_pos} positive tumours, got {actual_n_pos}"
        )

    assert int(sub["n_tumours"].iloc[0]) == spec["expected_n_tumours"]


def test_reproduce_model_v2_spillover_c1q_tam(reproduction_config):
    """Reproduce model-v2 spillover reference Q for C1Q_TAM (+0.004)."""
    root = Path(reproduction_config["repo_root"])
    spec = reproduction_config["benchmarks"]["v2_spillover_c1q_tam"]
    file_path = root / spec["file_path"]
    assert file_path.exists(), f"Benchmark file '{file_path}' does not exist."

    df = pd.read_csv(file_path)
    row = df[
        (df["target"] == spec["target"]) &
        (df["model"] == spec["model"]) &
        (df["group"] == spec["group"]) &
        (df["stratum"] == "all") &
        (df["level"] == "all")
    ]

    assert len(row) == 1, "myeloid:C1Q_TAM reference row not found in phase1_summary.csv"
    actual_q = float(row["Q_caf_expr"].iloc[0])
    expected_q = spec["expected_q"]
    tol = spec["tolerance_q"]

    assert abs(actual_q - expected_q) <= tol, (
        f"C1Q_TAM spillover reference Q: expected {expected_q} +- {tol}, got {actual_q}"
    )


def test_reproduce_model_v2_xgboost_cross_type(reproduction_config):
    """Reproduce model-v2 XGBoost cross-type share (-0.005)."""
    root = Path(reproduction_config["repo_root"])
    spec = reproduction_config["benchmarks"]["v2_xgboost_cross_type"]
    file_path = root / spec["file_path"]
    assert file_path.exists(), f"Benchmark file '{file_path}' does not exist."

    df = pd.read_csv(file_path)
    row = df[
        (df["receiver"] == spec["receiver"]) &
        (df["target"] == spec["target"]) &
        (df["model"] == spec["model"]) &
        (df["quantity"] == spec["quantity"])
    ]
    assert len(row) == 1, "myeloid C1Q_TAM cross_type row not found in phase2_quantities.csv"
    actual_gain = float(row["value"].iloc[0])
    expected_gain = spec["expected_cross_type_gain"]
    tol = spec["tolerance_gain"]

    assert abs(actual_gain - expected_gain) <= tol, (
        f"XGBoost cross-type gain: expected {expected_gain} +- {tol}, got {actual_gain}"
    )


def test_reproduce_model_v1_ridge_neighbour_gains(reproduction_config):
    """Reproduce model-v1 ridge neighbour gains (+0.033)."""
    root = Path(reproduction_config["repo_root"])
    spec = reproduction_config["benchmarks"]["v1_ridge_neighbour_gains"]
    file_path = root / spec["file_path"]
    assert file_path.exists(), f"Benchmark file '{file_path}' does not exist."

    df = pd.read_csv(file_path)
    row = df[
        (df["receiver"] == spec["receiver"]) &
        (df["model"] == spec["model"]) &
        (df["target"] == spec["target"])
    ]
    assert len(row) == 1, "myeloid ridge C1Q_TAM row not found in reproduction_check_vs_v1.csv"
    actual_v1 = float(row["v1_r2"].iloc[0])
    actual_v2 = float(row["v2_r2"].iloc[0])

    expected_v1 = spec["expected_v1_r2"]
    expected_v2 = spec["expected_v2_r2"]
    tol = spec["tolerance_r2"]

    assert abs(actual_v1 - expected_v1) <= tol, (
        f"Ridge v1 R2: expected {expected_v1} +- {tol}, got {actual_v1}"
    )
    assert abs(actual_v2 - expected_v2) <= tol, (
        f"Ridge v2 R2: expected {expected_v2} +- {tol}, got {actual_v2}"
    )
    assert bool(row["pass"].iloc[0]) is True, "Reproduction pass flag is False"

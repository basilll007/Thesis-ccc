"""Machine learning neighbourhood models and same-type vs cross-type feature ablations.

What this module does:
Fits predictive models to test whether cell state (gene or program expression) is predicted
by neighbouring cells beyond cell-intrinsic covariates and regional zonation:
1. Feature groups:
   - Own covariates: cell-intrinsic baseline and coarse regional tile mean.
   - Same-type neighbour features: aggregated features from neighbouring cells of the SAME type.
   - Cross-type neighbour features: aggregated features from neighbouring cells of OTHER types.
2. Grouped Cross-Validation:
   - Strict grouping by patient_id: cells from the same patient NEVER split across train and test.
3. Feature Ablation & Attribution:
   - Evaluates held-out R2 under full neighbours vs same-type-only neighbours.
   - Cross-type contribution = R2_full - R2_same_type.
4. Tumour-Level Paired Bootstrap:
   - Computes 95% confidence intervals by resampling patient/tumour units (2,000 bootstrap draws).

Why:
Apparent neighbourhood predictive power is often driven almost entirely by same-type homotypic
clustering or shared regional gradients rather than heterotypic cross-type signaling.
Ablating same-type vs cross-type neighbours tests how much of the predictive gain is truly
attributable to other cell types.

Adapted from:
`model_dev/src/v2engine.py`, `model_dev/src/baselines.py`, and `model_dev/src/v2_phase2.py`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score


@dataclass
class ModelAblationResults:
    """Held-out R2 metrics, cross-type shares, and patient-level bootstrap CIs."""
    r2_own: float
    r2_all_neighbours: float
    r2_same_type_only: float
    neighbour_gain: float
    cross_type_contribution: float
    cross_type_share_pct: float
    ci_neighbour_gain: tuple[float, float]
    ci_cross_type: tuple[float, float]
    per_patient_r2: pd.DataFrame


def build_neighbourhood_features(
    X_features: np.ndarray,
    cell_types: np.ndarray,
    receiver_mask: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate neighbourhood features into same-type and cross-type averages.

    Source file: `model_dev/src/baselines.py`.
    """
    n_cells, n_feat = X_features.shape
    receiver_indices = np.where(receiver_mask)[0]
    n_rec = len(receiver_indices)

    rec_map = {idx: i for i, idx in enumerate(receiver_indices)}

    same_sums = np.zeros((n_rec, n_feat), dtype=float)
    same_counts = np.zeros(n_rec, dtype=float)
    cross_sums = np.zeros((n_rec, n_feat), dtype=float)
    cross_counts = np.zeros(n_rec, dtype=float)

    # In directed edge u -> v, v is the receiver
    for u, v in zip(rows, cols):
        if v in rec_map:
            v_i = rec_map[v]
            if cell_types[u] == cell_types[v]:
                same_sums[v_i] += X_features[u]
                same_counts[v_i] += 1.0
            else:
                cross_sums[v_i] += X_features[u]
                cross_counts[v_i] += 1.0

    with np.errstate(invalid="ignore", divide="ignore"):
        same_means = np.where(same_counts[:, None] > 0, same_sums / same_counts[:, None], 0.0)
        cross_means = np.where(cross_counts[:, None] > 0, cross_sums / cross_counts[:, None], 0.0)

    return same_means, cross_means


def patient_bootstrap_ci(
    y_true: np.ndarray,
    y_pred_base: np.ndarray,
    y_pred_full: np.ndarray,
    patient_ids: np.ndarray,
    n_bootstraps: int = 2000,
    seed: int = 42,
) -> tuple[float, float]:
    """Compute paired patient-level bootstrap 95% confidence interval for delta R2.

    Source file: `model_dev/src/v2_phase2.py`.
    """
    unique_patients = np.unique(patient_ids)
    if len(unique_patients) <= 1:
        # Cannot bootstrap on single patient
        diff = float(r2_score(y_true, y_pred_full) - r2_score(y_true, y_pred_base))
        return (diff, diff)

    rng = np.random.default_rng(seed)
    deltas = []

    # Map patients to cell indices
    patient_to_indices = {p: np.where(patient_ids == p)[0] for p in unique_patients}
    m = len(unique_patients)

    for _ in range(n_bootstraps):
        boot_patients = rng.choice(unique_patients, size=m, replace=True)
        boot_idx = np.concatenate([patient_to_indices[p] for p in boot_patients])

        r2_b = r2_score(y_true[boot_idx], y_pred_base[boot_idx])
        r2_f = r2_score(y_true[boot_idx], y_pred_full[boot_idx])
        deltas.append(r2_f - r2_b)

    ci_low = float(np.percentile(deltas, 2.5))
    ci_high = float(np.percentile(deltas, 97.5))
    return (ci_low, ci_high)


def evaluate_neighbour_prediction(
    X_own: np.ndarray,
    X_features: np.ndarray,
    y: np.ndarray,
    cell_types: np.ndarray,
    patient_ids: np.ndarray,
    receiver_mask: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    n_splits: int = 5,
    seed: int = 42,
) -> ModelAblationResults:
    """Fit own vs same-type vs cross-type models using grouped cross-validation.

    Source file: `model_dev/src/v2engine.py`.
    """
    y_rec = y[receiver_mask]
    patients_rec = patient_ids[receiver_mask]
    X_own_rec = X_own[receiver_mask]

    same_nbrs, cross_nbrs = build_neighbourhood_features(
        X_features, cell_types, receiver_mask, rows, cols
    )

    X_all_nbrs = np.column_stack([X_own_rec, same_nbrs, cross_nbrs])
    X_same_only = np.column_stack([X_own_rec, same_nbrs])

    n_unique_patients = len(np.unique(patients_rec))
    actual_splits = min(n_splits, n_unique_patients)

    preds_own = np.zeros_like(y_rec)
    preds_all = np.zeros_like(y_rec)
    preds_same = np.zeros_like(y_rec)

    if actual_splits > 1:
        gkf = GroupKFold(n_splits=actual_splits)
        split_gen = gkf.split(X_own_rec, y_rec, groups=patients_rec)
    else:
        # Fallback if only 1 patient: train/test on self (with warning)
        split_gen = [(np.arange(len(y_rec)), np.arange(len(y_rec)))]

    for train_idx, test_idx in split_gen:
        # Fit Own model
        m_own = RidgeCV(alphas=np.logspace(-2, 4, 7))
        m_own.fit(X_own_rec[train_idx], y_rec[train_idx])
        preds_own[test_idx] = m_own.predict(X_own_rec[test_idx])

        # Fit Same-type only model
        m_same = RidgeCV(alphas=np.logspace(-2, 4, 7))
        m_same.fit(X_same_only[train_idx], y_rec[train_idx])
        preds_same[test_idx] = m_same.predict(X_same_only[test_idx])

        # Fit All-neighbours model
        m_all = RidgeCV(alphas=np.logspace(-2, 4, 7))
        m_all.fit(X_all_nbrs[train_idx], y_rec[train_idx])
        preds_all[test_idx] = m_all.predict(X_all_nbrs[test_idx])

    r2_own = float(r2_score(y_rec, preds_own))
    r2_all = float(r2_score(y_rec, preds_all))
    r2_same = float(r2_score(y_rec, preds_same))

    nbr_gain = r2_all - r2_own
    cross_contrib = r2_all - r2_same
    cross_share = (cross_contrib / max(nbr_gain, 1e-6)) * 100.0

    ci_nbr = patient_bootstrap_ci(y_rec, preds_own, preds_all, patients_rec, seed=seed)
    ci_cross = patient_bootstrap_ci(y_rec, preds_same, preds_all, patients_rec, seed=seed)

    # Per-patient metrics
    patient_df = []
    for p in np.unique(patients_rec):
        p_idx = np.where(patients_rec == p)[0]
        if len(p_idx) > 2:
            p_r2_own = r2_score(y_rec[p_idx], preds_own[p_idx])
            p_r2_all = r2_score(y_rec[p_idx], preds_all[p_idx])
            patient_df.append({"patient_id": p, "n_cells": len(p_idx), "r2_own": p_r2_own, "r2_all": p_r2_all})

    return ModelAblationResults(
        r2_own=r2_own,
        r2_all_neighbours=r2_all,
        r2_same_type_only=r2_same,
        neighbour_gain=nbr_gain,
        cross_type_contribution=cross_contrib,
        cross_type_share_pct=cross_share,
        ci_neighbour_gain=ci_nbr,
        ci_cross_type=ci_cross,
        per_patient_r2=pd.DataFrame(patient_df),
    )

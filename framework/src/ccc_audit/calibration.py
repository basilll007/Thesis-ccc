"""Ground-truth simulation calibration for spatial communication audits.

What this module does:
Simulates synthetic planted perturbations directly onto realistic or empirical tissue geometries
to establish empirical power and false-positive rates (FPR):
1. Planted Proximity: Receptor expression in receiver cells that directly touch ligand-positive
   sender cells is scaled by factor f >= 1.0 via Poisson resampling.
2. Planted Zonation: Both ligand and receptor expression are scaled by factor f across all cells
   within random regional tiles (e.g., 200 um) regardless of contact or cell type.

Why:
Without empirical calibration, an audit cannot distinguish whether a non-significant result
stems from a true biological null or an underpowered dataset (e.g. too few receiver cells).
Conversely, calibration reveals whether a test is fooled by regional expression zonation
(reporting high false-positive rates when no contact effect exists).

Adapted from:
`src/pdac/v3_phase2_calibration.py`.
"""

from __future__ import annotations

from typing import Sequence
import numpy as np
import pandas as pd

from ccc_audit.graphs import assign_tiles
from ccc_audit.nulls import global_null, block_null


def simulate_planted_proximity(
    raw_counts_L: np.ndarray,
    raw_counts_R: np.ndarray,
    sender_mask: np.ndarray,
    receiver_mask: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    effect_factor: float,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Plant a contact-specific effect by scaling receptor in touching receivers.

    Source file: `src/pdac/v3_phase2_calibration.py`.
    """
    cL = np.array(raw_counts_L, copy=True, dtype=float)
    cR = np.array(raw_counts_R, copy=True, dtype=float)
    if effect_factor <= 1.0:
        return cL, cR

    # Senders with non-zero ligand
    pos_senders = set(np.flatnonzero(sender_mask & (cL > 0)))

    # Receivers in contact with at least one positive sender
    touching_receivers = set()
    for u, v in zip(rows, cols):
        if u in pos_senders and receiver_mask[v]:
            touching_receivers.add(v)

    if not touching_receivers:
        return cL, cR

    rng = np.random.default_rng(seed)
    touching_arr = np.array(list(touching_receivers), dtype=np.int64)
    # Resample Poisson(f * c)
    lam = effect_factor * cR[touching_arr]
    cR[touching_arr] = rng.poisson(lam)

    return cL, cR


def simulate_planted_zonation(
    raw_counts_L: np.ndarray,
    raw_counts_R: np.ndarray,
    xy: np.ndarray,
    tile_size_um: float,
    effect_factor: float,
    high_tile_fraction: float = 0.5,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Plant regional expression zonation without spatial contact dependence.

    Source file: `src/pdac/v3_phase2_calibration.py`.
    """
    cL = np.array(raw_counts_L, copy=True, dtype=float)
    cR = np.array(raw_counts_R, copy=True, dtype=float)
    if effect_factor <= 1.0:
        return cL, cR

    codes, info = assign_tiles(xy, tile_size_um, random_offset=True, seed=seed)
    unique_tiles = np.unique(codes)
    rng = np.random.default_rng(seed)

    n_high = max(1, int(len(unique_tiles) * high_tile_fraction))
    high_tiles = set(rng.choice(unique_tiles, size=n_high, replace=False))

    high_cells = np.array([i for i, c in enumerate(codes) if c in high_tiles], dtype=np.int64)
    if len(high_cells) > 0:
        cL[high_cells] = rng.poisson(effect_factor * cL[high_cells])
        cR[high_cells] = rng.poisson(effect_factor * cR[high_cells])

    return cL, cR


def run_calibration_simulation(
    xy: np.ndarray,
    labels: np.ndarray,
    raw_counts_L: np.ndarray,
    raw_counts_R: np.ndarray,
    sender_type: int,
    receiver_type: int,
    rows: np.ndarray,
    cols: np.ndarray,
    effects: Sequence[float] = (1.0, 1.1, 1.25, 1.5, 2.0),
    n_sims: int = 50,
    n_perms: int = 1000,
    tile_size_um: float = 100.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Run calibration simulations comparing power and false-positive rates.

    Source file: `src/pdac/v3_phase2_calibration.py`.
    """
    sender_mask = labels == sender_type
    receiver_mask = labels == receiver_type
    total_counts = np.maximum(raw_counts_L + raw_counts_R, 1.0)

    # Directed edges connecting sender -> receiver
    e_obs_mask = (labels[rows] == sender_type) & (labels[cols] == receiver_type)
    e_indices = np.flatnonzero(e_obs_mask)

    results = []

    # Pre-generate global and block permutations
    g_perms = global_null(labels, n_perms=n_perms, seed=seed)
    b_perms, _ = block_null(labels, xy, tile_size_um=tile_size_um, n_perms=n_perms, seed=seed)

    null_dict = {"global": g_perms, "block": b_perms}

    for scenario in ("proximity", "zonation"):
        for f in effects:
            sig_counts = {"global": 0, "block": 0}
            valid_sims = 0

            for sim_i in range(n_sims):
                sim_seed = seed + sim_i * 1000 + int(f * 10)
                if scenario == "proximity":
                    cL, cR = simulate_planted_proximity(
                        raw_counts_L, raw_counts_R, sender_mask, receiver_mask,
                        rows, cols, effect_factor=f, seed=sim_seed,
                    )
                else:
                    cL, cR = simulate_planted_zonation(
                        raw_counts_L, raw_counts_R, xy, tile_size_um=tile_size_um,
                        effect_factor=f, seed=sim_seed,
                    )

                # Log-normalise
                with np.errstate(invalid="ignore", divide="ignore"):
                    L = np.log1p(cL / total_counts * 1e4)
                    R = np.log1p(cR / total_counts * 1e4)
                L = np.nan_to_num(L)
                R = np.nan_to_num(R)

                # Observed Q
                Lm = L[sender_mask].mean() if sender_mask.any() else 0.0
                Rm = R[receiver_mask].mean() if receiver_mask.any() else 0.0
                M = Lm * Rm
                if len(e_indices) == 0 or M <= 0:
                    continue

                E = np.mean(L[rows[e_indices]] * R[cols[e_indices]])
                Q_obs = E / M
                valid_sims += 1

                # Evaluate under nulls
                for n_name, perms in null_dict.items():
                    null_Q = []
                    for pl in perms:
                        s_m = pl == sender_type
                        r_m = pl == receiver_type
                        s_mean = L[s_m].mean() if s_m.any() else 0.0
                        r_mean = R[r_m].mean() if r_m.any() else 0.0
                        null_M = s_mean * r_mean

                        edge_m = (pl[rows] == sender_type) & (pl[cols] == receiver_type)
                        if not edge_m.any() or null_M <= 0:
                            null_Q.append(0.0)
                        else:
                            null_E = np.mean(L[rows[edge_m]] * R[cols[edge_m]])
                            null_Q.append(null_E / null_M)

                    null_Q = np.array(null_Q)
                    # Add-one p-value
                    p_val = (np.sum(null_Q >= Q_obs) + 1) / (len(null_Q) + 1)
                    if p_val < 0.05:
                        sig_counts[n_name] += 1

            for n_name in ("global", "block"):
                rate = sig_counts[n_name] / max(valid_sims, 1)
                results.append({
                    "scenario": scenario,
                    "effect_factor": f,
                    "null_model": n_name,
                    "rejection_rate": rate,
                    "metric_type": "power" if scenario == "proximity" else "false_positive_rate",
                    "n_sims_valid": valid_sims,
                })

    return pd.DataFrame(results)

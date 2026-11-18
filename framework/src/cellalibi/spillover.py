"""Segmentation spillover diagnostics and contamination controls.

What this module does:
Diagnoses optical/segmentation transcript spillover between adjacent cell types:
1. Identifies contamination marker genes (M) predominantly expressed in the sender population.
2. Partitions receiver cells into distance bins based on median cell diameter d (touching <= d,
   near d < dist <= 2d, separated dist > 2d).
3. Computes a per-cell contamination covariate s_c.
4. Identifies shared genes (CAF-shared) that should be removed from receiver targets.
5. Identifies spillover-free negative-control receiver genes (S) that are absent or negligible in senders.
6. Computes a mixing expectation model comparing observed transcript levels to pure optical mixing.

Why:
In high-plex spatial imaging (CosMx, Xenium, MERSCOPE), cell boundaries are imperfectly segmented.
A receiver cell touching a sender cell frequently absorbs transcripts from the sender's cytoplasm.
This produces a spurious correlation between sender and receiver state that mimics genuine juxtacrine
signaling. If an apparent communication effect is confined to touching cells, vanishes upon including
the contamination covariate, or matches negative control spillover rates, it is an optical artefact.

Adapted from:
`model_dev/src/v2_phase1.py` and `model_dev/reports/MODEL_REPORT_v2.md`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


@dataclass
class SpilloverDiagnostics:
    """Spillover diagnostics summary and calibrated gene sets."""
    median_cell_diameter_um: float
    distance_bins: pd.Series  # 'touching', 'near', 'separated'
    distances_to_sender_um: np.ndarray
    contamination_markers: list[str]
    contamination_covariate: np.ndarray
    shared_genes_to_remove: list[str]
    spillover_free_genes: list[str]
    mixing_excess_pct: float
    is_explained_by_spillover: bool
    diagnostic_details: dict[str, float]


def estimate_median_cell_diameter(xy: np.ndarray) -> float:
    """Estimate median cell diameter from positive nearest-neighbour distances.

    Source file: `model_dev/src/v2_phase1.py`.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.shape[0] <= 1:
        return 12.0  # default biological cell diameter in um

    nn = NearestNeighbors(n_neighbors=2, metric="euclidean").fit(xy)
    dists, _ = nn.kneighbors(xy)
    # 1st column is distance to self (0), 2nd is nearest neighbour
    nn_dist = dists[:, 1]
    pos_dists = nn_dist[nn_dist > 0]
    if len(pos_dists) == 0:
        return 12.0
    return float(np.median(pos_dists))


def compute_distances_to_senders(
    xy_receivers: np.ndarray,
    xy_senders: np.ndarray,
) -> np.ndarray:
    """Compute Euclidean distance in micrometres from each receiver to the nearest sender.

    Source file: `model_dev/src/v2_phase1.py`.
    """
    if len(xy_senders) == 0:
        return np.full(len(xy_receivers), np.inf, dtype=float)

    nn = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(xy_senders)
    dists, _ = nn.kneighbors(xy_receivers)
    return dists.ravel().astype(float)


def bin_by_distance(
    distances_um: np.ndarray,
    cell_diameter_um: float,
) -> pd.Series:
    """Classify receiver cells into touching, near, and separated bins.

    Bins:
    - touching: dist <= cell_diameter_um (1 median diameter)
    - near: cell_diameter_um < dist <= 2 * cell_diameter_um
    - separated: dist > 2 * cell_diameter_um
    """
    d = np.asarray(distances_um, dtype=float)
    bins = pd.Series(index=range(len(d)), dtype="object")
    bins[d <= cell_diameter_um] = "touching"
    bins[(d > cell_diameter_um) & (d <= 2 * cell_diameter_um)] = "near"
    bins[d > 2 * cell_diameter_um] = "separated"
    return bins


def identify_contamination_markers(
    expression_df: pd.DataFrame,
    sender_mask: np.ndarray,
    isolated_receiver_mask: np.ndarray,
    min_ratio: float = 10.0,
    fallback_ratio: float = 5.0,
    top_n: int = 20,
) -> list[str]:
    """Identify marker genes with high expression in senders relative to isolated receivers.

    Source file: `model_dev/src/v2_phase1.py` (gene_sets_and_adjusted_targets).
    """
    sender_means = expression_df.iloc[sender_mask].mean(axis=0)
    receiver_means = expression_df.iloc[isolated_receiver_mask].mean(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = sender_means / np.maximum(receiver_means, 1e-6)

    # Filter by min_ratio, fallback if fewer than 5 genes qualify
    candidates = ratios[ratios >= min_ratio].sort_values(ascending=False)
    if len(candidates) < 5:
        candidates = ratios[ratios >= fallback_ratio].sort_values(ascending=False)

    return list(candidates.head(top_n).index)


def identify_spillover_free_genes(
    expression_df: pd.DataFrame,
    sender_mask: np.ndarray,
    isolated_receiver_mask: np.ndarray,
    min_receiver_to_sender_ratio: float = 10.0,
    top_n: int = 5,
) -> list[str]:
    """Identify negative-control genes with high expression in receivers and near-zero in senders.

    Source file: `model_dev/src/v2_phase1.py`.
    """
    sender_means = expression_df.iloc[sender_mask].mean(axis=0)
    receiver_means = expression_df.iloc[isolated_receiver_mask].mean(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = receiver_means / np.maximum(sender_means, 1e-6)

    candidates = ratios[ratios >= min_receiver_to_sender_ratio].sort_values(ascending=False)
    return list(candidates.head(top_n).index)


def identify_shared_genes_to_remove(
    expression_df: pd.DataFrame,
    target_genes: Sequence[str],
    sender_mask: np.ndarray,
    isolated_receiver_mask: np.ndarray,
    threshold_ratio: float = 0.5,
) -> list[str]:
    """Identify target program genes that have substantial sender expression (CAF-shared).

    Source file: `model_dev/src/v2_phase1.py`.
    """
    valid_targets = [g for g in target_genes if g in expression_df.columns]
    sender_means = expression_df[valid_targets].iloc[sender_mask].mean(axis=0)
    receiver_means = expression_df[valid_targets].iloc[isolated_receiver_mask].mean(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = sender_means / np.maximum(receiver_means, 1e-6)

    shared = ratios[ratios >= threshold_ratio].index.tolist()
    return shared


def compute_mixing_expectation(
    marker_touching_mean: float,
    marker_isolated_mean: float,
    marker_sender_mean: float,
) -> dict[str, float]:
    """Compute physical transcript mixing expectation and excess contamination rate.

    Source file: `model_dev/reports/MODEL_REPORT_v2.md` and `b_spillover_rate.csv`.
    """
    excess = max(0.0, marker_touching_mean - marker_isolated_mean)
    sender_spread = max(1e-6, marker_sender_mean - marker_isolated_mean)
    excess_pct = (excess / sender_spread) * 100.0

    return {
        "marker_touching_mean": float(marker_touching_mean),
        "marker_isolated_mean": float(marker_isolated_mean),
        "marker_sender_mean": float(marker_sender_mean),
        "excess_level": float(excess),
        "excess_pct": float(excess_pct),
    }


def audit_spillover(
    xy: np.ndarray,
    expression_df: pd.DataFrame,
    sender_mask: np.ndarray,
    receiver_mask: np.ndarray,
    target_genes: Sequence[str],
) -> SpilloverDiagnostics:
    """Run complete spillover audit across distance bins and contamination models.

    Source file: `model_dev/src/v2_phase1.py`.
    """
    diameter = estimate_median_cell_diameter(xy)
    xy_rec = xy[receiver_mask]
    xy_send = xy[sender_mask]

    dists = compute_distances_to_senders(xy_rec, xy_send)
    bins = bin_by_distance(dists, diameter)

    # Isolated receivers are in separated bin (> 2 cell diameters)
    separated_idx = np.where(bins == "separated")[0]
    if len(separated_idx) == 0:
        # Fallback to the upper quartile of distances
        cutoff = np.percentile(dists, 75)
        isolated_rec_mask_full = receiver_mask.copy()
        rec_indices = np.where(receiver_mask)[0]
        isolated_rec_mask_full[rec_indices[dists < cutoff]] = False
    else:
        isolated_rec_mask_full = np.zeros(len(xy), dtype=bool)
        rec_indices = np.where(receiver_mask)[0]
        isolated_rec_mask_full[rec_indices[separated_idx]] = True

    # Identify contamination markers M and spillover-free genes S
    M_genes = identify_contamination_markers(expression_df, sender_mask, isolated_rec_mask_full)
    S_genes = identify_spillover_free_genes(expression_df, sender_mask, isolated_rec_mask_full)
    shared_to_remove = identify_shared_genes_to_remove(expression_df, target_genes, sender_mask, isolated_rec_mask_full)

    # Contamination covariate s_c
    if M_genes:
        s_c = expression_df[M_genes].iloc[receiver_mask].mean(axis=1).to_numpy()
    else:
        s_c = np.zeros(int(receiver_mask.sum()), dtype=float)

    # Mixing expectation
    if M_genes:
        touching_idx = np.where(bins == "touching")[0]
        touching_mean = float(expression_df[M_genes].iloc[receiver_mask].iloc[touching_idx].mean().mean()) if len(touching_idx) else 0.0
        isolated_mean = float(expression_df[M_genes].iloc[isolated_rec_mask_full].mean().mean())
        sender_mean = float(expression_df[M_genes].iloc[sender_mask].mean().mean())
        mix = compute_mixing_expectation(touching_mean, isolated_mean, sender_mean)
    else:
        mix = {"excess_pct": 0.0}

    # Spillover rule: if excess > 10% and touching shows elevation over separated
    is_spillover = mix["excess_pct"] >= 10.0

    return SpilloverDiagnostics(
        median_cell_diameter_um=diameter,
        distance_bins=bins,
        distances_to_sender_um=dists,
        contamination_markers=M_genes,
        contamination_covariate=s_c,
        shared_genes_to_remove=shared_to_remove,
        spillover_free_genes=S_genes,
        mixing_excess_pct=mix["excess_pct"],
        is_explained_by_spillover=is_spillover,
        diagnostic_details=mix,
    )

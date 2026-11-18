"""Input/Output and input contract validation for spatial transcriptomics data.

What this module does:
Validates and standardises AnnData objects representing spatial single-cell transcriptomics
data against the strict contract required by the cellalibi framework.

Why:
Spatial audits require unambiguous coordinates (in physical units, micrometres) to compute
accurate contact graphs and spatial bins, plus explicit hierarchy columns (patient_id, section_id)
so that grouping and cross-validation never leak patient-level dependencies across splits.
Failing fast with clear, actionable diagnostics prevents silent corruption of spatial analyses.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import anndata as ad
import numpy as np
import pandas as pd

REQUIRED_OBS_COLUMNS: tuple[str, ...] = (
    "cell_type",
    "patient_id",
    "section_id",
    "x_um",
    "y_um",
)


class DataContractError(ValueError):
    """Raised when an input AnnData violates the cellalibi input contract."""
    pass


def validate_spatial_data(adata: ad.AnnData) -> dict[str, Any]:
    """Validate that an AnnData meets the cellalibi input contract.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix containing single cells.

    Returns
    -------
    dict[str, Any]
        Metadata summary of validated dataset, including detected optional modalities.

    Raises
    ------
    DataContractError
        If required columns are missing, data contains invalid coordinates, or cells are zero.
    """
    if adata.n_obs == 0:
        raise DataContractError("Input AnnData contains 0 cells (n_obs == 0).")
    if adata.n_vars == 0:
        raise DataContractError("Input AnnData contains 0 variables/genes (n_vars == 0).")

    # Check required obs columns
    missing_cols = [c for c in REQUIRED_OBS_COLUMNS if c not in adata.obs.columns]
    if missing_cols:
        raise DataContractError(
            f"Input AnnData.obs is missing required columns: {missing_cols}. "
            f"The cellalibi contract requires {list(REQUIRED_OBS_COLUMNS)}."
        )

    # Check physical coordinates
    for coord in ("x_um", "y_um"):
        s = adata.obs[coord]
        if not pd.api.types.is_numeric_dtype(s):
            raise DataContractError(f"Coordinate column '{coord}' must be numeric, got {s.dtype}.")
        if s.isna().any():
            n_nan = int(s.isna().sum())
            raise DataContractError(f"Coordinate column '{coord}' contains {n_nan} NaN/missing values.")
        if np.isinf(s.to_numpy()).any():
            raise DataContractError(f"Coordinate column '{coord}' contains infinite values.")

    # Check cell_type, patient_id, section_id
    for col in ("cell_type", "patient_id", "section_id"):
        s = adata.obs[col]
        if s.isna().any():
            n_nan = int(s.isna().sum())
            raise DataContractError(f"Metadata column '{col}' contains {n_nan} missing values.")
        if (s.astype(str).str.strip() == "").any():
            raise DataContractError(f"Metadata column '{col}' contains empty strings.")

    # Detect optional negative control probes
    var_names = [str(v) for v in adata.var_names]
    neg_probe_prefixes = ("NegPrb", "Negative", "neg_probe", "Blank-", "SystemControl")
    detected_neg_probes = [
        v for v in var_names
        if any(v.startswith(prefix) or prefix.lower() in v.lower() for prefix in neg_probe_prefixes)
    ]

    # Detect optional protein channels in obs or obsm
    protein_channels: list[str] = []
    if "protein" in adata.obsm:
        protein_channels = [f"obsm['protein']:{i}" for i in range(adata.obsm["protein"].shape[1])]
    for col in adata.obs.columns:
        if any(marker in col.lower() for marker in ("cd45", "dapi", "panck", "cd3", "cd68", "sma", "protein")):
            protein_channels.append(col)

    metadata = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "cell_types": sorted(set(adata.obs["cell_type"].astype(str))),
        "n_cell_types": int(adata.obs["cell_type"].nunique()),
        "patients": sorted(set(adata.obs["patient_id"].astype(str))),
        "n_patients": int(adata.obs["patient_id"].nunique()),
        "sections": sorted(set(adata.obs["section_id"].astype(str))),
        "n_sections": int(adata.obs["section_id"].nunique()),
        "negative_control_probes": detected_neg_probes,
        "n_negative_control_probes": len(detected_neg_probes),
        "protein_channels": protein_channels,
        "x_range_um": (float(adata.obs["x_um"].min()), float(adata.obs["x_um"].max())),
        "y_range_um": (float(adata.obs["y_um"].min()), float(adata.obs["y_um"].max())),
    }
    return metadata


def load_spatial_dataset(
    file_path: str | Path,
    obs_mapping: dict[str, str] | None = None,
) -> ad.AnnData:
    """Load and validate an AnnData from disk.

    Parameters
    ----------
    file_path : str or Path
        Path to an `.h5ad` file.
    obs_mapping : dict[str, str], optional
        Optional rename mapping for `.obs` columns to meet required names
        (e.g., `{'Center_X': 'x_um', 'Center_Y': 'y_um'}`).

    Returns
    -------
    ad.AnnData
        Standardised, validated AnnData.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Spatial dataset not found: {path}")

    adata = ad.read_h5ad(path)
    if obs_mapping:
        adata.obs = adata.obs.rename(columns=obs_mapping)

    validate_spatial_data(adata)
    return adata

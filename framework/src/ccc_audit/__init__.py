"""ccc_audit: Spatial cell-cell communication audit and confounding diagnosis framework.

This package turns the thesis spatial analyses into a reusable, modular audit tool to determine
whether putative cell-cell communication signals are genuine contact effects or artefacts
of abundance, regional zonation, or segmentation spillover.
"""

__version__ = "0.1.0"

from ccc_audit.io import validate_spatial_data, load_spatial_dataset
from ccc_audit.graphs import knn_graph, radius_graph, assign_tiles
from ccc_audit.scores import compute_ccc_scores, LRSetup, ScoreResults
from ccc_audit.nulls import global_null, block_null, expression_tile_null
from ccc_audit.calibration import run_calibration_simulation
from ccc_audit.spillover import SpilloverDiagnostics
from ccc_audit.neighbour_models import evaluate_neighbour_prediction
from ccc_audit.verdicts import AuditVerdict, decide_verdict, generate_claims_register

__all__ = [
    "validate_spatial_data",
    "load_spatial_dataset",
    "knn_graph",
    "radius_graph",
    "assign_tiles",
    "compute_ccc_scores",
    "LRSetup",
    "ScoreResults",
    "global_null",
    "block_null",
    "expression_tile_null",
    "run_calibration_simulation",
    "SpilloverDiagnostics",
    "evaluate_neighbour_prediction",
    "AuditVerdict",
    "decide_verdict",
    "generate_claims_register",
]

"""cellalibi: Spatial cell-cell communication audit and confounding diagnosis framework.

This package turns the thesis spatial analyses into a reusable, modular audit tool to determine
whether putative cell-cell communication signals are genuine contact effects or artefacts
of abundance, regional zonation, or segmentation spillover.

Previously named `ccc_audit` (v0.1.0, branch framework-v1). Renamed to `cellalibi` in v0.2.0 (framework-v2).
"""

__version__ = "0.2.1"

from cellalibi.io import validate_spatial_data, load_spatial_dataset
from cellalibi.graphs import knn_graph, radius_graph, assign_tiles
from cellalibi.scores import compute_ccc_scores, LRSetup, ScoreResults
from cellalibi.nulls import global_null, block_null, expression_tile_null
from cellalibi.calibration import run_calibration_simulation
from cellalibi.spillover import SpilloverDiagnostics
from cellalibi.neighbour_models import evaluate_neighbour_prediction
from cellalibi.verdicts import AuditVerdict, decide_verdict, generate_claims_register
from cellalibi.verdicts import CHECKS, PRIORITY_RULE, ClaimAudit, audit_claim, claims_table
from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext, contact_statistics

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
    "CHECKS",
    "PRIORITY_RULE",
    "ClaimAudit",
    "audit_claim",
    "claims_table",
    "AuditConfig",
    "ExpressionData",
    "SpatialContext",
    "contact_statistics",
]

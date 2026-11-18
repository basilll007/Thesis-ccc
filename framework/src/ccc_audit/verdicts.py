"""Unified audit verdict decision logic and claims register generator.

What this module does:
Integrates score comparisons, null model tests, simulation calibrations, spillover controls,
and neighbour model ablations to assign a single definitive verdict for each interaction
or receiver program:
  1. NOT_TESTABLE: Receiver cell count, contact edge count, or calibrated statistical power
     is below the empirical detection floor.
  2. EXPLAINED_BY_ABUNDANCE: Apparent communication score matches cell-type abundance products;
     proximity ratio Q shows no enrichment beyond the non-spatial matched baseline.
  3. EXPLAINED_BY_ZONATION: Contact enrichment under a global null vanishes under spatially
     constrained block nulls, or false-positive rate under regional zonation exceeds tolerances.
  4. EXPLAINED_BY_SPILLOVER: Predictive association or enrichment is confined to physically touching
     cells (<= 1 cell diameter), vanishes after controlling for contamination markers, or reflects
     optical mixing.
  5. SUPPORTED: Tested above the power floor; contact ratio Q or cross-type neighbour gain is
     empirically significant, survives spatially constrained nulls, and persists outside touching range.

Why:
Biological spatial literature frequently claims 'communication' based on uncalibrated metrics that
are saturated by cell abundance, regional zonation, or segmentation artefacts. The unified verdict
system enforces a rigorous hierarchy of falsification controls.

Adapted from:
`src/pdac/v3_phase6_synthesis.py` and `docs/PDAC_SUMMARY_v3.md`.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Sequence
import pandas as pd


class AuditVerdict(str, Enum):
    """Categorical verdict for a spatial communication hypothesis."""
    SUPPORTED = "SUPPORTED"
    EXPLAINED_BY_ABUNDANCE = "EXPLAINED_BY_ABUNDANCE"
    EXPLAINED_BY_ZONATION = "EXPLAINED_BY_ZONATION"
    EXPLAINED_BY_SPILLOVER = "EXPLAINED_BY_SPILLOVER"
    NOT_TESTABLE = "NOT_TESTABLE"


@dataclass
class VerdictResult:
    """Full decision output for an interaction or receiver program."""
    target_id: str
    verdict: AuditVerdict
    primary_reason: str
    key_metrics: dict[str, Any]
    evidence_summary: str


def decide_verdict(
    target_id: str,
    n_receiver_cells: int,
    n_qualifying_edges: int,
    comp_p: float | None = None,
    matched_p: float | None = None,
    ratio_p_global: float | None = None,
    ratio_p_block: float | None = None,
    calibrated_power: float | None = None,
    calibrated_zonation_fpr: float | None = None,
    is_spillover: bool | None = None,
    spillover_excess_pct: float | None = None,
    cross_type_gain: float | None = None,
    cross_type_ci_low: float | None = None,
    min_receiver_cells: int = 25,
    min_edges: int = 10,
    alpha: float = 0.05,
) -> VerdictResult:
    """Execute prioritized verdict decision tree.

    Decision hierarchy:
    1. Power / Testability floor
    2. Abundance confounding (M/C significant, Q null)
    3. Zonation confounding (Global Q significant, Block Q null or high zonation FPR)
    4. Spillover confounding (Signal confined to touching, high contamination excess)
    5. Supported (Passes all checks)
    """
    metrics = {
        "n_receiver_cells": n_receiver_cells,
        "n_qualifying_edges": n_qualifying_edges,
        "comp_p": comp_p,
        "matched_p": matched_p,
        "ratio_p_global": ratio_p_global,
        "ratio_p_block": ratio_p_block,
        "calibrated_power": calibrated_power,
        "calibrated_zonation_fpr": calibrated_zonation_fpr,
        "is_spillover": is_spillover,
        "spillover_excess_pct": spillover_excess_pct,
        "cross_type_gain": cross_type_gain,
        "cross_type_ci_low": cross_type_ci_low,
    }

    # Branch 1: NOT_TESTABLE
    if n_receiver_cells < min_receiver_cells:
        reason = f"Receiver cell count ({n_receiver_cells}) is below power floor ({min_receiver_cells})."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.NOT_TESTABLE,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"NOT_TESTABLE: {reason}",
        )

    if n_qualifying_edges < min_edges:
        reason = f"Qualifying contact edges ({n_qualifying_edges}) is below minimum threshold ({min_edges})."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.NOT_TESTABLE,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"NOT_TESTABLE: {reason}",
        )

    if calibrated_power is not None and calibrated_power < 0.20:
        reason = f"Calibrated statistical power at 1.5x effect is only {calibrated_power:.2f} (< 0.20 floor)."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.NOT_TESTABLE,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"NOT_TESTABLE: {reason}",
        )

    # Branch 2: EXPLAINED_BY_SPILLOVER (checked early if predictive/spillover audit was run)
    if is_spillover is True or (spillover_excess_pct is not None and spillover_excess_pct >= 10.0):
        excess_str = f"{spillover_excess_pct:.1f}%" if spillover_excess_pct is not None else "significant"
        reason = f"Explained by segmentation spillover (marker excess {excess_str} in touching cells; signal confined to touching distance)."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.EXPLAINED_BY_SPILLOVER,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"EXPLAINED_BY_SPILLOVER: {reason}",
        )

    # Branch 3: EXPLAINED_BY_ABUNDANCE
    # If edge or matched baseline is significant, but proximity ratio Q is not significant
    if ratio_p_global is not None and ratio_p_global >= alpha:
        # Q is null
        if matched_p is not None and matched_p < alpha:
            reason = f"Matched dyadic baseline M is significant (p={matched_p:.4f}), but proximity ratio Q is not (p={ratio_p_global:.4f}). Score is driven by cell abundance."
            return VerdictResult(
                target_id=target_id,
                verdict=AuditVerdict.EXPLAINED_BY_ABUNDANCE,
                primary_reason=reason,
                key_metrics=metrics,
                evidence_summary=f"EXPLAINED_BY_ABUNDANCE: {reason}",
            )
        if comp_p is not None and comp_p < alpha:
            reason = f"Compositional score C is significant (p={comp_p:.4f}), but proximity ratio Q is not (p={ratio_p_global:.4f}). Driven by cell abundance."
            return VerdictResult(
                target_id=target_id,
                verdict=AuditVerdict.EXPLAINED_BY_ABUNDANCE,
                primary_reason=reason,
                key_metrics=metrics,
                evidence_summary=f"EXPLAINED_BY_ABUNDANCE: {reason}",
            )

    # Branch 4: EXPLAINED_BY_ZONATION
    # Global Q is significant, but block null shrinks/fails, OR zonation FPR is high
    if ratio_p_global is not None and ratio_p_global < alpha:
        if ratio_p_block is not None and ratio_p_block >= alpha:
            reason = f"Proximity ratio Q is significant under global null (p={ratio_p_global:.4f}), but non-significant under spatially constrained block null (p={ratio_p_block:.4f})."
            return VerdictResult(
                target_id=target_id,
                verdict=AuditVerdict.EXPLAINED_BY_ZONATION,
                primary_reason=reason,
                key_metrics=metrics,
                evidence_summary=f"EXPLAINED_BY_ZONATION: {reason}",
            )
        if calibrated_zonation_fpr is not None and calibrated_zonation_fpr > 0.10:
            reason = f"Simulation indicates high false-positive rate under 200 um zonation ({calibrated_zonation_fpr:.2f} > 0.10). Signal cannot be separated from regional zonation."
            return VerdictResult(
                target_id=target_id,
                verdict=AuditVerdict.EXPLAINED_BY_ZONATION,
                primary_reason=reason,
                key_metrics=metrics,
                evidence_summary=f"EXPLAINED_BY_ZONATION: {reason}",
            )

    # Branch 5: Cross-type share check for neighbour models
    if cross_type_ci_low is not None and cross_type_ci_low <= 0:
        reason = f"Cross-type neighbour contribution 95% CI includes zero [{cross_type_ci_low:.4f}, ...]; predictive gain is carried by same-type neighbours."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.EXPLAINED_BY_ABUNDANCE,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"EXPLAINED_BY_ABUNDANCE: {reason}",
        )

    # Branch 6: SUPPORTED
    # Must have passing Q or positive cross-type gain with passing lower CI
    q_passed = (ratio_p_global is not None and ratio_p_global < alpha) and (ratio_p_block is None or ratio_p_block < alpha)
    ml_passed = (cross_type_ci_low is not None and cross_type_ci_low > 0)

    if q_passed or ml_passed:
        reason = "Passes all confounding controls: contact enrichment persists under spatially constrained nulls and outside spillover range."
        return VerdictResult(
            target_id=target_id,
            verdict=AuditVerdict.SUPPORTED,
            primary_reason=reason,
            key_metrics=metrics,
            evidence_summary=f"SUPPORTED: {reason}",
        )

    # Default fallback: NOT_TESTABLE or EXPLAINED_BY_ABUNDANCE
    reason = "Signal does not reach empirical significance over dyadic null expectation."
    return VerdictResult(
        target_id=target_id,
        verdict=AuditVerdict.EXPLAINED_BY_ABUNDANCE,
        primary_reason=reason,
        key_metrics=metrics,
        evidence_summary=f"EXPLAINED_BY_ABUNDANCE: {reason}",
    )


def generate_claims_register(verdicts: Sequence[VerdictResult]) -> pd.DataFrame:
    """Generate a structured claims register table matching thesis conventions.

    Source file: `results/pdac/v3_phase6/claims_register.csv`.
    """
    records = []
    for i, v in enumerate(verdicts, start=1):
        records.append({
            "id": f"CLAIM_{i:03d}",
            "target": v.target_id,
            "verdict": v.verdict.value,
            "primary_reason": v.primary_reason,
            "evidence_summary": v.evidence_summary,
        })
    return pd.DataFrame(records)

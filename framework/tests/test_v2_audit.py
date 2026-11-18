"""Unit tests for the v2 claim audit (multi-label checks + pre-declared priority rule)."""
import numpy as np
import pytest

from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
from cellalibi.simulate import synthetic_tissue
from cellalibi.verdicts import AuditVerdict, CHECKS, audit_claim, claims_table, primary_verdict

CFG = AuditConfig(n_perms=1000)


def _checks(**kw):
    c = {k: "PASS" for k in CHECKS}
    c.update(kw)
    return c


@pytest.mark.parametrize("checks,m,expected", [
    (_checks(power="FAIL", spillover="FAIL", abundance="FAIL"), {"p_M_label": 0.001}, AuditVerdict.NOT_TESTABLE),
    (_checks(spillover="FAIL", abundance="FAIL", zonation="FAIL"), {"p_M_label": 0.001}, AuditVerdict.EXPLAINED_BY_SPILLOVER),
    (_checks(abundance="FAIL", zonation="FAIL"), {"p_M_label": 0.001}, AuditVerdict.EXPLAINED_BY_ABUNDANCE),
    (_checks(abundance="FAIL"), {"p_M_label": 0.5}, AuditVerdict.NOT_SUPPORTED),
    (_checks(zonation="FAIL"), {"p_M_label": 0.001}, AuditVerdict.EXPLAINED_BY_ZONATION),
    (_checks(zonation="NOT_RUN"), {"p_M_label": 0.001}, AuditVerdict.NOT_TESTABLE),
    (_checks(spillover="NOT_RUN"), {"p_M_label": 0.001}, AuditVerdict.NOT_TESTABLE),
    (_checks(), {"p_M_label": 0.001}, AuditVerdict.SUPPORTED),
])
def test_priority_rule_every_branch(checks, m, expected):
    v, _ = primary_verdict(checks, m, CFG)
    assert v == expected


ENGINES = ["0.2.0", "0.2.1"]


def _ctx(seed, engine="0.2.1", **kw):
    xy, labels, counts, genes = synthetic_tissue(seed=seed, n_cells=8000, size_um=1800.0, **kw)
    return SpatialContext(xy, labels, AuditConfig(n_perms=1000, seed=seed, engine=engine)), ExpressionData(counts, genes)


@pytest.mark.parametrize("engine", ENGINES)
def test_planted_contact_is_supported(engine):
    ctx, ex = _ctx(7, engine=engine, contact_effect=3.0)
    au = audit_claim(ctx, ex, "Sender", "Receiver", "L", "R")
    assert au.checks == {"power": "PASS", "abundance": "PASS", "zonation": "PASS", "spillover": "PASS"}, au.check_reasons
    assert au.primary == AuditVerdict.SUPPORTED


@pytest.mark.parametrize("engine", ENGINES)
def test_null_is_not_supported_and_all_checks_reported(engine):
    ctx, ex = _ctx(8, engine=engine)
    au = audit_claim(ctx, ex, "Sender", "Receiver", "L", "R")
    assert au.primary != AuditVerdict.SUPPORTED
    assert set(au.checks) == set(CHECKS)
    assert all(v in ("PASS", "FAIL", "NOT_RUN") for v in au.checks.values())


@pytest.mark.parametrize("engine", ENGINES)
def test_spillover_flagged_and_claim_specific(engine):
    ctx, ex = _ctx(9, engine=engine, spill_frac=0.3)
    au = audit_claim(ctx, ex, "Sender", "Receiver", "L", "R")
    assert au.checks["spillover"] == "FAIL", au.check_reasons["spillover"]
    assert au.primary == AuditVerdict.EXPLAINED_BY_SPILLOVER
    # a claim whose sender does not spill must not inherit the spillover flag (v1 bug)
    au2 = audit_claim(ctx, ex, "Other", "Receiver", "L", "R")
    assert au2.checks["spillover"] != "FAIL", au2.check_reasons["spillover"]
    tab = claims_table([au, au2])
    assert list(tab["check_spillover"]) == [au.checks["spillover"], au2.checks["spillover"]]


def test_min_attainable_p_reported():
    ctx, ex = _ctx(10)
    au = audit_claim(ctx, ex, "Sender", "Receiver", "L", "R")
    assert au.metrics["n_perms"] == 1000
    assert abs(au.metrics["min_attainable_p"] - 1 / 1001) < 1e-12


def test_knn_graph_keeps_nearest_neighbour():
    """Regression test for the v1 bug: knn_graph dropped each cell's true nearest neighbour."""
    from cellalibi.graphs import knn_graph
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0], [6.0, 0.0], [10.0, 0.0]])
    r, c = knn_graph(xy, k=1)
    assert list(c) == [1, 0, 1, 2, 3]
    assert not np.any(r == c)



def test_v021_power_plant_reaches_zero_count_cells():
    """v0.2.0 bug: multiplying observed counts leaves zero-count cells at zero. v0.2.1 plants on the mean."""
    from cellalibi.contact import _planted_receptor_logs
    rng = np.random.default_rng(0)
    counts = np.zeros((400, 2))
    counts[:, 0] = 50
    counts[:200, 1] = rng.poisson(0.3, 200)   # sparse receptor
    counts[200:, 1] = rng.poisson(0.3, 200)
    ex = ExpressionData(counts, ["X", "R"])
    R = np.arange(400)
    exposed = np.zeros(400, bool)
    exposed[:200] = True
    zero_exposed = exposed & (counts[:, 1] == 0)
    plants = _planted_receptor_logs(ex, "R", R, exposed, 2.0, AuditConfig(n_plants=20), rng)
    gained = np.mean([np.mean(p[zero_exposed] > 0) for p in plants])
    assert gained > 0.1
    base = ex.gene("R")[0]
    ratio = np.mean([np.expm1(p[exposed]).mean() for p in plants]) / base[exposed].mean()
    assert 1.6 < ratio < 2.4

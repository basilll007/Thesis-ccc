# v2 Phase 6 plan, amendment 1 (committed before any Phase 6 audit)

Date 2026-11-18. Amends `v2_phase6.md` ("engine = frozen Phase 4 engine v0.2.0").

## Reason

The pre-registered Phase 4 benchmark showed that v0.2.0 cannot return a decision on CosMx-like data:
- its 10x marker rule left the spillover check NOT_RUN in 100 % of T6-b and U5-a replicates;
- its power estimate was biased low for sparse genes;
- across 4,200 audits it gave 0 SUPPORTED, and 67-100 % NOT_TESTABLE per setting.

Running v0.2.0 alone on real data would mostly return NOT_TESTABLE because of these defects, not because of the data.

## Change

- **Primary engine: v0.2.1** (plan `v2_phase4b.md`, commit 2720a9e), whatever its Phase 4b benchmark result. That
  result is quoted next to the real-data verdicts so readers can judge them.
- **Secondary: v0.2.0** (`AuditConfig(engine="0.2.0")`) on a headline subset, to show what the frozen engine would
  have said. The subset is:
  - CAF_fibroblast->T_cell CXCL12-CXCR4;
  - CAF_fibroblast->epithelial_ductal_tumor CXCL12-CXCR4;
  - CAF_fibroblast->myeloid FN1-ITGB1;
  - myeloid->epithelial_ductal_tumor SPP1-CD44;
  - all 8 receiver-program claims.

Everything else in `v2_phase6.md` is unchanged: data, claims, aggregation, comparisons, runtime logging. Runtime and
peak memory are recorded for the whole section job, which includes both engines.

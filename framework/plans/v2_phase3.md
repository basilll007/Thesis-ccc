# v2 Phase 3 plan: rename `ccc_audit` -> `cellalibi`

- Import name `ccc_audit` -> `cellalibi`; distribution name `ccc-audit` -> `cellalibi`; CLI `ccc-audit` -> `cellalibi`
  (`[project.scripts] cellalibi = "cellalibi.cli:main"`). Version 0.2.0 (v2 engine).
- Renamed in: `framework/src` (git mv), `framework/tests`, `framework/scripts/v2`, `framework/pyproject.toml`,
  `framework/README.md`, `framework/docs/design_notes.md`, and the not-yet-executed `framework/plans/v2_phase4.md`.
- Not renamed (historical records, left untouched): earlier plans (`phaseA*.md`, `phaseB.md`, `v2_phase1-2*.md`),
  `framework/reports/FRAMEWORK_REPORT_v1.md`, v1/v2 result folders, progress logs, and the unrelated thesis module
  `src/pdac/ccc_audit.py`.
- The tracked build artefact `framework/src/ccc_audit.egg-info/` is removed from git and `*.egg-info/` ignored.
- README and package docstring carry the note "previously named `ccc_audit` (v0.1, framework-v1)".
- Availability checks (logged in `framework/docs/NAMING_LOG.md`): https://pypi.org/pypi/cellalibi/json (404 = free)
  and the GitHub repository/user search API. Nothing is published.
- Done when: `python -m pytest framework/tests` passes with the new import name and the `cellalibi` CLI runs the
  U5-a crop config.

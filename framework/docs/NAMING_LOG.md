# Naming log: `ccc_audit` -> `cellalibi`

Checks run on 2026-11-18 04:28 UTC (2026-11-17 local). Raw responses are cached outside git in
`F:/Thesis/data/v2_cache/naming/`.

| Check | Query | Result |
|---|---|---|
| PyPI JSON API | https://pypi.org/pypi/cellalibi/json | HTTP 404, so the name is not registered |
| PyPI JSON API | https://pypi.org/pypi/cell-alibi/json | HTTP 404, so the name is not registered |
| GitHub repository search (unauthenticated REST) | https://api.github.com/search/repositories?q=cellalibi | `total_count` 0 |
| GitHub user/org search (unauthenticated REST) | https://api.github.com/search/users?q=cellalibi | `total_count` 0 |
| GitHub code search | `gh api search/code -f q=cellalibi` | **not run**: `gh` is not authenticated in this shell, and code search requires authentication. UNVERIFIED |

Nothing was published or registered: no PyPI upload and no GitHub repository creation.

Rename scope:
- Import name `ccc_audit` -> `cellalibi`.
- Distribution `ccc-audit` -> `cellalibi`.
- CLI `ccc-audit` -> `cellalibi`.
- Version 0.1.0 -> 0.2.0.

Historical files that keep the old name on purpose: framework-v1 plans and report, v1 result folders, and the
progress log. `src/pdac/ccc_audit.py` is an unrelated thesis module with the same stem and is unchanged.

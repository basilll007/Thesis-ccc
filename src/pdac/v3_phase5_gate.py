"""v3 Phase 5 gate: Visium gene programs on GSE277783 run only if a GEO SOFT / series-matrix file for GSE277783 exists in data/
with verifiable per-GSM anatomic-site labels. Plan (with pre-declared gene lists): docs/v3_plans/phase5_visium_programs_gse277783.md.
Writes a status file either way; never downloads anything.
"""
from __future__ import annotations

import gzip
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, REPO_ROOT, finish_run, save_json, start_run  # noqa: E402

PHASE = "v3_phase5"
PATTERNS = ["*GSE277783*soft*", "*GSE277783*series_matrix*", "*GSE277783*.soft", "*GSE277783*.soft.gz", "GSE277783_family*"]


def site_labels_in(path: Path) -> dict:
    """Count '!Sample_characteristics' / '!Sample_source_name' lines mentioning a site; only reads text, never executes."""
    opener = gzip.open if path.suffix == ".gz" else open
    n_char, n_site = 0, 0
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_characteristics") or line.startswith("!Sample_source_name"):
                    n_char += 1
                    if re.search(r"primary|metast|liver|lung|periton|pancrea|site|tissue", line, re.I):
                        n_site += 1
    except Exception as exc:  # pragma: no cover
        return {"error": str(exc)}
    return {"n_characteristic_lines": n_char, "n_site_like_lines": n_site}


def main() -> None:
    log = start_run("v3_phase5_gate")
    found = []
    for pat in PATTERNS:
        found += [p for p in DATA_DIR.rglob(pat) if p.is_file()]
    found = sorted(set(found))
    details = {str(p.relative_to(REPO_ROOT).as_posix()): site_labels_in(p) for p in found}
    verifiable = any(d.get("n_site_like_lines", 0) > 0 for d in details.values())
    eda = REPO_ROOT / "eda_gse277783" / "tables" / "sample_map.csv"
    status = {"checked_at": datetime.now().isoformat(timespec="seconds"), "patterns": PATTERNS, "soft_or_series_matrix_files": details,
              "raw_tar_on_disk": (DATA_DIR / "GSE277783_RAW.tar").exists(), "eda_sample_map_exists": eda.exists(),
              "status": "RUN" if (found and verifiable) else "BLOCKED",
              "reason": None if (found and verifiable) else
              ("no GEO SOFT/series-matrix file for GSE277783 in data/; site labels unverifiable (eda_gse277783/EDA_REPORT.md: organ site UNKNOWN for 49/57 slides)"
               if not found else "SOFT/series-matrix file present but no site-like characteristic lines found")}
    log.info("PHASE 5 GATE: %s", status)
    save_json(status, PHASE, "phase5_gate_status.json", "Phase 5 gate check (GEO SOFT/series matrix with site labels)", log)
    if status["status"] == "BLOCKED":
        log.warning("Phase 5 BLOCKED: %s", status["reason"])
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

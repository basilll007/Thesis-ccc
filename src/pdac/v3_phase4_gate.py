"""v3 Phase 4 gate: replication on GSE310352 runs only if LICENSE_CONFIRMED_GSE310352.txt exists at the repo root.
Plan: docs/v3_plans/phase4_replication_gse310352.md. Writes a status file either way; never downloads anything.
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, REPO_ROOT, finish_run, save_json, start_run  # noqa: E402

PHASE = "v3_phase4"


def main() -> None:
    log = start_run("v3_phase4_gate")
    gate = REPO_ROOT / "LICENSE_CONFIRMED_GSE310352.txt"
    raw = DATA_DIR / "GSE310352_RAW.tar"
    extracted = sorted((DATA_DIR / "pdac" / "geo_gse310352").glob("*.csv.gz")) if (DATA_DIR / "pdac" / "geo_gse310352").exists() else []
    status = {"checked_at": datetime.now().isoformat(timespec="seconds"), "gate_file": gate.name, "gate_file_exists": gate.exists(),
              "raw_tar_on_disk": raw.exists(), "raw_tar_bytes": raw.stat().st_size if raw.exists() else None,
              "n_extracted_csv_gz": len(extracted),
              "status": "RUN" if gate.exists() else "BLOCKED",
              "reason": None if gate.exists() else "licence unconfirmed: LICENSE_CONFIRMED_GSE310352.txt absent (user creates it after checking the data licence)",
              "also_outstanding": ["um/px for this CosMx run UNVERIFIED", "no author cell labels (marker + Leiden labels would be used)",
                                   "patient overlap with the Mendeley set UNVERIFIED (needs GEO sample metadata)"]}
    log.info("PHASE 4 GATE: %s", status)
    save_json(status, PHASE, "phase4_gate_status.json", "Phase 4 gate check (licence file) and outstanding blockers", log)
    if not gate.exists():
        log.warning("Phase 4 BLOCKED (licence unconfirmed); skipping replication run")
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

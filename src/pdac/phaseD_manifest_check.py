"""Phase D (v2): manifest completeness check.

Every file under results/pdac/phase{A,B,C,D}*, figures/pdac/phase{A,B,C,D}* and every v2 run log (logs/phase[A-D]*)
must have a row in results/pdac/MANIFEST.csv, and every manifest row must point to an existing file.
Rows missing because of concurrent read-modify-write of the manifest (two scripts registering at the same time)
are re-registered here with an explicit note. Writes results/pdac/phaseD/manifest_check.csv.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import MANIFEST, REPO_ROOT, finish_run, register, rel, save_csv, start_run  # noqa: E402

import pandas as pd  # noqa: E402

PHASE = "phaseD"
V2 = re.compile(r"^(results/pdac/phase[A-D]|figures/pdac/phase[A-D]|logs/phase[A-D])")
SCRIPT_GUESS = [("phaseB/ccc_sections", "src/pdac/phaseB3_ccc_rerun.py"),
                ("consensus", "src/pdac/phaseB2_consensus_labels.py"),
                ("phaseB/sections", "src/pdac/phaseB1_protein_labels.py")]


def main() -> None:
    log = start_run("phaseD_manifest_check")
    with MANIFEST.open(newline="", encoding="utf-8") as f:
        man = {r["path"] for r in csv.DictReader(f)}
    files = []
    for root in ("results/pdac", "figures/pdac", "logs"):
        for p in (REPO_ROOT / root).rglob("*"):
            if p.is_file() and V2.match(rel(p)):
                files.append(rel(p))
    rows = []
    for f in sorted(files):
        if f in man or f == rel(log.log_path):
            continue
        script = next((s for k, s in SCRIPT_GUESS if k in f), "")
        register(REPO_ROOT / f, f.split("/")[2] if f.startswith(("results", "figures")) else PHASE,
                 "re-registered by phaseD_manifest_check (row lost to concurrent manifest write or not registered)", script)
        rows.append({"path": f, "issue": "file not in manifest", "action": "re-registered"})
    for m in sorted(man):
        if V2.match(m) and not (REPO_ROOT / m).exists():
            rows.append({"path": m, "issue": "manifest row without file", "action": "reported only"})
    log.info("v2 files: %d; issues: %d", len(files), len(rows))
    for r in rows:
        log.info("%s", r)
    save_csv(pd.DataFrame(rows, columns=["path", "issue", "action"]), PHASE, "manifest_check.csv",
             "Manifest completeness check for v2 outputs (issues found and actions)", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

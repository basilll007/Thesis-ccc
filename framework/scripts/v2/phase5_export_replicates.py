"""v2 Phase 5 step 1 (thesis-model-v2 env): export the Phase 4 replicates used for the tool comparison.

Regenerates replicates with the identical code and seeds as phase4_benchmark.py and saves panel counts
(int32 npz) + per-section geometry. Verifies the export by re-auditing 3 replicates with cellalibi and comparing
with the Phase 4 chunk rows. Plan: framework/plans/v2_phase5.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import phase4_benchmark as pb  # noqa: E402

OUT = Path("F:/Thesis/data/v2_cache/phase5_replicates")
TOOL_SETTINGS = [st for st in pb.SETTINGS if st[0] != "S5"]


def n_reps_for(sec: str, sel: dict) -> int:
    return 20 if sec == sel["roles"]["smallest"] else 10


def main():
    sel = json.loads((pb.OUT / "section_selection.json").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for sec in sel["selected"]:
        sim = pb.load_section(sec)
        np.savez_compressed(OUT / f"{sec}__geometry.npz", xy=sim.xy, labels=sim.labels.astype(str),
                            genes=np.array(sim.genes))
        for scen, level in TOOL_SETTINGS:
            for rep in range(n_reps_for(sec, sel)):
                s = pb.seed_of(sec, scen, level, rep)
                f = OUT / f"{sec}__{scen}__{level}__{rep}.npz"
                if not f.exists():
                    counts = sim.draw(scen, level, s + 1)
                    np.savez_compressed(f, counts=counts.astype(np.int32))
                manifest.append({"section": sec, "scenario": scen, "level": level, "rep": rep, "seed": s, "file": f.name})
        print("exported", sec, flush=True)
    pd.DataFrame(manifest).to_csv(OUT / "manifest.csv", index=False)
    # verification: re-audit 3 exported replicates and compare with Phase 4 chunks
    from cellalibi.contact import AuditConfig, ExpressionData, SpatialContext
    from cellalibi.verdicts import audit_claim
    checks = []
    sec = sel["roles"]["smallest"]
    sim = pb.load_section(sec)
    ctx = SpatialContext(sim.xy, sim.labels, AuditConfig(n_perms=1000, seed=pb.seed_of(sec, "geom"), engine="0.2.0"))
    for scen, level, rep in (("S0", None, 0), ("S1", 2.0, 1), ("S4", 0.3, 2)):
        s = pb.seed_of(sec, scen, level, rep)
        counts = np.load(OUT / f"{sec}__{scen}__{level}__{rep}.npz")["counts"]
        au = audit_claim(ctx, ExpressionData(counts, sim.genes), cfg=AuditConfig(n_perms=1000, engine="0.2.0"), seed=s + 3, **pb.CLAIM)
        ch = pd.read_csv(pb.CHUNKS / f"{sec}__{scen}__{level}.csv")
        row = ch[ch.rep == rep].iloc[0]
        checks.append({"section": sec, "scenario": scen, "level": level, "rep": rep,
                       "phase4_verdict": row.verdict, "reaudit_verdict": au.primary.value,
                       "phase4_p_global": row.p_contact_global, "reaudit_p_global": au.metrics["p_contact_global"],
                       "match": bool(row.verdict == au.primary.value and np.isclose(row.p_contact_global, au.metrics["p_contact_global"]))})
    res = pd.DataFrame(checks)
    res.to_csv(pb.ROOT / "framework/results/v2_phase5/export_verification.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    (pb.ROOT / "framework/results/v2_phase5").mkdir(parents=True, exist_ok=True)
    main()

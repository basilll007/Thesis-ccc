"""v2 Phase 4: fixed section-selection rule (plan framework/plans/v2_phase4.md)."""
import json
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[3]
LAB = ROOT / "results/pdac/phaseB/sections"
rows = []
for f in sorted(LAB.glob("*_consensus_labels.csv.gz")):
    sec = f.name.replace("_consensus_labels.csv.gz", "")
    lab = pd.read_csv(f, usecols=["consensus_T"])["consensus_T"].astype(str)
    rows.append({"section": sec, "n_cells": len(lab), "n_myeloid": int((lab == "myeloid").sum()),
                 "n_CAF_fibroblast": int((lab == "CAF_fibroblast").sum())})
df = pd.DataFrame(rows).sort_values("n_cells").reset_index(drop=True)
df["eligible"] = (df.n_myeloid >= 100) & (df.n_CAF_fibroblast >= 100)
el = df[df.eligible].reset_index(drop=True)
k = len(el)
sel = {"smallest": el.section.iloc[0], "median_lower": el.section.iloc[(k - 1) // 2], "largest": el.section.iloc[k - 1]}
out = ROOT / "framework/results/v2_phase4"
out.mkdir(parents=True, exist_ok=True)
df.to_csv(out / "section_eligibility.csv", index=False)
(out / "section_selection.json").write_text(json.dumps({"rule": "eligible = >=100 myeloid and >=100 CAF_fibroblast "
    "(consensus_T, after QC); smallest / lower-median / largest n_cells among eligible", "n_eligible": k,
    "selected": [sel["smallest"], sel["median_lower"], sel["largest"]], "roles": sel}, indent=2))
print(df.to_string(index=False)); print(sel)

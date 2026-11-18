"""v2 Phase 5 analysis: per tool x setting, the fraction of replicates in which the planted claim
(CAF_fibroblast -> myeloid via FN1 -> ITGB1) is reported as significant / SUPPORTED, with Wilson 95 % CIs.
cellalibi v0.2.0 values come from the Phase 4 chunks for the same (section, setting, replicate); v0.2.1 values from
the GPU run on the exported replicates (development set)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from phase4_analyse import setting_label, wilson  # noqa: E402

ROOT = HERE.parents[2]
OUT = ROOT / "framework/results/v2_phase5"
P4 = ROOT / "framework/results/v2_phase4/chunks"


def main():
    tools = pd.concat([pd.read_csv(f) for f in sorted((OUT / "chunks").glob("*.csv"))], ignore_index=True)
    tools["level"] = tools["level"].astype(str).replace({"nan": "None"})
    used = tools[["section", "scenario", "level", "rep"]].drop_duplicates()
    p4 = pd.concat([pd.read_csv(f) for f in sorted(P4.glob("*.csv"))], ignore_index=True)
    p4["level"] = p4["level"].astype(str).replace({"nan": "None"})
    p4 = p4.merge(used, on=["section", "scenario", "level", "rep"])
    p4 = p4.assign(tool="cellalibi_v020", significant=p4.verdict == "SUPPORTED", status="ok")
    allr = pd.concat([tools, p4[["section", "scenario", "level", "rep", "tool", "significant", "status", "verdict"]]],
                     ignore_index=True)
    allr["setting"] = allr.apply(lambda r: setting_label({"scenario": r.scenario,
                                                          "level": None if r.level == "None" else r.level}), axis=1)
    allr.to_csv(OUT / "tool_calls_all.csv.gz", index=False)
    rows = []
    for (tool, st, sec), g in list(allr.groupby(["tool", "setting", "section"])) + \
            [((t, s, "pooled"), g) for (t, s), g in allr.groupby(["tool", "setting"])]:
        n = len(g)
        k = int(g.significant.sum())
        lo, hi = wilson(k, n)
        rows.append({"tool": tool, "setting": st, "section": sec, "n": n, "n_errors": int((g.status != "ok").sum()),
                     "n_nan_p": int(g["p"].isna().sum()) if ("p" in g and not tool.startswith("cellalibi")) else 0,
                     "reported_rate": k / n, "lo": lo, "hi": hi})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tool_detection_rates.csv", index=False)
    pooled = res[res.section == "pooled"].pivot(index="setting", columns="tool", values="reported_rate")
    order = ["S0 null", "S1 contact 1.25x", "S1 contact 1.5x", "S1 contact 2x", "S2 abundance 2x", "S3 zonation 50 um",
             "S3 zonation 200 um", "S3 zonation 400 um", "S4 spillover 0.1", "S4 spillover 0.2", "S4 spillover 0.3"]
    pooled = pooled.reindex(order)
    pooled.to_csv(OUT / "tool_detection_rates_pooled_wide.csv")
    ns = res[res.section == "pooled"].pivot(index="setting", columns="tool", values="n").reindex(order)
    print(pooled.round(3).to_string())
    print(ns.to_string())
    secw = res[(res.section != "pooled")].pivot_table(index=["section", "setting"], columns="tool", values="reported_rate")
    secw.to_csv(OUT / "tool_detection_rates_by_section_wide.csv")
    rt = tools.groupby(["tool", "section"]).seconds.median().unstack()
    rt.to_csv(OUT / "tool_runtime_median_seconds.csv")
    print(rt.round(1).to_string())


if __name__ == "__main__":
    main()

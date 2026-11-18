"""v2 Phase 6 analysis: per-tumour and cohort-level verdicts, artefact frequencies, runtime/memory.
Aggregation rules pre-declared in framework/plans/v2_phase6.md."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "framework/results/v2_phase6"
PRIORITY = ["NOT_TESTABLE", "EXPLAINED_BY_SPILLOVER", "EXPLAINED_BY_ABUNDANCE", "NOT_SUPPORTED", "EXPLAINED_BY_ZONATION",
            "SUPPORTED", "ERROR"]


def mode_conservative(vals):
    c = Counter(vals)
    top = max(c.values())
    tied = [v for v in PRIORITY if c.get(v, 0) == top]
    return tied[0]


def cohort_verdict(unit_verdicts):
    testable = [v for v in unit_verdicts if v not in ("NOT_TESTABLE", "ERROR")]
    if len(testable) < 3:
        return "NOT_TESTABLE", len(testable), sum(v == "SUPPORTED" for v in testable)
    ns = sum(v == "SUPPORTED" for v in testable)
    if ns > len(testable) / 2:
        return "SUPPORTED", len(testable), ns
    rest = [v for v in testable if v != "SUPPORTED"]
    return mode_conservative(rest), len(testable), ns


def main():
    df = pd.concat([pd.read_csv(f) for f in sorted((OUT / "chunks").glob("*.csv"))], ignore_index=True)
    df["claim"] = df.sender + "->" + df.receiver + " " + df.ligand.str.replace("ANY_CONTACT", "contact") + "-" + \
        df.receptor.str.replace("PROGRAM:", "program:")
    df.to_csv(OUT / "all_audits.csv.gz", index=False)
    # per unit (tumour / slide)
    per_unit = (df.groupby(["dataset", "kind", "claim", "unit"]).verdict
                .agg(lambda s: mode_conservative(list(s))).rename("unit_verdict").reset_index())
    nsec = df.groupby(["dataset", "kind", "claim", "unit"]).section.nunique().rename("n_sections").reset_index()
    per_unit = per_unit.merge(nsec)
    per_unit.to_csv(OUT / "verdicts_per_tumour.csv", index=False)
    rows = []
    for (ds, kind, claim), g in per_unit.groupby(["dataset", "kind", "claim"]):
        v, nt, ns = cohort_verdict(list(g.unit_verdict))
        cnt = Counter(g.unit_verdict)
        rows.append({"dataset": ds, "kind": kind, "claim": claim, "n_units": len(g), "n_testable_units": nt,
                     "n_supported_units": ns, "cohort_verdict": v,
                     **{f"units_{k}": cnt.get(k, 0) for k in PRIORITY[:-1]}})
    coh = pd.DataFrame(rows)
    coh.to_csv(OUT / "verdicts_per_claim_cohort.csv", index=False)
    # artefact frequencies over all claim x section audits
    freq = []
    for (ds, kind), g in df.groupby(["dataset", "kind"]):
        r = {"dataset": ds, "kind": kind, "n_audits": len(g)}
        for v in PRIORITY[:-1]:
            r[f"primary_{v}"] = float((g.verdict == v).mean())
        for c in ("power", "abundance", "zonation", "spillover"):
            r[f"check_{c}_FAIL"] = float((g[f"check_{c}"] == "FAIL").mean())
            r[f"check_{c}_NOT_RUN"] = float((g[f"check_{c}"] == "NOT_RUN").mean())
        tst = g[g.check_power == "PASS"]
        r["n_powered"] = len(tst)
        for c in ("abundance", "zonation", "spillover"):
            r[f"powered_check_{c}_FAIL"] = float((tst[f"check_{c}"] == "FAIL").mean()) if len(tst) else np.nan
        freq.append(r)
    pd.DataFrame(freq).to_csv(OUT / "artefact_frequencies.csv", index=False)
    # runtime / memory per section
    rt = df.groupby(["dataset", "section"]).agg(n_cells=("n_cells", "first"), n_claims=("claim", "size"),
                                                load_s=("section_load_seconds", "first"),
                                                total_s=("section_total_seconds", "first"),
                                                peak_rss_mb=("section_peak_rss_mb", "first"),
                                                median_claim_s=("seconds", "median"),
                                                tile_um=("tile_um_used", "first")).reset_index()
    rt.to_csv(OUT / "runtime_memory_per_section.csv", index=False)
    summ = {"n_audits": int(len(df)), "n_errors": int((df.verdict == "ERROR").sum()),
            "cohort_verdict_counts": {f"{ds}/{k}": dict(Counter(g.cohort_verdict)) for (ds, k), g in coh.groupby(["dataset", "kind"])},
            "supported_cohort_claims": coh[coh.cohort_verdict == "SUPPORTED"][["dataset", "claim", "n_testable_units", "n_supported_units"]].to_dict(orient="records")}
    (OUT / "phase6_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2)[:3000])
    print(pd.DataFrame(freq).round(3).T.to_string())
    print(rt.round(1).to_string(index=False))


if __name__ == "__main__":
    main()

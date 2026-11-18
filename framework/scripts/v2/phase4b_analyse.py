"""v2 Phase 4b analysis (engine v0.2.1, fresh seeds, two claims). Same definitions as phase4_analyse.py,
evaluated per claim (pooled over its 3 sections) and over both claims."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from phase4_analyse import ARTEFACT, VERDICTS, acceptable, correct, setting_label, wilson  # noqa: E402

ROOT = HERE.parents[2]
OUT = ROOT / "framework/results/v2_phase4b"


def metrics(g):
    n = len(g)
    ks = int((g.verdict == "SUPPORTED").sum())
    lo, hi = wilson(ks, n)
    pw = g[g.check_power == "PASS"]
    kp = int((pw.verdict == "SUPPORTED").sum())
    plo, phi = wilson(kp, len(pw))
    return {"n": n, "accuracy": float(g.correct.mean()), "acceptable_rate": float(g.acceptable.mean()),
            "supported_rate": ks / n, "supported_lo": lo, "supported_hi": hi,
            "not_testable_rate": float((g.verdict == "NOT_TESTABLE").mean()),
            "n_power_pass": len(pw), "supported_given_power_pass": kp / len(pw) if len(pw) else np.nan,
            "sgp_lo": plo, "sgp_hi": phi,
            **{f"fail_{c}": float((g[f"check_{c}"] == "FAIL").mean()) for c in ("power", "abundance", "zonation", "spillover")},
            **{f"notrun_{c}": float((g[f"check_{c}"] == "NOT_RUN").mean()) for c in ("abundance", "zonation", "spillover")},
            "median_audit_seconds": float(g.audit_seconds.median())}


def main():
    df = pd.concat([pd.read_csv(f) for f in sorted((OUT / "chunks").glob("*.csv"))], ignore_index=True)
    df["setting"] = df.apply(setting_label, axis=1)
    order = list(dict.fromkeys(df.sort_values(["scenario", "level"]).setting))
    df["correct"] = [correct(s, v) for s, v in zip(df.scenario, df.verdict)]
    df["acceptable"] = [acceptable(s, v) for s, v in zip(df.scenario, df.verdict)]
    df.to_csv(OUT / "benchmark_all_replicates.csv.gz", index=False)
    rows = []
    groups = [("both", "pooled", df)] + [(c, "pooled", g) for c, g in df.groupby("claim")] + \
             [(c, s, g) for (c, s), g in df.groupby(["claim", "section"])]
    for claim, sec, g in groups:
        pd.crosstab(g.setting, g.verdict).reindex(index=order, columns=VERDICTS, fill_value=0).to_csv(
            OUT / f"confusion_matrix_{claim}_{sec}.csv")
        for st in order:
            h = g[g.setting == st]
            if len(h):
                rows.append({"claim": claim, "section": sec, "setting": st, "scenario": h.scenario.iloc[0], **metrics(h)})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "metrics_by_setting.csv", index=False)
    crit = {}
    for claim in ["C1", "C2", "both"]:
        p = met[(met.claim == claim) & (met.section == "pooled")].set_index("setting")
        art = p[p.scenario.isin(ARTEFACT)]
        fooled = art[art.supported_rate > 0.10]
        s1 = p.loc["S1 contact 2x"]
        sec = met[(met.claim == claim) & (met.section != "pooled")]
        crit[claim] = {
            "criterion_1_false_supported_le_0.10": bool(len(fooled) == 0),
            "criterion_2_S1_2x_sensitivity_ge_0.80": bool(s1.supported_rate >= 0.80),
            "PASS": bool(len(fooled) == 0 and s1.supported_rate >= 0.80),
            "S1_2x_sensitivity": float(s1.supported_rate),
            "S1_2x_sensitivity_given_power_pass": float(s1.supported_given_power_pass),
            "S1_2x_fraction_testable": float(s1.n_power_pass / s1.n),
            "max_false_supported_rate": float(art.supported_rate.max()),
            "settings_that_fool_it": fooled.index.tolist(),
            "section_specific_fooled": [f"{a} / {b}" for a, b in zip(sec[sec.scenario.isin(ARTEFACT) & (sec.supported_rate > 0.10)].section,
                                                                    sec[sec.scenario.isin(ARTEFACT) & (sec.supported_rate > 0.10)].setting)],
            "S1_2x_by_section": dict(zip(sec[sec.setting == "S1 contact 2x"].section, sec[sec.setting == "S1 contact 2x"].supported_rate.round(3))),
        }
    (OUT / "pass_criterion.json").write_text(json.dumps(crit, indent=2))
    for claim in ["C1", "C2"]:
        print(f"== {claim}")
        print(pd.read_csv(OUT / f"confusion_matrix_{claim}_pooled.csv", index_col=0).to_string())
        p = met[(met.claim == claim) & (met.section == "pooled")].set_index("setting")
        print(p[["n", "accuracy", "supported_rate", "supported_hi", "not_testable_rate", "supported_given_power_pass",
                 "fail_power", "fail_abundance", "fail_zonation", "fail_spillover", "notrun_spillover"]].round(3).to_string())
    print(json.dumps(crit, indent=2))


if __name__ == "__main__":
    main()

"""v2 Phase 4 analysis: confusion matrix, accuracy, sensitivity, false-SUPPORTED rates, pass criterion.
Definitions are those pre-declared in framework/plans/v2_phase4.md."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "framework/results/v2_phase4"
VERDICTS = ["SUPPORTED", "EXPLAINED_BY_ABUNDANCE", "EXPLAINED_BY_ZONATION", "EXPLAINED_BY_SPILLOVER", "NOT_SUPPORTED",
            "NOT_TESTABLE"]
EXPECTED = {"S0": None, "S1": "SUPPORTED", "S2": "EXPLAINED_BY_ABUNDANCE", "S3": "EXPLAINED_BY_ZONATION",
            "S4": "EXPLAINED_BY_SPILLOVER", "S5": None}
ARTEFACT = ["S0", "S2", "S3", "S4"]


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def correct(scen, v):
    if scen == "S0":
        return v != "SUPPORTED"
    if scen == "S5":
        return v in ("NOT_TESTABLE", "SUPPORTED")
    return v == EXPECTED[scen]


def acceptable(scen, v):
    if scen == "S1":
        return v in ("SUPPORTED", "NOT_TESTABLE")
    return correct(scen, v)


def setting_label(r):
    lv = r["level"]
    if r["scenario"] == "S0":
        return "S0 null"
    if r["scenario"] == "S1":
        return f"S1 contact {float(lv):g}x"
    if r["scenario"] == "S2":
        return f"S2 abundance {float(lv):g}x"
    if r["scenario"] == "S3":
        return f"S3 zonation {float(lv):g} um"
    if r["scenario"] == "S4":
        return f"S4 spillover {float(lv):g}"
    return f"S5 n={int(float(lv))} receivers"


def main():
    files = sorted((OUT / "chunks").glob("*.csv"))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df["setting"] = df.apply(setting_label, axis=1)
    order = list(dict.fromkeys(df.sort_values(["scenario", "level"]).setting))
    df["correct"] = [correct(s, v) for s, v in zip(df.scenario, df.verdict)]
    df["acceptable"] = [acceptable(s, v) for s, v in zip(df.scenario, df.verdict)]
    df.to_csv(OUT / "benchmark_all_replicates.csv.gz", index=False)
    # confusion matrices
    cm = pd.crosstab(df.setting, df.verdict).reindex(index=order, columns=VERDICTS, fill_value=0)
    cm.to_csv(OUT / "confusion_matrix_pooled.csv")
    for sec, g in df.groupby("section"):
        pd.crosstab(g.setting, g.verdict).reindex(index=order, columns=VERDICTS, fill_value=0).to_csv(
            OUT / f"confusion_matrix_{sec}.csv")
    # per-setting metrics (pooled and per section)
    rows = []
    for (sec, st), g in [(("pooled", s), df[df.setting == s]) for s in order] + \
            [((sec, s), df[(df.section == sec) & (df.setting == s)]) for sec in sorted(df.section.unique()) for s in order]:
        n = len(g)
        if n == 0:
            continue
        scen = g.scenario.iloc[0]
        ks = int((g.verdict == "SUPPORTED").sum())
        kc = int(g.correct.sum())
        lo, hi = wilson(ks, n)
        clo, chi = wilson(kc, n)
        pw = g[g.check_power == "PASS"]
        rows.append({"section": sec, "setting": st, "scenario": scen, "n": n,
                     "accuracy": kc / n, "accuracy_lo": clo, "accuracy_hi": chi,
                     "acceptable_rate": float(g.acceptable.mean()),
                     "supported_rate": ks / n, "supported_lo": lo, "supported_hi": hi,
                     "not_testable_rate": float((g.verdict == "NOT_TESTABLE").mean()),
                     "supported_given_power_pass": float((pw.verdict == "SUPPORTED").mean()) if len(pw) else np.nan,
                     "n_power_pass": len(pw),
                     **{f"fail_{c}": float((g[f"check_{c}"] == "FAIL").mean()) for c in ("power", "abundance", "zonation", "spillover")},
                     **{f"notrun_{c}": float((g[f"check_{c}"] == "NOT_RUN").mean()) for c in ("zonation", "spillover")},
                     "median_audit_seconds": float(g.audit_seconds.median())})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "metrics_by_setting.csv", index=False)
    pooled = met[met.section == "pooled"].set_index("setting")
    # pass criterion
    art = pooled[pooled.scenario.isin(ARTEFACT)]
    fooled = art[art.supported_rate > 0.10]
    s1_2x = pooled.loc["S1 contact 2x"]
    crit = {
        "criterion_1_false_supported_le_0.10_every_artefact_setting": bool(len(fooled) == 0),
        "criterion_2_S1_2x_sensitivity_ge_0.80": bool(s1_2x.supported_rate >= 0.80),
        "S1_2x_sensitivity": float(s1_2x.supported_rate),
        "settings_that_fool_it": fooled.index.tolist(),
        "max_false_supported_rate": float(art.supported_rate.max()),
    }
    crit["BENCHMARK_PASS"] = crit["criterion_1_false_supported_le_0.10_every_artefact_setting"] and crit["criterion_2_S1_2x_sensitivity_ge_0.80"]
    sec_weak = met[(met.section != "pooled") & met.scenario.isin(ARTEFACT) & (met.supported_rate > 0.10)]
    crit["section_specific_fooled"] = [f"{a} / {b}" for a, b in zip(sec_weak.section, sec_weak.setting)]
    sec_s1 = met[(met.section != "pooled") & (met.setting == "S1 contact 2x")]
    crit["S1_2x_sensitivity_by_section"] = dict(zip(sec_s1.section, sec_s1.supported_rate.round(3)))
    (OUT / "pass_criterion.json").write_text(json.dumps(crit, indent=2))
    print(cm.to_string())
    print(pooled[["n", "accuracy", "supported_rate", "supported_lo", "supported_hi", "not_testable_rate",
                  "fail_power", "fail_abundance", "fail_zonation", "fail_spillover"]].round(3).to_string())
    print(json.dumps(crit, indent=2))


if __name__ == "__main__":
    main()

"""v2 Phase 1c/1d: pre-registered GSE57495 confirmatory hypoxia test (plan: framework/plans/v2_phase1.md).

Steps: GPL15048 probe map -> coverage gate (>= 10/12) -> overlap screen vs other cohorts ->
score -> Cox (stage-adjusted primary) -> one-sided Wald test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import CACHE, read_series_matrix, char_value  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/translational/phaseB/v2_correction"
OUT.mkdir(parents=True, exist_ok=True)
GENES = ["VEGFA", "SLC2A1", "CA9", "LDHA", "PGK1", "ADM", "NDRG1", "BNIP3", "ANKRD37", "P4HA1", "ENO1", "HK2"]
log: dict = {}

# ---------- GPL15048 annotation ----------
rows, hdr, on = [], None, False
with open(CACHE / "geo_matrix/GPL15048_annot_head.txt", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if line.startswith("!platform_table_begin"):
            on = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if on:
            parts = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = parts
            else:
                rows.append(parts + [""] * (len(hdr) - len(parts)))
gpl = pd.DataFrame(rows, columns=hdr)
log["gpl15048_rows"] = int(len(gpl))
log["gpl15048_rows_with_symbol"] = int((gpl["GeneSymbol"].str.strip() != "").sum())

meta, X = read_series_matrix(CACHE / "geo_matrix/GSE57495_series_matrix.txt.gz")
log["series_matrix_shape"] = list(X.shape)
log["matrix_max"] = float(np.nanmax(X.to_numpy()))
log["matrix_min"] = float(np.nanmin(X.to_numpy()))
if log["matrix_max"] > 100:
    X = np.log2(X + 1)
    log["log2_transformed"] = True
else:
    log["log2_transformed"] = False

pm = []
for g in GENES:
    conf = gpl[gpl["GeneSymbol"].str.strip() == g]["ID"].tolist()
    multi = gpl[gpl["GeneSymbol"].str.contains(rf"(^|\s){g}(\s|$)", regex=True)
                & gpl["GeneSymbol"].str.contains("///")]["ID"].tolist()
    present = [p for p in conf if p in X.index and X.loc[p].notna().all()]
    chosen = None
    if present:
        means = X.loc[present].mean(1)
        chosen = means.idxmax()
    pm.append({"gene": g, "n_confident_probes": len(conf), "n_present_in_matrix": len(present),
               "n_multigene_probes_excluded": len(multi), "chosen_probe_maxmean": chosen,
               "all_present_probes": ";".join(present)})
pm = pd.DataFrame(pm)
pm.to_csv(OUT / "gse57495_gpl15048_probe_map.csv", index=False)
cov = int(pm["chosen_probe_maxmean"].notna().sum())
log["coverage_genes"] = cov
log["coverage_pct"] = 100 * cov / len(GENES)
print(pm.to_string(index=False))
if cov < 10:
    log["status"] = "BLOCKED: coverage below 10/12"
    (OUT / "gse57495_run_log.json").write_text(json.dumps(log, indent=2))
    sys.exit("BLOCKED: coverage < 80%")

# ---------- clinical ----------
os_m = pd.to_numeric(char_value(meta, "overall survival (month)"), errors="coerce")
vs = char_value(meta, "vital.status").str.upper()
stage = char_value(meta, "Stage").str.upper().str.strip()
clin = pd.DataFrame({"gsm": meta["Sample_geo_accession"], "title": meta["Sample_title"],
                     "os_months": os_m.to_numpy(), "vital": vs.to_numpy(), "stage": stage.to_numpy()})
clin["event"] = clin["vital"].map({"DEAD": 1, "ALIVE": 0})
log["stage_values"] = clin["stage"].value_counts(dropna=False).to_dict()
log["vital_values"] = clin["vital"].value_counts(dropna=False).to_dict()


def stage_adv(s):
    if s is None or (isinstance(s, float) and np.isnan(s)) or s in ("", "NA", "NONE"):
        return np.nan
    s = str(s).replace("STAGE", "").strip()
    if s.startswith("1") or s.startswith("I") and not s.startswith("II") and not s.startswith("IV"):
        return 0.0
    if s in ("2A", "IIA"):
        return 0.0
    if s.startswith("2") or s.startswith("3") or s.startswith("4") or s.startswith("II") or s.startswith("IV"):
        return 1.0
    return np.nan


clin["stage_advanced"] = clin["stage"].map(stage_adv)
clin = clin[(clin["os_months"] > 0) & clin["event"].notna()].copy()
log["n_samples_valid_os"] = int(len(clin))
log["n_events"] = int(clin["event"].sum())

# ---------- 1d overlap screen ----------
ov = {"identifier_check": {}}
others = {}
for acc in ("GSE62452", "GSE28735"):
    m2, x2 = read_series_matrix(CACHE / f"geo_matrix/{acc}_series_matrix.txt.gz")
    ov["identifier_check"][acc] = {
        "shared_gsm": int(len(set(m2["Sample_geo_accession"]) & set(meta["Sample_geo_accession"]))),
        "shared_titles": int(len(set(m2["Sample_title"]) & set(meta["Sample_title"]))),
        "shared_cel_names": int(len({str(f).split("/")[-1] for f in m2["Sample_supplementary_file"]}
                                    & {str(f).split("/")[-1] for f in meta["Sample_supplementary_file"]})),
    }
    if acc == "GSE62452":
        # GPL6244 gene_assignment not cached; use the Phase B probe map only for a gene-level screen is too small.
        others[acc] = (m2, x2)

# gene-level GSE57495 matrix (MaxMean probe per symbol, confident probes only)
gsym = gpl[(gpl["GeneSymbol"].str.strip() != "") & ~gpl["GeneSymbol"].str.contains("///")][["ID", "GeneSymbol"]]
gsym = gsym[gsym["ID"].isin(X.index)]
Xg = X.loc[gsym["ID"]].copy()
Xg["sym"] = gsym["GeneSymbol"].str.strip().to_numpy()
Xg["mean"] = Xg.drop(columns="sym").mean(1)
Xg = Xg.sort_values("mean", ascending=False).drop_duplicates("sym").set_index("sym").drop(columns="mean")
Xg = Xg[clin["gsm"]]


def load_cbio(path):
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    df = df.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    df = df.drop(columns=[c for c in ("Entrez_Gene_Id",) if c in df.columns])
    return np.log2(df.apply(pd.to_numeric, errors="coerce") + 1)


screen = []
for name, path in (("TCGA_PAAD", "F:/Thesis/data/tcga/paad_pancan/lfs/data_mrna_seq_v2_rsem.txt"),
                   ("CPTAC_PDA", "F:/Thesis/data/tcga/paad_cptac/lfs/data_mrna_seq_v2_rsem.txt")):
    try:
        E = load_cbio(path)
    except Exception as e:  # pragma: no cover
        screen.append({"cohort": name, "status": f"not loaded: {e}"})
        continue
    genes = Xg.index.intersection(E.index)
    A = Xg.loc[genes].astype(float)
    Bm = E.loc[genes].dropna(axis=1, how="any")
    keep = (A.std(1) > 0) & (Bm.std(1) > 0)
    A, Bm = A[keep], Bm[keep]
    A = A.sub(A.mean(1), axis=0).div(A.std(1), axis=0)
    Bm = Bm.sub(Bm.mean(1), axis=0).div(Bm.std(1), axis=0)
    Ra = A.rank().to_numpy()
    Rb = Bm.rank().to_numpy()
    Ra = (Ra - Ra.mean(0)) / Ra.std(0)
    Rb = (Rb - Rb.mean(0)) / Rb.std(0)
    R = Ra.T @ Rb / Ra.shape[0]
    thr = max(0.8, float(np.quantile(R, 0.999)))
    flagged = int((R > thr).sum()) if np.quantile(R, 0.999) > 0.8 else int((R > 0.8).sum())
    screen.append({"cohort": name, "n_common_genes": int(len(A)), "n_gse57495": int(A.shape[1]),
                   "n_other": int(Bm.shape[1]), "max_spearman": float(R.max()),
                   "p99_9_spearman": float(np.quantile(R, 0.999)), "n_pairs_flagged_gt_0.8": int((R > 0.8).sum())})
ov["expression_screen"] = screen
ov["note_gse62452_gse28735"] = ("GPL6244 gene-level annotation not cached; cross-platform screen vs GSE62452/GSE28735 "
                               "not run. Identifier check and institution (Moffitt vs NCI) reported instead.")
(OUT / "gse57495_overlap_check.json").write_text(json.dumps(ov, indent=2))
print(json.dumps(ov, indent=2))

# ---------- scoring ----------
Z = []
for _, r in pm.dropna(subset=["chosen_probe_maxmean"]).iterrows():
    v = X.loc[r["chosen_probe_maxmean"], clin["gsm"]].astype(float)
    Z.append((v - v.mean()) / v.std(ddof=1))
raw = pd.concat(Z, axis=1).mean(1)
clin["hypoxia_z"] = ((raw - raw.mean()) / raw.std(ddof=1)).to_numpy()
# sensitivity: mean of all confident probes per gene
Zs = []
for _, r in pm.dropna(subset=["chosen_probe_maxmean"]).iterrows():
    probes = r["all_present_probes"].split(";")
    zz = []
    for p in probes:
        v = X.loc[p, clin["gsm"]].astype(float)
        zz.append((v - v.mean()) / v.std(ddof=1))
    Zs.append(pd.concat(zz, axis=1).mean(1))
raw_s = pd.concat(Zs, axis=1).mean(1)
clin["hypoxia_z_allprobes"] = ((raw_s - raw_s.mean()) / raw_s.std(ddof=1)).to_numpy()
clin.to_csv(OUT / "gse57495_cohort_scored.csv", index=False)

# ---------- Cox ----------
from lifelines import CoxPHFitter  # noqa: E402
from lifelines.statistics import proportional_hazard_test  # noqa: E402


def fit(df, cols, label):
    d = df[["os_months", "event"] + cols].dropna()
    cph = CoxPHFitter()
    cph.fit(d, duration_col="os_months", event_col="event")
    s = cph.summary.loc["hypoxia_z" if "hypoxia_z" in cols else cols[0]]
    b, se = float(s["coef"]), float(s["se(coef)"])
    rec = {"model": label, "n": int(len(d)), "events": int(d["event"].sum()), "covariates": "+".join(cols),
           "beta": b, "se": se, "hr": float(np.exp(b)), "hr_lo95": float(np.exp(b - 1.959964 * se)),
           "hr_hi95": float(np.exp(b + 1.959964 * se)), "p_two_sided": float(2 * norm.sf(abs(b / se))),
           "p_one_sided": float(norm.sf(b / se))}
    try:
        ph = proportional_hazard_test(cph, d, time_transform="rank")
        rec["ph_test_p_hypoxia"] = float(ph.summary.loc[cols[0], "p"])
    except Exception as e:  # pragma: no cover
        rec["ph_test_p_hypoxia"] = f"failed: {e}"
    return rec


n_stage_missing = int(clin["stage_advanced"].isna().sum())
ev_by_stage = clin.groupby("stage_advanced")["event"].sum().to_dict()
log["stage_missing"] = n_stage_missing
log["events_by_stage_advanced"] = {str(k): int(v) for k, v in ev_by_stage.items()}
use_adjusted = (n_stage_missing <= 0.2 * len(clin)) and len(ev_by_stage) == 2 and min(ev_by_stage.values()) >= 5
log["primary_model"] = "stage-adjusted" if use_adjusted else "unadjusted (stage rule failed)"

res = []
res.append(fit(clin, ["hypoxia_z", "stage_advanced"], "stage-adjusted" + (" [PRIMARY]" if use_adjusted else "")))
res.append(fit(clin, ["hypoxia_z"], "unadjusted" + ("" if use_adjusted else " [PRIMARY]")))
c2 = clin.rename(columns={"hypoxia_z": "_h", "hypoxia_z_allprobes": "hypoxia_z"})
r_s = fit(c2, ["hypoxia_z", "stage_advanced"], "stage-adjusted, all-probe sensitivity")
res.append(r_s)
res = pd.DataFrame(res)
res.to_csv(OUT / "gse57495_cox_results.csv", index=False)
print(res.to_string(index=False))

prim = res[res["model"].str.contains("PRIMARY")].iloc[0]
D = int(prim["events"])
power = {f"power_at_HR_{h}": float(norm.cdf(np.log(h) * np.sqrt(D) - 1.645)) for h in (1.32, 1.43, 1.48)}
log["power_with_observed_events"] = power
if prim["p_one_sided"] < 0.05 and prim["hr"] > 1:
    decision = "CONFIRMED in an independent cohort"
elif prim["hr"] > 1:
    decision = "NOT CONFIRMED (direction consistent, underpowered)"
else:
    decision = "NOT CONFIRMED (direction inconsistent)"
log["decision"] = decision

# secondary descriptive pooling with GSE62452 (post-hoc)
hm = pd.read_csv(ROOT / "results/translational/phaseB/hypoxia_models.csv")
r62 = hm[(hm["cohort"] == "GSE62452") & hm["model_type"].str.contains("stage-adjusted")].iloc[0]
b1, s1 = float(r62["beta"]), float(r62["se"])
log["gse62452_row_used"] = r62.to_dict()
b2, s2 = float(prim["beta"]), float(prim["se"])
w = np.array([1 / s1**2, 1 / s2**2])
bb = np.array([b1, b2])
bf = float((w * bb).sum() / w.sum())
sf = float(np.sqrt(1 / w.sum()))
Q = float((w * (bb - bf) ** 2).sum())
tau2 = max(0.0, (Q - 1) / (w.sum() - (w**2).sum() / w.sum()))
wr = 1 / (np.array([s1, s2]) ** 2 + tau2)
br = float((wr * bb).sum() / wr.sum())
sr = float(np.sqrt(1 / wr.sum()))
log["secondary_pooled_gse62452_gse57495_POSTHOC"] = {
    "fixed_hr": float(np.exp(bf)), "fixed_lo": float(np.exp(bf - 1.96 * sf)), "fixed_hi": float(np.exp(bf + 1.96 * sf)),
    "fixed_p_one_sided": float(norm.sf(bf / sf)), "Q": Q, "I2_pct": float(max(0, (Q - 1) / Q) * 100) if Q > 0 else 0.0,
    "tau2": float(tau2), "dl_hr": float(np.exp(br)), "dl_lo": float(np.exp(br - 1.96 * sr)), "dl_hi": float(np.exp(br + 1.96 * sr)),
    "dl_p_one_sided": float(norm.sf(br / sr)),
}
(OUT / "gse57495_run_log.json").write_text(json.dumps(log, indent=2, default=str))
print(json.dumps(log, indent=2, default=str))

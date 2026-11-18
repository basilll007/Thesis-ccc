"""v2 Phase 1a: map GSE28735 samples into GSE62452 (plan: framework/plans/v2_phase1.md)."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import CACHE, read_series_matrix, char_value  # noqa: E402

OUT = Path(__file__).resolve().parents[3] / "results/translational/phaseB/v2_correction"
OUT.mkdir(parents=True, exist_ok=True)

m62, x62 = read_series_matrix(CACHE / "geo_matrix/GSE62452_series_matrix.txt.gz")
m28, x28 = read_series_matrix(CACHE / "geo_matrix/GSE28735_series_matrix.txt.gz")
print("GSE62452", x62.shape, "GSE28735", x28.shape)

common = x62.index.intersection(x28.index)
a = x28.loc[common].dropna()
b = x62.loc[a.index].dropna()
a = a.loc[b.index]
print("common probes (no NA):", len(a))


def corr_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = (A - A.mean(0)) / A.std(0)
    B = (B - B.mean(0)) / B.std(0)
    return A.T @ B / A.shape[0]


r_p = corr_matrix(a.to_numpy(float), b.to_numpy(float))
ra = a.rank().to_numpy(float)
rb = b.rank().to_numpy(float)
r_s = corr_matrix(ra, rb)

gsm28 = list(a.columns)
gsm62 = list(b.columns)
best_62_for_28 = r_p.argmax(1)
best_28_for_62 = r_p.argmax(0)

# GSE62452 running sample number from CEL filename "GSMxxxx_<n>_..."
suppl = m62.set_index("Sample_geo_accession")["Sample_supplementary_file"]
num62 = {}
for g, f in suppl.items():
    mm = re.search(r"GSM\d+_(\d+)_", str(f))
    num62[g] = int(mm.group(1)) if mm else None

t28 = m28.set_index("Sample_geo_accession")
t62 = m62.set_index("Sample_geo_accession")
tis28 = char_value(m28, "tissue").to_numpy()
mon28 = pd.to_numeric(char_value(m28, "survival_month"), errors="coerce").to_numpy()
ev28 = pd.to_numeric(char_value(m28, "cancer_death"), errors="coerce").to_numpy()
tis62 = char_value(m62, "tissue").to_numpy()
mon62 = pd.to_numeric(char_value(m62, "survival months"), errors="coerce").to_numpy()
ev62 = pd.to_numeric(char_value(m62, "survival status"), errors="coerce").to_numpy()
idx28 = {g: i for i, g in enumerate(m28["Sample_geo_accession"])}
idx62 = {g: i for i, g in enumerate(m62["Sample_geo_accession"])}

rows = []
for k, g28 in enumerate(gsm28):
    j = best_62_for_28[k]
    g62 = gsm62[j]
    srt = np.sort(r_p[k])[::-1]
    i28, i62 = idx28[g28], idx62[g62]
    rows.append({
        "gse28735_gsm": g28,
        "gse28735_order": k + 1,
        "gse28735_title": t28.loc[g28, "Sample_title"],
        "gse28735_tissue": tis28[i28],
        "gse28735_os_months": mon28[i28],
        "gse28735_event": ev28[i28],
        "best_gse62452_gsm": g62,
        "gse62452_sample_number": num62.get(g62),
        "gse62452_title": t62.loc[g62, "Sample_title"],
        "gse62452_tissue": tis62[i62],
        "gse62452_os_months": mon62[i62],
        "gse62452_event": ev62[i62],
        "pearson_best": float(srt[0]),
        "pearson_second": float(srt[1]),
        "spearman_best_pair": float(r_s[k, j]),
        "mutual_best": bool(best_28_for_62[j] == k),
        "rule_E_identity": bool(srt[0] > 0.99 and best_28_for_62[j] == k),
        "rule_T_order_match": num62.get(g62) == k + 1,
        "rule_S_os_match": bool(np.isfinite(mon28[i28]) and np.isfinite(mon62[i62])
                                 and abs(round(mon62[i62]) - mon28[i28]) <= 1
                                 and ev28[i28] == ev62[i62]),
    })
df = pd.DataFrame(rows)
df.to_csv(OUT / "gse28735_in_gse62452_map.csv", index=False)

is_t = df["gse28735_tissue"].astype(str).str.upper().eq("T")
summary = {
    "n_common_probes": int(len(a)),
    "gse28735_samples": int(len(df)),
    "gse28735_tumours": int(is_t.sum()),
    "gse28735_tumours_with_os": int((is_t & df["gse28735_os_months"].gt(0)).sum()),
    "identity_matches_all": int(df["rule_E_identity"].sum()),
    "identity_matches_tumours": int((df["rule_E_identity"] & is_t).sum()),
    "identity_tumours_with_os": int((df["rule_E_identity"] & is_t & df["gse28735_os_months"].gt(0)).sum()),
    "order_matches_all": int(df["rule_T_order_match"].sum()),
    "os_matches_tumours": int((df["rule_S_os_match"] & is_t).sum()),
    "pearson_best_min": float(df["pearson_best"].min()),
    "pearson_best_median": float(df["pearson_best"].median()),
    "pearson_second_max": float(df["pearson_second"].max()),
    "spearman_best_min": float(df["spearman_best_pair"].min()),
    "matched_gse62452_tumour_tissue": int((df["rule_E_identity"] & is_t
                                           & df["gse62452_tissue"].str.contains("tumor", case=False)
                                           & ~df["gse62452_tissue"].str.contains("non-tumor", case=False)).sum()),
    "gse62452_tumours_total": int(pd.Series(tis62).str.contains("tumor", case=False).sum()
                                  - pd.Series(tis62).str.contains("non-tumor", case=False).sum()),
}
(OUT / "overlap_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

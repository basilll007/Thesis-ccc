"""v2 Phase 1a supplementary diagnostic (POST-HOC, added after the pre-declared rule E was seen).

Rule E (raw Pearson > 0.99) is dominated by the shared tissue profile, so it cannot separate
"same array re-normalised" from "different patient". Here each probe is centred within its own
series before correlating; distinct samples then correlate near 0 and the same array strongly.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import CACHE, read_series_matrix  # noqa: E402

OUT = Path(__file__).resolve().parents[3] / "results/translational/phaseB/v2_correction"
m62, x62 = read_series_matrix(CACHE / "geo_matrix/GSE62452_series_matrix.txt.gz")
m28, x28 = read_series_matrix(CACHE / "geo_matrix/GSE28735_series_matrix.txt.gz")
a = x28.dropna()
b = x62.loc[a.index].dropna()
a = a.loc[b.index]
A = a.to_numpy(float)
B = b.to_numpy(float)
A = A - A.mean(1, keepdims=True)
B = B - B.mean(1, keepdims=True)
A = (A - A.mean(0)) / A.std(0)
B = (B - B.mean(0)) / B.std(0)
r = A.T @ B / A.shape[0]

mp = pd.read_csv(OUT / "gse28735_in_gse62452_map.csv")
g62 = list(b.columns)
best = r.argmax(1)
srt = np.sort(r, 1)[:, ::-1]
mp["centred_r_best"] = srt[:, 0]
mp["centred_r_second"] = srt[:, 1]
mp["centred_best_gsm"] = [g62[j] for j in best]
mp["centred_best_agrees_raw"] = mp["centred_best_gsm"] == mp["best_gse62452_gsm"]
mp["centred_mutual_best"] = [bool(r[:, j].argmax() == k) for k, j in enumerate(best)]
# null: correlation of each GSE28735 sample with non-matched GSE62452 samples
mask = np.ones_like(r, bool)
mask[np.arange(len(best)), best] = False
null = r[mask]
mp.to_csv(OUT / "gse28735_in_gse62452_map_centred_posthoc.csv", index=False)

is_t = mp["gse28735_tissue"].astype(str).str.upper().eq("T")
summ = {
    "note": "POST-HOC diagnostic added after rule E results were seen",
    "centred_best_min": float(mp["centred_r_best"].min()),
    "centred_best_median": float(mp["centred_r_best"].median()),
    "centred_second_max": float(mp["centred_r_second"].max()),
    "centred_nonmatched_pairs_p99_9": float(np.quantile(null, 0.999)),
    "centred_nonmatched_pairs_max": float(null.max()),
    "centred_best_agrees_raw_all": int(mp["centred_best_agrees_raw"].sum()),
    "centred_mutual_best_all": int(mp["centred_mutual_best"].sum()),
    "tumours_centred_mutual_best_same_number": int((is_t & mp["centred_mutual_best"] & mp["rule_T_order_match"]).sum()),
    "tumours_with_os_centred_mutual_best": int((is_t & mp["centred_mutual_best"] & mp["gse28735_os_months"].gt(0)).sum()),
    "raw_mutual_best_all": int(mp["mutual_best"].sum()),
}
(OUT / "overlap_summary_centred_posthoc.json").write_text(json.dumps(summ, indent=2))
print(json.dumps(summ, indent=2))
bad = mp[is_t & ~mp["rule_S_os_match"]][["gse28735_gsm", "gse28735_order", "gse28735_os_months", "gse28735_event",
                                          "gse62452_sample_number", "gse62452_os_months", "gse62452_event"]]
print("tumours with OS mismatch:\n", bad.to_string(index=False))

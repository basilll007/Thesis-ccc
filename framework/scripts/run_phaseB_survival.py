import gzip
import io
import urllib.request
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import norm

PROBE_MAP = {
    "VEGFA": "8119898",
    "SLC2A1": "7915472",
    "CA9": "8155083",
    "LDHA": "7938777",
    "PGK1": "8168500",
    "ADM": "7938390",
    "NDRG1": "8153002",
    "BNIP3": "7937079",
    "ANKRD37": "8098604",
    "P4HA1": "7934278",
    "ENO1": "7912198",
    "HK2": "8042942",
}

def load_geo_expression_and_clinical(acc, tissue_filter, month_idx, status_idx, stage_idx=None):
    prefix = acc[:5] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    print(f"Downloading {acc} from {url}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        content = resp.read()
    gz = gzip.GzipFile(fileobj=io.BytesIO(content))
    
    chars = []
    sample_ids = []
    probe_data = {g: None for g in PROBE_MAP}
    target_probes = {p: g for g, p in PROBE_MAP.items()}
    
    table_started = False
    for line in gz:
        l = line.decode("utf-8", errors="ignore")
        if not table_started:
            if l.startswith("!Sample_geo_accession"):
                sample_ids = [p.strip('"\n ') for p in l.split("\t")[1:]]
            if l.startswith("!Sample_characteristics"):
                chars.append([p.strip('"\n ') for p in l.split("\t")[1:]])
            if l.startswith("!series_matrix_table_begin"):
                table_started = True
            continue
            
        if l.startswith("!series_matrix_table_end"):
            break
            
        parts = l.strip().split("\t")
        probe_id = parts[0].strip('"\n ')
        if probe_id in target_probes:
            gene = target_probes[probe_id]
            vals = [float(x) if x != "null" and x != "" else np.nan for x in parts[1:]]
            probe_data[gene] = np.array(vals)

    # Filter to tumours with valid OS
    tissue = chars[0]
    tumour_indices = [i for i, t in enumerate(tissue) if tissue_filter(t)]
    
    records = []
    for idx in tumour_indices:
        sid = sample_ids[idx]
        m_str = chars[month_idx][idx].split(":")[-1].strip()
        s_str = chars[status_idx][idx].split(":")[-1].strip()
        stage_val = chars[stage_idx][idx].split(":")[-1].strip() if stage_idx is not None else None
        try:
            m = float(m_str)
            s = int(float(s_str)) if s_str.replace('.', '', 1).isdigit() else (1 if s_str.upper() in ["DEAD", "1", "DECEASED"] else 0)
            if m > 0:
                rec = {
                    "sample_id": sid,
                    "os_months": m,
                    "event": s,
                    "stage": stage_val,
                }
                for g in PROBE_MAP:
                    rec[g] = probe_data[g][idx] if probe_data[g] is not None else np.nan
                records.append(rec)
        except Exception:
            pass
            
    df = pd.DataFrame(records)
    return df

df_gse62 = load_geo_expression_and_clinical(
    "GSE62452",
    lambda t: "tumor" in t.lower() and "non-tumor" not in t.lower(),
    month_idx=3,
    status_idx=4,
    stage_idx=2,
)
print("GSE62452 loaded valid patients:", len(df_gse62), "events:", df_gse62["event"].sum())

df_gse28 = load_geo_expression_and_clinical(
    "GSE28735",
    lambda t: "tissue: T" in t or t == "T",
    month_idx=1,
    status_idx=2,
)
print("GSE28735 loaded valid patients:", len(df_gse28), "events:", df_gse28["event"].sum())

import urllib.request
import gzip
import io
import pandas as pd
import numpy as np

def vet_geo_cohort(acc, title, pubmed_citation, tissue_filter, month_idx, status_idx, stage_idx=None, grade_idx=None):
    prefix = acc[:5] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        content = resp.read()
    gz = gzip.GzipFile(fileobj=io.BytesIO(content))
    chars = []
    sample_ids = []
    platform = "Unknown"
    for line in gz:
        l = line.decode("utf-8", errors="ignore")
        if l.startswith("!series_matrix_table_begin"):
            break
        if l.startswith("!Series_platform_id"):
            platform = l.split("\t")[1].strip('"\n ')
        if l.startswith("!Sample_geo_accession"):
            sample_ids = [p.strip('"\n ') for p in l.split("\t")[1:]]
        if l.startswith("!Sample_characteristics"):
            chars.append([p.strip('"\n ') for p in l.split("\t")[1:]])
            
    n_samples = len(sample_ids)
    # Filter
    tissue = chars[0]
    tumour_idx = [i for i, t in enumerate(tissue) if tissue_filter(t)]
    
    os_months = []
    events = []
    for i in tumour_idx:
        m_str = chars[month_idx][i].split(":")[-1].strip()
        s_str = chars[status_idx][i].split(":")[-1].strip()
        try:
            m = float(m_str)
            s = int(float(s_str)) if s_str.replace('.', '', 1).isdigit() else (1 if s_str.upper() in ["DEAD", "1", "DECEASED"] else 0)
            os_months.append(m)
            events.append(s)
        except Exception:
            pass
            
    return {
        "cohort_id": acc,
        "platform": platform,
        "citation": pubmed_citation,
        "n_total": n_samples,
        "n_tumours": len(tumour_idx),
        "n_with_os": len(os_months),
        "n_events": sum(events),
        "status": "OBTAINED",
        "access_terms": "Public Domain (NCBI GEO open access)",
        "source_url": url,
    }

print("Vetting GSE62452...")
res62 = vet_geo_cohort("GSE62452", "Yang et al. Cancer Cell 2016", "Yang et al. Cancer Cell 2016 (PMID: 26859458)",
                       lambda t: "tumor" in t.lower() and "non-tumor" not in t.lower(), 3, 4, stage_idx=2, grade_idx=1)
print(res62)

print("Vetting GSE28735...")
res28 = vet_geo_cohort("GSE28735", "Zhang et al. PLoS ONE 2012", "Zhang et al. PLoS ONE 2012 (PMID: 22363747)",
                       lambda t: "tissue: T" in t or t == "T", 1, 2)
print(res28)

print("Vetting GSE57495...")
res57 = vet_geo_cohort("GSE57495", "Chen et al. CCR 2015", "Chen et al. Clin Cancer Res 2015 (PMID: 26085527)",
                       lambda t: True, 0, 1, stage_idx=2)
print(res57)

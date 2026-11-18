import urllib.request
import gzip
import io

HYPOXIA_GENES = [
    "VEGFA", "SLC2A1", "CA9", "LDHA", "PGK1", "ADM",
    "NDRG1", "BNIP3", "ANKRD37", "P4HA1", "ENO1", "HK2"
]

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    content = resp.read()
gz = gzip.GzipFile(fileobj=io.BytesIO(content))

table_started = False
gene_map = {}
for line in gz:
    l = line.decode("utf-8", errors="ignore")
    if l.startswith("!platform_table_begin"):
        table_started = True
        header = next(gz).decode("utf-8", errors="ignore").strip().split("\t")
        id_col = header.index("ID")
        sym_col = header.index("Gene symbol")
        continue
    if table_started:
        if l.startswith("!platform_table_end"):
            break
        parts = l.strip().split("\t")
        if len(parts) > sym_col:
            sym = parts[sym_col].strip()
            # Multiple symbols can be separated by ///
            for s in sym.split("///"):
                s = s.strip()
                if s in HYPOXIA_GENES:
                    if s not in gene_map:
                        gene_map[s] = []
                    gene_map[s].append(parts[id_col])

print(f"Mapped {len(gene_map)} of {len(HYPOXIA_GENES)} genes:")
for g in HYPOXIA_GENES:
    print(f"  {g}: {gene_map.get(g)}")

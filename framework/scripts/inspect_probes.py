import urllib.request
import gzip
import io

def check_matrix_rows(acc):
    prefix = acc[:5] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        content = resp.read()
    gz = gzip.GzipFile(fileobj=io.BytesIO(content))
    table_started = False
    first_rows = []
    for line in gz:
        l = line.decode("utf-8", errors="ignore")
        if l.startswith("!series_matrix_table_begin"):
            table_started = True
            continue
        if table_started:
            first_rows.append(l.split("\t")[0].strip('"\n '))
            if len(first_rows) >= 10:
                break
    return first_rows

print("GSE62452 rows:", check_matrix_rows("GSE62452"))
print("GSE28735 rows:", check_matrix_rows("GSE28735"))
print("GSE57495 rows:", check_matrix_rows("GSE57495"))

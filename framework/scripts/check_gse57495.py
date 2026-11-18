import urllib.request
import gzip
import io

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57495/matrix/GSE57495_series_matrix.txt.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    content = resp.read()
gz = gzip.GzipFile(fileobj=io.BytesIO(content))
lines = []
for line in gz:
    l = line.decode("utf-8", errors="ignore")
    if l.startswith("!series_matrix_table_begin"):
        # read next 20 lines
        for _ in range(20):
            lines.append(next(gz).decode("utf-8", errors="ignore").strip().split("\t")[:3])
        break
for r in lines:
    print(r)

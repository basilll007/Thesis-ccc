import urllib.request
import gzip
import io

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    content = resp.read()
gz = gzip.GzipFile(fileobj=io.BytesIO(content))
for line in gz:
    l = line.decode("utf-8", errors="ignore")
    if l.startswith("!platform_table_begin"):
        # next line is header
        header = next(gz).decode("utf-8", errors="ignore")
        print("Table columns:\n", header.split("\t")[:10])
        first_row = next(gz).decode("utf-8", errors="ignore")
        print("First row:\n", first_row.split("\t")[:10])
        break

import urllib.request
import gzip
import io

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
print("Checking GPL6244 annotation URL:", url)
try:
    with urllib.request.urlopen(req, timeout=30) as resp:
        print("Status:", resp.status)
        content = resp.read()
        print("Downloaded size:", len(content))
        gz = gzip.GzipFile(fileobj=io.BytesIO(content))
        header_lines = []
        for line in gz:
            l = line.decode("utf-8", errors="ignore")
            if not l.startswith("#"):
                header_lines.append(l)
                if len(header_lines) >= 5:
                    break
        print("Columns / header:\n", "".join(header_lines))
except Exception as e:
    print("Error:", e)

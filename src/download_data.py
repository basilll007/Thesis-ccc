"""Download and unpack the selected fully public Xenium breast dataset."""
from __future__ import annotations

import hashlib
import shutil
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

from config import CONTEXT_PATH, RAW_DIR

URL = "https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip"
TARGET = RAW_DIR / "10x_xenium_breast"
ARCHIVE = TARGET / "Xenium_V1_human_Breast_2fov_outs.zip"


def download() -> Path:
    TARGET.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists():
        request = Request(URL, headers={"User-Agent": "thesis-spatial-ccc/0.1"})
        with urlopen(request) as response, ARCHIVE.open("wb") as output:
            shutil.copyfileobj(response, output)
    return ARCHIVE


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def unpack(path: Path) -> None:
    marker = TARGET / ".unpacked"
    if marker.exists():
        return
    with zipfile.ZipFile(path) as archive:
        archive.extractall(TARGET)
    matrix_archive = TARGET / "cell_feature_matrix.tar.gz"
    if matrix_archive.exists():
        with tarfile.open(matrix_archive, "r:gz") as archive:
            archive.extractall(TARGET)
    marker.write_text("unpacked\n", encoding="utf-8")


def update_context(path: Path) -> None:
    text = CONTEXT_PATH.read_text(encoding="utf-8")
    stamp = datetime.now(timezone.utc).date().isoformat()
    size = path.stat().st_size
    provenance = (
        "## Data Provenance\n\n"
        f"- **Dataset**: 10x Genomics Xenium V1 Human Breast, 2 FOV public output bundle.\n"
        f"- **Source URL/accession**: {URL} (no accession; official 10x public dataset).\n"
        f"- **Download date**: {stamp} (UTC).\n"
        f"- **Size**: {size:,} bytes compressed; SHA-256 `{sha256(path)}`.\n"
        "- **License/access terms**: Creative Commons Attribution 4.0 International (CC BY 4.0), as stated on the official 10x dataset/support page; no access application required.\n"
    )
    before = text.split("## Data Provenance", 1)[0]
    after = text.split("## Reproducibility", 1)[1]
    CONTEXT_PATH.write_text(before + provenance + "\n## Reproducibility" + after, encoding="utf-8")


if __name__ == "__main__":
    archive = download()
    unpack(archive)
    update_context(archive)
    print(f"Downloaded and unpacked: {archive}")

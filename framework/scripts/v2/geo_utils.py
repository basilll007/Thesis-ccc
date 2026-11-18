"""Small helpers to parse GEO series-matrix files (v2 Phase 1)."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("F:/Thesis/data/v2_cache")


def read_series_matrix(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (sample_meta, expression) from a GEO series matrix (.txt.gz).

    sample_meta: one row per sample, columns = !Sample_* keys (characteristics numbered).
    expression: probes x samples (float), index = ID_REF, columns = GSM accession.
    """
    meta: dict[str, list[str]] = {}
    char_n = 0
    rows = []
    header = None
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        in_table = False
        for line in fh:
            if not in_table:
                if line.startswith("!Sample_"):
                    parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
                    key = parts[0][1:]
                    if key == "Sample_characteristics_ch1":
                        key = f"char_{char_n}"
                        char_n += 1
                    elif key in meta:
                        k2 = 1
                        while f"{key}_{k2}" in meta:
                            k2 += 1
                        key = f"{key}_{k2}"
                    meta[key] = parts[1:]
                elif line.startswith("!series_matrix_table_begin"):
                    in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if header is None:
                header = parts
                continue
            rows.append(parts)
    n = len(meta["Sample_geo_accession"])
    meta_df = pd.DataFrame({k: (v + [""] * n)[:n] for k, v in meta.items()})
    expr = pd.DataFrame(rows, columns=header).set_index("ID_REF")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return meta_df, expr


def char_value(meta: pd.DataFrame, prefix: str) -> pd.Series:
    """Find the characteristic column whose values start with `prefix:` and return the value part."""
    out = pd.Series([None] * len(meta), index=meta.index, dtype=object)
    for c in [c for c in meta.columns if c.startswith("char_")]:
        vals = meta[c].astype(str)
        hit = vals.str.lower().str.startswith(prefix.lower() + ":")
        out[hit] = vals[hit].str.split(":", n=1).str[1].str.strip()
    return out

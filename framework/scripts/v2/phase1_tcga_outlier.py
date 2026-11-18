"""v2 Phase 1d follow-up (POST-HOC): inspect the top GSE57495-TCGA cross-platform correlations."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import CACHE, read_series_matrix
OUT = Path(__file__).resolve().parents[3] / "results/translational/phaseB/v2_correction"
gpl_rows, hdr, on = [], None, False
for line in open(CACHE / "geo_matrix/GPL15048_annot_head.txt", encoding="utf-8", errors="ignore"):
    if line.startswith("!platform_table_begin"): on = True; continue
    if line.startswith("!platform_table_end"): break
    if on:
        p = line.rstrip("\n").split("\t")
        if hdr is None: hdr = p
        else: gpl_rows.append(p + [""] * (len(hdr) - len(p)))
gpl = pd.DataFrame(gpl_rows, columns=hdr)
meta, X = read_series_matrix(CACHE / "geo_matrix/GSE57495_series_matrix.txt.gz")
g = gpl[(gpl.GeneSymbol.str.strip() != "") & ~gpl.GeneSymbol.str.contains("///")]
g = g[g.ID.isin(X.index)]
Xg = X.loc[g.ID].copy(); Xg["sym"] = g.GeneSymbol.str.strip().values
Xg["m"] = Xg.drop(columns="sym").mean(1)
Xg = Xg.sort_values("m", ascending=False).drop_duplicates("sym").set_index("sym").drop(columns="m")
E = pd.read_csv("F:/Thesis/data/tcga/paad_pancan/lfs/data_mrna_seq_v2_rsem.txt", sep="\t", low_memory=False)
E = E.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
E = np.log2(E.apply(pd.to_numeric, errors="coerce") + 1)
genes = Xg.index.intersection(E.index)
A = Xg.loc[genes].astype(float); B = E.loc[genes].dropna(axis=1)
k = (A.std(1) > 0) & (B.std(1) > 0); A, B = A[k], B[k]
A = A.sub(A.mean(1), axis=0).div(A.std(1), axis=0); B = B.sub(B.mean(1), axis=0).div(B.std(1), axis=0)
Ra = A.rank().to_numpy(); Rb = B.rank().to_numpy()
Ra = (Ra - Ra.mean(0)) / Ra.std(0); Rb = (Rb - Rb.mean(0)) / Rb.std(0)
R = Ra.T @ Rb / Ra.shape[0]
order = np.dstack(np.unravel_index(np.argsort(-R.ravel())[:10], R.shape))[0]
rows = []
for i, j in order:
    rows.append({"gse57495_gsm": A.columns[i], "gse57495_title": meta.set_index("Sample_geo_accession").loc[A.columns[i], "Sample_title"],
                 "tcga_sample": B.columns[j], "tcga_tss_code": B.columns[j].split("-")[1] if "-" in B.columns[j] else "",
                 "spearman": float(R[i, j]),
                 "second_best_for_gse_sample": float(np.sort(R[i])[::-1][1]),
                 "second_best_for_tcga_sample": float(np.sort(R[:, j])[::-1][1])})
df = pd.DataFrame(rows); df.to_csv(OUT / "gse57495_vs_tcga_top_pairs_posthoc.csv", index=False)
print(df.to_string(index=False))
row_max = R.max(1); print("per-GSE57495-sample best r: median", np.median(row_max), "n>0.6", int((row_max > 0.6).sum()))
print("TSS codes of TCGA samples:", pd.Series([c.split('-')[1] for c in B.columns]).value_counts().head(20).to_dict())

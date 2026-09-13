"""Task A: Detection-rate sweep across all 20 curated LR pairs.

This script evaluates single-cell and per-cell-type detection rates for all 37 unique genes
across the 20 curated ligand-receptor pairs in src/pipeline.py, identifying which targets
are absent from the 10x Xenium 280-gene panel and determining the rate-limiting gene for each pair.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

# Add src to path for direct imports
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from pipeline import LR_PAIRS


def main() -> None:
    t0 = time.time()
    h5ad_path = Path("data/processed/xenium_breast_baseline.h5ad")
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Processed AnnData missing: {h5ad_path}")

    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    n_total_cells = adata.n_obs
    var_names_set = set(adata.var_names)

    # 1. Identify roles and unique genes
    ligands = set([lig for lig, _ in LR_PAIRS])
    receptors = set([rec for _, rec in LR_PAIRS])
    all_genes = sorted(list(ligands | receptors))
    print(f"Total unique genes across 20 LR pairs: {len(all_genes)}")

    roles = {}
    for g in all_genes:
        if g in ligands and g in receptors:
            roles[g] = "both"
        elif g in ligands:
            roles[g] = "ligand"
        else:
            roles[g] = "receptor"

    # Cell type groups
    cell_types = ["epithelial", "immune", "endothelial", "fibroblast"]
    ct_masks = {ct: (adata.obs["cell_type"] == ct).values for ct in cell_types}
    ct_counts = {ct: int(np.sum(ct_masks[ct])) for ct in cell_types}

    # 2. Gene-level detection diagnostics
    gene_records = []
    gene_det_rates = {}

    for gene in all_genes:
        in_panel = gene in var_names_set
        if in_panel:
            idx = adata.var_names.get_loc(gene)
            col = adata.X[:, idx]
            expr = np.asarray(col.toarray() if sp.issparse(adata.X) else col).ravel()

            n_det_global = int(np.sum(expr > 0))
            rate_global = n_det_global / n_total_cells
            mean_global = float(np.mean(expr))
            gene_det_rates[gene] = rate_global

            ct_rates = {}
            for ct in cell_types:
                m = ct_masks[ct]
                ct_expr = expr[m]
                ct_det = int(np.sum(ct_expr > 0))
                ct_rates[f"{ct}_detection_rate"] = ct_det / len(ct_expr) if len(ct_expr) > 0 else 0.0
        else:
            n_det_global = 0
            rate_global = 0.0
            mean_global = 0.0
            gene_det_rates[gene] = 0.0
            ct_rates = {f"{ct}_detection_rate": 0.0 for ct in cell_types}

        rec = {
            "gene": gene,
            "role": roles[gene],
            "in_panel": in_panel,
            "global_n_cells": n_total_cells,
            "global_n_detected": n_det_global,
            "global_detection_rate": rate_global,
            "global_mean_expression": mean_global,
        }
        rec.update(ct_rates)
        gene_records.append(rec)

    df_gene_diag = pd.DataFrame(gene_records)
    out_gene_csv = Path("results/gene_detection_diagnostic.csv")
    out_gene_csv.parent.mkdir(parents=True, exist_ok=True)
    df_gene_diag.to_csv(out_gene_csv, index=False)
    print(f"Saved {len(df_gene_diag)} gene diagnostics to {out_gene_csv}")

    # 3. LR pair coverage summary
    # Check current null_model_comparison.csv
    comp_csv = Path("results/null_model_comparison.csv")
    ranked_pairs = set()
    if comp_csv.exists():
        df_comp = pd.read_csv(comp_csv)
        for _, r in df_comp.iterrows():
            ranked_pairs.add((r["ligand"], r["receptor"]))

    pair_records = []
    for lig, rec in LR_PAIRS:
        lig_in = lig in var_names_set
        rec_in = rec in var_names_set
        in_ranked = (lig, rec) in ranked_pairs

        lig_rate = gene_det_rates[lig]
        rec_rate = gene_det_rates[rec]

        if not lig_in and not rec_in:
            limiting = "both (neither in panel)"
            lim_rate = 0.0
            reason = "Both ligand and receptor omitted from 10x Xenium 280-gene panel"
        elif not lig_in:
            limiting = f"{lig} (ligand not in panel)"
            lim_rate = 0.0
            reason = f"Ligand {lig} omitted from 10x Xenium 280-gene panel"
        elif not rec_in:
            limiting = f"{rec} (receptor not in panel)"
            lim_rate = 0.0
            reason = f"Receptor {rec} omitted from 10x Xenium 280-gene panel"
        else:
            if lig_rate <= rec_rate:
                limiting = f"{lig} (ligand)"
                lim_rate = lig_rate
            else:
                limiting = f"{rec} (receptor)"
                lim_rate = rec_rate

            if in_ranked:
                reason = "Active in ranked baseline"
            else:
                reason = "Both genes present in panel but zero co-expression across sender-receiver pairs"

        pair_records.append({
            "ligand": lig,
            "receptor": rec,
            "in_ranked_output": in_ranked,
            "ligand_in_panel": lig_in,
            "receptor_in_panel": rec_in,
            "ligand_detection_rate": lig_rate,
            "receptor_detection_rate": rec_rate,
            "limiting_gene": limiting,
            "limiting_rate": lim_rate,
            "reason": reason,
        })

    df_pair_summary = pd.DataFrame(pair_records)
    out_pair_csv = Path("results/lr_pair_coverage_summary.csv")
    df_pair_summary.to_csv(out_pair_csv, index=False)
    print(f"Saved {len(df_pair_summary)} pair coverage summaries to {out_pair_csv}")

    runtime = time.time() - t0
    print(f"\nTask A completed in {runtime:.2f} seconds.")
    print("\nSummary of 20 LR pairs:")
    print(df_pair_summary[["ligand", "receptor", "in_ranked_output", "ligand_in_panel", "receptor_in_panel", "limiting_gene", "limiting_rate"]].to_string())


if __name__ == "__main__":
    main()

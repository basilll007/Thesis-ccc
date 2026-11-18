"""Command-line interface for the ccc_audit framework.

What this module does:
Provides a unified entry point: `ccc-audit run config.yaml` that loads spatial data,
constructs contact graphs, evaluates CCC scores (C, E, M, Q) under unconstrained and block nulls,
audits segmentation spillover, assigns evidence verdicts, and exports summary tables,
diagnostic figures, and an audit markdown report.

Why:
Enables reproducible, end-to-end execution of spatial communication audits from a single
configuration file with zero manual intervention.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from ccc_audit.io import load_spatial_dataset, validate_spatial_data
from ccc_audit.graphs import knn_graph, assign_tiles
from ccc_audit.scores import compute_ccc_scores
from ccc_audit.nulls import global_null, block_null
from ccc_audit.spillover import audit_spillover
from ccc_audit.verdicts import decide_verdict, generate_claims_register, VerdictResult


def run_audit(config_path: str | Path, output_dir_override: str | None = None) -> int:
    """Execute end-to-end audit workflow from YAML configuration."""
    t0 = time.time()
    cfg_file = Path(config_path)
    if not cfg_file.exists():
        print(f"Error: configuration file '{cfg_file}' not found.", file=sys.stderr)
        return 1

    with open(cfg_file, "r") as f:
        cfg: dict[str, Any] = yaml.safe_load(f)

    out_dir = Path(output_dir_override) if output_dir_override else Path(cfg.get("output_dir", "results/audit"))
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Starting CCC Audit: {cfg.get('title', 'Spatial Communication Audit')} ===")
    print(f"Config: {cfg_file}")
    print(f"Output directory: {out_dir}")

    # 1. Load data
    data_cfg = cfg["data"]
    adata = load_spatial_dataset(data_cfg["path"], obs_mapping=data_cfg.get("obs_mapping"))
    meta = validate_spatial_data(adata)
    print(f"Loaded dataset: {meta['n_cells']} cells, {meta['n_genes']} genes across {meta['n_sections']} sections.")

    # 2. Build spatial graph
    xy = np.column_stack([adata.obs["x_um"].to_numpy(), adata.obs["y_um"].to_numpy()])
    k_nn = int(cfg.get("graph", {}).get("k", 6))
    rows, cols = knn_graph(xy, k=k_nn)
    print(f"Constructed {k_nn}-NN graph with {len(rows)} directed edges.")

    # 3. Interactions to audit
    interactions = [
        (item["sender"], item["receiver"], item["ligand"], item["receptor"])
        for item in cfg["interactions"]
    ]
    lr_pairs = sorted({(l, r) for _, _, l, r in interactions})
    print(f"Auditing {len(interactions)} interactions across {len(lr_pairs)} LR pairs.")

    # 4. Null permutations
    n_perms = int(cfg.get("nulls", {}).get("n_perms", 500))
    tile_size = float(cfg.get("nulls", {}).get("tile_size_um", 100.0))
    seed = int(cfg.get("seed", 42))

    labels_str = adata.obs["cell_type"].astype(str).to_numpy()
    types = sorted(set(labels_str))
    type_map = {t: i for i, t in enumerate(types)}
    lab_encoded = np.array([type_map[x] for x in labels_str], dtype=np.int64)

    print(f"Generating {n_perms} global and {tile_size}um block permutations...")
    g_perms = global_null(lab_encoded, n_perms=n_perms, seed=seed)
    b_perms, tile_info = block_null(lab_encoded, xy, tile_size_um=tile_size, n_perms=n_perms, seed=seed)

    # 5. Compute scores under both nulls
    print("Computing scores under global null...")
    res_global = compute_ccc_scores(adata, lr_pairs, interactions, rows, cols, g_perms, alpha=0.05)
    print("Computing scores under block null...")
    res_block = compute_ccc_scores(adata, lr_pairs, interactions, rows, cols, b_perms, alpha=0.05)

    df_g = res_global.dataframe
    df_b = res_block.dataframe

    # Save detailed scores
    df_g.to_csv(out_dir / "scores_global_null.csv", index=False)
    df_b.to_csv(out_dir / "scores_block_null.csv", index=False)

    # 6. Optional spillover audit
    spillover_results: dict[str, Any] = {}
    if cfg.get("spillover", {}).get("enabled", True):
        print("Running segmentation spillover audit...")
        # Build expression DataFrame for receiver/sender types
        target_genes = list({g for _, _, l, r in interactions for g in (l, r)})
        var_names = list(adata.var_names)
        t_indices = [var_names.index(g) for g in target_genes if g in var_names]
        X_sub = adata.X[:, t_indices]
        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()
        exp_df = pd.DataFrame(X_sub, columns=[var_names[i] for i in t_indices])

        # Pick primary sender and receiver from first interaction
        s_type, r_type = interactions[0][0], interactions[0][1]
        s_mask = labels_str == s_type
        r_mask = labels_str == r_type

        diag = audit_spillover(xy, exp_df, s_mask, r_mask, target_genes)
        spillover_results = {
            "median_cell_diameter_um": diag.median_cell_diameter_um,
            "mixing_excess_pct": diag.mixing_excess_pct,
            "is_explained_by_spillover": diag.is_explained_by_spillover,
            "n_contamination_markers": len(diag.contamination_markers),
            "shared_genes_removed": diag.shared_genes_to_remove,
        }
        pd.DataFrame([spillover_results]).to_csv(out_dir / "spillover_summary.csv", index=False)

    # 7. Evaluate verdicts
    verdicts: list[VerdictResult] = []
    for i, row_g in df_g.iterrows():
        row_b = df_b.iloc[i]
        target_label = f"{row_g['sender']}->{row_g['receiver']} ({row_g['ligand']}-{row_g['receptor']})"
        rec_count = int((labels_str == row_g['receiver']).sum())

        v = decide_verdict(
            target_id=target_label,
            n_receiver_cells=rec_count,
            n_qualifying_edges=int(row_g["n_qualifying_edges"]),
            comp_p=float(row_g["comp_p_plus1"]),
            matched_p=float(row_g["matched_p_plus1"]),
            ratio_p_global=float(row_g["ratio_p_plus1"]),
            ratio_p_block=float(row_b["ratio_p_plus1"]),
            is_spillover=spillover_results.get("is_explained_by_spillover"),
            spillover_excess_pct=spillover_results.get("mixing_excess_pct"),
        )
        verdicts.append(v)

    claims_df = generate_claims_register(verdicts)
    claims_df.to_csv(out_dir / "claims_register.csv", index=False)

    # 8. Generate summary figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Panel A: Ratio Q vs Matched baseline
    ax = axes[0]
    ax.scatter(df_g["matched_obs"], df_g["ratio_obs"], c="tab:blue", edgecolors="k", s=60, alpha=0.8)
    ax.axhline(1.0, color="gray", linestyle="--", label="Matched baseline (Q=1.0)")
    ax.set_xlabel("Matched Dyadic Baseline M (Lbar * Rbar)")
    ax.set_ylabel("Proximity Ratio Q (E / M)")
    ax.set_title("Edge Enrichment Q vs Non-Spatial Baseline M")
    ax.legend()

    # Panel B: Verdict distribution
    ax2 = axes[1]
    v_counts = claims_df["verdict"].value_counts()
    v_counts.plot(kind="bar", ax=ax2, color="tab:orange", edgecolor="black")
    ax2.set_title("Audit Verdict Distribution")
    ax2.set_ylabel("Count")
    ax2.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    fig_path = figures_dir / "audit_summary_figure.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()

    # 9. Write Markdown Report
    elapsed = time.time() - t0
    report_lines = [
        f"# Spatial Communication Audit Report: {cfg.get('title', 'CCC Audit')}",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Execution Runtime:** {elapsed:.2f} seconds",
        f"**Dataset:** `{data_cfg['path']}` ({meta['n_cells']} cells, {meta['n_genes']} genes)",
        "",
        "## Executive Summary & Verdicts",
        "",
        "| ID | Target Interaction | Verdict | Primary Reason |",
        "|---|---|---|---|",
    ]
    for _, r in claims_df.iterrows():
        report_lines.append(f"| {r['id']} | `{r['target']}` | **{r['verdict']}** | {r['primary_reason']} |")

    report_lines.extend([
        "",
        "## Diagnostic Details",
        f"- **Graph connectivity:** {k_nn}-NN with {len(rows)} directed edges.",
        f"- **Null model calibrations:** {n_perms} permutations; block tile size = {tile_size} um.",
        f"- **Segmentation spillover excess:** {spillover_results.get('mixing_excess_pct', 0.0):.1f}% (explained by spillover: {spillover_results.get('is_explained_by_spillover', False)}).",
        "",
        "## Generated Deliverables",
        f"- Detailed global scores: `{out_dir / 'scores_global_null.csv'}`",
        f"- Detailed block scores: `{out_dir / 'scores_block_null.csv'}`",
        f"- Claims register table: `{out_dir / 'claims_register.csv'}`",
        f"- Diagnostic figure: `{fig_path}`",
    ])

    report_path = out_dir / "AUDIT_REPORT.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"\nAudit completed successfully in {elapsed:.1f}s.")
    print(f"Report written to: {report_path}")
    return 0


def main() -> None:
    """CLI entrypoint."""
    parser = argparse.ArgumentParser(
        prog="ccc-audit",
        description="Run spatial cell-cell communication audit and confounding diagnosis.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run audit pipeline from config file.")
    run_parser.add_argument("config", type=str, help="Path to config.yaml")
    run_parser.add_argument("--output-dir", type=str, default=None, help="Directory to save audit deliverables")

    args = parser.parse_args()
    if args.command == "run":
        sys.exit(run_audit(args.config, args.output_dir))


if __name__ == "__main__":
    main()

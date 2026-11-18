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
from ccc_audit.contact import AuditConfig, ExpressionData, SpatialContext
from ccc_audit.verdicts import PRIORITY_RULE, CHECKS, ClaimAudit, audit_claim, claims_table


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

    # 4. Thesis-style C/E/M/Q tables under label-permutation nulls (kept for continuity with v1/thesis)
    n_perms = int(cfg.get("nulls", {}).get("n_perms", 1000))
    tile_size = float(cfg.get("nulls", {}).get("tile_size_um", 100.0))
    seed = int(cfg.get("seed", 42))

    labels_str = adata.obs["cell_type"].astype(str).to_numpy()
    types = sorted(set(labels_str))
    type_map = {t: i for i, t in enumerate(types)}
    lab_encoded = np.array([type_map[x] for x in labels_str], dtype=np.int64)

    print(f"Generating {n_perms} global and {tile_size}um block label permutations (thesis C/E/M/Q tables)...")
    g_perms = global_null(lab_encoded, n_perms=n_perms, seed=seed)
    b_perms, tile_info = block_null(lab_encoded, xy, tile_size_um=tile_size, n_perms=n_perms, seed=seed)
    res_global = compute_ccc_scores(adata, lr_pairs, interactions, rows, cols, g_perms, alpha=0.05)
    del g_perms
    res_block = compute_ccc_scores(adata, lr_pairs, interactions, rows, cols, b_perms, alpha=0.05)
    del b_perms
    df_g = res_global.dataframe
    df_b = res_block.dataframe
    df_g.to_csv(out_dir / "scores_global_null.csv", index=False)
    df_b.to_csv(out_dir / "scores_block_null.csv", index=False)

    # 5-7. v2 claim audit: every check computed per claim, spillover computed for each claim's own pair
    acfg = AuditConfig(k=k_nn, n_perms=n_perms, tile_um=tile_size, seed=seed,
                       **{k: v for k, v in cfg.get("audit", {}).items()})
    ctx = SpatialContext(xy, labels_str, acfg)
    if "counts" in adata.layers:
        ex = ExpressionData(adata.layers["counts"], adata.var_names)
    else:
        ex = ExpressionData(adata.X, adata.var_names, is_log1p_cp10k=True)
    print(f"Running v2 claim audit ({n_perms} permutations; minimum attainable p = {acfg.min_p:.6f})...")
    audits: list[ClaimAudit] = []
    for (s_t, r_t, l_g, r_g) in interactions:
        audits.append(audit_claim(ctx, ex, s_t, r_t, l_g, r_g, cfg=acfg))
    claims_df = claims_table(audits)
    claims_df = claims_df.rename(columns={"primary_verdict": "verdict"})
    claims_df.to_csv(out_dir / "claims_register.csv", index=False)
    pd.DataFrame([au.metrics for au in audits]).to_csv(out_dir / "claim_metrics.csv", index=False)

    # 8. Generate summary figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Panel A: Ratio Q vs Matched baseline
    ax = axes[0]
    ax.scatter(claims_df["p_M_label"].clip(lower=1e-4), claims_df["Q"], c="tab:blue", edgecolors="k", s=60, alpha=0.8)
    ax.axhline(1.0, color="gray", linestyle="--", label="Matched baseline (Q=1.0)")
    ax.set_xscale("log")
    ax.set_xlabel("naive M label-null p")
    ax.set_ylabel("Proximity Ratio Q (E / M)")
    ax.set_title("Contact ratio Q vs naive-score p")
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
        "Every check is computed for every claim (PASS / FAIL / NOT_RUN). The primary verdict follows the",
        f"pre-declared priority rule: {PRIORITY_RULE}.",
        "",
        "| ID | Claim | power | abundance | zonation | spillover | Primary verdict | Failing checks |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, r in claims_df.iterrows():
        report_lines.append(f"| {r['id']} | `{r['target']}` | {r['check_power']} | {r['check_abundance']} | "
                            f"{r['check_zonation']} | {r['check_spillover']} | **{r['verdict']}** | {r['failing_checks'] or '-'} |")
    report_lines += ["", "### Check reasons", ""]
    for _, r in claims_df.iterrows():
        report_lines.append(f"- `{r['target']}`: " + "; ".join(f"{c}: {r['reason_' + c]}" for c in CHECKS))

    report_lines.extend([
        "",
        "## Diagnostic Details",
        f"- **Graph connectivity:** {k_nn}-NN with {len(rows)} directed edges.",
        f"- **Permutations:** {n_perms} per null; minimum attainable p = 1/{n_perms + 1} = {1 / (n_perms + 1):.6f}; block tile = {tile_size} um.",
        f"- **Spillover:** computed separately for each claim's sender/receiver pair on all {adata.n_vars} genes (see claims_register.csv).",
        "",
        "## Generated Deliverables",
        f"- Detailed global scores: `{out_dir / 'scores_global_null.csv'}`",
        f"- Detailed block scores: `{out_dir / 'scores_block_null.csv'}`",
        f"- Claims register table (multi-label): `{out_dir / 'claims_register.csv'}`",
        f"- Per-claim metrics: `{out_dir / 'claim_metrics.csv'}`",
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

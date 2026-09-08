"""Run the complete v0 spatial CCC prototype."""
from __future__ import annotations

import logging
import random
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURES_DIR, LOGS_DIR, PROCESSED_DIR, RAW_DIR, RESULTS_DIR, SEED
from figures import make_figures
from pipeline import add_spatial_graph, annotate_and_cluster, baseline_lr, load_data, quality_control


def main() -> None:
    random.seed(SEED)
    np.random.seed(SEED)
    sc.settings.seed = SEED
    for directory in (LOGS_DIR, RESULTS_DIR, PROCESSED_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / "prototype.log"
    logging.basicConfig(filename=log_path, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    started = datetime.now(timezone.utc).isoformat()
    logging.info("prototype start=%s seed=%d", started, SEED)
    adata = load_data(RAW_DIR)
    messages = [f"Loaded {adata.n_obs} cells and {adata.n_vars} genes"]
    quality_control(adata, messages)
    annotate_and_cluster(adata)
    add_spatial_graph(adata)
    lr = baseline_lr(adata)
    lr.to_csv(RESULTS_DIR / "baseline_lr_ranked.csv", index=False)
    adata.write_h5ad(PROCESSED_DIR / "xenium_breast_baseline.h5ad", compression="gzip")
    make_figures(adata, lr, FIGURES_DIR)
    for message in messages:
        logging.info(message)
    logging.info("spatial graph parameters=%s", adata.uns["spatial_graph_parameters"])
    logging.info("baseline interactions=%d", len(lr))
    print("Prototype complete")
    print("; ".join(messages))
    print(f"Results: {RESULTS_DIR / 'baseline_lr_ranked.csv'}")
    print(f"Processed: {PROCESSED_DIR / 'xenium_breast_baseline.h5ad'}")
    print(f"Figures: {FIGURES_DIR}")


if __name__ == "__main__":
    main()

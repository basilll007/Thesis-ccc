"""Shared configuration for deterministic prototype runs."""
from pathlib import Path

SEED = 42
ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
LOGS_DIR = ROOT / "logs"
CONTEXT_PATH = ROOT / "docs" / "PROJECT_CONTEXT.md"
PROGRESS_PATH = ROOT / "docs" / "PROGRESS_LOG.md"

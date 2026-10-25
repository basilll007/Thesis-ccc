"""Shared helpers for the PDAC-extension scripts.

Every script in src/pdac/ uses these helpers so that:
  * a run log is written to logs/<script>_<timestamp>.log with package versions and git commit,
  * every produced file is registered in results/pdac/MANIFEST.csv,
  * figures are saved as PNG + PDF under figures/pdac/<phase>/ and are never overwritten,
  * all paths are relative to the repository root (no absolute drive paths).
"""
from __future__ import annotations

import csv
import importlib
import json
import logging
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
RESULTS_PDAC = REPO_ROOT / "results" / "pdac"
FIGURES_PDAC = REPO_ROOT / "figures" / "pdac"
LOGS_DIR = REPO_ROOT / "logs"
DATA_DIR = REPO_ROOT / "data"
MANIFEST = RESULTS_PDAC / "MANIFEST.csv"
SEED = 42

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

_PACKAGES = [
    "numpy", "pandas", "scipy", "scanpy", "squidpy", "anndata", "sklearn",
    "matplotlib", "statsmodels", "leidenalg", "igraph",
]

_CURRENT_SCRIPT = None


def rel(path: Path | str) -> str:
    """Repository-relative POSIX path string (used in logs and the manifest)."""
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True, check=True)
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
        return out.stdout.strip() + ("+dirty" if dirty else "")
    except Exception as exc:  # pragma: no cover
        return f"unknown ({exc})"


def package_versions(extra: list[str] | None = None) -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for name in _PACKAGES + (extra or []):
        try:
            mod = importlib.import_module(name)
            versions[name] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[name] = "not installed"
    return versions


def start_run(script_name: str, extra_packages: list[str] | None = None) -> logging.Logger:
    """Create logs/<script>_<timestamp>.log and log environment provenance."""
    global _CURRENT_SCRIPT
    _CURRENT_SCRIPT = script_name
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    log_path = LOGS_DIR / f"{script_name}_{stamp}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info("script=%s log=%s", script_name, rel(log_path))
    logger.info("git_commit=%s", git_commit())
    logger.info("seed=%d", SEED)
    logger.info("package_versions=%s", json.dumps(package_versions(extra_packages)))
    logger.info("platform=%s", platform.platform())
    logger.log_path = log_path  # type: ignore[attr-defined]
    return logger


def register(path: Path | str, phase: str, description: str, script: str | None = None) -> None:
    """Add/update a row in results/pdac/MANIFEST.csv (keyed by path)."""
    RESULTS_PDAC.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    header = ["path", "phase", "description", "script", "registered_at"]
    if MANIFEST.exists():
        with MANIFEST.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
    key = rel(path)
    rows = [r for r in rows if r["path"] != key]
    rows.append({
        "path": key,
        "phase": phase,
        "description": description,
        "script": script or (f"src/pdac/{_CURRENT_SCRIPT}.py" if _CURRENT_SCRIPT else ""),
        "registered_at": datetime.now().isoformat(timespec="seconds"),
    })
    with MANIFEST.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def out_dir(phase: str) -> Path:
    d = RESULTS_PDAC / phase
    d.mkdir(parents=True, exist_ok=True)
    return d


def fig_dir(phase: str) -> Path:
    d = FIGURES_PDAC / phase
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_figure(fig, phase: str, name: str, description: str, logger: logging.Logger | None = None) -> Path:
    """Save PNG + PDF under figures/pdac/<phase>/; never overwrite (adds a timestamp suffix)."""
    import matplotlib.pyplot as plt

    d = fig_dir(phase)
    base = d / name
    if base.with_suffix(".png").exists() or base.with_suffix(".pdf").exists():
        base = d / f"{name}_{datetime.now().strftime('%Y%m%dT%H%M%S')}"
        if logger:
            logger.warning("Figure %s already existed; saving as %s instead (no overwrite).", name, base.name)
    fig.savefig(base.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    for suf in (".png", ".pdf"):
        register(base.with_suffix(suf), phase, description)
    if logger:
        logger.info("Saved figure %s(.png/.pdf)", rel(base))
    return base


def save_csv(df, phase: str, name: str, description: str, logger: logging.Logger | None = None, index: bool = False) -> Path:
    path = out_dir(phase) / name
    df.to_csv(path, index=index)
    register(path, phase, description)
    if logger:
        logger.info("Saved %s (%d rows)", rel(path), len(df))
    return path


def save_json(obj, phase: str, name: str, description: str, logger: logging.Logger | None = None) -> Path:
    path = out_dir(phase) / name
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, default=str)
    register(path, phase, description)
    if logger:
        logger.info("Saved %s", rel(path))
    return path


def save_npz(phase: str, name: str, description: str, logger: logging.Logger | None = None, **arrays) -> Path:
    import numpy as np

    path = out_dir(phase) / name
    np.savez_compressed(path, **arrays)
    register(path, phase, description)
    if logger:
        logger.info("Saved %s (arrays: %s)", rel(path), ", ".join(f"{k}{getattr(v, 'shape', '')}" for k, v in arrays.items()))
    return path


def finish_run(logger: logging.Logger, phase: str) -> None:
    """Register the run log itself in the manifest."""
    path = getattr(logger, "log_path", None)
    if path is not None:
        register(path, phase, "run log (package versions, git commit, parameters, key numbers)")
    logger.info("run finished")
    for h in list(logger.handlers):
        h.flush()

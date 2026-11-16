"""Shared utilities for model_dev: paths, config loading, run logs, checksums, seeding, bootstrap.

Why a separate module: every phase script must write the same provenance (package versions, GPU, commit hash,
seed) and must never write large artefacts into git. Keeping that in one place keeps the phase scripts short.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import logging
import os
import platform
import random
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_DEV = REPO_ROOT / "model_dev"
RESULTS = MODEL_DEV / "results"
CACHE = REPO_ROOT / "model_cache"  # git-ignored: datasets, checkpoints, temp
DATA_DIR = REPO_ROOT / "data"
SEED = 42

# Existing project code that we reuse (loaders, QC, program scoring, type-aware flags). Never edited here.
for _p in (REPO_ROOT / "src", REPO_ROOT / "src" / "pdac", REPO_ROOT / "src" / "chain"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_PACKAGES = ["numpy", "pandas", "scipy", "sklearn", "scanpy", "anndata", "torch", "torch_geometric", "yaml", "pytest"]


def load_config(path: str | Path) -> dict[str, Any]:
    """Read a YAML config; relative paths inside are resolved against the repository root by callers."""
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def git_commit() -> str:
    """Current commit hash, with '+dirty' if the working tree has uncommitted changes."""
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True,
                           check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "model_dev"], cwd=REPO_ROOT,
                               capture_output=True, text=True).stdout.strip()
        return h + ("+dirty" if dirty else "")
    except Exception as exc:  # pragma: no cover
        return f"unknown ({exc})"


def package_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for name in _PACKAGES:
        try:
            out[name] = getattr(importlib.import_module(name), "__version__", "unknown")
        except Exception:
            out[name] = "not installed"
    return out


def gpu_info() -> str:
    try:
        import torch
        if torch.cuda.is_available():
            p = torch.cuda.get_device_properties(0)
            return f"{p.name} ({p.total_memory / 2**20:.0f} MiB); torch {torch.__version__}; cuda {torch.version.cuda}"
        return f"no CUDA device; torch {torch.__version__}"
    except Exception:
        return "torch not installed"


def start_run(phase: str, name: str, config: dict | None = None) -> logging.Logger:
    """Create model_dev/results/<phase>/run_<name>_<timestamp>.log and log provenance."""
    d = RESULTS / phase
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"run_{name}_{datetime.now().strftime('%Y%m%dT%H%M%S')}.log"
    log = logging.getLogger(f"model_dev.{name}.{os.getpid()}")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    log.propagate = False
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for h in (logging.FileHandler(path, mode="w", encoding="utf-8"), logging.StreamHandler(sys.stdout)):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.info("phase=%s script=%s log=%s", phase, name, rel(path))
    log.info("git_commit=%s", git_commit())
    log.info("seed=%d", SEED)
    log.info("package_versions=%s", json.dumps(package_versions()))
    log.info("platform=%s gpu=%s", platform.platform(), gpu_info())
    if config is not None:
        log.info("config=%s", json.dumps(config, default=str))
    log.log_path = path  # type: ignore[attr-defined]
    return log


def rel(path: Path | str) -> str:
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def sha256(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def register_artefact(path: Path, phase: str, description: str) -> None:
    """Record a large artefact (outside git) in model_dev/results/phase0/artefacts.csv with its checksum."""
    import pandas as pd
    reg = RESULTS / "phase0" / "artefacts.csv"
    reg.parent.mkdir(parents=True, exist_ok=True)
    row = {"path": rel(path), "phase": phase, "description": description, "bytes": Path(path).stat().st_size,
           "sha256": sha256(path), "registered_at": datetime.now().isoformat(timespec="seconds")}
    df = pd.read_csv(reg) if reg.exists() else pd.DataFrame(columns=list(row))
    df = df[df.path != row["path"]]
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    df.to_csv(reg, index=False)


def save_table(df, phase: str, name: str, log: logging.Logger | None = None, overwrite: bool = True) -> Path:
    """Save a small table under model_dev/results/<phase>/. Small tables are versioned by git, so overwriting is
    allowed (history is in git); large arrays never go here."""
    d = RESULTS / phase
    d.mkdir(parents=True, exist_ok=True)
    path = d / name
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    df.to_csv(path, index=False)
    if log:
        log.info("saved %s (%d rows)", rel(path), len(df))
    return path


def save_json(obj, phase: str, name: str, log: logging.Logger | None = None) -> Path:
    d = RESULTS / phase
    d.mkdir(parents=True, exist_ok=True)
    path = d / name
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
    if log:
        log.info("saved %s", rel(path))
    return path


def set_seed(seed: int, deterministic: bool = True) -> None:
    """Seed python, numpy and torch (if present); request deterministic kernels where torch allows it."""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
            try:
                torch.use_deterministic_algorithms(True, warn_only=True)
            except Exception:
                pass
    except ImportError:
        pass


def tumour_bootstrap(per_group: dict[str, tuple[float, float]], n_boot: int = 2000, seed: int = SEED):
    """Pooled R2 and its percentile CI by resampling groups (tumours), not cells.

    per_group maps group -> (SSE, SST). Pooled R2 = 1 - sum SSE / sum SST over the resampled groups.
    """
    keys = sorted(per_group)
    sse = np.array([per_group[k][0] for k in keys], float)
    sst = np.array([per_group[k][1] for k in keys], float)
    point = 1 - sse.sum() / sst.sum()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), size=(n_boot, len(keys)))
    boots = 1 - sse[idx].sum(1) / sst[idx].sum(1)
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return float(point), float(lo), float(hi)


def tumour_bootstrap_diff(a: dict[str, tuple[float, float]], b: dict[str, tuple[float, float]], n_boot: int = 2000,
                          seed: int = SEED):
    """R2(a) - R2(b) in pooled R2 with a paired tumour bootstrap; a and b share SST per group.

    Pass the model of interest FIRST: tumour_bootstrap_diff(model_pairs, reference_pairs) is the gain of the model
    over the reference (tests/test_metrics.py pins the sign)."""
    keys = sorted(set(a) & set(b))
    sa = np.array([a[k][0] for k in keys]); sb = np.array([b[k][0] for k in keys])
    sst = np.array([a[k][1] for k in keys])
    point = (sb.sum() - sa.sum()) / sst.sum()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), size=(n_boot, len(keys)))
    boots = (sb[idx].sum(1) - sa[idx].sum(1)) / sst[idx].sum(1)
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return float(point), float(lo), float(hi)

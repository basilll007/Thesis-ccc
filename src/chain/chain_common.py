"""Shared helpers for the chain-v1 scripts (src/chain/).

Same contract as src/pdac/common.py but writes to the chain locations:
  * run log logs/chain_<script>_<timestamp>_<pid>.log with package versions and git commit,
  * every produced file registered in results/chain/MANIFEST.csv,
  * figures as PNG + PDF under figures/chain/<phase>/ (never overwritten),
  * results never overwritten: save_* refuse to replace an existing file,
  * repository-relative paths only. Seed 42.
"""
from __future__ import annotations

import csv
import importlib
import json
import logging
import os
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
RESULTS = REPO_ROOT / "results" / "chain"
FIGURES = REPO_ROOT / "figures" / "chain"
LOGS_DIR = REPO_ROOT / "logs"
DATA_DIR = REPO_ROOT / "data"
CACHE = DATA_DIR / "chain_cache"   # git-ignored cell-level intermediates (deterministic, re-derivable)
MANIFEST = RESULTS / "MANIFEST.csv"
SEED = 42
os.environ.setdefault("NUMBA_CACHE_DIR", str(DATA_DIR / "pdac" / "numba_cache"))

for _p in (SRC_DIR, SRC_DIR / "pdac", SRC_DIR / "chain"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_PACKAGES = ["numpy", "pandas", "scipy", "scanpy", "anndata", "sklearn", "matplotlib", "statsmodels"]
_CURRENT = None


def rel(path) -> str:
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def git_commit() -> str:
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True,
                           check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT, capture_output=True,
                               text=True).stdout.strip()
        return h + ("+dirty" if dirty else "")
    except Exception as exc:  # pragma: no cover
        return f"unknown ({exc})"


def versions() -> dict:
    v = {"python": platform.python_version()}
    for n in _PACKAGES:
        try:
            v[n] = getattr(importlib.import_module(n), "__version__", "unknown")
        except Exception:
            v[n] = "not installed"
    return v


def start_run(name: str, params: dict | None = None) -> logging.Logger:
    global _CURRENT
    _CURRENT = name
    LOGS_DIR.mkdir(exist_ok=True)
    path = LOGS_DIR / f"chain_{name}_{datetime.now().strftime('%Y%m%dT%H%M%S')}_{os.getpid()}.log"
    log = logging.getLogger(f"chain.{name}.{os.getpid()}")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    log.propagate = False
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for h in (logging.FileHandler(path, mode="w", encoding="utf-8"), logging.StreamHandler(sys.stdout)):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.info("script=src/chain/%s.py log=%s", name, rel(path))
    log.info("git_commit=%s", git_commit())
    log.info("seed=%d", SEED)
    log.info("package_versions=%s", json.dumps(versions()))
    log.info("platform=%s", platform.platform())
    if params:
        log.info("params=%s", json.dumps(params, default=str))
    log.log_path = path  # type: ignore[attr-defined]
    return log


def register(path, phase: str, description: str) -> None:
    """Add/update a row of results/chain/MANIFEST.csv (keyed by path). File-locked for parallel workers."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    lock = RESULTS / ".manifest.lock"
    fd = None
    for _ in range(2400):
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except FileExistsError:
            time.sleep(0.05)
    if fd is None:
        raise RuntimeError("manifest lock timeout")
    try:
        header = ["path", "phase", "description", "script", "registered_at"]
        rows = []
        if MANIFEST.exists():
            with MANIFEST.open(newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
        key = rel(path)
        rows = [r for r in rows if r["path"] != key]
        rows.append({"path": key, "phase": phase, "description": description,
                     "script": f"src/chain/{_CURRENT}.py" if _CURRENT else "",
                     "registered_at": datetime.now().isoformat(timespec="seconds")})
        with MANIFEST.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            w.writerows(rows)
    finally:
        os.close(fd)
        os.remove(lock)


def out_dir(phase: str) -> Path:
    d = RESULTS / phase
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fresh(path: Path, overwrite_ok: bool) -> None:
    if path.exists() and not overwrite_ok:
        raise FileExistsError(f"{rel(path)} exists; earlier outputs are never overwritten (use a new name)")


def save_csv(df, phase, name, description, log=None, index=False, overwrite_ok=False) -> Path:
    path = out_dir(phase) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    _fresh(path, overwrite_ok)
    df.to_csv(path, index=index)
    register(path, phase, description)
    if log:
        log.info("Saved %s (%d rows)", rel(path), len(df))
    return path


def save_json(obj, phase, name, description, log=None, overwrite_ok=False) -> Path:
    path = out_dir(phase) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    _fresh(path, overwrite_ok)
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
    register(path, phase, description)
    if log:
        log.info("Saved %s", rel(path))
    return path


def save_figure(fig, phase, name, description, log=None) -> Path:
    import matplotlib.pyplot as plt
    d = FIGURES / phase
    d.mkdir(parents=True, exist_ok=True)
    base = d / name
    if base.with_suffix(".png").exists() or base.with_suffix(".pdf").exists():
        base = d / f"{name}_{datetime.now().strftime('%Y%m%dT%H%M%S')}"
    fig.savefig(base.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    for s in (".png", ".pdf"):
        register(base.with_suffix(s), phase, description)
    if log:
        log.info("Saved figure %s(.png/.pdf)", rel(base))
    return base


def finish_run(log, phase: str) -> None:
    p = getattr(log, "log_path", None)
    if p is not None:
        register(p, phase, "run log (package versions, git commit, parameters, key numbers)")
    log.info("run finished")
    for h in list(log.handlers):
        h.flush()
        h.close()


def mem_status() -> str:
    from common import mem_status as m  # src/pdac/common.py
    return m()


def add_one_p(null, obs: float, two_sided: bool = True) -> float:
    """Add-one empirical p = (k+1)/(n+1). Two-sided compares |null| >= |obs| (statistics are centred at 0)."""
    import numpy as np
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    k = int((np.abs(null) >= abs(obs)).sum()) if two_sided else int((null >= obs).sum())
    return float((k + 1) / (null.size + 1))


def bh(p):
    import numpy as np
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    if v.size:
        o = np.argsort(v)
        adj = np.minimum.accumulate((v[o] * v.size / np.arange(1, v.size + 1))[::-1])[::-1]
        out = np.empty_like(v)
        out[o] = np.minimum(adj, 1)
        q[ok] = out
    return q

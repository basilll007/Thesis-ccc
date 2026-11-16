"""Assemble receiver-level feature matrices and targets across sections, with caching.

What: for a split ('development' or 'external') and a receiver, builds the Phase 1 feature blocks for every
section with `features.build_section_features`, concatenates them, records section/group/row indices, and caches
the result in model_cache/features/ (git-ignored, registered with a checksum). Also loads the fold table.

Why: baselines, the MLP, ablations and the external test must all see byte-identical inputs; building once and
caching removes any chance of drift between phases.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import features as F  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DATASET_DIR = {"development": "mendeley_pdac_smi", "external": "geo_gse310352"}


def section_files(split: str) -> list[Path]:
    return sorted((C.CACHE / "datasets" / DATASET_DIR[split]).glob("*.npz"))


def load_section(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


def targets_spec() -> dict:
    return json.loads((C.RESULTS / "phase1" / "targets.json").read_text())


def all_classes(split: str) -> list[str]:
    """Annotation classes present in the split (stable, sorted); the same list is used for every section."""
    counts = pd.read_csv(C.RESULTS / "phase0" / "cell_counts.csv").query("split == @split")
    return sorted(counts.cell_type.unique())


@dataclass
class ReceiverData:
    X: dict[str, np.ndarray]          # block name -> matrix
    names: dict[str, list[str]]
    prov: dict[str, list[F.Provenance]]
    Y: np.ndarray                      # targets (programs then genes)
    target_names: list[str]
    section: np.ndarray
    group: np.ndarray                  # tumour (development) or unit (external)
    row: np.ndarray                    # index of the cell within its section file

    def matrix(self, blocks: list[str]) -> np.ndarray:
        return np.column_stack([self.X[b] for b in blocks])

    def feature_set(self, blocks: list[str]) -> F.FeatureSet:
        fs = F.FeatureSet(rows=np.arange(len(self.Y)))
        for b in blocks:
            fs.add(b, self.X[b], self.names[b], self.prov[b])
        return fs


def build_receiver(split: str, receiver: str, log=None, use_cache: bool = True) -> ReceiverData:
    spec = targets_spec()[receiver]
    classes = all_classes(split)
    cache = C.CACHE / "features" / f"{split}__{receiver}.npz"
    meta = cache.with_suffix(".json")
    if use_cache and cache.exists() and meta.exists():
        z = np.load(cache, allow_pickle=False)
        m = json.loads(meta.read_text())
        X = {b: z[f"X__{b}"] for b in m["blocks"]}
        prov = {b: [F.Provenance(p["source"], frozenset(p["genes"]), frozenset(p["programs"])) for p in m["prov"][b]]
                for b in m["blocks"]}
        return ReceiverData(X, m["names"], prov, z["Y"], m["target_names"], z["section"], z["group"], z["row"])
    parts = []
    for path in section_files(split):
        d = load_section(path)
        fs, y = F.build_section_features(d, spec["class"], spec["programs"], spec["genes"], classes)
        if len(fs.rows) == 0:
            continue
        F.assert_no_target_leak(fs, None, y, set(spec["genes"]), set(spec["programs"]))
        parts.append((str(d["section"]), str(d["group"]), fs, y))
        if log:
            log.info("features %s %s: %d receiver cells, %d columns", split, d["section"], len(fs.rows),
                     fs.matrix().shape[1])
    blocks = list(parts[0][2].blocks)
    X = {b: np.vstack([p[2].blocks[b] for p in parts]) for b in blocks}
    Y = np.vstack([p[3] for p in parts]).astype(np.float32)
    section = np.concatenate([[p[0]] * len(p[2].rows) for p in parts])
    group = np.concatenate([[p[1]] * len(p[2].rows) for p in parts])
    row = np.concatenate([p[2].rows for p in parts])
    names = parts[0][2].names
    prov = parts[0][2].prov
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, Y=Y, section=section.astype("U"), group=group.astype("U"), row=row,
                        **{f"X__{b}": X[b] for b in blocks})
    meta.write_text(json.dumps({"blocks": blocks, "names": names, "target_names": spec["programs"] + spec["genes"],
                                "prov": {b: [{"source": p.source, "genes": sorted(p.genes), "programs": sorted(p.programs)}
                                             for p in prov[b]] for b in blocks}}))
    C.register_artefact(cache, "phase2", f"{split} {receiver} receiver feature blocks and targets")
    return ReceiverData(X, names, prov, Y, spec["programs"] + spec["genes"], section, group, row)


def folds(kind: str = "grouped") -> pd.DataFrame:
    name = "folds.csv" if kind == "grouped" else "folds_loto_reproduction_only.csv"
    return pd.read_csv(C.RESULTS / "phase1" / name)


def fold_masks(data: ReceiverData, fold_table: pd.DataFrame, k: int, n_folds: int):
    """Boolean masks (train_inner, val_inner, train_all, test) for outer fold k; inner val = fold (k+1) mod n."""
    f = data.section.copy()
    sec2fold = dict(zip(fold_table.section, fold_table.fold))
    fold_of = np.array([sec2fold[s] for s in f])
    test = fold_of == k
    val = fold_of == ((k + 1) % n_folds)
    train_all = ~test
    train_inner = train_all & ~val
    return train_inner, val, train_all, test

"""Phase 1: patient-level fold assignment and target declarations.

What: assigns the 13 development tumours to 5 grouped folds (all sections of a tumour together), balanced on
myeloid (primary receiver) cell count, and writes model_dev/results/phase1/folds.csv plus targets.json. Also writes
the per-fold receiver counts against the declared floor and the leave-one-tumour-out table for the reproduction check.

Why: one fold file, written once and reused by every later phase, is what makes the comparisons between
baselines, graph models and ablations paired and leakage-free at the patient level.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase1"


def grouped_folds(counts: pd.Series, n_folds: int) -> dict[str, int]:
    """Greedy balanced assignment: largest group first, into the fold with the smallest running total."""
    order = sorted(counts.index, key=lambda t: (-counts[t], t))
    totals = np.zeros(n_folds)
    out = {}
    for t in order:
        k = int(np.argmin(totals))
        out[t] = k
        totals[k] += counts[t]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "splits.yaml"))
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "splits", cfg)
    counts = pd.read_csv(C.RESULTS / "phase0" / "cell_counts.csv").query("split == 'development'")
    sec_info = pd.read_csv(C.RESULTS / "phase0" / "section_info.csv").query("split == 'development'")
    prim_cls = cfg["receivers"][cfg["primary_receiver"]]["class"]
    per_tumour = counts[counts.cell_type == prim_cls].groupby("group").n_cells.sum()
    per_tumour = per_tumour.reindex(sorted(sec_info.group.unique()), fill_value=0)
    assign = grouped_folds(per_tumour, cfg["n_folds"])
    folds = sec_info[["section", "group"]].rename(columns={"group": "tumour"}).copy()
    folds["fold"] = folds.tumour.map(assign)
    folds["inner_val_fold_when_test_is"] = ""  # documented rule: inner val = (k + 1) mod n_folds of the training folds
    folds = folds.sort_values(["fold", "tumour", "section"])
    C.save_table(folds.drop(columns="inner_val_fold_when_test_is"), PHASE, "folds.csv", log)
    loto = folds[["section", "tumour"]].copy()
    loto["fold"] = loto.tumour.map({t: i for i, t in enumerate(sorted(loto.tumour.unique()))})
    C.save_table(loto, PHASE, "folds_loto_reproduction_only.csv", log)
    # floor check per receiver and fold
    rows = []
    for rname, r in cfg["receivers"].items():
        for k in range(cfg["n_folds"]):
            tums = [t for t, f in assign.items() if f == k]
            n = int(counts[(counts.cell_type == r["class"]) & counts.group.isin(tums)].n_cells.sum())
            rows.append({"receiver": rname, "class": r["class"], "fold": k, "tumours": ";".join(sorted(tums)),
                         "n_tumours": len(tums), "n_receiver_cells": n,
                         "passes_floor": bool(n >= cfg["floor_cells_per_fold"] and len(tums) >= cfg["floor_tumours_per_fold"])})
    fl = pd.DataFrame(rows)
    C.save_table(fl, PHASE, "fold_floor_check.csv", log)
    tcells = int(counts[counts.cell_type == "T_cell"].n_cells.sum())
    programs = json.loads((C.REPO_ROOT / cfg["programs"]).read_text())["programs"]
    targets = {}
    for rname, r in cfg["receivers"].items():
        genes = [g for p in r["programs"] for g in programs[p]]
        targets[rname] = {"class": r["class"], "role": r["role"], "programs": r["programs"], "genes": genes,
                          "n_genes": len(genes), "testable": bool(fl[fl.receiver == rname].passes_floor.all())}
    targets["T"] = {"class": "T_cell", "role": "not testable by declaration (too few cells)", "n_cells_total": tcells,
                    "testable": False}
    targets["_primary_endpoint"] = cfg["primary_endpoint"]
    targets["_inner_validation_rule"] = "inner val = training fold (k + 1) mod n_folds; inner train = other 3 folds"
    C.save_json(targets, PHASE, "targets.json", log)
    log.info("fold assignment: %s", json.dumps(assign))
    log.info("floor check:\n%s", fl.to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()

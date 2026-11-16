"""Phase 3 ablations: same-type only, other-type only, composition only, shuffled neighbours, no graph.

What: reruns the selected architecture/learning rate (results/phase3/architecture_selection.json) on the same
folds and seeds with each ablation mask, and reports the cross-type contribution = R2(full) - R2(same-type only)
with a paired tumour bootstrap, plus the other declared differences.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import train as T  # noqa: E402

import pandas as pd  # noqa: E402

PHASE = "phase3"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "ablations.yaml"))
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "ablations", cfg)
    sel = json.loads((C.RESULTS / PHASE / "architecture_selection.json").read_text())
    arch, lr = sel["selected_arch"], float(sel["selected_lr"])
    log.info("selected architecture %s, lr %g", arch, lr)
    modes = sorted({m for v in cfg["modes_by_receiver"].values() for m in v})
    res, pairs = T.run_experiments(cfg, modes, [arch], [lr], PHASE, "ablation", log, receivers_modes=cfg["modes_by_receiver"])
    C.save_table(res, PHASE, "ablation_results.csv", log)
    # full-model pairs from the Phase 3 main run (same arch/lr/seeds/folds)
    rows = []
    for receiver in cfg["receivers"]:
        fp = pd.read_csv(C.RESULTS / PHASE / f"per_tumour_sse_graph_{receiver}.csv")
        fp = fp[(fp["mode"] == "full") & (fp.arch == arch) & (fp.lr == lr)]
        for t in sorted(res[res.receiver == receiver].target.unique()):
            full = {g: (float(x.sse.mean()), float(x.sst.iloc[0])) for g, x in fp[fp.target == t].groupby("tumour")}
            r2_full = C.tumour_bootstrap(full, cfg["n_boot"])
            for mode in cfg["modes_by_receiver"][receiver]:
                ab = T.seed_mean_pairs(pairs, receiver, mode, arch, lr, cfg["seeds"], t)
                if not ab:
                    continue
                r2, lo, hi = C.tumour_bootstrap(ab, cfg["n_boot"])
                diff, dlo, dhi = C.tumour_bootstrap_diff(ab, full, cfg["n_boot"])   # full - ablation
                rows.append({"receiver": receiver, "target": t, "arch": arch, "lr": lr, "mode": mode, "r2_mode": r2,
                             "ci_low": lo, "ci_high": hi, "r2_full": r2_full[0], "full_minus_mode": diff,
                             "diff_ci_low": dlo, "diff_ci_high": dhi, "diff_ci_excludes_0": bool(dlo > 0 or dhi < 0),
                             "quantity": {"same_type": "cross-type contribution (full - same-type only)",
                                          "other_type": "full - other-type only", "composition": "full - composition only",
                                          "shuffled": "full - shuffled neighbour expression (zonation control)",
                                          "no_graph": "full - no graph"}[mode]})
    summ = pd.DataFrame(rows)
    C.save_table(summ, PHASE, "ablation_summary.csv", log)
    log.info("ablation summary:\n%s", summ.round(4).to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()

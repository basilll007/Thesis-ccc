"""Chain Phase 5b: de-duplicated view of results/chain/ranked_hypotheses.csv.

ranked_hypotheses.csv lists each hop once per chain program context (e.g. a B<-A hop repeated for every C program),
without a column naming that context, so rows look duplicated. The original file is kept unmodified (never overwrite);
this writes ranked_hypotheses_dedup.csv with one row per (chain, hop, sender/receiver program, ligand, receptor).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402

import pandas as pd  # noqa: E402


def main():
    log = cc.start_run("chain_phase5b_dedup")
    df = pd.read_csv(cc.RESULTS / "ranked_hypotheses.csv")
    key = ["chain", "hop", "sender_program", "receiver_program", "candidate_ligand", "candidate_receptor"]
    d = df.drop_duplicates(key).copy()
    d = d.sort_values(["hop_support", "hop_effect"], key=lambda s: s.abs() if s.name == "hop_effect" else s,
                      ascending=[False, False])
    path = cc.RESULTS / "ranked_hypotheses_dedup.csv"
    if path.exists():
        raise FileExistsError(path)
    d.to_csv(path, index=False)
    cc.register(path, "phase5", "De-duplicated ranked_hypotheses (one row per hop x LR pair); original kept unmodified")
    log.info("rows: original %d -> dedup %d; supported hops: %d", len(df), len(d), int((d.hop_support == "yes").sum()))
    cc.finish_run(log, "phase5")


if __name__ == "__main__":
    main()

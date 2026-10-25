"""Phase 1b: detailed gate breakdown computed from the Phase 1 tables (no new permutations).

Reports, per family:
  * C/E disagreements split by ligand-receptor pair, sparse-gene flag and low edge count (< 30 qualifying edges),
  * for every edge-significant interaction, whether the edge/matched ratio Q is also significant
    (E = M x Q: if Q is within its null range, E's significance is carried by the matched non-spatial baseline M),
  * where Q is significant, the matched-baseline z (to expose small-denominator ratio behaviour).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import finish_run, out_dir, save_csv, save_figure, save_json, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase1"
LOW_EDGES = 30


def main() -> None:
    log = start_run("phase1b_gate_details")
    df = pd.read_csv(out_dir(PHASE) / "breast_matched_baseline_all_families.csv")
    out = {}
    rows = []
    for fam, d in df.groupby("family", sort=False):
        ev = d[d["edge_status"] == "evaluated"]
        dis = ev[ev["comp_sig"] != ev["edge_sig"]]
        esig = ev[ev["edge_sig"]]
        qsig = ev[ev["ratio_sig"]]
        info = {
            "n_interactions": int(len(d)),
            "n_C_E_disagree": int(len(dis)),
            "disagree_by_pair": dis.groupby(dis["ligand"] + "-" + dis["receptor"]).size().to_dict(),
            "disagree_score_form_by_pair": dis[dis["decomposition"].str.contains("SCORE-FORM")]
                .groupby(lambda i: dis.loc[i, "ligand"] + "-" + dis.loc[i, "receptor"]).size().to_dict(),
            "disagree_sparse_gene": int(dis["sparse_gene_diagnostic"].sum()),
            "disagree_nonsparse": int((~dis["sparse_gene_diagnostic"]).sum()),
            "disagree_nonsparse_score_form": int(((~dis["sparse_gene_diagnostic"]) & dis["decomposition"].str.contains("SCORE-FORM")).sum()),
            "disagree_nonsparse_proximity": int(((~dis["sparse_gene_diagnostic"]) & dis["decomposition"].str.contains("PROXIMITY")).sum()),
            "disagree_proximity_with_lt30_edges": int((dis["decomposition"].str.contains("PROXIMITY") & (dis["n_qualifying_edges"] < LOW_EDGES)).sum()),
            "n_edge_sig": int(len(esig)),
            "n_edge_sig_with_ratio_sig": int(esig["ratio_sig"].sum()),
            "n_edge_sig_with_matched_sig": int(esig["matched_sig"].sum()),
            "n_ratio_sig": int(len(qsig)),
            "ratio_sig_rows": [f"{r.ligand}-{r.receptor} {r.sender}->{r.receiver} (Q={r.ratio_obs:.2f}, z_Q={r.ratio_z:.2f}, "
                               f"z_M={r.matched_z:.2f}, z_E={r.edge_z:.2f})" for r in qsig.itertuples()],
            "n_ratio_sig_with_negative_matched_z": int((qsig["matched_z"] < 0).sum()),
        }
        out[fam] = info
        log.info("[%s] %s", fam, info)
        for r in dis.itertuples():
            rows.append({"family": fam, "interaction": f"{r.ligand}-{r.receptor} {r.sender}->{r.receiver}",
                         "n_qualifying_edges": r.n_qualifying_edges, "comp_z": r.comp_z, "matched_z": r.matched_z,
                         "edge_z": r.edge_z, "ratio_z": r.ratio_z, "comp_sig": r.comp_sig, "matched_sig": r.matched_sig,
                         "edge_sig": r.edge_sig, "decomposition": r.decomposition,
                         "sparse_gene_diagnostic": r.sparse_gene_diagnostic, "lt30_edges": r.n_qualifying_edges < LOW_EDGES})
    dis_df = pd.DataFrame(rows)
    save_csv(dis_df, PHASE, "disagreements_detail.csv",
             "Every C/E verdict disagreement with matched/ratio z, decomposition, sparse-gene and <30-edge flags", log)
    save_json(out, PHASE, "phase1_gate_details.json",
              "Gate breakdown: disagreements by pair/sparsity/edge count; edge-significant rows vs ratio significance", log)

    # Figure: for edge-significant rows, matched z vs ratio z
    esig = df[(df["edge_status"] == "evaluated") & df["edge_sig"]].copy()
    esig["lab"] = esig["family"] + ": " + esig["ligand"] + "→" + esig["receptor"] + " " + esig["sender"] + "→" + esig["receiver"]
    esig = esig.sort_values("edge_z")
    y = np.arange(len(esig))
    fig, ax = plt.subplots(figsize=(10, 0.35 * len(esig) + 1.5))
    ax.barh(y - 0.2, esig["matched_z"].clip(upper=40), height=0.4, color="#1f77b4", label="matched baseline M z (clipped at 40)")
    ax.barh(y + 0.2, esig["ratio_z"], height=0.4, color="#2ca02c", label="ratio Q = E/M z")
    ax.axvline(1.645, color="#999", ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(esig["lab"], fontsize=8)
    ax.set_xlabel("z vs label-permutation null")
    ax.set_title("Edge-significant interactions: is the edge signal carried by M (composition) or by Q (proximity coupling)?",
                 fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, PHASE, "phase1b_edge_significant_matched_vs_ratio_z",
                "For edge-significant interactions, matched-baseline z vs edge/matched-ratio z", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

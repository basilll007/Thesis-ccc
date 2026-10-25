"""Phase 5: breast vs PDAC per-axis synthesis table (reads Phase 1 and Phase 3 outputs only).

Axis mapping (declared here; breast labels are the 4 coarse breast types, PDAC uses marker labels for the
primary section and marker/author labels across 16 sections):
  CXCL12->CXCR4 stroma->immune : breast fibroblast->immune      | PDAC CAF_fibroblast->T_cell
  CXCL12->CXCR4 stroma->tumour : breast fibroblast->epithelial  | PDAC CAF_fibroblast->epithelial_ductal_tumor
  CD274->PDCD1 tumour->immune  : breast epithelial->immune      | PDAC epithelial_ductal_tumor->T_cell
  CD86->CTLA4 APC/immune->T    : breast immune->immune          | PDAC myeloid->T_cell
  PTN->SDC4 stroma->tumour     : breast fibroblast->epithelial  | PDAC (PTN, SDC4 not on panel)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO_ROOT, finish_run, out_dir, save_csv, save_figure, start_run  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase5"
AXES = [
    ("CXCL12", "CXCR4", "stroma->immune", ("fibroblast", "immune"), ("CAF_fibroblast", "T_cell")),
    ("CXCL12", "CXCR4", "stroma->tumour", ("fibroblast", "epithelial"), ("CAF_fibroblast", "epithelial_ductal_tumor")),
    ("CD274", "PDCD1", "tumour->immune", ("epithelial", "immune"), ("epithelial_ductal_tumor", "T_cell")),
    ("CD86", "CTLA4", "APC/immune->T", ("immune", "immune"), ("myeloid", "T_cell")),
    ("PTN", "SDC4", "stroma->tumour", ("fibroblast", "epithelial"), (None, None)),
]
COLS = ["n_qualifying_edges", "comp_z", "matched_z", "edge_z", "ratio_obs", "ratio_z", "comp_sig", "matched_sig", "edge_sig", "ratio_sig"]


def main() -> None:
    log = start_run("phase5_synthesis_tables")
    p1 = pd.concat([pd.read_csv(out_dir("phase1") / "breast_matched_baseline_paper22.csv"),
                    pd.read_csv(out_dir("phase1") / "breast_matched_baseline_extra_CD86CTLA4_PTNSDC4.csv")])
    p3 = pd.read_csv(out_dir("phase3") / "mendeley_pdac_smi__U5-a_6NN_scores_both_annotations.csv")
    p3 = p3[p3.annotation == "cell_type_marker"]
    xs = pd.read_csv(out_dir("phase3") / "cross_section_summary_by_interaction.csv")
    xs = xs[xs.annotation == "cell_type_marker"]
    au = pd.read_csv(out_dir("phase3") / "author_annotation_6NN_all_sections.csv",
                     usecols=["section", "ligand", "receptor", "sender", "receiver", "edge_sig", "ratio_sig"])
    diag = pd.read_csv(REPO_ROOT / "results" / "cd274_pdcd1_detection_diagnostic.csv")
    rows = []
    for l, r, kind, (bs, br), (ps, pr) in AXES:
        rec = {"axis": f"{l}->{r}", "kind": kind, "breast_sender_receiver": f"{bs}->{br}",
               "pdac_sender_receiver": f"{ps}->{pr}" if ps else "not on panel"}
        b = p1[(p1.ligand == l) & (p1.receptor == r) & (p1.sender == bs) & (p1.receiver == br)]
        if len(b):
            for c in COLS:
                rec[f"breast_{c}"] = b.iloc[0][c]
            rec["breast_note"] = "6-NN, 500 perms (Phase 1)"
        else:
            note = "not in tested set"
            if (l, r) == ("CD274", "PDCD1"):
                n = int(diag.query("gene == 'PDCD1' and cell_type == 'immune'")["n_detected"].iloc[0])
                tot = int(diag.query("gene == 'PDCD1' and cell_type == 'immune'")["n_cells"].iloc[0])
                note = f"not evaluable: PDCD1 detected in {n}/{tot} breast immune cells (compositional score 0, excluded by baseline)"
            rec["breast_note"] = note
        if ps:
            pp = p3[(p3.ligand == l) & (p3.receptor == r) & (p3.sender == ps) & (p3.receiver == pr)]
            if len(pp):
                for c in COLS:
                    rec[f"pdac_U5a_{c}"] = pp.iloc[0][c]
                rec["pdac_U5a_edge_q_bh_all_rows"] = pp.iloc[0]["edge_q_bh_all_rows"]
                rec["pdac_U5a_ratio_q_bh_all_rows"] = pp.iloc[0]["ratio_q_bh_all_rows"]
                rec["pdac_ligand_n_cells_U5a"] = pp.iloc[0]["ligand_n_cells_detected"]
                rec["pdac_receptor_n_cells_U5a"] = pp.iloc[0]["receptor_n_cells_detected"]
            x = xs[(xs.ligand == l) & (xs.receptor == r) & (xs.sender == ps) & (xs.receiver == pr)]
            if len(x):
                for c in ("n_sections", "n_comp_sig", "n_matched_sig", "n_edge_sig", "n_ratio_sig", "n_edge_sig_bh", "n_ratio_sig_bh",
                          "median_edge_z", "median_ratio_z", "median_ratio_obs"):
                    rec[f"pdac_16sec_marker_{c}"] = x.iloc[0][c]
            a = au[(au.ligand == l) & (au.receptor == r) & (au.sender == ps) & (au.receiver == pr)]
            rec["pdac_16sec_author_n_edge_sig"] = int(a.edge_sig.sum()) if len(a) else np.nan
            rec["pdac_16sec_author_n_ratio_sig"] = int(a.ratio_sig.sum()) if len(a) else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    save_csv(t, PHASE, "breast_vs_pdac_axes.csv", "Per-axis breast (Phase 1) vs PDAC (Phase 3) synthesis table", log)

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    lab = [f"{a} ({k})" for a, k in zip(t.axis, t.kind)]
    y = np.arange(len(t))
    for ax, pre, title in ((axes[0], "breast_", "Breast Xenium (6-NN)"), (axes[1], "pdac_U5a_", "PDAC CosMx U5-a (6-NN, marker)")):
        for key, c, mk, off in (("comp", "#7f7f7f", "o", -0.27), ("matched", "#1f77b4", "D", -0.09),
                                ("edge", "#d62728", "s", 0.09), ("ratio", "#2ca02c", "^", 0.27)):
            col = f"{pre}{key}_z"
            if col in t:
                ax.scatter(t[col], y + off, c=c, marker=mk, s=40, label=key)
        ax.axvline(1.645, color="#aaa", ls=":")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xscale("symlog", linthresh=2)
        ax.set_yticks(y)
        ax.set_yticklabels(lab, fontsize=8)
        ax.set_title(title)
        ax.legend(fontsize=7)
    fig.suptitle("Breast vs PDAC per-axis z (missing marker = axis not evaluable / not on panel)")
    fig.tight_layout()
    save_figure(fig, PHASE, "phase5_breast_vs_pdac_axes_z", "Per-axis z comparison breast vs PDAC", log)
    finish_run(log, PHASE)


if __name__ == "__main__":
    main()

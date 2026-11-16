"""Chain Phase 5: name candidate ligand-receptor signals for chain hops (prior knowledge; 'candidate', not evidence).

Plan: docs/chain_plans/phase5_candidate_signals.md (committed before running).
Deliverable: results/chain/ranked_hypotheses.csv.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_common as cc  # noqa: E402
import chain_features as cf  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PHASE = "phase5"
KEY = {v: k for k, v in cf.TYPES.items()}


def main():
    log = cc.start_run("chain_phase5_candidates")
    from lr_catalog import catalogue
    tl = pd.read_csv(cc.RESULTS / "phase3" / "tumour_level_tests_author_mapped.csv").query("variant == 'forward'")
    calls = pd.read_csv(cc.RESULTS / "phase3" / "chain_calls_author_mapped.csv")
    sec = pd.read_csv(cc.RESULTS / "phase3" / "section_level_stats_author_mapped.csv")
    cal = json.loads((cc.RESULTS / "phase4" / "calibrated_range.json").read_text())
    ta = pd.read_csv(cc.RESULTS / "phase1b" / "gene_background_typeaware_flags.csv")
    usable = {t: set(ta[(ta.type == t) & ~ta.background_level_for_type].gene) for t in cf.TYPES}
    progs = cf.programs_used()
    rows = []
    for st, hop, s_role, r_role in (("beta_BA", "B<-A", "A", "B"), ("beta_CB", "C<-B", "B", "C")):
        d = tl[tl.stat == st]
        for r in d.drop_duplicates(["chain", "prog_A", "prog_B", "prog_C"]).itertuples():
            ch = r.chain.split("->")
            sender = ch[0] if s_role == "A" else ch[1]
            receiver = ch[1] if r_role == "B" else ch[2]
            sprog = r.prog_A if s_role == "A" else r.prog_B
            rprog = r.prog_B if r_role == "B" else r.prog_C
            support = bool(r.q_expr < .05 and r.p_label < .05)
            # tumours supporting: tumours whose estimate has the tumour-level sign
            ss = sec[(sec.variant == "forward") & (sec.chain == r.chain) & (sec.prog_A == r.prog_A) &
                     (sec.prog_B == r.prog_B) & (sec.prog_C == r.prog_C) & (sec.stat == st)]
            te = ss.groupby("tumour").obs.mean()
            n_same = int((np.sign(te) == np.sign(r.estimate)).sum())
            cl = calls[(calls.chain == r.chain) & (calls.prog_A == r.prog_A) & (calls.prog_B == r.prog_B) &
                       (calls.prog_C == r.prog_C)]
            chain_label = cl.label.iloc[0] if len(cl) else "NA"
            if not cal.get("trustworthy_setting_found"):
                cstat = "outside calibrated range (no trustworthy setting in Phase 4)"
            else:
                cstat = ("within calibrated range" if (len(cl) and bool(cl.c5_calibrated.iloc[0])) else
                         "outside calibrated range")
                if r.chain != "tumour->CAF->T":
                    cstat += " (extrapolated from primary-chain calibration)"
            for lig, rec, src in catalogue():
                if lig not in usable.get(sender, set()) or rec not in usable.get(receiver, set()):
                    continue
                in_prog = lig in progs.get(sprog, [])
                rows.append({"chain": r.chain, "hop": hop, "sender": sender, "receiver": receiver,
                             "sender_program": sprog, "receiver_program": rprog, "candidate_ligand": lig,
                             "candidate_receptor": rec, "catalogue_source": src, "ligand_in_sender_program": in_prog,
                             "hop_effect": r.estimate, "ci_low": r.ci_low, "ci_high": r.ci_high,
                             "q_expr": r.q_expr, "p_label": r.p_label,
                             "tumours_supporting": f"{n_same}/{len(te)}", "n_tumours_same_sign": n_same,
                             "hop_support": "yes" if support else "no", "chain_label": chain_label,
                             "calibration_status": cstat,
                             "suggested_perturbation": f"reduce {lig} in {sender} cells (e.g. knockdown or {rec} blockade "
                                                       f"in co-culture), measure the {rprog} program score in adjacent "
                                                       f"{receiver} cells",
                             "caveat": "association in expression data; candidate from prior LR knowledge; not causal"})
    df = pd.DataFrame(rows)
    if len(df):
        sup = df.hop_support == "yes"
        df["_abs"] = df.hop_effect.abs()
        a = df[sup].sort_values(["_abs", "n_tumours_same_sign", "ligand_in_sender_program"], ascending=[False, False, False])
        b = df[~sup].sort_values(["chain", "hop"])
        a["rank"] = np.arange(1, len(a) + 1)
        b["rank"] = np.nan
        df = pd.concat([a, b]).drop(columns="_abs")
    path = cc.RESULTS / "ranked_hypotheses.csv"
    if path.exists():
        raise FileExistsError("ranked_hypotheses.csv exists; never overwrite")
    df.to_csv(path, index=False)
    cc.register(path, PHASE, "Ranked candidate hypotheses (chain hop x catalogue LR; candidate, not evidence)")
    # hop-level eligibility summary (also when no LR pair is eligible)
    hops = []
    for st, hop in (("beta_BA", "B<-A"), ("beta_CB", "C<-B")):
        for r in tl[tl.stat == st].itertuples():
            hops.append({"chain": r.chain, "hop": hop, "prog_A": r.prog_A, "prog_B": r.prog_B, "prog_C": r.prog_C,
                         "estimate": r.estimate, "q_expr": r.q_expr, "p_label": r.p_label,
                         "hop_support": bool(r.q_expr < .05 and r.p_label < .05)})
    hops = pd.DataFrame(hops).drop_duplicates()
    cc.save_csv(hops, PHASE, "hop_support.csv", "Hop support per plan (q_expr < 0.05 and p_label < 0.05)", log)
    elig = []
    for t in cf.TYPES:
        for lig, rec, _ in catalogue():
            elig.append({"type": t, "gene": lig, "role": "ligand", "usable_for_type": lig in usable[t]})
            elig.append({"type": t, "gene": rec, "role": "receptor", "usable_for_type": rec in usable[t]})
    cc.save_csv(pd.DataFrame(elig).drop_duplicates(), PHASE, "lr_gene_usability_by_type.csv",
                "Catalogue genes usable (type-aware non-background) per cell type", log)
    summ = {"n_rows": len(df), "n_supported_rows": int((df.hop_support == "yes").sum()) if len(df) else 0,
            "n_hops_supported": int(hops.hop_support.sum()), "n_hops": len(hops)}
    cc.save_json(summ, PHASE, "phase5_summary.json", "Phase 5 key numbers", log)
    log.info("%s", summ)
    if len(df):
        log.info("top:\n%s", df.head(15)[["chain", "hop", "candidate_ligand", "candidate_receptor", "hop_effect",
                                          "hop_support", "tumours_supporting", "calibration_status"]].to_string(index=False))
    cc.finish_run(log, PHASE)


if __name__ == "__main__":
    main()

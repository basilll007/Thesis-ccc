"""v2 figures (PNG + PDF) from saved result tables only. Palette: reference categorical slots 1-3 (blue, orange,
aqua; validated all-pairs) with neutral text ink and recessive grid.

  python model_dev/src/v2_figures.py phase1|phase2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
FIG = C.RESULTS / "figures"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 9, "axes.titleweight": "bold", "text.color": INK})


def style(ax, grid_axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=170 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("saved", FIG / f"{name}.png")


def ci_points(ax, xs, rows, color, label, val="Q_caf_expr", lo="Q_ci_low", hi="Q_ci_high", off=0.0, annotate=True,
              fs=5.5, fmt="{:+.3f}", side=False):
    """Point + 95 % CI; with annotate=True the estimate and CI are written next to each point (above the CI, or to the
    right of the point when side=True, for crowded panels)."""
    v = np.array([r[val] for r in rows]); l_ = np.array([r[lo] for r in rows]); h = np.array([r[hi] for r in rows])
    x = np.array(xs, float) + off
    ax.errorbar(x, v, yerr=[v - l_, h - v], fmt="o", ms=5, color=color, lw=1.6, capsize=0, label=label)
    if annotate:
        for xi, vi, li, hi_ in zip(x, v, l_, h):
            txt = f"{fmt.format(vi)}\n[{fmt.format(li)}, {fmt.format(hi_)}]"
            if side:
                ax.annotate(txt, (xi, vi), xytext=(4, 0), textcoords="offset points", ha="left", va="center", fontsize=fs, color=INK2)
            else:
                ax.annotate(txt, (xi, hi_), xytext=(0, 2), textcoords="offset points", ha="center", va="bottom", fontsize=fs, color=INK2)


INDEX: list[dict] = []


def record(name: str, description: str, sources: list[str]) -> None:
    """Add a saved figure to the v2 figure index (written by write_index)."""
    import hashlib
    for ext in ("png", "pdf"):
        p = FIG / f"{name}.{ext}"
        INDEX.append({"figure": f"model_dev/results/figures/{name}.{ext}", "description": description,
                      "sources": "; ".join(sources), "bytes": p.stat().st_size,
                      "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})


def write_index() -> None:
    """Merge this run's figures into model_dev/results/figures/v2_figure_index.csv (one row per file)."""
    p = FIG / "v2_figure_index.csv"
    old = pd.read_csv(p) if p.exists() else pd.DataFrame(columns=["figure"])
    new = pd.DataFrame(INDEX)
    pd.concat([old[~old.figure.isin(new.figure)], new], ignore_index=True).sort_values("figure").to_csv(p, index=False)
    print("index", p)


def phase1():
    d = pd.read_csv(C.RESULTS / "v2_phase1" / "phase1_summary.csv")
    dec = json.loads((C.RESULTS / "v2_phase1" / "phase1_decision.json").read_text())
    b = pd.read_csv(C.RESULTS / "v2_phase1" / "b_contamination_by_adjacency.csv")
    f = pd.read_csv(C.RESULTS / "v2_phase1" / "f_target_genes_caf_vs_macrophage.csv")
    fig, axs = plt.subplots(2, 2, figsize=(13, 10))
    pf = lambda ok: "PASS" if ok else "FAIL"  # noqa: E731
    # A: Q by distance bin
    ax = axs[0, 0]
    bins = ["touching", "near", "separated"]
    for t, color, lab, off in (("myeloid:C1Q_TAM", BLUE, "C1Q_TAM (primary)", -0.28),
                               ("S_mean", AQUA, "spillover-free macrophage genes S", 0.0),
                               ("M_mean", ORANGE, "CAF marker genes M (pure spillover)", 0.28)):
        rows = [d[(d.group == "ref") & (d.stratum == "bin") & (d.level == lv) & (d.target == t) & (d.model == "xgb")].iloc[0] for lv in bins]
        ci_points(ax, range(3), rows, color, lab, off=off, side=True)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(3), ["touching\n(<= 12.1 um)", "near\n(12.1-24.3 um)", "separated\n(> 24.3 um)"])
    ax.set_ylabel("CAF-expression contribution Q (R2)")
    ax.set_title("A  CAF-expression contribution by distance to nearest CAF")
    ax.text(0.02, 0.03, f"check (a): C1Q_TAM CI > 0 beyond touching bin -> {pf(dec['check_a_pass'])}\n"
            f"stricter (separated bin only) -> {pf(dec['a_separated_bin_strict_pass'])}",
            transform=ax.transAxes, fontsize=7, color=INK, va="bottom", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GRID))
    ax.legend(frameon=False, fontsize=7, loc="upper right"); style(ax)
    # B: contamination by bin
    ax = axs[0, 1]
    bb = b[b.stratifier == "bin"]
    x = np.arange(3)
    ax.bar(x, bb.marker_cp10k, color=BLUE, width=0.55)
    ax.errorbar(x, bb.marker_cp10k, yerr=[bb.marker_cp10k - bb.ci_low, bb.ci_high - bb.marker_cp10k], fmt="none", ecolor=INK, lw=1)
    caf = float(b[b.stratifier == "CAF cells"].marker_cp10k.iloc[0])
    ax.axhline(caf, color=ORANGE, lw=1.6, ls="--")
    ax.text(2.35, caf, f"CAFs {caf:,.0f}", color=INK, va="bottom", ha="right", fontsize=7)
    for xi, v in zip(x, bb.marker_cp10k):
        ax.text(xi, v + 60, f"{v:,.0f}", ha="center", fontsize=7, color=INK)
    rate = pd.read_csv(C.RESULTS / "v2_phase1" / "b_spillover_rate.csv")
    rt = rate[rate.bin == "touching"].iloc[0]
    ax.set_xticks(x, ["touching", "near", "separated"])
    ax.set_ylabel("CP10k, sum over the 20 marker genes M")
    ax.set_title("B  CAF-marker transcripts in macrophages\n"
                 f"(touching: {rt.spillover_rate:.0%} [{rt.ci_low:.0%}, {rt.ci_high:.0%}] of the CAF excess)")
    style(ax)
    # C: checks
    ax = axs[1, 0]
    items = [("reference\nC1Q_TAM", ("ref", "myeloid:C1Q_TAM")), ("adjusted\nC1Q_adj+s_c", ("sc", "C1Q_TAM_adj")),
             ("S genes\nspill.-free", ("ref", "S_mean")), ("M genes\nCAF markers", ("ref", "M_mean")),
             ("CD45+ q75\nC1Q_TAM", ("cd45_gt_caf75", "myeloid:C1Q_TAM")), ("CD45+ q90\nC1Q_TAM", ("cd45_gt_caf90", "myeloid:C1Q_TAM"))]
    rows = [d[(d.group == g) & (d.stratum == "all") & (d.target == t) & (d.model == "xgb")].iloc[0] for _, (g, t) in items]
    colors = [BLUE, BLUE, AQUA, ORANGE, INK2, INK2]
    checks = ["reference signal", f"check (c): {pf(dec['check_c_pass'])}", f"check (d): {pf(dec['check_d_pass'])}",
              "spillover positive control", "identity check", "identity (strict)"]
    for i, (r, c) in enumerate(zip(rows, colors)):
        ci_points(ax, [i], [r], c, None, fs=6)
        ax.text(i, -0.0043, checks[i], ha="center", va="top", fontsize=6.5, color=INK,
                fontweight="bold" if checks[i].startswith("check") else "normal")
    ax.set_ylim(-0.0055, 0.0165)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(items)), [i[0] for i in items], fontsize=7)
    ax.set_ylabel("Q (R2), all receivers")
    ax.set_title(f"C  Checks (XGBoost). Verdict: {dec['verdict']}")
    style(ax)
    # D: target genes CAF vs isolated macrophages
    ax = axs[1, 1]
    f = f.sort_values("ratio_caf_over_iso")
    y = np.arange(len(f))
    ax.scatter(f.caf_mean_cp10k_median, y, color=ORANGE, s=30, label="CAFs", zorder=3)
    ax.scatter(f.iso_mac_mean_cp10k_median, y, color=BLUE, s=30, label="isolated macrophages", zorder=3)
    ax.scatter(f.touch_mac_mean_cp10k_median, y, color=AQUA, s=30, marker="D", label="touching macrophages", zorder=3)
    for yi, r in zip(y, f.itertuples()):
        ax.plot([r.caf_mean_cp10k_median, r.iso_mac_mean_cp10k_median], [yi, yi], color=GRID, lw=2, zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(y, [f"{g}{' *' if h else ''}  ({r:.2f})" for g, h, r in zip(f.gene, f.hypothesis_gene, f.ratio_caf_over_iso)], fontsize=7)
    ax.set_xlabel("mean CP10k (median over 13 tumours; log scale)")
    ax.set_title("D  Target genes in CAFs vs macrophages\n(label: ratio CAF / isolated macrophage; * hypothesis gene)")
    ax.legend(frameon=False, fontsize=7, loc="upper right"); style(ax, "x")
    fig.suptitle("v2 Phase 1: segmentation-spillover control for the CAF -> macrophage signal (development set, 13 tumours; tumour-bootstrap 95% CIs)", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig1_spillover")
    record("v2_fig1_spillover", "Phase 1 spillover control: Q by distance bin, contamination by adjacency, decision checks "
           "(PASS/FAIL), target genes CAF vs macrophage; estimates with tumour-bootstrap 95% CIs printed on the plot",
           ["model_dev/results/v2_phase1/phase1_summary.csv", "phase1_decision.json", "b_contamination_by_adjacency.csv",
            "b_spillover_rate.csv", "f_target_genes_caf_vs_macrophage.csv"])


def phase2():
    q = pd.read_csv(C.RESULTS / "v2_phase2" / "phase2_quantities.csv")
    show = [("myeloid_orig", "myeloid:C1Q_TAM"), ("myeloid_orig", "myeloid:SPP1_TAM"), ("tumour_orig", "tumour:classical"),
            ("tumour_orig", "tumour:basal_like"), ("myeloid_adj_sc", "C1Q_TAM_adj"), ("myeloid_adj_sc", "SPP1_TAM_adj")]
    show = [s for s in show if ((q.group == s[0]) & (q.target == s[1])).any()]
    qs = [("neighbour gain (all - own)", "neighbour gain", INK2), ("cross-type contribution (all - same)", "cross-type", BLUE),
          ("who is nearby: other-type composition (same_plus_other_comp - same)", "other-type composition", AQUA),
          ("what they express: other-type expression (all - same_plus_other_comp)", "other-type expression", ORANGE)]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axs[0]
    for j, (qq, lab, color) in enumerate(qs):
        rows = []
        for g, t in show:
            r = q[(q.group == g) & (q.target == t) & (q.model == "xgb (genes set)") & (q.quantity == qq)]
            rows.append(r.iloc[0] if len(r) else None)
        xs = [i for i, r in enumerate(rows) if r is not None]
        ci_points(ax, xs, [r for r in rows if r is not None], color, lab, "value", "ci_low", "ci_high", off=-0.24 + 0.16 * j)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(show)), [f"{t.split(':')[-1]}\n({g.split('_')[0]}{', adjusted' if 'adj' in g else ''})" for g, t in show], fontsize=7)
    ax.set_ylabel("held-out R2 difference")
    ax.set_title("A  XGBoost (genes set): neighbour gain and its other-type parts")
    ax.legend(frameon=False, fontsize=7); style(ax)
    ax = axs[1]
    rows = []
    for g, t in show:
        r = q[(q.group == g) & (q.target == t) & (q.model == "xgb (genes set)") & q.quantity.str.startswith("unique contribution")]
        rows.append(r)
    types = ["CAF", "tumour", "myeloid", "T", "other"]
    cols = {"CAF": ORANGE, "tumour": BLUE, "myeloid": BLUE, "T": AQUA, "other": INK2}
    for i, r in enumerate(rows):
        for k, t in enumerate(types):
            rr = r[r.quantity.str.contains(f" of {t} ")]
            if rr.empty:
                continue
            rr = rr.iloc[0]
            ax.errorbar(rr.value, i + (k - 2) * 0.15, xerr=[[rr.value - rr.ci_low], [rr.ci_high - rr.value]], fmt="o", ms=4,
                        color=cols[t], lw=1.3, label=t if i == 0 or (t == "myeloid" and i == 2) else None)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(range(len(show)), [f"{t.split(':')[-1]} ({g.split('_')[0]})" for g, t in show], fontsize=7)
    ax.set_xlabel("unique contribution (all - minus type), R2")
    ax.set_title("B  Unique contribution of each other cell type")
    h, l_ = ax.get_legend_handles_labels()
    ax.legend(h, l_, frameon=False, fontsize=7, title="neighbour type (other-type only)"); style(ax, "x")
    fig.suptitle("v2 Phase 2: cross-type share on the aggregated baselines (development set; paired tumour-bootstrap 95% CIs)", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig2_crosstype")
    record("v2_fig2_crosstype", "Phase 2 cross-type share on aggregated XGBoost: neighbour gain, cross-type, composition vs "
           "expression, unique per-type contributions; estimates with paired tumour-bootstrap 95% CIs printed on the plot",
           ["model_dev/results/v2_phase2/phase2_quantities.csv", "agreement_with_graph_estimate.csv"])


def phase3():
    """Survival: HR per SD with 95% CI, BH q (9-test family) and PH test p printed per estimate; CPTAC replication table."""
    import hashlib
    root = C.REPO_ROOT / "results" / "translational" / "phase3"
    run = sorted(p for p in root.iterdir() if (p / "power_min_detectable_hr.csv").exists())[-1]
    m = pd.read_csv(run / "cox_models.csv")
    ph = pd.read_csv(run / "PH_diagnostics.csv")
    v = pd.read_csv(run / "external_validation_cptac.csv")
    pw = pd.read_csv(run / "power_min_detectable_hr.csv")
    sigs = list(dict.fromkeys(m.signature))
    fig, axs = plt.subplots(1, 3, figsize=(19, 7), gridspec_kw={"width_ratios": [1.5, 1.5, 1.2]})
    for ax, fam, title in ((axs[0], "clinical_purity", "A  Adjusted: age, sex, stage, purity (primary)"),
                           (axs[1], "unadjusted", "B  Unadjusted")):
        mdhr = pw[(pw.cohort == "TCGA_confirmed_PDAC") & (pw.model == fam) & (pw.alpha == "0.05")].min_detectable_hr_per_sd.median()
        ax.axvspan(1 / mdhr, mdhr, color=GRID, alpha=0.7, lw=0)
        for k, sg in enumerate(sigs):
            for ds, color, off in (("TCGA", BLUE, -0.18), ("CPTAC", ORANGE, 0.18)):
                r = m[(m.cohort == f"{ds}_confirmed_PDAC") & (m.signature == sg) & (m.model == fam)]
                if r.empty or r.status.iloc[0] != "DONE":
                    continue
                r = r.iloc[0]
                y = k + off
                lab = f"{ds} strict (n={int(r.n)}, events={int(r.events)})" if k == 0 else None
                ax.errorbar(r.hr, y, xerr=[[r.hr - r.ci_low], [r.ci_high - r.hr]], fmt="o", ms=5, color=color, lw=1.5, label=lab)
                pt = ph[(ph.cohort == f"{ds}_confirmed_PDAC") & (ph.signature == sg) & (ph.model == fam) & (ph.term == "signature")]
                php = float(pt.p.iloc[0]) if len(pt) else np.nan
                ax.text(2.3, y, f"{r.hr:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}] q={r.q_bh:.3f} PH p={php:.2f}",
                        va="center", fontsize=6, color=INK)
        ax.axvline(1, color=INK2, lw=0.8, ls="--")
        ax.set_xscale("log"); ax.set_xlim(0.45, 2.25)
        ax.set_xticks([0.5, 0.75, 1, 1.5, 2], ["0.5", "0.75", "1", "1.5", "2"])
        ax.set_yticks(range(len(sigs)), [s.replace("_", " ") for s in sigs], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel(f"HR per SD (log). Grey band: |HR| below the 80%-power detectable HR in TCGA ({mdhr:.2f}).\n"
                      "Text: HR [95% CI], BH q over 9 signatures, PH test p", fontsize=7)
        ax.set_title(title); ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2); style(ax, "x")
    ax = axs[2]
    vv = v[v.model == "clinical_purity"].set_index("signature").reindex(sigs)
    ax.axis("off")
    rows = [[s.replace("_", " "), f"{r.TCGA_HR:.2f}", f"{r.TCGA_q_bh9:.3f}", f"{r.CPTAC_HR:.2f} [{r.CPTAC_ci_low:.2f}, {r.CPTAC_ci_high:.2f}]",
             r.same_direction, "yes" if r.holds_in_CPTAC_A3 else "no"] for s, r in vv.iterrows()]
    t = ax.table(cellText=rows, colLabels=["signature", "TCGA HR", "TCGA q", "CPTAC HR [95% CI]", "same dir.", "holds*"],
                 loc="center", cellLoc="center", colLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(6.5); t.scale(1, 1.5)
    t.auto_set_column_width(list(range(6)))
    for (i, j), c in t.get_celld().items():
        c.set_edgecolor(GRID)
        if i == 0:
            c.set_text_props(fontweight="bold")
    n_pass = int(((m.cohort == "TCGA_confirmed_PDAC") & (m.model == "clinical_purity") & (m.q_bh < 0.05)).sum())
    ax.set_title(f"C  CPTAC replication (adjusted models)\nPrimary rule (adjusted TCGA q < 0.05): {n_pass} of {len(sigs)} pass",
                 fontsize=9, fontweight="bold")
    ax.text(0.5, 0.04, "* holds = same direction AND CPTAC CI excludes 1 (amendment A3).\n"
            "Purity: TCGA ABSOLUTE vs CPTAC pathology cellularity.", transform=ax.transAxes, ha="center", fontsize=6.5, color=INK2)
    fig.suptitle("v2 Phase 3: bulk survival associations (Cox OS, HR per SD; strict-histology TCGA and CPTAC). "
                 "Prognostic expression associations, not causal.", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig3_survival")
    rel = run.relative_to(C.REPO_ROOT).as_posix()
    record("v2_fig3_survival", "Phase 3 survival: Cox HR per SD with 95% CI, BH q (9-signature family), PH test p, power band, "
           "CPTAC replication table", [f"{rel}/{f}" for f in ["cox_models.csv", "PH_diagnostics.csv",
                                                               "external_validation_cptac.csv", "power_min_detectable_hr.csv"]])
    for ext in ("png", "pdf"):
        km = run / f"signature_survival_v2.{ext}"
        if km.exists():
            INDEX.append({"figure": km.relative_to(C.REPO_ROOT).as_posix(), "description": "Phase 3 composite: median-split "
                          "KM curves (descriptive), adjusted forest plot, correlation matrix", "sources": rel,
                          "bytes": km.stat().st_size, "sha256": hashlib.sha256(km.read_bytes()).hexdigest()})


def tests():
    """Unit-test results (pytest JUnit XML) and the pre-declared decision checks, as one figure."""
    import xml.etree.ElementTree as ET
    root = ET.parse(C.RESULTS / "v2_tests" / "pytest_junit.xml").getroot()
    cases = []
    for tc in root.iter("testcase"):
        bad = tc.find("failure") is not None or tc.find("error") is not None
        st = "FAIL" if bad else ("SKIP" if tc.find("skipped") is not None else "PASS")
        cases.append([tc.get("classname").split(".")[-1], tc.get("name"), st, f"{float(tc.get('time', 0)):.2f}"])
    dec = json.loads((C.RESULTS / "v2_phase1" / "phase1_decision.json").read_text())
    checks = [["Phase 1", "(a) Q C1Q_TAM CI > 0 beyond touching bin", dec["check_a_pass"]],
              ["Phase 1", "(c) Q on spillover-adjusted target CI > 0", dec["check_c_pass"]],
              ["Phase 1", "(d) Q on spillover-free genes S CI > 0", dec["check_d_pass"]],
              ["Phase 1", f"verdict: {dec['verdict']}", None]]
    rc_p = C.RESULTS / "v2_phase2" / "reproduction_check_vs_v1.csv"
    if rc_p.exists():
        rc = pd.read_csv(rc_p)
        checks.append(["Phase 2", f"reproduces v1 R2 ({int(rc['pass'].sum())}/{len(rc)} within tolerance)", bool(rc["pass"].all())])
    sv = sorted((C.REPO_ROOT / "results" / "translational" / "phase3").glob("*/summary.json"))
    if sv:
        s = json.loads(sv[-1].read_text())
        checks.append(["Phase 3", "PH diagnostic synthetic validation", s["PH_validation"]["status"] == "PASS"])
        checks.append(["Phase 3", f"any signature adjusted TCGA q < 0.05 ({len(s['primary_adjusted_q_lt_0.05'])} of 9)",
                       bool(s["primary_adjusted_q_lt_0.05"])])
    checks = [[a, b, "PASS" if r is True else ("FAIL" if r is False else "-")] for a, b, r in checks]
    checks = [[a, b, ("YES" if r == "PASS" else "NO") if b.startswith("any signature") else r] for a, b, r in checks]
    n = {k: sum(c[2] == k for c in cases) for k in ("PASS", "FAIL", "SKIP")}
    fig, axs = plt.subplots(1, 2, figsize=(16, 0.3 * max(len(cases), 12) + 1.6), gridspec_kw={"width_ratios": [1.5, 1]})
    for ax, data, cols, title in ((axs[0], cases, ["module", "test", "result", "seconds"],
                                   f"A  Unit tests (pytest, all model_dev tests): {n['PASS']} passed, {n['FAIL']} failed, {n['SKIP']} skipped"),
                                  (axs[1], checks, ["phase", "pre-declared check", "result"], "B  Pre-declared decision checks")):
        ax.axis("off")
        t = ax.table(cellText=data, colLabels=cols, loc="upper center", cellLoc="left", colLoc="left")
        t.auto_set_font_size(False); t.set_fontsize(7); t.scale(1, 1.3)
        t.auto_set_column_width(list(range(len(cols))))
        for (i, j), c in t.get_celld().items():
            c.set_edgecolor(GRID)
            txt = c.get_text().get_text()
            if i == 0:
                c.set_text_props(fontweight="bold")
            elif txt in ("PASS", "FAIL"):
                c.set_text_props(fontweight="bold", color="#008300" if txt == "PASS" else "#e34948")
        ax.set_title(title, loc="left")
    ts = root.find("testsuite")
    fig.suptitle(f"v2 unit tests and decision checks (pytest run {ts.get('timestamp') if ts is not None else ''})", fontsize=10)
    fig.tight_layout()
    save(fig, "v2_fig0_tests_and_checks")
    record("v2_fig0_tests_and_checks", "Unit-test results (all model_dev tests) and every pre-declared decision check with its result",
           ["model_dev/results/v2_tests/pytest_junit.xml", "model_dev/results/v2_phase1/phase1_decision.json",
            "model_dev/results/v2_phase2/reproduction_check_vs_v1.csv", "results/translational/phase3/*/summary.json"])


if __name__ == "__main__":
    for step in sys.argv[1:]:
        {"phase1": phase1, "phase2": phase2, "phase3": phase3, "tests": tests}[step]()
    write_index()

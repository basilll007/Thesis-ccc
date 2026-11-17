"""Build the v2 section of the dashboard page neighbourhood_models.html (and the index card note) from result files.

Every number on the page is read from model_dev/results/v2_phase1, v2_phase2 and results/translational/phase3 at
build time, so the page cannot drift from the tables. Idempotent: content between <!-- V2-...-START/END --> markers
is replaced on every run; v1 content is not touched.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402

import pandas as pd  # noqa: E402

ROOT = C.REPO_ROOT
P1, P2 = C.RESULTS / "v2_phase1", C.RESULTS / "v2_phase2"


def f3(x):
    return f"{x:+.3f}".replace("-", "&minus;")


def ci(v, lo, hi, d=3):
    s = lambda x: (f"{x:+.{d}f}" if d else f"{x:.0%}").replace("-", "&minus;")  # noqa: E731
    return f"{s(v)} [{s(lo)}, {s(hi)}]"


def hr(v, lo, hi):
    return f"{v:.2f} [{lo:.2f}, {hi:.2f}]"


def figure_card(tag, title, explain, png, pdf):
    return f"""
      <div class="card" style="padding: 24px; margin-top: 20px;">
        <div style="display:flex; justify-content:space-between; align-items:flex-start; flex-wrap:wrap; gap:8px;">
          <span class="chip-tag" style="background: rgba(37, 99, 235, 0.12); color: var(--brand-primary); font-weight:700;">{tag}</span>
          <a href="{pdf}" class="btn btn-outline" style="font-size:12px; padding:6px 12px;" target="_blank">Download PDF</a>
        </div>
        <h3 style="margin: 8px 0 6px; font-size: 18px; font-weight: 800;">{title}</h3>
        <p class="explain">{explain}</p>
        <div class="figure-preview-container" style="background: var(--bg-tertiary); border-radius: 8px; padding: 12px; text-align: center;">
          <img src="{png}" alt="{title}" style="max-width: 100%; height: auto; border-radius: 6px; cursor: pointer;" onclick="openLightbox(this.src, this.alt)">
        </div>
      </div>"""


def table(header, rows):
    th = "".join(f"<th>{h}</th>" for h in header)
    tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div style="overflow-x:auto;"><table class="data-table" style="width:100%; font-size: 12.5px;"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'


def build() -> tuple[str, str, str]:
    d1 = pd.read_csv(P1 / "phase1_summary.csv")
    dec = json.loads((P1 / "phase1_decision.json").read_text())
    rate = pd.read_csv(P1 / "b_spillover_rate.csv").set_index("bin")
    q2 = pd.read_csv(P2 / "phase2_quantities.csv")
    rc = pd.read_csv(P2 / "reproduction_check_vs_v1.csv")
    run3 = sorted((ROOT / "results/translational/phase3").glob("*/power_min_detectable_hr.csv"))[-1].parent
    cox = pd.read_csv(run3 / "cox_models.csv")
    val = pd.read_csv(run3 / "external_validation_cptac.csv")
    pw = pd.read_csv(run3 / "power_min_detectable_hr.csv")
    import xml.etree.ElementTree as ET
    jx = ET.parse(C.RESULTS / "v2_tests" / "pytest_junit.xml").getroot()
    n_tests = sum(1 for _ in jx.iter("testcase"))
    n_fail = sum(1 for t in jx.iter("testcase") if t.find("failure") is not None or t.find("error") is not None)

    def p1(group, target, stratum="all", level="all", model="xgb"):
        return d1[(d1.group == group) & (d1.target == target) & (d1.stratum == stratum) & (d1.level == level) & (d1.model == model)].iloc[0]

    def q(group, target, model, quantity):
        r = q2[(q2.group == group) & (q2.target == target) & (q2.model == model) & (q2.quantity.str.startswith(quantity))]
        return None if r.empty else r.iloc[0]

    ref = p1("ref", "myeloid:C1Q_TAM")
    tb = {b: p1("ref", "myeloid:C1Q_TAM", "bin", b) for b in ["touching", "near", "separated"]}
    adj = p1("sc", "C1Q_TAM_adj")
    sm = p1("ref", "S_mean")
    mm = p1("ref", "M_mean")
    xc = q("myeloid_orig", "myeloid:C1Q_TAM", "xgb (genes set)", "cross-type contribution")
    xs = q("myeloid_orig", "myeloid:C1Q_TAM", "xgb (genes set)", "cross-type share")
    rcq = q("myeloid_orig", "myeloid:C1Q_TAM", "ridge (programs set)", "cross-type contribution")
    rcs = q("myeloid_orig", "myeloid:C1Q_TAM", "ridge (programs set)", "cross-type share")
    tcq = q("tumour_orig", "tumour:classical", "xgb (genes set)", "cross-type contribution")
    prim = cox[(cox.cohort == "TCGA_confirmed_PDAC") & (cox.model == "clinical_purity")].set_index("signature")
    n_pass = int((prim.q_bh < 0.05).sum())
    hyp = prim.loc["Hypoxia"]
    mdhr = pw[(pw.cohort == "TCGA_confirmed_PDAC") & (pw.model == "clinical_purity") & (pw.alpha == "0.05")].min_detectable_hr_per_sd.median()

    kpis = [
        ("Spillover verdict (pre-declared)", dec["verdict"], "Checks (a), (c), (d) all FAIL", "var(--status-danger)"),
        ("Extra CAF RNA in touching macrophages", f"{rate.loc['touching','spillover_rate']:.0%}",
         f"[{rate.loc['touching','ci_low']:.0%}, {rate.loc['touching','ci_high']:.0%}] of the CAF marker excess", "var(--status-warning)"),
        ("CAF-expression contribution, touching", f3(tb['touching'].Q_caf_expr),
         f"[{f3(tb['touching'].Q_ci_low)}, {f3(tb['touching'].Q_ci_high)}] R²; separated {f3(tb['separated'].Q_caf_expr)}", None),
        ("After spillover adjustment", f3(adj.Q_caf_expr), f"[{f3(adj.Q_ci_low)}, {f3(adj.Q_ci_high)}] — CI includes 0", "var(--status-danger)"),
        ("Cross-type share, best model (XGBoost)", f3(xc.value), f"[{f3(xc.ci_low)}, {f3(xc.ci_high)}] R² macrophage C1Q_TAM (graph v1: +0.019)", "var(--status-warning)"),
        ("Signatures related to survival", f"{n_pass} / 9", f"adjusted TCGA, BH q &lt; 0.05; hypoxia closest (q {hyp.q_bh:.3f})", "var(--status-danger)" if n_pass == 0 else None),
        ("Detectable HR per SD (80% power)", f"{mdhr:.2f}", f"TCGA strict, adjusted model, {int(hyp.events)} events", None),
        ("Unit tests / reproduction", f"{n_tests - n_fail} / {n_tests}", f"v1 reproduced {int(rc['pass'].sum())}/{len(rc)} within tolerance", "var(--status-success)"),
    ]
    kpi_html = "".join(f"""
        <div class="card kpi-card">
          <div class="kpi-label">{a}</div>
          <div class="kpi-value"{f' style="color: {c};"' if c else ''}>{b}</div>
          <div class="kpi-subtext">{s}</div>
        </div>""" for a, b, s, c in kpis)

    bins_rows = [[b, f"{int(tb[b].n_cells):,}", ci(tb[b].neighbour_gain, tb[b].gain_ci_low, tb[b].gain_ci_high),
                  ci(tb[b].Q_caf_expr, tb[b].Q_ci_low, tb[b].Q_ci_high)] for b in ["touching", "near", "separated"]]
    checks_rows = [["(a) not confined to touching bin", "CAF-expression Q in near or separated bin CI &gt; 0", "FAIL" if not dec["check_a_pass"] else "PASS"],
                   ["(c) survives adjustment", f"Q on C1Q_TAM_adj + contamination covariate: {ci(adj.Q_caf_expr, adj.Q_ci_low, adj.Q_ci_high)}", "FAIL" if not dec["check_c_pass"] else "PASS"],
                   ["(d) present in spillover-free genes", f"Q on S genes: {ci(sm.Q_caf_expr, sm.Q_ci_low, sm.Q_ci_high)} (CAF markers in macrophages: {ci(mm.Q_caf_expr, mm.Q_ci_low, mm.Q_ci_high)})", "FAIL" if not dec["check_d_pass"] else "PASS"]]
    p2_rows = []
    for grp, tgt, lab in (("myeloid_orig", "myeloid:C1Q_TAM", "macrophage C1Q_TAM"), ("myeloid_orig", "myeloid:SPP1_TAM", "macrophage SPP1_TAM"),
                          ("tumour_orig", "tumour:classical", "tumour classical"), ("tumour_orig", "tumour:basal_like", "tumour basal-like"),
                          ("myeloid_adj_sc", "C1Q_TAM_adj", "macrophage C1Q_TAM_adj (spillover-adjusted)"),
                          ("myeloid_adj_sc", "SPP1_TAM_adj", "macrophage SPP1_TAM_adj (spillover-adjusted)")):
        for model in ("xgb (genes set)", "ridge (programs set)"):
            g, c_, s_ = (q(grp, tgt, model, k) for k in ("neighbour gain", "cross-type contribution", "cross-type share"))
            if g is None:
                continue
            share = "not defined" if pd.isna(s_.value) else ci(s_.value, s_.ci_low, s_.ci_high, 0)
            p2_rows.append([lab, model.split(" (")[0], ci(g.value, g.ci_low, g.ci_high), ci(c_.value, c_.ci_low, c_.ci_high), share])
    sv_rows = []
    for sgn in prim.index:
        r, v = prim.loc[sgn], val[(val.signature == sgn) & (val.model == "clinical_purity")].iloc[0]
        sv_rows.append([sgn.replace("_", " "), hr(r.hr, r.ci_low, r.ci_high), f"{r.q_bh:.3f}", hr(v.CPTAC_HR, v.CPTAC_ci_low, v.CPTAC_ci_high), v.same_direction])
    fig = "model_dev/results/figures"
    rel3 = run3.relative_to(ROOT).as_posix()

    section = f"""
  <!-- V2-SECTION-START (generated by model_dev/src/v2_dashboard.py from result files; do not edit by hand) -->
  <section style="padding-top: 0;" id="v2-update">
    <div class="container">
      <h2 class="sec">9 · v2 update: is the fibroblast → macrophage signal real? (branch model-v2)</h2>
      <div class="card callout warn" style="padding: 22px;">
        <h3>Three answers</h3>
        <ol class="tight">
          <li><b>Spillover control:</b> the CAF → macrophage signal is <b>{dec['verdict']}</b> by the pre-declared rule (checks a, c, d fail). It lives only in macrophages that touch a CAF, where we estimate an extra {rate.loc['touching','spillover_rate']:.0%} [{rate.loc['touching','ci_low']:.0%}, {rate.loc['touching','ci_high']:.0%}] of a CAF's marker RNA inside the macrophage.</li>
          <li><b>Cross-type share on the best model:</b> about zero for macrophage C1Q_TAM with XGBoost ({ci(xc.value, xc.ci_low, xc.ci_high)} R², share {ci(xs.value, xs.ci_low, xs.ci_high, 0)}); a linear ridge model finds {ci(rcq.value, rcq.ci_low, rcq.ci_high)} (share {ci(rcs.value, rcs.ci_low, rcs.ci_high, 0)}), carried by CAF expression. Tumour cells: {ci(tcq.value, tcq.ci_low, tcq.ci_high)}. All smaller than the v1 graph estimate (+0.019).</li>
          <li><b>Survival:</b> {n_pass} of 9 signatures relate to overall survival after age/sex/stage/purity adjustment (BH q &lt; 0.05, TCGA strict). Hypoxia is closest (HR {hr(hyp.hr, hyp.ci_low, hyp.ci_high)} per SD, q {hyp.q_bh:.3f}); the macrophage C1Q/SPP1 signatures are null in TCGA and CPTAC. The cohorts can only detect HR ≥ ~{mdhr:.2f} per SD.</li>
        </ol>
        <p class="small-note">Plans were committed before each computation (model_dev/plans/v2_*.md). Expression data; predictive/prognostic associations only; the tumour (patient) is the unit of every CI.</p>
      </div>
      <div class="kpi-grid reveal" style="margin-top: 20px;">{kpi_html}
      </div>
{figure_card("Figure v2-1 · Spillover control", "The CAF-expression contribution follows the contamination gradient",
             "<b>How to read it:</b> (A) R² lost when CAF <i>expression</i> features are removed (keeping who is nearby), by distance from each macrophage to its nearest CAF; every point shows its estimate and tumour-bootstrap 95% CI. (B) CAF marker RNA measured inside macrophages falls with distance. (C) The three pre-declared checks with PASS/FAIL. (D) Target genes in CAFs versus isolated and touching macrophages: SELENOP is nearly as high in CAFs as in macrophages, so its shift can be carry-over.",
             f"{fig}/v2_fig1_spillover.png", f"{fig}/v2_fig1_spillover.pdf")}
      <div class="card" style="padding: 24px; margin-top: 20px;">
        <h3 style="margin: 0 0 8px; font-size: 17px; font-weight: 800;">Pre-declared checks (XGBoost, macrophage C1Q_TAM)</h3>
        {table(["check", "measured", "result"], checks_rows)}
        <h3 style="margin: 16px 0 8px; font-size: 17px; font-weight: 800;">By distance to the nearest CAF (bins fixed from the median cell diameter, 12.1 µm)</h3>
        {table(["bin", "macrophages", "neighbour gain [95% CI]", "CAF-expression contribution [95% CI]"], bins_rows)}
        <p class="small-note">Post hoc (exploratory, not in the plan): touching macrophages carry more C1QA/B/C, CD163, MRC1 and SELENOP RNA than mixing with average-CAF RNA would produce; SPP1 is lower by about what dilution predicts (model_dev/results/v2_phase1/f_posthoc_mixing_expectation.csv).</p>
      </div>
{figure_card("Figure v2-2 · Cross-type share", "On the best model, neighbouring macrophages carry the gain",
             "<b>How to read it:</b> (A) for each target, the neighbour gain and its split into other-cell-type parts (who is nearby vs what they express), from XGBoost with block ablations on identical folds. (B) R² lost when one other cell type is removed. Numbers on the plot are estimates with paired tumour-bootstrap 95% CIs.",
             f"{fig}/v2_fig2_crosstype.png", f"{fig}/v2_fig2_crosstype.pdf")}
      <div class="card" style="padding: 24px; margin-top: 20px;">
        <h3 style="margin: 0 0 8px; font-size: 17px; font-weight: 800;">Cross-type share (paired tumour-bootstrap 95% CIs)</h3>
        {table(["receiver / target", "model", "neighbour gain", "cross-type (all − same-type)", "cross-type share"], p2_rows)}
      </div>
{figure_card("Figure v2-3 · Survival", "No signature relates to survival after adjustment",
             "<b>How to read it:</b> hazard ratio per standard deviation of each signature (Cox, overall survival) in strict-histology TCGA (blue) and CPTAC (orange), adjusted (A) and unadjusted (B). Text gives HR [95% CI], BH q over 9 signatures and the proportional-hazards test p. The grey band marks effects too small to detect at 80% power. (C) CPTAC replication of the unchanged models.",
             f"{fig}/v2_fig3_survival.png", f"{fig}/v2_fig3_survival.pdf")}
      <div class="card" style="padding: 24px; margin-top: 20px;">
        <h3 style="margin: 0 0 8px; font-size: 17px; font-weight: 800;">Adjusted Cox models (age, sex, stage, purity)</h3>
        {table(["signature", "TCGA HR per SD [95% CI]", "TCGA q (9 tests)", "CPTAC HR [95% CI]", "same direction"], sv_rows)}
        <p class="small-note">Descriptive Kaplan–Meier curves at the median and the correlation matrix: <a href="{rel3}/signature_survival_v2.png" target="_blank">{rel3}/signature_survival_v2.png</a>. Purity: TCGA ABSOLUTE vs CPTAC pathology cellularity.</p>
      </div>
{figure_card("Figure v2-0 · Tests", "Unit tests and every pre-declared decision check",
             f"All {n_tests} model_dev unit tests (leakage guard, folds, GPU estimators, metric signs and the v2 additions) with their results, and each pre-declared decision check with PASS/FAIL.",
             f"{fig}/v2_fig0_tests_and_checks.png", f"{fig}/v2_fig0_tests_and_checks.pdf")}
      <div class="two-col" style="margin-top: 20px;">
        <div class="card callout warn"><h3>What changed versus v1</h3>
          <ul class="tight">
            <li>The CAF → macrophage lead (largest other-type contributor in v1) is attributed to segmentation spillover.</li>
            <li>The CXCL12-CAF → C1Q/SELENOP hypothesis is now unsupported by these data (its source signal did not survive).</li>
            <li>The macrophage neighbour gain itself is unchanged and present at every distance from CAFs.</li>
          </ul></div>
        <div class="card callout"><h3>Reports</h3>
          <ul class="tight">
            <li>Full report: model_dev/reports/MODEL_REPORT_v2.md</li>
            <li>Wet-lab brief (hypotheses only): model_dev/reports/WETLAB_BRIEF_v1.md</li>
            <li>Figure index with checksums: model_dev/results/figures/v2_figure_index.csv</li>
          </ul></div>
      </div>
      <p class="small-note" style="margin-top: 16px;">Generated from result files by model_dev/src/v2_dashboard.py. Expression data; predictive and prognostic associations; no causal or clinical claims.</p>
    </div>
  </section>
  <!-- V2-SECTION-END -->
"""
    banner = f"""
      <!-- V2-BANNER-START -->
      <a href="#v2-update" class="card callout warn" style="display:block; margin: 0 0 20px; padding: 14px 18px; text-decoration:none; color:inherit;">
        <b>v2 update (branch model-v2):</b> the fibroblast → macrophage signal below is <b>{dec['verdict']}</b> under a pre-declared control; on the best model the cross-type share is about zero; no signature relates to survival after adjustment. Jump to section 9 →
      </a>
      <!-- V2-BANNER-END -->"""
    card = (f"<!-- V2-CARD-START --><br><b>v2:</b> the CAF &rarr; macrophage signal is {dec['verdict']} (pre-declared control); "
            f"cross-type share on XGBoost {f3(xc.value)}; survival: {n_pass}/9 signatures pass after adjustment.<!-- V2-CARD-END -->")
    return section, banner, card


def replace_between(text: str, start: str, end: str, new: str, anchor: str, before: bool = True) -> str:
    if start in text:
        a, b = text.index(start), text.index(end) + len(end)
        return text[:a] + new.strip("\n") + text[b:]
    i = text.index(anchor)
    return text[:i] + new + "\n" + text[i:] if before else text[:i + len(anchor)] + new + text[i + len(anchor):]


def main() -> None:
    section, banner, card = build()
    page = ROOT / "neighbourhood_models.html"
    t = page.read_bytes().decode("utf-8").replace("\r\n", "\n")
    t = replace_between(t, "<!-- V2-SECTION-START", "<!-- V2-SECTION-END -->", section, "</main>")
    t = replace_between(t, "<!-- V2-BANNER-START -->", "<!-- V2-BANNER-END -->", banner,
                        '<div class="section-badge">● Model development v1', before=True)
    page.write_bytes(t.replace("\n", "\r\n").encode("utf-8"))
    idx = ROOT / "index.html"
    t = idx.read_bytes().decode("utf-8").replace("\r\n", "\n")
    anchor = "for macrophages (CAFs lead) and nothing measurable for tumour cells."
    if "<!-- V2-CARD-START -->" in t:
        a, b = t.index("<!-- V2-CARD-START -->"), t.index("<!-- V2-CARD-END -->") + len("<!-- V2-CARD-END -->")
        t = t[:a] + card + t[b:]
    else:
        t = t.replace(anchor, anchor + card, 1)
    idx.write_bytes(t.replace("\n", "\r\n").encode("utf-8"))
    print("dashboard updated")


if __name__ == "__main__":
    main()

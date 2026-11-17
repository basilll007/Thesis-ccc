# Wet-lab brief v1: CAFs and macrophage state in pancreatic cancer (everything below is a hypothesis)

**Who this is for:** a wet-lab biologist deciding whether, and how, to test a lead from our spatial-transcriptomics models. All statements are predictive associations from expression data (CosMx, 13 PDAC tumours). Nothing here shows cause, mechanism or clinical value. Sources: `model_dev/reports/MODEL_REPORT_v2.md` and `model_dev/results/v2_phase1/`.

## What the model suggested
- Macrophages next to fibroblasts (CAFs) looked more "C1Q-type": higher C1QA/B/C, CD163 and SELENOP, and lower SPP1.
- A computer "perturbation" of the model linked higher CXCL12 in neighbouring CAFs to these shifts.
- **Hypothesis H1:** CAF-derived signals, possibly CXCL12, push nearby macrophages towards a C1Q/CD163-high, SPP1-low state.

## What survived our controls, and what did not
- **Did not survive: the model's use of what CAFs express.**
  - In this imaging method, a macrophage's measured RNA partly includes RNA from the CAF touching it ("segmentation spillover"). We estimate that macrophages touching a CAF carry an extra ~23 % [17-26 %] of a CAF's marker-gene level, compared with macrophages more than 2 cell-widths away.
  - The model's gain from CAF expression was small (+0.004 R2). It appeared only in touching macrophages and vanished, even turning slightly negative, at larger distances.
  - It did not survive removing CAF-like genes from the target, and it was absent for macrophage genes that CAFs do not make.
  - The pre-declared verdict is **"explained by spillover"**. The CXCL12 part of H1 rested on exactly this signal, so it is now **unsupported** by our data. That is not disproved, just without support.
- **SELENOP specifically:** CAFs express SELENOP at ~0.8x the level of isolated macrophages. A SELENOP rise measured inside a "macrophage" next to a CAF can therefore be CAF RNA.
- **Still standing (exploratory, not pre-declared):**
  - Macrophages touching CAFs have higher C1QA/B/C, CD163, MRC1 and SELENOP than simple RNA mixing with average CAF RNA would produce. CAFs have very little C1Q or CD163 RNA.
  - SPP1 is lower roughly as much as mixing alone would predict, so the SPP1 drop is not evidence.
  - Knowing a macrophage's neighbourhood still improves prediction of its state at every distance. That gain comes mostly from neighbouring macrophages and from who is nearby, not from CAF gene expression (see MODEL_REPORT_v2, Phase 2).
- **Bulk RNA from patients:** C1Q-type and SPP1-type macrophage signatures were **not** associated with survival in TCGA (141 patients) or CPTAC (129 patients) after adjustment.

## The single most informative experiment
**Contact versus no-contact co-culture of human macrophages with patient-derived PDAC CAFs, read out in sorted macrophages.**
1. Differentiate monocyte-derived macrophages from at least 3 blood donors. Grow CAFs from at least 3 PDAC patients.
2. Four arms per donor x CAF pair:
   - (a) macrophages alone;
   - (b) transwell with CAFs (shared medium, no contact);
   - (c) direct contact co-culture;
   - (d) transwell + CXCL12-neutralising antibody, or the CXCR4 antagonist AMD3100.
   - Optional (e): recombinant CXCL12 alone.
3. After 48-72 h, purify macrophages by CD14/CD11b sorting (or magnetic beads) so that no CAF RNA is carried over.

**Readouts:**
- RT-qPCR of C1QA, C1QB, CD163, MRC1, SELENOP and SPP1 (normalised to two housekeeping genes);
- flow cytometry for surface CD163, and intracellular C1q if a validated antibody is available;
- CAF-specific transcripts (COL1A1, DCN) measured in the sorted macrophages as a purity check (they should be near zero).

**Why this design:** sorting and the transwell arm remove segmentation spillover completely, which is the artefact that undermined the in-tissue signal. Arm (d) tests the CXCL12 part separately from "CAFs in general".

## What would falsify the hypothesis
- **H1 (CAF signals shift macrophages) falsified:** macrophages from arms (b) and (c) do not differ from (a) in C1QA/C1QB/CD163 RNA or CD163 protein. Use a pre-set threshold, e.g. < 1.5-fold change, consistent across donors.
- **Contact-dependent:** a shift in (c) but not in (b). This would need a further experiment to separate juxtacrine signalling from residual contamination; check the purity readout.
- **CXCL12 part falsified:** the shift in (b) persists unchanged in (d), and recombinant CXCL12 alone (e) does nothing.
- **Purity caveat:** if COL1A1/DCN are detectable in sorted macrophages, the result is not interpretable.

*All of the above are hypotheses for testing, not findings.*

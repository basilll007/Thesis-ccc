# Chain Phase 5 plan: name the candidate signals

Committed before Phase 3 tumour-level results were aggregated and before any Phase 5 computation. These are candidate signals from prior knowledge (the existing ligand-receptor catalogue, `src/pdac/lr_catalog.py`, 33 pairs). They give direction by known biology only and are labelled **"candidate", not evidence**.

## Hops considered
- Every forward hop of the primary and secondary chains, i.e. B <- A (sender A, receiver B) and C <- B (sender B, receiver C), for each program pair.
- **Hop support** (Phase 3, tumour level, primary annotation): BH q < 0.05 under the expression block null **and** add-one p < 0.05 under the label block null for that hop's beta. These are the same as decision-rule criteria 1 and 2 for that hop.
- Hops with support are ranked first. Hops without support are also listed, with status "hop not supported" and no ranking weight. The file therefore always shows what was examined. The hypotheses come from supported hops only.

## Ligand-receptor eligibility (non-background genes only; Phase 1 Amendment 1, type-aware)
- The **ligand** must be on the panel and not background-level for the *sender* type (results/chain/phase1b/gene_background_typeaware_flags.csv). It counts as "in the sender's program" if it is one of the sender program's genes, for example CXCL12 for CAF:CXCL12. Otherwise it is "expressed by the sender", meaning above the type-aware background.
- The **receptor** must be on the panel and not background-level for the *receiver* type.
- Sender and receiver types missing from the type-aware table (fewer than 50 cells in every section) are ineligible.

## Deliverable: results/chain/ranked_hypotheses.csv
Columns:
- chain, hop, sender, receiver, sender_program, receiver_program;
- candidate_ligand, candidate_receptor, ligand_in_sender_program;
- hop_effect (standardised beta, tumour level), ci_low, ci_high;
- tumours_supporting (number of tumours with the same sign as the tumour-level estimate / number testable);
- hop_support (yes / no), chain_label (from Phase 3), calibration_status (Phase 4: within / outside calibrated range, or extrapolated);
- suggested_perturbation, for example "reduce CXCL12 in CAF (knockdown or blocking in co-culture), measure the T:cytotoxic program score in adjacent T cells";
- caveat, always: "association in expression data; candidate from prior LR knowledge; not causal".

Ranking: supported hops first; within those, by |hop_effect| descending, then by tumours_supporting. A ligand in the sender's own program ranks before a ligand that is only expressed. Unsupported hops follow, unranked (rank = NA).

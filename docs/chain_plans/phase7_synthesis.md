# Chain Phase 7 plan: synthesis

Committed before the summary is written.

- **docs/CHAIN_SUMMARY_v1.md**:
  - Lead with: "Do neighbours predict cell state beyond the cell's own type, and is there support for any A -> B -> C chain within the calibrated range?"
  - Then, in order: results by phase; the calibration table (from results/chain/phase4/calibration_table.csv); ranked hypotheses (top 10 rows of results/chain/ranked_hypotheses_dedup.csv, with their "not supported" status shown if no hop is supported); the claims register; limitations; at most 5 next steps.
- **Claims register**: results/chain/phase7/claims_register.csv, produced by src/chain/chain_phase7_synthesis.py.
  - Columns: id, claim, status (supported / exploratory / not supported / not testable), evidence file, key number.
  - Every number is read from a generated file by the script (no hand-typed numbers).
- **Key numbers table**: results/chain/phase7/key_numbers.csv.
- **Manifest check**: every file under results/chain and figures/chain is listed in results/chain/MANIFEST.csv, and every manifest path exists. Output: results/chain/phase7/manifest_check.csv.
- **Status labels**:
  - "supported": meets its pre-declared decision rule, within the calibrated range where the rule requires that;
  - "exploratory": passes some criteria or rests on a post-hoc observation;
  - "not supported": the rule was tested and failed;
  - "not testable": the measurement or the floors prevented the test.
- **README**: a section pointing to the summary, plans and key files.
- **Wording**: expression, not mutation; associations and ranked hypotheses; no causal, mechanistic or clinical claims; the tumour is the unit of inference.

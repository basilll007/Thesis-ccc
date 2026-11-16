# Translational progress log

## 2026-11-15: isolated start

Created translational-v1 worktree from origin/thesis-v3 at 8021098. Original working directory and PDAC/chain outputs are read-only inputs. Plans for all five phases are committed before computation. Seed 42, max four threads, one Visium slide in memory. No source dataset properties assumed from accession names.

## Acquisition

Source downloads and access blocks: results/translational/acquisition/20261116T035450398108Z/downloads.csv. Native Git LFS media endpoints are used only when a successfully read raw pointer specifies the exact expected SHA-256; no denied endpoint is bypassed.

## Phase 1 DONE

NOT ROBUST. Failed checks: ['b_tumour_host', 'c_2000_UMI']. Source: results/translational/phase1/20261116T040752758116Z/summary.json. Baseline aggregate/test and every raw slide reproduced. Random signature tail p=0.00999000999000999.

## Phase 2 DONE

Two public cohorts acquired. Published TCGA curated supplement is BLOCKED by CAPTCHA; strict clinical histology fallback and all-primary sensitivity used. Cohort counts and follow-up: results/translational/phase2/20261116T043554395029Z/cohort_summary.csv.

## User-requested PR hand-off

The user requested stopping further analysis and creating a PR, then authorised proceeding without clarification. Phase 1 DONE (NOT ROBUST), Phase 2 DONE, Phase 3 NOT RUN (execution approval interrupted; no phase3 output directory), Phase 4 NOT RUN, Phase 5 partial draft hand-off only. Draft PR targets thesis-v3 from translational-v1. No outcome modelling, chain analysis, main-branch push or other-track file edits were performed after this instruction.

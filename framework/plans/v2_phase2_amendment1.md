# v2 Phase 2 plan, amendment 1 (committed before any Phase 2 deliverable result)

Date 2026-11-17. Amends `framework/plans/v2_phase2.md` section 2.2.

## What was seen

While writing unit tests on purely synthetic tissues from `ccc_audit.simulate`, I found a power problem. No real data
and no benchmark scenario were involved.

- Setup: 12,000 cells, true 2x contact effect, about 71 % of receivers touching an L+ sender.
- The **block-null Q test** (L and Rg permuted within 100 um tiles, statistic Q = E/M) missed the effect in 1 of 10
  seeds (p = 0.14).
- The ring-control test had p = 0.001 in all 10 seeds.
- Cause: Q only contrasts edges against the type mean. When most receivers are exposed, the within-tile permutation
  can barely move R between exposed and unexposed cells, so Q has little power.

## Change

The zonation check uses a **within-tile receiver contrast** instead of the block-null Q:
- Statistic: mean Rg_log(exposed receivers) - mean Rg_log(unexposed receivers), over all R cells.
- Exposed is defined as in the ring test.
- Null: Rg_log permuted among R cells within the same 100 um tile (random offset), 1000 permutations.

Zonation PASS = `p_block_contrast < 0.05` and `p_ring < 0.05`.

The block-null Q p-value is still computed and reported (`p_contact_block`) for continuity with the thesis. It no
longer enters the verdict.

Chain power = min(power_global, power_block_contrast, power_ring), with the same plant-and-shift method.

## What does not change

- The abundance check: global within-type permutation of Q. It still detects apparent contact enrichment the way a
  naive spatial score would, including enrichment created by zonation. That is what the zonation check then has to
  explain.
- The priority rule, every threshold, the spillover check, and the power floor (0.5 at 1.5x).
- Unit tests now use synthetic tissues of 8,000 cells (1,800 um square). The 3,000-cell tissues sat below the
  declared power floor, so the engine correctly returned NOT_TESTABLE for them.

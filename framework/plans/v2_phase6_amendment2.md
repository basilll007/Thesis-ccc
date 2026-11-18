# v2 Phase 6 plan, amendment 2 (committed before any Phase 6 audit)

Date 2026-11-18. **The user asked to use the GPU wherever possible.**

## What runs on the GPU

- The Phase 6 real-data audits and the Phase 5 secondary cellalibi v0.2.1 runs use `AuditConfig(device="cuda")`.
  The permutation nulls run on the RTX 5060 through PyTorch (`cellalibi/gpu.py`).
- The test statistics, null designs, thresholds and verdict logic are unchanged. Only the random-number stream
  differs from the CPU path.

## Validation before use

`framework/tests/test_gpu.py`:
- CPU/GPU agreement on planted, spillover and null synthetic tissues: identical statistics, p-values within 0.06,
  same verdict;
- null calibration on the GPU: 200 datasets x 1000 permutations, rejection count inside the binomial 95 % interval.

## What stays on the CPU

- The Phase 4 / 4b benchmarks: Phase 4b was already running on CPU, and switching mid-run would mix RNG streams.
- The external tools: squidpy, SpatialDM, COMMOT and LIANA+ have CPU-only default code paths.

## Worker count

Phase 6 uses at most 3 worker processes. Each CUDA process adds about 0.5 GB of host RAM for the torch runtime, and
free RAM is close to the 8 GB floor. The GPU is shared by these processes.

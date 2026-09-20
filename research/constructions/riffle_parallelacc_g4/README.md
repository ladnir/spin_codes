# Riffle ParallelAcc g=4

**Status:** `DIAGNOSTIC_EXPLORATION`

This construction retains the Riffle DP g=4 outer stage and uniform four-bit
packet permutation. It replaces the deterministic inner map with four
lane-parallel accumulators.

Goal 01 derives the exact constituent interface and evaluates reduced instances:

- construction: `CONSTRUCTION.md`;
- goal: `GOAL_01_ACCUMULATOR_STRESS_TEST.md`;
- report: `proof/GOAL_01_ACCUMULATOR_STRESS_TEST.md`;
- primary receipt: `receipts/goal01_reduced_exact.json`;
- independent audit: `receipts/goal01_audit.json`.

Goal 02 applies the interface to all 26 authenticated support-33 outer words:

- goal: `GOAL_02_SUPPORT33_MOMENT.md`;
- report: `proof/GOAL_02_SUPPORT33_MOMENT.md`;
- moment receipt: `receipts/goal02_support33_moment.json`;
- independent moment audit: `receipts/goal02_audit.json`;
- low-state lower-family probe: `receipts/goal02_low_state_probe.json`.

Riffle ParallelAcc g=4 is a proof diagnostic. It does not replace the active
Riffle PacketMul-WrapMul-2Lap g=4 candidate and has no full-size distance claim.

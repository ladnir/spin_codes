# Riffle DP-2Lap g=4

**Status:** `PAUSED_OPEN`

This candidate retains the first-lap terminal state, wraps to the first
logical node, and applies the same deterministic convolution to the
overwritten first-lap word. It samples no additional permutation or packet
randomizer between laps.

The exact inner map is

\[
G(x)=F(F(x))\mathbin\oplus J(L(x)),
\]

where \(F\) is the zero-state one-lap map, \(L(x)\) is the first-lap terminal
state, and \(J\) injects that state into the wrapped lap.

## Pause point

The autonomous second-lap certificate is rigorous. A nonzero terminal state
grows beyond the target distance in the relevant wrapped prefix. The open
problem is the probability that the structured first lap reaches terminal
state zero.

The attempted local-distance proof reached component dimension 33 and then
encountered support-case growth. A global block-moment inequality remains a
diagnostic conjecture, not a proof.

Primary records:

- `../../explorations/riffle_dp_2lap_g4_initial_diagnostic.md`
- `../../explorations/riffle_dp_2lap_g4_global_permutation_proof_route.md`
- `../../explorations/riffle_dp_2lap_g4_global_block_moment_route.md`

No PacketMul proof or benchmark should be recorded under this candidate.

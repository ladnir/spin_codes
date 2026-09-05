# Riffle PacketMul-WrapMul-2Lap g=4

**Status:** ACTIVE_EXPLORATION

This candidate was approved on 2026-08-21. It retains the packet-local
nonzero GF(16) multipliers and adds one public nonzero GF(2^64) multiplier at
the wrap between laps.

The construction is frozen in CONSTRUCTION.md. Goal 01 is complete. It proves
the exact conditional-uniformity interface and derives the bad-state list
budget. The proof is in proof/GOAL_01_WRAP_UNIFORMITY_PROOF.md.

Goal 02 closes every first-node stratum from node 7,060 onward and isolates a
response-code spectrum obstruction in the early region. Goal 03 proves that
the wrapped response is exactly a two-node turnoff shift of the old zero-state
problem. It also shows that, after the first genuine packet, the remaining
object is a 128-state lifted orbit rather than the old 64-state autonomous
orbit.

Goal 04 tests that lifted route. Exact short-window certificates and
counterexamples refute every universal gap-block argument through nine nodes.
The first viable replacement theorem is the open ten-node bound
`d(C_10) >= 59` for a binary `[640,128]` observation code.

No security, distance, or integrated performance result is currently
certified.

## Artifact locations

- Proof reports: proof/
- Machine-readable receipts: receipts/
- Candidate-specific source: ../../scripts/ with prefix
  riffle_packetmul_wrapmul_2lap_g4_

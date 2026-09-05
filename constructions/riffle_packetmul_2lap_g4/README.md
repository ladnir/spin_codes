# Riffle PacketMul-2Lap g=4

**Status:** `PAUSED_OPEN`

Exploration moved on 2026-08-21 to the distinct successor candidate
Riffle PacketMul-WrapMul-2Lap g=4. This folder retains the completed Goals 01
through 05 and the unresolved four-node distance interval.

This candidate was approved on 2026-08-21. It adds one independent nonzero
GF(16) multiplier to each input packet before the existing two-lap map.

The construction specification is in `CONSTRUCTION.md`. The initial proof and
refutation targets are in `PROOF_PLAN.md`.

The first medium-sized checkpoint, `GOAL_01_ZERO_SYMBOL_GATE.md`, is complete.
It reduces the authenticated support-33 terminal-zero row to the explicit
mixed-character lemma \(q_\chi\le5/8\). The exact checkpoint is in
`proof/GOAL_01_ZERO_SYMBOL_GATE_CHECKPOINT.md`.

The larger checkpoint `GOAL_02_GLOBAL_ORBIT_WEIGHT.md` is complete. An exact
four-state certificate proves the global six-nonzero-nibbles-per-node bound.
The proof is in `proof/GOAL_02_GLOBAL_ORBIT_WEIGHT_PROOF.md`.

Together, Goals 01 and 02 close the authenticated support-33 terminal-zero
row. They do not yet prove the complete construction.

The medium checkpoint `GOAL_03_TRANSFER_LEDGER.md` is complete. Its audited
partial ledger is in `proof/GOAL_03_TRANSFER_LEDGER_REPORT.md`. The first open
row is the support-33, nonzero-terminal, nonsuffix event.

The medium checkpoint `GOAL_04_BOUNDARY_PREFIX.md` is complete. Its exact
four-node certificate moves the sufficient zero-prefix cutoff from node 22,653
to node 20,976 and closes 1,677 additional first-node strata. The proof is in
`proof/GOAL_04_BOUNDARY_PREFIX_PROOF.md`.

The first open row is now the support-33, nonzero-terminal event with first
occupied node at most 20,975.

The medium checkpoint `GOAL_05_FOUR_NODE_DISTANCE.md` is complete with a
bounded solver obstruction. It authenticates a weight-42 word, but the exact
distance remains in the interval 36 through 42. The report is in
`proof/GOAL_05_FOUR_NODE_DISTANCE_REPORT.md`. Goal 05 does not move the Goal 04
cutoff or change the probability ledger.

The recommended next checkpoint is a finite boundary-tail lemma for first
occupied nodes 20,972 through 20,975, where a bad output would need tail weight
at most 17.

No security, distance, or performance claim is currently certified for this
candidate. Results for Riffle DP-2Lap g=4 transfer only after an explicit
invariance or reduction argument.

## Artifact locations

- Proof reports: `proof/`
- Machine-readable receipts: `receipts/`
- Candidate-specific source: `../../scripts/` with prefix
  `riffle_packetmul_2lap_g4_`

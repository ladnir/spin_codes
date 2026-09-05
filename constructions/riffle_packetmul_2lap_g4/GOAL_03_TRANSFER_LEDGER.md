# Goal 03: transfer audit and partial probability ledger

**Status:** Complete. The audited report is
`proof/GOAL_03_TRANSFER_LEDGER_REPORT.md`; the exact ledger and independent
audit are `receipts/goal03_transfer_ledger.json` and
`receipts/goal03_audit.json`.

## Objective

Determine which results for Riffle DP-2Lap g=4 remain valid after independent
nonzero packet multiplication. Insert the proved Goal 02 character cap into
the authenticated support-33 ledger. Identify the first open row that blocks a
construction-level distance claim.

## Transfer classes

Each inherited artifact receives one class.

- **Transferred:** its statement and probability are unchanged.
- **Transferred with substitution:** its deterministic statement is unchanged,
  but PacketMul supplies a new probability bound.
- **Superseded:** a new PacketMul certificate proves the required statement by
  a different route.
- **Invalidated:** the artifact fixes packet values or binary weights that
  PacketMul changes, and no stated reduction preserves its claim.
- **Open:** no available certificate controls the required PacketMul event.

An artifact is not transferred merely because the convolution map is
unchanged. The audit must identify every random variable and every invariant
used by the artifact.

## Work packages

1. Record the exact boundary between the unchanged outer stage, packet-local
   multiplication, packet permutation, and deterministic two-lap map.
2. Audit the autonomous three-node certificate and its independent verifier.
3. Audit the one-data outer population and final-10,119-node placement row.
4. Replace the paused candidate's fixed-value terminal-zero obligation with
   the Goal 01 and Goal 02 PacketMul proof.
5. Reinterpret known fixed-value turnoffs under the multiplier probability
   space. Search for any claimed deterministic transfer that they refute.
6. Produce an exact partial ledger with no overlapping event charges.
7. State one explicit next lemma for the first uncovered stratum.

## Completion criteria

The goal is complete when:

1. every relevant inherited artifact has a transfer class and justification;
2. all transferred hashes, counts, and rational probabilities are audited;
3. the support-33 terminal-zero row is marked rigorous;
4. the known authenticated turnoff is correctly interpreted under PacketMul;
5. the partial ledger separates covered events from their complements; and
6. one finite next lemma or exact transfer failure is stated.

## Excluded work

This goal does not prove the complement ledger row. It does not prove the full
construction. It does not benchmark packet multiplication or the encoder.

## Outcome

The support-33 terminal-zero row and the support-35-through-39 suffix rows give
a disjoint partial upper bound below (2^{-40.972830672895}). The first open
row is the support-33 event with nonzero first-lap terminal state and first
occupied node before node 22,653. Goal 03 states that event as a finite
first-node-stratified tail lemma.

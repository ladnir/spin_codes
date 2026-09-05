# BlockAccumulateOnce

**Status:** `PAUSED_TRANSFER_THEOREM_RECORDED`

BlockAccumulateOnce applies one uniform global bit permutation and one binary
accumulator to a binary outer block code. The intended outer code combines
scaled BCH constituents with a scaled double-parity layer. Its local
parameters may grow with the final length.

This directory records a modular distance theorem. The theorem accepts any
upper bound on the outer binary weight spectrum and returns a distance bound
for the permuted accumulated code.

The record does not claim that a specific scaled BCH and double-parity family
meets the required spectrum condition. Establishing that condition is the
remaining outer-code problem.

Files:

- `CONSTRUCTION.md`: exact encoder and probability space;
- `TRANSFER_THEOREM.md`: exact finite-length theorem and proof;
- `manifest.json`: name, status, and boundary of the variant.
- `transfer_audit.json`: exhaustive verification through length 12.

The verification script is
`../../scripts/verify_block_accumulate_once_transfer.py`.

This variant is paused. The active proof effort remains
`Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4`.

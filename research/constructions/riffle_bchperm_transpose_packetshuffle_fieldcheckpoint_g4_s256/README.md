# Riffle BCHPerm-TransposePacketShuffle-FieldCheckpoint g=4 s=256

This candidate retains the four-bit packet permutation of the `g4` parent and
increases the FieldCheckpoint state from 64 to 256 bits.  A checkpoint epoch
still contains 256 input bits.  Hence, every state lane is visited once per
epoch instead of four times.

The change removes the aligned-lane obstruction found for the 64-bit state.
For a live state, the output of one epoch is uniform over all 256-bit vectors
except one vector determined by the input.  The packet positions therefore do
not enter the live-state transfer.

The complete floating-point certificate has 55.0716 bits of modeled margin at
9% relative distance.  It covers all group-rank histograms, all checkpoint
reset trajectories, and all exceptional all-one outer words.  The proof uses
an entrywise suffix envelope over the exact nonzero-packet-count matrices; it
does not use the false claim that a packed profile directly maximizes every
matrix moment.

The remaining numerical task is outward rounding.  The remaining coding task
is to prove or replace the modeled \([256,128,38]\) outer spectrum.

The first end-to-end implementation benchmark is substantially slower than
the 64-bit-state prototype.  On Peach, two isolated 21-sample runs give an
aggregate median of 17.48 ms for the 256-bit-state construction, compared
with 10.17 ms for the otherwise identical 64-bit-state control.  The current
GF(2^256) kernel processes two field elements in parallel and uses a
nine-product carry-less Karatsuba multiplication.  Thus, this receipt is a
serious baseline rather than a deliberately scalar implementation, but it is
not yet an optimized lower bound on the cost of the wider state.

See `CONSTRUCTION.md`, `PROOF_STATUS.md`, and
`proof/FLOATING_POINT_CERTIFICATE.md`.

The benchmark source and receipt are in `benchmark/` and
`receipts/packet4_s256_performance_peach_7950x.json`.

# Riffle BCHPerm-TransposePacketShuffle-LongAcc-FieldChecksum g=4 s=64

This folder explores the four-long-accumulator construction defined in
`CONSTRUCTION.md`.  The first calculation uses the same modeled
\([256,128,38]\) regular outer spectrum as the packet4 FieldCheckpoint work.

The current exact diagnostics cover every maximally packed profile and every
maximally split layered profile. Both leave at least 55.96 bits of
regular-word margin at relative distance 0.09. A coarse support-only
reduction fails, so it must not be interpreted as a construction failure.

The remaining proof issue is a two-extreme scalar comparison for arbitrary
mixed packet-rank profiles. Simpler epoch-level, entrywise, and packed-only
claims are false; `PROOF_STATUS.md` records the counterexamples and the
refined target.

The calculation uses an exact random field checksum. It does not yet replace
that checksum with a cheaper structured parity circuit.

See `PROOF_STATUS.md`, `receipts/regular_coarse_delta09.json`, and
`receipts/two_extreme_sparse_delta09.json`.

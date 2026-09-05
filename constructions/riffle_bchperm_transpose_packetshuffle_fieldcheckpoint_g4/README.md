# Riffle BCHPerm-TransposePacketShuffle-FieldCheckpoint g=4

This folder records a performance candidate derived from Riffle ExactPerm
FieldCheckpoint v1.  The candidate retains the BCH outer code and the
FieldCheckpoint inner code.  It replaces each uniform 8192-position region
permutation with a uniform permutation of 2048 fixed packets of width four.

Two isolated runs of the narrow packet-aware evaluator have medians of 10.11
ms and 10.08 ms on Peach.  Two surrounding production runs have medians of
11.12 ms and 11.06 ms.  Across all 42 samples per construction, the packet
median is 10.10 ms and the production median is 11.10 ms.  The packet
evaluator is approximately 9.0% faster and uses 5.5 MiB of hot setup schedule
instead of 12 MiB.

The complete transposed output matches a flat evaluator for the packet
construction.  This correctness test does not compare packet output with
exact-permutation output because the two constructions sample different
permutation distributions.

The benchmark justifies a proof attempt.  It does not establish that the
existing bit-transpose certificate transfers.  The earlier
TransposePacketShuffle analysis used different outer and inner constituents.

The first packet-specific proof audit is positive.  The one-active-block
transfer is exactly the bit-shuffle transfer and has 81.5283 bits of modeled
aggregate margin at 9% relative distance.  For two weight-38 outer words, an
exact packet-position calculation gives 163.9366 bits of modeled aggregate
margin.  The corresponding bit-shuffle calculation gives 165.7100 bits.
Thus, the first nontrivial sparse shell loses only 1.7735 bits.

The later medium-density audit is negative.  Activating all 2048 blocks in one
fixed packet lane defeats the 9% first-moment calculation by approximately
132,900 bits.  The favorable sparse results therefore do not justify a full
certificate attempt for the 64-bit-state candidate.  The distinct `g4 s256`
variant removes this aligned-lane mechanism by visiting every state lane once
per checkpoint epoch.

See `CONSTRUCTION.md` and
`PROOF_STATUS.md`.

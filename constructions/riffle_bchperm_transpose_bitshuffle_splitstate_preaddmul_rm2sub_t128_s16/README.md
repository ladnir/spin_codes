# Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul-RM2Sub t=128 s=16

**Status: fixed-inner regular proof closed.** This variant reduces the state
from 32 bits to 16 bits while retaining 128-bit epochs. It replaces the
parent's sparse nested maps with a fixed subcode of \(\mathrm{RM}(2,7)\).

The fixed map `A` has exact distance 48. The construction sets
`B` equal to the transpose of `A`. The selected code is self-orthogonal, so
`BA=0`. Exact MacWilliams transformation proves that the kernel of `B` has
distance four.

At 9% relative distance, the exact one-active margin is 80.1685 bits. The
aggregate margin for occupations 2 through 8192 is 125.5200 bits. Combining
those terms leaves 80.1685 bits of regular margin.

The smaller state does not yet have a wall-clock implementation. It halves
the width of field multiplication and the per-epoch scalar randomness. The
maps `A` and `B` also admit transposed circuits derived from the same
Reed--Muller evaluation circuit.

The modeled outer spectrum, large mixed all-one tail, and outward-rounded
calculation remain open.

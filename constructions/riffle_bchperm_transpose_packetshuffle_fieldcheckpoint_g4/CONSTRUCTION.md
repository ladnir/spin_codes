# Construction

The message has (2^{20}) bits.  Partition it into 8192 blocks of 128 bits.
Encode each block with `ExtendedBch256x128-Eq3`, which produces 256 bits.

Setup independently samples a permutation of the 256 coordinates for every
outer block.  Setup also partitions the outer blocks into the fixed groups

\[
G_j=(4j,4j+1,4j+2,4j+3),
\qquad 0\le j<2048.
\]

After the coordinate permutations, transpose the outer-block matrix into 256
regions.  In each region, the four bits from one group (G_j) form one
ordered packet.  Setup independently samples a uniform permutation of the
2048 packets in every region.  A packet permutation moves all four packet
lanes together and preserves their order.

Concatenate the permuted regions.  Apply the same 64-lane FieldCheckpoint
inner code used by Riffle ExactPerm FieldCheckpoint v1.  The inner code uses
256-element epochs and an independent nonzero multiplier in
\(\mathrm{GF}(2^{64})\) at every epoch boundary.

The packet law is different from a uniform permutation of all 8192 positions
in a region.  In particular, the membership of every four-block group remains
fixed across all regions.  Therefore, a distance proof must analyze this
packet law explicitly.

## Narrow transposed evaluator

One packet contains four 128-bit transposed elements and occupies 64 bytes.
The evaluator finalizes the four FieldCheckpoint values in reverse order and
writes them consecutively into one bucket stream.  All four destination outer
blocks belong to the same 2048-block tile.

One route record stores the local four-block group and four BCH coordinates.
The immutable hot schedule occupies 5.5 MiB:

- 1.5 MiB for one 24-bit bucket slot per input packet;
- 4 MiB for one 64-bit destination record per routed packet.

The reusable workspace occupies 40 MiB.  It contains a 32 MiB bucket buffer
and one 8 MiB outer tile.  The benchmark path performs no hot-path allocation.

The implementation is a performance prototype.  It remains separate from
the production exact-permutation class until a packet-distance argument is
available.

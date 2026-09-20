# Construction

The message has \(2^{20}\) bits.  Partition it into 8192 blocks of 128 bits.
Encode every block with `ExtendedBch256x128-Eq3`.

Setup samples an independent permutation of the 256 coordinates of every
outer block.  The fixed four-block groups are

\[
G_j=(4j,4j+1,4j+2,4j+3),
\qquad 0\le j<2048.
\]

After the coordinate permutations, transpose the outer-block matrix into 256
regions.  In each region, the four bits from \(G_j\) form one ordered packet.
Setup samples an independent uniform permutation of the 2048 packets in every
region.  The permutation preserves the four packet lanes.

Concatenate the regions and apply FieldCheckpoint with a 256-bit state.  One
epoch contains 256 input bits.  For position \(p\) in an epoch, compute

\[
q_p\gets q_p+x_p,
\qquad
y_p\gets q_p.
\]

At every epoch boundary, setup supplies an independent nonzero element of
\(\mathrm{GF}(2^{256})\).  The forward encoder multiplies the state by the
transpose of the corresponding multiplication map.  The transposed evaluator
uses the multiplication map itself.

The state size is part of the construction name.  This construction is
distinct from the 64-bit-state packet candidate.

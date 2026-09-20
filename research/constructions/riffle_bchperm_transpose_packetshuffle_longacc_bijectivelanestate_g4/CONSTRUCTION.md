# Construction

The working name is
`Riffle FullLaneFieldConv g=4 t=s=256`.  This is distinct from
`Riffle BCHPerm-TransposePacketShuffle-LongAcc-FieldChecksum g=4 s=64`.

The outer encoder and permutations are unchanged.  Encode 128-bit messages
with the modeled `[256,128,38]` block code, independently permute the 256
coordinates of each outer block, transpose blocks into four-block packets,
and randomly permute the packet positions in each region.

Partition a packet region into epochs of 64 packets.  Write an epoch input
as four lane words

\[
X=(X_0,X_1,X_2,X_3)\in(\mathbb F_2^{64})^4.
\]

The state is

\[
Q=(Q_0,Q_1,Q_2,Q_3)\in(\mathbb F_2^{64})^4.
\]

Let `A` be the invertible length-64 prefix-sum accumulator.  The epoch emits

\[
Y_\ell=A(X_\ell+Q_\ell),\qquad 0\leq\ell<4.
\]

Pack `Y` as one element of `GF(2^256)`.  Setup samples an independent
nonzero scalar `alpha_i` for each epoch and sets

\[
Q' = \operatorname{unpack}(\alpha_i\operatorname{pack}(Y)).
\]

The first state is zero.  For fixed setup randomness, the encoder is binary
linear.  Multiplication by a nonzero scalar is invertible.  Consequently, a
nonzero epoch output gives a nonzero next state.  For a fixed nonzero state,
a uniform nonzero scalar makes the next nonzero state uniform over
`GF(2^256)^*`.

The packet size remains `g=4`.  The convolution step and state size are both
256 bits: `t=s=256`.

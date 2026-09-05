# Construction

This candidate retains the outer code, coordinate permutations, and packet
permutations of the packet4 parent.  It replaces the parallel short
accumulators and the 256-bit checkpoint with four long accumulators and one
random field checksum.

Partition each packet region into epochs of 64 packets.  Write the input to
one epoch as four lane words

\[
X=(X_0,X_1,X_2,X_3)\in(\mathbb F_2^{64})^4.
\]

The incoming state is a word \(q\in\mathbb F_2^{64}\).  Let
\(A:\mathbb F_2^{64}\to\mathbb F_2^{64}\) denote the ordinary prefix-sum
accumulator.  The epoch emits

\[
Y_\ell:=A(X_\ell+q),
\qquad 0\leq \ell<4.
\]

Setup samples four independent coefficients
\(\gamma_{i,0},\ldots,\gamma_{i,3}\gets\mathrm{GF}(2^{64})\) for epoch
\(i\).  Interpret each \(Y_\ell\) as a field element and set

\[
q':=\sum_{\ell=0}^3\gamma_{i,\ell}Y_\ell.
\]

The four words \(Y_0,\ldots,Y_3\) are the 256 output bits of the epoch.  The
word \(q'\) is the input state of the next epoch.  The first epoch starts with
\(q=0\).

For fixed setup randomness, every operation is binary linear.  The state is
repeated across the four accumulator lanes before accumulation.  The field
checksum gives every nonzero 256-bit epoch output 64 random parity checks.

This construction is named
`Riffle BCHPerm-TransposePacketShuffle-LongAcc-FieldChecksum g=4 s=64`.
It is distinct from both FieldCheckpoint state widths.

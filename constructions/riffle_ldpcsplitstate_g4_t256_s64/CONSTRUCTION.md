# Construction interface

Fix the step width `t=256` and state width `s=64`.  Let

\[
A:\mathbb F_2^{64}\longrightarrow\mathbb F_2^{256},
\qquad
B:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{64}
\]

be binary linear maps that satisfy

\[
BA=0.
\]

One epoch receives an input `X` and an incoming state `Q`.  It computes

\[
Y:=X+A(Q),
\qquad
U:=Q+B(Y).
\]

Setup samples an independent
`alpha_i` from `GF(2^64)^*` for epoch `i`.  The epoch sets

\[
Q':=\alpha_i U
\]

and emits `Y`.  The first epoch starts with `Q=0`.

For fixed setup randomness, the epoch is binary linear.  The map from
`(X,Q)` to `(Y,U)` is invertible for every choice of `A` and `B`.  Indeed,

\[
Q=U+B(Y),
\qquad
X=Y+A(Q).
\]

The additional relation `BA=0` gives

\[
U=Q+B(X).
\]

Thus zero input preserves the incoming state before field multiplication.
If `Y=0`, then `U=Q`.  Output cancellation therefore cannot cancel a
nonzero state.

The construction name reserves the interface, not a particular pair
`(A,B)`.  Every modification of the constituent maps must receive a distinct
suffix.

## First fixed nested pair

Split the 256 output coordinates into `(q,a,p)` of widths `(64,128,64)`.
Let `P` map 64 bits to 128 bits.  Let `L` apply an accumulator, a fixed
permutation, and a second accumulator.  Let `S` and `R` map 64 and 128 bits,
respectively, to 64 bits.  Let `J` be a permutation of the final 64
coordinates.  Define

\[
A(q)=(q,L(Pq),J^{-1}(Sq+RL(Pq)))
\]

and

\[
B(q,a,p)=Sq+Ra+Jp.
\]

This gives `BA=0` by direct substitution.  Equivalently, `B=[S R J]` is the
parity-check map of a 192-dimensional code, and the additional 128 checks
that define `a=L(Pq)` cut out the 64-dimensional image of `A`.

The preferred fixed instance uses column weight 7 in `P` and column weight 3
in both `S` and `R`.  It reuses the packet-support-4 compressor selected for
the degree-6 baseline.  Its matrices are recorded in
`receipts/fixed_nested_pair_depth2_p7_fixedB_search.json`.  The degree-6
baseline remains in `receipts/fixed_nested_pair_depth2_search.json`.

One epoch is stored as 64 packets of four bits.  Packet lane zero contains
the 64 systematic coordinates.  Lanes one and two contain the 128 auxiliary
coordinates.  Lane three contains the final parity coordinates.  Setup
chooses `J` so that its column at packet slot `k` avoids the supports of the
three `S/R` columns at that slot.  Among such choices, setup selects a wiring
with no cancellation supported on at most three packets.  The wiring changes
no codeword weight and requires no XOR.

Setup rejects a sampled pair if one of the 256 output coordinates of `A` is
the zero linear form.  For every accepted pair, each coordinate is balanced
over the `2^64` state words.  The average weight over nonzero codewords is
therefore

\[
\frac{256\cdot2^{63}}{2^{64}-1}.
\]

Setup also rejects a pair if two output coordinates define the same linear
form or if three output-coordinate forms sum to zero.  An accepted pair has
3-wise independent output coordinates under a uniform state.  The preferred
fixed instance passes this audit and has one dependency of size four.

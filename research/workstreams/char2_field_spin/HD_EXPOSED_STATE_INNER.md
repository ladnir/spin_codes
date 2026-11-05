# Exposing the refreshed state in the packet output

Status: a separate construction proposal. The packet map and its inverse below
are algebraic. No implementation, packet-cap calculation, benchmark, or distance
certificate has been produced. The current encoder and its receipts are unchanged.

The existing inner keeps its refreshed state hidden until later emissions. A
proof can consequently admit long live-state paths whose marked emissions are
too small. This proposal changes which coordinates are transmitted: it replaces
some current emissions by the refreshed state already computed by the encoder.
The output length and state dimension stay unchanged.

## Packet map and inverse

Work over GF(256). A packet has input X in F^64 and incoming state Q in F^16.
Let A be the current 64-by-16 emission map, consisting of four RS[16,4] blocks.
Let C=A^T be the feedback map, so CA=0. The refresh M is the current invertible
diagonal--MDS--diagonal map, sampled independently for this packet in the ideal
ensemble.

Choose a fixed set S of sixteen output positions, with four positions in each
RS block. The sixteen-by-sixteen matrix C_S is invertible because each block
has independent sets of four columns. Let R be the other forty-eight positions.
Fix an ordering of S to identify its output with the sixteen state coordinates.
Define the new packet by

\[
W=Q+CX,\qquad Q'=MW,
\qquad \widetilde Y_R=X_R+A_RQ,
\qquad \widetilde Y_S=Q'.
\]

Initialize Q=0 and retain Q' as the state for the next packet. For fixed keys
and incoming Q, the inverse is

\[
\begin{aligned}
X_R&=\widetilde Y_R-A_RQ,\\
X_S&=C_S^{-1}
\left(M^{-1}\widetilde Y_S-Q-C_RX_R\right),\\
Q'&=\widetilde Y_S.
\end{aligned}
\]

Thus the complete causal inner remains a bijection. It emits sixty-four field
coordinates per packet and uses a sixteen-coordinate state. This preserves the
concatenated rate when the outer encoder is unchanged. The inverse is a proof
and correctness interface; it adds no encoder operation.

## The intended proof gain

Every nonzero outgoing state now contributes its actual weight to the current
packet's output mark. That mark is not deferred to a future emission or bounded
only through an endpoint factor.

For a fixed z in (0,1), put

\[
Z_z=(1+255z)^{16}-1,
\qquad
\nu_z(q)=\frac{z^{\operatorname{wt}(q)}}{Z_z}
\quad(q\ne0).
\]

The established MDS-refresh atom bound gives, for every fixed nonzero W and
every nonzero q,

\[
\Pr[MW=q]\,z^{\operatorname{wt}(q)}
\leq255^{-16}Z_z\,\nu_z(q).
\]

This suggests a zero/weighted-live comparison whose live reference is nu_z,
not uniform on nonzero states. The retained-output mark
`z^wt(X_R+A_RQ)` is independent of the fresh refresh randomness after fixing
Q and X. It can therefore multiply this pointwise bound. Weighted histories
still require a complete new packet-operator proof, especially the W=0 branch.

The reference must not be immediately replaced by its worst density relative
to the old uniform-live measure. Such a replacement can cancel the intended
gain. The new reference also does not justify reusing old regional entries or
old distance receipts.

One illustrative case is an empty packet with Q supported in one four-state
block. Puncturing four coordinates leaves an RS[12,4] emission with at least
nine nonzeros. The packet also emits the refreshed state's weight. Its minimum
can be only ten, below the old thirteen; the proposed benefit is probabilistic
spreading by M, not a universal increase in this minimum distance.

## An explicit adverse case

The punctured coordinates can hide an input cancellation. Fix one block and a
retained position j in that block. Set Q=0, choose X_j nonzero, and solve

\[
C_{S\text{ in block}}X_{S\text{ in block}}=C_jX_j.
\]

The four solved coordinates are all nonzero. Otherwise at most four columns
of that block's C would be dependent. This gives an input of weight five,
with CX=0 and Q'=0, whose new output has weight one. The old output had
weight five. Under independent uniform nonzero amplitudes on that fixed
five-position support, this cancellation has probability exactly 255^-4.

Consequently no entrywise improvement over the old packet operator is claimed.
A viable comparison must show that the stronger marked live-state behavior
outweighs these cheaper zero-return emissions in the relevant ensemble.

## Native dataflow and cost boundary

The inspected implementation is `spin/experiments/char2_field_spin/HdGroupedInner.h`.
Its `innerPacket4` calls four `innerCoset` instances, then `finishStriped`.
Each coset emits four ZMM vectors, striped across the four RS blocks. The
persistent refreshed state also occupies four ZMM vectors.

Choosing S to be one complete coset aligns the sixteen replacement coordinates
with those four stores. With the native-output layout, the proposal may store
the four refreshed-state vectors directly at those positions after refresh.
It need not add an output pass or a second state. This is a layout proposal,
not verified assembly equivalence.

The retained cosets still perform their existing A emission and C feedback.
The replaced coset needs only its C feedback; its A-emission butterflies and
XORs with X can be omitted. The same refresh remains. At the scalar arithmetic
level there is therefore no additional field multiplication, and there may
be fewer fixed emission products. An exact native instruction count depends
on the chosen coset and existing constant simplifications.

There is an important implementation constraint. Modes that feed Y into C
rely on the old identity CY=CX. That shortcut is invalid for the changed
output. Feedback must consume the original X in every coset, or use a proved
compensation. The input-feedback structure already exists in the header,
but changing its use can alter register pressure and instruction dependencies.

Other unresolved costs are the delayed four stores, physical state-coordinate
ordering, streaming-store order, and live-register requirements. The exact
diagonal gauge must also be rechecked now that state coordinates are visible.
A non-diagonal state-basis change is not automatically weight-preserving.
The proposal is compatible with the current operation budget in principle,
but no runtime near 3.4 ms is asserted.

## First falsifiable test

Before native work, construct a small-field packet test with explicit A, C,
and invertible refresh maps. Check the inverse and all weighted transition
branches against enumeration, including the five-coordinate cancellation.
Then compare the old and new packet operators under one fixed mark, retaining
the proposed nu_z reference instead of bounding it by uniform-live density.

The decisive question is whether the resulting all-support operator improves
the live excursions without making the zero-return branch dominate. A failure
there should stop this construction before another regional or per-arity
implementation. A success would justify a matched regional comparison; it
would still not establish whole-code coverage.

This alternative is independent of the current outer output-shuffle proposal.
The latter keeps the measured inner and may simplify the complete arity sum
without changing native state evolution. It remains the nearer proof test.

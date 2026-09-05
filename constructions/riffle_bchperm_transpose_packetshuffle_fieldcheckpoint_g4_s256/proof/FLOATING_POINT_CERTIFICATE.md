# Floating-point distance certificate

## Claim and scope

Fix a message length of (2^{20}) bits and an output length of
(N=2^{21}) bits.  Let

\[
d=\lfloor 0.09N\rfloor=188743.
\]

Under the modeled regular spectrum of a binary \([256,128,38]\) outer code,
the expected number of nonzero messages whose encoded word has weight at most
(d) is at most

\[
2^{-55.0716220969}.
\]

The expectation is over the coordinate permutations, packet permutations,
and nonzero field multipliers sampled during setup.  Thus the same number
upper-bounds the probability that a sampled encoder has such a bad message.
In particular, the modeled ensemble contains an encoder of minimum distance
at least (188744).

The calculation uses nearest binary64 arithmetic.  It is a floating-point
certificate, not an outward-rounded certificate.  It also assumes the
displayed outer spectrum; proving or replacing that spectrum is a separate
obligation.

## Outer reduction

Call an outer word regular when it is neither zero nor the unique all-one
word.  After the independent coordinate permutation of a regular outer
block, the modeled spectrum is bounded pointwise by twice the uniform
256-bit spectrum.  Therefore, one active regular block contributes the mass

\[
m=2(2^{128}-1)
\]

to the uniform proxy calculation.  If (a) of the 8192 outer blocks are
regular and active, their placement and word count are bounded by

\[
\binom{8192}{a}m^a.
\]

The exceptional all-one word is not included in this density reduction.  It
is handled after the regular calculation.

## Epoch relaxation

Consider one 256-bit checkpoint epoch.  The state is visited once, so its
output vector is also its final state before the nonzero field
multiplication.

First suppose that the incoming state is zero.  Replace every nonzero input
packet by one output bit and discard the packet's remaining Hamming weight.
If the epoch contains (h) nonzero packets, this replacement emits (h)
bits and leaves a nonzero state.  The replacement can only decrease output
weight.

Now suppose that the incoming state is nonzero.  The preceding random
nonzero field multiplier makes it uniform over the (2^{256}-1) nonzero
states.  For every fixed input vector, the epoch output is a translate of
that state.  Define

\[
\kappa=\frac{2^{256}}{2^{256}-1}.
\]

The exact live-state transition is bounded entrywise by \(\kappa\) times the
following normalized transition: sample a uniform 256-bit output, move to
the zero state if that output is zero, and otherwise remain live.  There are
8192 epochs, so replacing every live transition costs at most

\[
\kappa^{8192}.
\]

Its base-two logarithm is approximately
(1.02067\mathbin{\cdot}10^{-73}), and the checker includes it explicitly.

For a Chernoff variable (0<z<1), the normalized epoch matrices conditioned
on (h) nonzero packets are

\[
E_0(z)=
\begin{pmatrix}
1&0\\
2^{-256}&((1+z)^{256}-1)2^{-256}
\end{pmatrix}
\]

and, for (h>0),

\[
E_h(z)=
\begin{pmatrix}
0&z^h\\
2^{-256}&((1+z)^{256}-1)2^{-256}
\end{pmatrix}.
\]

These matrices retain every zero-to-live, live-to-zero, and repeated-reset
trajectory.

## Packet-count envelope

One region has 2048 packet positions, divided into 32 epochs of 64 packets.
Condition on exactly (H) nonzero packets in the region.  The packet
permutation makes their positions a uniform (H)-subset.  The checker uses a
hypergeometric matrix recurrence to compute the exact averaged region matrix
(R_H(z)) for every (0\le H\le2048).

The matrices (R_H) need not be monotone in (H).  The checker therefore
defines the entrywise suffix envelope

\[
\overline R_H(z):=\max_{K\ge H}R_K(z).
\]

This is the step that avoids a false global packing claim.

If a packet group contains (r) active regular blocks, its proxy packet is
nonzero with probability (1-2^{-r}).  A pairwise packing move preserves the
probability that two affected packets are both zero and decreases the
probability that both are nonzero.  Repeated moves show that the maximally
packed rank profile has the stochastically smallest nonzero-packet count.
Let (H_{\rm pack}) and (H_{\rm any}) be coupled counts with
(H_{\rm pack}\le H_{\rm any}).  Then, entrywise,

\[
R_{H_{\rm any}}(z)
\le \overline R_{H_{\rm pack}}(z).
\]

For (a=4q+r), the packed count is a binomial count from (q) rank-four
groups, plus one independent rank-(r) packet when (r>0).  The checker
averages \(\overline R_H\) over this distribution exactly in the log
semiring.  Independent regions then give the 256th matrix power.

For each active-block count, the checker selects one displayed Chernoff
parameter and applies

\[
\Pr[W\le d]\le z^{-d}\mathbb E[z^W].
\]

Summing all 8192 nonzero regular counts gives

\[
U_{\rm regular}\le2^{-55.0716221062}.
\]

The one-active-block row dominates this sum.

## Exceptional all-one words

Suppose (a) blocks contain regular words and at least one other block
contains the all-one word.  Their outer multiplicity is bounded by

\[
\binom{8192}{a}m^a(2^{8192-a}-1).
\]

The checker uses three valid input deletions.

For (0\le a\le105), it deletes every regular input and keeps one
deterministically nonzero packet per region.  For
(106\le a\le136), it keeps one forced packet and

\[
\left\lfloor\frac{a-3}{4}\right\rfloor
\]

packed rank-four regular packets.  The subtraction by three permits all
regular blocks sharing the retained forced packet to be deleted.  For
(a\ge137), it deletes every all-one input and reuses the regular-only bound.

Input deletion lowers the nonzero-packet count.  The suffix envelope makes
each comparison valid without assuming monotonicity of (R_H).  Summing all
mixed configurations gives

\[
U_{\rm all-one}\le2^{-82.2820153825}.
\]

Finally,

\[
U_{\rm regular}+U_{\rm all-one}
\le2^{-55.0716220969}.
\]

## Reproduction

Run

```powershell
python scripts/certify_riffle_packet4_fieldcheckpoint_s256_support.py
python scripts/certify_riffle_packet4_fieldcheckpoint_s256_allone.py
python scripts/verify_riffle_packet4_fieldcheckpoint_s256_certificate.py
```

The first command writes the regular receipt.  The second reads that receipt,
adds every exceptional all-one configuration, and writes the final receipt.
The third recomputes both row sums, checks that every active-count range is
covered, and checks the 40-bit target.

The proof artifacts are:

- `receipts/all_profiles_regular_support_relaxation_delta09.json`;
- `receipts/all_one_completion_delta09.json`.

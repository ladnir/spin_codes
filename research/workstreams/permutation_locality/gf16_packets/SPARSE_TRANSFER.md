# Transferring the Unweighted Sparse Argument

The old sparse proof contains a shape-uniform case that applies to the
GF(16) construction. Its inner operators must use all-one penalty one
and input-weight tilt one. The GF16 entry point fixes both parameters;
it does not import old numerical certificates using other values.
The last section gives a separate, normalized change of measure that
allows shape penalties without changing the outer counting measure.

Fix a nonzero message before sampling setup. Its four BCH words in group
g have independently shuffled row supports. Let U_g be the number of
regions where at least one of those rows is nonzero. GF multiplication
preserves packet activity, so it leaves U_g and the group's active-region
set unchanged.

Conditional on any positive-probability event specifying the U_g values,
each group's active-region set is a uniform U_g-subset of the 256
regions. Different groups remain independent. This follows from
exchangeability of the independent row shuffles, as in the original
independent-row bridge. The same averaged outer union-support CDF applies.

Conditional on the row supports, the GF randomizers make the active
packet values independent and uniform nonzero. Equivalently, for each
active (group,region) pair, first sample a packet weight from
(4,6,4,1)/15, then sample a uniform binary packet of that weight. This
sampling can be performed before the independent regional packet
permutations, because the original independent multipliers give this
same conditional distribution at every routed position.

Now condition on the complete sampled weight assignment attached to
those group-region pairs. The remaining packet values are independent
uniform masks in their prescribed weight shells. The regional
permutations remain uniform and independent. This is precisely a
distribution covered by the old shape-uniform inner operators: each
operator bounds every assignment of nonzero packet weights, before
composing across steps and averaging regional placement.

With both penalties equal to one, those bounds contain no reciprocal
factor depending on the original or transformed packet weights. Average
over the GF weight assignment for each fixed active-region pattern to
obtain the same conditional inner upper bound for the new construction.
Then average over the active-region patterns conditional on the U_g.
We do not assume that these patterns remain uniform after conditioning
on a weight assignment that identifies active positions. The bound is
uniform over weights before the active-pattern average is taken.
The outer union-support count,
the support-cover argument, and the factor choosing the active group
locations are unchanged.

This transfer does not say that the two codes have equal output-weight
distributions. It transfers a universal upper bound. In particular, a
bound involving the original number of all-one packets would not
transfer by this argument, because GF mixing changes that number.

`sparse_cover.py` therefore regenerates the penalty-one operators and
complete support covers at the requested cutoff. It checks every
requested occupancy; no interpolation in occupancy is used.
`assemble.py` regenerates these covers without cached operators and
replays the dense cover at the same cutoff before summing both ranges.
The sparse search checks its stopping bound after each subdivision; its
outward evaluation and exact support-volume coverage check are unchanged.

## Paying for a Shape Tilt Under the GF Distribution

The direct bridge protects against every nonzero packet shape, including
all ones at every active position. Under independent GF randomization,
an active packet equals 1111 with probability 1/15. We can exploit this
frequency without rebuilding the state envelope.

Fix the message and its active group-region pairs. Let P be the uniform
distribution on the 15 nonzero four-bit packets. For positive a and
0 < rho <= 1, define

    W(x) = a^wt(x) rho^[x=1111],
    c = E_P[1/W(X)]
      = sum_{w=1}^4 binom(4,w) / (15 a^w rho^[w=4]),
    Q(x) = P(x) / (c W(x)).

Q is a probability distribution and remains uniform within each weight
shell. For j independent packets and every nonnegative function H,

    E_{P^j}[H(X)] = c^j E_{Q^j}[H(X) product_i W(X_i)].

Regional routing and inner setup retain their original distributions.
The identity also holds after averaging those independent random choices.
The weighted shape-uniform operators bound the right-hand weighted
expectation for every weight assignment. Their fresh-state coordinate
allows convex mixtures of the same shell distributions, so the changed
shell probabilities do not require a new state invariant.

A region with j active groups contains exactly j active packets. Every
placement term therefore has the same compensation c^j, irrespective of
how those packets occupy inner steps. Multiply its weighted region
operator by c^j to obtain a valid unweighted operator for the GF ensemble.
Average the active-region patterns only after applying this bound, as in
the direct bridge above. The outer CDF and its group-location factor are
unchanged. No factor involving the original BCH packet weights is used.

For a=1, c=(14+1/rho)/15. Thus rho=1/2 costs 16/15 per active packet,
not a factor two per packet. The benefit must exceed that cost in the
complete output moment; this identity does not itself close a distance
certificate.

`shape_tilt.py` applies the compensation to exact outward region
operators. `--inner gf-shape-tilt` selects this method in the sparse
driver; `--sparse-inner gf-shape-tilt` selects it in the full assembler.
The normalizer is rounded upward before multiplication. Exact rational
tests check the identity for zero, one, and two packets, several positive
tilts, and every nonzero packet value. A full claim still requires
complete support coverage and the final sparse/dense sum.

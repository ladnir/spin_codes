# Retaining BCH Row Parity in the GF16 Comparison

The earlier dense comparison replaces a shuffled BCH row by a positive
mixture of independent bits. It discards the fact that every row has
even weight. After GF16 packet randomization, part of that constraint
can be retained without adding a state coordinate to the inner proof.

This changes a counting bound, not the encoder or its setup distribution.
It applies to the independently shuffled four-row groups described in
[the workstream README](README.md).

## One Mixture Component

Fix a component with independent lane bits having probabilities
p_1,...,p_4 in each of 256 packet positions. Let

    e = product_i (1-p_i),     a = 1-e.

Thus a packet is active with probability a. If a=0, every packet is zero
and all row parities are even. Otherwise fix an activity support S of
size u. Before field multiplication, the packets at positions in S are
independent copies of V conditioned on V!=0; the other packets are zero.

For a binary four-vector b, define its conditional character

    phi_b = (product_{i:b_i=1}(1-2p_i) - e)/a.

This equals E[(-1)^(b dot V) | V!=0], so |phi_b|<=1. The four rows all
have even weight exactly when the XOR of the u active packets is zero.
Character orthogonality therefore gives

    Pr[all four rows even | S] = (1/16) sum_b phi_b^u.

Suppose u>=d>=1. The following rational constant is a uniform upper
bound, independent of the positions in S:

    rho_d(p_1,...,p_4) = (1/16) sum_b |phi_b|^d.

Each independently sampled nonzero GF16 multiplier maps a fixed nonzero
packet uniformly onto the 15 nonzero values. Consequently, conditional
on S, the randomized values are independent uniform nonzero packets,
even after restricting the original rows to even weight. For each fixed
randomized output with support S, the original component's parity-filtered
mass is thus at most

    rho_d * a^u (1-a)^(256-u) / 15^u.

The right side is a positive iid activity comparison with its coefficient
multiplied by rho_d. Although the implementation samples multipliers
after regional routing, this calculation may precede routing: a fixed
packet permutation merely reindexes the independent uniform multipliers.
Pushing the resulting measure domination through that routing preserves it.

## Why the Support Floor Is Available

The authenticated BCH spectrum caps exclude every odd weight and every
nonzero weight below 38. A nonzero four-row message contains at least
one nonzero row, so its union support has size at least d=38. Restricting
the existing row-mixture envelope to even weights therefore preserves
its pointwise domination. Expand the four-row product, apply the bound
above to each component, and sum its positive coefficients.

The all-zero group remains the separate inactive component with mass one.
Only components bounding nonzero groups receive the factor rho_38.
Their iid comparisons may assign positive mass to supports below 38;
that extra mass cannot hurt domination because the actual nonzero-group
measure vanishes there. The active-group labels retain their previous
meaning in the outer sum.

The gain is not always a factor of 16. With one nondeterministic row and
three zero rows, activity is that row's bit pattern, so rho_38=1. With
four unbiased rows,

    rho_38 = (1 + 15*(1/15)^38)/16.

The implementation keeps the exact rational factor for every component;
it does not substitute 1/16 or assume independent row parities after
conditioning on activity.

## Implementation and Current Evidence

`parity_mixture.py` reconstructs the authenticated caps and checks the
parity premise. `scalar_cover.py --row-parity` applies the factors before
forming the existing activity mixture. Its witnesses record this option;
resume, replay, and the full assembler reconstruct the same comparison.
Older witnesses default to the original comparison.

Exhaustive tests enumerate all packet inputs through three positions,
including biased, zero, and deterministic-one lanes. They compare the
even-row measure with the bound for every activity support. Additional
tests check the one-row exception and preservation of the inactive term.

At 8.5% and a dense handoff of q=129, floating diagnostics improve the
proposed log2 bounds as follows. These are comparison coordinates, not
message weights, and selected points are not a full cover.

| Comparison coordinate | Earlier bound | With row parity |
|---|---:|---:|
| 0.03 | +113.09 | -259.19 |
| 0.04 | +1337.92 | +824.36 |
| 0.05 | +2155.81 | +1525.18 |

The parity refinement helps but does not close 8.5% with the current
return bound. At 8.2%, the full dense search for q=105--2048 closes with
173 accepted cells and no unresolved cells. Its full replay is running.
A whole-code claim still requires fresh outward replay and a matching
sparse prefix. The completed 8% dense cover does not use this refinement.

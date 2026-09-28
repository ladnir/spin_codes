# Retaining the profile-dependent shuffle loss

This scalar inequality is implemented as an optional affine-moment
witness. Exact small-instance tests and outward selected-cell checks pass.
It contributes to the complete 9% certificate in
[FIRST_CLOSURE.md](FIRST_CLOSURE.md). At higher distance, the
[exact-composition refinement](EXACT_COMPOSITIONS.md) removes additional
counting slack near the sparse boundary.

The heterogeneous shuffle argument first gives the pointwise ratio

    R(a) = n^n/n! * product_j g(a_j),   g(k)=k!/k^k,

where the three category counts sum to n and g(0)=1. The earlier bound
replaces R(a) by its maximum over capped counts. That maximum
can occur far from the category counts favored by a low-output event.

Fix any nonnegative integer anchors k_j. Define

    r(k)=g(k+1)/g(k)=(k/(k+1))^k,   r(0)=1.

Because log g is concave on nonnegative integers, its forward difference
at k is a supporting slope there. For every integer a>=0,

    g(a) <= g(k) r(k)^(a-k).

For a>k, each subsequent ratio is at most r(k). For a<k, every intervening
ratio is at least r(k), so reversing their product gives the same bound.
No approximation to a factorial or asymptotic argument is needed.

Write r_j=r(k_j), tau_j=r_j/r_0, and

    B = n^n/n! * r_0^n * product_j (g(k_j)/r_j^(k_j)).

Then tau_0=1 and

    R(a) <= B product_j tau_j^(a_j).

All these constants are rational. For large n, evaluate log B with Arb
from factorials and logarithms instead of constructing its enormous
rational numerator. The anchors need not sum to n, although choosing
them near a likely category profile is a natural numerical proposal.

For the SPIN comparison, apply this inequality in each of L=256 regions,
with n=G=4096. Replace the constant factor R_C^L by B^L, and multiply the
three unnormalized input weights by tau_j before evaluating H_lambda.
The count-dependent factor is exactly a product of per-packet weights,
so it is absorbed into the same positive input polynomial. The outer
component coefficients D'_i are unchanged. The affine moment argument
continues to apply because each tau_j is fixed throughout its cell.

This does not condition on typical counts: the supporting line bounds
every integer profile. A witness should record only the three integer
anchors. Replay must derive B and tau again and round the resulting
vertex weights upward. Constant-loss witnesses remain a valid fallback.

Tests exhaust every profile and anchor for n<=8 and check domination of
a small heterogeneous shuffled distribution. Arb constants enclose their
exact rational counterparts. Replay reconstructs the constants from the
recorded integer anchors.

At the light-only endpoint for q_min=448 and theta=2/5, the new bound gives
outward log2 upper +52.01777, improving the earlier +292 but still failing.
This is a diagnostic point, not a complete occupancy class. The next
priority is full coverage at a relaxed distance, not further optimization
of this one endpoint.

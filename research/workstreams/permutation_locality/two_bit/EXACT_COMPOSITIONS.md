# Exact composition counts near the sparse boundary

At 9.5% distance, one unresolved mixture cell contains exactly one integer
composition: 400 light-only pairs, one heavy-heavy pair, and 3695 inactive
pairs. The affine cell bound gives log2 upper about +67.02. Evaluating
that composition directly, with the profile-dependent shuffle bound,
gives log2 upper below -181.66. This refinement changes the proof, not
the encoder or its setup distribution.

The composition refers to the positive comparison measure in
[MIXTURE_BOUND.md](MIXTURE_BOUND.md). It is not an identified BCH message.
The complete-code claim still requires every mixture cell and the sparse
occupancies; a successful selected composition is not a certificate.

## Enumerate every composition in a small cell

Retain the existing component features f_i, active labels b_i, G=4096,
and minimum active count q_min. In the production comparison, the only
inactive type has feature zero. There is a unique active type minimizing
the sum of its first two feature coordinates. Call its feature f_* and
write a=f_*,1+f_*,2=1/7.

For any active type, define its excess c_i=f_i,1+f_i,2-a. The minimizing
type has excess zero; every other active type has positive excess.
Fix a rectangular cell C and an integer composition n whose mean lies in C.
If q=sum_i b_i*n_i>=q_min, then

    sum_{i active} n_i*c_i
        = G*(m1+m2)-a*q
        <= G*(upper_C(m1)+upper_C(m2))-a*q_min.

This bounds the number of nonminimal pairs near the lower occupancy
boundary. Enumerate all their nonnegative integer counts within the excess
budget. For each choice, the remaining minimal-type count must satisfy
all three coordinate intervals and q_min<=q<=G. Their exact intersection
is an integer interval; enumerate all its integers. The inactive count
is then G minus the active count.

The third coordinate also bounds the number of pairs containing central
rows. Bound those pairs using their minimum positive third coordinate;
bound the remaining nonminimal pairs using their minimum positive excess.
The sum remains an upper bound, even though both classes consume the
same excess budget. Use the smaller of this bound and the excess-only bound.

`boundary_compositions` returns either the complete list or no result.
It returns no result if these bounds permit more than sixteen nonminimal
pairs, or the search exceeds 20,000 nodes or 32 resulting compositions.
A budget limit never produces a truncated list that can pass verification.
The empty list is distinct from an inconclusive search.

## Bound each composition without separating its coordinates

Fix a complete composition n. Component i has coefficient d_i and packet
category probabilities pi_i=(pi_i,0,pi_i,1,pi_i,2). Fix a positive input
tilt z and an output tilt lambda>0. Put

    Z_i = pi_i dot z,
    w_j = sum_i (n_i/G)*pi_i,j/Z_i.

The coefficient for all labeled assignments of these component counts is

    M(n) = G! / product_i n_i!.

Use the constants B and tau from [DENSITY_TANGENT.md](DENSITY_TANGENT.md),
which dominate every regional packet profile. Let H_lambda be the exact
positive input polynomial from [AFFINE_MOMENT.md](AFFINE_MOMENT.md).
The contribution of this composition to the first moment is at most

    exp(lambda*d) M(n) product_i d_i^(n_i)
        * (product_i Z_i^(n_i))^256 * B^256
        * H_lambda(w0*tau0, w1*tau1, w2*tau2).

Here d is the inclusive bad-output-weight cutoff. The polynomial has degree
G*256, the number of packets, and includes the expectation over independent
transvections with fixed expansion and feedback maps. Its existing
two-state envelope supplies the outward bound. Irrational tilted weights
are rounded upward before normalizing; positivity justifies that rounding.

For comparison, each composition also gives exact category caps
h_j=sum_{i:pi_i,j>0} n_i. The constant capped shuffle bound uses these caps.
For the troublesome composition, h=(4096,401,1). This removes more than
1345 bits of loss from the universal three-category bound, but still gives
log2 upper +311.91. Retaining the count-dependent shuffle factor is essential
to the successful -181.66 bound.

The successful reference counts are (3892,203,0). Replacing the final zero
by one gives only about 15.87 bits of margin for this composition. A zero
reference count does not exclude double-bit packets: the supporting-line
inequality bounds every actual profile. The search now tries zero reference
counts for rare categories, rather than only rounding their means.

The input means can still be poor reference counts after weighting outputs
by exp(-lambda*weight). A bounded numerical refinement estimates category
counts from logarithmic derivatives of the tilted envelope. It rounds
these estimates to integers and reoptimizes the input and output tilts.
The derivative calculation proposes witnesses only. Every returned count
is checked by the existing pointwise inequality, which permits arbitrary
integer reference counts. No conditioning or additional randomness enters
the proof.

This refinement recovers about 100 bits for selected nearby compositions,
but some 9.5% bounds remain above one. [HILL_CLIMB.md](HILL_CLIMB.md)
records both the improvements and remaining failures.

## Whole-cell verification

For each exact-composition leaf, replay reconstructs the full integer list
from the original split-tree cell. It rejects missing or duplicate entries.
It recomputes every required contribution outward and sums them. A parent's
list may be reused on a child; replay selects exactly the reconstructed
child list and ignores the other entries. No cached bound or floating
infeasibility claim is trusted.

Small exhaustive tests compare the returned list with every composition,
check both budget limits, and reject incomplete witnesses. Additional tests
check the normalization, the capped-loss factor, and the production cell's
unique integer composition. The ordinary split-tree verifier still proves
coverage of the full domain.

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --minimum-groups 401 --threshold 199229 --exact-boundary --probe-cell 3677/262144 1839/131072 15/262144 1/16384 0 1/262144
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --resume tmp/two-bit-affine-dense-401-d095-exact.json --max-cells 2500 --max-depth 54 --target-bits 60 --exact-boundary --polish-tilt --reuse-witnesses --output tmp/two-bit-affine-dense-401-d095-exact-rare.json
```

The second command is a full-domain attempt, not a completed certificate.
Its subsequent refinement reached 222 accepted leaves and 704 unresolved
cells before adding the tilted reference-count search.
The complete 9% result and the separate four-bit work are preserved.

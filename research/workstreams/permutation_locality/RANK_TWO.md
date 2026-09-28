# All rank-two messages in one group: 53.38 bits

2026-09-24. For the uniform-shared-coordinate g=4,c=4 distribution in
RANDOM_BLOCKS_GFNI.md, an outward calculation bounds all rank-two messages
confined to one group by 8.519283e-17. This is approximately 53.382 bits for
output weight at most 209715, at K=2^20 and N=2^21. The union includes every
possible active group among the 2048 groups.

The probability is over the ideal uniform route and independent IMT
transvections, with the BCH code and expansion/feedback maps fixed. These
are the same distribution and maps used by the preceding rank-one bound.
This result does not cover rank-three or rank-four tuples, or messages in
multiple groups. It makes no additional claim about speed: the measured
candidate remains about 6.49 ms versus 9.62 ms for production.

## What the type calculation taught us

A two-dimensional binary subspace has exactly three nonzero words x, y,
and x+y. Let their weights be a,b,c, respectively, and let u=(a+b+c)/2.
The column types (0,0), (1,0), (0,1), (1,1) occur

    (256-u, u-b, u-a, u-c)

times. A shared uniform coordinate shuffle therefore has an exact
composition, not four independent choices at every coordinate.

`rank_two_types.py` bounds the number of subspaces of each sorted weight
triple using the authenticated BCH shell caps. For one ordering (a,b,c),
the number of ordered pairs is at most both

    A_a (A_b - 1[a=b])

and

    A_a binomial(a,(a+b-c)/2) binomial(256-a,(b+c-a)/2).

The second expression fixes x and counts all binary candidates y, without
requiring them to belong to BCH. Minimize over the weight orderings and
divide by the product of factorials of repeated-weight multiplicities.
This gives an integer upper bound on the subspace count. Each subspace
has 210 ordered spanning four-row tuples. Exhaustive small-code checks
validate the composition identity and these multiplicity factors.

An initial moment screen sampled column types independently and then
conditioned on the required counts. That upper bound was vacuous for
sparse types, even though conditioning costs only about ten bits. The
unconditioned moment is dominated by atypically small supports, and the
conditioning inequality does not discard their contribution.

The next screen fixes union support exactly, using independent labels only
for its nonzero columns. Conditioning then restores the required nonzero
label counts. The smallest feasible type, (38,38,38), improves from a
vacuous upper bound to about 62.40 bits. That diagnostic already includes
all row assignments and group locations for this one type. Five other
tested types also clear 40 bits. These are binary64 screens, not certificates
and not a union over the 87044 feasible weight triples.

## A simpler envelope covers all rank-two types

The successful support-conditioned screen suggests keeping exact union
support while discarding the nonzero labels entirely. This relaxation is
enough for rank two; no type-by-type union is needed in the final bound.

A nonempty 16-bit bundle has four column weights in {0,1,2,3,4}. Up to
column order, there are 69 nonempty shapes. Conditional on a shape,
the column order is uniform, and the row labels in each column are
permuted independently. The active group occupies one uniform epoch
and one uniform aligned window within that epoch.

The local census enumerates all 69 orbits against the actual E and B maps.
For each orbit it computes exact weighted-output counts by expansion level,
feedback-syndrome multiplicities, and cancellation-output counts. These
construct the same seven-coordinate moment envelopes as GROUP_DISTANCE.md.
The only shapes permitting Bx=0 are (0,2,3,3) and (1,2,2,3); these events
remain in the transfer from zero rather than being discarded.

Let R_sigma(z) be the resulting macroregion transfer for shape sigma,
averaged over the group's uniform placement among the 256 epochs. For
b=1,2,3,4, define V_b(z) by taking the entrywise maximum of R_sigma(z)
over shapes having exactly b nonzero columns. Set V_0(z)=R_0(z), the
empty-macroregion transfer. Every entry is nonnegative.

Fix a message tuple whose union support has size u. Conditional on its
support placement, its nonzero column labels may be arbitrarily dependent.
Each macroregion transfer is nevertheless bounded entrywise by V_b for
its occupancy b. Products preserve this order. Thus the envelope does not
assume that labels are independent, nor that one shape maximizes every
entry. The entrywise maximum deliberately permits a larger process.

The shared coordinate permutation makes the union support a uniform
u-subset of 256 positions. Condition on the first occupied macroregion:
let b be its occupancy and l the number of following macroregions. Define

    V(x,z) = sum_{j=0}^4 binomial(4,j) x^j V_j(z).

The conditional output moment is bounded by

    e_0 V_b(z) [x^(u-b)] V(x,z)^l 1 / binomial(4l,u-b).

As in the rank-one calculation, multiply by z^(-209715), minimize over
seven fixed tilts, and cap at one. Average over the first occupied
macroregion with weights

    binomial(4,b) binomial(4l,u-b) / binomial(256,u).

Call the resulting upper probability p(u). For CDF summation, replace it
by the decreasing majorant q(u)=max_{v>=u} p(v). This preserves an upper
bound, whether or not the numerical p values are monotone.

Let T_h(u) be the previously verified upper CDF on rank-h four-row tuples
with union support at most u. It combines ordinary-spectrum basis counts
with shortened-code dimension bounds. Summation by parts gives

    2048 [T_h(256) q(256)
          + sum_{u=0}^{255} T_h(u) (q(u)-q(u+1))].

All CDF coefficients are nonnegative. We do not interpret differences
between upper CDF values as exact shell counts. No independence between
the message events or possible active groups is required.

## Results and remaining slack

Independent outward runs at 192 and 384 bits agree on the following values:

| One active group's rank, union over all group locations | Bound expressed as margin |
|---|---:|
| 1, conservative shape envelope | 44.236 bits |
| 2, all tuples | **53.382 bits** |
| 3, all tuples | 30.216 bits |
| 4, all tuples | Vacuous |

The more specialized rank-one calculation remains stronger at 45.238 bits.
Rank two is dominated by small-support CDF terms, near u=57 through 63.
Rank three's largest terms occur near u=70 through 74. Rank four is instead
dominated by the dense tail, where the entrywise shape envelope and the
tuple multiplicity are too expensive. These are limitations of this bound,
not codeword counterexamples.

The subsequent RANK_THREE.md calculation retains rank-dependent column
information and verifies 40.643 bits for rank three in outward replays at
192 and 384 bits. Rank four and multiple active groups still require analysis.
The 2x performance objective also remains unmet.

## Reproduction

All scripts are isolated research tools and write no experiment files.

    python -B research/workstreams/permutation_locality/rank_two_types.py
    python -B research/workstreams/permutation_locality/rank_two_moment.py --fixed-support
    python -B research/workstreams/permutation_locality/random_group_moment.py
    python -B research/workstreams/permutation_locality/random_group_verify.py --precision 192
    python -B research/workstreams/permutation_locality/random_group_verify.py --precision 384

The first three commands provide exact counting tests and binary64 diagnostics.
Only the last two use outward arithmetic for the final probability bound.
The verifier authenticates 163 BCH dependency files and verifies 67 exact
rational shortening witnesses. It reconstructs all local orbit counts,
checks direct overlap counts for selected input patterns, and compares
every rank-one orbit transfer against the preceding implementation.
It checks the rank-two inequality directly against 2^-50, rather than
accepting a rounded logarithm. It does not replay every historical BCH proof.

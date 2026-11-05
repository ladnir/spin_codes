# Sparse updates in one interpolation half

The alternating expansion gains 193 nonzero feedback symbols per unchanged
state pair. The next question is what a sparse state update can destroy.
Here the actual state64 maps give a lower bound of150 out of256 when the
intervening output has weight one and both states are nonzero. This is a
local theorem over F=GF(2^128), not a whole-code distance certificate.

Use the [alternating construction](ALTERNATING_MAPS.md) without changes.
The outer, route, and setup distributions remain fixed. The bounds below
hold for every nonzero compression scaling R, not only on average.

## The compression image has known roots

Write s=64 and let V be the degree-s subspace polynomial for the first
s Cantor points. Split Y_0 into C_0 and C_1, each of size s, where V equals
0 and1, respectively. On the state evaluation set Z, V equals alpha, with
alpha not in {0,1}. On the two s-point halves of Y_1 it equals alpha and
alpha+1. These are four disjoint cosets.

For a finite set S, define g_S(z)=product_(x in S)(z-x). Fix a nonempty
output support W contained in C_e, with w=|W|<=s. Column scaling by R is
invertible on W, so absorb it into the supported output y. Let d(z) be the
degree-<s polynomial representing B y. Then its possible values are exactly

    d(z)=g_(C_e minus W)(z) v(z),     deg(v)<w.          (1)

To prove this, interpolate y on Y_0. Its polynomial f has degree<2s and
vanishes on Y_0 minus W. Hence

    f(z)=V(z)(V(z)+1) u(z)/g_W(z),     deg(u)<w.

Evaluating on Z replaces V(V+1) by alpha(alpha+1). The polynomial

    [alpha(alpha+1)/(alpha+e)] [(V(z)+e)/g_W(z)] u(z)

has degree<s and agrees with f on Z. Uniqueness proves (1), including
surjectivity as u varies. This argument is valid over F, although every
fixed map coefficient lies in E=GF256.

One consequence is that A_1 B restricted to any columns in one half is an
[128,w,129-w] generalized RS generator. Any h selected rows have rank
min(h,w). Cross-half minors need not have that rank; the retained toy
counterexample uses one column from each half.

## Bound holes before and after the update

Let p and p'=p+d be the two state polynomials. Require both nonzero.
First observe p on Y_0 and p' on Y_1. Let H_0 be the exact holes of p on
Y_0, and let h_1 count the holes of p' on Y_1. Define

    h_0=|H_0|,
    r=|H_0 intersect (C_e minus W)|,
    i=|H_0 intersect W|.

Then 0<=h_0,h_1<=s-1. Equation (1) shows that p' retains r roots outside
Y_1. Polynomial root counting therefore gives

    h_1 <= s-1-r.                                     (2)

A second bound uses the other interpolation half. Put
U=H_0 intersect C_e and M=C_(1-e) minus H_0. Thus |U|=r+i and
|M|=s-h_0+r+i. Factor p=g_(H_0) u_0, where deg(u_0)<=s-1-h_0.
On either half of Y_1, replace V by its constant value eta. Clearing the
nonvanishing denominator g_M g_W in p'=p+d gives the polynomial

    (eta+1-e) g_U g_W u_0 + (eta+e) g_M v.             (3)

Its degree is at most

    D=s-1-h_0+r+i+w.

Polynomial (3) cannot vanish identically. Otherwise p' would vanish on
an entire s-point half of Y_1, contradicting p'!=0 and deg(p')<s.
Each half therefore contributes at most D roots. Together with (2),

    h_1 <= min(s-1-r, 2(s-1-h_0+r+i+w)).                (4)

All divisions in the proof are only used at points outside Y_0, where
their denominators are nonzero. Equation (3) need not equal one polynomial
across both halves: root counting is performed separately on each half.

For the opposite starting phase, exchange p and p'. Characteristic two
gives p=p'+d, so the same total pair-weight bound applies. The physical
update is unchanged. Zero endpoints are excluded locally, not imposed at
the end of the full encoder.

Enumerating the feasible integer profiles in (4) gives:

| Nonzero output weight w, confined to one half | Feedback pair weight at least |
|---:|---:|
| 1 | 150 / 256 |
| 2 | 149 / 256 |
| 4 | 146 / 256 |
| 8 | 141 / 256 |
| 12 | 135 / 256 |
| 16..64 | 130 / 256 from the separate MDS bounds |

The helper checks the displayed weights. The last row uses the universal
130 bound; no assertion of exact minimum or monotonic optimization is
needed. Every weight-one update is automatically confined to one half.
The unchanged-state bound remains193. No new encoder operation buys these
bounds: they analyze the existing alternating variant.

## Remove zero-state and zero-update fibers exactly

A rank-deficient constraint matrix need not admit a live transition. Fix
hole sets H_0 subset Y_0 and H_1 subset Y_1, and an allowed support W in
one half. Here the hole sets need not be exact. Assume h_0,h_1<s;
otherwise one observed state must be zero and there are no live-to-live
solutions. Let d_K be the dimension of
the linear kernel on (Q,y_W) of

    A_0[H_0,:] Q=0,
    A_1[H_1,:] (Q+B[:,W]y_W)=0.

Use h_0=|H_0|, h_1=|H_1| and the same r as above. The kernel dimensions
when the old state, new state, or update is zero are, respectively,

    a=max(0,w-h_1),
    b=max(0,w-h_0+r),
    c=max(0,s-h_0-h_1).                               (5)

For a and b, use (1) and count polynomial roots outside its fixed root
set. For c, the two observations evaluate the same polynomial on disjoint
sets. Every pairwise intersection of these three zero fibers is {0},
because B restricted to W is injective. Inclusion-exclusion therefore
gives exactly

    q^d_K - q^a - q^b - q^c + 2,      q=|F|,          (6)

solutions with both states and the update nonzero. Individual coordinates
of y_W may still be zero. Formula (6) counts allowed support, not exact
support; it is an upper bound when exact support W is required.

A useful rank-free upper bound for d_K is

    min(max(0,s-h_1-r)+max(0,w-h_0+r),
        s-h_0+max(0,w-h_1)).                           (7)

For the first term, choose p' with its h_1+r prescribed roots; the fiber
of p has dimension max(0,w-h_0+r). For the second, choose p with H_0
roots first; the p' fiber has dimension max(0,w-h_1). Inconsistent choices
only decrease the dimension. Replacing d_K by (7) in (6) preserves an
upper bound. Equation (4) supplies an additional impossibility test for
live endpoints. Neither bound asserts generic multitime rank.

These formulas provide a local split count for the global proof. They do
not permit independently resampling Q or multiplying counts at shared
state boundaries. That composition and the complete route sum remain open.

## Three complete restricted trajectory classes

The local weight theorem can already strengthen the retained
[output-interval proof](OUTPUT_INTERVAL.md), without independently counting
states. Fix exactly one episode of v batches, including its activation.
Let z<=32 outside batches have zero output and zero incoming state. Every
other outside batch is nonzero and quiet, so j=16384-v-z. Set L=v-1.

Remove the activation batch and partition the L live-incoming batches into
consecutive disjoint pairs, starting at the first such batch. Leave the
last batch unpaired if necessary. Both incoming states of each pair are
nonzero. Its first output batch is precisely the update between them.

If at least g pairs have weight-one first output, the feedback support S
satisfies |S|>=65L+20g, by (4). If at least g pairs have zero first output,
then |S|>=65L+63g by the unchanged-state theorem. These statements hold
regardless of the absolute starting phase.

For a fixed output candidate, inverse membership has probability at most
q^(d_C(A_live)-|S|). The route-uniform interval bound remains

    d_C(A_live)<=8192 max(0,floor((L+126)/64)-128),

with dimension zero for L=0. Quiet constraints leave at most
2^(128(v+j))*q^(W-64j) candidate outputs, where W=1029701. Consequently,
for a specified interval and outside-zero set, the log2 failure bound is

    128(v+j)+128[W-64j+d_C(A_live)-65L-bonus*g].        (8)

There is no extra union over qualifying pairs. The class is defined by
having at least g qualifying slots in the fixed pairing; every output
in that class obeys the same support lower bound. The output count in (8)
is a relaxation of that entire class, not a count for one selected witness.

An integer-only replay covers all feasible (v,z) in each row:

| Complete restricted class | Cells | Largest per-pattern log2 bound | Failure after all placements |
|---|---:|---:|---:|
| v<=8700; at least 90% of canonical pairs have weight-one update | 269742 | -706688 | <=2^-706206 |
| v<=8800; every canonical pair has weight-one update | 273042 | -770688 | <=2^-770206 |
| v<=10000; at least 75% of canonical pairs have zero update | 312642 | -73216 | <=2^-72734 |

Percentages mean ceilings of the corresponding fractions of floor(L/2).
The count excludes only cells with 65j>W. For each row all interval locations,
outside-zero positions, supports, and amplitudes are included. At most
16384^2 *33*16384^32 <2^482 placements suffice, including every v and z;
the cell count is not another union factor. The bounds hold for every
route and every nonzero R realization, and average only over D.

These classes overlap and do not cover all one-episode trajectories.
Multiple episodes and dense updates in both halves remain open. The large
margins in the table must not be reported as a whole-code setup margin.
The terminal state remains free throughout.

## Literature comparison and next gate

[Hutchinson, Rosenthal, and Smarandache, *Convolutional Codes with Maximum
Distance Profile*](https://arxiv.org/abs/math/0307196) characterize optimal
column distances through nontrivial minors of a truncated system matrix.
This motivates checking multitime constraints, not inferring them from the
block MDS property. Their generic-existence result does not certify our
particular maps, randomized recurrence, or termination convention.

[Alfarano, Napp, Neri, and Requena, *Weighted Reed-Solomon convolutional
codes*](https://arxiv.org/abs/2012.11417) give a weighted-Vandermonde
construction and sufficient field-extension conditions for MDP. Their
Theorem3.5 requires an extension degree exceeding a determinant-degree
parameter. We have not established those hypotheses for this inner.
Importing the theorem or its coefficients would require a different
construction/proof and a separate cost analysis.

The restricted classes above pass, but extending only their pair-weight
bonus cannot close the present envelope. For one episode occupying all
16384 batches, even assuming every canonical pair has zero update gives
only |S|>=1580928. Replacing d_C(S) by the full outer dimension 1048576
leaves W+d_C(A_live)-|S|=497349 field symbols; the bound is vacuous.
This is a failure of the relaxation, not an exhibited low-weight codeword.

The [actual-support refinement](FEEDBACK_SUPPORT.md) now closes this
full-length class when every canonical pair has zero output update,
at >49.1% with failure below 2^-61. It retains S in the regional moment
and charges the outputs creating its holes, with explicit mask counts.
Dense updates in both halves and multiple episodes remain open. Retain
the exact zero-fiber subtraction before relaxing local ranks, and do not
multiply local counts at shared state boundaries.

## Checks

Run serially from the workstream directory:

    python -B coded_state/test_half_update_bound.py
    python -B coded_state/half_update_bound.py
    python -B coded_state/alternating_interval_certificate.py

The first exhausts524272 projectively normalized GF16 one-update live pairs,
checks six zero-fiber counts by enumeration, and checks304 GF256 fixtures
with constructed holes at state sizes4,8,16,32,64. The toy minimum11 is
stronger than this general root bound and is not extrapolated. The second
replays the actual-size integer optimization for the table. No generated
files or native timings are involved. The proof, rather than random
fixtures, establishes (1)--(7) over the full field.
The last command checks every integer cell in the restricted-class table
and its full placement factor. A separate JavaScript integer calculation
gave the same counts and maxima before the Python replay.

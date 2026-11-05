# Coarse state labels with exact routing counts

The joint-cohort calculation retains exact input counts inside a narrow
band. Extending that matrix literally to all513 input counts would spend
most of its memory on states outside the difficult middle range. This
note describes a conservative alternative. It changes only the proof
representation, not the encoder or its sampled setup.

The construction below is a matrix domination argument. Its continuation,
[coupled population bounds](QUARTET_COUPLED_ROW_CAPS.md), implements the
all-count matrix but finds a large loss from independent endpoint choices.
No distance certificate is claimed here.

## Keep the quota label on the transition

Partition the incoming histogram space into cells i=(A,R,S). A,R,S are
integer intervals containing a,n1,n2, with n1+n2<=a. They need not have
the same widths. Every admissible histogram must belong to exactly one
cell, including the count tails.

An edge label e contains an EXACT next input count l_e and an outgoing
cohort cell. Let j(e) be its coarse destination cell. Several edge labels
may have the same destination, but their l_e values remain separate.

Let K(n,e) be the exact potential-weighted transition coefficient from
incoming histogram n to the histograms represented by edge e. Construct
a nonnegative rectangular matrix M such that

    K(n,e) <= M_(i,e)  for every n in cell i.

One valid construction takes coefficientwise maxima of the existing
exact-a joint-bucket bounds over all a in A, with the incoming boxes
clipped to n1+n2<=a. This bound can be loose. It is nevertheless uniform
over the incoming cell. Do not maximize over an empty clipped box.

For any nonnegative function f of the destination cell,

    sum_e K(n,e) f(j(e)) <= sum_e M_(i,e) f(j(e)).

Induction gives the corresponding upper bound for every labelled path.
This argument does not claim that the coarse cells form an exact Markov
chain or that the omitted counts are independently resampled.

## Count-aware positive-vector bound

Fix a positive vector v on cells and z>0. Write d_e=l_e-408, and define

    rho = max_i sum_e M_(i,e) z^d_e v_(j(e)) / v_i,
    P_(i,e) = M_(i,e) z^d_e v_(j(e)) / (rho v_i).

Then every row sum of P is at most1. For t>=0, let p_(t,s)(i) be the
sum of P-weights of t-edge paths from i with total count deviation s.
With p_(0,0)=1 and other initial coefficients zero, the recurrence is

    p_(t+1,s)(i) = sum_e P_(i,e) p_(t,s-d_e)(j(e)).

After sixteen updates, apply the same regional binomial normalization
and payment marker as in [the regional bound](QUARTET_REGIONAL_FRESH_TRANSFER.md).
The proof uses the exact l_e values, even when j(e) stores only an interval
containing l_e. Replacing d_e by a representative of that interval would
break the quota argument and is not allowed.

The actual initial row remains separate and is used once. The free
terminal factor is bounded uniformly over each entire destination cell
and is also applied once. No reset is inserted at a regional boundary.

For the terminal bound, its logarithm is affine in a,r,s after maximizing
over the remaining two old types. Maximize over the vertices of

    a in A, r in R, s in S, r+s<=a.

Vertices either lie on an endpoint of A or satisfy a=r+s with r,s at
interval endpoints. This gives a finite uniform terminal bound without
choosing a favorable representative histogram.

`quartet_fresh_regional.coefficient_rows` accepts an optional target_map
for this rectangular representation. Its square-matrix behavior is the
special case j(e)=e. The added exact test enumerates all four-edge paths
on a two-state, three-edge example, including two edges with the same
destination but different count labels. This prepares the counter layer;
the full-width coarse-state tail matrix is evaluated in the continuation.

## A cheaper exact routing-coefficient recurrence

The joint transition uses

    A_(a,l,u) = sum_v binom(w-a,v) J_(l,u+v)/binom(w,u+v)
                     * d0^v h1^(l-u-v).

Here d0=h1 by the definition of the histogram potential. Hence

    A_(a,l,u) = h1^(l-u)
               * sum_v binom(w-a,v) J_(l,u+v)/binom(w,u+v).

Pascal's identity therefore gives, for 0<=a<w and admissible u,

    A_(a,l,u) = A_(a+1,l,u) + h1 A_(a+1,l,u+1).

The base is A_(w,l,u)=h1^(l-u) J_(l,u)/binom(w,u), with zero outside
the support. The recurrence counts every original summand and introduces
no probability comparison. It obtains all requested incoming counts from
one vectorized pass, rather than repeating the binomial sum for each pair.

`quartet_joint_regional.allowed_tables` implements the recurrence in log
coordinates and now supplies the joint driver. The coefficient tilt used
by the joint kernel is applied only after this recurrence, so it does not
alter the displayed identity.

The test source compares exact toy coefficients for h1=1 and h1=3/2,
then compares the floating implementation with the original direct sum at
counts0,1,127,392,408,424,511,512. All680 exact recurrence comparisons,
3770 exact coarse-count geometry comparisons, and12023 FLOAT coefficient
regressions pass. This validates the evaluation change at the stated
scope, not the full-width tail probabilities.

The rectangular counter passes324 exact labelled-edge paths with two
terminal vectors. The existing64 square paths,8692 uniform joint-bucket
comparisons,81 two-region quota pairs, and1155 toy coefficient checks
still pass. The whole audit takes0.43s remotely and writes no files.

The revised evaluator and counter reproduce the width16 four-count screen
at-9577.538462675759, within1.4e-8 of the previous log2 expression. Its
runtime is16.00s. The fine full-band result in the joint-cohort note was
computed with the previous equivalent evaluator, not rerun here.

## Uniform geometry across an incoming count interval

Positivity in the recurrence implies A_(a,l,u)<=A_(a-1,l,u) whenever
u<=a-1. Thus, for any a in [a_lo,a_hi] and any feasible u<=a,

    A_(a,l,u) <= A_(max(a_lo,u),l,u).

Using a_lo without the maximum with u would be wrong: the coefficient
at u>a_lo lies outside that smaller state's support. The displayed bound
keeps those feasible contributions. It can supply the geometry factor
of a coarse incoming cell without enumerating every a separately.

The old-cohort polynomial needs a separate uniform bound over a. For a
fixed fresh count r, its remaining degree count n=a-r varies in an interval.
For a coefficient of degree d, ignore n<d only because its true coefficient
is zero. The feasible (n,s) set is

    max(a_lo-r,d,s_lo)<=n<=a_hi-r,
    s_lo<=s<=s_hi, s<=n.

At every positive coefficient-evaluation point z, the logarithm of
L2^s max(L3,L4)^(n-s) is affine on this polygon. Maximize over its vertices
BEFORE minimizing over z. A point-mass correction likewise needs a variance
lower bound uniform over the polygon. This supplies a route to coarse-a
local bounds; its full-width implementation and numerical effect are open.

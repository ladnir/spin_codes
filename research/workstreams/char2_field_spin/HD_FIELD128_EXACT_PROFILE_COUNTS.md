# Exact counts for the two alternating matching profiles

The current support bounds use separate matching marginals and cyclic run counts.
Their minimum need not equal the joint count. This note derives that joint count
without changing the field, inner maps, route, or transfer matrices.

Fix an integer P>=1. A support is an indexed binary word
`x=(x_0,...,x_(2P-1))`, with indices taken modulo 2P. Define

    r = sum_i x_i,
    a = sum_j x_(2j) x_(2j+1),
    b = sum_j x_(2j+1) x_(2j+2).

Let N(P,r,a,b) count the words with this profile. Rotations are not identified.
The construction uses P=4096, but the counting argument applies to every P>=1.
An infeasible profile has count zero.
The feasible matching marginals satisfy `0<=r<=2P` and
`max(0,r-P)<=a,b<=floor(r/2)`; these conditions alone need not make the joint
profile feasible.

## A positive count for each phase flow

Fix a profile and let u count the even pairs of type 10. Put
`s=r-2a` and `v=s-u`, so v counts the even pairs of type 01. Represent a word
as a closed walk on the phase vertices E0,E1,O0,O1. A vertex specifies the
index parity and its bit. The edge multiplicities are

    even-to-odd edges:                 odd-to-even edges:
             O0       O1                       E0          E1
       E0  P-r+a       v                 O0   P-r+b       a+u-b
       E1    u         a                 O1   a+v-b         b

All multiplicities must be nonnegative. Together with feasible matching
marginals, this condition gives

    max(0,b-a) <= u <= min(r-2a,r-a-b).

The graph is balanced: each vertex has equal incoming and outgoing degree.
Remove vertices of degree zero. Write d_w for the degree of an active vertex w,
and n_e for an edge multiplicity. Let tau be the number of directed spanning
trees oriented toward any fixed active vertex, counting parallel edge choices.
In a connected balanced graph this number is independent of the root.
Set tau=0 when the active graph is disconnected.

The number of indexed words with this fixed u is

    T_u = P tau * product_active w (d_w-1)! / product_edges e n_e!.

To prove the formula, first label all parallel edges. For a walk starting at
vertex w, the last outgoing edges of the other vertices form a tree toward w.
After choosing this tree, order each vertex's other outgoing edges arbitrarily.
There are `d_w tau product(d_i-1)!` such orders. Following these orders gives
an Euler walk: the chosen last-exit tree prevents a premature closed component.

Each unlabelled indexed word has exactly `product n_e!` edge labellings.
Divide by that factor, then sum over starting vertices E0 and E1. Their degrees
sum to P, giving the displayed formula. The walk starts at the fixed index 0;
there is no additional rotation quotient, even for a periodic word.
If the active graph is disconnected, no single walk can traverse every edge.

For four active vertices, a positive expression for the tree count rooted at E0 is

    tau = u(P-r+b)(a+v) + a(a+v-b)(P-r+a+u).

These are the four possible directed tree monomials, grouped into two terms.
Define

    e0=P-a-u,  e1=a+u,  o0=P-a-v,  o1=a+v,
    e=a+u-b,   f=a+v-b.

An equivalent formula avoids large factorial intermediates:

    T_u = P tau C(e0,v) C(e1,a) C(o0,e) C(o1,b) / (e0 e1 o0 o1).

The numerator is divisible by the denominator. This formula applies only when
all four degrees are positive. If a phase vertex is absent, that phase is
constant. Choosing the remaining phase's one positions gives a single binomial
coefficient. In particular, each all-constant word has count 1. When
`r=P,a=b=0`, the two endpoint flows each give one alternating word; every
four-vertex interior flow is disconnected and gives zero.

Finally,

    N(P,r,a,b) = sum_u T_u

over the feasible integer interval above. Every term is a nonnegative integer.

## Exact recurrence and symmetries

For consecutive four-vertex flows, let
`Q_u=product(d_w-1)!/product n_e!`. With the quantities evaluated at u,

    Q_(u+1)/Q_u = e1 o0 v f / ((e0-1)(o1-1)(u+1)(e+1)).

Thus a nonzero T_u gives the next term by multiplying this ratio and
`tau_(u+1)/tau_u`. An implementation must not divide by a zero tree count.
It must also handle endpoints whose active vertex set changes directly.

Three bijections reduce the work. A one-position cyclic shift
exchanges a and b. Bit complementation maps the profile to
`(2P-r,P-r+a,P-r+b)`. The reflection `x'_i=x_(1-i)` preserves a and b and
exchanges u and v. Consequently the sum can retain only `u<=s/2`, doubling
each term except a possible central term. These are bijections of indexed
words, not identifications of words.

The isolated scalar helper uses these symmetries, four-binomial initialization,
and exact-divisibility checks. It supports P<=4096 and retains at most 4096
cached scalar results. Tiny exhaustive tests compare it with indexed words
and an independent tree enumeration, including disconnected and constant cases.

## Use as an optional support refinement

For any explicitly selected set of profiles, replace the previous count upper
bound by N on those profiles and retain the previous bound everywhere else.
All profiles remain in the transfer sum. This is valid even when a numerical
screen selected the profiles: the exact count identity holds for every profile.
The selection does not change the construction or the outer prefix.

This note and scalar helper supply no new coefficient computation, selected-weight
probability bound, whole-code certificate, or encoder performance measurement.

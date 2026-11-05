# A cheap upper bound for every joint matching profile

The [exact joint count](HD_FIELD128_EXACT_PROFILE_COUNTS.md) sums over the
orientation count of the even mixed pairs. Evaluating that sum separately
for every profile is expensive. This note bounds the sum using its largest
term. The bound changes only the proof calculation, not the construction.

Use the indexed words and fixed-flow counts `T_u` from the exact-count note.
Complement and shift the word, if necessary, to assume `r<=P` and `a>=b`.
Put `A=P-r`, `s=r-2a`, `v=s-u`, and `c=a-b`. All are nonnegative on the
feasible interval `0<=u<=s`. Reflection gives `T_u=T_(s-u)`.

Except when `r=P,a=b=0`, the terms are nondecreasing up to the midpoint.
Consequently, writing `N` for the exact joint count,

    N(P,r,a,b) <= (s+1) T_floor(s/2).

The exception consists of two alternating words, so its count is 2.
The constant profile has count 1. A profile with no words has bound zero.
For every nonempty profile, the midpoint bound is at most `(s+1)N`.
At P=4096 this factor is at most 4097, or about 12.001 bits.

## Why the midpoint is largest

First suppose `a>=1` and `A+b>=1`. All four phase vertices have positive
degree. Up to a positive factor independent of u, the factorial quotient
in `T_u` is `h(u)h(s-u)`, where

    h(t)=(a+t-1)! (A+a+t-1)! / (t! (c+t)!).

The ratio

    h(t+1)/h(t) = ((a+t)/(t+1)) ((A+a+t)/(c+t+1))

is nonincreasing in t: its two factors have this property because
`a>=1` and `A+b>=1`. Hence `h(u)h(s-u)` is nondecreasing for
`u<s/2`. The tree factor can be written

    tau_u=(A+a+b)u(s-u)+a(r-a-b)(A+a).

It is nonnegative and also nondecreasing up to the midpoint. Their
product therefore has the required monotonicity. A vanishing tree
factor counts no connected walk and causes no division by zero.

Two boundary cases remain. If `A=b=0` and `a>=1`, the tree factor is
`a(a+u)(a+v)`. Cancellation in the consecutive-term formula gives

    T_(u+1)/T_u = (a+u)v / ((a+v-1)(u+1)).

For `u<v-1`, this ratio is at least 1 because
`(a+u)v-(a+v-1)(u+1)=(a-1)(v-u-1)>=0`.
If `a=b=0` and `A>=1`, the interior ratio instead equals

    T_(u+1)/T_u = (A+u)v / ((A+v-1)(u+1)),

which is at least 1 by the same argument. The endpoint has a missing
phase vertex. For `r>=2`, direct counting gives
`T_0=C(P,r)` and `T_1/T_0=Ar/(A+r-1)>=1`, since
`Ar-(A+r-1)=(A-1)(r-1)>=0`. For r=0 or r=1 the midpoint is already
an endpoint; no ratio is needed. Finally, `A=a=b=0` is precisely the
alternating exception stated above.

## Constant work after one-dimensional preparation

Put `L=floor(r/2)` and `H=ceil(r/2)`. At the midpoint,
`u=L-a` and `v=H-a`. The four degrees are `P-L,L,P-H,H`, independent
of a and b. When r>=2, define

    U_a=C(P-L,H-a) C(L,a),
    V_b=C(P-H,L-b) C(H,b),
    D=(P-L)L(P-H)H,
    tau=(L-a)(P-r+b)H+a(H-b)(P-H).

Then `T_floor(s/2)=P tau U_a V_b / D`. The quotient is an integer.
Two one-dimensional arrays suffice; each subsequent profile needs only
integer products and an exact division. The r=0, r=1, and alternating
cases are handled directly. Complement and matching exchange preserve
the count, so the formula covers the full weight range.

The implementation takes the minimum with the previous valid profile
cap. Its floating-point interface is only a numerical screen. A directed
certificate must use the integer interface and retain every profile.

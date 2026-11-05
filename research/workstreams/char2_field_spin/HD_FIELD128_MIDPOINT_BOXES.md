# Grouping matching profiles without enumerating every pair

At one fixed active-group count, the exhaustive transfer calculation has
quadratically many matching profiles. We can instead partition those
profiles into rectangles. This note gives a cheap support-count bound
for each rectangle and states how that bound enters the transfer.
The construction and its setup distribution are unchanged.

Use the indexed binary supports and midpoint cap from
[the midpoint-count note](HD_FIELD128_MIDPOINT_COUNTS.md).
Let `P` be the number of pairs in each matching and `r` the support weight.
First complement supports when `r>P`; this replaces `r` by `2P-r` and
subtracts the original `r-P` from both full-pair counts.
Henceforth assume `r<=P`. Write

    L=floor(r/2), H=ceil(r/2), A=P-r,
    D=(P-L)L(P-H)H,
    U_a=C(P-L,H-a)C(L,a), V_b=C(P-H,L-b)C(H,b).

For `r>=2` and `a>=b`, the raw midpoint cap is

    C(a,b) = (P/D) (r+1-2a) U_a V_b
             * [LH(A+a+b)-Pab].

The formula is symmetric after exchanging a and b: `U_i/V_i` is
constant in i. It equals 1 for even r and `(P-L)/H` for odd r.
The exceptional profile `r=P,a=b=0` has count 2; its formula above is
zero, so we add 2 separately. The cases r=0 and r=1 are handled directly.
All feasible full-pair counts lie in `0..L` in these coordinates.

## Exact prefix sums of an upper bound

For `a>=b`, put

    E_a=(r+1-2a)U_a LH(A+a),
    F_a=(r+1-2a)U_a (LH-Pa).

Then the unscaled numerator is `E_a V_b+F_a bV_b`.
Prepare prefixes of V, bV, E, and F, and also prepare

    Q(x)=sum_(a<=x) [E_a sum_(b<=a)V_b + F_a sum_(b<=a)bV_b],
    J(x)=sum_(a<=x) V_a(E_a+aF_a).

Negative prefix indices have value zero. If `T(x,y)` denotes the
numerator sum over `0<=a<=x, 0<=b<=min(a,y)`, then

    T(x,y)=Q(x),                                      if x<=y;
    T(x,y)=Q(y)+Vprefix(y)(Eprefix(x)-Eprefix(y))
                +bVprefix(y)(Fprefix(x)-Fprefix(y)),  if x>y.

Inclusion-exclusion gives T on any rectangle intersected with `a>=b`.
Add the same calculation on the transposed rectangle and subtract J
on their common diagonal. Multiply the result by P and divide by D.
The division is exact because each summand is an integer flow-count cap.
Signed integer arithmetic is necessary: F can be negative.
Add the exceptional count 2 only when the rectangle contains `(0,0)`.

Preparation uses six linear-size prefix arrays. A query uses a constant
number of integer operations; their bit complexity depends on P.
No table indexed by both a and b is required.
This computes the sum of the raw midpoint caps, not the sum of their
pointwise minimum with another count bound.

## Directed transfer on a rectangle

In this section, r denotes the original active-group count, before
any complementation. The intervals I and J use the original full-count
coordinates. Only the count helper applies the complement shift;
regional matrices and the outer prefix remain at the original weight.
Fix this weight and the positive proof parameters. Let `R_a` be the
directed regional matrix after exact histogram normalization and
entrywise minimization over the available markers. For a rectangle
`I x J`, take componentwise maxima

    M_I=max_(a in I) R_a,  M_J=max_(b in J) R_b.

All matrix entries are nonnegative. If C is any valid upper bound on
the number of supports in this rectangle, its contribution is at most

    2^(constant+slope*r) C e_0 (M_I M_J)^128 1.

Use the smaller of the midpoint rectangle sum and any existing valid
support-count bound. The count is symmetric, but the matrix product
need not be: retain the ordered even/odd rectangles.

The implementation refines a freshly computed bound `2^u` from the
existing evaluator. If that evaluator used count `C_old>0`, replacing
it by `C_new<=C_old` gives the valid bound
`2^u C_new/C_old`. The ratio and logarithm use directed arithmetic.
A zero new count denotes an empty support event, even when the rectangle
contains marginally feasible profiles.

Only a disjoint, complete rectangle cover supplies a selected-weight
bound. Every integer weight must still be covered separately or by a
valid band argument. No interpolation between selected weights is used.

# A dimension obstruction and its cut interface

Fix the route and every quartet partition.  A state mask can certify a bad
codeword independently of the field maps.  The corresponding cut objective
has a compact graph representation, but its unconstrained scalar versions
are defeated by the empty and full masks.  This note records both facts;
it does not assert that a qualifying mask exists in the fixed construction.

There are `n=512` state coordinates, `T=4096` updates, and
`N=nT=2097152` positions.  There are `M=8192` outer RS `[256,128]` rows, and
the target bad-weight cutoff is `W=943718`.  Let `S_t` be any allowed state
mask, with `S_0=empty`.  Define

```
V_t(S) = S_t union Gamma_t(S_(t-1)),
L(S)   = sum_t |S_t|,
D(S)   = sum_t |S_(t-1) \ S_t|,
R(V)   = sum_outer_rows min(128, number of V positions in that row).
```

The terminal is free: there is no death charge after update `T`.

## A map-independent sufficient obstruction

If a mask satisfies

```
L(S) > R(V(S)),    R(V(S)) + D(S) + 1 <= W,                 (1)
```

then every choice of field maps has a nonzero codeword of output weight at
most `W` for this route and these partitions.

To prove this, restrict the state variables `Q` to the `L(S)` allowed
positions.  The induced input `X_t=Q_t-A_t Q_(t-1)` is supported in `V(S)`.
The outer parity-check equations restricted to these positions have rank
at most `R(V(S))`.  Their solution space therefore has dimension at least
`L(S)-R(V(S))>0`.  The causal maps from `Q` to `X`, and from `Q` to
`Y_t=Q_t+delta Q_(t-1)`, are injective.  Thus these solutions give a nonzero
linear output code of dimension at least `L-R`, supported on at most

```
sum_t |S_t union S_(t-1)| = L+D
```

positions.  Singleton gives an output of weight at most
`(L+D)-(L-R)+1=R+D+1`, proving (1).  The mask need only allow the state
support; the resulting word may use a proper subset of it.  The argument
does not assume that its nonzero amplitudes can be prescribed independently.

A qualifying mask would be an integer-checkable geometric obstruction to
the desired distance.  Conversely, absence of such masks would only exclude
this dimension argument.  It would not prove the distance bound: a special
amplitude assignment can still give a bad word when (1) fails.

## An exact directed-cut representation

For `eta,tau>=0`, the objective `R(V(S))+eta D(S)-tau L(S)` is represented
by a directed minimum cut.  A vertex on the source side means that its
binary variable is one.  Use state vertices `s_(t,i)`, footprint vertices
`v_(t,i)`, one vertex for each quartet at updates `2,...,T`, and one auxiliary
vertex for each outer row.  Add these arcs:

* source to every state vertex, capacity `tau`;
* each state vertex to its same-position footprint vertex, infinite capacity;
* each state vertex at time `t-1` to its quartet vertex at time `t`, infinite
  capacity;
* each quartet vertex to its four footprint vertices, infinite capacity;
* each footprint vertex to its routed outer-row auxiliary, capacity one;
* each outer-row auxiliary to the sink, capacity 128;
* each state vertex at time `t` to the same coordinate at time `t+1`, capacity
  `eta`, for `t<T`.

For fixed `S`, minimizing over a row auxiliary costs the smaller of its
covered-position count and 128.  The implication arcs force the canonical
footprint; adding further footprint positions cannot reduce its minimized
cost.  Consequently the minimum cut equals

```
tau*N + min_S [R(V(S))+eta D(S)-tau L(S)].                 (2)
```

Rational capacities can be scaled to integers, and an infinite capacity can
be replaced by an integer larger than a known finite feasible-cut cost.
The graph has about 4.73 million vertices and 12.59 million forward arcs.
No flow implementation or numerical search is included here.

## Why unconstrained scalarization cannot find the relevant masks

The empty mask has `(L,R,D)=(0,0,0)`.  The full mask has
`(L,R,D)=(N,N/2,0)`: every outer row has 256 covered positions, and there is
no final death charge.  For any target mask with `0<=L<=N`, the convex
combination of these two points with full-mask weight `L/N` has the same
`L` and weighted cost `R+eta D=L/2`.

In particular, whenever the target satisfies `R+eta D>L/2`,

```
min(0, N/2-tau*N)
  <= (L/N)*(N/2-tau*N)
  < R+eta*D-tau*L.                                      (3)
```

Thus the target cannot minimize (2), for any `tau`.  This applies to the
profiles under discussion, where `R` is close to `L`.  A finer scalar grid
cannot repair the obstruction.

The issue also rules out a Lagrangian certificate based only on the two
inequalities in (1).  The empty/full convex hull contains the point with
full-mask weight `2(W-1)/N`, whose coordinates are

```
L=2(W-1),    R=W-1,    D=0.
```

This fractional mixture satisfies both `L-R>=1` and `R+D<=W-1`.
Therefore no nonnegative linear combination of those two resource
constraints, minimized by the unconstrained cut oracle, can separate all
cuts from their joint feasible region.  The mixture is not asserted to be
an actual state mask.

A useful obstruction search must retain a genuine volume/rank constraint
or pinned geometry that removes this empty/full relaxation.  Scalar
minimum cuts then provide bounds for the restricted search, and any
returned mask can be checked directly against (1).  Such restrictions
would need their own coverage argument before they could support a global
theorem.  The next proof step remains the output-aware quotient with
closed-run or actual-death information, rather than an unconstrained
minimum-cut parameter search.

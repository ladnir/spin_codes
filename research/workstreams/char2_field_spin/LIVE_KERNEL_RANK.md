# Retaining the joint rank in a packet bound

A packet's state and active inputs can depend on few source forms, even
when the packet has many live coordinates. A bound that uses only the
numerical rank of its equations can lose this information. The following
lemma retains the dimension of the joint form space.

This is a local rank-potential bound, not a whole-code distance theorem.
It applies to the joint packet equations before their diagonal randomness
is exposed. It must not be applied after conditioning on some of those
same equations.

## Fixed interface and remaining randomness

Let F be a field of size q. Fix a subspace U of F^m with dimension d and
a matrix B over F with m columns and rank r. Assume 1<=d<=m<q. Sample
the m diagonal entries of Delta independently from F minus zero. Fix a
real mark tau in [1,q-1].

For each Delta, let K_Delta be the kernel of B Delta restricted to U,
and write j_Delta for its dimension. Let E_Delta mean that K_Delta
contains a vector with all m coordinates nonzero. The quantity to bound is

    E[ tau^(j_Delta-d) 1_{E_Delta} ].

If U is the image of a map from a larger source quotient, d-j_Delta is
still the rank gained by imposing these equations. The hidden kernel of
that source map does not change the rank gain.

## Bound

For j>=1, put

    L_j = q^j - 1 - m(q^(j-1) - 1).

These numbers are positive because m<q. If r=m, the quantity above is
zero. Otherwise it is at most

    min(1, (q^d-1)/(q-1)^r *
            max_{1<=j<=min(d,m-r)} tau^(j-d)/L_j).

A simpler, weaker version replaces the maximum by
`tau^(1-d)/(q-m)`.

To prove the bound, fix a full-support vector v in U. The vector Delta v
is uniform on (F minus zero)^m. After fixing the nonpivot coordinates,
the rank-r system B Delta v=0 determines its r pivot coordinates.
Therefore its probability is at most (q-1)^(-r). Summing over the
full-support vectors of U bounds the expected number of full-support
vectors in K_Delta by (q^d-1)/(q-1)^r.

Now suppose E_Delta holds and j_Delta=j. Every coordinate functional is
nonzero on K_Delta. Its zero set is a hyperplane with q^(j-1) vectors.
All m coordinate hyperplanes contain zero, so their union has size at
most 1+m(q^(j-1)-1). Thus K_Delta has at least L_j full-support vectors.
Divide tau^(j-d) by this lower bound, maximize over possible j, and use
the expectation above. Since Delta is invertible, j<=m-r. Also j<=d.

The simpler version uses at least q^(j-1)(q-m) full-support vectors and
the inequality (tau/q)^(j-1)<=1.

## Upper bounds on the dimensions

For fixed j, `(q^d-1) tau^(j-d)` increases with d when tau<=q-1.
Increasing d also enlarges the allowed range of j. Increasing m decreases
L_j and enlarges that range. Consequently upper bounds on d and m can be
substituted in the displayed formula.

For a packet with s state coordinates, incoming state-form rank at most
c, and w active input coordinates, first impose a canonical incoming
state-zero flat of rank k. On a feasible branch,

    d <= c-k+w,       m <= s-k+w.

The flat's rank costs tau^(-k), once. Apply the new bound to the full
return equation matrix, or to the return and selected-output equations
jointly. Their current input and state diagonals must remain independent
and unexposed at this point. A separately proved bound can be combined
with this one by taking their minimum.

For example, a width-eight/state-eight return with c=1,w=1,k=0 has
d<=2,m<=9,r=8. The refined factor is

    (256^2-1) / (255^9 tau).

At tau=6 this is much smaller than tau^(-8). This local improvement does
not imply a comparable global improvement: dominant paths may already
have saturated the state-form rank.

## Independent checks and scope

Two independent reviews checked the quotient interpretation, chronology,
rank gain, and monotonicity argument. The exact helper is
`theory/live_kernel_rank.py`.

`theory/live_kernel_rank_tests.py` checks every source subspace and every
equation row space in GF4^3, every nonzero diagonal, and three marks.
All 5547 moment comparisons pass. Additional tests cover monotonicity in
d and m, the singleton return formula, and the full-column-rank case.

This lemma counts local full-support vectors to bound a rank potential.
It does not restart a union bound over all messages of the global code.
Global source and routing accounting, class coverage, and outward
arithmetic remain separate requirements for a distance certificate.

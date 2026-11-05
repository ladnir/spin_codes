# Retain the quiet constraints in the live-state count

The [column-scaled inner](QUIET_BRANCH.md) has two useful partial proofs:
small source supports and trajectories that never enter a live state.
This note supplies a valid interface between quiet constraints and the
live-state proof. It also identifies a dense class where that forward
interface is insufficient, even with ideal rank. No construction changes
or new whole-code certificate are claimed.

## A common message space

Keep the outer, route, A, B, D, and R of QUIET_BRANCH.md. Write q=2^128
and e=256. First fix the route, a specified union V of live episodes,
and a disjoint set J of nonzero quiet batches. An episode includes its
activation batch. Let U be all coordinates in V and J. Input is zero
outside U. Fix these sets before exposing R or D; a final union must
cover every possible set.

Let C_U be an E-generator matrix, with N coordinate rows, for the routed
outer shortened to U. Its number of columns is

    d = sum_i d_i,    d_i=max(0,t_i-128),

where t_i is the number of coordinates of outer block i in U. Form L_J
by applying B separately to the coordinates in each batch of J and
retaining only those batches. Quiet input constraints give the matrix

    H_J(R) = L_J diag(R) C_U.

The remaining message space has dimension d_R=d-rank(H_J(R)). This
dimension is identical over E and F. It is a space containing all inputs
of the proposed mixed trajectory, not a claim that every element realizes
the trajectory. In particular, the live-state constraints are still absent.

## Rank from a common independent coordinate set

Choose a set S of r quiet coordinates such that L_J[:,S] has column rank r
and C_U[S,:] has row rank r. The choice depends on the fixed route and
sets, not on R. Choose r rows and r columns that make the corresponding
two square submatrices invertible. Condition on every R entry outside S.
Multiplying this r-by-r minor of H_J(R) by the two inverse matrices gives

    diag(R_S) + C,

where C is fixed and independent of R_S. These are invertible changes of
coordinates, not an assumption that the entire matrix is random dense.

Expose the diagonal entries in order. If the current diagonal pivot is
nonzero, eliminate its row and column. If it is zero, discard that row
and column without claiming a rank gain. The latter event has conditional
probability at most 1/(e-1), since the pivot is its fresh nonzero R entry
plus a value determined by earlier entries. A Schur complement preserves
this form for each unexposed diagonal. Discarding a failed pivot may lose
useful off-diagonal rank; that only weakens the bound.

After r steps, the rank deficit of the chosen minor is at most the number
of failed pivots. Conditional failure probabilities bound its upper tails
by those of Z distributed as Bin(r,1/(e-1)). Thus

    d_R is stochastically bounded above by d-r+Z,
    E_R[tau^d_R] <= tau^(d-r) (1+(tau-1)/(e-1))^r, tau>=1.  (1)

The expectation includes the earlier conditioning on entries outside S.
The bound holds uniformly in those entries. It does not assert that
quiet checks are independent or that the realized matrix has rank r.

## Determine r from the actual outer and batch structure

The shortened RS block on its t_i retained coordinates is [t_i,d_i] MDS.
Hence any at most d_i rows of its generator are independent. The matrix
L_J is block diagonal, and any at most 64 columns in one batch are independent.
Consequently a common independent set is exactly a set of quiet coordinates
using at most d_i from outer block i and at most 64 from each quiet batch.

Build a bipartite graph between outer blocks and quiet batches. Each retained
coordinate is one edge. Give source vertex i capacity d_i and every batch
vertex capacity 64. An integral maximum flow gives the largest possible r.
Parallel edges are retained through integer capacities. This is a small-field
matrix argument for the actual route, not a replacement route distribution.

In particular, when U consists entirely of quiet batches, assigning flow
d_i/t_i to each edge of block i saturates every source capacity. Since
d_i/t_i<=1/2, no 128-coordinate batch exceeds capacity 64. Integrality then
gives r=d. This is generic rank, not full rank for every sampled R.

## Combine with the live-state probability, in the correct order

Let v=|V|, j=|J|, and c be the number of episodes. Let a count their
prescribed resets; an episode reaching the final batch has no prescribed
terminal reset. Nonzero quiet batches contribute at least 65j output symbols.
For output threshold W, the retained live-state proof therefore applies with

    m = ceil(max(0,128v-191c-128a+65j-W)/2),
    T = 64a+(107/128)m.

Fix all R first. The proof over the independent D setup bounds each fixed
input's joint trajectory/low-weight event by 2q^-T. Inputs that do not
realize the specified quiet nonzero pattern contribute zero. Counting
projective representatives in the common kernel gives the bound

    min(1,4q^(d_R-1-T)).

For 1<tau<=q, apply the fractional power log_q(tau), then average R.
Together with (1), the conditional probability for this fixed route and
fixed trajectory pattern is at most

    tau^(1/64-1-T+d-r) (1+(tau-1)/255)^r.                  (2)

A final bound must average the actual route and sum all trajectory patterns.
Equation (2) is not such a certificate. It does, however, preserve the common
message space, conditions in the right order, and permits setup-adaptive
messages. It does not multiply independently averaged partial certificates.

## A verified obstruction to using this forward bound alone

Consider v=8192, j=8192, c=1, a=0, and U equal to the complete coordinate
set. The single episode can occupy the final half of the scan; quiet
batches precede it. This is a class allowed by the counting interface,
not an exhibited codeword. Here d=K=1048576. There are only 64j quiet
equations, so even ideal full row rank leaves at least 524288 dimensions.

At W=floor(.491N), the defect witness is m=275582 and T=14743637/64.
Even replacing all rank randomness by perfect rank, and omitting every
route and pattern-counting penalty, the exponent in (2) is at least

    524288-1-T+1/64 = 4702683/16 > 0.

Thus this forward envelope cannot close that class at any tau>1. Better
quiet rank, larger small-field coefficients, or tighter route averaging
alone cannot repair this particular deficit. It is a failure of the
bound, not an upper bound on the code's distance.
The subsequent [output-interval proof](OUTPUT_INTERVAL.md) closes this
example and the surrounding dense one-episode class without changing the
encoder. The failed forward envelope remains useful evidence about which
counting direction to use, not an unresolved instance of that class.

The next direction is output-side enumeration for dense mixed classes.
Quiet output batches impose 64 equations directly on their output variables.
Expanded-state holes also constrain outputs, and zero-output runs preserve
the state. Those constraints may reduce the output count much more than
the forward argument reduces its message count. Retain (2) for sparse
allowed-coordinate sets where the shortened outer dimension is already small.

## Exact checks

[test_quiet_rank.py](test_quiet_rank.py) checks the following over GF4:

- Every 2-by-2 constant matrix and 64 sampled 3-by-3 constant matrices,
  each under every nonzero random diagonal: 704 exact tail comparisons.
- 64 sampled rectangular products L diag(R) C under every diagonal:
  5184 exact matrix ranks, with rational moment comparisons at 3/2, 4, and 16.
- All 16 routes and all 256 coordinate restrictions of two RS[4,2] blocks:
  4096 capacity-flow comparisons against exhaustive coordinate-set selection.
- The K20 forward-envelope obstruction above, using exact rational arithmetic.

All pass. The first family exhausts the specified constants; the second
samples constants but exhausts their diagonals. These are tests of local
interfaces and a bound diagnosis, not large-code distance experiments.
The flow helper is for research-scale tests, not an encoding hot path.
Run serially, with no generated files:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_quiet_rank.py
```

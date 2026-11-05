# All sparse outer inputs through 31 blocks

For the unscaled [post-diagonal construction](POST_DIAGONAL.md), at K=2^20,
the probability that any nonzero message supported on at most31 outer blocks
has output weight at most floor(.491N) is less than2^-50. The probability
is over the sampled diagonal setup. The bound holds for every realization
of the specified regional route, so it also holds after averaging the route.
It includes messages chosen after setup. Inputs using32 or more outer blocks
are not covered by this theorem.

## Fixed parameters and source classes

Write q=2^128, N=2^21, M=N/256=8192 outer blocks, n=N/128=16384 batches,
and g=M/128=64 batches per region. The state dimension is s=64.
The maps A and B are those defined in DESIGN.md; no row scaling is used.

Fix a route and then a nonzero input before exposing the diagonals. Let p
be its number of nonzero outer blocks, and let h be the number of regions
containing at least one nonzero routed input. Every nonzero RS[256,128]
word has weight at least129, so 129<=h<=256.

For specified sets of p blocks and h regions, all such outer words lie in
a space of dimension d=p(h-128): each block is shortened to those h positions.
The number of one-dimensional subspaces is less than2q^(d-1).
Scalar multiples have identical output supports, state supports, and defect
events for each setup. Hence projective counting is valid. The total number
of representatives needed for this (p,h) class is at most

    M^p * 2^256 * 2q^(d-1) = 2^A,
    A = 13p+257+128(d-1).

This overcounts both block/region choices and words with smaller supports.
It does not depend on the route realization.

## Three exceptional events for one fixed input

A live batch has nonzero incoming state. A reset changes a live state to
zero. An activation changes an empty state to nonzero. For p<=64, a nonzero
input batch has weight at most p: each selected outer block contributes
at most one symbol to a region, and batches do not cross regional boundaries.
Since ker(B) has minimum distance65, every nonzero input batch activates
an empty state. This restriction is essential to the interval argument.

The exact refresh lemma gives uniform Q_t in F^64 whenever Q_(t-1)!=0,
conditional on all preceding diagonals. In particular, a prescribed reset
at such a step has probability q^-64. Iterating this conditional bound and
union bounding over steps gives

    Pr[at least r+1 resets] <= n^(r+1) q^(-64(r+1)).       (1)

Define E as the total number of zero coordinates of A Q_t, summed only over
steps with both Q_(t-1)!=0 and Q_t!=0. These are holes in nonzero refreshed
states, not holes after activations. Such a state has at most63 holes.
For any prescribed ell<=63 coordinates, the corresponding rows of A have
rank ell. Conditional uniformity therefore gives probability at most q^-ell
that these coordinates vanish and Q_t remains nonzero.

If E>=m, select exactly m actual zero coordinates across these refreshed
states. There are at most N^m choices. For a fixed choice, expose the selected
steps chronologically. At each selected step, incoming liveness is already
determined; if absent, the event is impossible. Otherwise the conditional
probability is at most q^(-ell) for that step's selected coordinates.
No selected step needs more than63 coordinates. Consequently

    Pr[E>=m] <= N^m q^-m.                                (2)

This charges total holes, not63 times the number of defective batches.

Define Z as the number of zero output coordinates at positions where
(A Q_(t-1))[j]!=0. Such a position has one fresh uniform coefficient D_t[j],
so its output is zero with probability1/q, conditional on preceding choices.
The set of eligible positions is known before exposing that batch's diagonal.
Exposing entries within a batch one at a time and selecting m witnesses gives

    Pr[Z>=m] <= N^m q^-m.                                (3)

The three events need not be independent. We only take their union.

## From few exceptions to many nonzero outputs

Suppose there are at most r resets and E,Z<m. There are at most r+1
activations. A state created by activation may not be uniform, but its
expansion has at most63 holes. Every other incoming nonzero state arose
from a live predecessor and its holes are counted in E.

Let L be the total number of live batches. Counting only their outputs gives

    wt(Y) >=128L -63(r+1) -2(m-1).                        (4)

We next lower-bound L without assuming anything about positions inside regions.
An episode starts with an activation batch and includes every following live
batch until a reset, or until the end. If there are c episodes, their total
length including activations is L+c. Every nonzero input batch belongs to
one of these episodes, because an empty state cannot ignore a nonzero input.

An interval of v batches meets at most v/g+2 regions. If it reaches the end
of the scan, it meets at most v/g+1 regions. If the last episode is closed,
c is at most r, so h<= (L+c)/g+2c. If the last episode is open, c<=r+1
and h<=(L+c)/g+2c-1. Both cases imply the safe common bound

    L >= max(0, g(h-2r-1)-(r+1)).                        (5)

For r=0, use the sharper L>=g(h-1). The first occupied region has index at
most256-h, and its activation batch occurs before the end of that region.
With no reset, all subsequent batches are live through the end of the scan.
No terminal-state constraint is imposed in either argument.

## Exact finite coverage

For each p,h, the verifier chooses

    r = ceil((A+64)/8178)-1,
    m = ceil((A+64)/107).

Here 8178=64*128-log2(n), and107=128-log2(N). Multiplying each of (1)--(3)
by the projective count2^A gives an upper bound at most2^-64 for each event.
Equations (4)--(5) then give a deterministic weight lower bound outside
their union. The verifier checks all128 values of h for every p=1,...,31.

All3968 cells have weight lower bound at least1,039,505, exceeding both
floor(.49N)=1,027,604 and floor(.491N)=1,029,701. There are three exceptional
events per cell, and 3*3968<2^14. The total failure probability is less than
2^14*2^-64=2^-50. The smallest weight bound occurs at p=31,h=132:
d=124, r=2, m=154, L>=8125.

Everything in this finite replay is integer arithmetic; no numerical
logarithm or floating comparison is used. The source proof supplies the
probability inequalities, while the script checks complete coverage and
their arithmetic. The small-timeline test separately enumerates possible
activation/reset histories and checks (5), including both final-state cases.

## Boundary and next refinement

At p=32,h=256 this envelope gives d=4096, r=64, m=4906, L>=8063, and
weight lower bound1,018,159, below the target. That is a failed lower bound,
not a low-weight codeword. Equation (5) lets each reset hide almost two
whole regions, regardless of how many active coordinates those regions
contain. Actual random regional positions should restrict those gaps.

The next refinement should retain regional input occupancy or use shortened
dimensions on unions of inactive intervals. Beyond p=64, nonzero batches
can also lie in ker(B), so activation itself needs a different argument.
Neither issue can be bypassed by extrapolating this sparse theorem.

Reproduce serially with Python -B:

    sparse_certificate.py --numerator 491 --denominator 1000
    test_sparse_certificate.py

The scripts write no data files and need no external numerical libraries.

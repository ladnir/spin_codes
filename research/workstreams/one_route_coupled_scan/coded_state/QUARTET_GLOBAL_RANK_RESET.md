# Global proof reset: count output masks and retained RS rank

2026-10-31. This note concerns the current field-quartet construction, not
the older width-four recurrence in the root CONSTRUCTION.md. The encoder,
outer, routing distribution, and setup acceptance rule are unchanged.
The whole45%/40-bit certificate is still open. The next
[global low-rank checkpoint](QUARTET_GLOBAL_LOW_RANK.md) groups optional
outputs before taking moments and certifies every history with
|V|>=53N/64, with total failure below 2^-6359.

The recent regional calculation refined a sufficient geometry event after
discarding several field constraints. Improving that event is not necessarily
the best route to minimum distance. Here we return to a direct count of
code outputs. We keep the RS parity rank and state-death equations together,
without first choosing an active outer-row set or maximizing a support score.

Three independent reviews agreed on the joint-rank argument below. The
bounded checks are in `quartet_inverse_rank_audit.py`. The rank and support
lemmas, not the small-instance tests, justify the general statements.

## Current construction and probability space

Work over J=GF(2^32), with q=2^32, K=2^20, N=2^21, and M=8192 outer rows.
Each row uses RS[256,128], extended from its subfield to J. Its parity-check
matrix has rank min(128,k) on every set of k columns. Coordinate j from
each row goes into region j; the M coordinates within each region are
independently uniformly permuted. There are256 regions.

The inner has512 state coordinates and4096 updates. Each update partitions
the state uniformly into quartets. On each quartet, A_t multiplies by an
independent raw uniform element of GF(q^4), in a fixed J-basis. Setup rejects
unless both local matrices A_t and A_t+delta I are superregular. The existing
global acceptance-density bound is gamma<=1024/989<2. All setup is sampled
once and fixed for every message. Delta is a fixed nonzero subfield scalar.

The recurrence and its inverse are

    Q_0=0,  Q_t=A_t Q_(t-1)+X_t,  Y_t=Q_t+delta Q_(t-1),
    Q_t=Y_t+delta Q_(t-1),  X_t=Q_t+A_t Q_(t-1).

Thus fixing an output vector Y fixes Q independently of every random map.
This is an enumeration device: it does not assume an adversarial message
is independent of setup. The terminal state remains free.

## Joint parity and death rank

Fix the route, quartet partitions, an output support J0 of size j, and an
exact state-support history S compatible with J0. Let d count transitions
from nonzero state to zero state, excluding any hypothetical terminal reset.
Let L contain all input coordinates in quartets with nonzero incoming state.
Put B=J0 outside L and V=L union B. Outside L the entire incoming quartet
is zero, so X=Y and B also equals the new state support outside L. Therefore
V depends only on S and the partitions, not on optional output zeros at
surviving state coordinates.

In outer row i, write ell_i=|L intersect row i| and b_i=|B intersect row i|.
Define

    R_L=sum_i min(128,ell_i),
    R_B=sum_i [min(128,ell_i+b_i)-min(128,ell_i)],
    R(V)=R_L+R_B=sum_i min(128,|V intersect row i|).

For each fixed candidate Y, the raw map on a live quartet makes its four
X coordinates jointly uniform over J^4: c->c Q_prev is a bijection when
Q_prev is nonzero. Different live quartets use independent coefficients.
If the fixed X outside L satisfies the residual parity equations, RS membership has
probability exactly q^(-R_L).

Compatibility of the fixed outside coordinates gives R_B independent
linear equations on the Y variables in B. Separately, each state death
gives one zero-state equation. Use the excursion interval ending at that
death; its last output variable is nonzero and lies in L. These intervals
give d distinct live pivot columns. No nonzero combination of their
equations is supported only on B, whereas every residual RS equation is.
Consequently the two systems have joint rank d+R_B.

The number of candidate output amplitudes is at most q^(j-d-R_B), after
relaxing the required nonzero values and any other state constraints.
For accepted setup, this proves the conditional bound

    Pr_field[there is an output with masks (S,J0)]
        <= min(1, 2 q^(j-d-R(V))).                         (1)

This counts all messages in the mask class, not one fixed message.

There is also a necessary compatibility exclusion. If some row has
b_i>0 and ell_i+b_i<=128, its residual equations have full rank b_i and
force those declared-nonzero outside outputs to zero. Such a mask is
impossible. This exclusion is important for late activations; live rank
R_L alone does not suffice. If every row has |V intersect row|<=128, the
shortened outer space is zero, so the entire nonzero mask is impossible.

## A global support census, without pending-run types

An output support J0 and its subset of death positions uniquely determine
the state-support history: a nonzero output activates an empty coordinate;
an active coordinate remains active except at a declared death. Zero
outputs preserve whether that coordinate is active. Infeasible choices
only reduce the number of histories.

Hence the number of masks with output weight j and d deaths is at most

    binom(N,j) binom(j,d).                                (2)

There is no union over outer-row subsets. In particular, inactive-row
count is not a parameter of this census.

Let ceilLog2(m) be the least integer c with m<=2^c, and define

    g(j,d)=ceil((ceilLog2 binom(N,j)
                 +ceilLog2 binom(j,d)+108)/32).

For any fixed route and partition geometry, the probability that any
nonzero output with 1<=j<=943718 satisfies

    R(V)+d >= j+g(j,d)                                    (3)

is less than2^-63. Indeed, (1) and (2) give at most2^-107 for each (j,d).
There are at most(N+1)^2<2^44 profiles. This includes every output amplitude,
every support history in (3), and every free terminal state. Averaging the
geometry preserves the same bound. No numerical optimization is required.

This is a broad class theorem, not a minimum-distance certificate: a bad
word outside this class may have rank below (3).
The class theorem is stated over J. A complete minimum-distance theorem
would transfer to the128-bit extension by taking a nonzero J-basis
component; a restricted rank-class theorem does not automatically transfer
with the same class label.

## Simple concrete corollaries

The support-mask matrix for a single state coordinate, retaining q^-d, is

    [[1,1],[q^-1,2]].

The two live-to-live choices distinguish zero from nonzero output. Both
row sums are at most2+q^-1. Starting every coordinate at zero and accepting
every final state gives a total weighted mask count at most(2+q^-1)^N.
Therefore all bad outputs with

    R(V)>=1009256 = 943718+65538

have total probability below2^-62, conditional on any geometry. This
permits rank deficiency K-R(V) up to39320, not just full parity rank.

Leaving q^(-d/2) in the mask sum instead proves a second class bound:

    R(V)+d/2>=1009257

has total failure below2^-71. Here the correction to2^N is at most
(1+2^-17)^N<e^16<(11/4)^16<2^24. The exact factor-two calculation yields
2^-71; an earlier discussion's2^-72 rounding was too strong.

These corollaries need not be combined with (3). If their certified
classes are combined, include their failure allowances explicitly.

## The remaining global geometry problem

To close the original target, bound all compatible nonzero masks with

    R(V)+d < j+g(j,d),  1<=j<=943718.

The task is to show that covered coordinates leave too little RS parity
rank only with sufficiently small total probability after the support
census. It is not enough to show that a typical fixed support has high
rank, or to assume a sampled code inherits a fixed-message probability.
The incompatible short rows and zero shortened outer space must remain
excluded. Keeping those constraints avoids paying for some small-support
geometry minimizers that no nonzero outer message can realize.

A deterministic lower bound is useful for concentrated histories. If v_r
counts covered coordinates in region r, then

    R(V)>=sum of the128 largest v_r.

Choose those128 coordinate regions. Every row has at most128 selected
columns, so all selected parity-check columns are independent. This
handles short concentrated bursts more sharply than R(V)>=|V|/2, but
does not settle general distributed support histories.

## A quota-respecting route census

Fix S and the quartet partitions before exposing the route. In region r,
the v_r covered positions map to a uniform v_r-subset of the M rows.
Let P=K-R(V)=sum_i(128-v_i)_+ be the resulting rank deficiency.
Here v_i denotes a row count, distinct from regional v_r.

Set p_r=v_r/M. Start instead with independent Bernoulli incidences for
each row and region, then condition region r to contain exactly v_r rows.
This produces the actual route law on covered sets. Its conditioning
probability is pi=product_r Pr[Bin(M,p_r)=v_r]. For any t>=1 and rho>0,
positivity and independent rows give

    E_actual t^P
      <= pi^-1 rho^(-|V|)
         [E(t^(128-Z)_+ rho^Z)]^M,                        (4)

where Z is the sum of independent Bernoulli(p_r) variables. The rho marker
keeps the deterministic total coverage; omitting it can allow impossible
occupancy censuses. Equation (4) retains each exact regional quota through
conditioning, unlike the previously rejected global inactive-count model.

Each regional count is a binomial mode, so pi>=(M+1)^(-256). Also the
function k->t^(128-k)_+ rho^k is convex: it is the maximum of two positive
exponentials. Pairwise equalizing two Bernoulli probabilities, with their
sum fixed, increases its expectation by a nonnegative second difference.
Thus replacing Z by Bin(256,|V|/N) gives a profile-uniform upper bound in
(4). This optional relaxation can be loose. Keep the actual profile for
concentrated cases, or use the deterministic top128-region bound.

Equation (4) is a new candidate tool for the low-rank remainder, not an
evaluated ensemble certificate. A fractional moment must be applied to
the conditional per-mask count before summing masks; no unjustified
interchange of a fractional power with the whole sum is permitted.

## Checks, scope, and next action

The bounded replay enumerates all65535 nonzero GF4 output vectors in a
two-row RS[4,2], two-update model. It compares624 exact support-mask
classes with the joint-rank bound and rejects255 impossible short-row
outputs. It also checks1397 exact output/death support counts,36 exact
conditioned-census cases, and the three integer failure-budget inequalities.
The small model uses the raw live-quartet point law; it is not presented
as an accepted GF4 instance of the field-quartet construction.

The construction and earlier checkpoints are preserved. No performance
work or new sampled-code claim is involved. The previous seven-type work
is not erased, but it is no longer the principal next step.

Next evaluate a deterministic global bound for the low-rank remainder
using (1)--(4), together with forced-birth compatibility. The footprint
law can be counted by incoming weight, outgoing weight, deaths, and the
number of active quartets; it does not require canonical run ages. If
that bound loses too much, diagnose the specific rank/census relaxation
before returning to local population refinements or changing the code.

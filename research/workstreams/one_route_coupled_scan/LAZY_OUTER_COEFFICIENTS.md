# Use the outer coefficients to charge later input matches

The first occupied columns supply fresh source coordinates. Later columns
do not: their values are linear combinations of source coordinates already
consumed. The current inner proof then discards some input/state matches,
even though the outer still has unused independent randomness.

This note exposes that randomness only when it is needed. It gives local
lemmas and a finite interface for two adjacent groups on the same complete
columns. It changes no encoder operation and asserts no numerical closure.

## 1. Source class and exposure experiment

Use the [specified construction](CONSTRUCTION.md), with E=GF(q), q=256,
and the degree-16 extension F. Fix the shared set H of h occupied RS
columns in two adjacent groups, where 33<=h<=64. Every selected column
has all four mixed coordinates nonzero in both groups. All other source
coordinates are zero. Put b=h-32.

Fix the full route before revealing the active columns' coefficients.
Fix the absent-column mixers to arbitrary invertible realizations. Their
zero output constraints are equivalent to zero pre-mixer RS columns.
The shortened message space V therefore has dimension 8b over F.

In each group, sample the first b occupied columns using the actual
distribution of invertible Toeplitz matrices. Expose these matrices before
the inner proof. Their mixed coordinates form an information set for that
group's shortened space. This preserves the existing
[fresh-source prefix](TWO_GROUP_ROUTE.md#8-fresh-source-coordinates-in-the-first-occupied-columns).

For the remaining h-b=32 columns in each group, first sample all seven
Toeplitz diagonals independently from E, without conditioning on
invertibility. Reveal them lazily below. Call this the raw experiment.
All inner coefficients retain their original independent distributions.

Let G be the event that all 64 raw matrices are invertible. Its probability
is (1-1/q)^64. Conditioning the raw experiment on G gives the actual setup
distribution for these columns. Thus, for any event B in this source class,

    Pr_actual[B] = Pr_raw[B and G]/Pr_raw[G]
                 <= (q/(q-1))^64 Pr_raw[B].                (1)

The factor is per matrix, not per output component or coefficient charge.
The established probability that a square Toeplitz matrix is invertible
is used here; see the
[column-mixer lemma](../char2_field_spin/HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers).
For q=256 the factor in (1) is less than 1.285.

The first b selected matrices are not sampled raw. Applying source
freshness to a potentially singular prefix would be invalid. All source
forms in that prefix remain independent because its matrices are sampled
from the actual invertible distribution.

Keep one common space W of zero equations and all exposed nonzero
requirements throughout the proof. The global potential is

    Phi=1[branch retained] tau^(-dim W),       tau>=1.

Future mixed input forms can depend on unrevealed raw diagonals. Check
their required nonzeroness when those inputs are exposed; do not condition
the coefficient distribution on future source supports. The complete
event still requires every selected mixed coordinate to be nonzero.

## 2. Partition a pre-mixer vector by its first nonzero coordinate

Immediately before processing a raw column, its pre-mixer RS vector

    v=(v_0,v_1,v_2,v_3)

is a known vector of linear forms on V. It does not depend on that
column's Toeplitz diagonals. A witness for the selected source class has
v!=0, because a zero vector could not produce four nonzero mixed outputs.

Partition such witnesses by their first nonzero coordinate ell. On branch
ell, impose v_0=...=v_(ell-1)=0 in the same W and require v_ell!=0.
Discard a branch if its new zero equations kill any previous nonzero
requirement. Retain the requirement v_ell!=0 for all later steps.
The all-zero branch cannot contain a witness and is omitted.

The partition need not cost four independent copies of the potential.
Let W_ell=W+span(v_0,...,v_(ell-1)). A branch can survive only if
v_ell is outside W_ell. As ell increases, these possible pivots occur at
successive prefix ranks 0,1,...,d-1, where d<=4. Coordinates already in
the preceding span admit no surviving nonzero-pivot branch.
Therefore the sum of relative rank potentials is at most

    g_tau = 1+tau^(-1)+tau^(-2)+tau^(-3).                  (2)

The exponent is the actual prefix-rank gain, not the coordinate index ell.
No independence among v_i or between v_i and prior constraints is assumed.

### Retain the change in residual state rank

Let v_live(r,S) be the verified positive live-state potential used by the
empty-gap envelope. Here S is the current state-support mask, with
1<=r<=|S|. It is distinct from the pre-mixer vector v above.

A pivot branch that adds i independent zero forms can lower the residual
state rank from r to a value s with

    max(1,r-i)<=s<=r.

The retained nonzero requirements keep its state-support mask S unchanged.
If the state rank became zero, that live branch would be rejected.
The source partition therefore has a common two-state diagonal envelope

    P_source=diag(g_tau,kappa),

    kappa=max_(r,S) sum_(i=0)^3 tau^(-i)
              max_(max(1,r-i)<=s<=r) v_live(s,S)/v_live(r,S). (3)

On a quiet branch, the state remains quiet and (2) applies directly.
Equation (3) permits arbitrary correlations with the incoming state and
earlier constraints. Apply P_source once for each group's raw column,
before its first scalar region. Do not apply it once per region or event.
For the present two-group class, that is two source partitions per raw
occupied column, and 64 source partitions in total.

## 3. A new Toeplitz coefficient at each component

Index a raw matrix by T_(a,i)=d_(a-i), with a,i in {0,1,2,3} and seven
independent diagonals d_-3,...,d_3. On the pivot branch ell,

    (T v)_a = d_(a-ell) v_ell
              + sum_(i>ell) d_(a-i) v_i modulo W.         (4)

Process components a=0,1,2,3 in their actual scalar-region order. Before
component a, reveal any still-unrevealed diagonals with index below a-ell.
The diagonal d_(a-ell) is fresh: earlier components used only smaller
indices after the zero forms v_i, i<ell, were imposed.

The coefficient multiplies the retained nonzero form v_ell. For any known
form u, a required zero equation

    u+d_(a-ell) v_ell=0

has the usual rank/root alternatives in the current common quotient:

- If u is outside W+F v_ell, the equation adds one rank for every diagonal.
- If u=c v_ell modulo W, at most one of the q diagonal values can survive
  without forcing v_ell into W.

Thus the expected relative potential contracts by at most
1/psi, where psi=min(tau,q). The diagonal may be zero; this causes no
problem for the affine coefficient argument. Other required nonzero
forms can only discard more branches.

This is an argument in the raw experiment. Do not condition an individual
diagonal on its matrix being invertible. Use the single final factor (1).

## 4. A later active-step charge

Fix an inner step in a raw column's scalar region. Each active coordinate
comes from a distinct outer group, hence has its own independent fresh
diagonal from (4). Let I be its exact input-support mask, S its incoming
state-support mask, J its output-support mask, and S' its next-state mask.
The standard necessary conditions remain

    I symmetric_difference S subset J subset I union S,
    J symmetric_difference R(S) subset S' subset J union R(S).

The inner update is still Y=X+D Q and Q'=R Q+Y. Let the inner diagonal
at coordinate j be a_j in theta+E; it is always nonzero.

At j in I, write X_j=t_j v_pivot,j+k_j, where t_j is the fresh outer
diagonal. If either Y_j or Q'_j is required zero, charge t_j using (4).
If both are zero, first impose (RQ)_j=0; the two remaining equations then
coincide. There is only one outer-coefficient charge at this coordinate.

At j in S outside I, the input is zero. If either Y_j or Q'_j is required
zero, use the ordinary inner coefficient a_j. These coordinates are
disjoint from those whose outer coefficients are charged.

Consequently an admitted support triple has the upper weight

    z^|J| psi^(-c_late),

    c_late = |(I union S) minus (J intersection S')|.        (5)

Sum (5) over the admitted J masks to obtain a positive later-step matrix.
This replaces the old exponent |S minus (J intersection S')|. It charges
input/state matches that were previously dropped when no new source
coordinate remained.

For example, the paired reset with I={0,1}, S={0,2}, and S'=empty needs
X_1=Q_2 as well as its inner-coefficient equations. The outer diagonal
inside X_1 now charges that match. This is a local necessary pattern, not
a constructed low-weight codeword.

### Sequential exposure and no double charging

First impose coefficient-independent zero equations, retaining earlier
nonzero requirements. Reveal the current inner coefficients at active
coordinates without charging them. Process the charged inactive inner
coordinates next, in the common quotient. Finally expose and charge the
active coordinates' fresh outer diagonals in a fixed order. Their affine
equation constants are known at exposure time. A previous equation that
kills a retained pivot form rejects the branch.

Reveal the remaining coefficients needed for the current step and append
its exposed source, output, and next-state nonzero requirements. Keep every
earlier requirement. No branch survival probability is renormalized away.

At an active coordinate already charged through t_j, do not also charge
the inner coefficient a_j for the same equation. Nor does (5) justify
subtracting every charged coordinate from the outgoing state-rank bound:
some charges survive through an exceptional coefficient without adding a
rank. Without a separate rank-sensitive refinement, use the safe bound

    1<=r'<=min(r+|I|,|S'|)

for a live next state. The stronger r+|I|-c_fresh bound remains restricted
to the actual fresh-source prefix.

## 5. Composition interface and unresolved calculation

Use the existing fresh-source active envelopes in the first b occupied
columns. Use the later-step matrix from (5) in the remaining 32 columns.
Its reset/live terminal functions can enter the same
[endpoint-tail construction](ENDPOINT_TAIL.md). Freshness of the relevant
outer diagonal persists throughout the preceding empty gap, which consumes
no input from that group and component.

For each raw occupied column, apply P_source twice before its four
chronological scalar-region matrices. Average those regional envelopes
under the actual paired route: keep pair collisions, short gaps, both
adjacent-pair parities, and the free terminal. Empty columns neither use
P_source nor advance the fresh-source prefix.

The global source factor remains tau^(8(h-32)-1). After the complete raw
experiment bound is obtained, multiply once by (q/(q-1))^64. The h-column
selection recurrence still counts each shared support set once; the
conditioning factor does not introduce another column-choice count.

The endpoint-aware implementation has now been evaluated with directed
80/120-digit arithmetic at tau=1024, z=.998, and lambda=.993. Its broad
shared-column union is about 5.751e162, so the class remains open. The
[online route recurrence](ONLINE_ROUTE.md) additionally retains rank and
support within a region; its first floating-point screen gives log2 union
286.351, still above one. Neither failed bound is a low-weight witness.

Both arguments use ordinary chronological equations, not lookahead
constraints involving unconsumed source forms. Partial columns need a
separate source-space argument, and no full-code guarantee follows from
this two-group class. The next test retains residual rank at region
boundaries before expanding the source class.

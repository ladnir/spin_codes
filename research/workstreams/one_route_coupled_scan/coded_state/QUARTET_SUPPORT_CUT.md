# Support cuts separate structural obstructions from counting slack

2026-10-27. The construction is unchanged. This note develops two distinct
tools: a sufficient condition for a structural low-distance witness, and
an exact graph inequality for a fixed routing/partition realization.
Neither tool alone supplies the required ensemble45%/40-bit certificate.

## An exact rowwise constraint count

Fix the route, every quartet partition, and every field map. Let S_t be
an allowed support for state Q_t, with S_0 empty. The state coordinates
inside S_t are variables; they are not required to be nonzero. Define

    A=sum_t |S_t|,
    D=sum_t |S_(t-1) outside S_t|.

Let Gamma_t(S_(t-1)) be the union of the quartets meeting S_(t-1). The
input X_t=Q_t-A_t Q_(t-1) is supported inside

    H_t=S_t union Gamma_t(S_(t-1)).

For outer row i, let h_i count its routed positions in the covered set H.
The RS[256,128] words supported on those positions have dimension
max(0,h_i-128). Requiring an arbitrary vector on H to belong to the outer
code therefore imposes

    F(S)=sum_i min(128,h_i)

independent constraints on its input coordinates. Pulling these equations
back to the A allowed state coordinates gives rank at most F(S). This
inequality does not assume independent random field coefficients.

Outputs Y_t=Q_t+delta Q_(t-1) are supported inside the union of consecutive
allowed supports. That union has A+D coordinates over the complete code.
For a target integer W, select z=max(0,A+D-W) of those output positions
and require their values to be zero. These add at most z equations. Thus
the allowed-state solution space has dimension at least

    A-F(S)-z.

If this number is positive, a nonzero solution exists and produces output
weight at most W. A nonzero state trajectory gives nonzero input and
output because both recurrences have identity diagonal blocks and Q_0=0.
It follows that

    A>F(S) and F(S)+D<W                                      (1)

are sufficient for a low-distance witness, for every choice of the field
maps respecting the fixed quartet partitions. Equivalently, whenever
A>F(S), this argument gives d_min<=F(S)+D+1.

The selected output positions here may include births or deaths, not just
support intersections. This is valid because S_t are allowed supports:
the solution may have smaller actual supports. Do not apply this argument
as though it counted exact nonzero support patterns.

Failure to find (1) does not prove distance. The argument is a sufficient
obstruction test, not a characterization of all low-weight codewords.

## A cut objective for a fixed inactive-row set

Fix a set T of r inactive outer rows, and put n=8192-r. Let C_T(S) count
its positions in H. Since each row outside T contributes at most128 to F,

    F(S)<=128n+C_T(S).

For a rational slope lambda=a/b in(0,1), consider

    C_T(S)+D-lambda A.                                       (2)

The minimum of (2) over every allowed support history is an s-t min-cut.
The cut's source side denotes selected state coordinates. Use these arcs:

- source to each state coordinate, capacity a;
- previous state coordinate to the same next coordinate, capacity b;
- each previous coordinate to its quartet auxiliary node, infinite capacity;
- the quartet node to each of its inactive input-position nodes, infinite capacity;
- each next-state coordinate to its inactive input-position node, infinite capacity;
- each inactive input-position node to sink, capacity b.

A quartet auxiliary node is needed only when it has inactive positions.
The initial state has no nodes and is fixed to zero. There is no death
edge after the final state. Infinite capacity can be aN+1: the all-empty
state cut costs aN, so a minimum cut cannot cross such an arc.

The inactive-position node is on the source side whenever its position
belongs to H. Its unit cost is paid once even if both causes select it.
After minimizing over auxiliary nodes, the cut capacity is exactly

    aN+b(C_T+D)-aA.

If the minimum capacity is V, every state-support history satisfies

    C_T+D >= lambda A - (aN-V)/b.                            (3)

The source helper verifies a feasible integer flow and independently
recounts A,C_T,D from the resulting state cut. Equality of their objectives
proves optimality; a favorable diagnostic number is not accepted merely
because an optimizer reported it.

## Full-size graph diagnostic

The C++ helper uses K=2^20, N=2^21, width512, r=2100, and256 regions.
Each region contains a permutation of8192 row identities. Each update has
an independently generated quartet partition. The reproducible test uses
mt19937_64 seed20261027 and GCC15.2.0/libstdc++ shuffling. This is one valid
pseudorandomly generated geometry, not a proof of its sampling distribution
or a statistical estimate of ensemble failure.

The inactive set is the first2100 row identities. The graph has3,159,042
nodes and7,263,696 directed capacity arcs before reverse residual arcs.
The verified results are:

| Slope | Exact offset in(3) | Minimizing A | C_T | D | F(S) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 22/100 | 71222/100 | 16801 | 141 | 2843 | 64730 |
| 24/100 | 106364/100 | 18086 | 350 | 2927 | 67824 |
| 26/100 | 782652/100 | 2095802 | 537082 | 0 | 1048576 |

None of these minimizing supports satisfies(1) at W=943718. This does not
rule out another support pattern or another inactive-row set. Inequality(3),
however, is verified for every support history for this fixed T and geometry.
It includes resets, arbitrary varying weights, and the free terminal state.

At slope0.24, the artificial pair from the earlier routing bound would have
A about1,039,773 and C_T+D about147,038. The right-hand side of(3) is about
248,482. Thus those averaged statistics cannot describe a distribution of
support histories on this particular graph. The numerical envelope averages
other geometries too; the graph check does not justify discarding them from
the ensemble proof.

## What the weighted bridge actually counts

Marked local coefficients retain uncovered positions, deaths, covered
inactive positions, and all inactive labels. The two representative
midpoint classes are(0,348,489) and(489,332,0). At p=0.70594552, their
expected regional inactive counts are2686.63 and2692.64, rather than2100.
The repeated artificial pair has active coverage about1,144,477, already
above the necessary129*6092=785868. The aggregate outer-coverage floor is
therefore not what excludes that diagnostic distribution.

Reducing p gives the following FLOAT bounds for these fixed classes:

| p | Growing-class log2 bound | Its inactive count | Decaying-class log2 bound | Its inactive count |
| ---: | ---: | ---: | ---: | ---: |
| 0.70594552 | 3018.464 | 2686.631 | 3289.799 | 2692.636 |
| 0.67 | 2903.133 | 2462.535 | 3173.781 | 2462.711 |
| 0.64 | 2852.946 | 2254.741 | 3123.672 | 2255.073 |

These are neither conditional code-state distributions nor exact-quota
coefficient extractions. The derivative of the fixed-class log bound with
respect to log odds is its expected inactive count minus r. This identifies
the direction of improvement without mistaking a positive bound for a code
obstruction.

The inherited field and selected-zero parameters are also not globally
optimal. Two bounded full513-state midpoint tests give:

| theta | ell | FLOAT log2 screen at r=2100 |
| ---: | ---: | ---: |
| 0.175 | 0.00377 | +404889.64 |
| 0.1 | 0.1 | +426592.92 |
| 0.13 | 0.03 | +303965.75 |

The last improves the full bound, but remains far from a certificate.
Other path families become dominant when the proof parameters change.

## Next ensemble obligation

The missing step is a probabilistic analogue of(3), uniformly over row sets
and state histories. For b>1 and a proposed offset B, positivity gives

    Pr[exists S: C_T+D<lambda A-B]
      <= b^(-B) E[sum_S b^(lambda A-C_T-D)].                  (4)

The all-empty history may be included in the right-hand sum. To cover
every row set of size r, multiply by binom(8192,r). Regional quotas and
quartet partitions must remain their actual distributions. This is a
geometry-only calculation; it has no selected-output-zero polynomial and
does not assume field refresh over the whole state.

`quartet_cut_moment.py` evaluates a FLOAT full-weight envelope for(4).
For a quartet with a incoming positions, next support V, and inactive set
T_local, its cost is c+d, where c=|T_local| if a>0 and c=|V intersect T_local|
otherwise. Its weight is b^(lambda|V|-c-d). Sum all supports and labels,
normalize incoming support by the uniform-partition binomial coefficient,
and use the preceding regional midpoint/quota bounds. A successful directed
geometry bound must then be combined with the field count and a complete
row-count union. No setup rejection rule is introduced here.

## Further cut diagnostics

### Canonical runs in a minimizing cut

A nonterminal run consists of consecutive selected states at one coordinate,
followed by an unselected state. Deleting a whole run of length L reduces A
by L and D by one. It cannot increase C_T, because coverage is monotone in
the allowed supports. The objective C_T+D-lambda A therefore changes by at
most lambda L-1. Every minimizing cut has lambda L>=1 for each such run.
Terminal runs are exempt: their deletion does not remove a death.

At slope0.24, every nonterminal run in a minimizing cut has length at least5.
The tested optimum has2927 nonterminal runs, shortest length5, and493
terminal runs. Exhaustive toy tests check36864 run erasures. This pruning
may reduce a cut enumeration, but does not yet bound its setup probability.

### A constructive flow without global optimization

The helper also constructs a feasible partial flow in chronological order.
Use integer capacities a=24,b=100. Each state node receives a new units
from the source. At an inactive position, serve those units immediately;
the position's sink arc retains capacity b-a. Within each quartet, use
the pooled residual sink capacity to serve the four previous-state queues.
Serve their largest queues first. Carry remaining work along the same-coordinate
arcs, each of capacity b, and discard any overflow. Discard all work still
queued after the final update.

Each service uses the quartet auxiliary arcs or the current direct arc.
The pooled service can be assigned to inactive sink arcs without exceeding
their residual capacities. Queue carry respects the same-coordinate arcs.
Removing the source-to-discarded-work paths leaves a feasible source-to-sink
flow. If L is total discarded work, including the terminal queues, its value
is aN-L. Weak flow duality therefore proves, for every support history on
this fixed geometry,

    C_T+D >= (a/b) A-L/b.

The observed L is352670, giving offset3526.70. This is weaker than the exact
offset1063.64, but obtained by a simple forward recurrence. The helper checks
work conservation and that its flow value does not exceed the verified
maximum. It does not claim a tail bound for L across random setups or row
sets. Neither this proof flow nor its queues are encoder operations.

### Tests of whether expansion alone resolves the gap

The first ensemble cut-moment bounds are far weaker than the fixed-geometry
cuts. At slope0.24 and r=2100, the midpoint envelope in(4) needs offsets
138256.38 at log2(b)=8 and107999.73 at log2(b)=4 to give a FLOAT60-bit bound.
These positive-envelope tests do not establish useful small-offset expansion.

Before investing in a queue tail bound, test what the desired expansion would
buy in the distance count. On the event R-lambda A>=-B, where R=C_T+D, and
for gamma>=1, the factor gamma^(R-lambda A+B) is at least one. Multiplying
the positive counting kernel by this factor bounds the bad-code event
intersected with that expansion event. Its complement still requires a
separate setup bound; no conditioning or renormalization is performed.

The helper implements this test with explicit `--cut-marker`, `--cut-slope`,
and `--cut-offset` options. Their defaults leave the original count unchanged.
Three exact rational toy checks verify the inequality direction. At r=2100,
theta0.13, ell0.03, and all513 state weights, the FLOAT results are:

| Extra assumption or tilt | log2 bound |
| --- | ---: |
| None | +303965.75 |
| R>=0.24A-5000, gamma=2 | +394432.61 |
| Same weighted test with B=0 | +389432.61 |
| No cut assumption, activity marker1.5 | +462968.00 |

The B=0 row follows by removing the exact additive5000 term; it is not a
claim that zero-offset expansion holds. These are specific proof-parameter
tests, not optimized lower bounds on what the method can achieve. They show
that the current choice of weighting does not exploit the new cut inequality
successfully. Do not promote the forward-flow tail as the final missing
lemma without first finding a useful joint distance bound. Next identify
the histories dominating that joint count and the rowwise constraints they
still evade. The construction remains fixed and the full certificate open.

## A proved geometry tail and stronger activity exclusion

This continuation replaces the difficult general cut tail by one scalar
geometry statistic. Fix an inactive-row set T. A state coordinate Q_t[j]
touches its own input position X_t[j] and, unless t is terminal, the four
input positions in its quartet at update t+1. Call it free if none of these
positions belongs to T. Let Z_T be the number of free state coordinates.

Every selected nonfree state coordinate touches a covered inactive position.
Each inactive input position touches at most five state coordinates: its
own next-state coordinate and four previous-state coordinates. Counting
these incidences proves, simultaneously for every allowed support history,

    A <= Z_T+5C_T.                                           (5)

Moreover, the set of all free coordinates has C_T=0. It is the unique
maximal allowed support with that property. This statement concerns allowed
supports, not existence of nonzero field amplitudes on all those positions.
On the same diagnostic geometry, Z_T=475819. That support has D=330239 and
F(S)=779756, so it is not a witness for the earlier rank obstruction.
The C++ helper checks(5) on196608 toy supports and verifies maximality
whenever C_T=0.

### A setup bound uniform over all row sets of size2100

Use an auxiliary experiment in which every input position has an independent
inactive label of probability p. Write a=1-p and choose v>1. Fix all quartet
partitions. The number of free coordinates at time t<4096 depends only on
labels at times t and t+1. Split these contributions by the parity of t;
assign the terminal contribution to the parity class whose pairs do not
use the final time. Within each class the pairs use disjoint time slices.

Within one pair, a quartet contributes zero unless all four labels at the
later time are active. Otherwise it contributes the number of active labels
at its four earlier positions. Its v-moment is

    L(v)=1-a^4+a^4(p+a v)^4.

The terminal contribution has moment g(v)^512, where g(v)=p+a v.
Cauchy--Schwarz between the two parity classes therefore gives

    E[v^(Z_T/2)] <= L(v)^(128*4095/2) g(v)^256.              (6)

This holds for every fixed partition realization. No independence between
the two parity classes is asserted. The actual route is the auxiliary
label experiment conditioned on exactly r inactive labels in each region.
Put pi=binom(8192,r)p^r a^(8192-r). Positivity, Markov's inequality, and
the union over row sets give

    Pr[exists T of size r: Z_T>=z]
      <= binom(8192,r) pi^(-256) v^(-z/2)
         L(v)^(128*4095/2) g(v)^256.                        (7)

This conditioning concerns the route only. The accepted field-map setup
does not change its distribution. All code randomness is still sampled
once; the artificial p is a proof parameter, not new setup randomness.

At r=2100, p=29/100, v=39/32, and z=599529, the replay squares(7) and uses
directed rational logarithms. The upper logarithm of that squared bound is

    -33012067465555/274877906944 < -120.

Thus, except with probability at most2^-60, every inactive-row set of size
2100 simultaneously satisfies Z_T<=599528. Inequality(5) then holds with
that same upper bound for every support history.

### Improved low-activity certificate at this row count

The old zero-selection-free field count for fixed supports is
2q^(A-128n-C_T), where q=2^32 and n=6092. On the preceding good-geometry
event, put z0=599528 and use

    C_T >= max(0,ceil((A-z0)/5)).

The integer exponent A-max(0,ceil((A-z0)/5)) is nondecreasing in A.
For a threshold b<=N/2, summing all activities1..b and all labeled supports
therefore gives the conditional upper bound

    2 binom(8192,2100) b binom(N,b)
      q^[b-128*6092-max(0,ceil((b-z0)/5))].                 (8)

Here binom(N,b)<=N^N/[b^b(N-b)^(N-b)], obtained by bounding the corresponding
binomial probability by one. This entropy bound avoids large coefficient
tables. At b=747590, the directed rational upper logarithm of(8) is

    -97240570033871/1099511627776 < -60.

Combining the field and geometry allowances proves:

    Pr[exists nonzero J-message with exactly6092 active outer rows
       and total state activity A<=747590] <= 2^-59.        (9)

The union covers every row set and arbitrary histories, including resets
and the free terminal state. This improves the previous threshold694488
at this row count. It does not cover other row counts or prove minimum
distance. As with the earlier activity lemma, its statement is over J;
a restricted activity property is not silently transferred to extension
messages. No setup rejection rule or construction change is introduced.

`quartet_zero_neighbors.py` replays(7)--(9) with the frozen rational p,v.
Its optional bounded FLOAT screen only proposes parameters; every reported
certificate uses directed rational comparisons. Exact enumeration checks
4096 toy label patterns, the parity moment inequality and quota conditioning,
and279 monotonicity steps for the field exponent. The replay uses only
scalar arithmetic and creates no files.

Next extend the activity exclusion across row-count intervals and use the
same geometry statistic jointly with the output count. The complete
45%/40-bit certificate is still open.

## Continuation: every unresolved row count

The later [joint support/output count](QUARTET_JOINT_PRIOR.md) proves that
a weighted cut with death cost0.74 would suffice for the r=2100 distance
class. It generalizes the fixed-geometry diagnostic below to separate
death and coverage costs. The required ensemble tail remains open.

The interval replay covers105<=r<=8191. The existing high-turnover certificate
already covers every message with r<=104; r=8192 cannot describe a nonzero
message. This new replay is an activity exclusion, not a distance proof.

Partition the range into127 intervals, starting at105, of width64, with a
short final interval. For endpoints l,u, use the frozen formulas

    rho_mid=(l+u)/(2*8192),
    p=min(rho_mid+3/100, (1+rho_mid)/2),    v=39/32.

The full replay performs no FLOAT search. The squared logarithm in(7),
before subtracting z log2(v), is

    (2-512)log2 binom(8192,r) + an affine function of r.

Binomial log-concavity makes this expression convex in integer r. Its
endpoint upper bounds cover every intermediate row count. Choose one integer
z per interval that makes both endpoints at most-148. Then each row count
has geometry failure at most2^-74, including its complete row-set union.

For the field bound, choose an integer b and threshold A<=b(8192-r).
Use fixed Bernoulli tangents for the two binomial coefficients in(8):

    binom(m,k) <= x^(-k)(1-x)^(-(m-k)),    0<x<1.

For row sets take x=(l+u)/(2*8192); for supports take
x=b[8192-(l+u)/2]/N. Bound the number of activity terms by b(8192-l).
The field exponent is at most the minimum of two affine functions,

    A-128n,    (4/5)A-128n+(z-1)/5.

The resulting upper logarithm is piecewise affine in r. Check both endpoints
and the crossover A=z-1, if it is inside the interval. Directed rational
logs verify a field allowance at most2^-74 for every accepted b. Searching
for b only proposes a threshold; its final check establishes validity.

There are8087 row counts and two allowances per count, giving

    2*8087*2^-74 < 2^-60.

The thresholds improve the original ones at7712 row counts. Retain the
maximum of old and new thresholds, charging each global allowance once.
For2089..2152, the replay gives z-1=604080 and b=122. At r=2100 the
remaining activity is at least743225. The separately retained single-row
checkpoint gives747591 with its own2^-59 allowance.

The executable rule is `certified_row_bounds(r)`; the complete replay is
`quartet_zero_neighbors.py --intervals 105 8191 --width 64`. It checks496
small exact log-concavity inequalities and638 Bernoulli tangent inequalities.
All127 interval checks complete in about7 seconds remotely, without files.

### Coupled output-count tests

On the good geometry event, gamma^[C_T-A/5+(z-1)/5]>=1 for gamma>=1.
The full helper implements this together with the stronger activity floor.
At r=2100, theta0.13, ell0.03, tau0.05, and the17-point p grid:

| Added weights | FLOAT log2 screen |
| --- | ---: |
| None | +271673.74 |
| Geometry gamma1.25 | +283239.96 |
| Geometry gamma1.25 and activity marker1.05 | +286010.29 |

The first trial shifts artificial state flow near weight93; the activity
marker moves it back to an alternating near-empty/near-full flow. These
weights do not improve the bound. The structural exclusion is proved, but
its current scalar insertion is not a closure argument.

The next [coverage-aware and trajectory-aware exclusions](QUARTET_LOW_ACTIVITY.md#coverage-aware-and-trajectory-aware-exclusions)
retain more of the field count. Stress tests assuming the small-offset cut
also remain positive. Proving its ensemble tail alone has not been shown
sufficient for the distance theorem.

## Validation summary

`quartet_support_cut.cpp` exhaustively checks196608 small support cuts and
three complete flow optimizations. Full runs verify integer capacities,
flow conservation, and matching independently recounted cut objectives.
The graph storage is explicitly budgeted below512MiB. Execution uses one
thread, a120s timeout, and a768MiB process address-space cap on Peach.
The three full runs took about7s,8s, and29s. A small temporary build directory
was validated under /tmp and removed after execution. No build or dump was
written locally or retained remotely.

`quartet_bridge_moments.py` checks an exact rational bridge,81 coverage
partition identities, and81 finite-difference label derivatives. The cut
moment helper checks25 local coefficients against an independent rational
binomial formula. All research jobs remain serial and file-free apart from
the transient C++ build. Existing certificates and production code are
unchanged.

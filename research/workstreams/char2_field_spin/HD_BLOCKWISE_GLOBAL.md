# Blockwise SPIN: structural identities and a capped routing bound

The selected blockwise map has an exact global identity: its inner encoder
is an involution. This does not imply a distance failure. Its quiet-state
subcode and the repeated-row pattern below also do not reject 40% distance.

There is nevertheless an actual obstruction to an uncapped codeword first
moment. A fixed set of only 39 outer rows gives expected low-weight
projective count greater than 2^780.17, through rare routes with large
message kernels. Counting those routes once instead succeeds: the entire
dimension-forced, full-row window family has probability below 2^-1511.82
for state16. This bound includes all source-row sets, packet-aligned window
starts and lengths, connections across region boundaries, and a free
terminal state.

The capped event is not all distance failure. Partial-row shortening and
general low-weight output supports remain unbounded. The results below give
one rigorous decomposition for the restarted proof, not a certificate.

## Exact construction under review

Use the [blockwise restart](HD_BLOCKWISE_RESTART.md): R=16384 independent
GRS[128,64] source rows, one independent uniform row permutation in each
coordinate region, and no mixing pre-stage or additional row shuffle.
There are K=2^20 message symbols and N=2^21 output symbols over
F=GF(65536). Shared scalar coefficients act on the eight field lanes of a
128-bit record. A one-lane witness therefore also witnesses record weight.

The packet width is t=64 and the state dimension is s=16. The terminal
state is free. Starting at Q_0=0, the inner emits and updates by

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The fixed maps in [HD_BLOCKWISE_MAPS.md](HD_BLOCKWISE_MAPS.md) use the
embedded GF256 sets U={0,...,63} and V={128,...,143}. Write
P(Z)=product_(v in V)(Z+v) and let D=P' be its constant nonzero derivative.
Then A_(u,j)=P(u)/(D*(u+v_j)), and C=A^T. These maps satisfy

    rank A = rank C = 16,       C A = 0,
    im A is an MDS [64,16,49] code,
    ker C is an MDS [64,48,17] code.

Each fresh M_i is the invertible rank-one update I+u_i v_i^T described in
[HD_BLOCKWISE_INNER.md](HD_BLOCKWISE_INNER.md). Its setup is independent
between packets and independent of the route. None of the routing bounds
below needs its exact distribution; they hold for every realization of
the causal linear inner. The distribution matters when proving the
remaining distance event.

## An involution, not an accumulator shortcut

Multiplying the emission equation by C gives CY_i=CX_i. Therefore the
inverse can use the same recurrence as the forward map:

    Q_(i+1) = M_i Q_i + C Y_i,
    X_i = Y_i + A Q_i.

Starting from the same zero state, induction shows that the second
application reproduces every first-pass state. It cancels the two copies
of A Q_i. Thus E^2=I for the full inner map E, with the same fixed M_i
schedule. The statement uses characteristic two and CA=0; it does not
replace either recurrence with an independent uniform-state process.

The identity supplies a useful proof interface. For a fixed proposed
output word Y, its inverse states depend only on the fixed syndromes CY_i
and the fresh packet mixers. They do not depend on the route. The remaining
test is whether the inverse word, after undoing the route, belongs to the
outer row code. Output weights can still bias the state measure, so the
local point/line/uniform law must retain that measure.

There are P=N/t=32768 packets. Every input whose first P-1 packets have
zero C syndrome is fixed pointwise by E. Its dimension is

    N - s(P-1) = 1572880.

Conversely, A is injective, so E(X)=X forces all emitted states to vanish
and hence forces those first P-1 syndromes to vanish. This is exactly the
fixed subspace; the last syndrome is unrestricted because there is no
subsequent emission. Its intersection with any routed K-dimensional outer
code has dimension at least

    K - s(P-1) = 524304.

The large fixed subcode is an important global constraint. Dimension alone
does not bound its minimum weight below 40%. No such inference is made here.

The most obvious invariant source pattern is also benign for this target.
Repeat one outer codeword identically in every source row. Each region is
then constant, so every regional permutation leaves the input unchanged.
The constant packet vector belongs to im A and hence to ker C. The state
stays zero and Y=X. The minimum nonzero weight among these repeated words
is exactly 65R, because the outer has distance65. Its relative weight is
65/128=0.5078125, above both 40% and 45%.

## A direct first-moment obstruction

Fix a set J of a outer rows. Use the balanced 40% packet windows from
[HD_BLOCKWISE_GEOMETRY.md](HD_BLOCKWISE_GEOMETRY.md): 77 regions have length
6528 and 51 have length6592. Their total length is 838848<0.4N. Use prefix
windows in the first127 regions and a suffix window in the last region.

Let E_J be the event that every coordinate of every row in J routes into
these windows. For independent uniform regional permutations,

    Pr_route[E_J] = product_c C(W_c,a)/C(R,a).

On E_J, all 64a message dimensions in the selected rows are supported
inside the windows. Setting the state to zero after the first127 windows
costs at most127s linear equations. Hence, for every inner setup, a message
space of dimension at least d=64a-127s has nonzero output only within the
windows. Injectivity makes its nonzero messages nonzero codewords.

Let Z be the number of output-code scalar lines of weight at most838860.
The actual code, with all outer parity checks enforced, satisfies

    E[Z] >= Pr_route[E_J] * (65536^d-1)/65535.

At s=16, the exact event and conservative directed lower bounds are:

| Fixed a | Dimension lower d | log2 Pr_route[E_J] | log2 E[Z] lower |
| ---: | ---: | ---: | ---: |
| 39 | 464 | -6611.8201838145 | 780.1798161854 |
| 40 | 528 | -6781.6940633279 | 1634.3059366721 |
| 64 | 2064 | -10863.7909982012 | 22128.2090017988 |

The helper uses exact binomial fractions and directed160-bit logarithms.
The final column conservatively replaces the projective count by
65536^(d-2). There is no union factor over J, so repeated counting of
codewords with smaller active row sets does not enter the calculation.

This is specific to the direct blockwise construction, not an imported
accumulator argument. Increasing the number of refresh rounds cannot fix
this first moment: the dimension bound holds for every causal inner setup.
It also does not reject the 40-bit distance goal. The displayed route
events have probabilities thousands of bits below2^-40.

## Counting the full-row traps once

Define the event B_full as follows. There is a nonempty row set J of
size a and one nonempty packet-aligned interval in every region, with
total length at most838860, such that:

- all b*a coordinates of J route into the chosen intervals;
- if the intervals form c connected components in the global stream,
  then 64a>s(c-nu), where nu=1 when the last component reaches the end.

Thus B_full contains every dimension-forced window trap that uses all
message dimensions of a chosen set of complete outer rows. It includes
the preceding first-moment witnesses. It does not include every trap
formed by partially shortened rows or a rank-deficient exit system.

We bound B_full over the route alone, uniformly in the mixer setup.
For fixed J and fixed windows, sampling without replacement gives

    Pr[all J coordinates contained]
      = product_j C(W_j,a)/C(R,a)
      <= product_j (W_j/R)^a.

Here b=128 and each region has B=R/t=256 packets. The total packet budget
is F=13107. Use the connected-window counting in
[WINDOW_CAPACITY_MOMENT.md](WINDOW_CAPACITY_MOMENT.md), specialized to all
b regions being used. The required formulas are given next so the finite
bound does not depend on treating disconnected windows as independent.

Fix c components, of which v occupy one region. There are

    H(c,v)=C(c,v)*C(b-c-1,c-v-1)

ordered component-size patterns when v<c. The all-singleton case has
H=1 only when c=v=b. The components force f=b-2c+v full regions and
leave d=2c-v=b-f variable window lengths. Put L=F-fB. Their lengths have
at most C(L,d) positive integer assignments. Singleton components have
at most B starts; the other endpoints are anchored at regional boundaries.
Thus the number of such fixed families is at most H(c,v) C(L,d) B^v.

The terminal-free cases have the following adjustments:

- For a final singleton, replace C(c,v) by C(c-1,v-1), replace B^v by
  B^(v-1), and charge c-1 exits.
- For a final component spanning several regions, replace C(c,v) by
  C(c-1,v), increase f by one, decrease d by one, and charge c-1 exits.

The second adjustment accounts for the final endpoint region being full.
Infeasible length counts are omitted. Counts that allow a variable length
above B are retained as upper bounds. No placements among unused regions
are needed because every region contains a coordinate of every full row.

In a group with f forced-full regions and d variable regions, put
p=L/(dB). Arithmetic-geometric mean bounds the per-row containment factor
by r=p^d. The smallest a that satisfies the dimensional condition is

    a0 = floor(s*(c-nu)/64)+1.

Sum over every source set and every larger cardinality. Since consecutive
terms have ratio at most eta=((R-a0)/(a0+1))*r<1,

    sum_(a=a0)^R C(R,a) r^a <= C(R,a0) r^a0/(1-eta).

Multiply by the appropriate window-family count, then sum all c,v and
the three terminal cases. Directed160-bit arithmetic evaluates7779 groups
and proves

    Pr_route[B_full] < 2^-1511.82711066       at s=16,
    Pr_route[B_full] < 2^-4051.75198256       at s=32.

Unlike the raw first moment, this bound does not multiply a rare route by
the number of codewords in its bad kernel. Boundary connections are material:
the largest s16 summand uses64 components and a free final singleton,
not128 separately charged exits.

## What the capped bound does and does not compose with

Let Bad be failure of the desired distance. The valid decomposition is

    Pr[Bad] <= Pr[B_full] + Pr[Bad and not B_full].

The first term now has ample40-bit slack. Since B_full depends only on the
route, conditioning on its complement leaves the independent mixer setup
law unchanged. The regional permutations themselves are no longer independent
uniform permutations under that conditioning. A subsequent route calculation
must retain the conditioning or use an explicitly dominating unconditioned
quantity; it cannot reuse the original route law as an equality.

Partially shortened rows remain a concrete gap. For one fixed window family,
the full shortened dimension is S=sum_i(t_i-64)_+. Negative association of
regional fixed-size samples and convex smoothing give the known bound

    E[z^S] <= m(z)^R,
    m(z)=E[z^((Binomial(128,0.4)-64)_+)],       z>=1.

For the conservative crossing S>=127s+1, the rational mark z=721/500 at
s16 gives a directed fixed-profile upper below2^-612.20530922. Multiplying
by the crude count [256*257/2]^128 of all one-interval-per-region families
makes the log2 upper positive1308.5146. This diagnostic does not close that
family, and it does not establish actual failure. Connected families can
need fewer exits, so they are not repaired by this fixed127-exit diagnostic.

The useful next decomposition is therefore source-aware. Retain the exact
shortened dimension of partial rows, the connected exit count, and the
actual rank or existence event. Use the local emission/syndrome transfer
from HD_BLOCKWISE_MAPS.md for the remaining output histories. A bound that
returns to the uncapped actual-code first moment over all routes cannot
succeed, even with perfect local state calculations.

## Validation

The new `theory/hd_blockwise_global.py` supplies the directed receipts and
`theory/hd_blockwise_global_tests.py` has nine passing tests. Small GF4 cases
check CA=0, the involution, packet syndromes, the pointwise fixed subspace,
and repeated outer rows. Exact combinatorial tests check containment,
source-set geometric tails, and connected-family counts against exhaustive
small window layouts. The remaining tests replay the actual first-moment
obstruction, the s16 capped gate, and the failed naive partial-row union.

```text
python -B research/workstreams/char2_field_spin/theory/hd_blockwise_global.py
python -B research/workstreams/char2_field_spin/theory/hd_blockwise_global_tests.py
```

Only this note and the two new helpers were written. No native build,
remote call, benchmark, existing-file edit, or commit was performed.

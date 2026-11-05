# A rank interface for full-alphabet messages

This rank argument applies to the
[cheap refresh construction](CHEAP_LANE_MIXING_20261016.md), including
messages with full GF256 lane span. The general continuation bound remains
open. The one-group section closes the class with one actual active source
group at 40% distance, including all lane ranks. Independent reviews checked
the local lemma, its reveal order, and the one-group continuation bound.

The approach follows the common-quotient rank accounting in the adjacent
`gen-BAA/field-universality/paper/distance-proof.tex`. It does not identify
the limited coefficient set with the full field.

## Track equations instead of counting field-valued messages

Fix an outer setup, route, and source group set S. Let V be the F-linear
space of messages supported on S, of dimension D=128|S|. Fix a sequence
of state and output supports. Every imposed zero coordinate is a linear
form on V; let W be the span of the forms imposed so far. Let R be a
family of forms that the support sequence requires to be nonzero.

Retain a branch only when every form in R lies outside W. This is a
necessary condition for a witness, not an assertion that a common
witness exists. For 1<=tau<=65280, define

    Phi(W,R) = 1[branch retained] tau^(-dim_F W).

A nonzero witness has dim W<=D-1. Therefore, pointwise in the setup,

    1[the fixed support sequence has a nonzero witness]
        <= tau^(D-1) Phi_final.

Summing support sequences and then averaging setup gives an existence
bound. No factor |F|^D is introduced. The parameter tau weights rank;
it is not a replacement for |F| in a message-counting union bound.

## One output packet

Fix the pre-refresh history. Let x be the 64 input forms on V and w
the 16 pre-diagonal state forms. Current state and output forms are

    q_j = lambda_j w_j,
    y   = x + A q,

where the independent lambda_j are uniform on the actual set
S_lambda={c+b theta:c in E, b in E minus zero}, of size 65280.
Fix a proposed state support H and output support J. Put Z=J^c.

First impose w_j=0 for j outside H, and append w_j for j in H to R.
These conditions can be imposed before revealing lambda because every
lambda is nonzero. Discard the branch if a required form lies in W.

Set M=A[Z,H], and use its **actual** matrix rank r. Let L have rows
spanning the left kernel of M. Impose L x_Z=0, again discarding a dead
branch. Let d be the total rank increase from these preliminary equations,
including the inactive-state equations. All of these operations precede
revealing lambda.

Choose r pivot columns P of M and reveal the coefficients indexed by
H minus P. Numeric row elimination leaves exactly r equations

    lambda_p w_p + u_p = 0,  p in P.

Each u_p depends only on known forms and nonpivot coefficients. It is
independent of every unrevealed pivot coefficient. Reveal the pivots in
any fixed order, appending each resulting equation to W and retaining
every previous nonzero requirement in R.

On a live branch, w_p is outside W. If u_p is outside W+F w_p, the new
equation adds one rank for every lambda_p. Otherwise write
u_p=kappa w_p modulo W. All lambda_p except -kappa force the required
w_p into W and kill the branch. The one possible surviving coefficient
has probability at most 1/65280. Thus each pivot contracts expected
potential by at least tau^(-1), and

    E[Phi_after | pre-diagonal history] <= tau^(-d-r) Phi_before.

After the pivots, active output forms may be appended to R; doing so
cannot increase the potential. Omitting those nonzero requirements is
also safe for an upper bound. They must not be used to condition on
unrevealed pivot coefficients earlier in the argument.

Earlier pivots can change W and the dependence classification of later
forms. They do not reveal later pivot coefficients, so the same
conditional argument still applies.

## The actual emission matrix matters

Each four-coordinate state block is encoded with the normalized basis

    1, x, x(x+1)/6, x^2(x+1)/6,

evaluated at the 16 GF256 points 0 through 15. This is an RS[16,4]
generator, but it is not superregular in this coefficient basis. Row
x=0 is (1,0,0,0). Selecting that row and state column 1 has rank zero,
not min(1,1). The independent audit rejected that tempting shortcut.

Since A has four diagonal blocks over E,

    rank_F A[Z,H] = sum_b rank_E A_b[Z_b,H_b].

For z=|Z_b| and h=|H_b|, valid cardinality bounds are

    max(0, h-max(0,4-z)) <= rank A_b[Z_b,H_b] <= min(z,h).

The lower bound follows by deleting 4-h columns from a matrix of rank
min(z,4). In particular, z>=4 gives rank h. A useful finite operator
should retain actual masks or use the conservative lower bound; it must
not use the rejected min formula as an equality.

The scalar `basis()` in `hd_grouped_forward.cpp` and `innerCoset()` in
`HdGroupedInner.h` specify the actual matrix. H16, the separate refresh
interpolation matrix, is superregular; do not confuse these two matrices.

## Chronology and composition

At output step t>=1, the pre-diagonal forms are

    w_t = H16 D1[t-1] (Q[t-1] + C X[t-1]).

Output t-1 does not depend on lambda[t-1], so that diagonal may remain
unrevealed until output t. The outer and route fix X[t]. This supplies
the fresh randomness required by the lemma. Step zero has zero state and
only deterministic output equations. The refresh after the last emitted
output cannot be charged to the distance bound.

The increments d are history-dependent. They cannot be pulled out of
expectations without a proved lower bound. For a fixed support sequence,
the weighted potential

    Psi_t = tau^(sum_{s<=t}(d_s+r_s)) Phi_t

is a nonnegative supermartingale. If every successful realization has
accumulated cost at least L, its probability is at most
tau^(D-1-L) Phi_0. Alternatively, a positive continuation operator can
retain the costs inside its conditional sum.

There are only 16(32768-1)=524272 refresh scalars affecting outputs at
K20, versus K=1048576 source dimensions over F. Charging every refresh
pivot but discarding deterministic ranks still leaves at least
1048575-524272=524303 unpaid dimensions. The deterministic left-kernel
equations and actual outer/route structure are therefore essential.

The next task is to bound accumulated costs jointly with output support
and routing constraints. This local lemma alone does not supply that
bound, nor does it justify multiplying marginal packet probabilities.

## A global reduction to scalar constraint counts

The next audit supplies a global upper-bound framework, including inputs
with full E-lane span. It does not yet give a useful numerical certificate.
The resulting first operator loses too much on sparse auxiliary inputs.

Use the full message space V=F^K, not a space restricted to a proposed
actual source support. Fix the outer, route, and all D1 matrices. Fix
the state-support path H_t and sets Z_t of requested output zeros.
Here H_0 is empty. Write R_t=H16 D1[t]. For a nonfinal packet define

    M_t = [ A[Z_t,H_t]          ],
          [ R_t[H_(t+1)^c,H_t] ]

    N_t = [ I[Z_t,:]             ].
          [ R_t[H_(t+1)^c,:] C  ]

The corresponding equations are N_t X_t + M_t Q_t[H_t]=0. The lower
rows impose the next state's inactive coordinates: the following diagonal
is nonzero and therefore does not change that support. For the final
packet omit these lower rows. At packet zero, H_0 is empty and M_0 has
no columns.

Let L_t span the left kernel of M_t. All equations L_t N_t X_t=0 are
over E and independent of every lambda. Collect them into a matrix B on
the original K source coordinates, and impose them before any lambda
is revealed. This includes future equations, but does not condition on
future random coefficients. Put r_t=rank_E M_t.

For a surviving quotient of dimension d over F, replace the earlier
power potential by

    P_q(d) = (q^d-1)/(q-1),  q=256.

Set the potential to zero when a required nonzero form enters W. A
nonzero witness leaves d>=1 and hence potential at least one. For an
independent new pivot equation,

    P_q(d-1) = (P_q(d)-1)/q <= P_q(d)/q.

For a dependent pivot, at most one of the 65280 actual coefficients
retains the branch. Its contribution is at most P_q(d)/65280, also at
most P_q(d)/q. Retain all earlier active-state requirements, and append
the next active-state requirements after the current pivots. Consequently,

    Pr_lambda[the fixed path has a nonzero witness | outer, route, D1]
       <= q^(-sum_t r_t) (|ker_E B|-1)/(q-1).

The initial potential has this exact scalar interpretation because B is
over E. These scalar kernel elements are auxiliary witnesses; they are
not projections promised to preserve the actual F-message support.
The subtraction of one excludes the zero scalar message before any
outer-density relaxation. Two independent reviews checked this reduction.

### Use positive zero-subset weights

A direct union over output supports wastes a completion factor. Instead,
for 0<z<1 and any output Y, use

    z^wt(Y) = sum_{Z subset zeros(Y)} z^(N-|Z|) (1-z)^|Z|.

If some nonzero word has weight at most L, the weighted sum of existence
indicators for zero sets is at least z^L. Apply the preceding path bound
inside this positive sum. This step concerns existence indicators, not a
sum over all F-valued messages.

Expand |ker_E B|-1 as a sum over nonzero scalar messages. Partition that
sum by the scalar messages' active source groups, denoted S. Apply the
preserved scalar outer point-density bound to each group in S. It
dominates their routed outputs by independent uniform E values on a
set I_t(S) in packet t. No uniform-F outer bound is asserted.

For auxiliary x sampled uniformly from E^I, elementary elimination gives

    q^(-r_t) Pr[L_t N_(t,I) x=0]
       = q^(-rank_E [M_t, N_(t,I)]).

Thus the local matrix entry for the first bound is

    K_(H,H')(I,z,R)
      = sum_{Z subset [64]} z^(64-|Z|) (1-z)^|Z|
                           q^(-rank_E [M, N_I]).

Equivalently, it is the expectation of

    z^wt(x+A_H u) * 1[(R(Cx+u))_(H'^c)=0]

for independent uniform x in E^I and u in E^H, with zeros elsewhere.
The auxiliary u is allowed to contain zeros. This is an upper-bound
counting device, not the actual distribution conditioned on state mask H.
The final packet uses the same expression without the next-state event.

The fractional-moment order remains available. First average lambda with
the alpha=1 path bound, then cap the bad-event bound by one. Partition
scalar witnesses by S, use subadditivity, and apply Jensen to the outer
setup conditional on route and D1. Apply the scalar density bound there.
Only then split state paths and average their independent D1 matrices.
The positive output-zero sum stays inside each K before taking K^alpha.
The resulting factor is z^(-alpha L) times a product of entrywise
fractional powers, for 0<alpha<=1. Parameters may be selected per scalar
source-count class, but not retrospectively per unexposed random outcome.

### Why this first operator is insufficient

The scalar support S need not equal the support of the actual F witness.
In particular, I_t(S) empty does not permit forbidding a zero-to-nonzero
actual state transition. That restriction would discard possible F
witnesses. Local feasibility tests using I_t(S) as the actual input
support are invalid in this ambient reduction.

There is a quantitative diagnostic that does not require constructing the
large mask matrix. Since [M,N_I] has at most |H|+|I| columns,

    K_(H,H')(I,z,R) >= q^(-|H|-|I|).

The positive zero-subset weights sum to one. Summing the lower envelope
over intermediate state masks in R packets gives

    q^(-alpha sum_t |I_t|) (1+q^(-alpha))^(16(R-1)).

This is a lower bound on the proposed upper-bound expression, not on the
actual bad-event probability. It exposes artificial state histories in
the relaxation. At K20, R=32768 and a one-group scalar support has
sum_t |I_t|=256. Even with the sharper density beta=256^160/255^32,
the one-support expression is at least

    [beta / (255 * 256^256)]^alpha
       * (1+256^(-alpha))^(524272),

where dropping z^(-alpha L)>=1 only decreases this floor. Both terms
decrease with alpha on (0,1], so alpha=1 gives its minimum. It is still
approximately 2^1916.98. A coarse exact check already proves it exceeds
2^1015: (257/256)^256>2, floor(524272/256)=2047, and
beta/(255*256^256)>256^(-129)=2^(-1032). Increasing numerical precision or optimizing markers
cannot fix this operator. This is not evidence of a low-distance word.

The next refinement should retain active-state linear forms or their
rank in the common quotient. Another valid interface uses a two-byte
scalar shadow: P_65536(d) contracts by at most 1/65280 per pivot and
counts the E-defined kernel over GF(256^2). A direct density argument would
require a new rank-two outer bound; the old scalar density cannot be
squared. Alternatively, Cauchy--Schwarz can reduce a two-shadow moment to
two scalar moments, with fractional exponent at most 1/2. That route must
retain the correlated union of the two supports and is not evaluated here.

The exact [emission-rank checker](theory/cheap_lane_emission_ranks.py)
now covers all 1048576 row/state-mask restrictions of one RS block.
It explicitly checks every four-row minor and enumerates the exceptions
with at most three zero rows. It supplies exact local data, not the
missing state-aware continuation bound.

### A valid support-preserving sparse refinement

For sparse actual messages, one scalar shadow can suffice with a controlled
loss. Fix the actual active-group set T, with r=|T|<256, and restrict
the message space to V_T=F^(128r). Let U=ker_E B for the preimposed
E-valued equations, and let d=dim_E U. Discard a branch if any group
projection is zero on U. Scalar extension then makes an actual witness
with group set T impossible.

Otherwise each group projection has a proper kernel in U. At most
q^(d-1) scalar vectors vanish on a specified group. A union bound shows
that the number n_T of vectors in U nonzero on every group satisfies

    n_T >= (1-r/q) q^d,
    P_q(d) <= n_T / ((q-1)(1-r/q)).

For r=1 the sharper identity is n_T=q^d-1, so there is no extra loss.
Use this comparison only for the initial E-defined quotient, before
averaging the outer or applying its point-density bound. It is not a
replacement potential after later equations become F-valued.

Now the actual witness and the counted scalar shadows have the same
group set T. The route positions I_t(T) are legitimate allowed-input
sets for local feasibility checks. In particular an empty actual packet
cannot activate a zero state. The multiplier 1/(1-r/256) costs at most
eight bits for r<=255. This supplies a concrete next sparse gate;
the refined local operator and its complete count remain to be evaluated.

For larger r, a bounded number of missing scalar groups gives another
valid initial comparison. Sample a uniform nonzero u from U. For each
group, the probability that its projection vanishes is strictly below
1/q. Therefore the expected number of missing groups is at most r/q.
For any integer h>=0 with c_h=1-r/(q(h+1))>0, Markov's inequality gives

    P_q(d) <= sum_(S subset T, |T minus S|<=h) n_S / ((q-1)c_h),

where n_S counts nonzero vectors of U whose exact group support is S.
The choice h=ceil(2r/q) has c_h>1/2 and allows at most 64 missing groups
at K20. Smaller h may give a better tradeoff between this prefactor and
the number of missing-group choices.

This comparison preserves the actual allowed positions I_t(T), while
the scalar reference kernel uses I_t(S). An exact route calculation must
retain their correlation. Alternatively, deleting input columns gives
the conservative comparison used in the dense-band certificate below.
Neither approach silently sets S=T or treats the sets as independent.

## A closed bound for one actual source group

The scalar-shadow relaxation pays for states that cannot occur. For one
actual source group, we can instead retain its allowed input positions
and count deterministic source constraints directly. This needs neither
an outer point-density bound nor a scalar shadow.

Fix one of the 8192 source groups. Its message space is V=F^128. Fix all
invertible outer mixers and all D1 diagonals. Each of the 256 regions has
one allowed input coordinate. Its packet is uniform among 128 packets,
independently across regions. These facts hold for either adjacent-pair
matching: the group's partner is inactive, and the fair swap gives a
uniform position inside the chosen pair.

### Deterministic constraints from quiet columns

Fix a state-support path, including the outgoing support of the last
packet. Call an outer column quiet on this path if the state is zero
throughout all four of its regions. At an allowed input coordinate i,
zero incoming and outgoing states imply C_i X_i=0. The column C_i has
a nonzero constant-basis entry, so X_i=0.

Thus a quiet outer column forces all four mixed RS evaluations to zero.
The column mixer is invertible. For c distinct quiet columns, these
equations have rank 4 min(c,32) over F, since each of the four constituent
RS[64,32] codes is MDS. They depend only on the fixed outer and route,
not on any fresh lambda. Impose them before revealing lambda. A nonzero
message is impossible when c>=32.

For any real 1<tau<=65280, use P_tau(d)=(tau^d-1)/(tau-1). Its independent
rank-step contraction is at most 1/tau, and its dependent-pivot contraction
is at most 1/65280<=1/tau. No finite field of size tau is needed. Before
the remaining pivots, the potential is bounded by

    tau^128/(tau-1) * tau^(-4c).

This bound is also conservative for c>=32, when the actual potential is
zero. We may either retain those impossible paths or discard them.

### A seventeen-state upper-bound operator

Let h count actual nonzero incoming state coordinates, and let j count
nonzero outgoing coordinates. Put w(h)=13 ceil(h/4). At least ceil(h/4)
of the four RS[16,4] emission blocks are nonzero. An empty packet emits
at least w(h) nonzero records. A packet with one allowed input coordinate
emits at least w(h)-1 when h>0.

For an empty packet, a live state cannot become zero, and h+j>=17 by
superregularity of H16. For a packet with one allowed input coordinate:
activation requires j>=13; extinction requires h<=4; a live-to-live
transition requires h+j>=13. These are necessary conditions only.

For a fixed pair of support masks, the equations for the next inactive
coordinates have matrix M=(H16 D1)[H'^c,H] in the incoming fresh diagonal.
Its rank is min(16-j,h). Choose that many pivot equations, dropping other
equations if necessary. The preceding rank argument charges
tau^(-min(16-j,h)). All active-state nonzero requirements remain in the
potential. The last outgoing support may be classified without revealing
or charging the unused last refresh diagonal.

For 0<z<1, define nonnegative 17 by 17 matrices E and I. Entries not
listed below are zero:

    E[0,0] = I[0,0] = 1,
    E[h,j] = binom(16,j) tau^(-(16-j)) z^w(h)
               for h,j>0 and h+j>=17,
    I[0,j] = binom(16,j) z                     for j>=13,
    I[h,0] = tau^(-h) z^(w(h)-1)               for 1<=h<=4,
    I[h,j] = binom(16,j) tau^(-min(16-j,h)) z^(w(h)-1)
               for h,j>0 and h+j>=13.

The destination-mask multiplicity is binom(16,j). Each row is a uniform
upper bound for every incoming mask of size h. In particular, aggregating
masks does not assume a symmetric distribution of actual states.

The regional operator is

    G = (1/128) sum_(k=0)^127 E^k I E^(127-k).

Let P have its sole nonzero entry P[0,0]=1. With only one possible input
per region, a path cannot activate and then return to zero in that same
region. Thus G[0,0]=1. Four consecutive regions form one outer column.
The column operator, including the deterministic quiet-column rank, is

    T = G^4 - P + tau^(-4) P.

The subtraction removes only the unique all-quiet subpath. Other paths
from zero to zero remain. For stable directed arithmetic, form G^4-P
by the positive recurrence R_1=G-P and R_(k+1)=R_k G+P(G-P).

For L=floor(delta*N), the union bound over all one-group supports is

    Pr[exists nonzero one-group message with output weight <= L]
      <= 8192 tau^128/(tau-1) z^(-L) e_0^T T^64 1.

For a fixed path, multiply the rank bound by z to its deterministic
output-weight floor. If a word of weight at most L realizes that path,
this marker is at least z^L. Sum paths and then average the independent
regional routes. This proves the displayed bound without counting all
F-valued messages individually. It holds pointwise in the fixed outer
mixers and D1 values, so averaging those setup choices preserves it.

### Directed arithmetic checkpoint

At tau=65280, z=19967/20000, and delta=7/20, the bound has base-two
logarithm approximately -183.693306321449. Directed evaluation at 192
bits proves it is strictly below 2^-180. The bound includes all 8192
choices of the active group and every nonzero F-message within it.
It does not cover messages with two or more active groups.

The bounded helper is [cheap_lane_sparse_gate.py](theory/cheap_lane_sparse_gate.py).
The older scalar-shadow mode remains a failed diagnostic, not a certificate.
Reproduce this checkpoint from the repository root with:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_sparse_gate.py --method density-free --verify --tau 65280 --z 0.99835 --distance 0.35 --precision 192 --target-bits 180
```

### Refinements close the same class at 40%

The dominant loss at 40% came from late activation. After 31 wholly quiet
columns, the coarse potential still leaves four source dimensions. Three
quiet regions in the next column leave only one. Ignoring those three
equations wastes almost 48 bits at tau=65280.

For each column, count its wholly quiet regions. Let c count columns with
four, and let m be the maximum count among the remaining columns; set
m=0 when none has a quiet region. If c<=31, the independent source rank
is at least 4c+m. Indeed, after imposing the c full columns, evaluation
at one further RS column is still surjective in all four constituents.
Any m rows of its invertible mixer are independent. This argument uses
one partial column; it does not add partial ranks across arbitrary columns.
Discard c>=32, which forces the source message to zero.

Two local refinements also reduce the bound. First, at an activation
packet the outgoing state is a nonzero scalar multiple of H16 D1 C_i.
For a fixed destination mask of size j, the event has probability at most
255^(-(16-j)) over the current D1. Superregularity proves this by fixing
nonpivot diagonal entries. Replace I[0,j] by

    binom(16,j) 255^(-(16-j)) z,  j>=13.

Second, extinction requires Q=C_i X_i with X_i nonzero, so the actual
incoming support must equal the support of C_i. In each sixteen-position
block, C_i has support {0} at point zero, {0,1} at point one, and all four
coordinates at the other fourteen points. Conditional on the occupied
packet, its coordinate is uniform among 64 positions. Thus multiply
I[h,0] by 1/64 for h=1,2, by 14/64 for h=4, and by zero otherwise.
These are maxima over incoming masks, not averages assuming random masks.

For these refinements, first fix the complete support path and occupied
packet indices. The rank and output-floor bounds are uniform in the
within-packet positions and D1 values. Activation and extinction events
belong to distinct regions, since each region has only one possible input.
Average those necessary setup events before aggregating masks. No
conditioning on the future message or quiet-kernel space is used.

Write G=P+B with B nonnegative and B[0,0]=0. The coefficients T_j in

    (B+vP)^4 = sum_(j=0)^4 v^j T_j

count column paths with exactly j wholly quiet regions; T_4=P. A positive
dynamic program tracks c in {0,...,31} and m in {0,...,3}. For j<=3,
replace m by max(m,j); for j=4, increment c. Charge tau^-4 on the latter
transition and tau^(-(max(m,j)-m)) on the former. The charges telescope
to tau^(-4c-m). At the end also multiply by

    1 - tau^(-(128-4c-m)),

so the initial quotient potential is bounded by the exact expression
P_tau(128-4c-m), rather than just its leading power. The prefactor remains
8192 tau^128/(tau-1) z^-L.

With tau=65280, z=623/625, and delta=2/5, directed arithmetic gives a
base-two bound approximately -44.29481729861079, strictly below -44.
This includes all nonzero messages in any one of the 8192 source groups.
It is not a certificate for the full code. Reproduce it with:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_sparse_gate.py --method density-free --verify --tau 65280 --z 0.9968 --distance 0.4 --activation-d1 --extinction-position --partial-quiet --precision 256 --target-bits 44
```

Floating-point proposals slightly underestimated this bound because small
intermediate contributions underflowed. Only the directed replay is used
for the stated margin; the helper's floating mode remains diagnostic.

## Extending the direct rank count to two through four groups

Fix an actual source set of size r<=4. A packet then contains at most
r allowed input positions. For l<=4 distinct positions, C_I has column
rank l: in each RS block, any at most four evaluation columns are
independent. Thus a zero-to-zero state transition forces every allowed
input value in that packet to zero.

Let B_l be the following seventeen-state cap, with B_0=E. For l>0,
B_l[0,0]=1. Activation requires j>=max(1,17-4l), and its output-weight
floor is one. Extinction requires h<=4l. A live-to-live transition
requires h+j>=17-4l. A live packet has output floor
max(0,13 ceil(h/4)-l), and pivot cost tau^-min(16-j,h). As before,
multiply by binom(16,j) for destination masks. Unlisted entries are zero.

There is a useful D1-averaged activation cap. Write P_256(l)=(256^l-1)/255.
If (H16 D1 C_I)[H'^c,:] has a nonzero kernel over F, its E-valued matrix
also has a nonzero kernel over E. Union over the P_256(l) projective
E inputs. For each input, C_I is injective and the usual diagonal-pivot
argument gives probability at most 255^(-(16-j)). Hence the activation
entry may be multiplied by

    min(1, P_256(l) * 255^(-(16-j))).

This event and the uniform rank-potential bound are combined before
averaging D1. It is not a claim that actual F input values are uniform
or belong to E.

### Exact averaging of the paired route

In one regional matching, let s count singleton active pairs and f count
fully active pairs, so s+2f=r. There are 4096 pair slots, thirty-two per
packet. Define a matrix polynomial

    A(x,y) = sum_(a<=s,b<=f) multinomial(32;32-a-b,a,b) x^a y^b B_(a+2b).

The regional operator is the coefficient of x^s y^f in A(x,y)^128,
divided by multinomial(4096;4096-s-f,s,f). Matrix multiplication retains
packet order. Pair orientations need not be counted separately because
the caps depend only on the allowed input count. Singleton orientation
factors cancel between numerator and denominator.

Each source support has fixed full-pair counts f_0 and f_1 in the two
alternating matchings. They satisfy f_0+f_1<r: the chosen vertices induce
a proper subgraph of the 8192-cycle, with at most r-1 edges. Maximize the
bound over all included profiles, then multiply by binom(8192,r). Do
not change the profile between columns.

Write each regional operator G_i=P+R_i, retaining R_i>=0 explicitly.
For r>=2, R_i[0,0] can be positive because a region can contain a live
excursion. Only P represents its wholly quiet path. Form the column
polynomial

    (R_0+vP)(R_1+vP)(R_0+vP)(R_1+vP) = sum_(q=0)^4 v^q T_q.

The helper uses positive quiet/nonquiet products throughout. It never
subtracts the full regional zero-to-zero entry.

### Retain constraints from the best thirty-two columns

Let q_j be the number of wholly quiet regions in column j. Each forces
r input values to zero. Any thirty-two outer columns can independently
realize arbitrary evaluations of the four RS constituents in each group.
Consequently the source rank is at least r times the sum of the largest
thirty-two q_j. The exact identity

    sum_largest32(q_j) = min_(v=0,...,4) [32v + sum_j max(q_j-v,0)]

gives four useful path classes. A nonzero witness has source rank below
128r, so a minimizing threshold exists with v<=3. Partition paths by
one such threshold before enlarging each class to all paths satisfying
e=sum_j max(q_j-v,0)<128-32v. Dropping this last cutoff reintroduces
zero-message paths and gives a useless bound.

For fixed v, propagate column coefficient T_q with counter increment
max(q-v,0) and factor tau^(-r max(q-v,0)). Start at state zero and
counter zero. Retain counters e<128-32v. The prefactor is
tau^(128r-32rv)/(tau-1); multiply terminal entries by
1-tau^(-r(128-32v-e)). Sum the four threshold bounds. Each threshold
may use its own fixed z. Finally apply z^-L within its threshold class
and the source-support union described above.

The bounded implementation is
[cheap_lane_few_group_gate.py](theory/cheap_lane_few_group_gate.py).
The next ledger checkpoint records its directed results; floating proposals
alone are not certificates.

### Further rank that the current operator discards

The quiet-region debit is not the whole deterministic rank. For any
packet, put M=R[H'^c,H], N=R[H'^c,:] C_I. Its left-kernel input equations
have rank d_t=rank[M,N]-rank M. They are E-valued and independent of all
fresh lambda coefficients. Within an outer column, different packets
and regions use disjoint input coordinates, so these ranks add.

Thus, with d_j the sum within column j, the same surjectivity argument
gives a source-rank lower bound sum_largest32(d_j). At h=0 this already
retains constraints imposed during activation, not only wholly quiet
regions. The five-group certificate below evaluates a deterministic lower
bound on these ranks. Using an averaged rank distribution would require
joint accounting with activation feasibility; multiplying their marginal
bounds is invalid.

## Five groups: retain rank even when quiet inputs cancel

Five allowed inputs can share an RS feedback block. A zero-to-zero
transition then need not have zero input. Let ell_b count allowed inputs
in block b, and put c=sum_b min(4,ell_b). The actual feedback matrix C_I
has rank c. For incoming and outgoing state sizes h and j, respectively,

    rank[M,N]-rank M >= max(0,c-h-j) =: d.

Indeed, deleting j rows from R C_I loses at most j ranks; adjoining the
h columns of M can account for at most h more. This lower bound holds
for every D1 realization. At h=j=0 it retains the full rank c, including
quiet packets inside otherwise active regions. The source-rank debit d
is therefore compatible with the separately averaged activation test.

For activation, union over projective vectors in im_E C_I, of dimension
c, rather than the whole input space. This excludes the feedback kernel
from the activation witness and gives min(1,P_256(c)255^(j-16)). At a
live packet, choose ceil(h/4) active emission blocks. Their output weight
is at least max(0,13 ceil(h/4)-sum_largest_ceil(h/4)(ell_b)).
These are conservative caps over actual state masks, not assumptions of
uniform state values.

The route calculation now retains the four occupancies. A half packet
has sixteen pair slots. With s singleton and f full pairs, u singletons
enter its first block, giving occupancies (f+u,f+s-u) and weight
multinomial(16;16-s-f,s,f) binom(s,u). Convolve two halves, and normalize
the regional coefficient by multinomial(4096;4096-s-f,s,f) 2^s. Thus
both pair positions and fair orientations follow the implemented route.

Sum d within each outer column to obtain d_j in {0,...,20}. The best
thirty-two columns give the source-rank lower bound. Partition paths by
a minimizing threshold v in {0,...,19}. Put D_v=640-32v and
e=sum_j(d_j-v)_+. A nonzero witness requires e<D_v. For v<19, the fixed
positive marker x=1/sqrt(2) bounds the truncated sum by

    sum_(e<D_v) w_e tau^-e
      <= x^(-(D_v-1)) sum_e w_e (x/tau)^e.

For v=19, only full-rank columns increment e; keep the exact counters
0,...,31 instead. The final factor 1-tau^(-(32-e)) retains P_tau exactly.
All twenty threshold bounds are added, and their Chernoff marks are
fixed before setup. Maximize over the nine possible fixed pair profiles
and multiply by binom(8192,5); a second profile union factor is unnecessary.

The [rank-debit helper](theory/cheap_lane_rank_debit_gate.py) implements
this operator. Its optimizer-free directed replay proves a union bound
below 2^-50 at 35% distance. The exact selected parameters and source
hash are in the goal ledger. This closes exactly five active groups,
not the intermediate-support range.

## A closed dense band by input-column deletion

For dense actual supports, losing a few source groups in the scalar
shadow can be charged without an exact joint route count. Fix an actual
group set T of size r>=7688. Use the initial comparison above with
h=32 missing groups and c_h>=1/33. If S omits j groups of T, deleting
their 256j input columns lowers the combined local matrix ranks by at
most 256j in total. Consequently the scalar kernels satisfy, pathwise,

    product_t K_t(I_t(S)) <= q^(256j) product_t K_t(I_t(T)).

After applying the outer density beta=q^160/255^32, the sum over missing
groups is at most beta^r C32, where

    C32 = sum_(j=0)^32 binom(8192,j) (q^256/beta)^j.

Padding I_t(T) to all sixty-four input positions costs another
q^(N-256r). Both comparisons hold for every route and D1 realization.
They do not replace the actual input distribution or require scalar
shadows to preserve state supports.

The full-input kernel has a simple bound. Substitute Y=X+A_H u in its
auxiliary scalar expectation. Since CA=0, the next pre-refresh state is
u+CY. Eliminating the requested output zeros gives

    rank[M,N_all] = |Z| + rank R_(H'^c)(C_(Z^c),I_H)
                 >= |Z| + min(16-j,h),

where h=|H| and j=|H'|. Superregularity is used only for H16, not the
coefficient-basis emission matrix A. With g(z)=z+(1-z)/q, summing the
positive zero-subset weights gives K_full<=g(z)^64 q^-min(16-j,h).
The final packet contributes g(z)^64 without a transition constraint.

For B[h,j]=binom(16,j)q^-min(16-j,h), the positive vector
v[h]=q^-h+q^-8 satisfies Bv<=rho v, where

    rho=(1+1/q)^16+65536 q^-8,
    v[0]/min_h v[h]=q^8.

This follows by bounding q^-min(a,h) by q^-a+q^-h and summing binomial
coefficients. Thus the sum over state paths is at most
g(z)^N q^8 rho^32767. The unused final refresh is not charged.

Set r0=7688, N=2097152, L=734003, and z=7/3315. Since beta/q^256<1,
union over all actual supports of size at least r0 gives

    U = 2^8192 * 33/255 * q^N * (beta/q^256)^r0 * C32
        * q^8 * rho^32767 * g(z)^N / z^L < 2^-98.

The [dense-band helper](theory/cheap_lane_dense_band.py) evaluates this
fixed expression with directed arithmetic. It covers arbitrary nonzero
F messages in all source-count classes 7688 through 8192, including full
lane rank. The intermediate classes still require a tighter calculation.

## Next local refinements, not yet a continuation certificate

The six-group diagnostic loses ranks during clustered activation and
extinction. The following activation bound charges rank and feasibility
together, instead of multiplying separate marginal bounds.

Fix an allowed input set I and an outgoing support H'. Put c=rank C_I>0,
m=16-|H'|, and let G be the restriction of H16 D1 to im_E C_I, followed
by projection onto H'^c. Write d=rank G. A nonzero activation requires
d<c; this condition is the same over E and F. For each projective E
vector in im C_I, the independent nonzero D1 entries give probability
at most 255^-m that its image under G vanishes. Therefore

    E_D1 P_q(c-d) <= P_q(c) 255^-m.

For tau>=q+1, the ratios P_q(n)/P_q(n-1)<=q+1, n>=2, imply

    E_D1[tau^-d 1[d<c]] <= 255^-m.

This replaces the earlier P_q(c) union factor only when the same local
rank d is charged in the source-rank weight. It cannot be added as an
independent improvement to the already averaged activation cap. For
c=0, nonzero activation is impossible.

A possible composition uses sum_largest32(d_j)>=sum_j d_j/2. With
s=sqrt(65280)<q and c>=1, the analogous bound is

    E_D1[s^-d 1[d<c]] <= Gamma_c(s) 255^-m,
    Gamma_c(s)=P_q(c)/s^(c-1).

The maximum ratio occurs at d=c-1. In particular Gamma_c<1.014 for c<=6
and <1.034 for c<=16. First apply the lambda potential pointwise in all
D1 values, including the deterministic source-rank bound. Only then
average the complete rank/feasibility factors over independent D1 values.
The exact source-rank cutoff still needs to be retained or replaced by
a justified positive bound. This local lemma alone does not close r=6.

At extinction, H' is empty, and invertibility of R gives the exact debit

    rank[M,N]-rank M = rank C_(H^c,I).

The current cardinality cap can miss this rank. The existing emission
minor checker supplies the corresponding transposed local matrix data;
a future operator must respect the paired route's common RS points.

For a scalar auxiliary return calculation, a stronger block bound is
also available. Let ell be the allowed input count and h the allowed
state-coordinate count in one four-state RS block. Set c=min(4,ell),
w=max(0,13-ell), and let k_ell be the exact quiet kernel from the ledger.
For independent uniform scalar x in E^I and u in E^H, define
K0=E[z^wt(x+A u) 1[u+Cx=0]]. Then

    K0 <= q^-h k_ell(z) [1+(q^min(h,c)-1)z^w].

To prove this, fix u=Cx. The occupied output coordinates form an affine
coset of ker C_I. Expanding the weighted mass by requested zero sets
bounds each affine coset by the homogeneous coset, whose normalized
mass is k_ell(z). There are at most q^min(h,c)-1 nonzero reachable u.
Each such u contributes at least w nonzero outputs outside I. Finally
average u, giving the factor q^-h. This argument concerns the auxiliary
scalar kernel, not the distribution of an actual field-valued state.

Both local lemmas passed independent mathematical review. The next section
integrates activation ranks into a closed six-group bound. The scalar
return refinement is not part of that certificate.

## Six groups: preserve activation ranks and the source cutoff

Charging half the total input rank gives a valid bound, but loses an
important exclusion: thirty-two full-rank outer columns already force
the entire source to zero. Even after removing physically quiet paths,
the half-rank screen retains too many nearly quiet histories. We therefore
return to the best-thirty-two-column cutoff and retain local rank tags.

For activation, keep the notation c,m,d above. On the event rank G=d,
the projective kernel count equals P_q(c-d). Consequently

    Pr_D1[rank G=d] <= min(1, P_q(c)255^-m/P_q(c-d)),  d<c.

Tag that same d as deterministic input rank in the packet polynomial.
There is no additional activation-union factor. The separate rank caps
may sum above one; replacing each rank-class mass in a positive sum is
still an upper bound. Rank c cannot activate and is omitted.

The four emission blocks strengthen the deterministic rank floor. Write
c_b=min(4,ell_b). A q_b-coordinate subset of block b intersects its
feedback image in dimension at least max(0,c_b-4+q_b). The MDS check
R[H'^c,:] is injective on any at most m state coordinates. Thus

    d >= max_(q_b in {0,...,4}, sum q_b<=m)
                  sum_b max(0,c_b-4+q_b).

This bound depends only on occupancies, so it preserves the exact paired
route average. For example, ell=(2,2,1,1) and m=5 give d>=2.

For live transitions with h<m, eliminate the h incoming-coordinate
columns from R[H'^c,:]. Their span does not depend on D1, because its
nonzero entries only rescale those columns. A fixed left-kernel basis
gives an (m-h)-row MDS check matrix on 16-h coordinates. Choose a fixed
k0=max(0,c-h)-dimensional subspace of im C[H^c,I] before revealing D1.
The projected rank of this subspace is a lower bound on the actual
input debit. Put a=m-h and n=16-h. Its rank-d class has mass at most

    min(1,P_q(k0)255^-a/P_q(k0-d))  for d<k0,

and mass at most one for d=k0. Restrict d to
max(0,a+k0-n)<=d<=min(a,k0). If h>=m or k0=0, use debit zero.
At extinction the projected check is square; the selected subspace
contributes exactly k0 ranks. These caps are averaged jointly with the
tagged ranks, after the pointwise lambda argument. They are not products
of marginal rank and feasibility estimates.

Compose these positive packet polynomials using the actual paired route.
Column debits range from zero to twenty-four. Apply all twenty-four
threshold classes, keeping the exact final thirty-two-counter cutoff.
The [joint-rank helper](theory/cheap_lane_joint_activation_gate.py) proves

    Pr[there is a message supported on exactly six outer groups
       whose encoded word has weight <=734003] < 2^-45.

The probability is over ideal setup, and messages range over the full
GF128 alphabet. The bound includes quiet histories, all fifteen fixed
matching profiles, and all binom(8192,6) source supports. Fixed 192- and
256-bit replays agree; an independent mathematical and source audit
passed. The selected parameters and hash are in the goal ledger.

### Next refinements retain geometry rather than add mixing

At extinction, rank and feasibility can be computed blockwise from the
actual RS point sets. Let I_b be the input points and H_b the incoming
state mask in one block. Put d_b=rank C[H_b^c,I_b]. A state with exact
support H_b is possible only if adding each row indexed by H_b increases
this rank. This tests whether every coordinate projection of
C[H_b,I_b] ker C[H_b^c,I_b] is nonzero. Since |E|>16, these individual
conditions also permit a vector with all the required coordinates nonzero.

A bounded enumeration illustrates why shared points matter. For three
full pairs split between the two half packets, ell=(2,2,1,1), and h=11,
the worst exact rank/feasibility factor is

    (112/120)(1/16) tau^-1.

Its state masks are (15,14,3,3), and all four emission blocks are active,
giving output floor 46 rather than the cardinality-only floor 34.
At tau=65280, the factor improves the old extinction cap by about 20.094 bits. A
four-way Holder relaxation improves it by only 2.006 bits because it
forgets that the two coordinates of a full pair use the same RS point.
This is a local finite calculation, not an additional class certificate.
Integration must preserve rank tags through the global cutoff.

For larger supports, output constraints offer a separate refinement.
For requested zeros Z, set p=|Z intersect I| and e=rank A[Z minus I,H].
The full output-and-feedback matrix satisfies rank[M,N_I]>=p+e and
rank M>=e. If s^2=tau and d=rank[M,N_I]-rank M, then

    tau^(-rank M) s^-d <= s^-p tau^-e.

This charges cancellations at occupied output positions instead of
using only an output-weight floor. It can be combined with a feedback
moment by taking the minimum of two upper bounds, not by multiplying
them. A useful whole-range operator still needs compatible source-rank
cutoff accounting. No claim for intermediate source counts follows
from this local inequality alone.

## A stronger dense bound from the output image

The previous dense bound padded every missing input coordinate at cost q.
A direct dimension argument avoids that padding. This argument concerns
the scalar auxiliary kernel; it does not change the actual alphabet or
the refresh distribution.

Fix I, H, the outgoing zero set J, and R=H16 D1. Put ell=|I|, h=|H|,
r_f=rank R[J,H]=min(|J|,h), and

    b=rank[R[J,H], R[J,:]C_I]=r_f+d.

For independent uniform x in E^I and u in E^H, the feedback kernel V
has dimension ell+h-b. The output is Y=x+A_H u. If Y=0, then x=A_H u
and CA=0 implies R[J,H]u=0. Therefore the kernel of the output map on V
has dimension at most h-r_f. Its output image has dimension at least
ell-d.

For any t-dimensional linear subspace W of E^64, choose t coordinates
on which projection is bijective. Uniform Y in W has independent uniform
entries on those coordinates. Dropping the other coordinates increases
z^wt(Y) for 0<z<1. With g=z+(1-z)/q, its normalized moment is at most g^t.
Consequently

    K_I(z) <= q^-b g^(ell-d)
           = q^-r_f g^ell (qg)^-d <= q^-r_f g^ell.

The final packet has no feedback condition and contributes at most g^ell.
This does not assert K_I<=g^(-missing)K_full for the exact kernels. That
inequality is false: on a quiet packet a single allowed input is forced
to zero by C, so deleting it can cost q rather than 1/g.

Apply the direct bound to each scalar shadow S of an actual source set T.
Every route has total allowed-input count 256|S|. With the preserved
outer density beta, write A0=beta*g^256. The source factor is A0^|S|.
The earlier comparison with at most thirty-two missing scalar groups
therefore gives, simultaneously for all actual supports of size r>=r0,

    U <= 2^8192 * 33/255 * A0^r0
         * sum_(j=0)^32 binom(8192,j) A0^-j
         * q^8 rho^32767 z^-L,

provided A0<1. The state-path bound and final boundary are unchanged.
The [fixed directed helper](theory/cheap_lane_dense_gpadding.py) uses
r0=7584 and z=3/1250. It proves U<2^-768 at L=734003, covering arbitrary
F messages in every actual source-count class 7584--8192. Independent
proof/source review and 192-/256-bit replays passed. The earlier dense
certificate remains frozen.

There is a stronger local inequality available for later use. Define

    e=rank([R[J,H]; A[I^c,H]])-r_f.

The zero-output kernel on V consists exactly of u satisfying both displayed
maps. Thus the output image has dimension exactly ell-d+e, and

    K_I(z) <= q^-r_f g^ell (qg)^-d g^e.

The current certificate discards the last two factors. Retaining them
could charge state-generated outputs at unoccupied positions. No new
class claim is based on this sharper expression yet.

## A compact way to retain the full-column exclusion

The intermediate-support calculation must exclude paths whose source
constraints already force the message to zero. The following bound keeps
one useful exclusion without storing every possible column rank.

Fix an actual source set T of size r. Each outer column contributes 4r
allowed input coordinates. Let d_j in [0,4r] be the certified source
rank in column j, for j=1,...,64, and put e_j=4r-d_j. The best thirty-two
columns give a source-rank lower bound. Since the source dimension is
128r, the remaining dimension is at most the sum of the smallest
thirty-two e_j.

A nonzero source therefore requires at least thirty-three columns with
e_j>0. Also, the sum of the smallest thirty-two deficits is at most half
their total. For tau>1 and any fixed marker y>=1, a valid continuation
upper bound is

    P_tau(remaining dimension)
      <= tau^(sum_j e_j/2)/(tau-1) * y^(#{j:e_j>0}-33).

For paths already excluded by source rank, the left side is zero. For
the other paths, both the half-total-rank bound and the marker inequality
apply. Increasing any source-rank tag can only decrease this right side.
Consequently positive lower-rank tags justified by weighted cumulative
bounds remain valid here.

For a packet with ell allowed inputs, source-rank tag d, and fresh-pivot
rank k, use the local factor tau^((ell-d)/2-k). Within each outer column,
retain a flag recording whether any packet has d<ell. Reset that flag
between columns and multiply deficient columns by y. Compute the two
flagged components with positive operations, not by subtracting a nearly
equal full-rank component from a total.

If weighted rank tags already contain tau^-k, multiply only by
tau^((ell-d)/2); do not charge the pivot rank twice.

This marker incorporates a necessary full-column condition, but does not
enforce it exactly at finite y. Even an exact deficient-column counter
fails the bounded six-group diagnostic. At z=.995, its floating log2
bound is about +5392. For example, thirty-one full columns and thirty-three
zero-rank columns have best-thirty-two rank 744, but half-total rank 372.
The missing partial-rank information is therefore material.

A sharper variant conditions on F full-rank columns, F=0,...,31. The
smallest thirty-two deficits sum to at most

    alpha_F sum_j e_j,  alpha_F=(32-F)/(64-F).

This follows by averaging the smallest 32-F of the 64-F positive deficits.
The combined bound is monotone under lowering ranks: a formerly zero
deficit increases both the total and alpha_F. Substitute one common tag
measure before partitioning by F. In particular, the Chernoff mark must
remain common across counter classes unless a new argument permits
different choices.

The adaptive variant also fails the six-group test: its best tested
floating log2 bound is about +5382. It repairs the thirty-one-full-column
example, but the dominant class has no fully ranked columns. Nearly full
columns still matter. Both diagnostics are failures of upper bounds, not
low-weight codewords. They do not justify changing the encoder.

A fixed middle threshold is another valid, but insufficient, shortcut.
Let c columns have rank at least v, with rank sums S_hi and S_lo above
and below v. Sorting gives the lower bound

    L_v = S_hi + (32-c)S_lo/(64-c)  if c<=32;
          32 S_hi/c                if c>=32.

The two expressions agree at c=32. L_v need not be monotone when a
column crosses the threshold. Therefore first substitute the common
weighted rank-tag measure against the exact monotone top-thirty-two
potential. Only then replace top-thirty-two rank by L_v, pointwise, and
partition by c. No counter class may retrospectively choose its own
local tags or output mark.

Tests at v=12 and 18 on the six-group profile (0,3) also fail. At
z=.995 their floating log2 bounds are about +5396; the zero-high-column
class alone contributes about +5381. Stable per-counter scaling was
used. This identifies a common-mark and rank-resolution limitation;
removing only the fully ranked ghost histories would not fix it.

## Positive output-rank tags for the intermediate regime

The next local operator charges requested output zeros as well as
feedback constraints. Fix I, H, H', and Z. Let k_f=rank M_fb, and let
d_f be a lower bound on the feedback input rank. Define

    p=|Z intersect I|,
    e=rank A[Z minus I,H],
    e_all=rank A[Z,H].

Write k=rank M and d=rank[M,N_I]-k for the full output-and-feedback
system. Adding rows cannot decrease this input-rank difference: any
increase in rank M also increases the joint rank. Thus d>=d_f. Eliminating
the occupied output identities gives output joint rank p+e. Consequently

    k>=k0=max(k_f,e_all),
    k+d>=R0=max(k_f+d_f,p+e),
    k<=h=|H|.

For each integer t in [0,|I|], the weighted point mass tau^-k at rank d
has cumulative mass at most

    F(t)=0                         if t<d_f or max(k0,R0-t)>h;
         tau^-max(k0,R0-t)         otherwise.

F is nondecreasing. Its nonnegative successive differences therefore
dominate the actual mass against every nonnegative, nonincreasing
continuation in d. This follows by expanding the continuation as a
positive sum of prefix indicators. The tags are weighted masses, not
probabilities; do not normalize them or charge tau^-k again.

F decreases pointwise as d_f increases. A lower feedback-rank measure
may therefore be averaged through it. The output ranks are independent
of fresh D1. Activation feasibility must remain in the activation measure.
The implementation currently tests only live incoming states.

Use the actual four-block emission ranks, or their conservative lower
bounds max(0,h_b-max(0,4-v_b)) for v_b requested zeros in block b.
These bounds are exact for h_b in {0,4}; they do not assume a superregular
coefficient basis. Sum output-zero subsets with their positive weights.
If an active block has four requested zeros outside I, omit that subset:
those equations would force the active state block to zero. Do not
renormalize the retained weights. This is a restriction on genuine state
supports, not on auxiliary scalar states in the dense kernel.

The old output-floor bound cannot be combined with these tags by a local
minimum. The two expressions bound different relaxed objects. For example,
an empty input and one active state coordinate can leave a tau^-1 residue
in the raw zero-subset sum; the floor-only expression has no such residue.
An early diagnostic merger was rejected in audit and removed before
freezing the helper. No class certificate used it.

A valid interpolation splits one fixed global mark z=a*b. If the actual
output has weight w>=f, then

    (ab)^w <= a^f b^w.

Use a^f for the known output floor and expand b^w by requested zero sets.
The global prefactor remains z^-L. The helper uses a=z^gamma and
b=z^(1-gamma); gamma=1 recovers the old tags, while gamma=0 gives the
new tags. A local choice may depend on I,H,H' before D1 and zero-set
averaging. Do not retrospectively optimize z for output-rank classes:
unlike the older feedback-only classes, these classes depend on the
requested zero sets themselves.

The [local helper](theory/cheap_lane_output_rank_tags.py) passed independent
mathematical/source review and exact tests on 1572 synthetic matrix/zero-set
cases. The tests include 17340 rational cumulative inequalities and 30156
monotone-continuation comparisons. They test algebra, not global distance.

For the unique h=16 state mask and j=16, the floating half-rank moment
at z=.5 improves by approximately 14, 16, 18, and 44 bits per packet for
occupancies (1,1,0,0), (2,2,0,0), (2,2,1,1), and (8,8,8,8), respectively.
At z=.8 the gains are approximately 4.507, 5.150, 5.794, and 14.164 bits.
Both screens use gamma=0 and the valid feasibility filter. Other outgoing
masks can favor the old gamma=1 endpoint. These are local comparisons,
not additive whole-code margins. The common-quotient source cutoff,
actual route average, support union, and terminal boundary still need
to be combined into a directed intermediate-range certificate.

### Factorized output constraints at the two extreme state supports

The intermediate pilot can retain output constraints cheaply when all
sixteen state coordinates are active, or when the incoming state is zero.
The following formulas use one common global mark z. They do not retune
that mark after observing a source-rank class.

First fix h=j=16. In one sixteen-coordinate emission block, let ell
coordinates belong to I. A requested zero subset chooses p of those
coordinates and v outside I. The emission map is RS[16,4], so its state
rank is k=min(4,p+v), and its joint rank is p+min(4,v). The input-rank
increment is therefore d=p+min(4,v)-k. An active state block cannot have
four requested zeros outside I. Its positive rank polynomial is

    F_ell(t) = sum_{p=0}^ell sum_{v=0}^{min(3,16-ell)}
      C(ell,p) C(16-ell,v) z^(16-p-v) (1-z)^(p+v)
      tau^-min(4,p+v) t^[p+min(4,v)-min(4,p+v)].

Multiply the four block polynomials. In the half-deficit moment, evaluate
at t=tau^-1/2 and multiply by tau^(sum ell/2). Every nonempty block has
d<ell, so a packet with any allowed input is deficient under these tags.
This is an exact local rank calculation after the stated feasibility
restriction; it is not a whole-code distance bound.

Next fix a quiet-to-quiet transition. The state is zero on both sides,
so the constraints are x_P=0 and Cx=0. In one block, p selected input
zeros give rank d=p+min(4,ell-p). Output zeros outside I impose nothing
and sum out. The half-deficit moment of that block is

    sum_{p=0}^ell C(ell,p) z^(ell-p) (1-z)^p
      tau^[(ell-p-min(4,ell-p))/2].

Its full-rank component is the same sum restricted to ell-p<=4. The
four block factors multiply. This retains the heavy-occupancy terms
without replacing their rank by min(4,ell).

For activation from h=0 to j>0, first fix the requested input-zero set
P. Eliminating its p coordinates leaves the feedback map restricted to
I minus P. Write ell'_b=ell_b-|P intersect block b| and
c'=sum_b min(4,ell'_b). Apply the existing joint activation/rank bound
to this restricted map. If its rank is d', total input rank is p+d'.
Activation requires d'<c'; terms with c'=0 vanish. Thus no activation
term has full input rank. Its half-deficit weight is
tau^[(sum_b ell'_b-d')/2]. The projective rank cap and structural floor
must both use the remaining occupancies ell', not the original ell.

The zero-set sum can be performed inside the route polynomial. Let u
mark a singleton pair and v mark a full pair. One pair initially has
polynomial 1+(u/2)(x+y)+v*x*y, where x and y mark the two residual input
positions. Replace each active input by (1-z)+z*x, or by (1-z)+z*y on
the other side. The polynomial becomes a+b*(x+y)+c*x*y, with

    a=1+u*(1-z)+v*(1-z)^2,
    b=u*z/2+v*z*(1-z),
    c=v*z^2.

Its sixteenth power gives one half packet. This substitution keeps the
original source markers and sums every requested-zero subset positively.
It does not introduce independent randomness or multiply separate
activation and rank probabilities. Two independent mathematical reviews
checked this activation reduction. Its effect on the global bound still
requires an actual-route calculation.

### An online rank-selection shortcut loses the exact cutoff

Selecting any thirty-two of the sixty-four column ranks gives a lower
bound on their largest-thirty-two sum. An online selection policy therefore
gives a valid upper bound on the remaining-dimension potential. However,
it need not select only full-rank columns even when thirty-two are present.
It may already have selected a deficient column before seeing later ranks.

For a small example, select two of four deficits. The first deficit is
one; the other three are independently zero or two with equal probability.
At tau=2, selecting the first deficit and then selecting a future zero
when available has expected exponential cost 2.75. Skipping the first
deficit has cost at least 4, even if the two best later deficits could be
chosen retrospectively. Thus a minimum-cost online policy selects the
first deficit. On half the paths, at least two later deficits are zero;
the best-two deficit is zero, but the online sum is positive.

Consequently an online rank-selection bound cannot replace the exact
full-column cutoff without an additional counter or another argument.
The earlier suggestion that selecting full columns whenever encountered
would preserve that cutoff was rejected before implementation.

## Actual-support scalar bounds at intermediate source counts

The failed half-rank pilot motivates a different local bound, not an
encoder change. Return to the normalized scalar kernel K defined above.
For one packet, let ell_b be the allowed inputs in emission block b,
ell=sum_b ell_b, h=|H|, j=|H'|, and c=sum_b min(4,ell_b). Define

    g=z+(1-z)/q,        r_f=min(16-j,h),
    d0=max(0,c-h-j),    a0=max(0,h-sum_b max(0,ell_b-12)),
    e0=max(0,a0-r_f).

Fix the actual D1 realization. The feedback constraints have joint rank
r_f+d. Conditional on those constraints, the output image has dimension
ell-d+e, where

    e=rank([R[H'^c,H]; A[I^c,H]])-r_f.

Deleting j feedback rows and then quotienting the h state columns gives
d>=d0. Each RS[16,4] block loses at most max(0,ell_b-12) rank after
restricting to output positions outside I, so e>=e0. A uniform vector
in a t-dimensional E-linear subspace has weight moment at most g^t:
choose t independent output coordinates and drop the remaining weight.
Consequently

    K <= q^-r_f g^ell (qg)^-d g^e
      <= q^-r_f g^ell (qg)^-d0 g^e0.

This bound is uniform over D1 and state masks of the specified sizes.
It does not require the auxiliary scalar state to be nonzero. The exact
quiet-to-quiet kernel and the exact full-state/no-feedback output kernel
can strengthen the corresponding entries. Bounds on this same scalar
kernel may be combined by a minimum; this does not license taking minima
of the different rank-tag measures discussed earlier.

### Preserve all 256 groups without missing-group padding

The elementary sparse comparison used 1-r/q and therefore stopped at
r=q. A hyperplane-union bound removes that endpoint loss. Fix the initial
E-defined kernel U of dimension d and an actual group set T of size
1<=r<=q. Discard a branch if any group projection vanishes identically
on U. Every remaining group-zero kernel is a proper subspace of U.

For d>=2, enlarge each such subspace to a hyperplane. Fix the first
hyperplane. It has q^(d-1) vectors; each further distinct hyperplane
contributes at most q^(d-1)-q^(d-2) vectors outside the first. Repeated
hyperplanes contribute nothing new. Therefore the number n_T of scalar
vectors nonzero on every group satisfies

    n_T >= q^(d-2) (q-1)(q-r+1),
    P_q(d) <= q^2 / [(q-1)^2(q-r+1)] * n_T.

For d=1, every group-zero kernel is {0}, and the same potential comparison
holds directly. For d=0 there is no retained nonzero witness. In
particular, at r=q=256 the multiplicative cost is only (256/255)^2.
The scalar support is exactly T; no missing-group padding is necessary.
Two independent reviews checked the hyperplane argument. Apply it only
to the initial E-defined kernel, before the outer-density relaxation.

For r>q, retain the earlier bounded-missing-group argument. If at most k
groups are missing, deleting their routed input columns gives the safe
comparison K(I(S))<=q^(256k) K(I(T)) over the whole stream. Thus, with
beta=q^160/255^32 and c_h=1-r/[q(h+1)]>0, a conservative source factor is

    beta^r / [(q-1)c_h]
      * sum_{k=0}^h C(r,k) (q^256/beta)^k.

The output kernel is evaluated on the actual set T, and its routing
correlations remain intact. At r=256 this padding bound is valid but
inferior to the hyperplane comparison. In either case, the outer factor,
source-support multiplicity, routing average, state paths, and final
packet must all be included before claiming a new covered class.

### Beyond256 groups: retain the ranks of initial group projections

The threshold256 is a limitation of the unrestricted scalar comparison,
not a demonstrated low-distance word. For example, let U=E^2 and use the
257 group projections x-ty for t in E, together with y. Place each scalar
in one coordinate of its128-coordinate group. Every E-valued vector
vanishes in at least one group. The F-valued vector(theta,1) does not.
Thus n_T=0 although the extended kernel has an exact-T witness. No finite
constant can give P_q(d)<=constant*n_T for every such initial kernel.
This is an abstract-kernel counterexample to the comparison, not a claim
that this kernel occurs in the proposed encoder's low-weight analysis.

A rank-sensitive comparison still covers many initial kernels. Fix one
fine-mask/output-zero term, its actual support T, and the initial E-kernel
U before revealing any refresh coefficients. Discard zero group
projections as above. Write e_g=dim_E pi_g(U), and let m count projections
with e_g=1. A uniform u in U vanishes on group g with probability q^-e_g.
Consequently,

    n_T >= q^d [1-m/q-(r-m)/q^2].

Whenever the bracket is positive, this compares P_q(d) to the same
exact-T scalar count. In particular, if m=0 the multiplier is
1/[(q-1)(1-r/q^2)]. No rank-two outer distribution is asserted.

If m>0, let c be the number of distinct kernels among the rank-one
group projections. Choose one such hyperplane H and restrict to U\H.
A different hyperplane occupies a1/q fraction of U\H. A subspace of
codimension at least two occupies at most a1/q^2 fraction: it either
lies in H or its intersection with H has codimension one inside it.
Repeated hyperplanes need be counted only once. The resulting bound is

    n_T >= q^d (1-1/q) [1-(c-1)/q-(r-m)/q^2],
    P_q(d) <= q^3 n_T /
              [(q-1)^2 (q^2+q-q*c-r+m)].

For q=256, r<=7583, and c<=228, the denominator's last factor is at
least69, since m>=c. Thus a uniform valid multiplier is
q^3/[(q-1)^2*69]. This statement is pointwise in the initial kernel;
it does not assume independent group projections.

To use it, put all m=0 terms and all m>0,c<=228 terms in the good class.
The same uniform multiplier also dominates the preceding m=0 bound.
Partition the positive potential/Chernoff sum into good and exceptional
terms before applying the fractional power.
Subadditivity permits separate bounds. On the first part, apply the
uniform comparison and then drop its rank-class indicator before the
existing outer-density and cone bounds. The local cone need not change.
This is a proved interface, not a computed certificate for any new weight.

Exceptional terms have at least229 distinct rank-one kernels and therefore
m>=229. Injectivity into the direct sum of group projections gives
d<=sum_g e_g<=128r-127m. This dimension debit alone does not close the
exceptional sum. A257-group pencil as above can be direct-summed with
unrestricted additional groups, leaving both large d and n_T=0.
Do not multiply q^(-127m) onto the old scalar bound: the rank condition
is correlated with the same source count already used in that bound.

Two scalar shadows provide another valid initial comparison. If ordered
pairs(u,v) in U^2 are counted whenever their group supports jointly equal
T, at least(1-r/q^2)q^(2d) pairs qualify. Hence P_q(d) is bounded by
the square root of that pair count divided by(q-1)sqrt(1-r/q^2).
This works for all8192 groups, but does not license squaring the scalar
outer density. After partitioning into S_1 union S_2=T, Cauchy--Schwarz
can bound the shared-outer expectation using two separate scalar moments.
The local actual-input filters must still use I(T); their two reference
kernels use I(S_1) and I(S_2). Applying the old exact-T cone separately
to the two shadows would discard possible actual paths. A joint
two-reference-support interface remains to be constructed and evaluated.

### Remove impossible actual paths before the scalar relaxation

Preserving T permits filters on the actual F-valued path before introducing
auxiliary scalar vectors. An empty actual packet cannot activate or
extinguish the state. In an active RS state block, four requested output
zeros outside I(T) are impossible. Remove those actual cases first; the
remaining positive scalar-kernel sum need not be an unfiltered scalar
weight moment. This does not assert that auxiliary scalar states are
nonzero. If shadows omit groups, retain I(T) in these filters while
padding their kernels; using I(S) instead would be invalid.

Activation has a useful joint filter. Fix h=0, j>0 and a requested
input-zero set P of size p. Write c' for the dimension of C on the
remaining input coordinates, m=16-j, and d' for the rank of R[H'^c,:]
on that image. Actual activation requires d'<c'. The projective D1
kernel bound gives

    E_D1[q^-d' * 1(d'<c')]
      <= min(1, [q/(q-1)](1-q^-c')255^-m).

Indeed, q^-d'1(d'<c') <= [q/(q-1)](q^-d'-q^-c'), and the expected
projective kernel size bounds the expectation of the difference. Terms
with c'=0 vanish. The eliminated coordinates contribute q^-p. Thus the
route substitution uses (1-z)/q for an eliminated input and z for a
retained input. This jointly charges activation and input constraints;
it does not multiply their marginal probabilities.

These filters should be included before investing in a larger coefficient
calculation. If scalar occupation markers still overcount crowded live
packets, a bounded next diagnostic can retain total input count exactly:
use the slot polynomial 1+(t/2)(x+y)+w*t^2*x*y, extract t^r over one
region, and leave the full-pair count as a positive w marker. This would
enforce the regional input budget with a single degree-r coefficient axis.
For a prescribed full-pair count f, divide by
multinomial(4096; r-2f,f,4096-r+f)*w^f. The slot polynomial already
averages the two singleton orientations; do not add a factor 2^(r-2f).
It is a proposed computation, not a completed bound.

### Filtered output moments and deterministic output floors

The new scalar pilot implements the actual-path filters before counting
auxiliary vectors. For one emission block, fix an actual state support of
size h_b>0 and an allowed input set of size ell. Put n=16-ell and
g=z+(1-z)/q. At most three positions outside that input set can be zero
on an actual nonzero RS[16,4] word. Define the retained zero-set mass

    R_n(z)=sum_{v=0}^{min(3,n)} C(n,v) z^(n-v)(1-z)^v.

The filtered output-only scalar kernel is at most g^ell times each of

    sum_{v=0}^{min(3,n)} C(n,v)z^(n-v)(1-z)^v
                         q^(-max(0,h_b-4+v)),
    q^(-h_b)R_n(z)+(1-q^(-h_b))z^max(0,n-3).

For the first bound, deleting 4-h_b columns from an independent set of
v<=3 rows loses at most 4-h_b rank. For the second, separate the zero
auxiliary state, of probability q^(-h_b), from nonzero auxiliary states.
The latter give nonzero RS words and therefore at most three roots.
These arguments do not assume that the auxiliary state has its claimed
actual support.

For h_b=1, an actual nonzero state can vanish only at roots of its basis
column. The implemented basis columns have at most two roots. Retaining
only those zero subsets gives the additional bound
g^ell*z^max(0,n-2), without an auxiliary-zero residue. This is a filter
on actual paths, not an assertion about an unfiltered scalar moment.

Take products over the four emission blocks and maximize over their
state-support sizes with sum h. The same retained weights can multiply
a feedback-only bound: each requested-zero kernel is at most the
feedback kernel before averaging D1. The pilot takes the minimum of
this bound, the filtered output-only bound, and the output-image bound.
All three bound the same filtered scalar object.

It is also important to retain the deterministic output floor. Fix
0<a,b<=1 and the single global mark z=a*b. On an actual path with packet weight w>=f,
z^w<=a^f*b^w. Apply this inequality before the scalar relaxation, and
expand only b^w into zero subsets. The pilot uses a=z^gamma and
b=z^(1-gamma). An active block contributes floor max(0,13-ell); a
zero-to-active transition contributes floor one. At gamma=1, the local
bound therefore charges the actual output floor together with the
feedback bound. Multiplying that floor into an already expanded
z-kernel would instead charge output weight twice and is not used.
The global output prefactor remains z^(-alpha L).

### Separate coefficient bounds and activity counts

The normalized route coefficient is fixed by its matching profile,
incoming and outgoing state supports, and activity flag. Different
positive coefficient markers can upper-bound this same coefficient.
Their pointwise minimum is valid before taking a profile envelope.
In particular, markers chosen for quiet paths need not equal markers
chosen for long active paths. These are different upper bounds on one
coefficient, not minima of incompatible rank measures.

A further disjoint split counts packets with nonzero incoming state.
Let k be that count, including the final packet. The initial state is
zero, so 0<=k<=32767. For an integer cutoff B, bound k<=B and k>=B+1
separately; each class may choose its own fixed z and alpha. After the
fractional-moment step and state-path subadditivity, introduce a positive
activity marker v. Multiply each live incoming row and the live terminal
entry by v, outside the fractional power of the scalar kernel. Then

    1{k<=B}   <= v^(k-B)       for 0<v<1,
    1{k>=B+1} <= v^(k-B-1)     for v>1.

The resulting bounds pay v^(-B) and v^(-B-1), respectively. For the
few-active class, z=1 is allowed: it bounds existence of such paths
regardless of output weight. No count axis is required for this
relaxation. A diagnostic using a few matching profiles is still not a
bound for all profiles or a certificate for an activity class.

Late first activity admits another audited comparison. Number outer
columns from zero. A first nonzero incoming packet state in column 63
requires C X_t=0 through packet 32254. Activation can occur in the
previous column's last update. Relative to the normalized scalar quiet
kernel imposing feedback at all 32768 packets, deleting the remaining
513 packet constraints costs at most q^(16*513), before alpha. If the
reference already omits final feedback, only 512 deletions are needed.
Apply the source comparison to the prefix event first, then use this
inequality between positive scalar kernels. It is not valid to multiply
an actual quiet-event probability by that loss: the full quiet kernel
can have no nonzero witness while the shortened prefix does.

### Next diagnostic: exact coefficients for pair-only profiles

The marker pilots still fail after these refinements. A bounded next
test can remove the input-count relaxation for one important profile,
without changing any local kernel. Suppose r is even and all active
groups are paired, so the full-pair count is f=r/2 and there are no
singleton pairs. One packet has occupancies (a,a,b,b), with 0<=a,b<=16.
Let M_(a,a,b,b) be its current local matrix majorant, including fractional
powers, destination multiplicities, and the positive activity flag.
Let T_(a,a,b,b) be its emission-only terminal column. Define

    P_d = sum_{a+b=d} C(16,a)C(16,b) M_(a,a,b,b),
    T_d = sum_{a+b=d} C(16,a)C(16,b) T_(a,a,b,b),
    P(t)=sum_d P_d t^d,          T(t)=sum_d T_d t^d.

The exact regional route averages of these majorants are

    [t^f] P(t)^128 / C(4096,f),
    [t^f] P(t)^127 T(t) / C(4096,f)

for an ordinary and final region, respectively. No swap multiplicity
is present: swapping a full pair does not change its occupancies.
The coefficient degree counts full pairs, not records. This computes
an exact route average of an upper bound, not an exact encoder failure
probability. It covers only the chosen pair-only profile.

Use positive coefficient convolution, truncate degrees above f, and
preserve small coefficients with per-degree scaling or logarithmic
arithmetic. FFT cancellation or underflow of quiet terms could falsely
suggest closure. Seven doubling steps can compute both expressions:
starting with A_1=P and B_1=T, use A_(2m)=A_m A_m and
B_(2m)=A_m B_m. The
[pair-only helper](theory/cheap_lane_paired_coefficient.py) now implements
this positive logarithmic convolution and checks it against small exact
slot enumerations. Its arithmetic is floating, not a directed certificate.

There is an important joint-profile limitation. A set fully paired in
both adjacent matchings must contain either every group or no group:
the two matchings together form a connected cycle. For 0<r<8192,
repeating f=r/2 through both matchings is therefore an artificial
diagnostic, not a realizable source class or a proved worst-case bound.
At r=256, if the first matching has 128 full pairs, the second has at
most 127. The profile (128,127) consists of 4096 aligned contiguous
supports; reversing the matchings gives (127,128). The second matching has two
singletons; its coefficient can be evaluated later with singleton degree
truncated at two. The current helper does not evaluate that profile.

### A full-column exclusion does not fix the tested restart terms

A regular packet with H=H'=empty and rank C_I=|I| forces X=0. Call
an outer column good when all its packets satisfy those conditions.
For the output-only terminal packet, require H=empty and I=empty
instead. Thirty-two good outer columns give thirty-two zero evaluations
of every RS[64,32] constituent, forcing the source to zero. Thus a
nonzero path has at least 33 columns that are not good. This is a valid
actual-path exclusion before scalar relaxation.

It is not sufficient for the tested scalar bound. Let N_00 be the
positive regional coefficient for paths starting and ending at zero
state but containing a nonzero incoming state. At r=256, z=.85,
alpha=.523, gamma=1, the exact-degree calculation gives
log2 N_00 approximately -429.319469. The source/output prefactor has
log2 approximately 228771.063269. The output-only terminal coefficient
is at least N_00: every regular last-step entry ending at zero is at
most the corresponding output-only entry, and all operations are positive.
Consequently the artificial repeated-marginal upper expression contains
a contribution with log2 at least

    228771.063269 + 256*(-429.319469) = 118865.279205.

Every region in that contribution contains activity, so all 64 columns
are not good. The proposed exclusion retains the contribution, even
with an exact column counter. These numerical values diagnose the
positive upper expression; they are not probabilities of code failures.
Independent review checked the terminal inequality and the exclusion.

### Joint output and feedback bounds for every nonzero outgoing state

One remaining local loss is avoidable. The auxiliary scalar feedback
kernel includes vectors with zero pre-refresh state, even on a path
whose actual outgoing state is nonzero. Extend the earlier activation
filter to every such transition, including live-to-live transitions.

Fix I, H, H', and a retained requested-zero set Z. Let j=|H'|>0,
J=H'^c, and m=|J|. Before sampling the current D1, define the E-linear
output-equation space

    U_Z = {(x,u) in E^I x E^H : (x+A_H u)_Z=0}.

Let r_out be the rank of these equations. Include u in E^16 by zero
padding outside H, and let G be the image of U_Z under V=Cx+u.
Write g=dim_E G and d=rank_E((H16 D1)[J,:]|G). Both U_Z and G are
independent of the current D1. Uniform auxiliary (x,u) conditioned on
U_Z gives uniform V in G, so the joint scalar mass is

    q^(-r_out) q^(-d).

An actual nonzero outgoing state requires a nonzero V in the scalar
extension of G that is killed by the rows J. Therefore d<g is necessary;
rank is unchanged by extension from E to F. If g=0, the actual branch
is impossible. Remove these impossible D1 outcomes before applying the
scalar relaxation.

For a fixed nonzero v in G, superregularity of H16 gives
Pr_D1[(H16 D1 v)_J=0]<=255^(-m). If wt(v)<=m, the probability is zero;
otherwise condition on nonpivot diagonal entries and solve the m pivot
entries. Summing over nonzero vectors of G gives

    E_D1[q^(-d)-q^(-g)] <= (1-q^(-g))255^(-m).

For g>0, q^(-d)1{d<g} is at most
q/(q-1)*(q^(-d)-q^(-g)). Thus the newly filtered joint kernel satisfies

    K_Z := q^(-r_out) q^(-d) 1{d<g},
    E_D1[K_Z] <= q^(-r_out) F_j,
    F_j := min(1, [q/(q-1)]255^(-(16-j))).

Set K_Z=0 when g=0.

Sum over retained Z with their positive output weights. Their earlier
actual-state filters and deterministic output floors are independent of
D1. If K_out denotes the resulting output-only scalar kernel, then

    K_joint_filtered <= F_j K_out,             j>0.

This is a joint rank calculation, not multiplication of independent
marginal probabilities. Take its minimum with the earlier upper bounds
on the same further-filtered kernel. Do not apply F_j to extinction
j=0 or to the output-only terminal. At j=16 the coarse factor is one.
Two independent mathematical audits checked this argument.

The terminal table upper-bounds K_out^alpha and includes the incoming
activity mark. Multiply its entry by F_j^alpha and C(16,j), keeping the
destination multiplicity outside alpha and the activity mark only once.
For h=0 the terminal entry omits the activation floor; this is weaker
but valid, and the earlier activation cap remains available. Choosing gamma
deterministically from {0,.5,1} for each I,h,j before D1 and the zero-set
expansion is also valid. It retains one common z and its original global
prefactor. It does not permit choices based on sampled D1 or Z, or
minima between incompatible rank-tag measures.

### Enforcing auxiliary support with an explicit local loss

The joint live filter removes a zero pre-refresh image. A stronger
comparison can require every prescribed auxiliary state coordinate to
be nonzero. This does not follow by conditioning the auxiliary vector:
the normalization changes, and an explicit comparison factor is needed.

Fix I,H,H',D1,Z. In the exact-source-support reduction, I is the allowed
input set of the actual source class. Let U be the kernel over E of
the combined output-zero and inactive-feedback equations on (x,u), with
x in E^I and u in E^H. Let rho be their joint rank. Choose the h
incoming coordinate forms u_i, i in H, and the j outgoing coordinate
forms (H16 D1(Cx+u))_k, k in H'. Put p=h+j. Define A to mean that
none of these p forms vanishes identically on U.

An actual F-valued witness implies A. All these equations and forms are
over E, and a form vanishing on U also vanishes on its scalar extension
to F. Therefore A is a valid actual-path filter before the scalar
relaxation. This use of I would not be valid in the earlier ambient
reduction, where auxiliary and actual source supports could differ.

When A holds, each of the p zero sets is a proper hyperplane of U.
For p<q, their union occupies at most p/q of U. Consequently

    1_A q^(-rho)
      <= q/(q-p) Pr_{x,u uniform in E^I x E^H}[
           combined equations hold,
           supp(u)=H, supp(H16 D1(Cx+u))=H'].

Coordinates outside H are zero by the domain, and inactive outgoing
coordinates are zero by the equations. The probability is normalized
by q^(|I|+h), not q^|I| (q-1)^h. This accounts for the rank-potential
normalization without an additional factor. If U={0} and p>0, A is
false and both sides vanish. If p=0, no filter or loss is needed.

Here p<=32 and q=256, so the worst local factor is 8/7 before the
fractional power. Extinction j=0 and the output-only terminal use h
forms, giving factor q/(q-h). For the terminal, U contains only output
equations. The comparison holds pointwise in D1 and Z. Apply the
feasibility filter, then the comparison; sum positive zero-subset weights,
average D1, and only then take the fractional power as above. Auxiliary
vectors in successive packets remain separate counting variables.
No assertion about the actual conditional state distribution is made.

A cheaper variant chooses one incoming coordinate from each active
emission block. When j>0, also choose one outgoing coordinate from H';
when j=0, choose none. Make these choices deterministically from the
prescribed masks before averaging D1 or Z. At most five forms then
force every active auxiliary block, and any prescribed nonzero outgoing
image, to be nonzero, at loss at most 256/251. Other coordinates may vanish.
This variant and the exact-support variant have different kernels and
must not be conflated. Both comparisons passed independent mathematical
review. Neither has yet been integrated into a whole-code calculation.

The [support-filter checker](theory/cheap_lane_support_filter.py) tests
883 exact comparisons over all specified small-field kernels and form
multisets, plus extinction, terminal, and duplicate-form examples. It
checks the zero-dimensional case and rejects p>=q. These tests confirm
finite arithmetic examples; the hyperplane argument proves the general
statement. The helper is not a whole-code certificate.

There is also a precise normalization bridge to a baseline-style scalar
transition. Let B_(I,H,H') use the same uniform x in E^I, but uniform u
in (E*)^H, with the same output weight and exact outgoing support event.
The exact-support mass above equals ((q-1)/q)^h B_(I,H,H'). Thus the
comparison factor relative to B is

    [q/(q-h-j)] ((q-1)/q)^h,

before alpha; replace j by zero for the terminal. This identity does not
transfer the old certificate. Its input distribution, outer comparison,
and route conditioning must match as well. If a baseline transition
requires exact nonzero input support S, first mix over S subset I with
weights q^(-|I|)(q-1)^|S|. The partial block-support variant does not
have this exact-support normalization.

### A cheap block-support cap to test next

The partial-support comparison permits a joint output/feedback cap without
enumerating all nonzero scalar states. Fix the actual block supports and
let b be the number of active emission blocks. Let ell_b be the allowed
input count in block b, ell their total, and

    f_out := sum over active blocks of max(0,13-ell_b),
    g(z) := z+(1-z)/q.

Select one incoming coordinate in each active block. For j>0, also
select one outgoing coordinate. For any auxiliary u satisfying the
incoming selections, each active block emits a nonzero RS word. Outside
I these words have total weight at least f_out. Conditional on u, uniform
x in E^I makes Y_I uniform, so its output moment is g(z)^ell. For every
fixed x,u, the D1 probability of the inactive outgoing equations and
the selected nonzero outgoing coordinate is at most 255^(-(16-j)).
It is zero when Cx+u=0. This bound is uniform in x,u and can be applied
inside their weighted expectation.

The selected incoming coordinates are independent under the auxiliary
uniform law. Applying the support comparison gives, for j>0,

    K_filtered <= [q/(q-b-1)] ((q-1)/q)^b
                   *255^(-(16-j))*g(z)^ell*z^f_out.

For the terminal, omit the outgoing selection and feedback probability,
and replace q/(q-b-1) by q/(q-b). This is not an extinction bound.
The floor counts only positions outside I; g(z)^ell counts positions
inside I. Thus output weight is not charged twice. No gamma split is
needed for this bound. Maximize over compatible block supports before
using a cardinality-only table, and apply alpha only after averaging D1.
Two independent mathematical reviews checked the joint probability and
the separation of inside-I and outside-I output weight.
The [bounded wrapper](theory/cheap_lane_block_support.py) tested this cap
after the previous gamma and joint-live refinements. At r=256, z=.85,
alpha=.523, its largest local improvement is only 0.0000435 bits and
the regional/chain diagnostics are unchanged to displayed precision.
This coarse cap is effectively dominated by the existing bounds. The
tenth-round ledger records the subsequent extinction test and its limits.

For extinction, exact local geometry can also help. With one allowed
input point a in one block, write v(a)=A[a,:], so Cx=v(a)x. Nonzero
extinction requires H=supp(v(a)), and the emitted word is
(I+AC)e_a*x. Using the implemented normalized basis
[1,a,a(a+1)/6,a^2(a+1)/6] over AES GF256 gives:

| Points a | Incoming mask | Weight of (I+AC)e_a |
|---|---:|---:|
| 0 | 1 | 15 |
| 1 | 3 | 16 |
| 2,3 | 15 | 15 |
| 4,...,15 | 15 | 16 |

Masks use the four basis coordinates as bits. The main agent computed
these sixteen columns exactly from the frozen emission helper. The
points are the integers 0,...,15 in AES GF256, not a GF16 subfield.
An independent algebraic check uses
(AC)_(t,a)=(1+at)(1+b(a)b(t)), where b(a)=a(a+1)/6. The resulting
I+AC has exactly three zero entries: (0,0), (2,3), and (3,2).
For the uniform auxiliary x,u, the exact-support extinction mass at a
fixed point is (1-1/q)q^(-|H|)z^weight when its mask matches, and zero
otherwise. Under a uniform point choice, average the displayed terms
once. Two blocks of a full pair share their point; multiplying independent
point averages would be invalid. The regional integration below uses the
sharper filtered formula, rather than this exact-support auxiliary mass.

A sharper form avoids even the comparison-factor loss here. For active
H, extinction forces x!=0. A requested zero at a nonzero coordinate
of (I+AC)e_a is impossible for the actual path. Remove such requests
before relaxation. The surviving requested-zero sum then equals the
fixed monomial z^weight, and the scalar feedback mass is q^(-|H|).
The resulting local cap is

    1{H=supp(v(a))} q^(-|H|) z^wt((I+AC)e_a).

For H empty, extinction instead forces x=0 and the scalar cap is q^-1.
The active formula has no auxiliary-zero output residue. Integrating it
requires the actual shared point geometry before mask aggregation.

### Integrating singleton extinction without breaking shared points

Fix the complete state-mask path as a summation index before averaging
routing. Do not condition the route on observing that path. Conditional
on the full-pair counts in each half-packet, their occupied point subsets
are independent and uniform. The two blocks within one half share the
same subset. Thus average each fixed-mask product raised to alpha over
the shared subset before maximizing over masks of a given cardinality.

Write Q=q^(-alpha) and p=z^alpha. For a half containing one full pair,
let T_1[s] be the maximum over its two fixed incoming masks whose total
size is s. The exact shared-point calculation gives

    T_1[0] = Q^2,
    T_1[1] = Q^2 p^15 /16,
    T_1[2] = max(Q^2 p^30, Q^3 p^16)/16,
    T_1[4] = max(Q^4 p^32, Q^5(2p^15+12p^16))/16,
    T_1[8] = Q^8(2p^30+12p^32)/16,

with all other entries zero. A full pair shares a position, not a
nonzero amplitude: one block may be quiet while the other is active.
There is no factor two for exchanging the active block, since this
operation takes a maximum over fixed masks, not another mask union.

For ell>=2, the first implementation keeps the coarse active-block cap
q^(-max(h,min(4,ell))) z^max(0,13-ell). A quiet block uses the exact
positive moment

    sum_{k=0}^ell C(ell,k) z^k ((1-z)/q)^(ell-k) q^(-min(4,k)).

The latter follows by expanding requested input zeros and using the MDS
rank of the remaining columns of C. For ell=0 only a quiet block can
extinguish. These expressions define T_ell[s] by maximizing the powered
product of two block caps. Independence between the two half-packets
then gives packet cap max_s T_a[s]T_b[h-s].

Take its minimum with the old uniform cardinality cap after the conditional
point average. Apply incoming activity once; no destination multiplicity
is needed for the single zero destination. Rebuild the live-to-zero
coefficient positively with weights C(16,a)C(16,b). Do not subtract an
old coefficient or average the two paired point sets independently.
The quiet-to-quiet entries, live destinations, and terminal are unchanged.
Two independent mathematical reviews checked this order.

The [singleton wrapper](theory/cheap_lane_singleton_extinction.py) tests
the implemented AES-GF256 basis, the shared-point formula at alpha=1
with exact rationals, and the powered formula at alpha=.523. The selected
regional diagnostic improves, but remains far from closure; see the ledger.

### Two-input extinction can retain output directions

For two distinct input points I, eliminate the state using u=C_I x.
For a fixed incoming mask H, put U=ker C[H^c,I] within E^2, and let
d=dim U. Discard a branch if any required coordinate of C[H,I] vanishes
identically on U. This is necessary for an actual F-valued witness.
The zero incoming mask forces x=0 and has scalar cap q^-2.

For active H and d=0, the branch is impossible. For d=1, take a nonzero
E-vector v spanning U. All nonzero actual inputs are F multiples of v,
so the exact filtered rank cap is q^(-h-1) z^wt((I+AC)_I v).

For d=2, the two output columns have no common zero row in the implemented
basis. Group their sixteen nonzero row forms by projective kernel line v.
Let n_v count rows in a class. Retain a line only when every required
coordinate of C[H,I]v is nonzero. Rank-two requested-zero subsets force
x=0 and are impossible; rank-one subsets all use a single retained class.
The positive filtered kernel is therefore

    q^-h [z^16 + q^-1 sum_{retained v} sum_{s=1}^{n_v}
                       C(n_v,s) z^(16-s)(1-z)^s].

Only H=15 with I different from {0,1}, or H=3 with I={0,1}, can have
a feasible d=2 branch. Other feasible masks use the d=1 formula. A
prospective calculation needs only 120 point pairs and at most sixteen
row-direction classes per pair. It must average both blocks over the
same I, after assembling and powering their complete positive kernels.
This extension is derived but not yet implemented.

### A two-packet filter from the degree of theta

Independent scalar states discard the relation between an activation and
its subsequent return to zero. That relation gives a different local
exclusion, without assuming the input lies in the byte field.

Consider two consecutive interior packets with states zero, Q, zero,
where Q is nonzero. Let W_i=im_E C_(I_i) for their allowed input sets,
and let B_0 be a basis matrix for W_0, of dimension c_0. Choose an
E-matrix P with kernel W_1. For the activation refresh, write D_b for
the nonzero b diagonal and D_a for the byte a diagonal. Define

    L := P D_b H16 D1 B_0,
    T(Tvar) := P D_b (Tvar I+D_a) H16 D1 B_0.

If L has rank k<16, a nonzero k-minor is the leading coefficient of
the matching minor of T(Tvar). That polynomial has degree k, so it
cannot vanish at theta, whose degree over E is sixteen. Thus T(theta)
has rank at least k over F for every a. In particular, if k=c_0,
no nonzero activated image can lie in W_1 extended to F, and the
two-packet excursion is impossible. The c_0=k=16 case requires W_1=0
and is separately excluded by invertibility of the refresh.

This filter depends on b. Reveal b and D1 before averaging a; conditional
on b, each b(theta+a) still has 256 distinct possible values. The current
q=256 rank potential remains valid, but a 65280-point contraction would
not survive this conditioning. Keep the filter's D1 dependence joint
with the local kernel, rather than multiplying independently averaged
bounds. A partial rank k is not an extra source debit to add to existing
pivot charges; that would require joint elimination. The full-rank case
is an actual-path exclusion and avoids this double counting.

Apply it only to interior returns already in the tracked state path,
not to an imposed return after the final output. Two independent audits
checked this argument. It is not yet a numerical class bound.

For fixed route images, a coarse bound on failure of this full-rank
test is available. Let c_1=dim W_1 and P_q(c_0)=(q^c_0-1)/(q-1).
Union over the E-projective directions of W_0. For a fixed nonzero
direction u, put v=H16 D1 u. The probability of any exact support S
of size g is at most 255^(-(16-g)). Conditional on v, D_b v is uniform
over the nonzero coordinates of S. At most q^min(c_1,g)-1 vectors of
W_1 lie there. Summing supports and directions gives

    Pr[rank L<c_0]
      <= min(1, P_q(c_0) 255^-16
                    sum_{g=1}^{16} C(16,g)(q^min(c_1,g)-1)).

The sum can be shortened when the input image has a uniform support
bound, but that improvement is not used here. This is a probability
over the activation packet's D1 and b, conditional on fixed route images.
It is not a free multiplier for a previously averaged scalar kernel.
For example, Holder gives the compatible joint bound
E[1_B K^alpha] <= Pr(B)^(1-alpha)(E K)^alpha for a defect event B.
Other uses must retain the same conditional dependencies.

### Retaining nonzero RS output moments on live paths

The simpler block-support cap uses only a minimum output weight. Replace
that floor by a positive moment that excludes the auxiliary zero block.
For a block with h>0 incoming coordinates, let n be its number of output
positions outside I. For a fixed requested-zero set Z in those positions,
the normalized mass of nonzero auxiliary states is

    q^(-rank A[Z,H])-q^-h.

At least four distinct requested zeros force the state to zero. For
v=|Z|<=3, deleting 4-h columns from the full RS map gives rank at least
r_lower(h,v)=max(0,h-4+v). Consequently the nonzero-block moment is at
most

    M_h(n,z) := q^-h sum_{v=0}^{min(3,n)} C(n,v)
                          z^(n-v)(1-z)^v (q^(h-r_lower(h,v))-1).

All terms are positive. For h=4, this is the exact moment for every
choice of the n outside positions. For h=1, also take the minimum with
(1-q^-1)z^max(0,n-2), since the implemented basis columns have at most
two roots. For h=2,3, the valid additional cap is
(1-q^-h)z^max(0,n-3). Define M_0(n,z)=1.

To justify this replacement in the F-valued proof, use zero-block
subspaces in the local support comparison. Each actual active block has
a nontrivial projection on the constraint kernel U. The kernel of that
projection is proper, with mass at most 1/q. For b active blocks and one
selected outgoing coordinate, the same union comparison therefore costs
q/(q-b-1). It no longer singles out an incoming scalar coordinate.

For a fixed block-mask allocation h_1,...,h_4 with total h, the resulting
joint live cap is

    [q/(q-b-1)] 255^(-(16-j)) g(z)^ell
                       product_{i=1}^4 M_(h_i)(16-ell_i,z),   j>0.

The conditional D1 bound is uniform in x,u, and dropping that event
restores the uniform translated output on I. The M terms already include
the mass of nonzero blocks: do not multiply by another such factor.
At the terminal, omit feedback and use q/(q-b). For j=16 the same
output-only comparison is also an upper bound and is slightly sharper.
Maximize the complete expression over compatible block-mask allocations,
then apply alpha and destination multiplicities. Do not reuse the new
terminal cap inside the older joint-live calculation, whose output
kernel has a different restriction. Two independent reviews checked
this replacement, including the conditional averaging order.

The bounded implementation in `theory/cheap_lane_nonzero_moment.py`
passes exact GF8 comparisons for all nonempty coordinate masks and all
puncturings of RS[8,4]. At r=256, z=.85, alpha=.523, the old local caps
already dominate almost all of this information. The regional and chain
diagnostics are unchanged at displayed precision. The tenth-round ledger
records the frozen source and numerical replay; this local lemma does
not by itself improve the certified support coverage.

### Two live packets: keep the state subspaces at the cut

The leading-theta exclusion need not wait for a return to zero. Requested
output zeros constrain the possible states on each side of a live cut.
Keeping these subspaces can exclude a sustained low-weight path that a
support-cardinality table permits.

Fix the route, state-mask path, requested output-zero sets, and D1
matrices. These are indices or fixed realizations, not conditioning on
an observed path. Let H be the nonzero state support at a cut, with
h=|H|. In the preceding packet's E-linear local constraint kernel,
include its requested output zeros and inactive outgoing coordinates.
Let G_minus be the image of its pre-diagonal state H16 D1(Q+CX),
restricted to H. Let G_plus be the projection onto the incoming state
coordinates H of the following packet's local constraint kernel. Both
are E-subspaces of E^H. They are independent of the current diagonal
b(theta+a). At the final emitted packet, use its output constraints only;
there is no additional final refresh or forced return.

Write g=dim G_minus and v=dim G_plus. Take a basis matrix B for
G_minus and a matrix P with kernel G_plus. Reveal the nonzero byte
diagonal D_b, but retain the independent a bytes. Put

    L := P D_b B,
    T(X) := P D_b (X I+D_a) B.

If rank L=g<16, a corresponding g-minor of T(X) has nonzero leading
coefficient and degree g. Since theta has degree sixteen over E, this
minor cannot vanish at theta. Hence T(theta) has full column rank.
No nonzero vector in G_minus extended to F can refresh into G_plus
extended to F. The branch is impossible for every a and every F-valued
message consistent with these local constraints. If g=16 and L has
full rank, then h=16 and G_plus=0; invertibility excludes the branch
directly. Full column rank requires g+v<=h. A rank-deficient L does not
establish feasibility; this is only a sufficient exclusion test.

Apply the exclusion before aggregating the two requested-zero sets.
Conditional on the fixed route, D1, masks, and zero sets, the b variables
are independent at distinct cuts. Conditional on b, each refreshed
coefficient still has 256 distinct a choices, preserving the existing
q=256 rank potential. Do not add rank L again as a source-rank debit.
Do not multiply separately averaged cut probabilities and packet kernels.
The next task is a positive two-packet kernel with this dependency order.
Start with disjoint packet pairs. Alternatively, a transfer state must
retain all shared data. Independently averaging overlapping pair kernels
and multiplying them would reuse shared D1 and zero-set data incorrectly.

A cheaper probability bound can instead use larger spaces independent
of the current D1. In the preceding packet retain only its output
constraints, and set W=im_E(Q+CX). In the next packet retain only output
constraints and incoming support H; its state space is

    V = {u in E^H : A[Z_next outside I_next,H] u=0},

embedded in E^16. For fixed W,V of dimensions g,v, the earlier joint
D1,b defect bound applies to P_V D_b H16 D1 B_W. It remains conditional
on the route and requested zero sets. It does not apply after making W
depend on current-D1 feedback restrictions. Fresh diagonal b alone has
no dimension-only mixing guarantee: it preserves coordinate subspaces,
including span(e_1).

For the conditional b-only calculation, actual subspace geometry gives
the projective union bound

    Pr_b[rank(P D_b B)<g]
      <= min(1, sum_{lines ell in P_E(G_minus)}
           |{v in G_plus : supp(v)=supp(ell)}|
                         /255^|supp(ell)|).

For one fixed line, a representative multiplied by D_b is uniform over
nonzero vectors on that support. Union over E-projective lines detects
rank failure of the E-matrix L. This is not a union over F-valued messages.
It is valid even when G_minus depends on the already fixed D1, but its
value must remain joint with the local scalar kernels and zero sets.

Two independent reviews checked this live-cut extension and its
conditioning restriction. No numerical gain or new certificate is yet
claimed. A product-of-linear-factors atom estimate is not a substitute:
H16 coupling produces sums of products, and arbitrary F-valued forcing
can exploit the allowed input images. The subspace exclusion handles
those inputs through a rank statement instead of a fixed-input atom count.

### Pool internal state masks before taking a fractional power

The output-only live-cut filter has a quantified limitation at the
current tilt. For two full incoming states and no allowed inputs, a
retained outside-zero set has at most three points in each RS block.
Its rank is its size. Full-column-rank exclusion requires at least
sixteen such requested zeros across the pair. Their total positive
scalar mass is at most q^-16, whereas the empty-zero term contributes
z^128. At z=.85, the relative removable mass is at most 2^-97.9884.
This comparison concerns the output-only kernel; additional feedback
filters can change its denominator. It nevertheless rejects this local
filter as the main way to improve that full-state output channel.

The present scalar operator splits every state-mask path before applying
the fractional power. This pays a separate concavity loss at every packet.
Instead, retain the intermediate-mask sum within a disjoint pair of packets.
This changes the proof operator, not the encoder or its randomness.

Fix the actual source set T, with |T|<=256, and use its existing exact-support
comparison factor eta and scalar outer density beta^|T|. After averaging
the a coefficients, the positive path bound is capped by one. For
0<alpha<=1, min(1,x)<=x^alpha. Jensen over the outer setup, conditional
on route and D1, bounds this power by the corresponding scalar-density
expression. The uniform eta and beta factors can be taken outside the
entire path sum. This leaves the fractional power of a positive sum of
products of local scalar kernels.

Partition the packet sequence into disjoint pairs. First sum over the
intermediate masks and both requested-zero sets in each pair. Split only
the pair-boundary mask paths by subadditivity of x^alpha. Conditional on
the route, distinct pairs use independent D1 variables. Their expectations
therefore factor; Jensen inside each pair bounds its fractional moment
by the fractional power of its raw expected kernel. The output-only last
packet is part of the final pair. No final refresh is introduced.

For fixed allowed-input sets, let K_i[h,j] be a uniform bound on the raw
expected one-packet kernel for one fixed incoming mask of size h and one
fixed destination mask of size j. Define D_j=binom(16,j). The raw pair
cap and its transfer entry are

    B[h,k] = sum_{j=0}^{16} D_j K_0[h,j] K_1[j,k],
    T[h,k] = D_k B[h,k]^alpha.

The intermediate multiplicity D_j is inside alpha; the boundary
multiplicity D_k remains outside. If an imported matrix M_i already
contains destination multiplicities, use

    T[h,k] = D_k ((M_0 M_1)[h,k]/D_k)^alpha.

Taking a second power of already powered one-packet entries is invalid.
The grouped calculation must start with raw alpha=1 kernels. At the
terminal, replace K_1 by its output-only column and omit D_k.

Quiet and nonquiet histories can still be separated positively. A pair
is quiet when its initial and intermediate masks are both empty. For an
empty initial mask, the nonquiet sum uses j>=1; for a nonempty initial
mask it uses every j. Apply the fractional power separately to these two
positive sums. This adds a category split but never subtracts close
quantities. The first implementation uses activity marker one.

The existing gamma choices can be retained. For each complete fixed
mask-path index and occupancy pattern, choose each packet's gamma from
its deterministic averaged caps before evaluating the random setup.
Each choice gives a valid positive expansion with the same global z,
eta, and beta. Summing these path-specific expansions is therefore valid.
It does not amount to choosing gamma retrospectively from the codeword.

For the paired route, condition on all half-packet occupancies. The
occupied point subsets in different halves are independent uniform
subsets. Jensen can average these subsets inside the raw pair kernel.
Within one full-pair half, both emission blocks share a point subset;
the frozen alpha=1 singleton formula retains that dependence. Masks
remain summation indices, not observations used to condition the route.

Two independent reviews checked the grouped-alpha order, including the
uniform source factor, raw-kernel normalization, and terminal. This
permits a new finite operator; it does not itself establish a distance
certificate or remove the need to cover every feasible route profile.

### Average route placements within a tile only after normalizing

For the pair-only marginal, p packets contain 32p routing slots. Let
c_d[h,k] be the degree-d coefficient after multiplying p raw packet
polynomials. Each packet polynomial includes its route-placement count
and destination-mask multiplicity. Thus c_d includes the endpoint factor
D_k as well as all interior-mask sums.

Conditional on a tile containing d full pairs, its binom(32p,d) occupied
slot subsets are equally likely. Jensen over those subsets gives the
powered tile coefficient

    binom(32p,d) D_k
      [c_d[h,k]/(binom(32p,d) D_k)]^alpha.

The normalization applies separately to the positive quiet and nonquiet
coefficients. For the terminal column, set D_k=1. Multiplying 128/p tile
polynomials and extracting degree f then counts all regional placements.
Divide that final coefficient by binom(4096,f), outside alpha.

This route Jensen step may lose more than explicitly averaging the
fractional moments of all p individual occupancy patterns. It permits
larger tiles without exponentially many occupancy axes. It is never
valid to raise an unnormalized route coefficient to alpha: that would
fractionalize the route multiplicity itself. Independent review checked
this normalization and the fixed-mask uniformity needed by the raw caps.

For orientation only, a simplified zero-coordinate penalty gives a
mask-sum factor (1+255^-alpha)^16 when split packet by packet, but
(1+255^-1)^(16alpha) at a pooled internal boundary. At alpha=.523,
their log2 values are 1.23871 and 0.04725. This calculation identifies a
potential source of loss; it is not the gain of the complete operator.

### Retain an outer-input weight floor before scalar density domination

The exact-T scalar comparison preserves a restriction that the current
local kernels discard. Fix a scalar source with exactly r active groups.
In each active group, some RS[64,32] constituent is nonzero, so at least
33 columns are nonzero. Each invertible four-coordinate column mixer
preserves column nonzeroness. Therefore the routed scalar outer word X
has at least 33r nonzero records, regardless of the route.

For v>=1, insert 1<=v^(wt(X)-33r) before replacing the scalar outer
distribution by its density bound. The marker is nonnegative, so the
same density domination applies. The resulting raw local kernel is

    K_v = E_{x,u}[v^wt(x) z^wt(x+A_H u) 1[inactive feedback is zero]],

where x and u have the same uniform E-valued auxiliary distributions
as the unmarked kernel. The global fractional-moment bound acquires
v^(-33r alpha). This factor is separate from the unchanged eta, beta,
and z factors. The restriction concerns exact-T scalar shadows; it
cannot be reused unchanged after allowing missing source groups.

A cheap first diagnostic multiplies every unmarked raw cap by v^ell,
where ell is the number of allowed input coordinates. This follows from
wt(x)<=ell before any support or rank relaxation. It remains valid for
each deterministic gamma choice. Quiet-to-quiet kernels permit an exact
replacement: incoming state and feedback are zero, so C_I x=0 and y=x.

For one emission block with ell allowed points, let a_w count kernel
vectors supported exactly on a fixed set of w points. The MDS rank of
these columns gives

    a_0=1,
    a_w=q^max(w-4,0) - sum_{j=0}^{w-1} binom(w,j) a_j,
    W_ell(t)=sum_{w=0}^ell binom(ell,w) a_w t^w.

Its exact marked quiet kernel is q^-ell W_ell(vz). Multiply these kernels
over the four emission blocks. For ell<=4 the kernel equals q^-ell and
is independent of v. At the output-only terminal, the quiet kernel is
[(1+(q-1)vz)/q]^ell for the total packet occupancy ell.

The coefficients a_w are nonnegative integers; this positive polynomial
also works when vz>1. Do not use a zero-subset expansion with negative
coefficients in that regime. Compare v=1 and v>1 using these same exact
quiet replacements in both runs, since replacing a loose unmarked quiet
cap can itself improve the bound. Two independent reviews checked the
weight floor, the density step, and these marked local kernels.

## Retain the number of active RS state blocks

The cardinality-only interface applies the weakest placement of h state
coordinates to every mask of size h. Refine the actual F-state class to
(h,b), where b is the number of nonzero four-coordinate RS blocks. There
are 35 classes, including (0,0). Their exact multiplicities are

    D(h,b) = C(4,b) [x^h] ((1+x)^4-1)^b.

Summing over b recovers C(16,h). At h=8, the counts for b=2,3,4 are
6,1968,10896. These classes concern actual F-state masks; they do not
assert that an auxiliary E-valued state has the same support.

Fix allowed input occupancies ell_i in the four emission blocks. For a
class (h,b), consider each placement t_i in {0,...,4} with sum t_i=h and
exactly b positive entries. Maximize the local cap over these placements;
do not multiply the maximum by C(4,t_i). Mask counts enter only through
the destination multiplicity D.

For each placement, the actual output has the deterministic floor

    f(t,ell) = sum_{i:t_i>0} max(0,13-ell_i).

Two rank lower bounds sharpen the existing scalar caps:

    k(t,ell) = sum_i max(0,min(4,ell_i)-t_i),
    a(t,ell) = sum_i max(0,t_i-max(0,ell_i-12)).

The first follows by deleting t_i rows from the rank-min(4,ell_i)
feedback block. The second follows by deleting 4-t_i columns from the
rank-min(4,16-ell_i) expansion block outside the allowed input positions.
Neither argument assumes that the normalized expansion basis is
superregular.

Write j for the outgoing support size and m=16-j. The feedback-only cap
for a fixed placement is

    j=0:       256^(-h-k),
    j>0,h>=m:  256^(-m),
    j>0,h<m:   256^(-h) [256^(-k)+(1-256^(-k))*255^(-(m-h))].

For the last case, retain a fixed k-dimensional image of the input
outside the incoming support. A uniform auxiliary input has zero image
with probability 256^-k. On a nonzero image, the MDS quotient equations
charge the independent nonzero D1 coefficients. The lower-rank choice
enlarges the resulting cap. This is a joint feedback calculation, not a
product of marginal rank probabilities.

A return to zero requires t_i=0 whenever ell_i=0. For j>0, a further
necessary condition is

    j + 4*#{i:ell_i>0} + sum_{i:ell_i=0} t_i >= 17.

It follows from the superregular H16 map and an upper bound on the
support of Q+CX. No restriction on the outgoing block count follows
from this cardinality condition alone.

For each placement, combine the feedback cap with the existing output,
retained-zero, padding, and joint-live caps. Use a(t,ell) in the padding
bound and f(t,ell) in the deterministic weight split. Maximize only after
taking these compatible minima. The implementation minimizes over three
gamma choices for each occupancy, incoming class, and outgoing cardinality.
Each choice is fixed before averaging the remaining randomness. Finally,
intersect with the old uniform-h cap, including its shared-point singleton
bound. The terminal uses only the output cap.

The extra classes need not increase matrix-power dimension. Let
R[(h,b),j] be the per-fixed-destination cap and define
L[j,(j,c)]=D(j,c), with all other entries zero. The raw transfer is R L,
so a p-step product is R(L R)^(p-1)L. Its interior products have dimension
17. After endpoint-normalized Jensen, D(j,c) remains outside the power.
The p-packet terminal product is R(L R)^(p-2)L T for p>=2. Thus it has
p-1 updates and no terminal feedback or destination multiplicity.

Independent proof and source reviews passed. The fixed evaluated-route
test improves the log2 expression from +37827.01773 to +37145.25256.
The zero/full16 channel is unchanged at +17650.47668, as required by
its unique block placements. This is a modest diagnostic gain, not a
new certificate. Further subdivision of intermediate masks alone cannot
remove that positive control at these fixed parameters.

## A fixed cone can pool live masks before averaging routes

The route-first calculation needs to sum partial state supports without
introducing a separate fractional-moment penalty at every mask boundary.
Use one positive envelope for all nonzero actual masks throughout the
trajectory. The envelope is fixed before sampling any route or setup.

Let q=256, D_h=C(16,h), and let K_HJ(I) denote a valid raw local cap for
actual incoming and outgoing masks H,J and allowed input set I. The old
cardinality matrix M_hj includes D_j. Choose positive p_h and assign
weight f_h=p_h/D_h to each actual mask of size h. Define

    Q=M_00,  A=max_j M_0j/p_j,
    L=max_j (sum_h p_h M_hj)/p_j,  E=sum_h p_h M_h0,
    (T,U)=(t_0, sum_h p_h t_h),                 h,j in {1,...,16}.

An incoming row satisfying v_H<=s f_|H| maps to a live row bounded by
s L f_|J|. A birth enters this cone with coefficient A; a return and the
output-only terminal cost E and U. Thus the raw sum over all fine masks
of a fixed binary zero/live schedule is bounded by its two-state product.
The test uses

    f_h=min(1, [q/(q-1)] (q-1)^(-(16-h))).

The same vector is used for every packet and occupancy. Its normalization
does not affect a complete zero-initialized trajectory. Independent local
minimization over cone vectors would not have this compositional property.

Actual masks are exhaustive union indices, not independently sampled
states. First fix a binary schedule and sum its fine actual-mask paths.
Apply the cone bound to that positive raw sum, then take alpha. Only then
average the route occupancies and extract the exact regional coefficient.
This ordering retains quiet-route suppression. It introduces neither an
extra binary-state multiplicity nor another fractional power per region.

### Exact supports improve the return to zero

At a return, the pre-refresh equation is u=C_I x. The former cardinality
maximum may count the same auxiliary x under many incoming masks. Enforce
the exact incoming support before summing those masks.

Fix H, its size h, the allowed input set I, and a requested output-zero
set Z. If a required incoming-coordinate form vanishes identically on
the local E-linear constraint kernel, the actual F-valued branch is
impossible. Otherwise the hyperplane-union comparison preserves all h
nonzero forms at cost q/(q-h). Apply this comparison before the positive
sum over Z. Eliminating u=C_I x supplies a new per-mask cap

    Knew_H,0(I) := [q/(q-h)] q^-h
                  E_x[z^wt(S_I x) 1(supp(C_I x)=H)],
    S=Id_64+AC,  x uniform in E^ell,  ell=|I|.

The normalization is q^-h, not (q-1)^-h: the nonzero conditions are
indicators under uniform auxiliary sampling, not conditioned sampling.
The indicators for distinct H are disjoint. Hence the actual-mask cone
return sum is at most

    E_new(I) = kappa q^-16 E_x[z^wt(S_I x) 1(C_I x!=0)],
    kappa=max_(1<=h<=16) f_h q^(16-h) q/(q-h)
         = q^3/[(q-1)^2(q-15)] = 16777216/15671025.

The maximum occurs at h=15. There is no extra source-density, gamma, or
nonzero-coordinate factor. Since CA=0 in characteristic two, S^2=I.
Its restriction to ell input coordinates is injective. Selecting ell
independent output coordinates therefore bounds its moment by g(z)^ell,
where g(z)=z+(1-z)/q. On C_I x=0, S_I x=x. Subtracting this exact quiet
contribution gives

    E_new(I) <= kappa q^-16 [g(z)^ell-Q_I(z)].

Here Q_I is exact, not a quiet upper bound. For four block occupancies
ell_b it equals the product of

    Q_l(z)=sum_(k=0)^l C(l,k) z^k [(1-z)/q]^(l-k) q^-min(4,k).

The MDS rank of the remaining k columns of C proves this positive
requested-zero formula. A directed implementation must use a lower
enclosure for Q_I when subtracting it. The empty-input return is zero.

The new return expression need not bound the old cardinality-max
expression. Both independently bound the actual-mask cone return.
One may select the smaller proved cap deterministically from the allowed
occupancy, before the auxiliary averaging. This selects a proof method;
it does not multiply incompatible relaxed scalar measures.

### Retain the RS output floor at returns

For one block, put c=min(4,l) and d=max(0,13-l). For a fixed nonzero
syndrome C_I x, the inside-I output ranges over an affine coset of
ker C_I. Its q^-l-normalized weight sum is at most that of the zero coset: in the
positive requested-zero expansion, each affine system has either zero
solutions or the same number as its homogeneous system. The outside-I
output has at least d nonzeros, because it belongs to an RS[16,4] word.
Summing the q^c-1 nonzero syndromes gives the compatible bound

    D_l := min(g(z)^l-Q_l(z), (q^c-1) Q_l(z) z^d)

on the block's nonzero-syndrome moment. Combine four blocks positively,
starting with Qprod=1,R=0 and updating

    R <- R (Q_l+D_l)+Qprod D_l,   Qprod <- Qprod Q_l.

Then kappa q^-16 R bounds the packet return. This recurrence retains the
event that at least one block has nonzero feedback. It does not subtract
products of loose upper bounds. Exact GF4 affine-coset checks and two
independent mathematical reviews passed.

### Use the implemented return matrix for one or two allowed inputs

The preceding RS floor need not be tight for the actual return matrix.
Within one16-coordinate block, write S_b=I_16+A_b A_b^T. For an allowed
point set I, define its nonzero-syndrome moment under uniform x in E^I by

    D_I(z)=q^-|I| sum_(x:C_I x!=0) z^wt(S_b[:,I]x).

For |I|<=2, the columns of C_I are independent. The excluded event is
therefore exactly x=0, and the quiet moment is q^-|I|. The actual S_b
is invertible because A_b^T A_b=0 and S_b^2=I_16.

The normalized basis gives

    (A_b A_b^T)[a,t]=(1+at)(1+B(a)B(t)),
    B(a)=a(a+1)/6,                  a,t in{0,...,15} subset AES GF256.

The only zero entries of S_b are(0,0),(2,3),(3,2). Thus three columns
have weight15 and the remaining thirteen have weight16. For |I|=1,

    D_I(z) <= (q-1)z^15/q.

This is a pointwise bound for every allowed point. It improves the
generic floor12 without changing the encoder or its randomness.

For |I|=2, every row of S_b[:,I] is nonzero. Each row annihilates
exactly one projective direction of E^2. Let n_l count rows with kernel
line l, and let m count distinct such lines. Every projective direction
contains q-1 nonzero vectors of equal output weight. Consequently,

    D_I(z)=(q-1)/q^2 * [(q+1-m)z^16 + sum_(l:n_l>0) z^(16-n_l)].

Only sixteen row directions per point pair are needed. Maximizing this
positive polynomial over the120 possible pairs gives a uniform cap.
Intersect that cap with the previous D_2, and likewise intersect the
one-input cap with D_1. Keep every other D_l and every exact Q_l unchanged.
The same positive four-block recurrence then bounds a return, with the
same factor kappa q^-16. Take alpha only after constructing that cap.

This argument does not average point pairs independently. It is valid
when paired emission blocks share their allowed point set. Independent
source-level algebra review checked the column weights and projective
enumerator; the numerical impact requires a separate controlled replay.

### Remember the first live packet after a birth

A birth whose input occupies b of the four emission blocks has incoming
feedback support at most 4b. Superregularity of H16 makes its nonzero
outgoing support at least 17-4b. Keep that restriction for one packet
instead of replacing the newborn row by the whole live cone immediately.

Use states zero, newborn1,...,newborn4, and mature. A birth enters only
newborn_b, where b is fixed by the current occupancy. Its envelope is
f_h 1(h>=17-4b). On the next packet, a live transition enters mature
with the corresponding restricted sum in L. A return uses the restricted
sum in E, capped by the independently valid unrestricted return bound.
The terminal uses the restricted output-only U.

If the next packet occupies b' blocks, return requires h<=4b'. Therefore
b+b'<=4 excludes an immediate return. In particular, a birth in two
blocks activates at least nine state coordinates, while another
two-block packet can cancel at most eight. The birth label is determined
by the route and adds no multiplicity. The numerical six-state operator
retains this information across regional boundaries as well.

These refinements change the proof, not the implemented coefficient
distribution or encoder. The excursion ledger reports their numerical
scope. None currently supplies the missing whole-code support cover.

### A smaller independently valid cone

The directed check need not reproduce every historical local refinement.
Keep the same fixed mask envelope f_h and put p_h=C(16,h)f_h. For a
packet with block occupancies l_1,...,l_4, let ell=sum l_i and
g=z+(1-z)/q. The following simpler entries suffice:

    Q=product_i Q_(l_i),
    A=g^ell-((1-z)/q)^ell,              A=0 when ell=0,
    E=kappa q^-16 R,
    L=U=sum_(h=1)^16 p_h T_h,
    terminal=(g^ell,U).

Here R is the positive nonzero-syndrome recurrence above, and T_h is
the original output-only cap for one fixed incoming mask of size h.
These entries bound the same actual-mask sums, independently of the
older cardinality-return and feedback refinements.

For activation, retain only that at least one input coordinate remains
after expanding the requested output-zero set. All other terms are
impossible for an actual birth. Their positive scalar mass is
g^ell-((1-z)/q)^ell. For each fixed nonempty outgoing mask J, D1 mixing
contributes at most f_|J|. Dividing by this fixed envelope gives A.
This deletes an impossible requested-zero branch; it does not assume
that auxiliary scalar support equals actual extension-field support.

For completeness, T_h can be computed without the feedback table.
For each gamma in {0,1/2,1}, put a=z^gamma, b=z^(1-gamma), and
g_b=b+(1-b)/q. In one block with l allowed inputs and t incoming
state coordinates, use g_b^l if t=0. Otherwise set n=16-l and

    w_v=C(n,v)b^(n-v)(1-b)^v,           0<=v<=min(3,n),
    W=sum_v w_v,
    J=sum_v w_v q^-max(0,t-4+v),
    B=q^-t W+(1-q^-t)b^max(0,n-3).

The block cap is g_b^l a^max(0,13-l) min(J,B). For t=1, also
include b^max(0,n-2) in this minimum. Multiply the four block caps
and maximize over t_i in {0,...,4} with sum t_i=h. Then minimize
over the three gamma choices to obtain T_h. Each choice is fixed
from the occupancy and h before auxiliary averaging. No state-mask
multiplicity belongs inside this maximum.

For every fixed outgoing mask J, the original joint-live bound is
K_HJ<=f_|J| T_|H|. Thus the same T_h bounds both the live transition
and the terminal sum, giving L=U. This argument requires the original
output-only caps. A later, independently tightened terminal cap cannot
be inserted into the live transition without another joint proof.
All these simplified formulas are uniform over point subsets at the
given occupancies; no shared-point averaging refinement is used.

The fractional-moment order remains essential. Fix the actual source
set T and complete route first. For each binary zero/live schedule,
sum all refining fine-mask paths in its positive conditional event bound.
Apply min(1,x)<=x^alpha, then Jensen over outer setup and D1 for that
fixed schedule. The established scalar interface contributes the common
factor eta beta^|T| z^-L. The fixed cone bounds the remaining raw sum
by a product of its schedule entries. Only now take the packetwise
alpha powers and average the route. Finally sum the binary schedules.
This introduces no fractional power on a route-placement count and no
new alpha at a regional boundary. The resulting common prefactor is
(eta beta^|T| z^-L)^alpha, applied once to the complete chain.

### Count source runs without identifying rotations

The two alternating matchings retain a source support's two full-pair
counts. They do not retain the locations or spacings of its runs:
independent uniform pair permutations and fair flips erase those details.
For size256 support, write f_0=128-a and f_1=128-b. Its number of cyclic
runs is c=a+b, and its singleton counts are 2a and 2b. A run is a maximal
nonempty interval of active positions on the indexed8192-cycle.

Here is an exact count for each ordered profile. More generally let
the cycle have 2M positions, weight 2w, and v=M-w, with 0<w<M. Mark
one of the c boundaries at the start of a run. The parities of the 2c
alternating run-start and run-end boundaries determine which lengths
are odd. Exactly 2a of these boundaries lie in the first matching.
Let H_c(d,o,t) count these parity words with d first-matching boundaries,
o odd run lengths, and t odd gap lengths. It is the positive coefficient

    H_c(d,o,t)=[u^d x^o y^t] tr((D E(x) D E(y))^c),
    D=diag(u,1), E(x)=[[1,x],[x,1]].

The one-step matrix is [[u^2+uxy,u^2y+ux],[ux+y,uxy+1]]. For o=2i
and t=2j, the positive length compositions number
C(w+i-1,c-1) C(v+j-1,c-1). Odd o or t cannot give the even totals.
There are M starting positions of each parity; each indexed support
has c possible marked starts. Consequently its ordered profile count is

    N(a,b)=(M/c) sum_(i,j) H_c(2a,2i,2j)
                              C(w+i-1,c-1) C(v+j-1,c-1).

Division by c removes the marked run, not rotations. The total over
a+b=c is (2M/c) C(2w-1,c-1) C(2M-2w-1,c-1). The implementation uses
exact integers and checks final divisibility. Exhaustive small-cycle
enumeration checks both matching counts and this total.

For one ordered profile, extract regular matrices R_a,R_b and terminal
column V_b from the exact route coefficients. Its finite chain is
e_0 (R_a R_b)^127 R_a V_b. Multiply by N(a,b) and the common powered
source factor, then sum all selected ordered profiles. No additional
orientation factor is present. Singleton degree at most 2c_max retains
every required coefficient for all profiles with at most c_max runs.

### Count fragmented supports using a four-vertex graph

The boundary-parity recurrence becomes expensive when the support has
many runs. A second exact count uses only O(r) positive terms per
profile. It applies to odd as well as even source weights.

Fix an indexed cycle of length2M, weight r, matching counts f_0,f_1,
and t active even positions. Use vertices E0,E1,O0,O1, recording position
parity and whether that position is active. Consecutive cycle positions
give the following directed edge multiplicities:

    E->O: [[M-r+f_0, r-t-f_0], [t-f_0, f_0]] = [[A,B],[C,D]],
    O->E: [[M-r+f_1, t-f_1], [r-t-f_1, f_1]] = [[E,F],[G,H]].

The outdegrees are d=(M-t,t,M-r+t,r-t), equal to the indegrees.
Sum over feasible t, so all edge counts are nonnegative. Delete vertices
of degree zero. Let tau be the directed spanning-tree cofactor of the
remaining graph. Balanced indegrees and outdegrees make its root
cofactors equal. If all four vertices remain, direct expansion gives
tau=C E(G+H)+D G(E+F). A disconnected graph has tau=0.

To count walks from a fixed root v, first label parallel edges. Choose
a last-exit tree toward v. Order outgoing edges freely at the root;
at every other vertex put its tree edge last. Following these orders
gives exactly tau d_v product_u(d_u-1)! closed walks using every edge.
Each unlabelled position sequence has product_edges(n_edge!) labellings.
Sum over the allowed even roots; their degrees sum to M. Rewriting the
factorial quotient therefore gives the nonnegative integer

    N_t = M tau product_u C(d_u, n_(u,first outgoing class))
                      / product_u d_u.

Omit zero-degree vertices from both products. The total profile count
is sum_t N_t. Each division is exact. These walks start at the fixed
position0. We do not quotient indexed words by cyclic rotation;
periodic words require no correction. Four initial binomials and their
exact adjacent-t recurrences evaluate the sum without large factorials
or a parity-histogram cache.
Exhaustive all-weight small-cycle tests, all-root cofactor checks, and
comparison with the first counting method verify the implementation.

### Pack coefficient products without changing the bound

For target full degree at most F and singleton degree at most S, map
the bivariate monomial x^f y^s to t^(f+(2F+1)s). Both factors have
f<=F, so a product has full degree at most2F, strictly below the radix.
No forbidden full degree can carry into the singleton degree. Decode
after each product and discard f>F or s>S. At source weight256, also
discard2f+s>256. These are downward-closed truncations: nonnegative
future degrees cannot bring a discarded term back to a target.

Different rectangles can recover different marginal matching profiles.
Exact positive variable scales improve interval conditioning and cancel
in each coefficient's normalization before composing the regional
matrices. They do not change the local alpha powers or introduce a
Chernoff approximation. Small exact matrix products check the packing,
including terms whose full degrees must be discarded rather than carried.

### Use the refined cone after padding a missing scalar group

For257<=r<=511, the earlier missing-group argument can use h=1. Fix an
actual support T of size r and its initial E-kernel U. Discard the term
if U=0 or any group projection vanishes identically on U. A uniform nonzero
scalar kernel vector misses fewer than r/q groups on average. Markov's
inequality shows that at least c_1=1-r/(2q)>0 of these vectors miss at
most one group. Therefore

    P_q(dim U) <= [n_T+sum_(g in T) n_(T without g)] / [(q-1)c_1].

Fix the actual state masks, D1 values, and requested output zeros before
padding. Embed each reduced scalar input space into the actual-T input
space by setting the missing group's256 routed coordinates to zero.
The normalized positive kernel increases by at most q^256 when those
coordinates are removed. This comparison is pointwise; it does not
assume independent routes for the actual and scalar supports.

Perform this padding before applying the sharper return hyperplane
comparison or any other local actual-support refinement. All feasibility
filters still use I(T). A required state form may vanish on the reduced
shadow space while remaining possible on the actual-T space, so applying
the refined cone directly to I(T without g) would be invalid.

The resulting complete source factor, before its fractional power, is

    Omega_r = beta^r [1+r*q^256/beta] / [(q-1)(1-r/(2q))].

For r<=256 use Omega_r=eta_r*beta^r instead. In either case the global
factor is (Omega_r*z^-L)^alpha, applied once. The missing-group sum stays
inside this power. No new encoder randomness or full-field distribution
is used. This argument licenses a calculation; it does not assert that
the resulting numerical bound passes.

### Retain the missing groups' nonzero projections until their outer average

The q^256 padding loss above can be avoided. This requires a different
order of averaging, not a smaller assumed rank. Fix the actual support T,
route, D1 matrices, state masks, and requested output zeros. Eliminate
the auxiliary state columns M and let D=L_M N be the resulting E-valued
constraint map on ambient outer outputs. Write R=sum_t rank M_t.

Let J be a missing set, S=T without J, and 1<=j=|J|<=q. Condition on
the outer maps G_S. Let n_S count scalar source vectors nonzero on every
S group and satisfying D_S G_S x=0. Quotient the constraint codomain by
im(D_S G_S), and let L_J be the image of D_J in that quotient. Then

    projection_J ker(D_S G_S,D_J G_J) = ker(L_J G_J).

This identity is over E and survives extension to F. A retained actual-T
witness implies that every J-group projection of this E-kernel is nonzero.
A vector space over E cannot be covered by at most q proper subspaces:
enlarge them to hyperplanes and count their nonzero vectors, whose union
has size at most q*(q^(d-1)-1)<q^d-1. The one-dimensional case is immediate.
Thus the kernel contains a scalar vector nonzero on every J group, and
at least q-1 such vectors by scalar multiplication.

The outer counting measure applies to an arbitrary fixed linear map L:
the expected number of nonzero source vectors with L G x=0 is at most
beta*q^(-rank L). This follows by summing the point-density domination
over ker L; it does not assume an isotropic outer. Independence across
the J groups gives the joint beta^j bound. Explicitly, let n_J count
vectors x_J nonzero in every J group with L_J G_J x_J=0. Then
1_retained<=n_J/(q-1) and E_(G_J)n_J<=beta^j*q^(-rank L_J).
Since n_S is fixed under the G_J average,

    E_(G_J)[n_S * 1_retained]
      <= n_S * beta^j * q^(-rank L_J)/(q-1).

The coded image im(D_S G_S) is contained in the ambient image im(D_S).
Consequently rank L_J >= Delta_J := rank D_T-rank D_S. Define the
normalized ambient kernels K_S=q^-R*q^(-rank D_S) and likewise K_T.
They satisfy K_S=q^Delta_J*K_T exactly. Average G_S only now, using its
counting density, to obtain

    E_outer[1_retained * q^-R * n_S]
      <= beta^r * q^-Delta_J * K_S/(q-1)
       = beta^r * K_T/(q-1).

This is a simultaneous average over the missing maps, not a product of
separately conditioned group events. Retain the missing-projection event
until that average. Actual-T local feasibility filters remain fixed;
additional outer-dependent global filters may be dropped. Apply the
refined actual-T cone and fractional power only after this comparison.
Neither the local kernel nor the implemented randomness changes.

For h<r, h<=q, and c_h=1-r/[q(h+1)]>0, the earlier missing-count argument
now gives the complete source factor

    Omega_r = beta^r / [(q-1)c_h]
                * [1 + (1/(q-1))*sum_(j=1)^h C(r,j)].

Choosing h=floor(r/q) is valid through r=8192, where h<=32. For one
missing group, the simpler conservative factor beta^r*(r+1)/[(q-1)c_1]
is also valid. This joint derivation passed an independent proof audit.
It replaces the earlier ambient padding penalty; it does not by itself
supply numerical certificates for uncomputed source sizes.

The choice of h need not be the smallest admissible integer. For each
fixed actual support size r, let H_r contain the integers with
0<=h<r, h<=q, and c_h>0. Every h in H_r bounds the same actual-T
ambient kernel, with all feasibility filters still referring to T.
Consequently one may replace eta_r by

    min_(h in H_r) [1 + sum_(j=1)^h C(r,j)/(q-1)] / [(q-1)c_h].

This deterministic choice depends only on r. There is no union over h,
and monotonicity permits taking the same minimum before the positive
fractional power. For r<q the admissible h=0 expression is
q/[(q-1)(q-r)]; it can be below one because it normalizes a projective
count, not an event probability. It is unavailable at r=q. The current
numerical helpers retain h=floor(r/q); this optional refinement has not
been substituted into their frozen receipts.

### Bound a cell of source sizes and run counts

The regional coefficient tables can support a contiguous cover without
assuming monotonicity in source size. Fix z,alpha and a cell C of integer
pairs(r,c), where r is actual source weight and c is its cyclic run count.
Use 1<=r<8192 and 1<=c<=min(r,8192-r); the full-cycle support is separate.
For every size in the cell, let Omega_r be a proved complete source factor.
All required marginal coefficients must actually be computed.

Let R_(f,s) be its2-by-2 regular regional matrix and V_(f,s) its terminal
column, with r=2f+s. Fix positive rational marks rho,sigma. Scale each
marginal by

    d_(r,f)=Omega_r^(alpha/256)*rho^(r/256)*sigma^(f/128),
    Rhat_(f,s)=d_(r,f)*R_(f,s), Vhat_(f,s)=d_(r,f)*V_(f,s).

For every compatible ordered pair of matching profiles in C, form
B=Rhat_0*Rhat_1 and Cend=Rhat_0*Vhat_1. Compatibility requires the same
r in both profiles and f_0+f_1=r-c. Let Bbar and Cbar be entrywise upper
envelopes of these pair products over the cell. The finite chain has128
regions of each orientation, so the scales contribute exactly
Omega_r^alpha*rho^r*sigma^(r-c). Positivity then bounds the cell's union by

    z^(-alpha L) * e_0*Bbar^127*Cbar
      * sum_((r,c) in C) N_(r,c)*rho^(-r)*sigma^(-(r-c)),
    N_(r,c)=(8192/c)*C(r-1,c-1)*C(8191-r,c-1), e_0=(1,0).

Here N_(r,c) counts indexed source supports, independently of subsequent
routing. The envelopes enlarge products already valid at fixed r; they
do not allow actual support size to change along a path. The terminal
column preserves the output-only final packet. Sum disjoint cells to
cover their union. Directed maxima and powers are required in a verifier.
The envelope may be loose, but it requires neither interpolation nor
monotonicity of the marginal matrices. The directed cell verifier now
uses this bound for contiguous source-size bands; the goal ledger records
the exact bands and witnesses.

### A cheap large-support diagnostic with positive coefficient bounds

Exact regional coefficient tables become expensive at larger source sizes.
For a fixed matching profile (f,s), let P(u,v) be the positive packet
matrix, after raising each local cone cap to alpha and before averaging
route placements. Let T(u,v) be its output-only terminal column. The
markers u,v are positive rationals. A full pair has weight u; a singleton
has total weight v, divided equally between its two orientations.
Each packet has32 pair slots, and each region has4096 pair slots.

Define M(f,s)=C(4096,f) C(4096-f,s). Positivity of the coefficients gives
the entrywise bounds

    R_(f,s) <= P(u,v)^128 / [M(f,s) u^f v^s],
    V_(f,s) <= P(u,v)^127 T(u,v) / [M(f,s) u^f v^s].

The fractional power precedes route averaging, as in the exact coefficient
calculation. The coefficient inequality adds slack but changes neither
the local kernel nor the setup distribution. If the two matching profiles
coincide, the complete chain is bounded by

    H(u,v) / [M(f,s) u^f v^s]^256,
    H(u,v)=e_0 P(u,v)^32767 T(u,v).

Multiply by the proved source factor (Omega_r z^-L)^alpha and an upper
bound on the number of indexed supports of this profile, where r=2f+s.
For example, the exact number of all c-run supports with c=s is a valid
upper bound. This counts only the selected profile, not every profile
at source size r. The terminal contains route degrees but no refresh.

Marker calibration can be diagnosed without changing any other bound.
With z and alpha fixed, H is a positive polynomial in u,v. Under its
tilted monomial distribution, let F_total and S_total be the total full
and singleton degrees. Then

    d_ln(u) log2(bound) = [E(F_total)-256f]/ln(2),
    d_ln(v) log2(bound) = [E(S_total)-256s]/ln(2).

Positive dual matrix arithmetic computes these means for the finite chain.
No additional alpha multiplies these derivatives: alpha already appears
inside P and T. The log bound is convex in (ln u,ln v), since log H is a
log-sum-exp and the logarithmic denominator is affine. A gradient sign
justifies a local search direction, not an unverified finite-step gain.
Any final certificate still requires a directed evaluation at fixed markers.

### Reuse route markers across profile cells

To cover all profiles at one size r without a separate marker solve for
each, partition the feasible full counts into integer intervals I_i. The
feasible range is max(0,r-4096)<=f<=floor(r/2). Fix one positive marker
pair (u_i,v_i) per interval, with common z and alpha. Define

    B_i=P(u_i,v_i)^128, V_i=P(u_i,v_i)^127 T(u_i,v_i),
    H_ij=e_0(B_i B_j)^127 B_i V_j,
    d_i(f)=C(4096,f) C(4096-f,r-2f) u_i^f v_i^(r-2f).

For f_0 in I_i and f_1 in I_j, the complete conditional chain is at most
H_ij/[d_i(f_0)d_j(f_1)]^128. Each matching has128 regions, including the
terminal region. Keep the cell pair ordered: terminal asymmetry need not
give H_ij=H_ji.

Fix a feasible cyclic run count c, so f_1=r-c-f_0. The compatible values
of f_0 form an integer interval [a,b], obtained by intersecting the two
cells and the feasible full-count range. Empty intersections contribute
nothing. To bound its denominator, set M_r(f)=C(4096,f)C(4096-f,r-2f).
Its successive ratio is

    M_r(f+1)/M_r(f)
      = (r-2f)(r-2f-1) / [(f+1)(4096-r+f+1)].

This ratio decreases throughout the feasible range. Multiplication by
u_i^f v_i^(r-2f) preserves log concavity. Hence
log d_i(f)+log d_j(r-c-f) is concave, and its minimum on [a,b] occurs
at an endpoint. Let D_ijc be the smaller endpoint product. Then

    (Omega_r z^-L)^alpha
      * sum_c N_(r,c) max_(compatible ordered i,j) H_ij / D_ijc^128

bounds every source support of size r. The run count is charged once,
not once per cell. This uses only positive coefficient inequalities and
scalar envelopes; it adds no Jensen step or entrywise matrix maxima.
Zero-multiplicity profiles can safely enlarge an envelope.

If a rare cell dominates, its source multiplicity can be bounded more
sharply. The exact number of indexed supports with first-matching count f
is K_r(f)=M_r(f)2^(r-2f). For a compatible interval [a,b], its count is
at most the minimum of N_(r,c), sum_(f=a)^b K_r(f), and
sum_(f=r-c-b)^(r-c-a) K_r(f). These bounds define capacities for the
ordered cells. Actual cell counts sum to N_(r,c). Allocating that total
greedily to the largest cell-bound upper endpoints, subject to capacities,
therefore gives a valid relaxed upper bound. Check that total capacity
is at least N_(r,c) and that the allocation leaves no remainder.

Both formulas passed independent proof review. Their computational value
still requires a complete directed profile-cover calculation; the
selected-profile diagnostics do not supply that calculation.

### Retain regional degree residues in the coefficient bound

Ordinary Cauchy evaluation includes every regional degree. Some of that
slack can be removed with a small positive polynomial calculation. Fix
positive markers u,v and an integer m>=1. In each local packet entry,
assign the term with full/singleton counts (f,s) the additional monomial
t^(2f+s), and calculate in the quotient by t^m-1. Retain the same local
cap raised to alpha before route averaging. Call the resulting matrix
P_m(t) and the output-only terminal column T_m(t).

Let [t^a]_m denote the sum of coefficients whose ordinary degrees equal
a modulo m. For the target regional profile (f,s), put r=2f+s. Then

    R_(f,s) <= [t^(r mod m)]_m P_m(t)^128 / [M(f,s)u^f v^s],
    V_(f,s) <= [t^(r mod m)]_m P_m(t)^127 T_m(t) / [M(f,s)u^f v^s].

To prove these inequalities, expand the finite products before reducing
exponents. Their coefficients are nonnegative. Every term of the desired
bivariate coefficient has total degree2f+s=r, so the residue sum includes
it with weight u^f v^s. It may include additional terms, which only enlarge
the entry. Each residue sum is also at most the sum over all residues.
Thus these bounds are entrywise no larger than ordinary Cauchy bounds at
the same markers. No conditional-independence or anti-concentration claim
is used, and the source factor is unchanged.

Extract the residue separately for each completed region, before composing
the ordered regional chain. The terminal contributes its route degree but
no final refresh. With m=1 the calculation exactly recovers ordinary
Cauchy evaluation. Summing every residue must likewise recover P(u,v)^128
or P(u,v)^127 T(u,v); these are useful implementation controls. This
refinement passed independent proof review. At the calibrated size512
equal profile(f,s)=(223,66), moduli4,8,16 recover512,768,1024 bits
respectively. The modulus16 bound remains positive. These are directed
selected-profile diagnostics, not a certificate for the entire source size.

### Retain both regional occupancy residues

The same argument can retain full-pair and singleton counts separately.
Let P(X,Y) be the nonnegative packet matrix before marker evaluation,
and let T(X,Y) be its output-only terminal column. The coefficient of
X^i Y^j includes its exact route multiplicity and its local cap raised
to alpha. For positive markers u,v, compute P(ux,vy) and T(ux,vy) in
the positive residue representation x^a=y^b=1. Then the residue
(f mod a,s mod b) of the regional product contains the desired coefficient
times u^f v^s, together with only nonnegative aliases. Thus

    R_(f,s) <= [x^(f mod a)y^(s mod b)] P(ux,vy)^128
                   / [M(f,s)u^f v^s],
    V_(f,s) <= [x^(f mod a)y^(s mod b)] P(ux,vy)^127 T(ux,vy)
                   / [M(f,s)u^f v^s].

Extract each completed region's prescribed residue before composing the
ordered regional chain. Summing all residues recovers ordinary Cauchy
evaluation. If both new moduli are multiples of the respective old
moduli, the target entry cannot increase. With equal moduli m, summing
the residue pairs satisfying2i+j=r modulo m recovers the one-axis
total-degree filter; the selected pair is no larger than that sum.
These statements follow from nonnegative coefficient inclusion, not
from a uniform-residue model.

An alternative is exact extraction: after each positive polynomial
product retain only0<=i<=f and0<=j<=s. This downset truncation preserves
the exact target coefficient because exponents never decrease. It adds
no alias terms. Both refinements leave source factors, profile counts,
alpha placement, and the missing final state update unchanged. The
lemma passed independent proof review; numerical efficacy is separate.

### Reuse one residue table across marginal profiles

Exact extraction closes the selected clustered512-group profile, but a
complete source-size bound must include every ordered profile. A residue
table can serve several marginals without assuming interpolation. Fix
z,alpha and the local cone caps throughout this calculation. For a positive
marker pair(u_i,v_i), cache the two cyclic regional products

    B_i=P(u_i x,v_i y)^128,
    C_i=P(u_i x,v_i y)^127 T(u_i x,v_i y).

For any f in0,...,256, put s=512-2f and
d_i(f)=C(4096,f)C(4096-f,s)u_i^f v_i^s. Extract the residue pair
(f mod a,s mod b) from B_i and C_i, and divide by d_i(f). Each result
bounds the same exact regional matrix R_f or terminal column V_f,
regardless of where the marker pair was calibrated.

One may assign a marker pair deterministically to each f. Alternatively,
take entrywise minima across these independently valid bounds. Denote
the selected bounds by Rbar_f and Vbar_f. Positive matrix multiplication
preserves their upper-bound property. For each ordered actual profile
pair(f_0,f_1), use

    N(f_0,f_1) (Omega_512 z^-L)^alpha
       e_0 (Rbar_f0 Rbar_f1)^127 Rbar_f0 Vbar_f1,

then sum all ordered pairs. The terminal orientation is not interchangeable.
There is no additional union over marker choices, no new fractional power,
and no profile interpolation. Fixed marker cells are a computational choice,
not a claim that nearby profiles have nearby probabilities. This reuse
argument passed independent proof review. The numerical complete-size
calculation remains necessary.

The same packet polynomial can also serve neighboring source sizes when
z,alpha and the local caps stay fixed. For each feasible(r,f), extract
(f,r-2f) and divide by its own normalization and marker powers. The source
factor remains Omega_r. For r>4096, the feasible lower limit is f=r-4096;
it is not zero. A retained exact rectangle must contain every requested
target, while a positive residue table has no target-degree cutoff.

The earlier source-size/run-count cell lemma can combine such bounds.
Its envelopes must range over compatible pairs with the same r and retain
the output-only terminal. Before expanding a band, test the envelope loss
at one certified size: replacing exact profile counts by run-total maxima
can discard the margin. Neither this reuse argument nor successful endpoint
calculations certify an uncomputed interval of source sizes.

### Retain exact counts when grouping ordered profiles

A run-total maximum can charge a rare, unfavorable profile for every
support with that run count. To avoid this loss, partition actual ordered
triples p=(r,f_0,f_1) into disjoint finite cells. Let N_p be the exact
indexed support count. Each triple uses two matching profiles of the same
source size r; its terminal belongs to the second matching.

Fix common z,alpha and positive rational marks rho,sigma within a cell.
Use the scale d_(r,f) from the source-size cell lemma and define

    B_p=d_(r,f_0)d_(r,f_1) R_(r,f_0) R_(r,f_1),
    C_p=d_(r,f_0)d_(r,f_1) R_(r,f_0) V_(r,f_1).

Let Bbar and Cbar be entrywise upper bounds for these objects over the
cell. Positivity bounds its failure contribution by

    z^(-alpha L) e_0 Bbar^127 Cbar
      * sum_(p in cell) N_p rho^(-r) sigma^(-f_0-f_1).

Indeed, the 128 pair factors contribute Omega_r^alpha rho^r sigma^(f_0+f_1).
The scalar weights cancel rho^r sigma^(f_0+f_1), leaving Omega_r^alpha.
The remaining expression is the existing bound for each actual profile.
Taking the entrywise envelopes only increases each nonnegative chain.
The scalar sum uses exact counts, not the total count for its run class.

A singleton cell recovers the original profile bound. Reversing a profile
changes its terminal, so reversed triples remain separate even though
their counts agree. Directed arithmetic is required for numerical claims.
This lemma permits aggregation; it does not assert that a chosen cell
width preserves enough margin. Test that loss at size512 before using
cells to extend the certificate.

### Count high-contribution profiles first and bound the tail

Exact profile counts need not all be evaluated when a certified tail is
already negligible. Fix a source size r and run count c. List every
feasible ordered profile p in this class. Let H_p be its nonnegative
conditional-chain majorant and u_p a certified upper endpoint. Sort the
profiles by decreasing u_p, keeping the two terminal orientations separate.
Let N_c be the exact total count of indexed supports with r active groups
and c runs.

After evaluating the exact counts N_1,...,N_m, set

    R_m = N_c - sum_(i=1)^m N_i.

This is the remaining support count, so it must be nonnegative. Every
remaining profile has H_p<=u_(m+1). With S_r=(Omega_r z^-L)^alpha, the
class contribution is at most

    S_r * (sum_(i=1)^m N_i H_i + R_m u_(m+1)).

For an empty tail use zero in place of its final term and check R_m=0.
This bound exceeds the complete profile-weighted majorant by at most
S_r R_m u_(m+1). Hence a stopping budget must apply after multiplying by
S_r. Budgets for separate run classes add. A zero-count profile can be
removed without spending any count mass. Reversed profiles may reuse their
equal exact count, but retain their distinct chain values.

The sorting keys must be proved upper endpoints, not interval midpoints.
Sorting then remains sound even when interval enclosures overlap. This
lemma passed independent review. Its numerical value requires a control
against the complete profile sum before replacing that sum in a verifier.

### A route-independent test for broad source-size intervals

The next range-level diagnostic can avoid profile counts entirely. Fix
z,alpha and use the current local cone caps, not a new refresh law. For
an occupancy shape ell=(ell_1,...,ell_4), write d=sum ell_i. The frozen
local operator supplies

    M_ell = [[Q_ell^alpha, A_ell^alpha],
             [E_ell^alpha, U_ell^alpha]],
    T_ell = (g^(alpha*d), U_ell^alpha)^T.

The terminal emits the last packet without a state update. Fix t,w>0
and p=(1,w)^T. Choose kappa,kappa_T to satisfy, entrywise for every shape,

    M_ell p <= kappa t^d p,
    T_ell <= kappa_T t^d p.

The least such constants are the maxima of the corresponding two row
ratios over all 4,845 occupancy shapes. Positivity, applied backward
from the terminal, bounds every conditional chain at source size r by

    kappa^32767 kappa_T t^(256r).

Here the 32,767 updates and one terminal contain exactly256r allowed
input coordinates in total. There is no terminal update or residual
potential factor: the initial state selects p_0=1. Summing any finite
set I of source sizes therefore gives the upper bound

    z^(-alpha L) kappa^32767 kappa_T
      * sum_(r in I) C(8192,r) Omega_r^alpha t^(256r).

This sum uses each size's own proved source factor. It requires no
interpolation, route independence, or uniform full-field coefficients.
At fixed z,alpha, the logarithm of the displayed expression is convex
in(log t,log w): each row ratio is a positive sum of exponentials,
and the source sum has the same form. Directed evaluation is still
required for a certificate. A failed optimization identifies restrictive
occupancy shapes; it does not show that the code has low distance.
The formula and terminal chronology passed two independent reviews.

### Where the route-independent bound necessarily fails

The preceding bound cannot replace route averaging in the middle range.
This is a limitation of that majorant, not a counterexample to the code.
Fix any admissible z,alpha,t,w and constants kappa,kappa_T satisfying its
constraints. Write q=256, g=(1+255z)/256, and L=734003.

Empty occupancy forces kappa>=1 and kappa_T>=1. At full occupancy,
the frozen live cap is exactly U_full=P g^64, where P=sum_h p_h>=1
uses the existing state-support weights. To see this, set every outside
block size to zero in the frozen output table. Each compatible product
is g_b^64; its minimum over b in{z,sqrt(z),1} is g^64. Thus

    kappa t^64 >= P^alpha g^(64 alpha).

For r<=8191, the per-size expression in the route-independent bound is
consequently at least

    C(8192,r) [eta_r beta^r P^(4r) g^(256r) z^-L]^alpha.

Indeed, its remaining factor is kappa^(32767-4r) kappa_T>=1.
Here beta=q^160/255^32>=q^128 and Omega_r=eta_r beta^r. For r>=512,
the source factor includes C(r,2), so eta_r>=C(512,2)/255^2>1.
Put A=256r. Dropping these factors only lowers the expression, leaving
the bracket q^(A/2) g^A z^-L.

If A<=L, use g>=z and z<=1; the bracket is at least q^(A/2)>1.
Otherwise put x=L/A. Minimizing over all z>0 gives the lower bound

    (1/A) log2[q^(A/2) g^A z^-L]
        >= -4 + H_2(x) + x log2(255).

For 381/1000<=x<1 this is strictly positive. The right side is concave
in x. At x=1 it is log2(255/16)>0. At x=381/1000, positivity follows
from the exact integer comparison

    255^381 * 1000^1000 > 2^4000 * 381^381 * 619^619.

The main agent checked this inequality directly; its logarithm divided
by1000 is approximately0.004593889756726406 bits. Since
floor(1000L/(381*256))=7525, every individual size512<=r<=7525 has a
route-independent majorant above one, for every admissible witness.
No numerical optimizer can repair this relaxation in that range.

There is also a concrete conditional explanation for part of this loss.
Fix a contiguous set T of r groups, with 33<=r<=2815. Each alternating
matching has at most floor(r/2)+1<=1408 occupied pair slots. Thus an
allowed route can place every occupied slot in the first44 packets of
each region. Fix such a route and any remaining setup realization.
The encoder is F-linear. Requiring zero state after each window imposes
at most256*16=4096 equations on128r source coordinates. Their solution
space has dimension at least128r-4096>0.

At most32 whole-group projections vanish identically on this space.
Every other group-zero condition defines a proper subspace, of fraction
at most1/|F|. Because |F|>r, some solution is nonzero on every such group.
It uses between r-32 and r groups and has output supported within the
windows: at most44*64*256=720896 positions, or34.375% of N.
In particular, r=2048 prevents any uniform-over-route certificate for
the entire source-size band1024--2048. The witness depends on setup;
no claim is made about a fixed message chosen before setup.

These concentrated routes have positive but potentially tiny probability.
Their existence does not refute a distance guarantee over random setup.
It explains why the proof must retain route probabilities in the middle
range. The route-independent bound can still be tested above7525, but
the existing dense certificate already starts at7552.

### Pool one internal zero/live boundary

Return-marked diagnostics identify paths containing returns as the source
of the positive size1024 bounds. The current proof splits every binary
zero/live schedule before alpha. The fixed cone permits a smaller split.

Fix the actual source set and complete route. Let M_i be the raw two-state
cone matrix for packet i, before alpha. Its live envelope is identical
at every step, so M_0 M_1 bounds the raw sum over both intermediate
binary states and all refining fine masks. Split only at pair boundaries,
then use the entrywise power(M_0 M_1)^alpha. Subadditivity gives

    (M_0 M_1)[i,j]^alpha
        <= sum_k M_0[i,k]^alpha M_1[k,j]^alpha.

At the final pair use(M_0 T_1)[i]^alpha, without a final refresh. The
common source and output factor remains(Omega_r z^-L)^alpha once.
There is no additional binary-state multiplicity. In particular, the
quiet-to-quiet alternatives Q_0 Q_1 and A_0 E_1 now remain inside one
power. Independent review confirms this fixed-route comparison.

A bounded implementation can retain the pair's total route degrees.
Let M(x,y),T(x,y) be the raw route-count polynomials, including fair
singleton orientations but no alpha. Set H=M^2 and J=MT. A pair has64
routing slots. For total degrees(f,s), put D=binom(64,f)binom(64-f,s).
Conditional on those degrees, Jensen gives the pair coefficient caps

    D^(1-alpha) ([x^f y^s]H)^alpha,
    D^(1-alpha) ([x^f y^s]J)^alpha,

where the power applies entrywise to the extracted coefficient. This
notation does not mean extraction after powering a polynomial. Compare
these caps with coefficients of the old powered-packet product, which
bounds the same conditional pair transfer. Their entrywise minimum is
also a bound. The two-packet Jensen step can lose route suppression,
so its usefulness is not implied by the fixed-route inequality.

Compose64 pair operators per ordinary region and63 followed by the
terminal pair in the final region. Retain the exact regional target
degrees and original multinomial denominator. This preserves127 ordinary
updates and one terminal emission in the final region. Required controls
are equality at alpha=1, correct terminal order, and comparison at the
same route degree. The round25 replay confirms the comparison, but the
degree-Jensen bound alone is too loose. Its minimum with the old bound
gains about223/225 bits on the two selected profiles. This is a FLOAT
diagnostic, not a directed or whole-code certificate.

### A correction that does not average occupancies before alpha

The preceding conditional Jensen step can lose more than pooling gains.
An alternative retains the old route average and subtracts a proved
lower bound on its excess. For0<alpha<=1 and x,y>=0, put c=2-2^alpha.
Then

    (x+y)^alpha <= x^alpha+y^alpha-c min(x^alpha,y^alpha).

To prove this, divide by the larger powered term and write t for their
ratio. The function(1+t^(1/alpha))^alpha is a convex norm on0<=t<=1.
Its endpoint chord is1+(2^alpha-1)t. Zero inputs follow by continuity.

For a fixed pair entry(i,j), set a_k=M_0[i,k]^alpha and
b_k=M_1[k,j]^alpha. For every rho>0,

    min(a_0 b_0,a_1 b_1)
        >= min(rho a_0,a_1) min(b_0/rho,b_1).

The product on the right is no larger than either product on the left.
The factors depend on separate packet occupancies. Sum each against
its packet's route-count weights, then convolve the resulting scalar
polynomials. At every total route degree, that convolution L_rho is
a lower bound on the correction to the old powered-packet product.
Thus old_pair-c L_rho bounds the grouped entry without route Jensen.
The terminal uses the same argument with the second packet's column T.

For a finite fixed bank of positive rho values, max_rho L_rho is still
a lower bound on the same correction. Selection per entry and degree
adds no union factor and does not depend on sampled setup. One may also
intersect this upper bound with the preceding grouped-Jensen bound.
In directed arithmetic, subtract a lower enclosure for c L_rho from an
upper enclosure for old_pair. Floating subtraction alone is diagnostic.
The scalar inequality, separation, and selection order passed independent
proof review. With the exact one-/two-input return caps, this correction
improves the two selected log2 bounds from+1105.680/+812.934 to
-80.755/-423.995. Intersecting with degree Jensen gives-101.005/-443.497.
These are FLOAT diagnostics for selected profiles, not whole-size bounds.

### Count three-input returns exactly

The return cap also loses information inside each16-coordinate block.
Fix a three-coordinate allowed set I and put B=S[:,I], where S=I_16+AC
is the actual16-by16 return matrix over E=GF256. Every C_I has rank3,
so C_I x is nonzero exactly when x is nonzero. The desired block moment is

    D_I(z)=q^-3 sum_(x in E^3, x!=0) z^wt(Bx), q=256.

All16 rows of B are nonzero, and B has rank3. The zero-row set of Bx
is therefore a proper flat of the row configuration. A flat here is the
set of all rows contained in the span of some subset of rows. Enumerate
rank1 flats by proportional rows and rank2 flats by the planes generated
by row pairs. Their closures identify the flats uniquely.

For a rank2 flat F, exactly q-1 nonzero x have that exact zero-row set.
For a rank1 flat F, the number is q^2-1 minus the counts of rank2 flats
containing F. The empty flat receives the remaining q^3-1 vectors.
Each flat contributes its nonnegative integer count at weight16-|F|.
This produces a positive polynomial without enumerating q^3 inputs.

The helper checks all560 choices of I, including rank, distinct closures,
and the first two zero-incidence moments. A separate enumeration of all
65793 projective directions checks I={0,1,2}. Take max_I D_I(z), intersect
with the earlier D3 bound, and retain the existing D1/D2 improvements.
This uniform maximum is valid without independently resampling I in
correlated blocks. Only raw return E changes, before alpha; quiet Q,
the source factor, routing, and the output-only terminal remain unchanged.

At the fixed round24 witnesses, this refinement alone changes the two
selected log2 bounds to+30.712 and-288.298. It is independent of the
binary-pair correction. Recompute the correction polynomials and the
grouped-Jensen caps from the refined raw matrix; do not reuse the larger
old return-based correction. Their composition gives selected log2 bounds
-1068.237 and-1431.422, including the degree-Jensen minimum. These FLOAT
results still do not establish whole-size or whole-code coverage.

### Directed composition and complete run-class envelopes

For directed reconstruction, first choose a coherent exact rational
majorant of every raw local entry. Compute both the old powered cap and
all correction terms from that same table. Rounding the old cap upward
and the correction downward preserves the inequality; using a correction
from a different, larger raw table need not. The exported rational pair
coefficients are chosen upper bounds, not exact moments.

The regional polynomial has nonnegative degrees and coefficients. Thus
discarding terms outside a downward-closed rectangle preserves every
coefficient inside it. One may also retain one degree exactly and fold
only the other degree modulo a positive integer. Both refinements remove
positive aliases from a fully cyclic bound. Normalize every target by
its own multinomial and marker monomial; the markers cancel on an exact
coefficient, even if different targets use a common table.

For1<=r<=4096, let R_f and V_f denote the normalized regular and terminal
regional upper bounds at actual source size r. For matching degrees(a,b), retain
the ordered chain e0(R_a R_b)^127 R_a V_b. Its reverse uses V_a, not V_b.
For a cyclic support with c runs, a+b=r-c. Therefore the sum over all
indexed supports is at most

    sum_(c=1)^r N_c max_(a+b=r-c) e0(R_a R_b)^127 R_a V_b,
    N_c = (8192/c) binom(r-1,c-1) binom(8192-r-1,c-1).

Here the maximum includes every degree-feasible ordered pair; including
pairs of actual multiplicity zero is harmless. The exact identity
sum_c N_c=binom(8192,r) checks the complete support count. This envelope
needs no individual profile-count query and leaves no residual tail.
Multiply by the common source factor(Omega_r z^-L)^alpha once, outside
the sum. The support counts are not raised to alpha.

The positive regional arithmetic also admits an integer implementation.
Represent a polynomial matrix coefficientwise by U*2^e, where every
entry of U has nonnegative integer coefficients and e is shared by the
matrix. If this majorizes its target, exact integer multiplication with
another such matrix gives a majorant at exponent e+e'. Apply the same
positive residue fold or downward-closed truncation as before. To limit
operand size, replace each retained coefficient n by ceil(n/2^h) for a
common nonnegative h, and increase the exponent by h. Since
ceil(n/2^h)*2^h>=n, this preserves the upper bound, including exact zeros.
Use separate exponents for the regular matrix and terminal vector.
No floating convolution accuracy assumption is needed. Upward rounding
can inflate the eventual bound, so its usefulness still requires an
actual complete replay; equality of rounded and unrounded DC sums is
not a valid test.

For actual source size r>4096, the same coefficient calculation can use
holes instead of active positions. This is an algebraic reindexing, not
a change to the code or to any local probability distribution. For a
matrix polynomial H of slot capacity64, define

    H_h(u,v) = u^64 H(1/u,v/u).

Its coefficient at(f_h,s) is the original coefficient at(64-f_h-s,s).
Apply the same transformation to the terminal pair J. Nonnegative local
degrees make the transformed objects polynomials. Products preserve their
order and add capacities, so H_h^64 and H_h^63 J_h reverse degrees at
regional capacity4096. At actual regional degree(f,s), where2f+s=r, the
hole degrees are(f_h,s)=(4096-r+f,s), with2f_h+s=8192-r.
The multinomial denominator is unchanged by exchanging full and empty
positions. Complementing an indexed cyclic support preserves its run
count and both matching orientations. Hence its exact run mass is also
unchanged. Use the source factor for actual r, not for8192-r. Reindex the
complete local table before any degree truncation; do not evaluate the
local caps at hole occupancies.

### A further outside-coordinate cap for at least four allowed inputs

The same matrix admits a cheap bound without enumerating larger flats.
Fix an allowed set I of size ell>=4 and let n=16-ell. The rows of A
evaluate the basis(1,a,a(a+1)/6,a^2(a+1)/6) at16 distinct AES-GF256 points.
This is an invertible change of the cubic monomial basis. Vandermonde
rank therefore makes C_I surjective onto E^4.

For uniform x in E^ell, u=C_I x is uniform in E^4. Outside I, S[:,I]x
equals A_(I^c)u. Discarding the inside-I output weights increases z^wt
for0<z<1. Each u has q^(ell-4) preimages, so the nonzero-syndrome moment
is at most

    q^-4 sum_(u in E^4, u!=0) z^wt(A_(I^c)u)
      = sum_(v=0)^(min(3,n)) binom(n,v) z^(n-v) (1-z)^v
          * (q^-v-q^-4).

Expand each output factor as z+(1-z)1_[output=0]. Any v<=4 distinct
rows have rank v. For v>=4, only u=0 satisfies all selected equations,
and that vector is excluded. This proves the positive expression and
its uniformity over I. It may be intersected with the existing D_ell
before the four-block nonquiet recurrence. The ell4 expression strictly
improves its old distance-floor branch; for larger ell, discarding more
inside coordinates can erase the gain. This is a proved local candidate,
not yet included in the round25 numerical results.

The inside coordinates can also be retained cheaply. Fix I with ell>=4,
write K_I=ker C_I, and let

    Q_ell = q^-ell sum_(x in K_I) z^wt(x),
    B_n = q^-4 sum_(u in E^4, u!=0) z^wt(A_(I^c)u),  n=16-ell.

The preceding positive polynomial evaluates B_n. The MDS rank property
makes Q_ell independent of I; it is the existing quiet-block moment.
For x sampled from E^ell, C_Ix is uniform in E^4. Conditional on each
fixed syndrome u, the inside output x+A_Iu is uniform on an affine
coset of K_I. For every affine coset v+K_I,

    E_(y in v+K_I)[z^wt(y)] <= E_(y in K_I)[z^wt(y)] = q^4 Q_ell.

To prove the inequality, expand each output factor as
z+(1-z)1_[y_j=0]. Every coefficient is nonnegative. Each requested set
of zero coordinates has either no affine solutions or the same number
as its homogeneous system. Apply this comparison separately to every
zero set and sum.

The inequality holds for every fixed u. It therefore bounds the inside
moment while retaining the u-dependent outside weight, without assuming
independence between the two parts. Averaging over u!=0 gives

    D_I(z) <= q^4 Q_ell B_(16-ell).

Intersect this bound with the existing D_ell before composing the four
blocks. For ell4 its inside factor is1. For larger ell it retains weight
information that the outside-only bound loses. Moreover,
B_n<=(1-q^-4)z^max(n-3,0), so this bound is no worse than the previous
(q^4-1)Q_ell z^max(13-ell,0) branch. Only the nonquiet return cap changes;
quiet caps, routing, source counts, and the output-only terminal do not.
This local lemma is not a whole-code distance certificate.

### Retaining the syndrome constraint for almost all-zero inside outputs

The coset bound can be loose when it treats a nonzero-syndrome affine
coset as the homogeneous coset. The actual tensor rows rule out some
of those cases. Write a row as v(a)=(1,a,b,ab), where b=a(a+1)/6.
For any point subset T, put M_T=I_4+C_T A_T. Then rank(M_T)>=2.

Indeed, M_T is symmetric, M_00 is0 or1, M_ii=1+M_0i^2 for i>0,
and M_03=M_12. Rank zero contradicts the diagonal identities. Over
GF256, a symmetric rank-one matrix has the form xx^T because every
field element has a square root. If M_00=1, then x_0=1, contradicting
M_ii=1+M_0i^2. If M_00=0, then x_0=0 and the other three x_i equal1.
This contradicts M_03=M_12. The argument uses the implemented tensor
rows, not independently sampled rows.

Fix the allowed input set I, |I|=ell>=4. In the positive zero-subset
expansion of z^wt(S[:,I]x), let T be requested zeros inside I and V
requested zeros outside I. Put t=|T|, v=|V|, k=ell-t, and n=16-ell.
For syndrome u=C_Ix, the zero conditions imply

    x_T=A_Tu,  M_Tu=C_(I\T)x_(I\T),  A_Vu=0.

If v>=4, the outside MDS equations force u=0 and contribute nothing
to the nonzero-syndrome moment. If k>=4, C_(I\T) is onto, with exactly
q^(k-4) preimages per syndrome. The normalized number of solutions with
u!=0 is therefore exactly q^-t(q^-v-q^-4).

If k<4, C_(I\T) is injective. Each feasible syndrome determines at most
one input. Project M_T onto the quotient by im C_(I\T): its rank is at
least2-k. Together with the outside constraints, the feasible syndrome
space has dimension at most min(4-v,k+2). A third bound follows from
the invertibility of S: the constrained output is supported on at most
n+k-v coordinates, and determines the input and syndrome. Thus use

    d = min(4-v,k+2,n+k-v),
    q^-ell(q^d-1)

as the normalized nonzero-syndrome count when k<4. Multiply the applicable
count by binom(ell,t)binom(n,v) z^(16-t-v)(1-z)^(t+v), and sum over
0<=t<=ell and0<=v<=min(3,n). Every term is nonnegative and bounds every
actual choice of I,T,V. This gives a uniform block cap no larger than
q^4 Q_ell B_n. Intersect with the existing cap, then rebuild the return
and pair tables. No routing independence, source factor, or encoder
operation is introduced.

### Retaining the first packet in the pair correction

The separable correction need not use one rho for every first-packet
occupancy. Fix a matrix entry of the pair product. At left occupancy
omega and right occupancy nu, write its two powered path contributions
as ab and cd. The nonnegative values a,c depend only on omega; b,d depend
only on nu. For every rho>0,

    min(rho*a,c) min(b/rho,d) <= min(ab,cd).

If a,c>0, choosing rho=c/a gives equality for every b,d. Thus much of
the loss can arise from discarding the left occupancy before choosing rho.
This observation does not require changing the random route.

Partition the left occupancies into disjoint buckets for each incoming
state row. Fix one positive rho for each bucket. Form its restricted-left
polynomial from min(rho*a,c), with the original route weights and degrees.
Form the corresponding full-right polynomial from min(b/rho,d). Multiply
these two polynomials, then sum over buckets. Every ordered route pair
occurs once in that sum. The scalar inequality therefore bounds each
coefficient by the integrated minimum of the two powered paths.

The left partition may differ between the two incoming state rows.
The implementation shares each row's partition across the two outgoing
state columns and the output-only terminal. A zero left factor contributes zero correction.
No route count or fractional-power factor is added. Selecting buckets
from approximate ratios is permissible for a diagnostic: any fixed
partition and positive rho values satisfy the inequality. A directed
replay must still reconstruct lower corrections from the same coherent
raw table used for the upper pair cap.

A finite nearest-ratio partition need not improve every coefficient over
the existing global-rho choices. Take the coefficientwise maximum of
these independently valid lower corrections before subtracting
(2-2^alpha) times that maximum from the old pair cap. Intersect the result
with the independently valid degree-Jensen cap as before. This refinement
changes the proof operator, not any encoder operation or setup distribution.

The bucket spacing has a quantitative interpretation. Suppose the ratio
r=c/a lies between the grid endpoints. If adjacent positive grid points
differ by at most a factor kappa, nearest-log assignment gives
max(rho/r,r/rho)<=sqrt(kappa). For every nonnegative b,d,

    min(rho*a,c) min(b/rho,d) >= min(ab,cd)/sqrt(kappa).

To verify this, put t=rho/r. For t>=1, the left minimum is c and
the right minimum loses at most the factor t. For t<=1, the left
minimum is tc and the right minimum is at least its value at t=1.
The two cases prove the claim. It applies pointwise before route averaging.
The dyadic grid has kappa=2. The rational grid
2^e times {1,5/4,3/2,7/4}, for e=-16,...,15, followed by2^16,
has kappa=5/4. Its guaranteed fraction is therefore sqrt(4/5).
Clipped ratios retain the original valid lower bound, but not this
quantitative guarantee. Zero factors contribute zero.

Finer nearest-log assignment is not automatically coefficientwise stronger.
Intersect the completed upper cap with the prior completed upper cap.
The round28 implementation performs this intersection, retains all right
occupancies, and changes no probability distribution or source count.
Independent source review and exact rational toy checks pass; its numerical
whole-chain replays remain FLOAT diagnostics until directed reconstruction.

### Exact correction between different ratio buckets

More of the same minimum can be retained without enumerating all pairs
of occupancies. For positive a,c,b,d, define r=c/a and s=b/d. Partition
both ratios using the same ordered intervals. If r lies in bucket i
and s in bucket j<i, then ab<=cd. If j>i, then cd<=ab. Hence the
minimum is exact on every off-diagonal bucket rectangle. Use the existing
rho_i separable lower bound only when i=j.

Let A_i,C_i be the route polynomials for a,c restricted to left bucket i.
Let B_j,D_j be those for b,d restricted to right bucket j. Put
L_i=min(rho_i*a,c) on that left bucket and R_i=min(b/rho_i,d) on the
same right bucket, retaining the route weights and occupancy degrees.
The correction polynomial is

    sum_i [ A_i * sum_(j<i) B_j + C_i * sum_(j>i) D_j + L_i * R_i ].

Products are ordinary positive polynomial convolutions. Prefix and suffix
sums avoid enumerating individual occupancy pairs. The bucket rectangles
are disjoint, so each ordered pair contributes exactly once. Its correction
is at least the previous one-sided rho witness and at most min(ab,cd).
Zero left or right minima contribute zero and may be omitted.

Use separate partitions for each of the two left rows and each of the
three right columns, including the output-only terminal. All partitions
use the same interval boundaries. Clipped endpoint buckets remain ordered;
no quantitative ratio-spacing guarantee is needed for the off-diagonal rule.
The argument introduces no independence assumption or extra source factor.
Independent mathematical review passes.

In a directed implementation, ordering must refer to the same coherent
powered caps used in the base upper bound. Midpoint comparisons alone do
not certify order for overlapping intervals. Use certified ratio intervals
or a coherent exact representation, and reconstruct every subtracted
correction with lower bounds. FLOAT tests of this lemma are not certificates.

One directed construction can avoid interval-order ambiguity. Freeze each
powered local cap upward to an exact nonnegative rational, and rebuild
both the uncorrected pair sum and the correction from that table. For
powered path products P,Q, define

    g(P,Q)=max(P,Q)+(2^alpha-1)min(P,Q),  0<alpha<=1.

The scalar concavity inequality bounds the pooled raw-path target by g.
This function is nondecreasing in both arguments. Upward freezing gives
Pbar>=P,Qbar>=Q. Any L<=min(Pbar,Qbar) then satisfies

    Pbar+Qbar-(2-2^alpha)L >= g(Pbar,Qbar) >= g(P,Q).

Consequently exact rational ordering is sound even though the frozen
atoms need not equal the original powers. Use a lower enclosure for the
subtracted product. Apply the argument to all six pair channels before
positive route averaging and composition. An independently valid Jensen
cap for the original pair target may still be intersected entrywise.
The directed local constructor now implements this plan with exact
integer route polynomials and upward2^-256 powered atoms. Its rational
tests and three independent source reviews pass. Selected regional
192-/256-bit controls also pass. A complete directed regional/run join
now certifies the entire size2048 class at35% with86.7425 bits. Other
source sizes remain open; neither local correctness nor this single-size
certificate establishes a whole-code bound. See round29 in the goal ledger.

The ordered correction also dominates every single-rho witness whose
rho belongs to the same nearest-log grid. For one cell centered at rho,
divide the witness by a*d and write r=c/a, s=b/d. The center witness is
min(r,s,rho,rs/rho). If rho_minus is the preceding grid point, the lower
cell boundary gives r,s>=sqrt(rho_minus*rho), hence rs/rho>=rho_minus.
Any smaller grid point yields a witness at most rho_minus. The upper
boundary similarly gives r,s<=sqrt(rho*rho_plus); any larger grid point
yields a witness at most rs/rho_plus, and each of the four center terms
is at least that value. At clipped endpoint cells only the applicable
side is needed. Different cells already give the exact minimum.

Apply this argument to each occupancy pair before route averaging.
It permits omitting redundant global-grid products only when every
candidate uses the same frozen powered atoms and a subset of this grid.
It does not compare coefficients built from independently rounded tables.

There is a further unimplemented refinement. For powered path ratio w>=1,
the exact deficit per unit of the smaller path is

    h_alpha(w)=w+1-(w^(1/alpha)+1)^alpha.

For0<alpha<1, its derivative is
1-[w^(1/alpha)/(w^(1/alpha)+1)]^(1-alpha)>=0. It increases from
2-2^alpha at w=1 toward1. An ordered bucket rectangle that certifies
w>=w_min may therefore use h_alpha(w_min), instead of the universal
coefficient2-2^alpha, on its minimum contribution. A few fixed separation
thresholds permit positive incremental prefix/suffix corrections. This
could recover loss remaining after exact off-diagonal minima.

If powered caps are frozen upward, justify this variant by monotonicity
of the exact pooled norm (P^(1/alpha)+Q^(1/alpha))^alpha. Do not assume
the bucket-dependent formula itself is monotone. Its ordering and deficit
coefficient must be certified for the same frozen table. This candidate
has an independently checked scalar argument but no numerical replay yet.

### Share the regional calculation across a source-size band

Fix one z, alpha, and output-weight threshold L. Let B be a finite set
of actual source sizes in1..4096. For each r in B, write
Omega_r=eta_r beta^r. Suppose R_(r,f) and V_(r,f) are valid normalized
regional bounds for actual full-pair count f, including that target's
own coefficient and multinomial/marker denominator. The ordered chain is

    Phi_r(a,b) = e0 (R_(r,a) R_(r,b))^127 R_(r,a) V_(r,b).

It contains255 regular matrices and one output-only terminal vector.
Set s_r=Omega_r^(alpha/256), and define entrywise envelopes

    Rbar_f = max_(r in B) s_r R_(r,f),
    Vbar_f = max_(r in B) s_r V_(r,f).

An invalid marginal index contributes a zero matrix/vector to these
maxima. Scaling the256 factors contributes exactly Omega_r^alpha.
Positivity therefore bounds Omega_r^alpha Phi_r(a,b) by the same ordered
chain formed from Rbar,Vbar. This does not assume a common maximizing
source size in different entries or factors.

Let N_r(c) count indexed cyclic supports of size r with c runs. Set it
to zero outside1<=c<=min(r,8192-r). Each actual support has
d=a+b=r-c. Define C_d=sum_(r in B) N_r(r-d). The resulting band bound is

    z^(-alpha L) sum_d C_d
      max_(a+b=d) e0 (Rbar_a Rbar_b)^127 Rbar_a Vbar_b.

The maximum ranges over ordered pairs; no additional factor of two is
needed. Counts are exact and stay outside alpha. In particular,
sum_d C_d=sum_(r in B) binom(8192,r). A singleton band recovers the
existing run-envelope bound, up to directed arithmetic.

At fixed z, alpha, and markers, the local and marked regional polynomials
do not depend on r. Thus all needed targets can be extracted while a
table is live, before its arrays
are discarded. Already-normalized size2048 entries alone cannot recover
neighboring singleton coefficients. Each neighbor still needs its own
extraction and denominator. Do not infer a band from the size2048 margin.

For hole coordinates, the regrouping changes: if
f_h=4096-r+f, then d_h=8192-r-c and C_d=sum_r N_r(8192-r-d).
The source factor remains Omega_r. Never combine hole-indexed matrices
with the actual-full-index formula d=r-c. The first proposed test avoids
this complication by restricting to actual sizes2047--2049.

The lemma has two independent mathematical reviews. Exact small examples,
singleton reproduction, and the count identity pass. A selected2047--2049
control loses9.797 bits relative to the sum of its individual bounds.
The complete directed band replay closes35% with74.2104 bits, including
the new(1024,1024) endpoint. The loss permits different source
sizes in different entries and factors. It can also combine one size's
large conditional bound with another size's large exact count.

### Keep run counts with their own source size

For a fixed d, let B_d contain exactly the sizes r with N_r(r-d)>0.
Choose positive weights p_(r,d) for r in B_d, fixed across all ordered
profile pairs with this d. Set

    q_(r,d) = (Omega_r^alpha N_r(r-d) / p_(r,d))^(1/256).

Form entrywise envelopes Rbar_(d,f),Vbar_(d,f) from q_(r,d) times each
own-r marginal for which both f and its partner d-f are valid. This
restriction retains every actual ordered chain. Each chain contains256 factors, so
Omega_r^alpha N_r(r-d) Phi_r(a,b) is at most p_(r,d) times its envelope
chain. Maximizing over a,b with a+b=d and then summing over r gives

    sum_(r in B_d) Omega_r^alpha N_r(r-d) max_(a+b=d) Phi_r(a,b)
      <= (sum_(r in B_d) p_(r,d))
           max_(a+b=d) e0(Rbar_(d,a) Rbar_(d,b))^127
                           Rbar_(d,a) Vbar_(d,b).

The final sum over d is multiplied by z^(-alpha L) once. Counts enter
q with exponent1/256, not alpha/256. Zero-count sizes are omitted;
invalid marginals are zero only for that source size. The profile maximum
retains both terminal orders. Weights may depend on d and the deterministic
proof tables, but cannot change separately with each candidate profile
and then be discarded under this maximum.

Taking every p_(r,d)=1 absorbs the counts before the matrix maxima and
leaves the factor |B_d|. Taking p_(r,d)=N_r(r-d) is no larger than the
previous band bound: it excludes zero-count sizes and invalid partner
indices. It reproduces that bound when both use the same eligible
region set. Neither weighting is assumed to dominate the other.
The minimum of their completed degree bounds is also valid. This lemma
has two independent mathematical reviews, and its exact small controls
pass. At the selected pair(895,895), the2040--2056 count-aware envelope
improves the directed log2 bound from+11.1082 to-65.2459. The individual
controls sum to-69.1517, so3.906 bits of envelope loss remain. That
selected control alone is not a complete degree or band bound. The
subsequent31-cell directed replay covers the entire2040--2056 band
jointly at35% with60.2491326 bits, including all1058840 ordered envelope
pairs and every state path. All own-size normalizers and exact count
identities pass. This remains a support-size band at fixed K=2^20,
not whole-code coverage. No encoder operation changes.

### Direct all-singleton boundary at actual size4096

Fix actual source size r=4096 and one matching with f=0 full pairs.
Then s=4096 pairs are singleton. Let H and J be the regular and
output-only terminal pair polynomials used by the current ordered bound.
Every monomial has nonnegative degrees f,s with f+s<=64. A region
contains64 pair factors. Its target degree(0,4096) therefore forces
degree(0,64) in every factor. Consequently its normalized matrices are

    R_0 = H_(0,64)^64,
    V_0 = H_(0,64)^63 J_(0,64).

The multinomial normalizer is binom(4096,0)binom(4096,4096)=1.
Marker powers cancel exactly. Singleton orientations have already been
averaged in the local coefficients; no factor2^4096 belongs here.
This extracts the boundary of the selected pair majorant exactly,
not the true encoder enumerator. Arithmetic used to evaluate that
majorant must still be identified as FLOAT or directed.

If both matchings have f=0, the source support alternates around the
8192-cycle. There are exactly two indexed supports. Their finite chain
reduces to e0 H_(0,64)^16383 J_(0,64), with no final refresh. The general
run envelope already counts those two supports at c4096; a separate
endpoint diagnostic must not be added to it again. This boundary
calculation removes one extreme marker/projection problem without
changing the code or substituting full-field randomness.

### Exact near-boundary coefficients in full/empty coordinates

The all-singleton argument extends to nearby profiles without a large
coefficient table. Replace the local singleton degree s by the empty-pair
degree e=64-f-s. For either regular or terminal pair polynomial, define

    G(X,Y) = sum_(f,s) H_(f,s) X^f Y^(64-f-s),
    K(X,Y) = sum_(f,s) J_(f,s) X^f Y^(64-f-s).

The regional coefficients of G^64 and G^63 K at (f,e) equal those of
H^64 and H^63 J at (f,4096-f-e). An actual source size r has
e=4096-r+f and s=r-2f. Divide these coefficients by their own exact
multinomial 4096!/(f!e!s!). No marker, new alpha, source factor, or
singleton-orientation factor is introduced.

All degrees are nonnegative. Therefore truncating every intermediate
product to a fixed rectangle 0<=f<=F,0<=e<=E preserves every coefficient
inside that rectangle exactly. Six positive matrix doublings suffice;
the terminal update remains left multiplication, giving G^63 K. At
r4096,f1, the four-coefficient rectangle F=E=1 includes both a single
factor with one full and one empty pair and two separate factors carrying
those defects. Their positions and noncommuting matrix order must remain.
At r4095,f0, the target is (0,1), not the r4096 all-singleton endpoint.

This is exact extraction of the chosen local majorant, not an assertion
that the majorant is the encoder's exact enumerator. Directed arithmetic
is still required for a numerical certificate. The bounded implementation
`cheap_lane_edge_regions_directed.py` has independent source audits;
numerical replays are recorded separately in the goal ledger. This
paragraph alone adds no support-size coverage.

### Select ordinary witnesses entrywise

For a fixed actual size and matching target, each positive marker gives
a valid bound H(u,v)^64/D and H(u,v)^63 J(u,v)/D, with that target's
own denominator D. Their entrywise minima are still valid bounds for
the same regional matrix and terminal vector. Positivity preserves
these inequalities through the finite ordered chain. The six entries
may choose different witnesses; no extra marker union factor is needed.
Witnesses are deterministic proof parameters, not setup randomness.

With directed arithmetic, retain the candidate enclosure having the
smallest certified upper endpoint. These are different majorants, not
enclosures of one unknown probability, so do not intersect their intervals.
At actual size4096, all denominators for one marker can be generated from

    D_0=v^4096,
    D_(f+1)/D_f=(u/v^2)(4096-2f)(4095-2f)/(f+1)^2.

Stop at f2048 before advancing through a zero numerator. An exact endpoint
projection may replace its own target afterward. An entrywise selection
result does not establish the performance of a whole-marker cell bound
over neighboring sizes; that is a separate envelope below.

### Ordinary marker cells over an actual-size band

This band lemma can exploit slack in an ordinary-Cauchy
bound. Fix z,alpha and a positive marker bank (u_i,v_i). Use either actual
coordinates k(r)=r for a band at most4096, or hole coordinates k(r)=8192-r
for a band above4096. Require1<=k(r)<=4096 and the established source-bound
domain throughout. This excludes the full-support endpoint r8192, which
has no runs and requires its separate dense-support bound.
In the latter coordinate convention reverse the complete local pair
table first. Do not mix the coordinate conventions inside one band.
The source factor is always Omega_r=eta_r beta^r for the actual size r.

Let B_i=H(u_i,v_i)^64 and V_i=H(u_i,v_i)^63 J(u_i,v_i). These ordinary
numerators do not depend on the regional target. For the ordered marker
pair(i,j), define

    h_ij = e0 (B_i B_j)^127 B_i V_j,
    D_i(r,a) = 4096!/[a!(k(r)-2a)!(4096-k(r)+a)!]
                 * u_i^a v_i^(k(r)-2a),
    Dtilde_i(r,a) = D_i(r,a) beta^(-alpha r/256).

The256 regional factors contain128 copies of each matching denominator.
Thus the source-scaled conditional chain is bounded by
eta_r^alpha h_ij/[Dtilde_i(r,a)Dtilde_j(r,b)]^128. The output factor
z^(-alpha L) remains outside. Both ordered terminal numerators matter.

Assign every valid integer marginal to a marker interval I_i. For an
actual-size integer interval B and each interval pair, take the polytope
of(r,a,b) with r in the closed real interval spanned by B,
a in I_i, b in I_j, 2a<=k(r),2b<=k(r), and
a+b<=k(r)-1. The last inequality excludes zero-run profiles. The
nonnegative empty counts are automatic when k(r)<=4096 and a,b>=0.
There are nine bounding planes. Every vertex can be found by solving
independent triples of planes exactly and checking all inequalities.
Empty cells contribute nothing.

Extend factorials by Gamma. The logarithm of each multinomial is jointly
concave: it is a constant minus three convex log-Gamma functions of
affine nonnegative counts. Marker and source-tilt terms are affine.
Hence the minimum of log(Dtilde_i Dtilde_j) occurs at a polytope vertex.
Let m_ij be a directed positive lower bound for the product at all these
vertices. Fractional vertices must not be rounded. Directed log-Gamma
evaluation suffices; alternatively the convexity bound

    log Gamma(n+theta+1) <= log(n!) + theta log(n+1), 0<=theta<=1,

gives a lower denominator when substituted with its negative sign.
The entire band then has the valid bound

    z^(-alpha L) * [sum_(r in B) eta_r^alpha binom(8192,r)]
                 * max_(i,j with nonempty cell) h_ij/m_ij^128.

The binomial already counts every indexed support; there is no additional
profile count, orientation factor, or band length. This bound drops the
correlation between run count and chain value, so it need not be tight.
Split disjoint size/run cells if necessary, retaining their exact run
counts. A target-specific endpoint replacement is not a constant
numerator and cannot be inserted into h_ij without a separate argument.
The lemma has independent mathematical reviews. The round33 fixed-cell
replay closes actual sizes6144--6160 jointly, but its wider intervals fail.
The ordinary size6144 result does not justify interpolating other sizes.

### Keep marker cells aligned with relative occupancy

A fixed full-pair interval can change meaning as the source size changes.
For example, a1023 is almost entirely full at hole size2048, but leaves258
singletons at hole size2304. A near-full marker is unsuitable for the
latter target. This mismatch changes the bound, not the construction.

For hole size k>0 and full-pair count a, define x=2048a/k. Then x lies
in [0,1024]. The reference marker bank was calibrated at k2048, where
x=a. Replace each internal integer-cell boundary h by the rational cut
h+1/2; use0 and1024 as the outer boundaries. The resulting closed
intervals cover every real x, including those arising from odd k.
Their shared endpoints cause no extra count: the band bound uses a
maximum over cells and counts each source support only in its outer sum.

For normalized cell [ell,u], replace its fixed full-count constraints by

    ell*k <= 2048a <= u*k.

Do the same for b. These are linear inequalities with rational
coefficients. Together with the k bounds,2a<=k,2b<=k, and a+b<=k-1,
they define a bounded nine-plane polytope. The preceding log-concavity
argument and exact vertex enumeration therefore still apply. The
denominator is evaluated at the actual rational vertex(k,a,b), not at
the reference value x. Neither the marker values nor the ordinary
numerators change. Source factors still use actual size r=8192-k.

This alignment is a proof-witness choice. It adds no setup randomness,
encoder operation, coefficient approximation, or transferred certificate.
Numerical band coverage must be recorded separately after directed replay.

### Actual coordinates and cell-specific run counts

The same normalized geometry applies below the midpoint. For an actual
support size 1<=r<=4096, set k=r and retain the unreversed local H/J
tables. Every matching marginal a satisfies 0<=a<=r/2. The two
matchings have c=r-a-b positive cyclic runs, so a+b<=r-1. Empty
matching slots number 4096-r+a and are automatically nonnegative.
Thus no additional capacity plane is needed. Evaluate the source tilt
as -alpha*r*log(beta)/256, not as -alpha*(8192-r)*log(beta)/256.
The log-concavity and exact-vertex argument is otherwise unchanged.

The binomial-times-maximum band bound forgets that extreme matching
profiles may be rare. Retain part of this information as follows. Fix
one coordinate convention and a band of actual sizes. For every r,
set k=r in actual coordinates or k=8192-r in hole coordinates.
Let A_i(k) be the integer full-pair counts assigned to normalized cell i.
Use the existing tie-to-lower-cell convention, so the A_i(k) partition
0,...,floor(k/2). For a cell with cuts [ell_i,u_i], its upper endpoint
is min(floor(k/2),floor(k*u_i/2048)). Its lower endpoint is zero in
the first cell, and floor(k*ell_i/2048)+1 otherwise. Empty intervals
receive no assigned targets.

For nonempty A_i(k),A_j(k), define the integer interval

    C_ij(r) = [k-max A_i(k)-max A_j(k),
               k-min A_i(k)-min A_j(k)] intersect [1,k].

Every indexed source support assigned to ordered cell (i,j) has its run
count in C_ij(r). The same c applies to a support and its complement.
The exact number of indexed cyclic supports with actual size r and c
runs is

    N_r(c) = (8192/c) binom(r-1,c-1) binom(8192-r-1,c-1).

Consequently M_ij(r)=sum_(c in C_ij(r)) N_r(c) is an upper bound on
the number assigned to the cell. Set M_ij(r)=0 when either marginal
interval or the run interval is empty. It need not equal the cell's
actual count; several cells can charge the same run class.

Let U_ij=h_ij/m_ij^128 be the earlier beta-tilted conditional upper
bound for that ordered cell over the whole size band. Define

    W_ij = sum_(r in band) eta_r^alpha M_ij(r).

Partition actual supports by their assigned ordered cell, apply U_ij
inside each part, and enlarge each part's count to M_ij(r). Positivity
gives the complete band bound

    z^(-alpha L) sum_(i,j) W_ij U_ij.

No independence is needed. Off-diagonal cells have equal count bounds
W_ij=W_ji, but their two terminal orders must still be added separately.
Diagonal cells appear once. There is no extra profile factor or power
of alpha. The overlap among run intervals only enlarges the count.
The sum of these count bounds is at most 961 times the binomial count;
this is not a bound of ten bits on the final weighted loss against exact
joint-profile counting. A run class can be charged to a cell with a
larger conditional bound even when no support realizes that assignment.

Compute N_r(c) by the existing exact recurrence, using the complementary
weight when r>4096. Its prefix sums answer all interval queries; the
checksum remains sum_c N_r(c)=binom(8192,r). Stream one actual size at
a time, retaining its own eta_r. The minimum of this completed band
bound and the earlier completed binomial-times-maximum bound is valid.
Do not mix factors from their different aggregations. Independent
source reviews validate the refinement; numerical gains require a
separate directed replay.

### Keeping run counts and conditional bounds together

The preceding count can still be loose within a cell. Its largest run
mass and largest conditional bound may occur at different run counts.
Fix the local table, positive markers, alpha, and z. Partition the
integers1,...,k_hi into consecutive disjoint intervals C_s=[l_s,u_s].
An interval may contain just one integer. No setup randomness changes.

For each cell (i,j), intersect its existing polytope with

    l_s <= k-a-b <= u_s.

These two affine constraints preserve boundedness and concavity of the
log denominator. Enumerating all independent triples of the eleven
planes therefore finds every vertex, including lower-dimensional cases.
Use the resulting denominator lower bound to obtain U_ij,s for that
ordered cell and run interval. Keep fractional vertices exact.

For each actual size r, define

    M_ij,s(r) = sum_(c in C_ij(r) intersect C_s) N_r(c),
    W_ij,s = sum_(r in band) eta_r^alpha M_ij,s(r).

Each assigned support belongs to exactly one integer interval C_s.
Consequently sum_s W_ij,s U_ij,s bounds that cell's full contribution.
The count intervals can overlap between cells, as before. A positive
M_ij,s(r) requires a nonempty restricted polytope: an integer run count
in the interval can be realized by integer a,b in the assigned marginal
ranges. This does not assert that every such profile has a nonzero
indexed-support count.

For an unordered off-diagonal cell, add both terminal orders inside
every interval before summing. A diagonal cell contributes once. The
minimum of this completed cell bound and the preceding completed cell
bound is valid. Sum these minima over refined cells and keep all
unrefined cells in the sum. Finally apply z^(-alpha L) once; the earlier
complete binomial-times-maximum bound can also be minimized globally.
There is no additional interval-count factor or alpha power.

Source factors always use actual r, including in hole coordinates.
The exact interval masses must sum to the original cell mass for each r.
This identity checks the partition; it does not by itself validate the
conditional bounds. Numerical gains from this refinement remain a
separate replay obligation.

### Optional exact count for matching imbalance

The preceding run count forgets the difference between the two matching
occupancies. That difference can be counted without changing the code.
Fix an even cycle length n, a support size0<r<n, and a run count
1<=c<=min(r,n-r). Positions are indexed0,...,n-1; supports are not
identified up to rotation.
Let a count full edges starting at even positions and b count full edges
starting at odd positions. Set D=a-b, so a+b=r-c. Write N(r,c,D) for
the number of these indexed supports with the specified values.

Track the parity of the start of each positive ones run. A run of even
length contributes+1 to D from an even start, and-1 from an odd start.
It preserves parity. An odd-length run contributes0 and changes parity.
Its following positive zero gap changes parity exactly when its length
is odd. Mark ones by x, zeros by y, and imbalance by t. The two matrices
of formal power series are

    R(x,t) = (1-x^2)^(-1) [[x^2 t, x], [x, x^2/t]],
    G(y)   = (1-y^2)^(-1) [[y^2, y], [y, y^2]].

The row and column order is even,odd. Multiplying R G advances through
one run and its gap. Taking the trace after c such pairs enforces return
to the starting parity. Therefore

    N(r,c,D) = n/(2c) [x^r y^(n-r) t^D] tr((R(x,t) G(y))^c).

Indeed, there are n/2 possible indexed starting positions of each parity.
Every indexed support has exactly c choices of a designated ones-run
start, including periodic supports. Dividing this count by c proves the
identity. At t=1, the coefficient of x^r y^(n-r) in the trace is
2 binom(r-1,c-1) binom(n-r-1,c-1), recovering the previous run count.
The identity also gives symmetry under D->-D, by rotating once. The
constraint D congruent r-c modulo2 follows from a+b=r-c.

For real0<x,y,t<1 and an integer d>=0, positivity gives the count bound

    sum_(D<=-d) N(r,c,D)
      <= n/(2c) x^(-r) y^(-(n-r)) t^d tr((R(x,t) G(y))^c).

For such D, t^D>=t^(-d); evaluating the positive series also upper-bounds
the selected x,y coefficient. The minimum with the exact total run count
is valid. This is an independently audited counting identity, not a
numerical certificate. No implementation or numerical gain is claimed.
It remains available if the complete coefficient calculation needs a
sharper count of asymmetric profiles.

### A common numerator for all regional targets

The regional calculation can also supply a bound that is independent of
the target degrees. Fix the exact nonnegative local pair table and positive
rational markers u,v. Mark a degree(f,s) by u^f v^s before multiplication.
Let P be one entry of either H^64 or H^63 J. Reduce both degree coordinates
modulo256, using positive coefficient folding and upward integer rounding
at each multiplication. Write the resulting entry as

    P_hat(x,y) = 2^e sum_(0<=a,b<256) n_ab x^a y^b,
    M = 2^e max_(0<=a,b<256) n_ab.

Here n_ab are nonnegative integers. The exponent e belongs to this matrix;
the regular and terminal matrices may have different exponents. Every
exact marked coefficient of P is at most its corresponding coefficient
in P_hat, and hence at most M. This follows from positivity at folding
and from entrywise upward rounding, not from uniformity of the residues.
The maxima of different entries need not occur at the same residue.

For an actual target(f,s), divide M by its own denominator

    binom(4096,f) binom(4096-f,s) u^f v^s.

The four regular maxima and two terminal maxima therefore give an
entrywise upper bound for every normalized regional matrix. Positive
matrix multiplication preserves this inequality. They can replace the
ordinary unnormalized numerators in the existing cell and run-interval
proofs. Those proofs must still cover every actual target and retain its
source factor, support count, denominator, and terminal ordering. There
is no extra65536 factor, alpha power, or independence assumption.

The scan has a simple exact control: if D_hat is the represented cyclic
DC mass, then M<=D_hat<=65536 M. Thus this particular replacement gains
at most16 bits per regional factor relative to its represented DC mass.
Rounding can make D_hat exceed the unrounded ordinary DC mass, so their
comparison requires a separate directed check. A common numerator is a
proof tool, not a new encoder or a whole-code certificate. Its numerical
tightness and any resulting interval coverage are recorded separately.

### Selecting regional witnesses without changing support cells

The normalized cells assign supports to counting classes. They need not
also determine the only markers used to bound those supports. Fix the
local table, z, alpha, and a verified bank of positive rational markers.
For marker p, let R_p,V_p be its four regular and two terminal numerator
bounds. These are unnormalized universal coefficient bounds as above.

Fix an assigned cell pair(i,j), with its exact feasible polytope P_ij
of actual-size/full-count triples(r,a,b). For any candidate markers p,q,
define

    L_pq = min_(r,a,b in P_ij)
        [log Dtilde_p(r,a) + log Dtilde_q(r,b)],

where Dtilde includes the existing beta tilt. The concavity argument still
puts this minimum at an exact rational vertex. Let h_pq^forward and
h_pq^reverse be the two finite-chain numerators from(R_p,V_p),(R_q,V_q).
For either required terminal order, the completed conditional bound is

    U_pq^order = h_pq^order exp(-128 L_pq).

Every such bound covers the same assigned cell and retains the common r
in both denominators. Therefore the minimum over any finite candidate set
is another valid bound for that cell and order. This is selection among
upper bounds, not an intersection of enclosures of one unknown value.
Directed implementation selects by certified upper endpoints and retains
the selected enclosure. Including the previously assigned pair ensures
that its upper endpoint cannot increase.

Keep the support cells and indexed run-interval counts unchanged. Add both
terminal orders off diagonal and count diagonal cells once, as before.
Then sum the cell contributions and apply the output factor once. Candidate
markers require no union-bound multiplier: they do not define additional
events or change setup randomness. Never identify the two terminal
numerators merely because the support cell is diagonal; different markers
can yield different ordered chains. Numerical gains require a complete
replay with all cells retained.

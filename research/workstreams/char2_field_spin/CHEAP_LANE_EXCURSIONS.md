# Cheap lane mixing: quiet excursions and route concentration

This note studies the unchanged [cheap-lane candidate](CHEAP_LANE_MIXING_20261016.md).
It explains a concrete loss in the current proof bound, then separates a
rare routing event from the remaining outer-code problem. These results
do not supply a whole-code distance certificate. The certified 3.8-ms
baseline and its implementation remain unchanged.

## What the decomposition identifies

Fix the existing artificial 256-group diagnostic at z=.85, alpha=.523,
with four-packet pooling. Its repeated fully-paired regional profile is
not a feasible whole-chain source class. Compare one regional transfer
entry, with full 16-coordinate state support at both ends:

| Allowed intermediate supports | Log2 expression |
|---|---:|
| Full support throughout | -1000.287280 |
| Any nonempty support | -943.647442 |
| Including zero support | -525.633675 |

Restrictions are imposed before pooling and taking fractional powers.
These are diagnostic comparisons, not probabilities of complementary
events. Nevertheless, they locate the main loss: allowing the state to
stop and restart changes the expression much more than partial support
while the state stays nonzero.

The full-state local cap is already close to its exact filtered output
moment. For the tested sparse occupancies, it exceeds z^64 by only
0.063627 bits per packet. Further tuning of that local cap is therefore
unlikely to remove the much larger gap.

## A real, rare route that admits low-weight words

There are 8192 source groups, each with 128 F-valued source coordinates,
where F=GF(2^128) has the tower representation fixed in the construction.
Each of 256 regions contains 128 packets, or 4096 routing-pair slots.
Its permutation is independent and uniform. The two alternating pair
matchings each serve 128 regions.

Fix the consecutive group set T={0,...,255}. It occupies 128 full pairs
in the even matching. In the shifted matching it occupies 127 full pairs
and two singleton pairs. Consider the event that every occupied pair
lands in the first 44 packets of its region. Its exact probability is

    [C(1408,128)/C(4096,128)]^128
    [C(1408,129)/C(4096,129)]^128.

Conditional on that route, impose zero state after each 44-packet window.
These impose at most 256*16=4096 F-linear equations on the 32768 source
coordinates supported in T. Their kernel has dimension at least 28672.
Every word in this kernel has zero output outside the windows: the input
is zero there, and each preceding window leaves zero state. Its output
weight is at most 44*64*256=720896, or 34.375% of N.

Thus this route necessarily admits a nonzero word below 35%, for every
realization of the other setup variables. The word need not use all 256
groups. At most 32 group projections can vanish identically on the
kernel. Since |F|>256, a union of the remaining proper group-zero kernels
cannot cover it. Some witness consequently uses between 224 and 256 groups.

The route event is extremely rare. Floating evaluation gives log2
probability -52149.139093 for this fixed T and window placement. A union
over all 256-group sets and all independently located 44-packet windows
has log2 upper bound -48870.140249. The latter calculation uses the fact
that a proper group set cannot consist of full pairs in both matchings.
These floating values illustrate the mechanism; the next bound uses
exact integer arithmetic.

The analogous fixed-window term in the artificial scalar majorant is
positive, at +73742.613560. It retains the source-count prefactor and
counts a large space of potential witnesses where event probability
cannot exceed one. This comparison does not identify the artificial
mask path with the actual kernel witnesses. It demonstrates why rare
routes should be bounded as events, rather than charged through their
possibly enormous number of codewords.

## An exact bound for all concentrated allowed-position sets

For a source-group set T, define E_T to mean that all its allowed pair
positions fit within one contiguous 44-packet window in each region.
The window may differ between regions. This concerns every position
allowed by T, not just the nonzero entries of a particular outer word.

If E_T holds for |T|>=7, it also holds for every seven-group subset of T.
Fix one such subset. It occupies at least four pairs in each matching;
select four of those pairs before sampling the route. Their destinations
are a uniform four-subset of 4096 slots.

If these destinations fit within 1408 consecutive slots, one of them
anchors an interval containing the others among its next 1407 slots.
Allowing cyclic intervals enlarges the event. Union over the four anchors
and sample the remaining destinations without replacement. This gives

    p <= 4*C(1407,3)/C(4095,3) <= 4*(67/195)^3.

Independence across regions and a union over seven-group subsets give

    Pr_route[exists T, |T|>=7, E_T]
      <= C(8192,7) * [4*(67/195)^3]^256
      <= 2^-592.

The final comparison is exact. The helper verifies both integer
comparisons that place the displayed rational bound between 2^-593 and
2^-592. No logarithmic rounding enters this certificate.

This closes one routing event across all source sizes. It does not close
the nonquiet 7--7583-group distance classes. A word can vanish on allowed
positions outside its live intervals, or have nonzero inputs in the
feedback kernel. Such words do not imply E_T.

## Retain outer constraints across quiet intervals

The remaining issue is not only where inputs may land, but which of them
can satisfy the outer code while the inner stays quiet. The following
necessary-condition bound retains that information without analyzing
the live intervals.

Fix T with r=|T|<=256, a route, and a prescribed set Q of regular packets.
Require both endpoint states to be zero on every packet in Q. Do not
impose the unused update after the final output. For packet t, let I_t
be its allowed input positions and ell_tb the count in emission block b.
Invertibility of the refresh forces C X_t=0 on each such packet. Define

    D_Q := sum_(t in Q) rank_E C_(I_t)
         = sum_(t in Q,b) min(4,ell_tb),
    eta_r := 256^2 / [255^2*(257-r)],
    beta := 256^160 / 255^32.

The equality uses the four direct RS[16,4] maps. The equations concern
disjoint routed input coordinates, so their ambient ranks add.

For each sampled outer, the source kernel of these equations is defined
over E=GF256. An exact-T F-valued witness implies that every group
projection of this E-kernel is nonzero. The reviewed hyperplane bound in
the [rank interface](CHEAP_LANE_RANK_INTERFACE.md) compares its projective
count to the number n_T of scalar kernel vectors nonzero on every group.
The outer-density bound then gives

    Pr_outer[an exact-T witness satisfies the Q equations | route]
       <= min(1, eta_r * beta^r * 256^(-D_Q)).

This cap is taken before averaging the route. For any fixed 0<alpha<=1,
it implies the simpler bound

    (eta_r * beta^r)^alpha * E_route[256^(-alpha*D_Q)].

Here alpha can be selected for this necessary-condition event. It need
not equal the exponent used by the previous output-weight operator.
Support and window unions remain outside the expectation. Q is a class
index fixed before setup; it is not conditioning on a previously observed
quiet trajectory. A packet preceding a live-state interval may activate
the state and must not be declared quiet merely because its incoming
state is zero.

### Exact closure for consecutive 256-group supports

Fix T to be any consecutive set of 256 groups on the cycle. Consider
words with exact source support T for which all regular steps having
either endpoint state nonzero fit within one 44-packet window per region.
There are 8192 choices of T and 85 window choices per region. This class
does not require a bound on output weight; controlling it supplies a
necessary-condition bound for low-weight words in the same class.

For each window, take the first 83 of its 84 outside packets as Q. This
omits packet 127 whenever it is outside the window, so the final unused
update is never charged. There are 166 quiet half-packets per region.
The remaining 90 halves contain 1440 unconstrained routing-pair slots.

Put y=3/4. In a quiet half, let f count full pairs and s count singleton
pairs. The fair swaps assign j singletons to one emission block and s-j
to the other. Define the positive polynomial

    P(u,v;y) = sum_(f+s<=16) C(16,f) C(16-f,s) u^f v^s 2^-s
               * sum_(j=0)^s C(s,j)
                   y^[min(4,f+j)+min(4,f+s-j)].

For a regional profile (f,s), its rank moment is exactly

    R(f,s) = [u^f v^s] P(u,v;y)^166 (1+u+v)^1440
               / [C(4096,f) C(4096-f,s)].

The orientation average is already in P; do not multiply the denominator
by another power of two. Consecutive supports have profiles (128,0) and
(127,2), in either order, each in 128 regions.

Choose alpha=log(4/3)/log(256)<1/16. The global source factor satisfies

    (eta_256 beta^256)^alpha
       <= 2*(4/3)^32768*(256/255)^512.

The factor 2 bounds eta_256^alpha. The other exponent relaxation uses
32*256*alpha<=512. Therefore the complete class probability is at most

    8192 * 85^256 * 2*(4/3)^32768*(256/255)^512
        * [R(128,0) R(127,2)]^128
      <= 2^-579.

This comparison is exact integer arithmetic. Truncation to u-degree 128
and v-degree two preserves both requested coefficients. Multiplying
each retained coefficient of P by 2^18 makes it integral; positive
integer convolution computes P^166 without rounding. The helper checks
both inequalities locating the resulting rational bound between 2^-580
and 2^-579. Its 223 exact tests and independent source audit passed;
the main replay took 0.45 seconds and generated no files.

Unlike the allowed-position cover, this bound permits nonzero inputs in
the quiet feedback kernels and allowed coordinates outside the windows.
It is restricted to exact consecutive source supports and the specified
state trajectories. It does not cover every 256-group set, every path,
or the entire low-distance event.

## Extend the quiet-window calculation

The same exact coefficient table covers consecutive supports of every
size from 7 through 256. A size 2m interval has matching profiles (m,0)
and (m-1,2); a size 2m+1 interval has profile (m,1) in both matchings.
Each size has 8192 intervals. Keep the 44-packet windows and the selected
83 quiet packets per region unchanged.

The source factor is at most 2*y^(-128r)*(256/255)^(2r), with y=3/4.
Exact integer comparisons, followed by upward dyadic rounding and an
exact union, give a bound of 2^-45 for sizes 196 through 256 jointly.
The preceding tail, sizes 195 through 256, gives only 36 bits by this
calculation. All 250 sizes were checked; no monotonicity was assumed.
The helper passed 270 coefficient-normalization checks and 278 direct
cyclic-profile checks. This remains an interval-support/window-path
result, not a certificate for all supports in that size band.

### A sharper scalar witness bound for quiet equations

Quiet equations need less machinery than the live-state potential.
Fix a source set T, a route, and a prescribed quiet-packet set Q before
sampling the outer. For each realized outer, let K be the E-defined
source kernel of the quiet equations, restricted to T. An F-valued
witness with exact support T implies that every group projection of K
is nonzero. For a uniform vector of K, each group is zero with
probability at most 1/256. If r=|T|<=256, there is a vector nonzero on
every group: otherwise every vector would miss at least one group, and
the zero vector would make the expected number missed strictly greater
than r/256. Multiplying a nonzero scalar witness by E minus zero gives
at least 255 witnesses. Hence the quiet-only event indicator is at most
n_T/255, and its conditional probability is at most

    min(1, beta^r * 256^(-D_Q) / 255).

This replaces eta_r for the quiet event, not for the live-state
potential. The earlier exact interval certificates remain valid with
their more conservative factors; their frozen helpers are unchanged.

For r>256, the same expectation argument gives a nonzero scalar witness
missing at most ceil(r/256)-1 groups. Thus an F-support band [r0,r1]
can be covered by scalar supports of sizes
[r0-ceil(r0/256)+1,r1], provided Q is the same prescribed set.
Union over these scalar supports once; do not also count their possible
T supersets. This scalar witness satisfies the quiet equations, not
necessarily the F-valued trajectory or its output-weight constraint.

### Arbitrary source supports: count both matching profiles

For a size-256 source set, let f0,f1 count full pairs in the two
matchings. Their singleton counts are 256-2f0 and 256-2f1, and the
number of cyclic source runs is c=256-f0-f1. There are 16640 compatible
ordered profiles. Their multiplicities cannot be replaced by two
independent matching counts.

For 0<f<L, the indexed cycle count is
N(L,f,c)=(L/c)*C(f-1,c-1)*C(L-f-1,c-1) for feasible positive c.
The earlier profile bound takes the minimum of N(8192,256,c) and the
two marginal counts C(4096,f_i)*C(4096-f_i,256-2f_i)*2^(256-2f_i).
If either matching has no singletons, contracting its pairs gives the
exact count N(4096,128,c). No quotient by rotations is taken.

An additional count bound deletes singleton pairs from one matching.
Let M=4096, f=f0, s=256-2f, and L=M-s. Deleting a 01 or 10 pair
reduces the shifted matching's 11-edge count by zero or one. The
remaining indexed 00/11 cycle has weight f and a run count c' in
[f-f1,f-f1+s]. If N(L,f,c') counts indexed binary cycles with these
parameters, the original profile count is at most

    C(M,s)*2^s * sum_(c' in [f-f1,f-f1+s]) N(L,f,c').

Clamp the sum to feasible run counts; constant cycles have one choice
and zero runs. Selecting the deleted positions, their orientations,
and the compressed cycle reconstructs the original set. Index the
compressed cycle from the least remaining original position, so no
extra rotation factor occurs. Take the minimum with the previous bound
and with the matching-swapped version. Direct enumeration of small
cycles checks every weight/profile, including the constant cases.

The first floating-point screen used the positive Cauchy coefficient
bound for P^166(1+u+v)^1440, multiplied by the maximal atom of the
1440-slot multinomial distribution. With y=3/4, the complete profile
union has log2 upper expression +511.841. The sharper quiet-only
source factor and a choice among complete bounds at y=3/4,2/3,3/5
reduce this to +119.297. The deletion counts give +118.728. They remove
122.166 bits from the former worst profile (125,51), but the worst
case moves to (115,54), with 87 source runs. None of these positive
expressions closes this restricted event; they do not refute the code.
The marks are common to both matchings within each complete bound.
Minimize complete bounds per profile, then sum the profiles.

### Retain randomness from the quiet halves

The preceding Cauchy bound discards the coefficient's probability mass.
More of that mass can be bounded without enumerating the full bivariate
coefficient table. Fix positive coefficient markers u,v and 0<y<1;
a marker may be zero only when its requested degree is zero. Write

    U=1+u+v,  B=1+y^2*u+y*v,  theta=B^16/P(u,v;y),
    p=y^2*U/B.

Coefficientwise, P dominates B^16, because each half's quiet rank is
at most 2f+s. Its normalized coefficient distribution is therefore a
mixture. With probability theta it consists of 16 independent B-slots;
otherwise it has the nonnegative residual distribution. In a B-slot,
the identity B=y^2*U+(1-y^2)+(y-y^2)*v extracts a U-slot with probability
p. All mixture choices precede the categorical draws.

Let M count the U-slots extracted from the 166 quiet halves. Its exact
probability generating function is

    E[t^M] = [1-theta+theta*(1-p+p*t)^16]^166.

The slots within one half share its first mixture coin; M is not a
binomial variable. Conditional on every mixture choice, the extracted
slots and the original 1440 free slots are independent of the residual
draws. All use the same probabilities (u,v,1)/U. Hence the coefficient's
normalized mass is at most E[a_M], where a_m is the largest atom of the
multinomial distribution with 1440+m slots and these probabilities.
Convolution shows that a_m is nonincreasing. For an integer cutoff m>=1,

    E[a_M] <= a_m + delta*(a_0-a_m),
    delta = min(1, t^(-(m-1))*E[t^M]),  0<t<1.

The exponent m-1 bounds the strict event M<m. Every marker, cutoff,
and t is a fixed numerical witness for its source profile, not a choice
made after observing the route.

The three-mark screen tests only cutoffs 0,800,1000,1200, where zero
denotes the old a_0 bound. It selects 1200,1000,800 for y=3/4,2/3,3/5,
respectively. The complete size-256 profile union has floating log2
bound -76.03854856; all 16640 terms are negative. The largest term is
about -85.0600 at profile (115,55), with 86 source runs. This improves
the bound without changing the construction or sampling distribution.
The mixture proof and source passed independent audits. Exact tests
cover 459 coefficient inequalities, nine mixture identities, 36 toy
compound-tail inequalities, and 32 multinomial modes.

Directed verification at 192 and 256 bits proves that the complete
restricted-class bound is below 2^-76. The verifier treats optimized
markers only as witnesses: it converts their exact binary ratios and
recomputes all 387 regional bounds with Arb interval arithmetic. It also
checks 774 multinomial modes by exact integer exchange inequalities.
For mode counts k and integer category weights w, the inequalities are
k_i*w_j <= (k_j+1)*w_i for every possible transfer i to j. These suffice
for a global mode because successive transfers from an original surplus
to an original deficit only decrease their probability ratios.

All 16640 complete profile bounds, their source multiplicities, and the
window union are included. The 256-bit replay encloses log2 of the bound
at -76.038548563327705159412054256907766001994286345490256137478612416063778,
with radius at most 5.22e-70. This certifies the arbitrary size-256
quiet-window class, not all size-256 trajectories or the whole code.

## Other refinements tested in the preceding round

The actual coefficient set has S=65280 elements, rather than 256. Its
extra entropy helps a dependent pivot but not an ordinary rank increase.
Let d0 be the initial quotient dimension and R the number of charged
fresh-coefficient pivots on a retained path. At most d0-1 pivots can
increase rank if a nonzero witness survives. Applying the projective
potential with weights 256 for independent pivots and S for surviving
dependent pivots gives

    Pr[path] <= P_256(d0) * 256^(-R)
                * (1/255)^max(0,R-d0+1).

The independent-pivot contraction and the dependent-pivot probability
use the full coefficient without exposing its b part separately.
For 1/255<=t<=1, the final factor is at most t^(R-d0+1). The existing
feedback matrix supplies the rank floor min(16-j,h) for incoming and
outgoing support sizes h,j. Marking this floor before pooling is valid.
At t=.99, however, the artificial score worsens from +98082.918905 to
+98299.703806 after the source-dimension correction. No marker sweep is
justified: this path family offers too few charged pivots.

Simply enlarging the initial potential does not evade scalar counting.
For tau=65280, kappa=log_256(tau), and alpha*kappa<=1, collapsing the
moment of P_tau(d) back to the scalar first moment reproduces the old
scalar exponent alpha*kappa. A genuine two-shadow calculation must retain
their shared outer randomness. For a random 4-by-4 Toeplitz map, two
independent input directions have joint image rank only 5--7, not 8.
Adjacent basis directions give rank 5. Treating these two images as
independent would be invalid.

The derivation is elementary: for input polynomials f,g of degree at
most three, divide out their gcd and put m=max(deg f0,deg g0). Relations
Af+Bg=0 with deg A,deg B<=3 have dimension 4-m. The joint map therefore
has rank 8-(4-m)=4+m. A future two-shadow bound must retain this geometry
and its source multiplicities, rather than multiply two outer densities.

### A short-return lemma retained for later use

There is also a necessary condition for a two-packet excursion from zero
to nonzero state and back to zero. Fix its allowed input sets I0,I1 and
its intermediate support H. Put G=im_E C_I0, c=dim G, and
V=im_E C_I1 intersect E^H, with v=dim V. Assume 1<=c<16 and v>0.
Here H16 is the fixed superregular matrix. D1 and D2 are independent
diagonal matrices with independent uniform entries in E minus zero,
and D2=diag(b). The implication below holds for every a in E^16.
Let B_G be a basis matrix for G and let pi_V be a quotient map with
kernel V. A possible excursion requires

    rank_E(pi_V D2 H16 D1 B_G) < c.

Indeed, if this leading matrix has full column rank, some c-by-c minor
of pi_V D2(theta I+diag(a))H16 D1 B_G is a nonzero polynomial in theta
of degree c<16. It cannot vanish in the specified extension field.
The possible birth state therefore cannot belong to V extended to F.

For each nonzero fixed u in G and each nonzero target w in E^16,
Pr[D2 H16 D1 u=w]<=255^-16. Superregularity charges the zero positions
through D1; the independent D2 entries charge the nonzero target values.
Union over projective E-directions in G and nonzero elements of V gives

    Pr[the necessary rank condition]
       <= min(1, P_256(c)*(256^v-1)/255^16).

If c=0 or v=0, the excursion is impossible. The lemma is not integrated
into a path calculation. It exposes the b=D2 coefficients, so it cannot
be combined on the same refresh with an argument that keeps all 65280
coefficient values fresh. Longer excursions, rather than this short case,
remain the current priority.

## Joining quiet constraints to live paths

The window calculations do not partition all trajectories. A global
argument must also charge repeated activation and prolonged live state.
Let R=N/64=32768, and write b_t for the indicator that the incoming
state of output packet t is nonzero. There are R-1 regular updates,
and b_0=0. Put a=sum_t b_t, let k count 0-to-1 transitions, and set
epsilon=b_(R-1). The number q0 of regular 0-to-0 transitions is

    q0 = R-1-a-k+epsilon.

For a>0, the number of binary paths with a active states and k runs is
C(a-1,k-1)*C(R-a,k). Fixing epsilon=1 replaces the second factor by
C(R-1-a,k-1); epsilon=0 gives C(R-1-a,k). Treat the all-quiet path
separately. These counts concern the binary path, not source supports
or field-valued state assignments.

A large quiet count alone cannot prove distance. For a sparse scalar
outer word supported on s groups, at most 256s packets contain input.
One can declare those packets and neighboring transitions live, leaving
almost all other steps quiet, without requiring the declared live path
to be attainable. The quiet equations then hold automatically at empty
input packets. Dropping activation and return constraints therefore
admits deterministic witnesses, not just a slightly loose count.

For chosen thresholds k_star,q_star, a complete next test should partition
paths into three classes:

1. k<=k_star and q0>=q_star: use the quiet equations with the run count.
2. k>k_star: retain and charge actual activation/return transitions.
3. k<=k_star and q0<q_star: use persistent-live output bounds.

Raw markers for these counts must enter each atomic transfer before
pooling intermediate masks and taking fractional powers. For example,
multiply 0-to-1 entries by u and 0-to-0 entries by v; the unused
terminal update receives neither marker. If a class implies
u^k*v^q0>=c, its fractional-moment bound has correction c^(-alpha).
For the second class take u>1,v=1,c=u^(k_star+1). For the third,
take u=1,0<v<1,c=v^(q_star-1). Do not assign these markers to a pooled
tile after its internal transitions have been forgotten.
Alternatively, a quiet-only operator built directly from already
fractional rank moments may insert count marks afterward. Its correction
is then c^(-1), not c^(-alpha). Keep these two conventions separate.

This is a proof plan, not a computed cover. The three bounds must be
closed together for the same thresholds and actual source profiles.
Further restricted-window certificates alone do not close the global
7--7583-group gap.

### An all-binary-path quiet pilot

The next pilot removes the 44-packet-window restriction, but retains the
branch k<=K and q0>=Q. It tests four actual size-256 source profiles:
(128,127), (115,55), (64,64), and (0,0). The two coordinates count full
pairs in the alternating matchings. For matching i, put s_i=256-2f_i.
Fix positive route markers u_i,w_i, and write U_i=1+u_i+w_i. Use the
quiet-half polynomial P_i from above, evaluated at y=256^-alpha. Define

    S_i = [[v*P_i^2/U_i^32, b], [1,1]],
    H_i = 4096*ln(U_i)-f_i*ln(u_i)-s_i*ln(w_i)
          -ln(C(4096,f_i)*C(4096-f_i,s_i)).

A zero marker is allowed only for a zero target degree. Each matrix
entry bounds a fixed incoming and outgoing binary state. The complete
chain starts in state zero and sums over both terminal states:

    exp(128*(H_0+H_1))
      * e_0^T (S_0^128 S_1^128)^127 S_0^128 S_1^127 (1,1)^T.

Multiply by the source count, (beta^256/255)^alpha, and b^-K*v^-Q.
Here 0<b<=1 and v>=1 are post-fractional count marks. The terminal packet
has neither a state update nor a count mark. Matching order is retained.

The earlier independent-slot extraction also gives a uniform coefficient
mass bound for every binary path. A quiet half has extraction PGF J(t);
a free half has PGF t^16<=J(t). There are 256 halves per region. Hence
Pr[M<m]<=t^(-(m-1))*J(t)^256, and the coefficient mass is at most
a_m+delta*(1-a_m). Here a_m is the maximum multinomial atom with m slots;
there are no 1440 guaranteed free slots in this calculation. The bound
is uniform over paths, so it also applies to their positive matrix sum.

At Q=24576, the four-profile floating union has log2 bound -3141.784
for K=256 and -510.227 for K=512. At K=256,Q=20480 it is +3085.348.
The source and derivation passed independent reviews; exact toy tests
cover 146 path terms and 355 branch-mark inequalities. These are selected
floating witnesses, not a directed certificate or an all-profile cover.
Further profile expansion must wait until the complementary branches
work at compatible thresholds.

### Fixed paths test the complementary bound

A separate diagnostic extracts positive terms from the selected scalar
majorant. A positive term can reject that majorant at fixed z,alpha and
all admissible count marks. It is not a codeword, a lower bound on code
failure, or a rejection of every z,alpha choice. The diagnostic averages
raw-cap fractional powers over the exact oriented route placements.
It uses the actual consecutive-support profiles (128,0) and (127,2),
4096 supports per orientation, and the true output-only terminal.

Two tested paths lie outside the previously certified window class.
One has two separated excursions in each region, 512 births in total,
and 24319 quiet regular steps. The other straddles regional boundaries:
it has 256 births, 24320 quiet regular steps, and 8192 active incoming
states. Its first region starts at zero; later regions start and end
with full state support. All but the first region fail the nonwrapping
44-packet-window restriction. Five exact path checks verify these counts.

At Q=24576, their log2 contributions are:

| (z,alpha) | Two separated excursions | Boundary-straddling excursions |
|---|---:|---:|
| (.85,.523) | +65247.85 | +80227.72 |
| (.7,.25) | +21179.96 | +27146.15 |
| (.5,.125) | +14163.40 | +16331.04 |
| (.99,.05) | -9162.14 | -6938.10 |
| (.9,.05) | -7582.65 | -5476.13 |

Thus the failure at the first three anchors is parameter-dependent.
The negative entries do not close the complementary branch: they test
individual paths, not the sum over intermediate supports and paths.
The main replay passed 36 geometry checks, 256 direct slot assignments,
37 rational convolution checks, six terminal checks, 18 count-mark
inequalities, and the five out-of-window path checks.

The low-alpha anchors do not evade every fixed path. Consider a path
that activates at the first packet, has full state support thereafter,
and ends with the output-only terminal. Every selected full-state
continuation and terminal cap is at least z^64. The activation cap has
the same lower bound whenever the first packet contains input. These
are lower bounds on the chosen caps, not on transition probabilities.
They follow by retaining the zero-requested-output term; each gamma
split, joint-live minimum, and singleton refinement preserves them.

For the (128,0) first matching, its first paired slot is full with
probability 128/4096=1/32. One orientation of 4096 source intervals
therefore contributes at least

    128 * (256^32768 * z^1363149)^alpha.

Here 1363149=N-floor(.35*N). The expression drops factors greater than
one from eta and beta. This path has one birth and no quiet steps, so
the complementary branch's admissible count marks cannot reduce it.
At z>=.9 and alpha=1/20, the exact inequality .9^6>1/2 makes its log2
contribution greater than 1754.625. Thus those high-z anchors still fail
for this untruncated majorant. The actual all-live path is not a bad
codeword: its output floor already exceeds the target. This obstruction
does not apply after imposing that global floor.

For a useful additional constraint, let h_t count the nonzero coordinates
of the actual incoming state, and put S=sum_t ceil(h_t/4), including the
terminal packet. At least ceil(h_t/4) of the four RS state blocks are
nonzero. Each contributes at least 13 output records before adding X_t.
Thus wt(Y_t)>=13*ceil(h_t/4)-wt(X_t). An outer support of r groups has
total input weight at most 256r, so every output of weight at most L has

    S <= floor((L+256r)/13) = 61503  when r=256.

A raw marker w^ceil(h_t/4), with 0<w<=1, imposes this necessary condition
with correction w^(-alpha*61503). Mark the terminal emission as well.
These masks refer to the actual state before scalar relaxation, not the
support of an auxiliary E-valued witness. No independence from the output
weight event is assumed: both inequalities hold on each retained bad path.
The extra marker can exclude long full-state paths without maintaining
an explicit counter of length 61504. It is distinct from the earlier
lower bound on outer input weight. The next summed-path pilot includes
this marker; w=1 recovers the untruncated expression.

### Summing the complementary state paths

A bounded follow-up uses all 17 state-support sizes and the actual two
consecutive-support profiles. It places the birth, quiet, and global-floor
marks into raw packet matrices. It then pools the internal masks of each
16-packet tile before taking the fractional power. The terminal tile
contains 15 regular updates and one output-only emission.

This pilot also applies Jensen to route placements at fixed coefficient
markers. If Mbar is the normalized raw packet matrix, its tile entry is
D_j^(1-alpha)*(Mbar^16)_[h,j]^alpha, where D_j=C(16,j). The regional
Cauchy normalization stays outside alpha. Local singleton degrees above
two are omitted, because neither target coefficient can use them; the
full route mass remains in Jensen. No source class is replaced by an
artificial repeated marginal.

The four anchors z in {.8,.85}, alpha in {.05,.1} do not close either
complementary branch. The best selected log2 expressions are +37827.018
for q0<=24575 and +39277.686 for k>=513, both at z=.85,alpha=.05.
The selected global-floor mark is one. This does not show that the floor
constraint is useless: the bounded single-start search does not establish
a minimum over all witnesses. Independent source review passed, including
the terminal floor mark. Exact tests include 96 route normalizers and
2500 deterministic output-floor comparisons.

Do not continue a broad parameter or tile-size scan of this relaxation.
The next controlled comparison should retain full-pair and singleton
degrees inside each raw tile. Apply Jensen separately for each degree,
restoring the route multiplicity before extracting the actual regional
coefficient. Compare one anchor against the evaluated-route version.
This isolates the loss from averaging different occupancies together;
the positive pilot alone does not quantify that loss or refute the code.

### Conditioning on tile degrees: controlled comparison

The next comparison retains each raw tile's full-pair count f and
singleton count s. Fix the same consecutive source supports, z=.85,
alpha=.05, quiet mark .0044813232524843144, and floor mark one. The
branch remains q0<=24575. No encoder or setup distribution changes.

Let C_p[f,s,h,j] be the raw p-packet matrix coefficient, including
intermediate mask sums and destination multiplicity D_j=C(16,j). Put
V_p(f,s)=C(32p,f)C(32p-f,s). Degree-wise Jensen gives the tile coefficient

    V_p(f,s)^(1-alpha) D_j^(1-alpha) C_p[f,s,h,j]^alpha.

For the terminal tile, use the coefficient of P^(p-1)T and omit D_j.
After evaluating the coefficient markers, divide by U^(32p). Keep the
regional Cauchy correction outside alpha, as before. Truncating f>128
or s>2 inside a tile is safe for these two target regional profiles.
This also deletes impossible degrees that the earlier evaluated-route
relaxation retained; the comparison includes both improvements.

At the previous route witnesses, the log2 expression falls from
+37827.01773005103 to +36376.29710010885. Optimizing only the three
route markers gives +36327.86231064020; the evaluated-route expression
at those same markers is +37787.75257500522. Thus the observed gain is
about 1451 bits, not the speculative tens-of-thousands estimate. The
dominant path sum does not realize that estimate. The remaining positive
gap is too large to justify more route-marker tuning alone.

The main replay took 29.24 seconds. Independent source review passed,
as did 108 exact coefficient entries, four alpha-one invariances, six
Jensen comparisons, and 96 exact route normalizers. These floating
expressions are diagnostics, not directed certificates or evidence of
a low-distance codeword.

### Longer pooling and an exact regional control

At the same parameters and route markers, a second diagnostic pools
16,64,128 raw packets. It reports the full 17-state expression and two
restrictions: only tile endpoints in {0,16}, or every packet state in
{0,16}. Both restrictions are positive parts of the chosen majorant;
neither bounds the full source class.

| Packets per tile | All 17 states | Endpoints only 0/16 | All packets only 0/16 |
|---|---:|---:|---:|
| 16 | +36376.30 | +10629.79 | +9966.86 |
| 64 | +19350.27 | +12951.83 | +12649.21 |
| 128 | +16552.60 | +13361.62 | +13187.08 |

These are floating log2 expressions. Pooling changes three relaxations:
mask counting at boundaries, conditional route Jensen, and deletion of
impossible total degrees. A 16-packet tile permits total degree at most
128, but evaluating eight such tiles can still count regional degrees
above 128. A 128-packet tile deletes them before evaluation. Therefore
the improvement is not a measurement of mask entropy alone. In particular,
a proposed upper limit on the improvement based only on deleted mask
boundaries was withdrawn. The restricted scores need not improve with
tile length, because the route Jensen loss also changes.

The audited run took 56.15 seconds and peaked below 95 MiB. No output
files were generated. The longer pool is useful, but the positive
zero/full-state control prevents closure by intermediate-mask refinements
alone at this fixed witness.

The region-length raw polynomial also permits exact coefficient
extraction. For profile (f,s), let V=C(4096,f)C(4096-f,s). Replace Cauchy
evaluation by the entry

    D_j^(1-alpha) (C_128[f,s,h,j]/V)^alpha.

Use (E_128[f,s,h]/V)^alpha for the terminal, where E_128 is the coefficient
of P^127 T. No route markers or regional H factors remain. This still
uses conditional Jensen within the specified regional degree; it is not
an exact average of fractional powers over individual routes.

For the all-packet 0/full16 control, D_j=1. Exact extraction reduces
+13187.079472798836 to +11617.188069039783, a gain of 1569.8914 bits.
The main replay took 2.78 seconds and passed the pinned control, 12
coefficient-versus-Cauchy comparisons, and the prior exact tests.
Independent source review passed. The positive result obstructs this
selected majorant, not the code or every choice of marks and exponents.

The original persistent-live branch also requires at most 512 births.
These controls dropped that restriction. They retained the global floor
only with mark one. For zero/full states, let a count active incoming
packets, including the terminal. The intended branch satisfies

    q0=32767-a-k+epsilon <=24575,  k<=512,  4a<=61503,
    hence 7680<=a<=15375.

Before changing the encoder or launching another parameter grid, measure
the tilted birth, quiet, and active counts in the exact regional control.
They determine whether the omitted birth restriction or unused floor can
improve this fixed bound. A positive dual coefficient calculation can
obtain these derivatives without finite differences or a sweep.

### Branch sensitivity: the omitted birth cap is not the immediate loss

The positive dual calculation now measures those counts. Each coefficient
carries its mass and its nonnegative count-weighted mass. Products add
counts; sums average them with positive weights. Apply alpha to the log
derivative when taking the regional fractional power. The finite chain
then gives a hierarchical tilted measure, not a distribution of actual
codewords or observed encoder behavior.

The main replay obtains

    E[k]=148.0937357934,
    E[q0]=25579.9768024540,
    E[a]=7039.9238323517,
    E[epsilon]=0.9943705947.

Their identity E[q0+k+a-epsilon]=32767 agrees within 5e-9. It checks
the terminal count independently of the bound reproduction. The replay
also passes 432 exact-rational dual entries and four exact finite-chain
derivatives. Independent source review passed; runtime was 2.61 seconds.

For a raw mark exp(x) on a count C with upper limit c, the derivative
of the corrected log2 expression is alpha*(E[C]-c)/ln(2). At the fixed
witness, the derivatives for births, quiet steps, and the floor count
4a are respectively -26.25029, +72.49375, and -2405.21101. The log bound
is convex in these log marks. Since the birth and floor marks are at
their upper endpoint x=0, moving either alone to x<0 cannot improve it.
A stronger quiet penalty does improve it locally.

The next bounded test should solve only the quiet-mark stationarity
condition E[q0]=24575, keeping z,alpha and the raw caps fixed. Recheck
the birth and floor derivatives there; their signs here do not exclude
joint effects after changing the quiet mark. This diagnoses the remaining
Chernoff mismatch without a parameter grid. No new support class or
whole-code certificate follows from the sensitivity calculation.

### Quiet calibration does not repair the regional pool

The bounded stationarity check keeps z=.85, alpha=.05 and all raw caps
fixed. It changes only the raw quiet mark. The optimum is approximately
.0023312810999, with tilted quiet count 24575.00000004. The log2 control
falls from +11617.18806904 to +11592.38616256: only 24.80 bits.
The birth and floor slopes there remain negative, approximately -26.25
and -2115.22. Joint convexity in the three log marks gives the KKT
diagnosis: these marks cannot repair this control at the selected anchor.
These are floating calculations, not directed numerical certificates.

### Isolate the loss on one admissible path

The boundary-straddling path above has 256 births, 24320 quiet updates,
and 8192 active incoming packets. It satisfies the complementary budgets
k<=512, q0<=24575, and 4a<=61503. This is an admissible mask path; no
claim is made that some input realizes it.

Fix that path, the two actual consecutive-support profiles, z=.85 and
alpha=.05. Compare two orders of operations, using identical raw caps:

1. Average products of powered local caps over regional routes.
2. Average products of raw caps over those routes, then apply alpha.

There is no intermediate-mask sum in this comparison. Exact regional
coefficient extraction supplies the same route normalization in both
cases; support multiplicities and the source factor are unchanged.
The log2 expressions are respectively -4598.42402459 and +10679.55735626.
Thus route Jensen alone loses 15277.98138 bits on this selected term.

The second term remains positive after any of the three allowed raw
marks u,v,w<=1. Its correction multiplies it by

    u^(-256 alpha) v^(-255 alpha) w^(-28735 alpha) >= 1.

Consequently, neither those marks nor exact enforcement of their three
budgets can fix this particular route-pooled majorant at this anchor.
This does not reject the construction, other anchors, or other proofs.

### Retain routing while summing binary state paths

Form a two-state polynomial matrix P=[[Q,A],[E,L]]. Each entry sums its
local occupancy cap raised to alpha, weighted by the exact route-placement
count. Use terminal column (T,U). For profile (f,s), extract

    [x^f y^s] P^128 / V(f,s),
    [x^f y^s] P^127 (T,U)^T / V(f,s),
    V(f,s)=C(4096,f)C(4096-f,s).

Apply no further fractional power. The alternating finite regional chain
now sums every path whose actual masks are zero or full16. This is a
restricted part of the route-first operator, not an upper bound on all
source words. Its ordering against the earlier route-pooled operator is
not automatic: routing Jensen and path concavity are different losses.

With neutral marks, the log2 expression is +537.70899537. The all-quiet
path alone contributes -4480.15361210. Thus an already-covered quiet path
is not responsible for the positive result. The previously optimized
regional-pool quiet mark gives +8459.29800625 here; it is inappropriate
for this changed operator.

Positive dual convolution gives tilted counts k=873.53194, q0=15999.67744,
a=15894.71663, and terminal activity .92601266. Their path-count identity
agrees within 2e-8. Unlike the earlier pool, the new dominant terms exceed
the birth and output-floor budgets. The corrected slopes for birth, quiet,
and floor marks are +26.07902, -618.57877, and +149.74212. Restoring birth
or floor penalties can improve this expression locally; no optimized or
negative branch bound is claimed.

A single preselected comparison at z=.9, alpha=.05 gives +3546.74229505,
despite improving the one fixed path to -5476.12789464. Improving that
path therefore does not imply improving the sum. No further anchor scan
was performed.

### Pool live masks without averaging away quiet suppression

The fixed-path comparison suggests a weighted rather than uniform route
average. For each route, write the raw path product as B*V. Let B contain
the quiet, birth, and return factors; V contains sustained-live factors
and the full-state output-only terminal. For 0<alpha<1, weighted Jensen
gives

    E[(B V)^alpha] <= E[B^alpha]^(1-alpha) E[B^alpha V]^alpha.

Both expectations use the same regional route distribution. If E[B^alpha]
is zero, both sides are zero; otherwise normalize B^alpha to a probability
density and apply concavity. The highly variable B factors therefore
retain their fractional power before route averaging. The two coefficient
normalizers combine into one V(f,s), not its square or its alpha power.

For the selected boundary-straddling path, this bound is -4598.42402459,
within 3e-11 bits of the route-first value. It recovers the full diagnosed
15278-bit loss to floating precision. Six exact-rational weighted-Jensen
tests, 96 free-packet route counts, the terminal checks, and independent
source review passed. This near-equality concerns one diagnostic path,
not all occupancies, paths, or source supports.

The next test should retain a fixed binary quiet/live schedule and fixed
birth/return endpoint masks. Sum partial nonzero masks only inside each
sustained-live segment before applying alpha. In that test the Q/A/E
factors are common, so the weighted inequality remains applicable. For
the first comparison, also fix full16 masks at regional cuts and at the
terminal. The base expectation E[B^alpha] must use free route counts on
live packets, with no hidden-mask multiplicities. Only the mixed
expectation E[B^alpha V] contains those masks and their raw matrix products.
Normalize after combining all segments of a region, not separately per
segment. For a full binary-schedule aggregation with unspecified endpoint masks, only
the Q factors are automatically common: A and E depend on those masks.
Either retain the endpoints as indices or move their factors into V.
Do not silently reuse the fixed-path factorization after pooling them.

This is the next proof architecture to test, without changing the encoder.
No additional whole-source class closed in this round. The nonquiet
7--7583-group cover remains open.

## Seventeenth round: return supports and a fixed live cone

The first test follows the sixteenth-round plan. Keep the straddling
binary schedule, full16 birth/return endpoints, regional cuts, and final
incoming mask fixed. Sum all nonzero masks only inside sustained-live
segments. The weighted route calculation changes from -4598.42402459
to -4544.22732062, a cost of only 54.20 bits. The selected bound remains
negative. This is a restricted schedule with fixed endpoint masks,
not a complete source-class certificate. Independent source review and
432 explicit route/mask checks passed; the run took 23.42 seconds and
peaked at 95,035,392 bytes.

The next operator removes the fixed binary schedule and endpoint masks.
A single positive vector bounds every nonzero actual-mask row throughout
the chain. Birth enters this cone, a live step preserves it, and return
and terminal use the corresponding row sums. This yields a two-state
operator before alpha, without averaging regional routes away. Its exact
regional coefficients cover all actual masks and binary schedules for
the consecutive256-group source class. They do not cover every source
set of size256. The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) states
the actual-mask invariant and the order of averaging.

At z=.85, alpha=.05, neutral path marks, the controlled results are:

| All-mask bound for consecutive256-group sources | Log2 expression |
|---|---:|
| Fixed live cone with old return caps | +4312.54152383 |
| Enforce disjoint incoming supports at return | +3736.19672317 |
| Also retain nonzero-syndrome RS output floors | +3262.21324769 |
| Also retain the first live packet's support restriction | +3133.11957238 |

The old zero/full16-only reference remains +537.70899537. It covers
fewer paths and is not an alternative upper bound on this all-mask class.
The table compares proof refinements only; the encoder is unchanged.

The return improvement uses u=Cx before summing incoming supports H.
An explicit factor q/(q-|H|) permits the auxiliary scalar support to be
exactly H. Those supports are then disjoint, instead of overlapping
across a cardinality maximum. With the fixed cone, the complete factor
is kappa q^-16, where kappa=16777216/15671025. The output moment is that
of (Id+AC)x with Cx nonzero. The positive block-floor recurrence further
retains the RS distance outside the allowed input positions. Exact GF4
checks and independent proof/source reviews passed.

The one-step restriction uses six states: zero, four birth types, and
mature. A birth touching b input blocks activates at least 17-4b state
coordinates; a following packet touching b' blocks can cancel at most
4b'. Thus b+b'<=4 excludes an immediate return. This removes 172 local
return records in the test, but improves the score by only 129.09 bits.
Immediate returns alone are not the full remaining loss.

The two-state floor operator's tilted counts are 1665.46 births,
15063.80 quiet updates, and 16038.60 active incoming packets. They are
counts under the proof majorant, not encoder measurements. The typical
contribution has roughly 9.6 live packets per birth. Positive dual
convolution verifies the 32767-update identity and the honest terminal.

### Recalibrate alpha after changing the operator

The alpha=.05 witness was inherited from the older route-pooled proof.
For the new two-state floor operator, its score derivative is
-172170.31 bits per unit alpha. A single predeclared comparison at
alpha=.10 lowers the score to +137.11622265 and changes the derivative
to +16628.65. Thus the convex one-dimensional minimum lies between
these two anchors. This is a substantial parameter mismatch, not evidence
that another construction change is required.

The derivative calculation seeds each local powered mass with
-log(raw cap), averages these costs within each exact route coefficient,
and propagates them through the finite chain. The derivative is the
source rate minus the expected total raw cost in bits. All caps in this
diagnostic are at most one; the helper rejects larger caps rather than
misusing its nonnegative dual arithmetic. Source multiplicity has zero
derivative. Independent review, exact rational cost tests, and scalar/dual
agreement passed. These are floating diagnostics, not directed bounds.

The bounded solve replays both endpoints and performs six bisections.
Its best evaluated witness is alpha=119/1280=.09296875, with two-state
score +78.57605118. The final derivative bracket is [.09296875,.09375].
At that same witness, the newborn refinement gives +71.54683757. The
small 7.03-bit difference confirms that extra first-step memory is not
essential once the exponent is calibrated.

A final local check rebuilds the complete local caps at three nearby
output marks, keeping alpha=119/1280 and every path mark equal to one:

| Output mark z | Two-state floor-cone log2 expression |
|---|---:|
| .84 | -218.66953880 |
| .85 | +78.57605118 |
| .86 | +613.63398816 |

The newborn operator at the selected z=.84 gives -227.35829705. Both
negative expressions include all fine state masks, all binary schedules,
all 8192 consecutive source sets of size256, and the actual alternating
regional profiles. There is no window restriction. They concern arbitrary
F-valued messages with exact source support in that class and the
implemented limited refresh distribution, not uniform full-field refresh.

This is the first negative all-trajectory floating bound for this source
class in the present proof search. It is not a directed certificate and
does not close arbitrary size256 source supports or the entire 7--7583
gap. The simpler two-state witness already has about 219 bits of numerical
room; verifying it is preferable to first adding the six-state machinery.
The final local comparison took 9.75 seconds. The preceding root check
took 27.91 seconds. Both passed independent source review and old-score
controls, and generated no files.

## Eighteenth round: directed verification and a broader support class

The two-state witness survives a substantial simplification of its local
caps. Use exact quiet and activation entries, the positive return-floor
bound, and the original output-only bound for both live continuation and
terminal. Dropping the remaining feedback refinements costs only
0.00001178 bits in the selected floating comparison. Independent review
checked that each simplified entry bounds the same actual-mask sum.

A new verifier implements those formulas directly, without importing the
old17-state numerical kernels. It uses exact rational geometry and Arb
intervals throughout the fractional powers and polynomial products.
Both192- and256-bit replays put the consecutive-support log2 bound at
-218.669527023371, proving the explicit200-bit check at35% distance.
The256-bit run uses exact rational endpoint tests: testing containment of
a rounded rational ball was unnecessarily stronger than testing the
rational point itself. The adapter changes only this toy test, not the
frozen verification arithmetic.

Next, count actual source patterns by their cyclic runs. If the matching
full-pair counts are(128-a,128-b), the support has a+b runs and the two
singleton counts are(2a,2b). A positive boundary-parity recurrence counts
each ordered profile exactly, without identifying rotations. Its small
tests exhaust5440 indexed masks and748 profiles. The
[rank interface](CHEAP_LANE_RANK_INTERFACE.md) gives the counting formula
and the finite chain for each ordered profile.

The directed union bounds for at most2,8,16 runs are respectively
-208.301408590818,-166.933170386442,-133.814902217474 in log2 units.
The last includes all152 ordered profiles and every actual-state path.
Its192-/256-bit replays agree; no interval intersection with the
nonnegative half-line was needed. The largest contribution is the
balanced16-run profile, but this observation is not a monotonicity proof.

The next calculation removes the run restriction. Carry-safe polynomial
packing evaluates four coefficient rectangles, covering all129 marginal
matching profiles. A last-exit-tree count on a four-vertex graph replaces
the expensive boundary-parity histogram. Exhaustive tests cover5440 masks
and1370 all-weight profiles;279 cofactor checks and27 crosschecks with
the earlier counter also pass.

The384-bit full replay gives a directed log2 union of
-93.5574055068572906550312346546009957059259358406723846743,
with radius below3.3e-110. All16640 ordered matching profiles and all
state trajectories are included. Counts agree with every cyclic-run total
and sum to C(8192,256). The consecutive-source control agrees with its
earlier192-/256-bit replays. The explicit40-bit test passes; no positivity
intersection is needed. This is now a certificate for every size256 support
at35%, not merely selected or bounded-run patterns.

The largest run class is c=43, contributing log2 bound-96.981220741206.
Its intermediate fragmentation explains why the consecutive control alone
was insufficient. A conservative per-run maximum gives total-93.511494866948,
close to the exact profile-weighted count here. This numerical agreement
does not imply a general logarithmic bound on the conservative count's loss.

A256-bit follow-up retains all16641 coefficients with2f+s<=256 and
checks nine source weights using their own counts and comparison factor.
It reproduces the size256 interval and closes all size255 supports with
log2 bound-105.267775674949285. All48152 nonempty ordered chains pass the
exact count checks, including odd source weights. The same fixed witness
at sizes7,8,16,32,64,128,192 gives positive bounds; these are failures of
this witness, not low-distance examples. The remaining nonquiet source
weights are7--254 and257--7583.

Next retune the output mark toward1 for smaller sizes. Above256, use the
rank-sensitive initial-kernel split or the existing bounded-missing-group
comparison before investing in a new two-shadow operator. These are proof
refinements, not encoder changes. The
[rank interface](CHEAP_LANE_RANK_INTERFACE.md) states their scope and the
unresolved exceptional class; no source-size monotonicity is assumed.

## Size-dependent witnesses

Retuning proof parameters, without changing the encoder, closes all
source supports of sizes7--32 jointly at35% distance with at least178 bits.
The fixed schedule is z=1-r/1280 and alpha=1/2 for r<=15,else3/16.
Every size uses its own exact indexed counts; all26 directed bounds are
summed as exact conservative dyadic probabilities. No source-size
interpolation is used. The full run takes63.78seconds with stdout only.

At r=128,z=91/100,alpha=119/1280, the directed log2 union is
-196.40889813805544 over all4224 ordered profiles and all trajectories.
The largest run class has39runs. Positive bounds from earlier parameter
choices were proof slack, not low-weight codewords.

At the end of round19 the nonquiet gaps were33--127,129--254,and258--7583. The new
actual-support padding and paired-cell-envelope interfaces in the
[rank note](CHEAP_LANE_RANK_INTERFACE.md) are independently audited;
their usefulness beyond the existing coverage still needs computation.
The first h=1 check at257 groups gives directed log2 bound+0.895648,
not a certificate. Its exact count and arithmetic checks pass. The
missing-group penalty therefore remains a concrete proof bottleneck;
it is not evidence of a low-distance codeword.

Retaining the missing group's nonzero kernel projection until its outer
average removes that penalty. The new source factor is
beta^r*(r+1)/[(q-1)*(1-r/(2q))]. Applying its exact ratio to the preceding
completed257-group sum gives log2 bound-94.28703297309922, verified by
directed256-/384-bit algebraic adjustment. The coefficient sum is reused,
not freshly recomputed. A separately audited joint missing-set argument
now supplies a source comparison through all relevant source sizes;
those larger regional calculations were still pending. See the rank
interface for the conditioning order and exact joint factor.

The new shared-coefficient cell verifier now gives contiguous35% coverage
through224 active groups, with more than43 bits jointly. It computes every marginal in each stated band,
then sums exact indexed run counts with paired matrix envelopes. It does
not interpolate between endpoints. The joint missing-group factor extends dense coverage to7552--8192
with more than781 bits jointly. See the round20--21 results in the
[goal ledger](CHEAP_LANE_MIXING_GOAL.md) for fixed rational witnesses.

For larger supports, positive coefficient bounds avoid full regional
tables. Initial selected-profile witnesses are loose, but an exact
finite-chain derivative shows their natural markers are badly calibrated.
At the selected512-group profile, calibrating only those markers gives a
directed107.5447-bit bound at35%, replayed at256 and384 arithmetic bits.
The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) now gives a way to reuse
markers across profile intervals while covering all indexed supports.
The initial complete FLOAT size512 screen is positive, dominated by clustered
supports near66 runs. Calibrating the actual equal profile(223,66) leaves
log2 bound+2760.155. Retaining both occupancy residues modulo32 gives
+311.432, and a bounded output-mark adjustment gives+304.353. Convexity
rules out alpha-only closure of this fixed residue witness. Modulus64 then
gives+44.0843, but exact coefficient extraction gives-82.3388. A separate
384-bit replay at different positive route scales confirms this82-bit
selected-class certificate. A subsequent residue128 FLOAT calculation sums
all66048 ordered profiles and gives log2-59.9352. The directed hybrid
calculation then closes the complete source size with log2-60.8663747948,
including all matching profiles and every state trajectory. Exact high-count
coefficients and mixed exact/cyclic low-count tables remove positive aliases;
the encoder is unchanged. Neighboring source sizes remain to be covered.

## Reproduction and next step

Run the following helpers with `python -B`; each has a bounded runtime,
a 20-GiB free-space preflight, and stdout-only results:

- `theory/cheap_lane_path_decomposition.py`: restricted regional operators.
- `theory/cheap_lane_surplus_pivots.py`: two fixed surplus-pivot marks.
- `theory/cheap_lane_window_diagnostic.py`: real route and artificial term.
- `theory/cheap_lane_window_cover.py`: exact 592-bit routing-event bound.
- `theory/cheap_lane_quiet_window_cover.py`: exact 579-bit quiet-window class.
- `theory/cheap_lane_quiet_window_band.py`: exact 45-bit consecutive-size tail.
- `theory/cheap_lane_quiet_profile_screen.py`: first arbitrary-profile screen.
- `theory/cheap_lane_quiet_profile_marks.py`: three complete rank-mark bounds.
- `theory/cheap_lane_quiet_profile_counts.py`: singleton-deletion counts.
- `theory/cheap_lane_quiet_profile_smoothing.py`: compound-slot screen.
- `theory/cheap_lane_quiet_profile_verify.py --precision 256`: directed 76-bit class bound.
- `theory/cheap_lane_quiet_path_pilot.py`: selected all-binary-path quiet bounds.
- `theory/cheap_lane_path_marks.py --run`: fixed-path majorant diagnostics.
- `theory/cheap_lane_marked_state_pilot.py`: complementary 17-state tile pilot.
- `theory/cheap_lane_degree_conditioned_tile.py`: fixed tile-degree A/B comparison.
- `theory/cheap_lane_pooling_scale_pilot.py`: fixed 16/64/128-packet controls.
- `theory/cheap_lane_active_blocks.py`: actual (weight,active-block-count) classes.
- `theory/cheap_lane_exact_region_control.py`: exact regional coefficient control.
- `theory/cheap_lane_branch_sensitivity.py`: positive dual count derivatives.
- `theory/cheap_lane_quiet_stationarity.py`: one-dimensional quiet calibration.
- `theory/cheap_lane_fixed_path_jensen.py`: isolate route Jensen on one path.
- `theory/cheap_lane_route_first_control.py`: exact-route binary path sum.
- `theory/cheap_lane_route_first_sensitivity.py`: its counts and pure-quiet term.
- `theory/cheap_lane_route_first_anchor.py`: one preselected anchor comparison.
- `theory/cheap_lane_quiet_tilt_jensen.py`: weighted route Jensen on the fixed path.
- `theory/cheap_lane_live_segment_tilt.py`: partial masks inside fixed live segments.
- `theory/cheap_lane_live_cone.py`: one fixed cone for all live masks.
- `theory/cheap_lane_disjoint_return.py`: exact incoming supports at returns.
- `theory/cheap_lane_return_floor.py`: nonzero-syndrome output floors and path counts.
- `theory/cheap_lane_newborn_cone.py`: one-step birth-support memory.
- `theory/cheap_lane_cone_alpha.py`: alpha derivative and one predeclared comparison.
- `theory/cheap_lane_cone_stationarity.py`: bounded alpha calibration and same-alpha A/B.
- `theory/cheap_lane_cone_output_anchor.py`: three nearby output marks at the selected alpha.
- `theory/cheap_lane_cone_simple.py`: controlled floating simplification.
- `theory/cheap_lane_cone_interval_tests.py --precision 256`: directed consecutive-support replay.
- `theory/cheap_lane_run_profile_counts.py --max-runs 6`: exact indexed source counts.
- `theory/cheap_lane_cone_run_verify.py --max-runs 16 --precision 256`: directed bounded-run union.
- `theory/cheap_lane_cone_packed.py --precision 192`: carry-safe coefficient products and controls.
- `theory/cheap_lane_fast_profile_counts.py`: exact all-weight source counts using last-exit trees.
- `theory/cheap_lane_cone_profile_cover.py --precision 384`: every size256 support, all trajectories.
- `theory/cheap_lane_cone_weight_pilot.py --precision 256`: exact nine-weight pilot and size255 certificate.
- `theory/cheap_lane_cone_weight_witness.py --weight 128 --z 91/100`: size-dependent full-support witness.
- `theory/cheap_lane_cone_small_band.py --precision 256`: joint directed7--32-group certificate.
- `theory/cheap_lane_missing_group_receipt.py --precision 384`: audited algebraic adjustment closing257 groups.
- `theory/cheap_lane_cone_cells.py --rlo 169 --rhi 184 --z 331/375 --alpha 119/1280`: joint directed support-band bound.
- `theory/cheap_lane_dense_joint.py --precision 384`: joint dense endpoint7552--8192.
- `theory/cheap_lane_cone_cauchy_anchors.py`: four selected-profile coefficient-bound diagnostics.
- `theory/cheap_lane_cauchy_marker_gradient.py`: finite-chain derivative at the512-group anchor.
- `theory/cheap_lane_cauchy_marker_calibration.py --precision 384`: directed selected-profile recalibration, not a whole-size cover.
- `theory/cheap_lane_cauchy_profile_screen.py`: FLOAT complete size512 screen with exact source counts.
- `theory/cheap_lane_cauchy_profile_calibration.py --full 223`: directed clustered-profile marker calibration.
- `theory/cheap_lane_cauchy_modular.py --full 223 --u 29411/311125 --v 10441/460035 --include-16`: one-axis regional degree residues.
- `theory/cheap_lane_cauchy_bivariate_modular.py`: two-axis regional occupancy residues; see the goal ledger for fixed z/alpha diagnostics.
- `theory/cheap_lane_cauchy_fine.py --mode exact --precision 384 --u 3/32 --v 3/128`: exact clustered-profile coefficient with independent precision/scales.
- `theory/cheap_lane_cauchy_hybrid_directed.py --precision 256`: complete size512 certificate, all66048 ordered profiles and state paths.
- `theory/cheap_lane_cauchy_neighbors_directed.py --rlo 508 --rhi 516 --seconds 1200 --precision 256`: shared-table directed band; the508--515 subtotal has43.5679 bits, whereas including516 gives37.0975 bits.
- `theory/cheap_lane_profile_peeling.py --weights 512`: FLOAT diagnostic of hypothetical on-demand counting savings, not a faster directed verifier.
- `theory/cheap_lane_broad_peeling_screen.py --weight 512`: implemented on-demand counting control, with a complete residual envelope;3.3-second FLOAT replay, not a certificate.
- `theory/cheap_lane_broad_peeling_screen.py --weight 1024 --z 304397809/625000000 --ordinary-fallback`: first broad FLOAT diagnostic; explicitly retains unresolved boundary coefficients by a loose Cauchy bound.
- `theory/cheap_lane_broad_profile_calibration.py`: fixed-count marker and residue controls for two selected size1024 profiles; not a whole-size bound.
- `theory/cheap_lane_singleton_endpoint_directed.py`: exact directed all-singleton size1024 profile, more than4318 bits at35%; no claim about other profiles.
- `theory/cheap_lane_return_decomposition.py`: return-marked diagnostic for the same two selected profiles; marked values are not complete code bounds.
- `theory/cheap_lane_exact_small_return.py`: exact one-/two-input return moments and a controlled FLOAT replay;734/797-bit improvements, still positive.
- `theory/cheap_lane_exact_triple_return.py`: exact three-input return moments, including an independent projective enumeration; selected FLOAT replay only.
- `theory/cheap_lane_binary_pair_pooling.py`: pool the binary intermediate state, retaining total pair degrees; Jensen alone loses too much.
- `theory/cheap_lane_pair_concavity.py`: occupancy-level separable correction to the pair cap, with nine fixed rho values.
- `theory/cheap_lane_pair_concavity_triple.py`: compose exact D1/D2/D3 and pair corrections; selected FLOAT log2 bounds-1068.237/-1431.422, not a whole-size certificate.
- `theory/cheap_lane_pair_bank_screen.py --centers 32,96,160,224,288,352,416,448,464,472,480,488,496,504,508`: complete r1024 FLOAT screen; seven unresolved boundary marginals dominate its explicit fallback bound.
- `theory/cheap_lane_old_boundary_regions.py`: directed old-cone regional bounds for full degrees0..6, with singleton residue512 and an exact singleton control; source/count factors are omitted here.
- `theory/cheap_lane_pair_bank_boundary_screen.py`: complete mixed-region r1024 FLOAT sum, log2-688.3404 over all263168 profiles with no fallbacks; not a directed certificate.
- `theory/cheap_lane_pair_residue_control.py --modulus 256`: same complete FLOAT witness, log2-688.2132; modulus128 loses too much to aliasing.
- `theory/cheap_lane_pair_local_directed.py`: coherent exact rational upper pair coefficients, including only lower corrections and the Jensen minimum.
- `theory/cheap_lane_directed_run_mass.py`: exact tests for the complete indexed-run envelope and both terminal orientations.
- `theory/cheap_lane_integer_regions.py`: exact Fraction tests for integer convolution, positive folding/truncation, and upward dyadic rounding.
- `theory/cheap_lane_pair_regions_integer.py --seconds 300`: directed complete size1024 certificate, more than430 bits at35%, all263168 ordered profiles;154.5 seconds and no generated files.
- `theory/cheap_lane_integer_precision_control.py`: fixed128-/192-bit control of the high rectangle and two ordered profiles; explains258 bits of rounding loss but does not improve the whole-size certificate.
- `theory/cheap_lane_pair_sizes_screen.py --weight 6144`: complete FLOAT screen with audited hole-degree reindexing; retains the actual source size and both terminal orientations.
- `theory/cheap_lane_coset_sizes_screen.py --weight 2048`: matched complete screen with the inside-coordinate coset return bound; no new certificate.
- `theory/cheap_lane_pair_profile_calibration.py`: fixed size2048 profile(768,768), comparing inherited, own-profile, and fitted markers at residues256/512.
- `theory/cheap_lane_coset_profile_control.py`: same selected comparison with the coset return bound.
- `theory/cheap_lane_pair_parameter_control.py --z 23/100 --alpha 13/128`: selected coset-cap parameter control, FLOAT log2+1353.3423.
- `theory/cheap_lane_residual_rank_return.py --tests-only`: exact checks for the tensor-row rank lemma and positive residual counts.
- `theory/cheap_lane_residual_rank_control.py`: matched selected comparison with that rank bound; use `--complete --weight 2048` for a complete FLOAT screen.
- `theory/cheap_lane_rank_parameter_control.py --z 23/100 --alpha 13/128`: combined selected refinement, FLOAT log2+419.2413; still positive.
- `theory/cheap_lane_rho_bucket_control.py`: same-parameter/marker comparison of nine global choices,33 global choices, and33 first-occupancy buckets; FLOAT+419.24,+48.13,−337.41 respectively.
- `theory/cheap_lane_rho_bucket_sizes.py --fitted-center768 --modulus 512 --seconds 300`: complete size2048 FLOAT screen, all1050624 profiles; positive total and resolved-profile envelopes, explicit fallback retained.
- `theory/cheap_lane_rho_profile_calibration.py --full 898`: bounded actual-profile calibration and coherent rational parameter controls; not a whole-size certificate.
- `theory/cheap_lane_rho_calibrated_bank.py`: complete31-center size2048 FLOAT screen; zero fallbacks, log2+686.4416.
- `theory/cheap_lane_calibrated_bank_parameters.py --z 23/100 --alpha 3/32`: same complete coverage at an explicit proof witness, log2+243.5942.
- `theory/cheap_lane_rho_fine_buckets.py --full 915`: matched coarse/fine rational-bucket control at identical markers, FLOAT+239.0055 to+46.7391.
- `theory/cheap_lane_rho_fine_bank.py`: complete fine-grid size2048 FLOAT bound,+51.9694; all profiles, no fallbacks.
- `theory/cheap_lane_rho_grid_bank.py --subdivisions 8`: retain prior caps and refine the rational grid; complete FLOAT+13.6371. The16-subdivision option is unrun.
- `theory/cheap_lane_rho_ordered_buckets.py --full 920`: exact off-diagonal ratio minima, matched selected FLOAT+47.3975 to-34.5311.
- `theory/cheap_lane_rho_ordered_bank.py`: complete ordered-ratio FLOAT bound-29.9703 over all1050624 size2048 profiles, zero fallbacks; not a directed certificate.
- `theory/cheap_lane_rho_ordered_bank.py --alpha 23/256 --modulus 256`: selected complete witness for directed follow-up, FLOAT-86.7425; modulus512 gives-86.7463 at identical markers.
- `theory/cheap_lane_rho_weighted_ordered.py`: unrun optional prototype for a fixed5/4 ratio-separation deficit; not used in the86.7-bit witness or any certificate.
- `theory/cheap_lane_ordered_local_directed.py`: exact coherent powered-table freezing, integer ordered corrections, and independent Jensen caps for all six pair channels; local tests and three independent reviews pass.
- `theory/cheap_lane_ordered_regions_directed.py`: selected center893192-/256-bit control; `--center 1024` checks the high-full-degree endpoint without old boundary substitutions.
- `theory/cheap_lane_ordered_regions_directed.py --complete --seconds 1200`: complete directed size2048 certificate,86.7425 bits at35%, all1050624 profiles;848.8 seconds, no generated files.
- `theory/cheap_lane_band_regions.py`: one cyclic table with separate actual-weight coefficient extraction and exact normalizers; selected2048 endpoints reproduce the frozen engine exactly.
- `theory/cheap_lane_band_run_mass.py`: source-scaled band envelopes, exact indexed counts, and complete original-marginal checks; singleton and small exact controls pass.
- `theory/cheap_lane_ordered_band_directed.py --complete --seconds 1200`: directed joint certificate for sizes2047--2049,74.2104 bits at35%, including all1050625 ordered envelope pairs;856.3 seconds, zero generated data.
- `theory/cheap_lane_counted_band_run_mass.py`: degree-specific exact-count roots before envelope maxima, with own-r partner validity; exact tests and two independent reviews pass.
- `theory/cheap_lane_counted_band_directed.py`: selected2040--2056 A/B, directed log2+11.1082 to-65.2459 at(895,895). Its `--complete --seconds 1800` replay certifies the entire band jointly at35% with60.2491 bits, all1058840 ordered envelope pairs, in981.4 seconds.
- `theory/cheap_lane_ordered_large_screen.py`: actual size4096, current rank/ordered pair operator and exact all-singleton boundary; complete FLOAT log2-4621.1552 at modulus256 and-4621.1766 at512, all4198400 ordered profiles, no fallbacks. Not a directed certificate.
- `theory/cheap_lane_ordered_large_directed.py`: directed selected profile(2048,0) at actual size4096, log2-4621.6241621; complete mode now closes all4198400 profiles at35% with4621.1552 bits in941.6 seconds.
- `theory/cheap_lane_ordered_ordinary_directed.py`: same local caps and31 markers without coefficient projection, exactf0 retained; complete directed log2+262.5035, not a useful guarantee.
- `theory/cheap_lane_ordered_hole_screen.py`: actual6144, current operator, complete once-only hole reversal; all1050624 profiles, FLOAT log2-53312.9361, no fallback. Not a certificate.
- `theory/cheap_lane_ordinary_hole_directed.py`: frozen6144 markers, ordinaryH64/H63J, own hole denominators and actual source factor; complete directed35% bound49748.3343 bits in32.3 seconds.
- `theory/cheap_lane_ordinary_profile_screen.py`: complete FLOAT diagnosis plus selected directed checks; edge0 identifies(2013,36), not f1, as the coarse ordinary bound's leader. Optional endpoint refinement is not used for that diagnosis.
- `theory/cheap_lane_edge_regions_directed.py`: exact near-singleton FULL/EMPTY coefficient rectangle, with380 rational checks passing; reusable partial regional extraction, not a source-size certificate.
- `theory/cheap_lane_ordinary_marker_envelope.py`: all31 existing ordinary witnesses, own denominator recurrence, six independent entry selections; complete directed35% bound256.2334 bits at4096 in85.6 seconds. Preserve the stronger4621.1552-bit coefficient receipt.
- `theory/cheap_lane_ordinary_hole_band.py`: exact rational fixed-cell geometry and directed log-Gamma denominator bounds; complete6144--6160 joint bound42710.6165 bits. Wider fixed-cell trials fail because marker assignment does not track relative occupancy.
- `theory/cheap_lane_normalized_hole_band.py`: same marker bank and local operator, with cells indexed by2048a/k; complete joint35% bounds25386.7993 bits for5888--6400 and4134.3211 bits for5632--6528. All integer sizes, including odd sizes, are covered. The5120--6528 trial is positive and adds no coverage.
- `theory/cheap_lane_hole_marker_fit.py`: bounded FLOAT proposals of31 positive rational marker pairs, using a read-only copy of the already-built local table; not a certificate.
- `theory/cheap_lane_calibrated_hole_band.py`: independent directed replay of those witnesses with unchanged local coefficients and own-size source counts. Five new disjoint bands join the earlier bounds to cover4864--8192 jointly at35% with more than780 bits. The4352--4863 trial is positive and adds no coverage; exact parameters and hashes are in the goal ledger.
- `theory/cheap_lane_calibrated_actual_band.py`: audited unreversed-coordinate extension below4096, with the actual-r source tilt. The3072 controls remain positive; no coverage is inferred.
- `theory/cheap_lane_cell_run_counts.py` and `theory/cheap_lane_calibrated_run_band.py`: exact run-interval masses upper-bound cell populations, retaining both terminal orders and each actual size's source factor. Matched directed replay closes4352--4383 jointly at35% with221.7125 bits. The wider4352--4479 and singleton3072 checks improve but remain positive.
- `theory/cheap_lane_route_free_screen.py --bands 4301:6143,6144:7551 --z 3/1250`: historical FLOAT diagnostic; the stronger analytic obstruction now rules out this relaxation at every size512--7525.

The proof and source reviews were independent of the implementation;
the main agent also replayed the selected checks. These runs changed no
encoder. Rounds37--38 retain two selected regional-bound witnesses totaling
27539 bytes; earlier runs generated no data files. Together with the earlier small-support
certificates, the directed bounds now cover every size through751,
1024,2040--2056,3008--3392,4096,4352--4383,and4544--8192, at the40-bit working target. Round39 improves the old37.1139-bit
size516 bound to1205.5529 bits inside its new band replay. Joint whole-code coverage is still missing.
Extend to the remaining source sizes. Source multiplicity and singleton routing compete,
so no profile monotonicity is assumed.
Round24 rules out dropping route probabilities in most of that gap.
At size1024, exact extraction repairs the all-singleton numerical fallback,
but selected clustered profiles remain positive after marker fitting and
residue512 controls. Return suppression identifies paths with returns as
the remaining source of loss. The actual return-matrix refinement recovers
734/797 bits but does not close those profiles. Round25 adds exact D3
and pools a binary intermediate state without the costly occupancy-Jensen
relaxation. Both selected profiles now have negative FLOAT bounds with
large slack. The complete size1024 screen, with directed boundary regional
bounds replacing seven noisy FFT entries, gives log2-688.3404. It includes
every profile and residual count, with no unresolved fallbacks. Round26
reconstructs a complete directed bound using exact integer polynomial
products and upward coefficient rounding. It closes size1024 with430.0100
bits. A fixed128-/192-bit comparison isolates258 bits of leading-profile
loss from upward rounding of the tiny f512,s0 return coefficient. This
selected-profile control does not replace the complete430-bit bound.
Round27 tests actual sizes2048,4096,6144 and separates marker/projection
loss from local return-cap loss. A new affine-coset lemma and the actual
tensor-row rank constraint tighten return moments without encoder changes.
At the selected size2048 profile(768,768), calibrated markers and these
local bounds reduce the FLOAT log2 bound to+419.2413. Retaining first-packet
occupancy in the pair correction lowers the same comparison to−337.4096.
This remains FLOAT evidence and does not cover other profiles in the run.
The complete size2048 screen reproduces the anchor but remains positive:
the full envelope is+18329.88 with364 fallback marginals, and the separate
resolved-profile envelope is+546.18, led by(898,898). Round28 establishes
that the selected positive bound survives own-profile calibration. A fixed
calibrated31-center bank removes every fallback. At alpha3/32, finer pair
corrections and exact off-diagonal ratio minima reduce the complete
size2048 FLOAT bound to-29.9703. One midpoint proof-parameter control
improves this to-86.7463, with a matched smaller-grid replay at-86.7425.
Round29 reconstructs these fixed witnesses using exact integer products
and directed rounding. It closes the entire size2048 class at35% with
86.7424999 bits. The two selected precision controls and complete replay
use identical exact local coefficients. Round30 extracts sizes2047--2049
from shared tables and certifies the joint band with74.2104 bits at35%.
The new2049 endpoint is included; exact support counts and singleton
reproduction pass. Larger source-size gaps remain open. The count-aware
envelope recovers76.3541 bits in a selected17-size control, without
changing local caps or the encoder. Round31 verifies that entire band
jointly with60.2491 bits at35%. It then tests the current operator at
size4096: two complete FLOAT replays agree within0.022 bits, near
log2-4621.16, with no fallback coefficients. A directed selected-profile
control reproduces the leading(2048,0) bound and verifies the singleton
boundary. Round32 completes directed size4096 coverage at4621.1552 bits.
It closes actual size6144 with an ordinary-Cauchy bound of49748.3343 bits,
retaining the actual source factor after the hole-coordinate transformation.
The ordinary shortcut at4096 initially remains positive; a complete
profile diagnosis locates its worst loss at coarse marker-cell edges,
not at the all-singleton endpoint. Entrywise selection among the same31
witnesses repairs it, giving256.2334 bits with a much cheaper directed
verification. Preserve the existing construction; support-size bands,
not isolated checkpoints,
are the next coverage objective.
See the goal ledger for the remaining gaps. No encoder work was added.
Round35 adds the continuous4544--8192 cover and the4352--4383 band.
The new run-aware count reduces the3072 log2 bound from+7548.07 to+2015.56,
but does not close that size. Next retain run-count/conditional-cost
correlation inside the leading cells before changing local caps.
Round36 completes that run-interval refinement over all496 cells. The
size3072 log2 bound becomes+1546.85 and remains positive. A matched
directed selected-profile comparison identifies coefficient-bound slack:
at full counts(1186,1452),434 runs, residue256 projection changes the
full-run-charged bound from+1468.91 to-3379.95. Both terminal orders pass.
The same rational local table and markers are used throughout. This is
not a complete source-size certificate by itself. The subsequent complete
regional replay closes every3072-group support at35%, with3356.7491 bits,
including all2362368 ordered profiles and all state histories. The encoder
is unchanged. Extend shared-table coverage to neighboring sizes next;
do not infer those sizes from the singleton's margin.
Round37 verifies universal maxima over the positive cyclic256 tables,
separately for each regular and terminal entry. The same31-marker bank
then closes3008--3136 with197.1555 bits and3137--3168 with957.9652 bits.
These disjoint bands cover all161 sizes jointly with more than196 bits
at35%, using actual-size denominators, source factors, and indexed counts.
The compact selected witness has186 exact dyadics and a pinned canonical
hash. Cached replay reproduces the first band's full directed interval in
18.4 seconds without rebuilding regional polynomials; it remains conditional
on the verified generation of that witness. The wider2944--3007 and
3137--3264 trials are positive. Their outer endpoints remain positive
individually, so band splitting alone does not repair this fixed bank.
These are unresolved upper bounds, not low-distance counterexamples.
Recalibrate near those endpoints next. The encoder and whole-code status
remain unchanged; the goal ledger records exact intervals and provenance.
Do not cap individual transfer entries at one: they are not complete
conditional event probabilities.

Round38 recalibrates the proof witnesses at3264,z=1/12 and verifies all31
regional tables. This closes3169--3328 jointly with1308.6306 bits.
The new13,785-byte selected witness reproduces that complete interval
in19.6 seconds; its metadata and marker fingerprint receive an independent
audit. Cached replay also closes3329--3392 with2383.1518 bits. The four
disjoint bands now cover all385 sizes3008--3392 jointly with more than196
bits at35%. Generation, not a matching hash alone, establishes the saved
maxima. No encoder change is involved.

All31-by31 saved marker pairs at2944 leave its completed bound unchanged
at log2+1146.3154. Only minor cells improve. The3393--3456 band remains
positive at+635.5442. These failed bounds do not establish low distance.
Source review also identifies641 as a historical verifier limit, not a
hypothesis of the local proof. A separately tested lower-domain adapter
is the recommended next step; it has not yet been implemented. Preserve
the actual-size source factor, including eta_r<1 below255. Whole-code
coverage remains open; exact commands and pins are in the goal ledger.

Round39 checks the lower-domain extension, including odd source sizes
and the source-factor boundaries. Ordinary225 and516 bounds remain
positive. A new bounded helper instead uses the current ordered pair
caps, exact downward coefficient extraction, and upward integer rounding.
It covers225--232 jointly with334.5700 bits and every516-group support
with357.8913 bits, at35% distance. Retuned old witnesses close229--254;
separately complete contributions from the positive257--289 trial close
258--279 jointly with more than55 bits. The resulting1--279 cover has
more than43 bits jointly. No result is inferred for280--289 from that trial.

The new516 result is not a pure exact-versus-approximate extraction
comparison: its local ordered-pair operator is also stronger than the
historical one-packet cone. The source factor, route denominators, threshold,
and update chronology agree, and the new run-count envelope is conservative.
The encoder and its work remain unchanged. Whole-code coverage is still
the task; extend bounded neighboring bands without extrapolating their
endpoints.

The first516--548 extension has a positive aggregate bound at192 coefficient
bits. Its odd/even split and scale-dependent516 result point to large
relative upward-rounding inflation on rare coefficients. A matched256-bit
replay keeps all local caps, rectangle markers, and source factors unchanged
and closes the entire band with842.2922 bits jointly. It takes136.6 seconds,
peaks at264601600 bytes, and writes no data files. The508--548 cover now
exceeds43 bits jointly. This refinement changes proof arithmetic, not the
encoder or its setup distribution.

The next280--312 band closes jointly with198.4058 bits at z677/825 and
alpha119/1280, using256-bit coefficient rounding. The small-source union
now covers every size1--312 with more than43 bits jointly. This is complete
per-size verification, not interpolation. Next test313--345 at the proposed
z5449/6765 and keep bridging actual sizes; the larger whole-code gaps remain
open.

Round40 closes all195 sizes313--507 with the same frozen exact helper.
All runs retain alpha119/1280,256-bit coefficient rounding,384-bit Arb
arithmetic, actual-size source factors, indexed counts, and both terminal
orders. A first346--378 witness gives only26.8648 bits. Retuning z alone
closes the same33-size band with440.5422 bits; it does not trim the band or
change the code. Subsequent directed bands close379--507. No result is
inferred by interpolating bounds: interpolated witness proposals are each
checked on every included size.

The next549--581 band closes with86.9979 bits at the older z17447/25000
witness. This adds228 sizes overall. The resulting1--581 cover exceeds42
bits jointly; nested1--507 exceeds43.
Including all other currently certified classes and the true-quiet
bad-distance event still gives a union bound below218*2^-50<2^-42.
This accounting excludes the weak26.8648-bit trial and does not double-count
nested receipts. Six support-size gaps remain. The encoder, actual limited-
support refresh distribution, and certified baseline are unchanged.
Retune the proof witness for582--614 and continue exact verification before the
more expensive middle gaps. Commands and directed intervals are in the
[goal ledger](CHEAP_LANE_MIXING_GOAL.md#fortieth-round-continue-the-exact-cover-and-retune-proof-witnesses).

Round41 closes582--614 and615--647 with1151.0231 and1172.7312 bits jointly.
An independently reviewed wrapper reuses each exact rectangle for up to129
actual sizes, while retaining fresh target coefficients, normalizers,
source factors, and complete explicit-weight joins. Exact controls pass.
The full648--776 trial completes but has log2 bound+1978.4852. Its104
completed contributions648--751 are each below2^-76, verified by exact
rational comparison of all reported interval endpoints. Thus this prefix
alone has more than69 bits jointly. No later source size is admitted from
the failed aggregate. Continuous coverage now reaches751, with the full
covered-class union still below2^-42. Six gaps prevent whole-code closure.

The failed tail shows large odd/even jumps. The next matched control keeps
the local table,z,alpha,source factor,and run mass unchanged at r776,c55.
Both matching full counts then lie in333--388. A smaller exact singleton
downset retains all56 ordered profiles and both terminal orders while
potentially reducing upward-rounding inflation. This is proof arithmetic,
not additional encoding work; the control cannot certify all size776 paths.

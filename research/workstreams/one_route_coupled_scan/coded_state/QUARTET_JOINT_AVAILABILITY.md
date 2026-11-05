# Joint global availability and unavailable coverage

This checkpoint keeps the construction and full 45%/40-bit objective.
It develops the joint count specified in
[the availability census](QUARTET_AVAILABILITY_CENSUS.md), without
returning to local seven-type refinements. The full certificate is open.
All bounds retain the once-sampled accepted setup and free terminal.

## Retain a subset of unavailable covered positions

For a canonical history, write E_t=V_t\A_t. Every quartet meeting A_(t-1)
must avoid E_(t-1), and every quartet meeting E_(t-1) must be nonfull in
V_t. Keeping only the size of A dropped this joint restriction.

Choose fixed retained-count levels J, including zero. After reconstructing
E_t=Gamma_t(A_(t-1))\A_t, choose uniformly a subset of size
floor_J(|E_t|). Carry that subset to the next step and enforce its required
nonfull-quartet conditions. Canonical histories pass every retained
constraint. Auxiliary choices have total weight one, so no binomial
penalty is added. These are proof variables; encoder setup is unchanged.

The count state is (u,e), with u available positions and e retained
unavailable covered positions, where e is in J and u+e<=512. The sets
are disjoint. States with ceil(u/4)+ceil(e/4)>128 cannot have a separating
partition; their entire transition and terminal rows are zero. They may
be removed without dropping a complete history. The basis is fixed
independently of the routing quota tilt.

This rounding must retain an actual subset. The barycentric upward
rounding used for the earlier survival tags does not automatically apply
to this different constraint.

## Positive transition and terminal

Put n=512, M=128, a_birth=mp chi, r0=r chi, and

    c_(k,u)=[x^u]((1+x)^4-1)^k,
    Q(y)=1+4y+6y^2+4y^3,
    G_(s,e)(y)=[z^e]((1+z)^4 Q(y)+y^4)^s.

The last polynomial counts retained old holes in nonhit quartets: a
quartet with old holes cannot become a full fresh birth. If k quartets
are hit by old availability, their contribution is

    binom(M,k)c_(k,u)r0^(4k) / multinom(n;u,e,n-u-e)
      * (alpha w+b y)^u(w+b y)^(4k-u)G_(M-k,e)(a_birth y).

Extract y^v and sum w-degrees f whose retained level is floor_J(f)=g.
This is K_((u,e),(v,g)). The resulting subset is uniform conditional on
its size, so the next independent partition sees only these count states.
Fix complete candidate masks before averaging partitions; do not condition
a fresh partition on the realized availability set that it helps define.

The final update enforces ACTUAL new E empty, not merely retained g=0:

    t_(u,e)=sum_k binom(M,k)c_(k,u)(r0 b)^(4k)
       G_(M-k,e)(a_birth) / multinom(n;u,e,n-u-e).

Thus the complete count is e_(0,0) K^4095 t. This includes every A_T and
requires E_T=empty. Alpha is a positive closure weight, not a probability
of retaining the new E subset, so class markers may make alpha exceed one.

The implementation uses positive polynomial Horner recurrences. Exact
new-hole counts below max(J) are kept as polynomial rows; the final row
sums the tail at or above max(J). Multiplication by w+b y shifts exact
counts and grows the saturated tail positively. The old-availability
factor alpha w+b y is then applied with positive binomial convolutions.
No large nearly equal quantities are subtracted.

## Empty regions and exact quotas

The zero-coverage operator is Z=v e_(0,0)^tr with v_(u,e)=1[u=0], not
a projector supported only at (0,0). Old retained holes may disappear
before an empty region. Z^2=Z. For P=K-Z, positive doubling uses

    Q_(2r)=Q_r^2 + v(e_(0,0)^tr Q_r) + (Q_r v)e_(0,0)^tr,
    Q_r=K^r-Z.

The last-region positive vector is K^15 t-v, computed from t-v by the
recurrence d<-K d+P v. The occupied-region generating calculation uses
A_j=e_(0,0)^tr Q^j v and B_j=e_(0,0)^tr Q^j d. This preserves empty
regions entered from a retained-hole-only state.

Condition the eligible-label counter to exactly h positions in each
8192-position region. The actual conditioned operator has a common basis
for every Bernoulli tilt p. Therefore entrywise minima of quota-adjusted
REGIONAL positive operators are valid. Do not minimize one-step operators
before extracting the regional quota.

`quartet_joint_subset.py` passes 748 exact matrix comparisons, 60 terminal
checks, ten impossible-state checks, 59 explicit uniform-subset rows,
nine composed free-terminal comparisons, six general-zero checks and
18 occupied-region checks. Independent audits also cover normalization,
Horner saturation, conditioning, and endpoints.

## Two alternative positive comparisons

The exact two-count kernel can act on an exponential potential
phi_(u,e)=t^u s^e without constructing its dense matrix. Evaluating the
local polynomial at new A=t and new E=s leaves A(x)+B(z), where

    A_i=r0^4 binom(4,i)(s+b t)^(4-i)(alpha s+b t)^i, i=1..4,
    B_0=(1+a_birth t)^4,
    B_j=binom(4,j)Q(a_birth t), j=1..4.

The coefficient array of (A+B)^128 is a sum of 129 positive outer
products of univariate polynomial powers. Divide by the old multinomial
to obtain K phi. All 131841 count states can be checked without a dense
physical transition. Sixty explicit small-instance action checks pass.

A scalar supremum K phi/phi is valid but loose in the selected probes:
its worst row is full availability. A four-sector comparison distinguishes
(u,e)=(0,0), u=0<e, e=0<u, and u,e>0. Disjointly supported potentials
make KH<=HB a four-state upper comparison with an exact corresponding Z.
Its positive presence recurrences pass 240 action checks, 240 row
inequalities, 128 regional entries and eight occupied-tail checks.
The selected t=s=1 result remains positive (+936801 FLOAT log2 bound).
Neither potential comparison is preferred to the existing tagged control.

## Current numerical evidence

For h=2048, theta=1/16 and log2 z=-7/16, three quota tilts give:

| Retained E levels | Count states | FLOAT log2 upper bound |
|---|---:|---:|
| {0,1} | 1022 | +723667 |
| {0,16} | 1010 | +504681 |
| {0,64} | 962 | +642928 |

The last choice is not a stronger rounding of {0,16}: it discards every
actual unavailable set of size below 64. Its all-positive-region marginal
is concentrated near 450 available positions and zero retained holes.
This is a marginal of the relaxed weighted count, not a distribution of
actual code trajectories.

`quartet_joint_bins.py` adds a positive comparison with disjoint count
bins. Within each retained-E level, keep u=0 separate and bin positive u
in intervals. For a fixed potential phi=t^u and disjoint lift H, compute
each entry of B_plus as the maximum of (K-Z)H/phi over the source bin.
Then (K-Z)H<=HB_plus and ZH=HZ_bin exactly. The terminal-positive vector
is bounded the same way. Regional composition and quota minima remain
valid in this common basis. Source-bin maxima may choose different u at
different steps; this is an upper comparison, not exact lumping.

Exact toy checks pass 2542 row inequalities, 276 endpoint inequalities,
eight singleton-bin identities, 24 zero intertwinings and 24 regional
comparisons. At the same three-tilt point, width-four u bins give:

| Retained E levels | Compressed states | FLOAT log2 upper bound |
|---|---:|---:|
| {0,1,16} | 382 | +506781 |
| {0,1,16,64} | 495 | +379257 |

The second run takes about 204 seconds on one Peach thread, within its
350-second guard. A rank-cap marker of log2 weight -1/4 on the first
grid gives +486486. These controls remain weaker than the earlier
tagged control; none closes a row class.

These are complete-history upper-bound diagnostics, not certified margins
or sampled-code statistics. They show a material benefit from retaining
more unavailable coverage, but do not close the remainder. The previous
tagged control near +100811 used a different relaxation and nine tilts;
do not interpret the table as a same-parameter performance comparison.

All experiments are serial, guarded, single-threaded and file-free on
Peach. No production, paper, construction, benchmark, commit or push
changes accompany this checkpoint.

## A small norm control and its limitation

There is also a 15-state local norm bound for the exact two-count
interface. For n coordinates put M_n(u,e)=multinom(n;u,e,n-u-e),
S_n=diag(sqrt(M_n)) and D_n=diag(t^u s^e). The orthonormal count matrix
is S_n D_n^-1 K_n D_n S_n^-1. Averaging the uniform quartet partition
is an orthogonal symmetrization commuting with D_n. Consequently its
norm is at most sigma^128, where sigma is the norm of the corresponding
four-coordinate operator. Its 81-configuration norm equals the norm
on the 15 symmetric count states: the nonnegative Gram matrix has a
nonnegative maximizing eigenvector, which can be symmetrized without
vanishing. The free-terminal factor is (1+t^-2)^256.

This comparison preserves the exact availability interface but drops
the requirement of 129 occupied regions. The helper passes 4860 exact
physical/count entries, 12 norm comparisons, three zero-operator norms
and six terminal checks. A 144-point marked diagnostic remains above
1.1 million log2 units at its best tested point. It is not competitive.

There is a useful warning independent of that search. The local full-A
diagonal is (rb)^4 and is unchanged by the potential. The exact-h
binomial term implies r^8192/pi >= a^(-(8192-h)). At the selected
unmarked point this alone forces the complete norm bound above
19658.6875, even before the outer-row-set union. Yet full coverage makes
all 8192 rows eligible, not 2048. Thus this is an ineligible history
admitted by the relaxation, not a low-distance codeword. Coverage and
rank constraints must be retained before interpreting such a norm.

## Couple unavailable coverage to outer zero rows

The accepted map supplies an additional necessary condition. Let H be
the actual set of outer rows with covered degree at least 129. Every
row outside H is zero by RS shortening. In a quartet Q hit by old
availability, put r=|A_prev intersect Q|. Then

    |E_new intersect H^c intersect Q| <= r-1.          (1)

If the left side is zero this is immediate. Otherwise E_new is
nonempty in Q. A zero actual incoming state would give V_new=S_new
there, and S_new is contained in A_new, a contradiction. Let k>=1 be
the actual incoming support size, so k<=r. On an outside-H position,
X_new=0; an E_new position also has Q_new=0. It is therefore a zero
coordinate of the image of that k-supported incoming state. The
superregular quartet map permits at most k-1 such zeros. This proves
(1), without an independence assumption or an additional probability
penalty. Two independent audits agree.

Let E_out count all such outside-H unavailable positions. There are
(F-f)/4 old-hit quartets in total. Summing (1), with the free endpoint
retained, gives

    F-f+4E_out <= 4L-4|A_T| <= 4L.                   (2)

The weaker form fits the existing positive transfer. For eta>=1 insert
the single marker eta^(4L-F+f-4E_out). Put

    r_A=p+(1-p)/a,       r_E=p+(1-p)/(a eta^4),
    chi'=chi/eta,        m'=eta^5 m,
    b'=eta^4 b r_A/r_E,  r'=r_E.

Keep alpha unchanged. There is no global compensator. If an outside
rank-cap marker is also used, multiply both occurrences of (1-p)/a
by its weight. The ratio in b' is an evaluation shortcut: with the
eligible-label counter present, available and unavailable coordinates
have distinct positive route factors, and the denominator cancels.
Their quota-conditioned operator is common across p even though b'
itself depends on p. Apply (2) canonically before subset relaxation.

The stronger quartet-by-quartet filter (1) is under development. These
are constraints on the existing global count, not a new encoder or a
return to the parked local seven-type proof.

The aggregate marker is valid but numerically inadequate: log2 eta=1/16
raises the selected 382-state diagnostic from +506781 to +694211.
It should not be retuned as the main approach. Its 2766 exact labeled
monomial checks also verify the telescoping old/new availability term;
the marker need not be at least one at each individual update.

The exact local filter has a positive coefficient formula. Write
o=(1-p)/a and r_route=p+o. For i old available positions and f new
unavailable positions, define

    D_(i,f)=[w^f](1+alpha w)^i(1+w)^(4-i).

The coefficient of x^i y^(4-f) w^f in the hit-quartet polynomial is

    binom(4,i) D_(i,f) b^(4-f) chi^4 r_route^(4-f)
       * sum_(j=0)^min(f,i-1) binom(f,j) p^(f-j) o^j. (3)

The truncated sum imposes (1). At a quota counter z, use pz in covered
eligible factors, retain o on outside positions, and give every
uncovered position its raw weight 1-p+pz. Thus regional conditioning
still extracts exactly h eligible labels and cancels the Bernoulli
normalization. Formula (3) only removes impossible nonnegative terms;
it adds no field suppression factor.

`quartet_joint_mds.py` implements this local operator, positive hit-power
recurrences and its action on exponential potentials. Checks pass 4860
explicit physical/label entries, 900 numerical/quota comparisons, 150
unchanged terminal entries, 252 exact hit-power coefficients and eight
physical/count norm comparisons. An independent GF(8) Cauchy-matrix
check verifies all 69 superregular minors, 4095 nonzero input vectors
and 390369 outside-hole incidences. This toy matrix checks the MDS lemma;
it is not a sample from the actual accepted field-multiplier family.

The exact local filter leaves E_new=empty transitions unchanged, including
the full-availability diagonal. It therefore cannot repair that particular
norm obstruction by itself. Its intended use is in the global joint
count with quota conditioning and all unavailable positions retained.

## A second accepted-map constraint: outside output weight

The accepted map A+delta I is also superregular. In an old-hit quartet
with r old available positions and j positions outside H, the output
has at least max(0,j-r+1) nonzero positions. For j>0 the actual incoming
state must be nonzero: otherwise the full covered quartet would require
a nonzero new state at an outside-H position whose input is zero.
On H^c, Y=(A+delta I)Q_prev; a k-supported nonzero incoming state has
at most k-1 image zeros, and k<=r. Fresh births in unhit quartets each
produce one nonzero output. Consequently

    B=f+sum_(old-hit quartets) max(0,j-r+1) <= wt(Y) <= W. (4)

For eta>=1 the canonical multiplier eta^(W-B) is at least one. Its
global factor is eta^W. Locally discount births by eta^-1 and sum
outside available/unavailable labels jointly, multiplying each hit
term by eta^-max(0,j-r+1). Keep the unavailable-label cap (1) in the
same sum. This is one deterministic event marker, not another field
probability or another charge for the earlier closure penalty.

The helper's optional output_discount (default one) implements this
positive local factor; --outputgap supplies log2 eta. Unlike the
unmarked MDS exclusion, this marker changes E_new=empty terminal
transitions, so marked composition must use the marked endpoint.
Checks additionally verify a GF(256) pair A,A+I with all 138 minors,
1263 deterministic vectors and 138768 output-budget incidences.

## A full-count comparison without retaining a subset

The local polynomial also permits a positive comparison that retains
every actual E count. Fix disjoint bins of (u,e), with the zero values
of each coordinate separated, and a common lift phi(u,e)=t0^u.
Let P=K-Z. For a destination bin J and positive tau,sigma,

    (P H_J)(u,e)
       <= (P[(t0 tau)^v sigma^f])(u,e)
          / min_((v,f) in J) tau^v sigma^f.            (5)

For each physical source and target bin, minimize (5) over a fixed
finite family of tilts, divide by phi at the source, then maximize over
the source bin. This yields PH<=HB_plus. The order of the minimum and
maximum is the tighter valid choice; they need not be numerically equal
in the other order. Because phi=1 for every u=0 state, ZH=HZ_bin exactly.
The lift and bins stay fixed across quota tilts. Compose the regional
operators before taking the quota-adjusted entrywise minimum.

At a fixed new potential, the local polynomial is A(x)+B(z), so its
128th power is a sum of 129 positive outer products. The u=0 part of
P is computed as B^128-B0^128 using a positive difference-of-powers
recurrence. No nearly equal positive quantities are subtracted.
The actual free-terminal action is P applied to 1[E_new=0], evaluated
at new t=1,s=0, including the output marker when enabled.

`quartet_joint_chernoff.py` implements (5) on all 131841 physical counts.
The initial checks pass 128 exact comparison entries, 240 physical-bin
inequalities, 60 positive actions, 30 terminal entries and four exact
zero intertwinings, including marked outputs and rational 16-step
regional comparisons. Subsequent checks add independent width-eight
two-quartet convolutions. No dense physical transition is allocated.

The first width-32-bin probe uses 169 bins, nine destination tilts and
three quota tilts. It gives +920524 FLOAT log2 units in about 30 seconds.
Thus full E retention alone does not make this coarse comparison tight.
Expanding the destination family to 49 tilts gives +716079 in about
143 seconds. This still does not close the selected row class.

The implementation now also uses exact zero-coordinate slices: targets
with e=0 can use P(t0 tau,0), and targets with u=0 can use P(0,sigma).
The destination (0,0) has no positive transition and is kept exactly zero.
An optional piecewise potential phi=t0^u s0^e for u>0 and phi=1 for u=0
preserves the same zero intertwining. Positive-u targets use the shifted
interior actions; zero-u targets MUST use the unshifted exact axis alone.
Using the shifted interior action for those targets when s0<1 would be
an invalid underbound. Marked endpoints divide by the same source phi.

Final exact checks pass 256 comparison entries, 480 physical-bin
inequalities, 180 positive actions, 540 independently convolved width-eight
checks, 60 terminal entries and eight zero intertwinings. Underflow guards
distinguish structurally zero axis rows from required positive actions.

A piecewise probe with log2 t0=-1/8 and log2 s0=-2 gives +4805700.
The normalization is valid (independently audited), but this choice
greatly worsens the upper comparison. No conclusion about actual code
distance follows from that failure.
With both potentials zero, the exact-axis version gives +918902 in
about 44 seconds, versus +920524 before axis sharpening. This gain is
small. The full-count comparison is sound but has not improved the best
existing global bound.

## What would constitute a genuine structural obstruction

There is a separate deterministic test that must not be confused with
a positive upper bound. Fix any allowed state masks U_0=empty, not
necessarily canonical availability masks. Define

    V_t=Gamma_t(U_(t-1)) union U_t,
    L_U=sum|U_t|,
    D_U=sum_(t=1)^T |U_(t-1)\U_t|,
    R_U=sum_i min(128,|V intersect row i|).

If L_U>R_U and L_U+D_U<=W, the geometry itself forces a low-weight
codeword for every choice of local field maps. Indeed, states supported
on U form an L_U-dimensional space. The causal invertible map Q to X
places X inside V. The outer parity checks restricted to V have rank
R_U, so their pullback has a nonzero kernel when L_U>R_U. The resulting
nonzero output lies inside U_(t-1) union U_t, of total size L_U+D_U.
The terminal remains free. This criterion passed independent audit.

No such masks have been found here. A formal profile near R=L, a
positive transfer bound, or a failed proof does not supply them.
Useful necessary screens are F>=L_U+D_U, R_U>=F/2, and R_U at least
the sum of the 128 largest regional coverage counts. The last inequality
follows by choosing those same 128 regions in each outer row. This
formulation separates a potential matroid-expansion obstruction in the
geometry from loose field counting.

The obvious unconstrained cut relaxation cannot decide this question.
For lambda>=0, the objective R_U-L_U+lambda(L_U+D_U) has a directed
min-cut representation: availability nodes imply their covered input
neighborhood; each covered position has a capacity-one edge to its
outer-row node; row nodes have capacity 128 to the sink. Time edges
charge lambda D_U, and unary terms charge (lambda-1)L_U. There is no
terminal-reset edge.

However, its minimum is exactly min(0,N(lambda-1/2)) for every geometry.
Writing C=L_U+D_U, the screens above imply R_U>=C/2 and L_U<=C<=N,
so the objective is at least (lambda-1/2)C. Empty masks attain zero;
full masks attain N(lambda-1/2). Thus an unpinned parametric cut search
would merely rediscover those two extremes. It cannot establish or
refute the required size-constrained expansion at C<=W. This negative
control is proved analytically; no large graph was constructed.

## Next comparison and stopping point for this diagnostic round

Do not continue broad unconditioned potential searches on these results.
The next useful comparison is on common joint (R,L,f,F) classes. The
older global bound and the new full-E comparison may be dominated by
different relaxed profiles; one large positive total does not show
which bound controls a fixed class. The previously suggested box centered
at (900000,900000,20000,1200000), with half-width 1024, was subsequently
found EMPTY: F_min-R_max=297952>128h=262144. It must not drive further
proof work. The replacement center uses F=1150000; it passes the
aggregate necessary conditions without establishing realizability.
Apply each geometric count and its valid conditional field bound
to that SAME class, then take the smaller class bound. Do not take a
kernel minimum across different physical bases.

For this test the geometric interface is explicit. Set a=2^loga,
b=2^logb, m=b*2^logu, chi=2^logchi, alpha=1 and output_discount=1.
Use eligible weight p and outside weight (1-p)/a. The only external
geometric factor is log2 binom(8192,h)-128h*loga. The resulting moment
counts a^-R b^L (2^logu)^f chi^F. Subtract the minimum of
-loga R+logb L+logu f+logchi F over the chosen box, then add its
independently optimized conditional field bound. There is no theta,
output-weight prefactor, gamma/q, or field-endpoint correction inside
this geometric moment. Do not cap its expected count at one. The
existing markers (2, approximately 25/16, -7/16, 0) provide a fixed
first comparison against the previously positive hard box.

There is also a useful negative result for the output-budget marker.
Its deterministic budget obeys

    B <= 5F/4-L-f/4+|A_T|,

because max(0,j-r+1)<=5-r in each hit quartet. At the previously
recorded middle profile (F,L,f)=(1132763,891035,71426), even allowing
all 512 terminal coordinates gives B<508000, far below W=943718.
Thus the output-budget condition excludes no geometrically consistent
history in that class; improving a relaxation would be a separate issue.
Moving the marker into the older
tagged transfer would also require a common quota basis: its exact
old-occupancy-dependent route coefficient would otherwise make the
previous common factor depend on p. A naive reduced-kernel minimum
would not be justified. These are reasons to retain the marker as a
proved auxiliary constraint without treating it as the next closure plan.

The full certificate remains open. All jobs from this round have exited.
All prior class theorems and frozen helpers are preserved; no genuine
bad-geometry witness or low-distance codeword has been established.

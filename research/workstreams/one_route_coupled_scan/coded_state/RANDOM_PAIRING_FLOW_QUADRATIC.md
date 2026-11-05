# A quadratic histogram potential and a directed local bound

The [two-update diagnostic](RANDOM_PAIRING_FLOW_BLOCK.md) motivates keeping
auxiliary-state diversity between updates. This continuation obtains a
pointwise upper integrand for a quadratic potential. Its numerical screen
has slack for the selected1974-row structural tail. A directed analytical
envelope now certifies the integral at the worst histogram found by that
screen and on explicit neighboring regions. Uniform histogram coverage
remains open. This is not a distance certificate.

The encoder is unchanged. All histograms in this note describe the
five-state causal flow, not the encoder's field-valued state or its support.

The [joint integral continuation](RANDOM_PAIRING_FLOW_CONTROL_POINTS.md)
retains population weights across Gaussian annuli. With an outcome-based
global curvature bound, it closes a larger tight boundary region and a
4,884,615-histogram interior band. A disjoint partial partition and wider
failed controls are recorded there. Uniform coverage remains open.

## Candidate and rejection tests

Keep n=512, p=1974/8192, and define

    phi(h)=u dot h + kappa ||h||^2/n,
    theta=1/10, eta=1/5, kappa=2/25,
    u=(0,-84450914,-150416955,-196990597,-222300536)/10^9.

These log weights are now fixed rationals. They were selected by a spectral
screen; subsequent checks need not rely on that search. For one update,
let H' be the next histogram, L its overflow loss, and J its target-label
count. The required conditional moment is

    E exp(theta L + eta(J-np) + phi(H')-phi(h)).

First test lower bounds on this moment. Jensen's inequality uses the exact
mean of L, H', and ||H'||^2. The latter needs the joint distribution of two
distinct matched pairs, not independent partner states. The helper averages
the three pairings of four drawn coordinates without replacement. Under
regional sampling it also uses the joint tag-count distribution on those
four positions. Exact width4 comparisons cover all70 histograms under both
IID and fixed-count labels:140 tests pass.

The initial kappa=0.01 candidate fails. At h=(192,256,64,0,0), an exact
rational Jensen check makes the selected union expression exceed exp(35),
even before endpoint costs. It uses the exact row count >2^6519, the
regional tag probability <=2^(-6.5), and ln2>69/100. This rejects that
potential at the fixed markers and weights, not the code.

Larger penalties survive the mean check. A stronger lower bound changes
the pair-tag probabilities while keeping the matching uniform:

    log E_P exp(Z) >= E_Q Z - D(Q||P).

Conditional on the matching, Q assigns independent tag counts to its pairs.
The helper computes the finite-population quadratic mean and relative
entropy for any such proposal. Iteration merely finds a useful proposal;
validity of the inequality does not require convergence to an optimizer.

At kappa0.08, a495-point histogram grid followed by local refinement checks
591 histograms. Its largest sampled lower bound is about7.62493 at
h=(185,188,83,56,0). A4096-sample importance probe there estimates the log
moment at7.62810, with estimated log standard error0.00119. This sampling
result is a diagnostic, not an upper bound or certified confidence interval.

## A matching bound that retains rare colors

We next bound an actual uniform matching. Fix an input histogram h and
a symmetric real pair-score matrix f. Put

    U_f(M)=sum_({i,j} in M) f_(state_i,state_j).

For a set A of available colors, define

    Delta_A(f)=max_(a,b in A)
      [max_(x in A)(f_xa-f_xb)-min_(x in A)(f_xa-f_xb)].

Expose a matching by choosing one unmatched vertex, then revealing its
partner. Compare two possible partners a,b by coupling completions through
the switch (x,a),(b,c) -> (x,b),(a,c). This is a bijection between the two
conditional uniform completion spaces. The score difference is at most
Delta_A(f) in absolute value. Thus the next Doob increment has conditional
range at most Delta_A(f). The final pair is determined and costs no range.

A constant range over all n/2-1 exposures is wasteful when an extreme color
is rare. Fix a priority order of the colors and always choose a vertex of
the first remaining color. The phase for color a has at most h_a exposures.
After that phase, a is absent and Delta cannot increase.

For a priority order pi, form a deterministic squared-range budget by
allocating at most h_pi(i) exposures to each successive color, stopping
after n/2-1 exposures. Fill earlier phases first. This maximizes the budget
under these quotas because the ranges decrease when colors are removed.
Let B_pi(h,f) be one eighth of that budget. Conditional Hoeffding bounds
and their exponential supermartingale then imply

    E_M exp(U_f(M)) <= exp(E_M U_f(M) + min_pi B_pi(h,f)).

The priority order may be chosen separately for each fixed matrix f.
The mean E_M U_f is the exact average over unordered pairs drawn from h.
Exact tests check the deterministic budget on2520 matching/priority paths
at width8. The switch coupling and conditional-range argument establish
the bound at arbitrary width; the tests do not substitute for that argument.

## A Gaussian integral for the quadratic potential

Let P_j=binom(2,j)p^j(1-p)^(2-j). For input colors a,b and tag count j,
write o_abj for the outgoing count vector and ell_abj for its loss.
Choose any center t with sum t_a=n, and set beta=sqrt(2kappa/n). For
Z distributed as a standard five-dimensional Gaussian, define

    f_ab(Z)=ln sum_j P_j exp(theta ell_abj + eta j
                   + (u+2kappa t/n+beta Z) dot o_abj),
    C=-eta np-phi(h)-kappa ||t||^2/n.

The Gaussian square identity and conditional independence of tag counts
give the exact representation

    E exp(theta L+eta(J-np)+phi(H')-phi(h))
      = exp(C) E_Z [exp(-beta Z dot t) E_M exp(U_f(Z)(M))].

Insert the preceding matching bound. The resulting positive integrand

    exp(C - beta Z dot t + E_M U_f(Z) + min_pi B_pi(h,f(Z)))

is a pointwise upper integrand for the desired moment. The common Gaussian
direction cancels because each o_abj sums to two and t sums to n. Only four
orthonormal contrasts need integration. The chosen center is close to the
mean under the variational proposal, but no fixed-point assumption is used.

This interface preserves the matching distribution and the histogram.
It does not take an independent worst case for every product potential
inside the Gaussian mixture.

## Numerical screen, including endpoint costs

The original unstratified matching range was too loose. Excluding absent
colors and then exposing rare colors first reduces its loss substantially.
With the fixed rational weights above, an order3 four-dimensional Gaussian
quadrature screen checks592 grid/refined histograms. The worst tested
histogram remains (185,188,83,56,0).

At that histogram, higher orders give:

| Quadrature order | Log of the upper-integrand integral, numerically |
|---|---:|
| 5 | 7.6482438954 |
| 7 | 7.6482715785 |
| 9 | 7.6482769737 |

These are not directed upper bounds on the integral. Agreement between
orders does not certify integration error or Gaussian tails.

The endpoint potential ratio is bounded by minimizing
u dot x+kappa||x||^2 over the continuous simplex. Relaxing integer
histograms this way is valid. Its numerical cost is51.77773745 nats.
The ignored initial input batch costs1.93683432 nats in the centered tag
moment. Including both, the target-row union and regional conditioning
leaves a per-update log-moment allowance7.6791686911.

If the order9 value were established as a uniform upper bound over every
histogram, the selected structural union exponent would be approximately
-182.50 bits. Neither premise is established. In particular, the592 tested
histograms are not exhaustive coverage of the512-coordinate simplex.

## An analytical Gaussian tail envelope

The integral has a useful global tail bound. A gate's outgoing count
vectors have Euclidean diameter at most2sqrt(2). Consequently the Hessian
of its log moment, as a function of Z, is at most2beta^2 I. Summing n/2
gates bounds the Hessian of E_M U_f(Z) by2kappa I.

Write delta=Delta_A(f(0)), where A contains the colors present in h, and
L_s=4sqrt(2) beta. Each four-score switch has gradient norm at most L_s,
so Delta_A(f(Z))<=delta+L_s||Z||. Let q=n/2-1 and let

    epsilon = ||gradient_Z E_M U_f(0)-beta t||,
    A_0 = 1/2-kappa-q L_s^2/8,
    B_0 = epsilon+q delta L_s/4,
    D_0 = q delta^2/8,
    C_0 = C+E_M U_f(0).

Use the weaker unstratified range for this tail bound. For R with
A_R=A_0-B_0/R>0, radial integration in four dimensions yields

    integral_(||Z||>R) upper_integrand dGaussian(Z)
      <= exp(C_0+D_0) (1+A_R R^2) exp(-A_R R^2)/(4 A_R^2).

The derivation uses B_0 r <= (B_0/R) r^2 for r>=R. It does not assume
epsilon is zero. At the selected histogram A_0=0.10125. Floating evaluation
gives log absolute tail bounds about-8.8275 at R16 and-39.0334 at R24.
These evaluations still require outward rounding before certification.
The central integral and uniform histogram coverage remain unresolved.

## Directed analytical envelope at the selected histogram

Fix h=(185,188,83,56,0) and the rational potential above. The moment uses a
uniform matching and IID Bernoulli(p) target labels. The directed check proves

    E exp(theta L+eta(J-np)+phi(H')-phi(h)) < exp(7.675753).

This is a one-update statement at this particular histogram. It does not
assume that the process stays there, or bound transitions from other histograms.
The construction and accepted pair matrices are unchanged.

The first attempt bounded all Gaussian directions by their worst direction.
Its floating log bound exceeded7.74 and failed. Retaining directional
information gives an analytical envelope that can be integrated explicitly.
The following derivation replaces numerical quadrature at this histogram.

Work in the four orthonormal Gaussian contrasts already defined. For each
pair type, let pi_j be its normalized weights at Z=0. Let V_j be the projected
outgoing vector, mu its pi-weighted mean, and v_j=V_j-mu. Jensen's inequality
gives sum_j pi_j exp(v_j dot Z)>=1. Hence, on ||Z||<=R,

    Hessian f_ab(Z)
      <= sum_j pi_j exp(R||v_j||) v_j v_j^T =: H_ab,R.

The inequality is in positive-semidefinite order. To obtain it, bound the
tilted covariance by its second moment around mu, then bound each numerator.
Summing with the exact expected pair counts gives a matrix C_R. If G(Z)
denotes C-beta Z dot t+E_M U_f(Z), then

    G(Z) <= G(0)+gradient G(0) dot Z+Z^T C_R Z/2.

Next bound the matching switches without replacing all gradients by their
norms. A switch s is a signed sum of four pair scores. Write its affine
part as w_s dot (1,Z). For each available-color set A, choose a matrix M_A
and verify M_A-w_s w_s^T is positive definite for every switch on A.
This implies |w_s dot (1,Z)|<=sqrt((1,Z)^T M_A(1,Z)).

Let B_ab be the diameter of the projected outgoing vectors. Likelihood
ratios between the tilted and untilted three-point laws lie between
exp(-B_ab R) and exp(B_ab R). The variational characterization of covariance
therefore gives the same two-sided comparison for their covariance matrices.
For a switch with pair coefficients c_ab, a valid Hessian norm bound is

    D_s(R)=||Hessian s(0)||_F
           +sum_ab |c_ab|(exp(B_ab R)-1) trace(Hessian f_ab(0)).

Let D_A(R) be the maximum over switches on A. Taylor's theorem yields

    Delta_A(f(Z))
      <= sqrt((1,Z)^T M_A(1,Z))+D_A(R)||Z||^2/2.

For any epsilon>0, square this bound using
(a+b)^2<=(1+epsilon)a^2+(1+1/epsilon)b^2. On the radius-R ball,
replace ||Z||^4 by R^2||Z||^2. The resulting bound is quadratic in (1,Z).
Use the fixed color priority (3,2,0,1), whose exposure quotas are56,83,116.
Add these bounds with their quota/8 factors to obtain a5-by-5 matrix A_R.

Combining the two quadratic envelopes gives a Gaussian precision matrix
P_R=I-C_R-2(A_R)_zz and linear term b_R=gradient G(0)+2(A_R)_0z.
Whenever P_R is positive definite, integration over all four dimensions is

    exp(G(0)+(A_R)_00)
      det(P_R)^(-1/2) exp(b_R^T P_R^(-1)b_R/2).

This full integral bounds the contribution from ||Z||<=R. For a shell with
lower radius l, multiply by an upper bound on its Gaussian probability.
Put mu_R=P_R^(-1)b_R and verify P_R>=lambda_R I. The event ||Z||>=l implies
that a standard four-dimensional Gaussian has norm at least
sqrt(lambda_R) max(0,l-||mu_R||). Its exact squared-norm tail gives the factor
(1+x)exp(-x), where x=lambda_R max(0,l-||mu_R||)^2/2.
Ignoring the shell's upper radius only enlarges this bound.

The verifier uses the central ball of radius4 and shells ending at5,6,8,12,16.
The global radial formula above bounds the tail beyond16. Although each
local quadratic integral extends past its shell, summing them remains an
upper bound on the disjoint partition of the original integral.

[The floating envelope helper](random_pairing_flow_envelope.py) proposes
matrices. [The directed verifier](random_pairing_flow_envelope_verify.py)
converts them to rational decimal candidates and checks every matrix
inequality with interval LDL elimination. It uses outward-enclosed Decimal
operations for arithmetic, square roots, exponentials, logarithms, solves,
and determinants. Floating optimization and eigenvalues propose candidates
only; they are not accepted without the directed matrix checks. The center
is also a rational candidate, with its fifth coordinate fixed so its sum
is exactly512. No numerical integration or unbounded tail is left here.

The65- and80-digit replays give a log-moment upper bound7.675752586290.
The80-digit replay also checks183 arithmetic cases and recomputes the
row-union, regional-conditioning, and endpoint costs with directed bounds.
The resulting per-update allowance is greater than7.679168691101.
If the proved local bound held uniformly over all histograms, it would
leave more than20.18 bits in this selected structural calculation.
That uniform premise is not proved, so this is not a20-bit certificate.

The directed replay takes0.616 seconds on Peach. The earlier broad radial
envelope failed, while the directional envelope passed; this was a proof
bound improvement, not a construction or performance change.

## Directed bounds for neighboring histogram regions

The next continuation extends the verifier from a single histogram to
constrained boxes. All16 histograms obtained by moving one coordinate of
(185,188,83,56,0) into another color pass individually. The largest bound
among them is below7.675895, at (184,188,83,57,0). In particular, introducing
one coordinate of the fifth color does not break the local check.

Individual checks cannot provide scalable coverage. Define a box B(r,l,u)
by integer coordinates

    |h0-185|<=r, |h1-188|<=r, |h2-83|<=r, l<=h4<=u,
    h3=512-h0-h1-h2-h4.

Every box below has nonnegative coordinates. Its cardinality is
(2r+1)^3(u-l+1). Choose a reference histogram with h4=floor((l+u)/2),
and choose one rational Gaussian center for the whole box. Each histogram
in the box still uses exactly the fixed potential, markers, and transition
distribution specified above; only the analytical bound is recentered.

The expected matching score has the exact quadratic form

    E_M U_f = [h^T F h - sum_a h_a F_aa]/[2(n-1)],

where F_ab=f_ab. Expand this form in changes to h0,h1,h2,h4, with h3 absorbing
their negative sum. This preserves the cancellation with -kappa||h||^2/n.
Use the same expansion for the Gaussian gradient and curvature matrix.

Use color priority (4,3,2,0,1). In the verified boxes, the exposure quotas are

    h4, h3, h2, 255-h4-h3-h2.

They sum to255 exactly. The verifier retains their affine dependence instead
of independently maximizing each quota. Consequently the Gaussian exponent
has a quadratic core, an affine exposure contribution, and a quadratic
linear-term vector b(h). Interval matrix elimination encloses P(h)^(-1)
and log det(P(h)) throughout the box. Combining the core with
b(h)^T P(h)^(-1)b(h)/2 produces an interval polynomial of degree at most4.
Its coefficients are intervals because the inverse matrix varies with h.
For each fixed h, those intervals contain the coefficients of its expression.

Evaluate the combined polynomial using tensor Bernstein coefficients on the
box. The Bernstein basis is nonnegative and sums to one, so the greatest
coefficient upper endpoint bounds every value. This retains cancellations
that were lost when the core, gradient, and squared Gaussian shift were
maximized separately. The matrix and tail checks remain directed.

The80-digit replay gives these rounded-up log-moment bounds:

| Box (r,l,u) | Integer histograms | Log-moment upper bound | Below allowance? |
|---|---:|---:|---|
| (1,0,0) | 27 | 7.678556 | yes |
| (2,0,0) | 125 | 7.682042 | no |
| (2,2,2) | 125 | 7.660477 | yes |
| (4,4,8) | 3,645 | 7.650837 | yes |
| (8,8,16) | 44,217 | 7.668881 | yes |

The failed box is a failure of this upper bound, not a code counterexample.
These boxes overlap and do not cover the histogram simplex. In particular,
the large passing box does not remove the unresolved zero-fifth-color region.
No minimum-distance margin follows from these local moment bounds alone.

The progression identifies two avoidable losses. Plain interval evaluation
of separate terms failed even on27 histograms. Retaining polynomial
cancellation closed that box. Keeping the original Gaussian center for the
largest box gave a bound above8.35; recentering the analytical representation
reduced it below7.668881. Neither change modifies the encoder or its setup.

The serial box replay takes10.448s; the final combined neighbor/box replay
takes21.653s and passes399 arithmetic and polynomial checks. Polynomial tests compare the exact matching quadratic
and its square with direct rational evaluation at54 points; the region
proof uses the Bernstein identity, not point sampling. The verifier also
asserts positivity for every matrix inequality and checks the infinite tail.
An initial interval eigenvalue candidate was too aggressive for a region;
the verifier now pads candidates and accepts only those whose matrix
inequality passes. A zero-population guard was corrected to use the exact
integer domain rather than a tiny negative outward-rounding endpoint.

Next subdivide or tighten the failing zero-fifth-color region and generalize
the region argument away from this reference neighborhood. A full histogram
cover is still required before multiplying conditional bounds across updates.

## Closing the local boundary and handling exposure-phase changes

The formerly failing box B(2,0,0) is covered by eight radius1 boxes. Their
first-three-coordinate centers are (185,188,83)+(e0,e1,e2), with each e_i
in {-1,1}. Each coordinate interval is either [center-2,center] or
[center,center+2], so the eight boxes cover all125 integer histograms.
Every box passes. Before the shell refinement below, the largest directed
log bound was7.678944. This confirms that the failed radius2 bound was
region-evaluation slack, not evidence of a failing point.

There is no eightfold probability penalty. These boxes cover possible
conditioning histograms for the same one-update expectation; take the
largest conditional moment bound, rather than summing probabilities.

### Removing Gaussian mass outside each shell

The earlier shell bound integrated beyond the upper radius. Let P be its
Gaussian precision and mu its mean. Verify lambda_min I<=P<=lambda_max I.
For a standard four-dimensional Gaussian W, write Z=mu+P^(-1/2)W.
The triangle inequality gives

    Pr(||Z||>=l)
      <= T(lambda_min max(0,l-||mu||)^2/2),
    Pr(||Z||>R)
      >= T(lambda_max (R+||mu||)^2/2),
    T(x)=(1+x)exp(-x).

Subtract the second bound from the first to bound the shell l<=||Z||<=R.
Use one for the first term when l=0. Both spectral inequalities and both
probability terms are checked with outward rounding. The unbounded tail
beyond16 still uses the independent global radial bound.

At the original selected histogram, this reduces the directed log bound
from7.675753 to7.672287. If the latter bound held uniformly, the selected
structural expression would have40.66 bits of slack. It does not hold
uniformly by any result established here, so this is not a40-bit certificate.
The previous passing boxes remain passing. The unsplit radius2 boundary
box still exceeds the allowance, at7.679669; subdivision remains useful.

### A bound across changes of exposure phase

The box verifier now accepts different first-three-coordinate centers and
arbitrary color priorities. With a fixed final phase, its quotas are affine
functions of h. This permits regions outside the original neighborhood.
A region can, however, cross the point where the255 exposures are exhausted.
Restricting every region to one phase would force unnecessary subdivision.

For each phase i through a chosen reference final phase j, let M_i be its
quadratic upper bound on the per-exposure squared range divided by8. The
bound is evaluated on w=(1,Z). Enlarge these matrices so

    M_0 >= M_1 >= ... >= M_j

in positive-semidefinite order. The verifier works backward, adding a
nonnegative scalar multiple of the identity to each matrix. A Gershgorin
bound proposes the shift; interval LDL verifies each resulting difference.
The original per-phase bounds remain valid because the shifts are positive.

Let S_i=sum_(k<=i) h_(priority_k). The true squared ranges decrease as colors
are removed. Consequently phases after j may also use M_j. A quadratic
majorant for the capped exposure budget is

    255 M_j + sum_(i<j) min(255,S_i)(M_i-M_(i+1)).

Each difference is positive semidefinite. Replacing min(255,S_i) by S_i
therefore enlarges this matrix, giving one affine majorant valid on both
sides of every intervening phase boundary. Its final affine quota may be
negative outside the reference phase; validity comes from this matrix
inequality, not from interpreting that quota as a number of exposures.
The code never substitutes a majorant into a negative actual quota.

For boxes contained in a single phase, the verifier retains the sharper
original matrices. The nested construction is a fallback for crossing boxes.
It does not change the random matching or any field operation.

Selected directed checks away from the original neighborhood are:

| First-three center | Radius | Fixed h4 | Histograms | Log upper bound |
|---|---:|---:|---:|---:|
| (48,48,48) | 2 | 48 | 125 | 1.847460 |
| (64,320,64) | 2 | 32 | 125 | 2.657955 |
| (100,100,100) | 4 | 100 | 729 | 7.185916 |
| (64,64,64) | 4 | 64 | 729 | 4.992448 |

The last box crosses exposure phases and uses the nested majorant. The
others use their respective fixed-phase formulas. All are below the same
allowance, but they are selected regions, not a complete simplex cover.

The phase/shell replay passes5018 arithmetic, polynomial, phase and Gaussian
checks at80-digit precision in10.331s. Exact small-composition tests check
the affine-extension inequality for every candidate final phase. Isotropic
Gaussian tests compare the shell bound with its exact chi-square-four mass.
The algebra and spectral inequalities establish the bounds; the tests
independently check their implementations.

The final80-digit replay of the eight boundary boxes takes17.820s and
reproduces all passes with the sharper shell formula. Its uniform log bound
over the125-histogram local box is below7.676047.

Next construct a complete, explicitly accounted histogram cover. The current
tool can change its reference center and cross exposure-phase boundaries;
the remaining task is coverage of the entire simplex, including its faces,
not an additional local phase formula. Any uncovered region must remain
listed as open. Only after this one-update interface is uniform may it feed
the full selected routing tail, and that still does not close all distance
profiles or field-cancellation counts.

## Accounted coverage and a coarse radial envelope

The next obstacle is coverage, not another selected histogram. Let
H={h in Z_{≥0}^5: sum h=512}. It contains C(516,4)=2,919,618,945 points.
The entries of h count auxiliary flow states; they are not field symbols.
`random_pairing_flow_cover.py` represents a region by five inclusive integer
intervals intersected with H. It tightens their endpoints using sum h=512.

For lower endpoints l_i, upper endpoints u_i, residual r=512-sum l_i,
and widths w_i=u_i-l_i+1, the exact region size is

    sum_(S subset {0,...,4}) (-1)^|S| C(r-sum_(i in S) w_i+4,4).

A term is zero when its residual is negative. Splitting one interval into
two disjoint integer intervals gives an exact partition. The bounded search
maintains certified+pending+unclosed_singletons=total. Neither exact
accounting nor a partial search implies that all points pass.

The verifier may eliminate any one coordinate using sum h=512. This permits
faces where the previously eliminated coordinate is fixed at zero. Exact
small-domain tests check counts, splits, reference points and all five
coordinate choices. An absent-state3 face with h0,h1,h2,h4 in[124,132]
contains489 histograms and passes with log bound below5.782922. Two vertex
neighborhoods of70 histograms also pass.

### A cheap bound for broad regions

Write z for the four-dimensional auxiliary Gaussian variable, and r=||z||.
For exposure phase i, bound its switch range by d_i+L_i r. Here d_i is
computed with directed rounding at z=0 and L_i=4 sqrt(2) beta. A phase
with one remaining color has both terms zero. Backward maxima make d_i
nonincreasing; L_i is already nonincreasing.

The three coefficients of (d_i+L_i r)^2/8 are therefore nonincreasing
across phases. Apply the capped-prefix extension above to each coefficient.
This gives affine population expressions bounding the constant, linear,
and quadratic terms of the exposure correction. Bound these expressions
over the region, and cap each by its full-color255-exposure bound. Denote
the resulting nonnegative upper bounds by c0,c1,c2.

Let G(h,z) denote the Gaussian core, including -||z||^2/2 but excluding
the matching exposure correction. The core without that Gaussian density
has Hessian at most2 kappa I. Directed bounds on G(h,0) and its z-gradient
give

    G(h,z)+correction <= C-a r^2+b r,
    C=upper G(h,0)+c0,
    a=1/2-kappa-c2,
    b=upper ||gradient_z G(h,0)||+c1.

For any fixed0<epsilon<a, Young's inequality gives the normalized
four-dimensional Gaussian integral bound

    log moment <= C+b^2/(4 epsilon)-log(4(a-epsilon)^2).

The verifier proposes epsilon in floating point, fixes its decimal value,
and checks the displayed inequalities with outward rounding. No quadrature
or finite-radius truncation is used in this coarse bound.

Near a vertex, the region is a translated simplex rather than a box.
When every coordinate width equals512-sum l_i, use total-degree Bernstein
coefficients on that simplex to bound the population polynomials. For
a monomial u^k of total degree at most D, its coefficient at multi-index
alpha is product_i falling(alpha_i,k_i)/falling(D,sum k_i). The implicit
fifth barycentric coordinate is1-sum u_i. The convex-hull property then
gives the bound. General clipped boxes retain tensor Bernstein bounds.

The80-digit replay gives the following log upper bounds. Entries are
rounded upward here; the allowance is greater than7.679168691101.

| Maximum deficit from dominant color | Histograms per cap | Color0 | Color1 | Color2 | Color3 | Color4 |
|---|---:|---:|---:|---:|---:|---:|
|64|814,385|2.043869|0.117715|-2.767517|-1.202296|-0.923862|
|96|3,921,225|8.992811|6.009838|2.857256|4.152069|4.747205|
|128|12,082,785|16.777719|12.307711|9.369855|10.408171|11.167153|

Use deficit64 for color0 and deficit96 for each other color. These five
caps are disjoint because any two dominant thresholds sum to more than512.
They certify16,499,285 histograms. Their complementary region has
h0<=447 and h1,h2,h3,h4<=415, containing2,903,119,660 histograms.
The ladder replay completes in3.190s and passes696 exact implementation
checks. A bounded global continuation verifies7 regions, certifies
16,505,404 histograms, and leaves2,903,113,541 pending after20 attempts.
Its worst certified log bound is below7.249294. No singleton is rejected.
Its stopping point depends on the wall guard; the five initial caps provide
the fixed reproducible coverage statement. No search tree is saved.

### Why joint concavity is not a shortcut

One possible shortcut would maximize G(h,z) using joint concavity in the
histogram and Gaussian coordinates. The natural real extension of this
core is not jointly concave. Directed second-derivative checks give:

- At h=(185,188,83,54,2), a direction in both h and z has curvature
  greater than0.000264213334646.
- At h=(506,2,2,2,0), a direction keeping the fifth population fixed has
  curvature greater than0.000152164937913.

The script fixes rational directions before verifying each inequality.
These witnesses concern the core excluding the matching correction. They
reject this concavity argument, not the moment bound or minimum distance.
The remaining coverage needs sharper interior bounds or subdivision; the
current verified caps cover only about0.57% of H.

## Retaining population dependence through the Gaussian integral

The preceding coarse bound maximizes the core and the gradient separately.
That can be expensive: both depend on the same histogram. Keep the linear
Gaussian term until integration, rather than replacing it by its norm.
For a fixed region and a fixed Gaussian center, write C(h)=G(h,0) and
g(h)=gradient_z G(h,0). Both are quadratic population polynomials.
Let B0(h) be the affine upper bound on the constant exposure term. Let
c1,c2 be uniform nonnegative bounds on the other two exposure terms.

Suppose the core's Hessian without Gaussian density is at most H I in the
integration region. Set a=1/2-H/2-c2 and require a>0. For0<epsilon<a,
write gamma=a-epsilon. Young's inequality applies only to the exposure
term c1||z||. Completing the square for g(h) dot z then gives

    log integral <= upper_h [C(h)+B0(h)+||g(h)||^2/(4 gamma)]
                    +c1^2/(4 epsilon)-log(4 gamma^2).

The bracketed expression is a degree-four polynomial, bounded jointly by
Bernstein coefficients. This is not the sum of separate maxima. Scalar
epsilon is proposed numerically, then fixed and checked with intervals.
For a global envelope take H=2 kappa. For a ball-local envelope use the
local covariance bound already derived above. The latter cannot be used
outside its ball.

### Ball-local bounds and an independent tail

The local scalar bound computes the matrix of upper covariance expressions
over the population region. A Gershgorin bound gives H; take the smaller
of that bound and the global2 kappa bound. For each switch s, Taylor's
theorem on ||z||<=R gives

    |s(z)| <= |s(0)| + (||gradient s(0)||+R D_s(R)/2)||z||.

Here D_s(R) is the previously justified bound on its Hessian norm. Cap
the slope by the global Lipschitz bound and take backward maxima across
exposure phases. Applying the same affine quota extension gives c1,c2
and B0(h) for this ball. The integrated envelope may extend beyond the
ball because its integrand is nonnegative, but it bounds only the actual
contribution from inside the ball.

Outside the ball, use the independent global envelope. Completing its
square produces a Gaussian with mean mu(h)=g(h)/(2 gamma) and covariance
I/(2 gamma). If m bounds ||mu(h)|| over the population region, then

    Pr(||Z||>=l) <= (1+x)exp(-x),
    x=gamma max(0,l-m)^2.

Multiply the integrated global envelope by this factor at l=R. The
implementation bounds m using the population polynomial ||g(h)||^2.
It also retains the older radial tail formula and takes the smaller valid
upper bound. Retaining population dependence in this tail matters: a
separate core/norm bound can swamp an otherwise passing local integral.

The same lower-radius factor can bound a local annulus l<=||z||<=R.
Sum nonnegative bounds for disjoint annuli, then the independent global
tail. The optional scalar diagnostic uses radii4,8,16. It does not replace
the sharper anisotropic matrix method.

### Checked regions and the remaining directional loss

The global coupled bound alone certifies disjoint caps with dominant
thresholds(384,384,352,352,352), containing111,318,573 histograms. Adding
the local bound and coupled tail enlarges them as follows. Log bounds in
the table are rounded upward; all are below7.679168691101.

| Dominant color | Minimum population | Histograms | Log upper bound |
|---|---:|---:|---:|
|0|352|29,051,001|5.892285|
|1|352|29,051,001|5.672818|
|2|320|59,626,385|3.764036|
|3|320|59,626,385|5.803631|
|4|320|59,626,385|5.896333|

These five caps are disjoint and contain236,981,157 histograms. Their
complement contains2,682,637,788. The80-digit bounded cover rechecks each
cap, attempts12 regions in18.196s, and maintains the exact accounting
identity throughout. It certifies no additional region in that run.
Coverage is about8.12%, not a probability of success or proof completion.

The scalar annulus bound certifies the16,730-histogram region
96<=h_i<=108, sum h=512, below7.665537. The wider region88<=h_i<=116
still fails that bound. More importantly, the scalar bound at the single
histogram(185,188,83,56,0) is7.851141, above the allowance. The preserved
matrix bound at the same point is7.672287 and passes. Subdividing to a
point therefore does not repair this scalar loss: the anisotropic matrix
information is necessary for the present envelope near that boundary.
This is a distinction between proof bounds, not a code counterexample.

The cover passes711 domain, polynomial and Gaussian checks. New controls
include exact isotropic integrals and tails, plus a polynomial example
whose core cancels ||g(h)||^2/(4a) identically. The latter verifies that
the implementation retains the intended population cancellation. The
separate80-digit matrix regression passes5018 checks and reproduces its
selected-point bound. Numerical runs remain serial, bounded and file-free.

Next use cheap bounds for broad regions and preserve matrix directions
near tight cases. Most histograms remain open; other target-row counts,
support profiles and field-cancellation counting are not closed by this
one-update calculation. The construction and all prior controls remain
unchanged.

## Validation and next step

The finite-population first/second moments pass140 exact checks; the
exposure budget passes2520 exact path checks. The initial kappa rejection
is separately verified with rational arithmetic. The final frozen-parameter
upper-grid screen finishes in14.027s. Numerical runs are serial and file-free.

Local free space started at25.19GiB, then the20GiB guard stopped a later
run. Remaining runs used the configured Peach connection, after checking
63GiB free in /tmp and an idle account. Sources were streamed over pinned
stdin under Python-B, a30s timeout and256MiB address-space cap. No remote
files or packages were created. An overly strict absolute check on high odd
quadrature moments was replaced by a magnitude-scaled roundoff check; the
successful final runs include all three reported quadrature orders.

The subsequent directed-envelope replays cover the selected histogram and
the explicit regions above, including the local125-point boundary box.
Next extend coverage across the histogram simplex. Even a successful uniform
bound would close only this selected structural tail. Other target-row
counts, support profiles, and field-cancellation counting remain required
for the full minimum-distance proof. No native encoding benchmark, library
change, paper edit, commit, or promotion was performed.

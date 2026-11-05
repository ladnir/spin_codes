# Joint integral bounds over histogram regions

The [quadratic continuation](RANDOM_PAIRING_FLOW_QUADRATIC.md) obtains
directed local moment bounds, but needs uniform histogram coverage. Scalar
envelopes fail even at a difficult point that the matrix envelope certifies.
This continuation retains matrix directions and population dependence through
Gaussian integration. It does not change the encoder or close its distance
proof. Parameters remain the selected1974-row structural-tail parameters.

## Integrate control exponents before taking a maximum

Fix a region of the auxiliary five-state histogram h, with sum h=512.
Eliminate one coordinate and write x for the four remaining coordinates.
On Gaussian annulus S_r, the existing verified envelope has exponent

    psi_r(x,z)=C_r(x)+l_r(x) dot z-z^T P_r(x) z/2.

The Gaussian density is included in this exponent; the normalization is
(2 pi)^(-2). The functions C_r,l_r,P_r are polynomial in x of degree at
most two. They retain the matching mean and the affine exposure majorant.
Their coefficients depend on the fixed reference center and the annulus.

Expand all these polynomials in the same tensor Bernstein basis on the
population box. Its weights B_alpha(x) are nonnegative and sum to one.
Use one common degree for all coefficients and all annuli. Then

    psi_r(x,z)=sum_alpha B_alpha(x) psi_(r,alpha)(z).

Convexity of the exponential gives

    integral_(S_r) exp(psi_r(x,z)) dz
      <= sum_alpha B_alpha(x) integral_(S_r) exp(psi_(r,alpha)(z)) dz.

Each control exponent has a constant, linear vector and precision matrix.
The verifier checks that every control precision matrix is positive
definite. Completing the square gives its full log integral

    C_alpha-log det(P_alpha)/2
      +l_alpha^T P_alpha^(-1) l_alpha/2.

Multiply by the previously derived Gaussian annulus bound, using mean
P_alpha^(-1)l_alpha and covariance P_alpha^(-1). All operations use
directed intervals. Control coefficients need not correspond to actual
histograms: their role is an exact convex representation of the exponent.

The population weights B_alpha(x) are the same in every annulus. Therefore
sum each control's annulus bounds before maximizing over controls:

    sum_r I_r(x) <= max_alpha sum_r J_(r,alpha).

The weaker sum_r max_alpha J_(r,alpha) permits a different population
control in every annulus. That loss alone prevents one tested region from
passing. Degree elevation also tightens the Bernstein envelope without
changing the polynomial or introducing new encoder randomness.

## Retain the population polytope when a box includes impossible points

Eliminating one population and bounding the others independently can add
points outside the original five-coordinate bounds. A second control
representation avoids this relaxation. Enumerate the vertices v_i of
the five-coordinate box intersected with sum h=512. At every vertex,
at least four coordinates are at bounds; the helper enumerates these
possibilities and deduplicates them exactly.

For each h in this polytope, fix any convex representation h=sum_i lambda_i v_i.
Use that same representation in every annulus. For a quadratic polynomial q,
define its symmetric polar coefficients by

    q[v_i,v_j]=2q((v_i+v_j)/2)-(q(v_i)+q(v_j))/2.

Then

    q(h)=sum_i lambda_i^2 q[v_i,v_i]
          +sum_(i<j) 2 lambda_i lambda_j q[v_i,v_j].

The weights are nonnegative and sum to one. Apply the same integral argument
to these common polar coefficients of C,l,P. No triangulation or choice of
one preferred convex representation is required. This representation can
be sharper geometrically, but its controls can still overestimate curvature.

## The infinite tail must remain covered

Local matrix envelopes are valid only on their specified annuli. The region
verifier now also uses the preceding globally valid coupled radial tail,
rather than rejecting a region when the older scalar tail coefficient is
nonpositive. A failed old tail formula was an inconclusive bound, not a
successful certificate.

For the control representation, retain the affine exposure polynomials
B0(h) and B1(h), and a uniform bound c2 on its quadratic coefficient.
Choose a fixed0<epsilon<1/2-kappa-c2. Young's inequality bounds the linear
radial correction by epsilon||z||^2+B1(h)^2/(4 epsilon). The resulting
global exponent is still quadratic in h. Its precision is the constant
matrix2(1/2-kappa-c2-epsilon)I. Integrate each control only over ||z||>=16.
Add that tail to its own annulus contributions before taking the maximum.
The verifier also retains the independent uniform tail bound and takes the
smaller valid complete result. No local envelope is extrapolated to infinity.

## Results and failed controls

The allowed log moment is greater than7.679168691101. Define B_r by

    h0 in[185-r,185+r], h1 in[188-r,188+r],
    h2 in[83-r,83+r], h4=0, h3=512-h0-h1-h2.

Selected directed80-digit results are:

| Region | Histograms | Method | Log upper bound | Outcome |
|---|---:|---|---:|---|
|B_2|125|One common lower precision matrix|7.680537|Inconclusive|
|B_2|125|Degree2 integral controls, separate shell maxima|7.677500|Pass|
|B_4|729|Degree4 controls, separate shell maxima, before tail improvement|7.681221|Inconclusive|
|B_4|729|Degree4 controls, joint shell weights and improved tail|7.677137|Pass|
|88<=h_i<=116|421,750|Common lower precision matrix|7.584312|Pass|
|88<=h_i<=124|871,720|Polytope controls, joint shells and tail|7.485025|Pass|
|80<=h_i<=136|4,884,615|Same polytope method|8.618565|Inconclusive|
|64<=h_i<=128|8,250,385|Same polytope method|10.172558|Inconclusive|

All bands also impose sum h=512. Bounds are rounded upward. Earlier passing
bounds remain valid after the later tail refinement. These regions overlap;
do not sum their sizes without a disjoint partition. The two broad failures
are not code counterexamples. Their local integral bounds are also above
the allowance, and their tails grow loose because a single reference center
serves a broad population region. Fixing only the tail is insufficient.

The final B_4 replay takes11.936s. Its joint finite-annulus bound is below
7.677136630, compared with7.680695308 from separate shell maxima. The
independent coupled tail gives the final7.677136805 bound. The larger
interior passing band takes9.199s in its first replay.

The cover helper now passes778 exact domain, polynomial and Gaussian
checks. New checks reconstruct quadratic polynomials from barycentric
polar controls and degree-elevated tensor controls, using rational inputs.
Analytical Gaussian controls check constant and varying precision examples.
All numerical runs are serial, bounded, and file-free. Prior helpers and
proofs are preserved; no native encoding timing or production edit is made.

## Global curvature from the pair outcomes

The medium population band first failed only because of its infinite-tail
bound. Degree3 polar controls put its finite-annulus contribution below
7.672986, but the isotropic tail still prevented certification. A global
matrix curvature bound removes that particular loss.

Fix a pair type(a,b). Its Gaussian score is the logarithm of a weighted
sum of three exponentials. Let v_1,...,v_r be its distinct Gaussian slope
vectors, with r<=3. Duplicate output vectors are merged for this argument;
their positive exponential weights induce some categorical probabilities p.
The Hessian is the covariance of these vectors under p, at every Gaussian
point. The distribution p changes with the point; it is not fixed or uniform.

For every categorical probability vector p,

    diag(p)-p p^T <= I/2.

Indeed, the matrix on the left has row absolute sum2p_i(1-p_i)<=1/2.
Gershgorin's bound proves the claimed spectral inequality. If V has columns
v_i-bar(v), where bar(v) is their unweighted mean, covariance is unchanged
by this translation. Consequently the pair Hessian is globally bounded by

    H_(a,b) := V V^T/2.

This retains its actual directions instead of replacing it by a multiple
of the identity. It is zero for one distinct outcome. For two distinct
outcomes it is the familiar one-quarter outer product of their difference.
The bound does not assume independent coordinates or a fresh uniform state.

The Hessian of the mean matching score is bounded by

    H(h)=sum_(a<=b) E[e_(a,b) | h] H_(a,b).

The exact expected matching counts are quadratic polynomials in h. Thus
the tail precision matrix

    P_tail(h)=(1-2c2-2epsilon)I-H(h)

still has degree at most two. Retain the earlier Young bound on the
exposure term, expand all tail coefficients in the same control basis,
and verify positive definiteness at every control. Sum each control's tail
and annulus contributions before maximizing. This is an additional valid
complete bound, not an entrywise minimum of matrices or a per-control
minimum between unrelated envelopes. The older complete bounds remain
available and are compared only after their full contributions are summed.

For the band80<=h_i<=136, sum h=512, this reduces the maximum tail-control
log bound from8.237293 to-1.180483. Combined with degree3 controls, the
complete one-update moment bound is below7.672987 and passes. The band
contains4,884,615 histograms. This is still the selected1974-row moment,
not a minimum-distance certificate.

### Higher-degree polar controls and their limit in this pass

Polar controls now support degrees3 and4. For a multi-index alpha with
sum alpha=D, the elevated coefficient is

    [sum_i alpha_i(alpha_i-1) q[v_i,v_i]
      +sum_(i<j) 2 alpha_i alpha_j q[v_i,v_j]] / [D(D-1)].

This is the degree-D multinomial Bernstein representation of the same
quadratic polynomial in the convex weights. Exact rational reconstruction
tests cover the new degrees. The implementation rejects more than2000
polar controls, and checks a wall deadline during control integration.

The larger band64<=h_i<=128 remains inconclusive. With the new tail bound,
degree3 gives7.754862; degree4 gives7.712463. The latter bounded replay
takes46.049s under a75s limit. Its finite-annulus bound alone is7.712370,
so the infinite tail is no longer the relevant loss there.

A targeted alternate exposure order(4,3,2,0,1) reduces the degree3 result
to7.715887, but does not close it. Exposure order is a choice in the matching
proof, not a change to the code's random pairing. This test motivates
population-dependent partitions; it does not justify an unbounded search
over polynomial degrees or a claim that the actual moment fails.

## An explicitly accounted partial partition

The cover helper now subtracts certified integer boxes from the pending
domain using disjoint coordinate slabs. Exact small-domain tests check
point membership and the count identity, including empty intersections.
Its fixed seed partition contains the five previously certified caps,
the band80..136, and the tight boundary region B_4. These seven regions
are pairwise disjoint.

The full80-digit replay rechecks every region, rather than trusting stored
passing labels. It certifies241,866,501 histograms and leaves2,677,752,444
in25 explicit pending regions. Their counts sum exactly to2,919,618,945.
The band bound is7.672986463; the boundary bound is7.677136631. The replay
takes37.630s under a55s deadline. It reports both complete_histogram_cover
and whole_code_certificate as false.

The helper passes1040 checks after adding exact degree-elevation,
categorical-covariance, translated outcome-matrix, and partition tests.
Source specifies the finite partition; no search tree or raw timing dump
is saved. All numerical runs are serial, file-free and above the20GiB disk
guard. No construction, native kernel, paper or production code is changed.

Reproduce the combined partial certificate with:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --verify-partition --precision80 --deadline-seconds=55
```

## Next discriminating step

Refine the25 pending regions, using the matrix method near the tight
boundary and cheaper radial bounds away from it. Split according to
population differences and exposure order where the broad controls lose
slack; do not merely keep increasing degree. The full histogram cover is
still open. Other target-row counts, support profiles and field-cancellation
counting also remain before any whole-code distance certificate.

## Retuning the potential

The first broad refinement of the25 pending regions tested21 boxes in a
35.835s bounded run. It certified none. This does not refute the moment
inequality: the broad radial bounds can be far above the pointwise values.
Before investing in finer partitions, we revisited the proof potential.

Keep the same construction, auxiliary flow, target count1974, markers
theta=.1 and eta=.2, and quadratic coefficient kappa=.08. Only the linear
weights in phi(h)=u.h+kappa*||h||^2/512 change. The new candidate is

    u=(0,-79200914,-154166955,-195740597,-227550536)/10^9.

The original weights remain the verifier's default. Optional weights are
passed explicitly to the integrand and endpoint allowance. Each result
records the exact numerator tuple and denominator. Old and new regions
cannot simply be combined in one telescoping argument: that argument needs
one fixed potential, or a separately justified potential-switching bound.

### Faster proposal generation

For a pair type(a,b), let V be its random next-population vector under the
exponentially tilted three-outcome distribution. Each V has nonnegative
coordinates summing to two, and its diameter is at most sqrt(8). Thus
Cov(V)<=2I. With e_ab denoting the expected number of pairs of this type,
the map T from a proposed center to the tilted next mean satisfies

    D T = (2*kappa/512) sum_ab e_ab Cov(V) <= 2*kappa I.

Here sum_ab e_ab=256. Direct iteration is therefore contractive at the
fixed kappa=.08; damping is unnecessary. Pair counts are computed once.
The iteration needs only first moments, not the quadruple sums required
for the final second moment. Any resulting center remains only a proposal;
the certificate rechecks its rounded value using interval arithmetic.

The preserved damped solver and direct solver agree within3.32e-13 in
every pair-outcome probability on72 test histograms. Their scores differ
by less than1e-8. Timings were2.998s versus0.166s over that test set.
The140 exact finite-population moment tests also pass. These are proof
helper timings, not encoder measurements.

### Search and directed results

A bounded coordinate search used variational lower bounds, not certificates,
on a coarse simplex grid and neighbors of its largest values. On its final
1150 sampled points, the largest lower-bound-minus-allowance changes from
-0.054241 to-0.647274. The candidate sampled maximizer is
(189,168,115,40,0). Order-five Gaussian quadrature of the upper integrand
gives7.058405 there; numerical quadrature is not a directed upper bound.

The independent80-digit analytical verifier then gives log M<=7.083655
at that point, against allowance>7.679295226742. This local slack does not
establish a uniform moment bound or a whole-code margin.

Directed region checks under the same new potential give:

| Region | Integer histograms | Log upper bound | Method |
| --- | ---: | ---: | --- |
|64<=h_i<=128, sum h=512|8,250,385|7.659564|degree3 polar controls|
|h0=173..205,h1=152..184,h2=99..131,h4=0..16; h3 determined|593,589|7.464759|degree2 polar controls|
|Five previously used caps|236,981,157|7.016544 maximum|radial bounds|

The caps and the two new regions are pairwise disjoint, so these checks
cover245,825,131 histograms under the new potential. This is still partial
coverage of2,919,618,945 histograms. The old partition remains separately
available. The larger band fails with degree2 controls at8.020547; degree3
closes it in19.451s, below a55s deadline. The two-region probe takes35.808s;
the five-cap replay takes4.629s. All runs are serial and file-free.

The original-weight partial partition also replays successfully in37.006s.
The original point verifier passes5018 arithmetic/polynomial checks at80
digits and gives7.672286057. Six further interval checks confirm that adding
the same constant to every linear weight changes neither the core/gradient
nor the endpoint allowance. This tests consistent parameter propagation.
The frozen pair-defect diagnosis helper is unchanged.

Reproduce the proposal regression and candidate search with:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_quadratic.py --proposal-regression
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_quadratic.py --retune-potential
```

Reproduce the new directed checks with:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-probe --precision80 --deadline-seconds=55
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-probe --only-retuned-band --degree3 --precision80 --deadline-seconds=55
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-caps --precision80 --deadline-seconds=25
```

The next test is a broader upper-integrand search, including each boundary
face, to look for displaced maxima. Then construct coverage using one
fixed potential. Do not infer global closure from the lower-bound search
or spend unbounded effort subdividing the old tight envelope.

### Broader numerical check of the fixed candidate

The subsequent upper-integrand search examines863 distinct histograms:
the coarse simplex grid, followed by local searches in the full domain and
on each face h_i=0. The largest sampled order3 value is7.059943 at
(191,167,113,41,0). At that point, orders5 and7 give7.060040 and7.060025.
The allowance is7.679295. The sampled face maxima at order7 are6.850945
(h0=0),6.589195(h1=0),5.398263(h2=0),6.927041(h3=0), and7.060025(h4=0).
Agreement of quadrature orders is diagnostic, not a certified error bound.
The16.580s run found no displaced numerical peak; it does not prove absence.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_quadratic.py --retuned-upper-search
```

The next proposed refinement retains the actual population polytope in
the cheap radial bound, not just in the more expensive matrix controls.
One can use the same quadratic polar weights for the core and gradient,
then bound their integrated energy by convexity. This is a proposed next
step, not an implemented result or additional certified coverage.

## Population-aware scalar bounds

The next pass implemented that refinement, keeping the retuned potential
and the construction fixed. It sharpens some broad upper bounds but does
not replace the directional matrix method everywhere.

Write a scalar Gaussian envelope as

    (2*pi)^(-2) integral exp(C(h)+g(h).z-a||z||^2+b||z||) dz,

where a>0, b>=0, and the Gaussian contrast vector z has four coordinates.
The core C and each component of g are quadratic in the population h.
Choose a fixed epsilon in(0,a) and let gamma=a-epsilon. Young's inequality
and Gaussian integration give the logarithmic upper bound

    C(h)+||g(h)||^2/(4*gamma)+b^2/(4*epsilon)-log(4*gamma^2).

Previously, the polynomial part was bounded on a coordinate box. That box
can contain impossible populations because the eliminated fifth population
must also satisfy its bounds and sum h=512. Instead, express h as a convex
combination of the exact vertices of the clipped population box. The same
nonnegative polar weights express C and g. Convexity of the squared norm
then gives an upper bound from their common controls. The old complete
bound remains available; compare complete upper bounds, not matrix entries.

### Keeping signed cross terms at degree four

Squaring each gradient control separately introduces another convexity
loss. The optional degree-four version represents the squared polynomial
exactly. If b_alpha and c_beta are degree-two Bernstein coefficients in
barycentric coordinates, the degree-four product coefficient at nu is

    sum_(alpha+beta=nu) multinomial(2,alpha)*multinomial(2,beta)
        /multinomial(4,nu) * b_alpha*c_beta.

Apply this to each squared gradient component and sum. Elevate C to degree
four and combine before maximizing. Signed cross terms are retained; a
coefficient of a nonnegative polynomial need not itself be nonnegative.
The coefficient count is capped at2000. Larger cases are reported as
inconclusive rather than silently weakening or skipping their obligations.

Exact rational tests independently expand the polynomials in ordinary
barycentric monomials and check every product coefficient. Other tests
reconstruct values at rational convex combinations and check a Gaussian
integral on the segment x+y=1, where the enclosing rectangle loses slack.

### What the directed tests say

All values below use80-digit directed arithmetic and the same retuned
weights. The allowance is greater than7.679295226742.

- On the505,215,280-histogram low-h0 region, the scalar bound improves
  from49.451373 to31.358239, but still fails. Degree four exceeds the guard.
- On the64..128 band, the scalar bound improves from15.213360 to9.802749.
  Exact degree-four products improve this only to9.765471. Neither closes.
- Splitting Gaussian radii into shells improves the band's degree-two
  bound to9.574624. The directional matrix bound is still needed there.
- On the593,589-histogram boundary box, the shell-based scalar bound passes
  at7.549269, versus7.979229 without shells. Population controls do not
  materially improve this box; the shell split is responsible for closure.

The degree-two comparison takes13.774s, the guarded degree-four comparison
7.953s, and the shell comparison31.359s. These are proof-helper runs, not
encoder timings. They recertify an existing region more cheaply; they do
not establish a new distance certificate or complete histogram coverage.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-radial-controls --precision80 --deadline-seconds=45
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-radial-controls --degree4 --precision80 --deadline-seconds=55
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-radial-controls --radial-shells --precision80 --deadline-seconds=55
```

### A complete-domain subdivision pilot

The pilot starts from the entire histogram domain. It rechecks the five
caps under the retuned weights, then subdivides their exact disjoint
complement. It uses longest-coordinate splits and visits the larger child
first. Every successful region contributes its exact integer count; every
unfinished region remains in the pending sum. This does not import any
old-weight certificate or assume that an unvisited region passes.

The first short run checks11 regions after the caps and adds no coverage.
A24-node run then closes the following five disjoint regions. Bounds are
listed in coordinate order(h0,h1,h2,h3,h4), always with sum h=512.

| Bounds | Histograms | Directed log upper bound |
| --- | ---: | ---: |
|(0..43,44..87,80..119,40..79,240..319)|2,616,996|7.601821|
|(0..43,44..87,80..119,0..39,240..319)|2,177,804|7.088908|
|(0..43,44..87,120..159,0..39,240..319)|2,616,996|7.222525|
|(0..43,44..87,120..159,40..79,240..308)|915,920|7.411315|
|(0..43,0..43,80..159,0..79,240..319)|7,044,675|7.672885|

The new regions contain15,372,391 histograms. Including the caps rechecked
in this run,252,353,548 pass and2,667,265,397 remain pending in15 regions.
There are no failed singleton regions in this pilot. The exact counts sum
to2,919,618,945. The run takes83.655s under a100s deadline. Its source and
compact results are retained; no search tree is written.

These five regions are also disjoint from the previously retuned band and
boundary box, but the pilot reports only the regions it actually rechecks.
The optional newer scalar machinery leaves the old default verifier intact.
Its original-point replay still passes5018 checks at80 digits. The expanded
coverage helper passes1060 tests, including the independent rational
degree-four coefficient comparison. The frozen diagnosis hash is unchanged.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-cover-pilot --precision80 --deadline-seconds=100
```

Next use cheap scalar bounds where they suffice and directional matrix
bounds for difficult leaves. Share repeated switch/envelope preparation
to make the coverage test practical. The scalar comparisons do not justify
blindly increasing degree. Full histogram coverage and the subsequent
outer/field-cancellation obligations remain open.

## Triangulated population controls

The latest pass tests whether difficult regional bounds reflect difficult
populations or loss in their polynomial representation. The encoder and
retuned potential remain fixed. These are still one-update moment bounds
for the selected target-row count, not whole-code distance certificates.

### Retain only controls needed by a geometric cover

The rectangular representation uses few controls but includes infeasible
populations. The previous polytope representation instead includes every
pair of vertices, even when a smaller collection suffices to represent
every feasible population.

Triangulate the clipped population polytope. For a point h in one simplex,
write h=sum_i lambda_i v_i, with nonnegative weights summing to one. Every
quadratic polynomial p then has its exact polar representation

    p(h) = sum_i lambda_i^2 b_ii + 2 sum_(i<j) lambda_i lambda_j b_ij.

Only pairs of vertices in that simplex are needed. Take the union of these
pairs over the triangulation. Use this same collection for the constant,
linear and precision polynomials in every Gaussian shell and the tail.
For each h, its simplex and barycentric weights stay fixed across shells.
Thus the earlier Jensen argument and joint-shell summation remain valid.
There is no union-probability factor for the simplices: each h belongs to
at least one simplex, and the common bound applies there.

The implementation uses an exact pulling triangulation. Every face of a
coordinate box intersected with sum h=512 is defined by coordinate
equalities. Within each face, cone its least vertex over triangulated
facets not containing that vertex. This covers the face; using one vertex
order makes the facet triangulations compatible. All rank tests are rational.

Two optional refinements preserve the same domain. A strict interior
integer point can be coned over all triangulated boundary facets. Otherwise
the routine falls back to pulling. An integral edge midpoint can also
split every incident simplex into two. The latter removes a particular
long-edge control without dropping any populations. These controls currently
use degree two; the previous higher-degree method remains available.

### Directed evidence

All values use80-digit arithmetic and allowance greater than7.679295226742.

| Region or diagnostic | Controls | Log upper bound | Result |
| --- | ---: | ---: | --- |
| Merged region, rectangular basis |81|7.838282|Inconclusive|
| Merged region, pulling triangulation |157|6.654506|Pass|
| Extended region, pulling triangulation |157|6.713346|Pass|
| Balanced band64..128, pulling |45|7.895777|Inconclusive|
| Balanced band, interior-point cones |56|7.895777|Inconclusive|
| Balanced band, edge(0,9) split |55|7.823917|Inconclusive|
| Balanced band, then edge(0,7) split |63|7.764475|Inconclusive|

The merged region has bounds
((0,43),(44,87),(80,159),(0,79),(240,319)) and8,327,716 histograms.
It uses58 simplices and157 controls instead of378 unrestricted controls.
Its standalone check takes15.228s. The rectangular check takes9.964s but
fails; an unrestricted-polytope replay exceeds its28s budget. The earlier
successful unrestricted bound remains recorded separately.

The extended region replaces h0<=43 by h0<=63. It contains10,794,081
histograms, adding2,466,365 beyond the merged region. Its full command,
including the then-current5596 checks, takes20.398s. This does not claim
a newly replayed complete-domain accounting total.

The balanced-band failure is informative. After two edge splits, a worst
first-shell control joins (64,128,64,128,128) and (128,128,128,64,64).
Their midpoint (96,128,96,96,96) has a separately directed point bound
below6.999814. The regional bound remains above the allowance. This
demonstrates slack in this representation, not failure of that population.
It does not prove that every population in the remaining region passes.
The older degree-three certificate already covers this balanced band;
the experiment tests a cheaper representation, not a new obstruction.

### Validation and next step

The final suite passes7627 exact checks. New tests exhaust small integer
population domains, compare simplex volumes with an independent
inclusion-exclusion formula for a box slice, and reconstruct quadratic
polynomials at rational barycenters. They exercise pulling, interior-point
cones, edge bisection, and degenerate fallbacks. An initial zero-dimensional
star case failed and was corrected before the successful replay.

The earlier tensor and integer-polar arithmetic optimizations also pass
independent rational coefficient checks. The original-point verifier still
passes5018 checks and bounds its moment by exp(7.672286057). The old
seven-region partition replays successfully. The frozen diagnosis hash is
unchanged. The preceding hybrid run's final output was unavailable; its
missing process handle was confirmed, and no final coverage is inferred.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --simplex-probe=extended --precision80 --deadline-seconds=50
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --simplex-probe=balanced "--simplex-splits=0,9;0,7" --control-midpoint --precision80 --deadline-seconds=45
```

Next integrate triangulated controls into the explicitly accounted cover,
using a bounded number of targeted edge splits before ordinary region
subdivision. Retain the old degree-three bound where it is cheaper. Do not
spend an unbounded run re-solving the already certified balanced band.
Uniform coverage, other target-row counts, general histories, and the
field-cancellation composition remain open. All runs were serial and
file-free; none measures encoder performance.

## Accounted coverage with larger corner regions

These regions lie near vertices of the population simplex. The
construction, target-row count1974, and retuned potential remain fixed.

The first full-domain triangulated pilot rechecks the old five regions
and adds41,964,450 histograms. Its total is278,945,607, with2,640,673,338
pending. The earlier scalar-only pilot added15,372,391. This comparison
concerns proof coverage under bounded runs, not encoding performance.

### Larger simple regions

For coordinate i, the region h_i>=T_i contains binom(516-T_i,4) integer
histograms. Use thresholds (320,320,288,288,288). Any two thresholds sum
to more than512, so these five regions are disjoint. Their exact complement
is 0<=h0,h1<=319 and 0<=h2,h3,h4<=287, always with sum h=512.

All five regions pass directed verification. Higher degree below refers
to a polynomial representation in the proof; it changes no code operation.

| Coordinate | Threshold | Polynomial degree | Histograms | Log upper bound |
| --- | ---: | ---: | ---: | ---: |
|0|320|2|59,626,385|7.347898|
|1|320|2|59,626,385|5.789763|
|2|288|4|109,658,025|7.347218|
|3|288|3|109,658,025|6.727748|
|4|288|3|109,658,025|7.473030|

Together they cover448,226,845 histograms, compared with236,981,157 in
the earlier five regions. The rejected degree-two bounds for coordinates
2,3,4 were9.163016,8.270596,9.704433. Degree three still failed for
coordinate2 at7.834521. Retaining these failures distinguishes stronger
polynomial bounds from favorable point sampling.

### Spend matrix work after enough subdivision

The first pilot with the enlarged regions admitted matrix calculations
when every coordinate width was at most160. It added no interior coverage
within its roughly100-second work budget. Three matrix calculations
finished without passing; another timed out. All populations remained
accounted for, including the unfinished region.

Reducing the width trigger to96, with the same scalar cutoff14, closes
four interior regions containing47,851,881 histograms. The run reports
496,078,726 certified and2,423,540,219 pending. It takes107.466s including
tests, has one inconclusive deadline, and reports no failed singleton.
This supports avoiding expensive matrix calls on overly broad boxes;
it does not imply that those boxes contain counterexamples.

The source now preserves these passing regions as a compact partition,
not a search tree. Its pending-only mode explicitly separates previous
certificates not replayed in that call from newly checked coverage. A
one-node continuation closes another2,223,095 histograms below7.437052.
That region is included in the source checkpoint too.

The ten checkpoint regions are the five large regions above and these
five interior boxes, listed in coordinate order:

| Bounds | Histograms | Log upper bound |
| --- | ---: | ---: |
|(0..79,80..159,0..71,72..143,144..215)|15,958,619|7.497206|
|(0..79,80..159,0..71,72..143,216..287)|13,712,913|6.810244|
|(0..79,80..159,0..71,0..71,216..287)|15,958,619|6.620282|
|(0..79,80..159,0..71,0..71,144..215)|2,221,730|7.634662|
|(0..79,80..159,0..71,72..143,60..143)|2,223,095|7.437052|

Their union contains498,301,821 histograms. Exact subtraction accounts
for2,421,317,124 pending histograms without overlap or omission. The
checkpoint does not silently incorporate other previously certified boxes.
Every bound uses the same retuned potential; old-weight certificates stay
separate.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --verify-retuned-partition --precision80 --deadline-seconds=115
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-pending-cover --triangulated-cover --matrix-width=96 --cover-nodes=40 --precision80 --deadline-seconds=115
```

The updated exact suite passes7630 checks, including the new partition's
count identity. A fresh89.063s replay verifies all ten regions, reports
worst log upper bound7.634662, and leaves six pending regions with the
exact total above. The original-point control also replays5018 checks and
retains its7.672286057 bound. The controller reserves time for accounting, retains
inconclusive regions, and reports the largest checked log bound separately
from unproved regions. No partial maximum is a uniform bound or a whole-code
failure probability. Next work on the checkpoint's complement, using the
large mixed-population regions to choose further refinements. Uniform
coverage and subsequent outer/field-cancellation composition remain open.

## Continuation, precision, and a sharp global sensitivity bound

The next pending-only run certifies78,927,954 additional histograms at
80-digit precision. It checks22 regions in100.002s of search time, with
one unfinished deadline and no failed singleton. The previous498,301,821
are explicitly reported as not replayed during that search.

The source checkpoint now contains19 disjoint regions and577,229,775
histograms. It accounts for2,342,389,170 unproved histograms in its exact
complement. The added bounds and coordinates are retained in
`random_pairing_flow_cover.py:retuned_partition`, rather than a raw search
tree. They use the same retuned potential as the old regions.

### Optional40-digit directed arithmetic

Added an opt-in40-digit mode; the existing default and80-digit mode remain.
It uses the same outward-rounded Decimal operations. At40 digits, the7630
coverage checks and5018 original arithmetic/polynomial/Gaussian checks pass.
Lower precision is not treated as a reason to accept a failed inequality.

A serial A/B check on the coordinate2 corner region, threshold288 and
degree4, takes9.595s at80 digits and5.835s at40 digits. Both pass below
7.347218, and their reported upper bounds differ by less than10^-25.
This is one proof-tool timing comparison, not an encoder measurement or
a universal speedup claim.

A fresh40-digit replay verifies all19 checkpoint regions in104.881s.
Its worst log upper bound remains below7.634662, versus allowance above
7.679295. This is still partial coverage of the selected one-update moment,
not a uniform histogram theorem or a whole-code distance certificate.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --verify-retuned-partition --precision40 --deadline-seconds=115
python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --retuned-pending-cover --triangulated-cover --matrix-width=96 --precision40 --deadline-seconds=115
```

### Can common Gaussian variables reduce the global switch constant?

The existing global bound assigns a matching switch Lipschitz constant
4*sqrt(2)*beta. Bounding four terms separately suggests possible slack,
because all terms depend on the same Gaussian variables. We tested this
using the actual auxiliary pair transitions. A limiting example shows that
this constant is sharp over unrestricted Gaussian arguments.

Index the five auxiliary colors by0..4. For positive variables x_i, write

    F_ab(x) = sum_j P_j exp(theta*loss_abj + eta*j) x^(o_abj).

All coefficients are positive. The exponent vector o_abj is the outgoing
color population of the pair. The Gaussian-dependent score is log F_ab
after substituting x_i=exp(u_i+2*kappa*t_i/n+beta*Z_i).

Consider the switch

    S(x) = log F_22(x) - 2 log F_24(x) + log F_44(x).

Its actual monomial supports give positive coefficients a,...,g such that

    F_22 = a*x1^2 + b*x3^2 + c*x4^2,
    F_24 = d*x1*x3 + e*x4^2,
    F_44 = f*x3^2 + g*x4^2.

Set x=(1,epsilon^3,1,epsilon^-1,1) and let epsilon decrease to zero.
The unique leading terms are b*x3^2, e*x4^2, and f*x3^2. Therefore the
gradient of S with respect to log x converges to(0,0,0,4,-4), whose
squared Euclidean norm is32. Multiplication by beta restores the claimed
constant. The gradient has zero coordinate sum. Moreover, S is invariant
under common scaling of x, so removing the common Gaussian direction does
not exclude this path. The conclusion holds for every positive choice of
the coefficients, including the actual markers and label probabilities.

The helper `random_pairing_switch_lipschitz.py` checks these monomial
supports and unique leading exponents exactly. It also derives a useful
screening interface: if S=log(A/B), set D=AB and
N_i=(x_i partial_i A)B-A(x_i partial_i B). Coefficientwise nonnegativity
of c*D^2-sum_i N_i^2 suffices for a squared-gradient bound c. The FLOAT
screen covers55 distinct switches and returns maximum32; that screen is
not a directed certificate. Its190 exact rational identity tests and the
coefficient-independent limiting witness are separate evidence.

This is an obstruction to reducing the uniform global sensitivity
constant, not to the code or its distance. The limiting path moves into
an extreme Gaussian direction. A sharper proof must retain its Gaussian
mass, or derive smaller bounds on bounded radii, instead of replacing the
global constant by a smaller number. Next compare such radius-dependent
common-variable bounds against the existing local matrix envelopes while
continuing the accounted population cover. Construction and prior controls
remain unchanged; all runs were serial, bounded and file-free.

## Bounded sensitivity and the central-region bottleneck

The unrestricted sensitivity constant is sharp, but its limiting example
requires an extreme Gaussian tilt. We now bound the same switch on each
ball of Gaussian arguments. The encoder and auxiliary potential are unchanged.
This section uses the retuned weights and target-row count1974 only.

Let o_j be the five-coordinate output population of one auxiliary pair
when its tag count is j. Fix a nonnegative center t with sum512. Embed the
four-dimensional Gaussian Z isometrically into the zero-sum subspace of
R^5. With beta=sqrt(.16/512), its tilted outcome weights are

    w_j = P_j exp(.1*loss_j + .2*j
                  + <u + .16*t/512 + beta*Z, o_j>).

For outcomes j,k, put d=o_k-o_j. On ||Z||<=R, Cauchy--Schwarz gives

    log(w_k/w_j) in log(P_k/P_j) + .1*(loss_k-loss_j) + .2*(k-j)
                    + <u+.16*t/512,d> +/- beta*R*||d||.

For a bound uniform over all such centers, replace the center contribution
by [.16*min(d), .16*max(d)]. The identity
p_j=1/(1+sum_(k!=j) w_k/w_j) then gives probability intervals.
These intervals need not be independently attainable. Using them as a
superset is valid; no independence between pair scores is assumed.

Each coordinate of o_j is0,1,or2. Thus its tilted expectation equals
Pr[o_j[i]>=1]+Pr[o_j[i]>=2]. For each threshold, intersect the probability
interval obtained from its outcomes with the interval obtained from their
complement. This retains the probability-sum constraint. Adding the signed
expectation intervals for a matching switch bounds its gradient with
respect to log x. Squaring and summing encloses its squared norm.

The helper `random_pairing_switch_ball.py` checks all55 switches. Directed
upper bounds, rounded upward here, uniform over every center, are:

| Gaussian radius R | Squared gradient bound before the beta^2 factor |
| --- | --- |
| 0 | 5.286450 |
| 4 | 7.097414 |
| 8 | 9.200644 |
| 12 | 11.545935 |
| 16 | 14.073988 |

These replace32 only on the stated balls. For a switch S, integration
along the segment from0 to Z gives
|S(Z)-S(0)| <= beta*sqrt(B_R)*||Z||. The optional region verifier takes
the better of this bound and its existing local Taylor bound. It uses the
fixed center of the region, which can improve the table further.
It retains the old unrestricted bound on the infinite tail.

The new `shells32` option separates radii0,4,8,16,32 and then the infinite
tail. The finite-shell integration may integrate its positive majorant
beyond its domain to obtain an upper bound; it does not assert that the
majorant bounds the actual integrand there. Only the separate global
bound covers Gaussian arguments beyond32.

### Checks and discriminating cases

The helper passes2358 rational and interval checks. These include threshold
expectations under rational probability boxes, direct fixed-center
evaluation at radius zero, and the inherited arithmetic checks. An initial
test was stopped because converting a rounded-zero Decimal endpoint to a
Fraction constructed a million-digit denominator. The revised test compares
zero directly; it does not change interval arithmetic or any proof bound.

Using the same degree-two controls, the new scalar bound changes the
balanced region h_i in[64,128], sum(h)=512, from9.574624 to8.091269.
It still fails the7.679295 allowance. The stronger matrix bound for that
region is7.895776; its central-shell contribution alone is7.895776209.
Consequently further tail improvements cannot make this degree-two matrix
certificate pass. This diagnoses a central-envelope limitation, not a
construction counterexample. The previously established higher-degree
certificate for this region remains separate.

Two regions in the saved checkpoint's complement now have scalar proofs:

| Bounds on (h0,h1,h2,h3,h4), with sum512 | Histograms | Old scalar upper | New upper |
| --- | ---: | ---: | ---: |
| [0,79],[80,159],[180,233],[0,71],[144,252] | 5,623,240 | 9.511594 | 7.651482 |
| [0,79],[80,159],[180,233],[0,71],[0,71] | 4,423,190 | 7.541224 | 7.379158 |

The second region already passes the old scalar method, but had not been
included in the checkpoint. Both new proofs were replayed at80 digits;
their upper bounds agree with the40-digit runs to better than10^-30.
The neighboring26,133,720-histogram region remains unproved: its scalar
upper falls from20.878262 to8.836943, still above the allowance.

The source checkpoint now contains21 disjoint regions covering587,276,205
histograms, with2,332,342,740 pending. The exact partition accounting
includes both new regions. A fresh40-digit replay verifies all21 regions
in111.220seconds after7630 coverage checks; the worst upper is7.651482.
The original80-digit point control still passes5018 checks and retains
its7.672286056280 upper bound. Reproduce these serially with:

    python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_switch_ball.py --precision80
    python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --verify-retuned-partition --precision40 --deadline-seconds=120
    python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_envelope_verify.py --precision80

This remains a partial one-update certificate
at one target-row count, not a minimum-distance theorem. Next retain the
bounded-radius improvement as a cheap regional test, but focus sharper
central controls on regions whose central contribution already fails.
Do not tune the infinite-tail constant again or change the construction
on the basis of these failed upper bounds.

## Higher-degree controls on the existing triangulation

The balanced region's degree-two central bound fails even with a negligible
tail. This loss does not come from the encoder: its polynomial envelope
is represented by control exponents that need not describe actual
histograms. Degree elevation reduces that representation loss. Restricting
the controls to a triangulation also avoids unnecessary vertex mixtures.

Fix the existing triangulation of a population region. Let E contain
every unordered vertex pair belonging to a common simplex, including the
diagonal pairs. For each histogram h, choose one containing simplex and
fix its barycentric weights lambda. Use these same weights in all Gaussian
annuli, including the tail.

Write b_ij for a polynomial's degree-two polar controls. For degree d=3
or4 and a nonnegative integer multi-index alpha with sum d, define

    b_alpha = (sum_i alpha_i*(alpha_i-1)*b_ii
               +2*sum_(i<j) alpha_i*alpha_j*b_ij) / (d*(d-1)).

This is exact degree elevation for any polynomial of degree at most two.
Retain alpha only when every pair in its active vertex set belongs to E.
All multi-indices supported on the chosen simplex are retained. Every
omitted multi-index therefore has zero barycentric weight at h. The
retained weights still sum to one and represent the same polynomial.
Some retained cliques need not belong to a common simplex; their extra
controls only weaken the maximum bound. No histogram is excluded.

Apply this representation jointly to the exponent's constant, gradient,
and precision matrix. Jensen's inequality and Gaussian integration then
apply as before. Controls remain indexed identically across all annuli.
The implementation preserves the degree-two path and caches only the
small exact-integer elevation plans. Enumeration and retained-control
guards bound memory and work.

### The tail must use the same refined population controls

For the new mixed-region certificate, degree elevation fixes the central
bound but the old tail still fails. We combine the bounded switch
estimate on16<||Z||<=32 with the old global estimate beyond32. The
core covariance bound remains valid for every tilt. Both parts retain
the same population controls as the central shells, so their integrals
are added per control before maximizing. The finite-radius estimate is
never substituted into the infinite tail.

### Directed comparisons

All results below retain the retuned potential and target-row count1974.
The allowance exceeds7.679295; the displayed quantities are natural-log
moment bounds, not distance margins.

| Region and proof representation | Retained controls | Log upper |
| --- | ---: | ---: |
| Balanced h_i in[64,128], degree2 triangulated | 45 | 7.895777 |
| Same region, degree3 unrestricted | 220 | 7.659564 |
| Same region, degree3 restricted | 136 | 7.628991 |
| Mixed region below, degree3 restricted, old tail | 275 | 10.957919 |
| Mixed region, degree4 restricted, bounded annulus plus global tail | 670 | 7.597061 |

The mixed region is

    h0 in[0,79], h1 in[80,159], h2 in[180,233],
    h3 in[72,252], h4 in[0,180], sum h=512.

It contains26,133,720 integer histograms. Its degree-three central bound
alone is7.838647, so the tail refinement alone cannot close that bound.
Degree four reduces the central bound to7.597060107935. The old geometric
tail bound remains10.935476. With the bounded annulus, the16..32 contribution
is at most exp(-10.533933), and the beyond32 contribution is at most
exp(-37.469656). Retaining their common population controls yields the
final7.597060108312 upper bound. Thus both refinements matter here.

The mixed region would have4845 unrestricted degree-four controls, above
the existing2000-control guard. The restricted basis uses670. A40-digit
run completes in27.605seconds. Its80-digit replay, including the expanded
test suite, completes in57.299seconds and returns7.59706010831118765757.
These are proof-check timings, not encoder timings.

The additional5477 exact checks use an independent ordered-pair blossom
formula, barycentric reconstruction at interior and boundary points,
simplex-support coverage, and a negative test excluding a non-clique triple.
The complete coverage suite now passes13107 checks. The new region is
disjoint from the prior checkpoint and raises its total to613,409,925,
with2,306,209,020 histograms pending in five regions.

Replay supports explicit slices so each serial job remains bounded. A
slice reports only the histograms replayed in that call, separately from
the source checkpoint total. Reproduce the new case with:

    python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_flow_cover.py --verify-retuned-partition --partition-slice=21:22 --precision80 --deadline-seconds=60

Fresh40-digit replays of slices0:11 and11:21 verify514,260,440 and73,015,765
histograms in62.820and55.795seconds. Together with the80-digit new-region
replay, these disjoint slices recheck the entire22-region checkpoint.
The original80-digit point control passes5018 tests and retains its
7.672286056280 bound. The frozen fixed-pair diagnosis helper is unchanged.

This is still a partial one-update certificate. Neither the new control
basis nor the mixed-region check supplies the remaining target-row counts,
arbitrary-history coverage, or field-cancellation composition required for
whole-code distance. Next apply the same central/tail diagnosis to the
remaining mixed populations; do not treat every failed regional maximum
as evidence against the construction.

## The next intact region and a broader flow interface

The region

    h0 in[0,79], h1 in[80,159], h2 in[234,287],
    h3,h4 in[0,198], sum h=512

contains32,314,965 previously uncovered histograms. The bounded-annulus
method passes at degree four with upper6.818549. Degree three also passes,
at7.376753, using275 controls rather than670; the checkpoint retains this
cheaper certificate. Its80-digit replay, including13107 checks, completes
in32.747seconds and returns7.37675208771730657765. The source checkpoint
now covers645,724,890 histograms in23 disjoint regions, with2,273,894,055
pending. The prior22 certificates were not replayed in this continuation;
their source paths and verifier are unchanged.

The adjacent broad region
([0,79],[80,159],[72,179],[0,287],[0,287]), sum512, contains155,373,250
histograms. Its degree-four bound remains9.025866, essentially all from
the central contribution9.025865986442. The bounded-tail improvement has
already removed the tail as the limiting term. This is a failed envelope,
not a bad histogram or a construction counterexample.

Splitting that region at h3=143 produces two disjoint children. Both
degree-four checks hit their43-second deadlines; neither adds certified
coverage. The corrected serial probe ends after86.012seconds. Next inspect
the affine extension of the exposure-phase budget, since the broad region
crosses phase boundaries, before repeating the same axis-split approach.

Separately, the [varying-weight flow interface](RANDOM_PAIRING_FLOW.md#the-deterministic-interface-also-covers-varying-state-weights)
retains aggregate activity and the telescoping endpoint correction. It
shows how a completed setup-tail certificate would apply simultaneously
to many histories outside the constant-half-full family. The code,
auxiliary potential, and transition kernel are unchanged. The deterministic
extension does not complete the histogram cover or the field-count proof.

## Exposure cells and the remaining interpolation loss

The next diagnostic keeps the construction, potential, and target-row count1974
fixed. It studies the pending region

    h0 in[0,79], h1 in[80,159], h2 in[72,179],
    h3,h4 in[0,287], sum h=512.

Fix color priority(0,1,2,3,4). An exposure cell has final phase k when
the prefix before k is at most255 and the prefix through k is at least255.
Its greedy exposure quotas are h0,...,h_(k-1),255-sum_(i<k)h_i.
They are nonnegative on this cell. The actual squared switch ranges decrease
as colors disappear. Therefore these quotas upper-bound the exposure budget;
each range can then be replaced by its own matrix majorant. Unlike an affine
extension outside the cell, this argument needs no nested-matrix correction.

The helper enumerates vertices with exact rational elimination and uses all
box and prefix facet normals in a pulling triangulation. Closed real cells
overlap at boundaries. Integer accounting assigns each histogram to the first
prefix reaching255, giving105,823,530,45,894,905,and3,654,815 histograms in
phases2,3,and4. Their sum is155,373,250. These counts are not certificates.

Phase2 isolates the following losses. All listed bounds use degree-three
Gaussian controls, the same reference and center, and bounded outer tails.

| Proof representation | Directed log upper bound |
|---|---:|
| Force the nested-matrix correction on this cell | 9.646791297 |
| Exact phase quotas, pulling triangulation | 9.530488820 |
| Exact phase quotas, center-based triangulation | 9.530488820 |
| Retain integrated control weights, elevate to degree8 | 9.101206642 |
| Also subdivide two boundary edges | 9.080092587 |

Every entry fails the7.679295 allowance; no coverage is added. The40-digit
diagnostics use directed arithmetic. The last two checks take47.058and68.031
seconds, respectively. These are proof-check times, not encoding times.
These comparisons isolate degree-three machinery. The previously recorded
degree-four bound9.025866 on the entire parent box is still smaller and
remains preserved; this diagnostic does not replace it with a better bound.
A preceding one-edge check terminated at its55-second deadline without a
result. The later two-edge check is a separate, completed bounded probe.

The original worst control uses vertices
(0,80,175,0,257),(79,80,96,0,257),(79,104,72,257,0).
All lie on h0+h1+h2=255. Their centroid is
(158/3,88,343/3,257/3,514/3). Adding an interior center leaves this boundary
control unchanged. The nearby integer histogram(53,88,114,86,171) has a
directed upper bound6.988342 with the region's fixed Gaussian center and
priority. Thus that local test passes even without recentering. It does
not establish a maximum over the region or make a claim about code distance.
An80-digit replay agrees. Changing only the Gaussian center, with the same
priority and tail method, gives6.934139713 at80digits. The initial override
probe was rejected by the positive-definiteness check because its floating
matrix proposal still used the old center. Both proposal and directed check
now use the requested center; no result from the rejected probe is retained.

For the two-edge probe, the helper subdivides edges(0,15)and(10,15) at their
exact rational midpoints. An interval conversion encloses each rational;
no midpoint is rounded to a different histogram. Replacing each incident
simplex by its two cones preserves its entire real domain. Histogram counts
still refer to the original cell, not to the added proof vertices.

### Keep scalar control weights after integration

On one simplex let B_alpha be its degree-d barycentric Bernstein weights.
The common-control argument already proves M(h)<=sum_alpha B_alpha(h)J_alpha,
where J_alpha includes all central annuli, the16..32 annulus, and the far
tail for that same control. Previously the verifier immediately used max J.
Instead, elevate this scalar polynomial after integration. Its new coefficient
at multi-index beta, of total degree d+1, is

    J'_beta = sum_(i: beta_i>0) beta_i J_(beta-e_i)/(d+1).

These are convex combinations and preserve the polynomial exactly. Their
maximum is another valid upper bound. Repeat independently on every simplex
and take the largest bound. Controls not supported on any simplex need not
enter this maximum. All arithmetic is outward rounded. This operation does
not redo the costly Gaussian integrals or change the underlying code.

The scalar refinement reduces, but does not remove, the current failure.
After the two edge cuts, the worst control moves to a neighboring triangle
on the same face. The next discriminating test should split the full direction
that exchanges h3 and h4, rather than chase individual controls. If this is
still ineffective, inspect the Jensen loss of the Gaussian-control mixture
before adding further local subdivisions.

Reproduce each probe separately with `random_pairing_flow_phases.py` and
`--probe=nested`, `--probe=pulling`, `--probe=star`, `--probe=integrated`,
or `--probe=split`. The default deadline is95seconds, with an explicit upper
limit120. `--probe=point` checks the fixed-center integer histogram;
`--recenter` changes only its center, retaining the exposure priority.
Use `--precision80` for higher-precision replay. The default invocation checks
phase geometry and rational subdivisions; the separate
`random_pairing_flow_integrated_controls.py` checks scalar elevation exactly.

The final local suite passes5377 pulling-cell geometry checks,5381 star-cell
checks,451 native quota checks, and40 rational-subdivision checks in30.018s.
Scalar elevation passes1635 exact checks, including independent polynomial
reconstruction and a nonconstant example whose coefficient maximum is loose.
The unchanged original-u point control replays at80digits with5021 checks
and upper7.672286056280. The general coverage suite passes13107 checks.
The previous32,314,965-histogram high-middle region also replays at80digits,
returning the unchanged7.37675208771730657765 upper bound in36.818s including
the coverage tests. This is one checkpoint region, not a full replay.

The source checkpoint remains23 regions and645,724,890 covered histograms;
2,273,894,055 remain pending. No whole-code certificate or native random-pairing
timing is claimed. The current phase experiments do not change the old
checkpoint's verifier options.

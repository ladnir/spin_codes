# Short-length IMT bounds

These refinements address the finite BCH-256 calculation at K=2^16 and
K=2^18. They do not change either IMT map, the permutation distribution, or
the one-transvection update. Both covers are now complete: the accepted
K16 dense cover has 1,096 leaves and passed 512-bit replay. Its full union
has 41.8183267960 bits of margin. The historical point calculations below
are not themselves full certificates; see `PAPER_RESULTS.md` for the
current accepted inputs and `coupled_input_dense.py` for the last refinement.

The first useful retained point checks are in `FIXED_INPUT_POINTS_v1.json`.
Their 512-bit replay passed. At the previously failing density coordinate
21/128, the Q=218 contribution at K16 has margin 53.0405001895 bits;
Q=870 at K18 has margin 654.2693726381 bits. An isolated point does not
cover nearby densities or other occupancies.

## State activation

The seven-state representation is the one defined in
[`TRANSFER_ARGUMENT.md`](../asymmetric/TRANSFER_ARGUMENT.md): zero mass Z,
arbitrary nonzero mass D, and mass on five uniform image-weight shells.
Let a_v count the states in shell v. A coefficient C_v represents mass
C_v times the uniform probability measure on that shell.

Suppose the incoming state is zero and the epoch input is uniform among
weight-j vectors. Write k_j for the number with zero B-syndrome, c_j for an
upper bound on each syndrome fiber, and z=exp(-lambda). The output has
weight j. Every nonzero target has weighted mass at most

    c_j z^j / binom(128,j).

Thus a valid complete zero row is

    Z' = k_j z^j / binom(128,j),
    D' = 0,
    C'_v = a_v c_j z^j / binom(128,j).

`activation_density.py` uses this row only in the regime specified by its
integer threshold; elsewhere it retains the historical row. Each complete
row is valid. Taking entrywise minima between different representations
would not be justified.

For independent Bernoulli(r) input bits, put g=1-r+rz and
rho=|1-2rz/g|. Fourier inversion bounds each weighted B-syndrome mass by

    H = g^128 / 2^19 * (1 + sum_w b_w rho^w),

where b_w counts nonzero words of weight w in the image of B^T. The exact
zero-to-zero entry can be retained, with D'=0 and C'_v=a_v H.
`bernoulli_activation.py` evaluates the complete resulting transfer. It
compares complete moments with the historical transfer, not individual entries.

## Summing band assignments before taking a maximum

This refinement is implemented in `holder_split.py`. It improves one short
point but does not close either target by itself.

Let B=256 be the number of regions. Fix h all-one outer rows and d ordinary
active rows. Let R_j be the nonnegative region transfer conditioned on j
input ones. An ordinary band g has Bernoulli reference p_g and pointwise
row majorant Gamma_g. Define the positive linear operator

    (T_g V)(j) = Gamma_g^(1/B) [(1-p_g)V(j) + p_g V(j+1)].

For a band sequence sigma of length d, define C_sigma by applying its d
operators to R and evaluating at h. The existing row domination gives its
region contribution as e_0 C_sigma^B 1.

The objective is to sum these contributions over all band sequences.
Expand each matrix power into products along B-step state paths. For a fixed
path, Holder's inequality bounds the sum over sigma by the product of the
entrywise l_B norms of the B matrix entries. Summing paths then gives

    sum_sigma e_0 C_sigma^B 1 <= e_0 S_d(h)^B 1,

where S_d is the entrywise l_B norm over all length-d sequences.
Minkowski's inequality supplies a computable majorant:

    V_0(j) = R_j,
    V_(d+1)(j) =
      [sum_g Gamma_g ((1-p_g)V_d(j)+p_g V_d(j+1))^B]^(1/B).

All powers and roots in this recurrence are entrywise. Induction gives
S_d(j)<=V_d(j). This is an inequality between nonnegative matrices; it
does not assert independence between band assignments or messages.

Multiply the terminal bound by binom(L,Q)binom(Q,h), where Q=d+h, and
by the Chernoff factor exp(lambda H_cut). There is no extra factor 12^d:
the norm recurrence has already summed ordinary-band assignments. All bands
must remain in that sum, including bands discarded by a maximum-based hull.

## A composition-sensitive routing comparison

The old dense bound compares each region with iid input bits at cost
min(L+1,4 ceil(sqrt(L))). That uniform penalty is unnecessarily large when
many rows have nonconstant Bernoulli references.

Fix Q active rows. Their reference probabilities p_i belong to the finite
set P of ordinary-band references or equal one. The other L-Q probabilities
are zero. These are auxiliary reference inputs used to dominate outer rows,
not a new encoder setup. Let nu=(sum_i p_i)/Q and theta=Q nu/L.
After an independent region permutation, the input is uniform conditional
on its weight. It therefore suffices to bound

    Pr[sum_i X_i=j] / Pr[Bin(L,theta)=j],

where the X_i are independent Bernoulli(p_i) reference variables.

For 0<j<L, put y=j/L and x=y(1-theta)/(theta(1-y)). Exponential tilting
changes p_i to r_i=p_i x/(1-p_i+p_i x). The binomial reference tilts to
Bin(L,y). The likelihood-ratio identity separates the desired ratio into
a generating-function ratio and a ratio of tilted point masses.

Concavity of log gives

    product_i(1-p_i+p_i x) <= (1-nu+nu x)^Q.

For fixed theta and x, this upper bound increases with Q. For a box
Q in [q_lo,q_hi], nu in [nu_lo,nu_hi], set theta_lo=q_lo nu_lo/L and
theta_hi=q_hi nu_hi/L. An upper bound on the log generating-function ratio is

    q_hi log(1-y+(L/q_hi)(y-theta_*))
      + (L-q_hi)log(1-y) - L log(1-theta_*),

where theta_* is y clipped to [theta_lo,theta_hi]. Differentiation shows
that the expression increases up to theta=y and decreases afterwards.

To bound the tilted point mass, let V=sum_i r_i(1-r_i). Fourier inversion
and |1-r+re^(it)| <= exp(-2r(1-r)sin^2(t/2)) imply

    max_j Pr[sum_i X'_i=j] <= exp(-V) I_0(V).

Here I_0(V)=(1/pi) integral_0^pi exp(V cos(t)) dt. The function
exp(-V)I_0(V) decreases with V, directly from this integral representation.

For each p in P, the tilt x lies between the values computed from theta_hi
and theta_lo. The function p(1-p)x/(1-p+px)^2 is unimodal in x, so its
minimum over that interval is at an endpoint. Call this minimum v_p and
set v_1=0. The least mean variance at mean reference nu is the lower convex
envelope of the points (p,v_p), including (1,0). A two-support linear
program computes it exactly. This envelope decreases with nu: it is convex,
nonnegative, and zero at its right endpoint. Therefore

    V >= q_lo * envelope(nu_hi).

`density_comparison.py` evaluates these rational variance bounds, the Bessel
bound, and the binomial point mass with Arb. It maximizes over j=1..q_hi;
j>Q is impossible, and the endpoint ratios at 0 and L are at most one by
concavity. The result is uniform over the supplied Q/nu rectangle. Boxes
ending at nu=1 retain the historical factor instead.

The point factors are approximately 1.4156 at K16 and 1.4111 at K18,
compared with historical factors 92 and 184. The region comparison is paid
256 times. This explains the large change in the bound without any change
to the encoder. `test_density.py` checks the formula against exact
Poisson-binomial distributions, including all-one rows, and tests boxes
containing multiple Q values and means.

## Fixing the input tilt within a box

The old fixed-reference bound keeps an iid probability r fixed and changes
the input tilt x with the composition mean theta. Its proof consequently
pays for all possible compositions. At Q=218 that penalty is about 64.89
bits. It cannot simply be dropped from that bound.

`fixed_input.py` instead selects one positive x for the whole box. The
reference probability r(theta)=theta/(x+(1-x)theta) then varies over an
interval. Arb interval evaluation gives a positive transfer that dominates
every reference probability in that interval. The terminal moment therefore
has one uniform upper bound throughout the box.

For the outer sum, let Gamma'_g(x) be the input-tilted row majorant and
choose a real eta. Set

    S = sum_g Gamma'_g(x) exp(-eta p_g),

including the all-one band with p_g=1. The products of the normalized
weights Gamma'_g(x) exp(-eta p_g)/S sum to one over all labelled sequences.
Taking a uniform bound on the remaining mean-dependent scalar therefore
already sums all compositions. The existing scalar tangent bound controls

    Q log(S) + eta Q nu + 256 L log(1-theta+theta/x).

No composition-count factor is needed in this calculation. Exact sums of
binom(L,Q) count active positions across the Q interval. The implementation
takes the better of this entire bound and the old fixed-reference bound.
It does not remove a factor from the old expression while retaining its
composition-dependent choice of x.

## Coupling the normalizer and the moment

Bounding the scalar above separately from the moment can lose their shared
dependence on theta. `coupled_input_dense.py` keeps that dependence until
after it has formed an upper bound for their product.

Fix the output tilt lambda, input tilt x>0, and label tilt eta for a box.
For an N-bit input word u, let f(u) be its exponential output-weight moment,
averaged over the setup. Setup randomness is independent of the input law.
Write M(r)=E[f(U)] for independent Bernoulli(r) input bits. Then

```text
D(theta) = 1-theta+theta/x,
r(theta) = theta/(x+(1-x)theta),
D(theta)^N M(r(theta))
  = sum_u f(u) x^(-wt(u)) theta^wt(u) (1-theta)^(N-wt(u)).
```

The right side is a Bernoulli moment of a fixed nonnegative function.
Put z=logit(theta), h(z)=log(1+exp(z)), and c=eta L. The function

```text
F(z) = N log D(theta) + log M(r(theta)) + c theta
```

has second derivative at least -(N+abs(c))/4. To see this, write its first
two terms as the logarithm of a positive exponential sum minus N h(z).
The logarithm is convex, h''<=1/4, and abs(theta''(z))<=1/4.
Adding (N+abs(c)) z^2/8 therefore gives a convex function. Its secant bound
implies, for z between z0 and z1,

```text
F(z) <= max(F(z0),F(z1)) + (N+abs(c)) (z1-z0)^2/32.
```

Outward endpoint moment bounds can replace F(z0) and F(z1). For each integer
Q in the box, the calculation adds Q log S, log binom(L,Q), and lambda H.
It exponentiates and sums these Q bounds, then applies the uniform routing
density factor. The implementation uses this extra bound only for boxes
with at most eight integer Q values; larger boxes retain their earlier bound.

The high-band split uses the same endpoint calculation. Its ordinary branch
uses S_low and the restricted routing factor. Its exceptional branch uses
S_low+exp(-xi)S_high and adds xi. Both branches are formed directly and
summed, without subtracting separately rounded normalizers. Tests verify the
change-of-input identity exactly and check the curvature bound on positive
polynomials, including both signs of a rational c.

## Numerical status and reproduction

The l_256 split alone gave a best point margin of -408 bits at K16/Q218;
singleton bands and common reference-odds shifts did not improve it. The
routing comparison was the useful subsequent step. These negative margins
describe uninformative upper bounds, not lower bounds on actual failure.

```text
python -B workstreams/inner_design/finite_migration/fixed_input.py --output workstreams/inner_design/finite_migration/FIXED_INPUT_POINTS_v1.json --verify
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p "test_*.py"
```

Replay requires the locally retained producer record and authenticates its
source hashes. Generation uses the same command without `--verify` and a
fresh output path. The replay recomputes numerical bounds at 512 bits; it
is not an independent proof implementation. Tests supplement the derivation.

`short_dense.py` applies these bounds to a complete rectangular partition.
Unresolved leaves remain explicit. A passing point, a completed search
process, or a partial sparse range must not be reported as a full certificate.
No performance number is changed by this mathematical refinement.

## Total-budget acceptance

At K16, requiring 80 margin bits from each dense leaf is stronger than the
application needs. A retained bottleneck point currently has about 53 bits;
such a point cannot satisfy the old stopping rule even before interval loss.

`budget_dense.py` retains the same outward bound on each region. Let p_i be
the upward integer exponent stored for leaf i, so its contribution is at
most 2^p_i. The driver checks that the leaves cover every remaining
occupancy and density, then computes the exact rational sum

    U_dense = sum_i 2^p_i.

It refines the largest contributions until this sum falls below its chosen
budget or the bounded search ends. The initial K16 search allocates
U_dense < 2^-42. This is a search target, not the final distance assertion.
Higher-precision replay must recompute every retained leaf, including leaves
that do not individually attain 80 bits.

`verify_budget.py` separately authenticates the Q1, sparse, and dense
receipts. It reconstructs Q1 from the BCH inequalities, verifies every sparse
occupancy, and checks the complete dense partition and exact sum. It accepts
the full certificate only if

    U_Q1 + U_sparse + U_dense < 2^-40.

No number of failed or unexamined leaves is permitted. The driver records
an unmet budget explicitly and never labels its output as a full certificate.
The historical 80-bit-per-leaf producers and their accepted records are
unchanged. Tests include valid unequal leaf contributions, missing coverage,
altered sums, and equality at the final 2^-40 threshold.

## Separating all-one-rich band assignments

The density ratio above must allow many all-one rows. Such assignments
can have small reference variance, but their outer multiplicity is also small.
Applying the worst density ratio to every assignment loses this distinction.

Fix a box of active-row counts and reference means. Let h be the number
of all-one bands in an assignment, and put alpha=h/Q. The following split
bounds two subsets of the message sum, not two conditioned setup laws.
The encoder's setup randomness and each reference-input law are unchanged.
Use the fixed input tilt and transfer bound from the preceding section.

First consider assignments with alpha<=a, where a=1/32. At each Fourier
tilt, keep the ordinary variance lower bounds v_p defined above. Let c_p
be the fraction of active rows in ordinary band p. Their mean variance
is bounded below by the following finite linear program:

    minimize  sum_p c_p v_p
    subject to c_p>=0, 0<=alpha<=a,
               sum_p c_p + alpha = 1,
               nu_lo <= sum_p p c_p + alpha <= nu_hi.

`constant_density.py` enumerates its vertices with rational arithmetic.
On either alpha boundary, a vertex uses one ordinary reference, or two
ordinary references and a tight mean boundary. Away from both alpha
boundaries, a vertex uses one ordinary reference and a tight mean boundary.
The all-one reference contributes zero variance. Multiplying the minimum
by q_lo gives a uniform variance lower bound for the box.

Insert that lower bound into the same Bessel bound and generating-function
comparison. The resulting density factor applies to every assignment in
this first subset. The unrestricted outer band sum still upper-bounds its
restricted sum, so the scalar part of the fixed-input bound can be retained.
If nu_lo>a+(1-a)max(P), this subset is empty.

For the second subset, alpha>=a, retain the unrestricted density factor.
Choose a fixed xi<=0 and replace the band normalizer by

    S_xi = sum_g Gamma'_g(x) exp(-eta p_g - xi 1[g is all-one]).

For any labelled assignment, its product of tilted band costs equals its
product of normalized weights times

    S_xi^Q exp(eta Q nu + xi h).

The products of normalized weights sum to one over all assignments. Since
xi<=0 and h>=aQ, the additional scalar is at most xi*a*q_lo in log form.
The existing scalar tangent therefore bounds this subset after replacing S
by S_xi and adding that term. The implementation uses xi=-3335/16.
It adds the two complete subset bounds; it does not mix their transfer entries.
Assignments on the boundary can be counted twice, which is conservative.

At K16, Q=206, and density coordinate 49/256, optimizing the actual
fixed-input expression gives 31.4197780178 bits. Its density factor is
1.4786590073. The split lowers the first subset's factor to 1.4156227664,
giving 47.5100238751 bits. The second subset has 1958.19198917 bits.
Their sum retains 47.5100238751 bits to the displayed precision.
Both the unsplit and split point passed 512-bit replay.

The same derivation permits any fixed 0<a<1. A seven-cutoff comparison
in `M16_CONSTANT_THRESHOLDS_v1.json` gives 50.3999770231 bits at a=1/1024,
with 82.85 bits for the second subset. Its 512-bit replay passed.
At a=1/2048 the second subset becomes large enough to worsen the sum,
illustrating why the cutoff should be selected against the complete bound.
The first rectangular search retains a=1/32; these sharper point witnesses
do not change that search's acceptance rules or frozen producer.

The natural next simplification uses the integer count directly: split
h=0 from h>=1. The first subset uses the variance program with a=0.
For the second, xi*h<=xi replaces the fractional lower-count penalty.
A preliminary point evaluation gives about 50.49 bits for the first subset
and 323.07 bits for the second. This isolated point was not recorded in
a separate higher-precision replay receipt.

`zero_constant_dense.py` now implements this integer split uniformly over
rectangles. It also removes the all-one term from the h=0 scalar sum.
Its first bounded cover retained 702 leaves after 271 refinements, with
an unmet dense budget and margin -52.4435354395 bits. New center checks
found a Q=205, v=201/1024 point with only 42.5846958627 bits, so the
remaining loss is not solely interval width.

## Extending the split to high-reference bands

The same argument applies to a larger subset of bands, without changing
the encoder. Call a band high when its fixed Bernoulli reference exceeds
3/4; this includes the all-one band. Let h count high bands in an assignment.
Split the message sum into h=0 and h>=1. This is a partition of counted
messages, not conditioning on a setup event.

Write S_low and S_high for the sums of Gamma'_g(x) exp(-eta p_g) over
the two band classes. For h=0, remove high bands from the scalar normalizer
and from the variance program's allowed references. If the mean interval
lies above every low reference, this subset is empty. For h>=1, retain
the unrestricted routing factor, use the normalizer

    S_low + exp(-xi) S_high,

and add xi to its logarithmic bound. The inequality xi*h<=xi holds because
xi<=0 and h>=1. Both subsets retain the same complete positive transfer.
Their bounds are added before comparison with the unsplit bound.

The fixed cover implementation `high_band_dense.py` uses xi=-5123/64.
At the five retained difficult centers, `M16_HIGH_BAND_POINTS_v1.json`
gives margins 67.0411, 90.0779, 158.8202, 81.9355, and 68.6235 bits.
Every point passed 512-bit replay. Tests check the marked subset-counting
inequality exactly and compare the point and rectangle implementations.
These point bounds still require a complete dense cover for a K16 certificate.

An earlier cutoff of 1/2 was too aggressive: its high-band subset was not
small enough, and its complete bound was worse at four of these centers.
`M16_BAND_SPLIT_POINTS_v1.json` preserves that result and its replay.
Selecting the cutoff by the low-band bound alone would have missed this loss.

`combined_dense.py` retains the minimum of the complete all-one-only split
and the complete high-band split. Each bound covers the same box of
messages, so selecting their minimum is valid. Separate checker instances
keep their constrained variance programs distinct. Only their unrestricted
routing comparison cache is shared.

`fast_density.py` accelerates the same variance calculation. Before solving
the variance program for a weight j, it bounds the tilted point mass by one.
If this upper bound on the density ratio is already below the current
maximum, the more expensive variance calculation cannot increase that maximum.
The skipped term is still bounded, not omitted from the proof.

For the unrestricted comparison, the mean constraint already implies an
all-one fraction at most (nu_hi-p_min)/(1-p_min). Adding this redundant
constraint permits the same pruned evaluator. `fast_combined_dense.py`
uses these evaluations for the two complete subset bounds. Tests compare
both constrained and unrestricted factors against the original formulas.

`constant_dense.py` applies the same argument uniformly to rectangles and
uses the total-budget driver. A complete, replayed partition is still
required before reporting a K16 distance certificate. Tests check the
variance vertices and density ratios against exact toy distributions,
including rectangles, and compare the point and rectangle implementations.

## Endpoint bounds for a Bernoulli interval

`convex_input_dense.py` provides another uniform moment bound without
changing the IMT construction. Fix a nonnegative function f on N input
bits, and define M(p)=E[f(X)] for independent Bernoulli(p) input bits.
For the encoder calculation, f is the exponential output-weight moment,
averaged over setup randomness independently of the input.

With z=log(p/(1-p)), write

```text
log M(p) = log(sum_x f(x) exp(z wt(x))) - N log(1+exp(z)).
```

The first term is convex in z. The second derivative of h(z)=log(1+exp(z))
lies between zero and 1/4. On an interval [z0,z1], the difference between
the secant of h and h is therefore at most (z1-z0)^2/32. If U0 and U1
upper-bound log M at the endpoints, then throughout the interval

```text
log M(p) <= max(U0,U1) + N (z1-z0)^2 / 32.
```

The implementation evaluates both endpoints with the existing outward IMT
transfers and takes the minimum of this bound and the original interval
calculation. The input tilt remains fixed across the box. Positive-polynomial
tests exercise interior points and the zero-width limit. This backend has
not produced a completed K16 cover; it is not part of an accepted full
certificate yet.

## Joint-map audit: valid but negligible improvement here

`joint_overlap.py` separately audits the binary code generated by the rows
of A and B^T. Its rank is 38, so the two 19-dimensional images intersect
only at zero. The all-one word belongs to the B^T image. Reducing all
349633 supports of weight at most three modulo the joint generator basis
gives distinct syndromes. A collision would produce a nonzero joint word
of weight at most six; hence this joint code has minimum distance at least seven.

For q!=0, neither Aq+B^T z=0 nor Aq+B^T z=1 is possible. Complementing
within the joint code therefore also bounds this sum's weight by 121.
If v=wt(Aq), w=wt(B^T z), and h is their support overlap, the resulting
additional restrictions are ceil((v+w-121)/2)<=h<=floor((v+w-7)/2).
Intersect these with the usual support bounds before maximizing each
Fourier term. The zero-state row remains unchanged.

This audit barely improves the retained difficult point: about 0.00004 bits.
It is not used in the successful all-one-count point calculation or the
new cover search. The tests check its hypotheses, low-weight collision
rejection, and overlap inequalities independently on small examples.

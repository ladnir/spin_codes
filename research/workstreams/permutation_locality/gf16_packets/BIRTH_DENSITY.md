# Density of Newly Activated States

The existing three-coordinate bound sends every nonzero state created
from zero into the arbitrary-state component. This loses information
when several randomized packets contribute to its feedback. We keep
the rare low-occupancy inputs arbitrary and bound the remaining feedback
distribution pointwise.

This is a new bound for the same encoder. It applies to the dense proof's
iid packet comparison, not to the actual outer distribution directly.
The surrounding positive-measure domination and counting argument are
unchanged.

## The Weighted Feedback Measure

Fix the packet activity probability p and output weight z in (0,1].
The comparison input X consists of W=32 independent four-bit packets.
Each packet is zero with probability 1-p and is otherwise uniform among
the 15 nonzero values. Let J count its active packets. The fixed binary
feedback matrix C maps 128 input bits to s=19 state bits. Put S=2^s
and L=S-1.

When the entering state is zero, this epoch emits X and the next state
is CX. Fix an occupancy cutoff k. Define the nonnegative measure

    nu_k(b) = E[z^wt(X) 1_{CX=b, J>=k}].

We seek an upper bound on nu_k(b) for every b!=0. This measure includes
its occupancy probability; it is not a distribution conditioned on J>=k.

For a feedback character u, let r_i(u) be the binary weight of C_i^T u,
where C_i is the four-column restriction to packet i. Define

    a_r(z) = ((1+z)^(4-r) (1-z)^r - 1)/15.

This is the weighted character average of a uniform nonzero packet.
With a formal variable v, set

    P_u(v) = product_{i=1}^W (1 + v a_{r_i(u)}(z)).

Its coefficient of v^j sums the character products over j-subsets of
packet positions. Therefore the Fourier transform of nu_k is

    F_k(u) = product_i (1-p+p a_{r_i(u)}(z))
             - sum_{j=0}^{k-1} p^j (1-p)^(W-j) [v^j] P_u(v).

The expression is polynomial in p and includes the endpoints p=0 and
p=1. The factors and coefficients can be signed. The measure nu_k
itself is nonnegative.

## Centering the Fourier Sum

Fix b!=0. For any real constant c, character orthogonality gives

    nu_k(b) = (1/S) sum_u (-1)^(u dot b) (F_k(u)-c).

Consequently,

    nu_k(b) <= B_k := (1/S) sum_u |F_k(u)-c|.

The implementation chooses c near a weighted median of the character
values. Its exact optimality is irrelevant: every real c gives a valid
bound. Character profiles with equal counts of r_i values are grouped,
and their exact multiplicities are retained in the sum.

The total mass of nu_k supplies a second upper bound:

    H_k = sum_{j=k}^W binom(W,j) p^j (1-p)^(W-j) a_0(z)^j.

Use min(B_k,H_k) for the pointwise cap. In particular, the cap is exactly
zero when the retained occupancy range is empty.

## The New Zero-State Row

The state invariant uses zero mass Z, arbitrary nonzero mass M, and a
nonzero-state measure dominated pointwise by U/L. The high-occupancy
feedback measure fits the U coordinate with amplitude

    U_new = L min(B_k,H_k).

The low-occupancy nonzero feedback measure has total mass at most

    M_new = sum_{j=1}^{k-1} binom(W,j) p^j (1-p)^(W-j) a_0(z)^j.

The j=0 case cannot create a nonzero state. The displayed upper bound
may still count zero-feedback inputs at j>=1; that only enlarges its
bound on nonzero mass. The original zero-to-zero bound remains valid
for the entire input measure.

Thus replace the zero-source row by

    (original zero-to-zero bound, M_new, U_new).

Every nonzero-source row and every terminal weight stays unchanged.
No independence between feedback and emitted weight is assumed: the
factor z^wt(X) is inside nu_k throughout the calculation.

Each cutoff and the original row define a valid transfer matrix.
The search may choose the best candidate for its fixed comparison
probability and output tilt. That choice is a proof parameter, not a
change to setup or encoding. The selected matrix is then propagated
through the same sequence of epochs.

## Evaluation and Evidence

`birth_refresh.py` implements the `birth-refresh` dense kernel. Floating
arithmetic proposes a cutoff; outward evaluation recomputes its bound
with Arb and an exact dyadic centering constant. Newton identities
evaluate the truncated polynomial coefficients from their power sums.
They are algebraic identities; interval arithmetic encloses the signed
subtractions and divisions. No floating Fourier cap enters a certificate.

Tests compare coefficients against direct polynomial multiplication,
check every cutoff against exact small input enumerations, and include
zero and rank-one feedback matrices. Separate checks propagate the
state-domination invariant through repeated transitions.

At 8.5%, with row parity retained and a dense handoff at q=129, reoptimized
floating point checks are negative at every tested coordinate from
0.025 through 0.06. At coordinate 0.04, the proposed log2 bound is
-119.0705264457; fresh 256-bit evaluation gives
-119.0705264463. This is one point, not a domain cover. The complete dense
search and the matching sparse range still need independent replay
before any whole-code claim at 8.5%.

## Retaining the Expansion Weight

The next refinement keeps more information about low-occupancy births.
Let A be the fixed expansion map, and let E_l contain the nonzero states
b for which wt(Ab)=l. For the actual maps, l ranges over 48, 56, 64, 72,
and 80. An arbitrary measure supported on E_l need not be uniform.

The weighted mass entering E_l can be computed without enumerating
all 128-bit inputs. Let f_l be its indicator and hat(f_l) its unnormalized
Walsh transform. Character orthogonality gives

    E[z^wt(X) f_l(CX)] = (1/S) sum_u hat(f_l)(u) F_0(u).

For a restriction J<k, replace F_0 by the sum of its occupancy terms
with j<k. The integer Walsh coefficients are grouped by the feedback
character profiles already used above. This produces five exact signed
censuses, each with 9218 entries, before numerical evaluation.

The diagnostic adds a mass coordinate for each E_l. On empty input,
the lazy transvection branch preserves that class and emits weight l.
On nonempty input, that branch is bounded by arbitrary nonzero mass.
The refreshed branch retains the existing uniform-density bound. The
old arbitrary-state return bound also applies to each restricted class.
The hybrid keeps the Fourier density envelope for J>=k and uses these
classes only for J<k.

`birth_classes_probe.py` implements the integer census and floating
transfer proposals. Tests compare its census with direct input
enumeration and check eight transitions from every small initial state,
including zero and rank-one feedback maps. The separate `birth_classes.py`
now supplies outward evaluation for dense search and replay. No complete
whole-code certificate yet uses it.

At 9%, q>=129, and comparison coordinate 0.04, the optimized floating
score improves from approximately +2409.62 with `birth-refresh` to
+1985.49 with the hybrid. The gain is useful but does not close this
point. The pure class bound, without the high-occupancy density component,
is weaker than the hybrid here.

The same point's earlier diagnostic removed only the M and U return
entries, not the added class rows. It therefore did not measure the
effect of removing every return. `frontier_probe.py` now includes all
nonzero source rows and handles the variance-partition sum when measuring
shuffle loss. This correction affects diagnostics, not any certificate.

## Outward Class Transfer

The class coordinate F_l bounds the mass of an arbitrary measure
supported on E_l. It does not assert uniformity within E_l. Together
with Z, M, and U, these coordinates dominate the weighted state measure
by a sum of nonnegative measures of the stated types.

Fix an iid comparison activity p, output weight z in (0,1], and an
occupancy cutoff k. The verifier recomputes the exact weighted mass
of births into each E_l with J<k. Newton identities evaluate the
truncated occupancy polynomial. Signed integer Walsh coefficients and
Arb enclosures retain cancellations; only the final nonnegative class
mass is rounded upward. Births with J>=k use the preceding density
envelope. Another candidate retains all birth classes without a tail.

For a transition from F_l, define

    H_l = max_{a in E_l} E[z^wt(Aa+X)],
    h_l = Pr(X=0) z^l = (1-p)^W z^l.

The expansion histograms evaluate H_l as a finite maximum of polynomials,
with outward rounding. Put alpha=2^(-r), beta=1-alpha, and L=S-1. The
empty-input lazy branch stays in E_l and contributes alpha h_l to F_l.
The nonempty lazy branch contributes at most alpha(H_l-h_l) to M.
The refreshed branch contributes at most beta H_l to U. Returns to zero
are also included separately; this overlap only enlarges the bound.

Let R_{h,j}(z) be the existing rank/trimmed return bound for expansion
histogram h and j active packets. For class E_l, use

    R_l = sum_j Pr(J=j) max_{h of weight l} R_{h,j}(z).

The Z entry is at most alpha R_l + beta H_l/L. The verifier takes the
minimum of this bound, H_l, and the existing unrestricted-state Z entry.
All other source rows retain the previous three-coordinate bound.
No new uniformity assumption is made about states created by sparse input.

The search may select the original kernel, all birth classes, or a
class/density cutoff. Each choice defines a valid transfer matrix; the
outward evaluator reconstructs its entries rather than importing the
floating matrix. Small tests compare class masses with exact input
enumeration and check three transitions with rational arithmetic from
every state, including degenerate feedback maps. Separate floating
tests cover eight transitions and one, two, or three updates.

The `birth-classes` kernel is supported by saved dense witnesses and
the full assembler. With eight variance bins, selected mean intervals
at 8.7% and q>=97 pass at 256-bit precision. For example,
[0.029999,0.030001] has log2 bound -271.5351135632977, and
[0.039999,0.040001] has bound -226.2135159802855. The complete dense
cover now has 156 accepted cells and no unresolved cells. The full
assembler independently replayed this cover, giving dense log2 upper
-78.4855903358. Its sparse replay passed every occupancy except q=4,
where the coarse tilt grid did not meet the per-occupancy stopping
budget. A finer grid subsequently verifies q=4 even at 9%, with 89.2548
bits of margin. Fresh full assembly at 9% uses that grid; the partial
8.7% run is not a whole-code certificate.

The same class census also gives exact conditional birth masses at a
fixed packet occupancy. `occupancy_birth_classes.py` uses those masses
before the without-replacement regional placement polynomial. Its
derivation is in [the fixed-occupancy note](OCCUPANCY_BOUND.md). At 8.7%,
q=48, the resulting all-support bound has 175.1443 bits of margin, with
fresh 256- and 384-bit evaluations agreeing to the displayed precision.
The expanded tilt grid also lets the older rank bound close this
occupancy, at 63.0187 bits. Both experiments concern the same encoder.

At 9%, combining class-specific return bounds with a 32-bin variance
partition still leaves selected scores positive: about +657.10 at
x=0.04 and +562.54 at x=0.05. These diagnostics identify remaining
proof loss, not a counterexample to the construction.

Retuning the reference tilt changes that diagnosis. With base tilt
3/16 instead of 1/8, the 9% comparison has negative bounds on all five
sampled reference activities 0.15, 0.20, 0.25, 0.30, and 0.35. Narrow
mean intervals around activities 0.125, 0.15, 0.175, and 0.20 have also
passed outward checks at 256 and 384 bits. The weakest of these four
interval bounds is about -129.05. The complete dense search with the
new tilt closes in 116 cells. Fresh assembly is replaying this partition
and regenerating all sparse occupancies. This remains the same encoder:
only its proof changes, and a whole-code claim awaits the aggregate.

## Remaining Loss at 9.5%

The corrected diagnostic uses base tilt 3/16 and eight variance bins.
At reference activities 0.15 and 0.20, its log2 scores are +906.36 and
+993.78. Removing every nonzero-to-zero entry leaves +525.26 and +534.61.
Removing the shuffle comparison loss as well leaves +270.32 and +283.91.
These are deliberately optimistic counterfactuals at fixed witnesses,
not bounds on an altered encoder. Extra updates also leave these fixed
witnesses positive; retuning their witnesses remains a separate question.

The zero-state diagonal is close to the probability of entirely empty
input: 0.00551387 versus 0.00551322 at activity 0.15. At activity 0.20,
the values are 0.000792662 versus 0.000792282. This points to long empty
stretches in the comparison measure as another source of loss.

`activity_floor_probe.py` tests the constraint J>=38q, where J is total
packet activity and q counts active outer groups. It inserts a positive
activity mark and reoptimizes the output tilt. The tested marks do not
improve either point, even when the penalty is applied per active group
inside the outer count. The latter scan retains the original count duals,
so it does not rule out a jointly optimized bound. The regional-count
refinement now retains a separate shuffle density bound for each
packet count. Its whole-interval derivation is in
[the shuffle note](SHUFFLE_VARIANCE.md#retaining-the-regional-packet-count).
Three difficult intervals at 9.25% pass outward checks at 256 and 384
bits. The complete dense search subsequently covered q=97--2048 in
172 mean cells. Its independent 384-bit replay is running. The sparse
prefix remains separate: the initial q=96 search exhausted its budget,
and a finer output-tilt grid is being tried. No whole-code 9.25%
certificate is claimed.

The activity-floor diagnostic was also repeated with regional count
bounds and reoptimized outer-count duals. For a mark gamma>=0, the
constraint J>=38q gives

    1{J>=38q} <= exp(gamma J-38 gamma q).

The diagnostic multiplies the active-packet input weight by exp(gamma)
and the coefficient of every active comparison group by exp(-38 gamma).
It reoptimizes the mean, variance, and occupancy counting duals after
that change. The component labels preserve q, so the mark charges
each active group rather than only the lower bound q_min.

At 9.5%, reference activity 0.175, and eight variance parts, the
unmarked floating proposal is +135.5516. Every tested positive mark
from 1/256 through 1/8 makes it worse. Output-tilt scale 39/40 has
unmarked proposal +129.5785 and the same outcome for all tested marks.
This rules out that tested parameter route, not every possible use of
the support floor. No certificate uses this diagnostic.

## Density at a Fixed Packet Count

The regional bound conditions on the exact number of active packets.
Its class operator does not use the preceding iid density tail. To test
the resulting loss, `occupancy_birth_density_probe.py` bounds the weighted
feedback density separately for each occupancy j.

Retain the polynomial P_u(v) defined above. Conditional on J=j, its
normalized coefficient is the weighted feedback character moment

    F_j(u,z) = [v^j] P_u(v) / binom(W,j).

For every nonzero feedback target b and every real c, character
orthogonality therefore gives

    E[z^wt(X) 1{CX=b} | J=j]
        <= (1/S) sum_u |F_j(u,z)-c|.

The probe selects an exact dyadic centering constant and evaluates the
sum with outward arithmetic. The total moment a_0(z)^j is another cap.
Write B_j for the smaller cap. Replacing the nonzero birth-class entries
by U=L B_j gives a valid density envelope. The zero-state entry is
unchanged. A chosen occupancy cutoff selects where to use this replacement.
This changes a proof bound, not the encoder or its randomness.

Exact small-input tests check every nonzero feedback target, degenerate
feedback maps, and ordered products of mixed empty and active steps.
The first performance-independent screen uses three updates, activities
0.175 and 0.20, and output tilts 0.10, 0.12, and 0.14. Among the tested
cutoffs, the largest improvement in the iid moment is about 5.91 bits.
That screen does not evaluate the full regional bound. The improvement
is too small to prioritize over the remaining regional gaps, so the
candidate is not integrated into any certificate verifier.

## Three-Update Regional Diagnosis at 10%

`regional_frontier_probe.py` holds the saved counting witnesses fixed
and deletes selected transition contributions. The starting bounds use
64 variance parts and an output-tilt scale of 39/40. At reference
activities 0.175 and 0.20, the floating log2 bounds are +162.1393 and
+182.0288. Deleting lazy-return contributions changes them to -229.0484
and -188.2013. Deleting the contribution that moves uniform-density
mass into arbitrary mass changes them to -1660.9211 and -1707.3039.
These counterfactuals identify influential losses; they do not describe
a valid altered encoder or prove that all that slack is removable.

The actual refinement in `fixed_joint_return_probe.py` uses the exact
joint census through three active packets. It tightens only the return
entries, preserving the existing domination invariant. With the same
regional witnesses, the floating log2 bounds improve to +74.2901 and
+104.2342. Thus the gain is about 78--88 bits, but neither interval
closes. Extending the census through four packets checks all
1,820,475,000 choices at that occupancy against the independent Walsh
feedback marginal. It improves the scores to +63.8535 and +90.3329.
The fourth packet therefore adds only about 10--14 bits in this
comparison. These are selected-interval diagnostics; no complete
certificate uses the joint-return refinement.

The diagnostic normalizes each conditional count matrix after every
placement step. Tests cover tiny products and structurally zero matrices,
as well as agreement with direct products. This avoids an underflow
found in the first counterfactual run. The unchanged baseline is also
checked against its saved floating score before any comparison is used.

## Retaining Density Through Lazy Updates

The birth-class operator sends the lazy part of its uniform-density
coordinate into arbitrary mass. This loses information even when the
outgoing weighted distribution still has a useful density bound.
The diagnostic `lazy_density_probe.py` tests a different allocation
without changing the other rows or their return-to-zero bounds.

Fix an occupancy j and a nonzero outgoing state b. If the incoming
measure is bounded by U/L at each nonzero state, its lazy contribution
at b is at most alpha U g_j(b)/L, where

    g_j(b) = E[z^wt(A(b+CX)+X) 1{b+CX != 0} | J=j].

Here the j packet positions form a uniform subset, and their nonzero
values are independent and uniform. A common bound g_j on these
quantities therefore permits sending alpha g_j into U instead of
alpha Hbar_j into arbitrary mass. The existing zero-state bound is
retained separately. This is a density envelope, not a claim that
the conditioned state is uniform.

The probe reuses the density bounds from the fixed-occupancy analysis.
The initial cap averages z^max(0,d_A-wt(X)). The exact one-packet census
can sharpen it. For j>=2, conditioning on the other j-1 packets gives
the checked cap

    W/(W-j+1) * (((1+1/z)^4-1)/15)^(j-1) * G_1,

where G_1 also bounds the one-packet quantity at target zero.
The restricted location average costs W/(W-j+1); the other inputs
cost their averaged output-weight perturbation. See
[the one-packet argument](OCCUPANCY_BOUND.md#exact-one-packet-refinement).

An empty lazy step has a separate exact allocation. It leaves the
state unchanged, so each expansion-weight class E_l receives mass
alpha |E_l| z^l/L in its existing F_l coordinate. The refreshed U
entry remains unchanged. A selected occupancy cutoff determines which
nonempty lazy steps retain density. Each such choice gives a valid
local comparison, but need not improve the whole-length bound.

Small exact tests enumerate every eight-bit input and every three-bit
entering state. They check density caps and ordered transition products
for two and three updates, including degenerate feedback maps.
The first iid screen with the one-packet refinement, activity 0.20,
and output tilt 0.12 improves by about 460 bits when density is retained
through six active packets. At activity 0.175 and the same tilt, the
best tested gain is below one bit. These are floating diagnostics,
not regional or whole-code certificates.

The local bounds now live in `lazy_density.py`. The regional verifier
records an integer `regional_lazy_density_through` cutoff and recomputes
the one-packet census with outward dyadic weights. A separate
`regional_joint_return_through` flag selects the exact return census.
The latter is cross-checked against the independent integer Walsh
marginal once per fresh model; only its unweighted integer counts are
reused across output tilts. Neither flag accepts a saved numerical cap.

The regional calculation applies the return refinement first, then
the density refinement. The return formula uses the original refreshed
U entry, before it gains a lazy-density contribution. At 10% distance,
q>=97, and three updates, fresh 256-bit checks give:

| Reference activity | Existing bound | Density through six | Also exact returns through three |
|---|---:|---:|---:|
| 0.175 | +162.1393 | -158.6102437806 | -220.1614986164 |
| 0.20 | +182.0288 | -9.2073982920 | -68.7217789318 |

The entries are log2 bounds for complete mean intervals centered at
7/183 and 3/67, each of radius 1/1000000. They use 64 variance parts
and output-tilt scale 39/40. The combined bound clears 40 bits on both
intervals. Independent 384-bit evaluation agrees to the displayed digits. These two
intervals do not cover the full composition domain or the sparse prefix.
The scalar cover can now propose either refinement and keeps the old
bound whenever it is better. Previously saved witnesses retain their
original interpretation.

The fixed-occupancy sparse builder uses the same local refinements in
the same order, before averaging ordered placements. It regenerates and
cross-checks the return census once per invocation and recalculates the
density census for each output tilt. No new independence assumption is
introduced: conditional packet locations and values are exactly those
in the local bound above. The flags select a complete comparison
operator, not entrywise minima of incompatible mass allocations.
They default to disabled because a valid density allocation need not
improve every tilt or occupancy. Full assembly can select these flags
but must still regenerate and sum all sparse occupancies.

With two updates and 9.5% distance, the combined bound also verifies
three intervals centered at the tilted means corresponding to reference
activities 0.15, 0.175, and 0.20. The 256-bit log2 upper bounds are
-886.7218, -498.2154, and -381.6294. These checks use direct count masses
and otherwise the same interval radius, variance resolution, and tilt
scale as above. Independent 384-bit evaluation agrees to the displayed
digits. Complete dense and sparse coverage remain necessary.

```text
python -B research/workstreams/permutation_locality/gf16_packets/regional_count.py --distance 1/10 --updates 3 --minimum-groups 97 --activities .175 .2 --variance-bins 64 --tilted-atom --fine-tilts --tilted-variance --tilt-scales 39/40 --lazy-density-through 6 --joint-return-through 3 --precision 384 --output tmp/gf16-r3-d10-joint-lazy-selected-p384.json
```

## Density Supplied by the GF Feedback

The preceding refinement requires a density bound on the entering state.
The GF packet values can also supply density when the entering state is
fixed or otherwise arbitrary. This gives another possible replacement
for the lazy branch that currently enters the arbitrary-mass coordinate.

Fix a nonzero entering state a. Conditional on j active packets, let X
have the same uniform location and independent nonzero-label distribution
as above. For each expansion profile h, define a bound kappa_j(h) such that

    E[z^wt(Aa+X) 1{CX=c} | J=j] <= kappa_j(h)

for every a with profile h and every feedback target c, including zero.
Here h counts the four-bit packets of Aa at each weight from zero to four.
For an outgoing nonzero state b, the lazy branch has c=b+a. Therefore
an entering measure of mass M contributes at most alpha kappa_j M at b,
where kappa_j is the maximum over the allowed entering profiles.
The density convention U/L permits allocating this branch entirely to U,
with coefficient alpha L kappa_j. A birth class uses only profiles with
its specified expansion weight. The refreshed contribution and the
existing return-to-zero bound remain unchanged.

`fiber_density.py` computes two bounds and takes their minimum. The first
reuses the feedback-rank argument: a fiber on selected coordinates is
an affine binary space, regardless of its target. Dropping nonzero-packet
restrictions and selecting systematic coordinates bounds its output
enumerator by a power of 1+z, multiplied by the fixed output weight
outside those coordinates. The existing annihilator census bounds the
average over selected packet locations.

The second bound selects the lightest outputs up to an event-count
budget. Its budget uses the larger of the exact zero-feedback probability
and the bound on nonzero feedback atoms. The previous return-only budget
is insufficient here: a nonzero outgoing state can equal the entering
state, so its feedback target is zero.

The prototype applies this allocation only on a selected positive packet
count interval. It does not change the zero-input step or the U row.
Exhaustive eight-bit-input, three-bit-state tests check every feedback
target and every nonzero outgoing state for two and three updates.
The tests include identically zero and rank-deficient feedback maps.
The actual-size test at reference activity 0.30, three updates, and 9.9%
distance does not improve the saved bound. Replacing all eligible M/F
rows for packet counts 2--32, 3--32, 4--32, 6--32, or 8--32 gives
floating log2 scores +604.9566, +590.8621, +585.3699, +483.4942, and
+473.5934. The unchanged local allocation gives +290.1241 at that
witness. This prototype is not selected by the certificate verifier.

The density coefficient L kappa_j can exceed the original mass bound.
Replacing one with the other therefore need not improve later steps.
Entrywise minima between incompatible allocations would not be justified.
The negative experiment does not establish a distance limit of the code.

### Retaining Both Mass and Density

The state classes already used for births also partition every nonzero
outgoing state. Let E_l consist of states b for which wt(Ab)=l, and let
n_l be its exact size. For a fixed entering class and local occupancy j,
write H for the existing upper bound on the lazy branch's total weighted
mass, including its factor alpha. Let kappa be the corresponding bound
above, maximized over that entering class. For any entering measure of
mass M, the lazy outgoing measure nu satisfies

    nu(E_l) <= M min(H, alpha n_l kappa).

Both terms bound the same quantity: the first by the branch's total
mass, the second by summing its pointwise bound over E_l. Thus this
minimum is valid. Allocate the outgoing measure restricted to E_l to
the existing F_l coordinate with that coefficient, and remove the old
lazy M coefficient. Keep the return-to-zero and refreshed contributions
unchanged. This proves the local comparison for an arbitrary entering
measure by linearity; it does not assume independent entering-state bits.

The F_l coordinates now include later lazy transitions, not just births.
Their invariant remains unchanged: arbitrary nonnegative mass supported
on E_l. The class bounds may sum to more than H. This loses the joint
total-mass constraint, but remains an upper bound. No entrywise minimum
between alternative coordinate representations is taken.

`fiber_density.candidate(..., allocation='classes')` implements this
allocation on a selected positive occupancy interval. Exact small-state
tests enumerate every input and check each outgoing class, for two and
three updates, two injective expansion maps, and ordinary or identically
zero feedback. The regional verifier accepts it through the explicit
integer witness fields `regional_feedback_classes_from` and
`regional_feedback_classes_through`. It reconstructs the integer census
for each fresh model and the weighted fiber bounds at each output tilt
and precision; it does not import numerical caps from a saved receipt.
Old witnesses keep their original interpretation.

At reference activity 0.30, three updates, and 9.9% distance, the same
fixed witness now gives floating log2 scores -398.1013, -398.7713,
-395.9749, -371.7946, and -247.3122 for starts 2, 3, 4, 6, and 8,
respectively, all through 32. The old allocation gives +290.1241.
These are diagnostics, not a full-domain certificate. A fresh 384-bit
outward check for start 3 gives log2 upper -398.7713020749524 over the
entire tilted-mean interval 9/121 +/- 1/1000000, including all 16
variance parts and all q>=49. This is a checked interval, not coverage
of all means. Its receipt is
`tmp/gf16-r3-d099-feedback-classes-p03-p384.json`, SHA256
`a2c2dfcbee6c4a2356609701f39bc1c86e9231e32fb9d00529693667e53a4435`.
The scalar proposer can select this refinement after the earlier
regional alternatives fail; it retains the previous bound if better.
A separate fresh 256-bit calculation agrees to the displayed digits.

The same argument has a density-bounded entering version. An entering
measure satisfying mu(a)<=U/L is dominated by U times the uniform
measure on nonzero states. Replace the maximum kappa by

    kappa_bar = (1/L) sum_{a!=0} kappa_j(h(a)).

The exact profile multiplicities evaluate this sum. If the U row still
allocates its lazy branch to M, that branch can instead be divided
among the F_l classes using min(H, alpha n_l kappa_bar), where H is
the existing U-to-M coefficient. Branches already allocated to U or F
are left unchanged. This allocation is tested against exact
multistep positive terminal costs. At reference activity 0.35, the best previously tested output
tilt (scale 9/10) gives floating log2 score +517.2243. Allocating the
remaining U-row lazy mass to classes lowers this to -9.1928, without
changing that witness. This improves the bound but does not yet supply
the required 40-bit margin, even for that selected interval. At activity
0.40, the score improves from +3594.5604 to +3028.6192 and still fails.

Replacing complete U rows also allows revisiting the earlier lazy
density allocation. Rebuild the source rows before that allocation,
apply the class-mass argument, and select those complete rows for the
chosen packet counts. This avoids subtracting rounded density bounds
and leaves other rows unchanged. Both row choices bound the same
coordinate invariants, so either is valid independently at each count.

For counts 3--32, the activity-0.35 score improves to -74.4654. Fresh
384-bit outward evaluation verifies log2 upper -74.4654194057170 for
the whole tilted-mean interval 21/229 +/- 1/1000000 at 9.9% distance,
q>=49, three updates, and 16 variance parts. The output-tilt scale is
9/10. The receipt is `tmp/gf16-r3-d099-feedback-uniform-p035-p384.json`,
SHA256 `dbc5c5ffd02895c24e3bacbceef724ce0955d7fa9514cd0e26ef4a3e061883ab`.
A fresh 256-bit check agrees to the displayed digits. The optional witness fields
`regional_feedback_uniform_classes` and
`regional_feedback_uniform_replace` select these two refinements.
The scalar proposer keeps an earlier bound when it is better.

At activity 0.40, complete U-row replacement still gives floating
score +2980.9415. The next argument retains one shared mass budget
across the classes; see [MASS_DENSITY.md](MASS_DENSITY.md). Neither
selected interval nor the new diagnostic completes the whole-code proof.

# Four-bit proof attack: fresh history and stronger mixing

The target remains the four-bit independent-row construction at K=2^20,
with BCH[256,128] and IMT(128,19). The completed two-bit, two-update
certificate remains the control. This investigation changes the proof
state first; changing the number of updates is a separate ensemble.
No production encoder or published performance number is changed.

## What the new proof state remembers

The previous continuation experiment kept each current packet shape
together. Its fresh-state coordinate still forgot which packet created
that state. This can combine a cheap activation with an incompatible
worst-case continuation. It also charges the minimum expansion weight
again after each empty step.

Write E for the state-expansion map and B for the input-feedback map.
For a packet weight a in {1,2,3,4} and expansion weight v, define

    S[a,v] = { Bx : x occupies one four-bit window, |x|=a, |E(Bx)|=v }.

The current maps have nonzero, distinct feedbacks for all 480 single-packet
inputs. Let n[a,v] be the size of a nonempty S[a,v]. A new coordinate f[a,v]
means that the corresponding tilted submeasure is bounded pointwise by
f[a,v]/n[a,v] on this set. It does NOT assert that a conditioned state is
uniform. The coordinate has terminal weight one, unlike the auxiliary
mature-density and tail coordinates.

When a weight-a packet activates state zero, it creates this component
with coefficient exp(-theta*a) rho^(a=4) n[a,v]/(32 binom(4,a)). The current
packet shape must be retained until applying the continuation potential;
independently maximizing these new target entries would sum incompatible
activation choices.

For R updates, put alpha=2^(-R), beta=1-alpha, and S=2^19-1. At an empty
step the lazy component stays in S[a,v] and pays exactly
alpha exp(-theta*v). The refreshed component has pointwise coefficient
beta exp(-theta*v)/S on the nonzero state space. Thus successive empty
steps retain the expansion weight instead of repeatedly using 48.

## Nonempty steps

Fix the current packet-weight shape h. The input X is uniform over the
distinct-window and lane assignments for h. Conditional on the assignment
of packet labels to epochs, these choices are independent of the incoming
state. Let Q be uniform on S[a,v] for defining the following linear
domination calculation, and define

    A = E[exp(-theta*|E(Q)+X|)],
    B0 = E[exp(-theta*|E(Q)+X|) 1_{Q=BX}].

The expectations here average Q and X independently. They bound the actual
fresh component by the pointwise inequality defining f[a,v]. Include the
factor rho for each weight-four packet after computing A and B0.

The lazy return coefficient is alpha B0; the lazy nonzero mass coefficient
is alpha(A-B0). The refresh return coefficient is at most beta A/S.
For every expansion class of size A_v, the corresponding refreshed
uniform-envelope mass is at most beta A A_v/S. This last statement is a
pointwise bound after translating the uniform nonzero distribution by BX,
not an assertion that the translated distribution is uniform on nonzero
states.

Exact window histograms supply an outward upper for A through eight active
windows. The joint cancellation census now retains the incoming fresh
packet weight and expansion class through three active windows. It gives
B0 from exact integer histograms. When subtracting B0 from A, the code
subtracts its lower endpoint from the upper endpoint for A.

For the lazy mature-density coordinate, let W be the total current input
weight and p_max bound the largest feedback atom. Convolution gives the
pointwise upper

    alpha exp(-theta*max(0,v-W)) min(1/n[a,v], p_max).

The existing feedback census provides p_max through six windows. Above
that, the existing conditional atom bound is used. Unrefined occupancies
retain the old valid envelopes. The code can also forget a new class by
dominating it with the old fresh component, or with mature mass one,
density 1/n[a,v], and the applicable tail indicators. Those alternatives
bound the same outgoing lazy and refreshed components. Taking their
componentwise minima does not splice different mass decompositions.

`candidates/fresh_history.py` implements this extension. The optional
weight-only variant retains a but not v. The global evaluation uses the
same without-replacement placement and finite continuation witness as
[SHAPE_POTENTIAL.md](SHAPE_POTENTIAL.md). Global arithmetic is binary64;
these are witness searches, not certificates.

## Stronger-mixing control

An existing bound at lazy probability alpha0 can be converted conservatively
to any smaller alpha1. Scale pure lazy contributions by alpha1/alpha0 and
pure refresh contributions by (1-alpha1)/(1-alpha0). A zero-return entry
that combines both requires care: retain the scaled old upper and add
the extra refresh contribution using a separately valid refresh upper.
Do not subtract a loose estimate of the old lazy contribution.

The empty uniform-class self-loop has a known exact lazy term; that term
can be subtracted before rescaling the refresh contribution.
`candidates/mixing_attack.py` implements these rules. The alpha1=0 case
is an ideal uniform-refresh diagnostic, not a claimed implementation.
Finite R tests use the corresponding changed setup distribution.

The same driver can delete selected transition entries to measure their
influence. Those deletion experiments are explicitly labeled INVALID
ABLATION: their results are not probability bounds or certificates.

## Validation and reproducibility

The direct single-packet checks enumerate all incoming fresh classes and
all current packet inputs on the actual maps. They compare lazy nonzero
mass, zero returns, pointwise density, and tail quantities against exact
integer-weighted sums at exp(-theta)=7/8. Empty-step checks cover two,
four, and ideal-refresh updates. The refined census sums back to every
previous fresh-state histogram. Separate tests exercise update retargeting
and exact-endpoint memo serialization.

`candidates/attack_cache.py` binds the generator sources, fixed maps,
parameters, and precision. It stores local integer counts and exact upper
endpoints under ignored tmp/, so searches need not repeat the expensive
census. Omit its directory argument when calling build to regenerate it.
The cache is a local memo, not independent evidence of correctness.

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/candidates/attack_cache.py --tilts .056 .072 .088
python -B research/workstreams/permutation_locality/independent_rows/candidates/mixing_attack.py --ablate
python -B research/workstreams/permutation_locality/independent_rows/candidates/fresh_history.py --tilt .056 --groups 96 128
```

The selected global screens use union support 200 in every active group
and output cutoff 193986 (9.25%). Neither this support class nor the three
tilts exhausts the proof obligation. Closure still requires all support
vectors and occupancies, an outward global replay, and a summed failure
bound for the intended construction.

## Results of the expanded attack (2026-09-28)

The new outward result is a **selected-event bound**, not a whole-code
certificate. Fix three updates, q=96 active four-row groups, and union
support exactly 200 in every active group. At output cutoff 193986, the
expected number of bad messages in this class, summed over the choices
of active groups, is less than 2^-374.76. The row messages and packet
weights within those union supports are not fixed. The proof uses their
authenticated weighted outer count and the universal inner envelopes.
It concerns the independently sampled ideal setup described in this
workstream; the two-update sparse certificate cannot be added to it
without a separate argument for the changed update count.

`candidates/selected_replay.py` verifies this inequality at 192-bit global
precision, with a second global replay at 384 bits. It restores and checks
the actual Flint precision after the legacy BCH helper imports. It uses
regenerated, source-bound local bounds with six-window
feedback counts. It does not consume the historical ten-window records
discussed below. The binary64 finite-product proposal is -378.08845 bits;
the common-potential replay gives log2 upper -374.76513. The difference
is slack from the finite terminal factor, not numerical uncertainty.

For a positive column v normalized by v[Z]=1, the verifier computes the
regional action Pv using exact combinatorial placement coefficients and
outward arithmetic. It checks Pv <= lambda*v coordinatewise. If tau is
the terminal column and c=max_i tau[i]/v[i], then

    e_Z P^256 tau <= c lambda^256.

The final bound multiplies this by exp(theta*193986), binomial(2048,96),
and the 96th power of the weighted outer count divided by the positive
binomial conditioning probability. The numerical search supplies only
rational p and v; no numerical score is trusted in the replay. Exact
small-geometry tests enumerate every slot subset with noncommuting local
matrices and check the vector placement recurrence.

### What stronger mixing buys

Here are selected binary64 log2 first-moment bounds from the current
six-window local family, at theta=.072 and all-four penalty rho=.9.
All rows have common union support 200 and cutoff 193986. Positive values
mean that this upper bound is uninformative, not that bad codewords exist.

| Active groups | Two updates | Three updates | Four updates | Eight updates |
|---:|---:|---:|---:|---:|
| 96 | +1399.11 | -378.09 | -1028.16 | -1496.43 |
| 128 | +5467.50 | +3450.85 | +2788.41 | +2262.78 |

The tested homogeneous supports 160,176,192,200,208,216,224,240,256 all
have negative proposals at q=96 with three updates. Among those samples,
u=200 is worst. The same test at q=112 already has positive proposals
around u=192..224. This is not interpolation or coverage of the omitted
supports, and says nothing by itself about unequal support vectors.

Increasing the output tilt changes the comparison substantially. At
theta=.088, q=128, u=200, eight updates give +541.12, versus +412.01
in the ideal uniform-refresh limit. Ideal refresh is a diagnostic, not
an implemented encoder. These values do not prove that more updates are
useless, but they show that a fixed-tilt update-count sweep is insufficient.

Extending the tilt grid to .096, .104, and .12 did not improve the tested
q=128, u=200 bound. Even the ideal-refresh proposals rise to +1345.15,
+3188.79, and +7053.69, respectively. This rejects that simple retuning
direction for the current envelope; it does not bound the true distance.

The earlier measured four-bit implementation costs remain 6.442 ms for
two updates and 6.727 ms for three updates at K=2^20. No four- or
eight-update performance was measured in this investigation.

### Which proof refinements helped

* Remembering fresh packet and expansion classes did not yield a robust
  global improvement. At theta=.056 the initial split improved the tested
  bound by only about 10 bits. At theta=.072, even after exact one-packet
  density and tail refinements, its q=128 proposal was +5842.34, worse
  than the unsplit +5467.50. This branch remains experimental.
* The adjoint diagnostic puts substantial sensitivity at seven occupied
  windows, immediately above the exact-density cutoff of six. High
  occupancies above sixteen have negligible sensitivity at this witness.
  Thus extending an already small local census is more promising than
  unbounded enumeration of high occupancies.
* Historical exact-ten/joint-four/density-ten local records improve
  q=96 even with two updates, but do not close q=128 at theta=.072.
  One such record plus fixed coupled replacements gives about +2731.38
  for three updates at q=128. These records pass their stored checksums,
  but their source hashes no longer match today's generators. Their
  reuse is explicitly **diagnostic only**, not a new authenticated replay.
  `candidates/strong_mixing.py --rebuild` regenerates the base refinement
  independently when a promising witness warrants the cost.
* Averaging packet weights under the existing row-conditioned Bernoulli
  reference can help. With three updates and q=128, the class with all
  four rows of every group at weight 80 has a -659.42 proposal. The
  analogous weight-96 class remains at +4132.57 in this coarse screen.
  Two-active-row classes are much easier in the tested range. These are
  different classes from fixed union supports and must not be added to
  those covers as though they formed a partition.

The row screen removes the proof-only rho penalty before averaging. Up to
eight occupied windows this is shape-specific; above eight it pays the
conservative factor rho^-j. It does not yet reuse the stronger averaged
high-occupancy row model, so a poor row-screen result is not an obstruction
to that proof method.

### Post-feedback mixing: a tested alternative, not a replacement

We also tested a different inner recurrence:

    y = x + E(q),       q' = R_post(R_pre(q) + B(x)).

The current inner omits R_post. An independent post-feedback update fixes
zero and acts on each nonzero state as alpha times identity plus
1-alpha times uniform refresh. On the existing envelope, every nonzero
coordinate retains its alpha-scaled lazy contribution. Only the genuine
mass coordinates F,M,U contribute the new uniform component; C,L48,L56
are auxiliary constraints, not additional mass. This yields a positive
eleven-coordinate matrix H, and the epoch envelope becomes T_pre H.
`candidates/post_mix.py` implements this composition. Small exact-kernel
tests check zero preservation, the uniform classes, and mass accounting.

At the tested three- and four-update budgets, this ordering change did not
beat allocating all updates before feedback. For example, at theta=.072
and q=128,u=200:

| Pre / post updates | Selected log2 proposal |
|---|---:|
| 3 / 0, current ordering | +3450.85 |
| 2 / 1 | +3528.54 |
| 4 / 0, current ordering | +2788.41 |
| 3 / 1 | +2828.53 |
| 2 / 2 | +2916.99 |

No implementation or performance claim follows from this experiment.
It does not justify replacing the current inner. Its finite-update and
ideal-refresh screens are reproducible with:

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/post_mix.py
```

### Next proof work

Retain the two-bit full certificate and implementation as the control.
For four bits, prioritize a complete three-update support cover near
q=96, and the actual row-distribution/affine-mixture bound for the denser
regime. Reuse the ten-window density refinement where it matters. Avoid
spending the next iteration solely enlarging the fresh-state split or
increasing the update count at a fixed proof tilt. A whole-code result
still needs sparse coverage, dense coverage, and a verified sum for the
same construction and output cutoff.

Additional reproduction commands (run from the repository root):

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/selected_replay.py
python -B research/workstreams/permutation_locality/independent_rows/candidates/selected_replay.py --precision 384
python -B research/workstreams/permutation_locality/independent_rows/candidates/mixing_attack.py --tilts .056 .072 .088 --rounds 3 4 --groups 96 112 128 --supports 160 176 192 200 208 216 224 240 256
python -B research/workstreams/permutation_locality/independent_rows/candidates/row_mixing_screen.py
python -B research/workstreams/permutation_locality/independent_rows/candidates/strong_mixing.py --historical-directory tmp/independent-row-local-families
```

The last command requires the old local memos and is not an independently
reproducible certificate. The 384-bit command increases global replay
precision while retaining the same 192-bit local upper endpoints. Raw
screens and memos stay under ignored tmp/. No production encoder,
parameter default, published performance claim, or paper text is changed.

The candidate test suite passes 85 tests after these changes, including
exact fresh-state enumeration, update-count conversion, packet-law
averaging, the post-refresh kernel, and outward vector-placement checks.

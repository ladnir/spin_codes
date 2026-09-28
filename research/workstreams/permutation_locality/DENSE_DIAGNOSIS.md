# Diagnose the two-update intermediate-occupancy gap

The fourteen-group certificate does not establish a route to all occupancies.
The next diagnostic fixes 64 active groups and equal group support sizes.
Every score below includes group locations and the outer count cap. These
are binary64 proposals, not outward certificates or bad-word witnesses.

## Larger tilts alone do not close the tested points

With tilts .032, .064, and .1 and all-one penalties .75 and 1, the best
log2 contribution uppers at supports 80, 128, 160, 224, and 256 are
+2569.59, +3809.06, +5177.84, +4630.17, and +2079.69. Earlier smaller
tilts did better at support 80, but none of these points has closed.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --joint-witness --penalties .75 1 --tilts .032 .064 .1 --probe-supports 80 128 160 224 256

## Sensitivity is distributed across the state envelope

At support 128, tilt .032, and penalty .75, the baseline score is +3809.06.
Artificially setting one family of active-input transition coefficients
to zero gives the following scores. None of the altered operators is a
valid bound; they identify possible targets only.

| Removed contribution | Diagnostic log2 score |
|---|---:|
| Direct zero-to-zero | +3789.67 |
| Fresh-state cancellation | +3704.44 |
| Mature-density cancellation | +3396.10 |
| Uniform-class cancellation | +3675.83 |
| Mature refresh cancellation | +3794.19 |
| All active-input returns to zero | -5453.82 |

Removing any single listed cancellation mechanism does not close this
point. Even halving every return-to-zero coefficient only gives +3429.26.
The results do not justify attributing the gap to one loose coefficient.

occupancy_sensitivity.py uses a normalized hypergeometric recurrence for
binary64 placement diagnostics. It agrees with exact rational products
of noncommuting matrices on small examples. Its baseline reproduces the
outward-operator screen above within 1e-7 bits. It is separate from the
certificate's outward recurrence. Optional snapshot files contain only
diagnostic matrices and are never certificate inputs.

    python -B research/workstreams/permutation_locality/occupancy_sensitivity.py

## Keep total outer weight coupled to support

Fix a nonzero ordered tuple of four BCH words. Let u be its union support
size and w the sum of its four word weights. The permutation preserves w.
The existing envelope maximizes over the nonzero weights of arriving
columns. Counting many tuples while allowing sparse column shapes may
lose useful information; this is a proposed target, not a diagnosed cause.

For a>=1, multiply every input-weight-w_epoch local transition by
a^w_epoch before taking its shape maximum. The complete tilted moment
then carries a^w. Its probability upper is multiplied by a^-w. The outer
count must therefore use the same weighted count, not its ordinary count.

Let U(u) bound the number of nonzero four-tuples with union support at most
u. Let A(z) have authenticated upper coefficients for the BCH enumerator,
including the zero word. The coefficients of A(z)^4-1 bound tuple counts
by total weight. Let S(w) be their cumulative sum. For fixed u,

    C_u(w) = min(U(u), S(w))

bounds the number of tuples with support at most u and total weight at
most w. Such tuples have total weight at most 4u. Summation by parts,
using the decreasing weight a^-w, bounds their weighted count by

    sum_{w=0}^{4u} (C_u(w)-C_u(w-1)) a^-w,    C_u(-1)=0.

These differences are extremal coefficients of a cumulative upper bound,
not true weight-shell counts. total_weight_caps.py rounds the result upward
to an integer and checks monotonicity in u. Exhaustive four-tuple enumeration
on three small binary codes checks the bound for a in {1,1.02,1.1,2}.

The cover's optional --weight-tilt flag propagates the same multiplier
through fresh, multi-window, cancellation, uniform, and mature-tail terms.
The direct component regression includes it. It currently requires all-one
penalty 1; no unsupported joint count in (w,J) is assumed. Default a=1
retains the existing certificate path. The first 64-group screen with a=1.02,
all-one penalty 1, and tilts .016, .032, .05 gives log2 uppers +2340.24,
+5013.04, +6181.23, +5247.93, and +2922.30 at supports 80, 128, 160, 224,
and 256. It does not improve the best previous bounds. This single trial
does not rule out other tilts, but supplies no reason to prioritize a large
parameter sweep. The default fourteen-group certificate was replayed at
192 bits after these changes and reproduces its previous upper exactly.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --weight-tilt 1.02 --fresh --window-average --multi-average --fresh-collision --zero-moment --mature-tail 64 --pair-tail --joint-witness --penalties 1 --tilts .016 .032 .05 --probe-supports 80 128 160 224 256

## Exact five- and six-window zero-feedback counts

feedback_character_census.py computes an exact alternative to brute-force
input enumeration. For a dual character c in F_2^19, let r_i(c) be the
number of negative signs it induces on feedback columns in window i.
For a four-bit input weight b, the signed sum over its masks is

    K_b(r) = sum_v (-1)^v binomial(r,v) binomial(4-r,b-v).

The coefficient of x_1^n1 ... x_4^n4 in

    product_{i=1}^{32} (1 + sum_{b=1}^4 K_b(r_i(c)) x_b)

counts distinct-window assignments with those weight multiplicities,
signed by their total feedback. Averaging over all 2^19 characters removes
exactly the nonzero-feedback assignments. The coefficient denominator is
binomial(32,j) j!/product_b n_b! times product_b binomial(4,b)^n_b.

There are 9218 distinct histograms of r_i among the 524288 characters.
Grouping identical histograms makes the coefficient calculation small.
For j<=6, checked combinatorial bounds keep every coefficient, partial
sum, and weighted character sum within signed 64-bit range. The result is
exactly divisible by 2^19 and lies between zero and its denominator.

The worst zero-feedback probabilities over weight shapes are 1/201376
for five windows, and 17257/3711762432 for six. The three-/four-window
results agree with the independent direct-rank census on all 320 ordered
three-/four-window weight shapes, including zero-count shapes. These new counts
are not yet used in a complete occupancy certificate. The sensitivity
result above suggests that this direct-zero refinement alone is insufficient.

    python -B research/workstreams/permutation_locality/feedback_character_census.py

## Next structural question

Four binary BCH words can be identified with one word in the scalar
extension of the BCH code to F_16. Its symbol support is their union support.
Indeed, choose an F_2 basis of F_16 and expand each symbol in that basis.
The binary generator matrix acts independently on the four components.
Thus the grouped-support enumerator is an ordinary Hamming enumerator of
this extension code. This is an analysis viewpoint, not an encoder change.

The current rank/basis bounds do not use the full F_16 enumerator identities.
Investigate whether those identities or the binary dual structure improve
the intermediate-support counts. Separately, the state envelope maximizes
each matrix entry over input shapes; determine whether those maxima combine
incompatible shapes and how much a joint bound would save. Neither avenue
has been shown to close the gap. They are more informative next tests than
another performance change or an unbounded witness-grid search.

## The elementary F_16 bound is already too weak

The existing joint_oa_caps.py applies the q-ary Christoffel bound with
q=16 and orthogonal-array strength 29. It therefore already uses one
consequence of the scalar-extension interpretation. Conditional on that
strength, its shell log2 caps at supports 80, 96, 128, 160, 192, 224, and
256 are 399.13, 403.59, 414.37, 429.22, 453.75, 500.59, and 489.63.
At support 128 this is weaker than the existing rank/basis count. Merely
renaming the four-word distribution as an F_16 code does not improve it.
Stronger constraints from that interpretation remain untested.

Extending shortened_bound.dimension_caps from last_lp=104 to 192 accepted
74 exact rational witnesses and rejected 81 numerical proposals. The
dimension cap at length 128 improves from 50 to 49; at length 160 it stays
81. These results do not replace the current certificate inputs. Rejected
proposals provide no bound.

## A weighted fresh-state envelope saves little at the tested point

Let p_a be the feedback distribution of one uniformly placed four-bit
window of weight a, for a in {1,2,3,4}. The old fresh coordinate bounds a
mixture of these distributions. It forgets which a caused activation.
fresh_gauge.py tests generators p_a/kappa_a, where every kappa_a >= 1.
The fresh coordinate now bounds the sum of the coefficients in this
weighted representation; it is not necessarily probability mass.

Activation by weight b costs kappa_b. A bound on the contribution from
generator p_a is divided by kappa_a. Empty lazy steps preserve the fresh
type and retain their old pointwise multiplier. The old outgoing bounds
remain valid because p_a/kappa_a <= p_a. The old terminal coefficient 1
also remains valid. This changes only the proof representation, not the
encoder or its setup distribution.

The experiment replaces selected outgoing coefficients for at most four
active windows. It retains the old envelope elsewhere. For one window,
the feedback supports for different a are disjoint, so lazy cancellation
requires the new input weight to equal a. The calculation keeps this
restriction, the exact cancellation-output histogram, and the averaged
refresh contribution before maximizing over a and b. Independent direct
summation checks the single-input cancellation and mass bounds for all
16 pairs (a,b), under each of two gauges: 32 cases in total.

At the same 64-group, support-128, tilt-.032, penalty-.75 point:

| Gauge (kappa_1,...,kappa_4) | Diagnostic log2 upper |
|---|---:|
| Original envelope | +3809.0603 |
| (1,1,1,1) | +3809.0452 |
| (1,1,1,1.25) | +3794.6919 |
| (1,1,1,1.5) | +3790.0837 |
| (1,1,1,2) | +3849.4552 |
| (1,1.2,1.2,1.5) | +3821.4554 |
| (1,1.5,1.5,2) | +3868.1227 |

These are binary64 diagnostics, not outward certificates. The best tested
gauge saves about 19 bits, not enough to close the point. A genuinely
shape-resolved state analysis could do better; this restricted experiment
does not evaluate it. It supplies no reason for a large gauge sweep.

Reproduce the baseline cache with occupancy_sensitivity.py --snapshot
tmp/r2-q64-sensitivity.npz, then run:

    python -B research/workstreams/permutation_locality/fresh_gauge.py

At this point the all-one-weighted outer count cap has log2 value
352.5504 per group. Holding the inner moment fixed, a count improvement
of (3809.0603+40)/64 = 60.1416 bits per group would give a 40-bit upper
for this selected support vector. This is a target for count refinements,
not a claim that such a bound exists or that the full range would close.
In particular, the current cheap-basis argument does not enforce the
joint weight constraints of all fifteen nonzero combinations of a
rank-four tuple. Those constraints are a concrete next counting target.

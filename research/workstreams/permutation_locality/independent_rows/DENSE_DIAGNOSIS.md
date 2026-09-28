# Where the middle-occupancy bound loses strength

The complete sparse cover and selected row classes are different results.
This diagnostic concerns only q = 80 active groups, each containing four
nonzero rows of weights in [64,192]. It does not cover occupancy 80 in
full. The construction, threshold, and counting measure are specified in
[ROW_CONDITIONING.md](ROW_CONDITIONING.md).

Before the density refinements below, the screened log2 upper at
Bernoulli witness p = 1/2 was
+1351.499 at output tilt .052, and +1154.443 at tilt .06. Neither proves
a useful probability bound. The calculation includes class-specific tail
arrivals, high-occupancy output averaging, and region-boundary clamping.
No codeword search or performance experiment is involved.

## Transition ablations

`row_sensitivity.py` sets selected transition coefficients to zero and
recomputes the binary64 score. These modified matrices are **not valid
upper bounds**. They only measure how strongly the calculation depends on
particular estimates. Here C is a pointwise mature-state density bound;
it is not probability mass. Z denotes the zero state.

| Deleted contribution in nonempty epochs | Score at tilt .052 | Score at tilt .06 |
|---|---:|---:|
| None | +1351.499 | +1154.443 |
| C to Z | +1125.018 | +171.177 |
| C to C | +1196.611 | +285.880 |
| Uniform-class mass to Z | +1235.043 | +1134.349 |
| All returns to Z with at least four nonzero packets | +1122.965 | +527.737 |
| All returns to Z with at least nine nonzero packets | +1338.104 | +1091.322 |

Thus deleting the cancellation contributions beyond the exact three-window
census is insufficient by itself at either tested tilt. Extending that
census alone is not the next priority. The calculation is more sensitive
to density persistence at the larger tilt. These observations concern the
bound, not the frequency of any failure mechanism in the actual code.

```sh
python -B research/workstreams/permutation_locality/independent_rows/row_sensitivity.py --occupancies 64 80 128 --tilt .052
python -B research/workstreams/permutation_locality/independent_rows/row_sensitivity.py --occupancies 64 80 128 --tilt .06
```

## Retaining total input weight

The reference experiment includes row weights outside [64,192]. We tested
whether excluding those weights through their aggregate sum removes the
gap. Let n = 4q be the number of nonzero rows, and let W be their total
weight. The original class requires 64n <= W <= 192n.

For one shuffled row, define d_w = Abar_w/binomial(256,w). A particular
array of n rows has capped mass equal to the product of its row densities.
For any auxiliary Bernoulli witness s in (0,1), that product is at most

    gamma_I(s)^n s^W (1-s)^(256n-W).

On a band L <= W <= U, divide this expression by the Bernoulli(p)
reference density. The maximum ratio occurs at U when s >= p, and at L
otherwise. This gives a valid scalar factor for that band. Different
disjoint bands may use different reference probabilities, auxiliary
witnesses, and output tilts. Their bounds are then summed. The reference
moment for each band can include inputs outside that band.

`aggregate_weight.py` supplies exact small-instance checks and outward
scalar factors for this reduction. `weight_band_screen.py` is only a
binary64 search over those factors and reference moments.

At q = 80, we partitioned [20480,61440] into 32 consecutive bands. The
probability grid was .36, .40, .44, .48, .50, .52, .56, .60, .64, .70,
and .80, with output tilts .052 and .06. The resulting log2 sum remained
positive, about +1156.028. The largest terms include W = 33290 through
37132, well inside the permitted interval. This finite grid does not rule
out better witnesses, but it gives no evidence that the aggregate weight
constraint resolves the present gap.

```sh
python -B research/workstreams/permutation_locality/independent_rows/weight_band_screen.py --groups 80 --tilts .052 .06
```

## Next calculation

The immediate target is the lazy density column. Its exact expression
contains both the feedback Bx and the output E(t+Bx)+x. Bounding their
effects separately can overestimate how much tilted mass remains
concentrated at a target state t. A conditioned one-window bound is now
available as `--column-density`; it has exact toy checks and independent
review. The stronger bound below is enabled by `--feedback-density 6`.
Neither refinement changes the encoder or supplies the missing full-code
claim.

### Feedback distribution and expansion classes

Let E map a 19-bit state to its 128-bit expansion, and let B map an input
epoch to its 19-bit feedback. Fix a shape of nonzero four-bit packets,
with total input weight W. Let n_r count inputs of this shape with feedback
r, and set D = sum_r n_r. These counts include the uniform choice of
distinct windows and the uniform masks within those windows.

On the lazy branch, a nonzero target t comes from source s = t + Bx.
For each nonzero expansion weight v, define

    G_v(t) = sum_{s: weight(E(s))=v} n_(t+s).

The triangle inequality gives weight(E(s)+x) >= |v-W|. Consequently, an
incoming mature measure with pointwise density at most C contributes
at most

    C * alpha/D * max_(t != 0) sum_v G_v(t) exp(-lambda |v-W|)

to the outgoing mature density, where alpha = 1/4 for two transvections.
The source-zero term is absent because it belongs to a separate component.
For incoming uniform-class mass U_v, use the single-class expression and
divide by the class size A_v. These are lazy-branch bounds only; the
refresh contributions remain unchanged.

`feedback_density.py` computes each G_v by exact integer XOR convolution.
It reuses the feedback character polynomials and computes five transforms
of expansion-class indicators. Before each signed int64 transform, it
checks (2^19)D < 2^63. A partial Walsh sum has magnitude at most (2^19)D:
character orthogonality restricts the second primal index to a coset,
whose size cancels the number of summed characters. This bound applies
to every intermediate butterfly, not merely the final counts.

The implementation checks that the resulting counts are integral and
nonnegative, that sum_v G_v(t) = D-n_t, and that selected translated-class
sums agree with direct enumeration. It rounds each exponential upward
to a dyadic rational. The dyadic precision is chosen so that every
weighted integer sum also fits int64. Thus it can take the maximum over
all nonzero target states without floating-point maximizer selection.
Only the final rational coefficient is converted outward to Arb.

Independent review found no blocker in the derivation, integer-width
argument, or implementation. The driver applies the new scalar bounds
only to C-to-C and uniform-class-to-C entries of the unsplit density
column, after the two-update conversion. The other transitions are
unchanged. Existing bounds remain available and are intersected with
these bounds entry by entry where they bound the same contribution.

### First production replay

The six-window census checked all 209 nonempty packet shapes. The largest
integer-width bound was 394071581644554240, below 2^63. At q = 64,
tilt .052, and p = 1/2, the complete replay for the stated four-row
middle-weight class gives

    U < 2.485528803403e-122,
    -log2(U) > 403.961674 bits.

The replay used 192-bit outward arithmetic, all three earlier refinements,
and both new density refinements. Every class member and group location
is included, but other row patterns at occupancy 64 remain outside this
result. At tilt .06 the same class has a weaker verified margin exceeding
380.232978 bits. These overlapping witnesses are alternatives, not terms
to add together.

The q = 80 screen improves from +1154.443 to +309.284 at tilt .06.
At tilt .064 it reaches +201.724, still insufficient. The following table
reports the best **binary64** score among tilts .052, .06, .064, .072,
and .08 for each larger selected class:

| Active groups | Best tested tilt | Log2 score |
|---|---:|---:|
| 80 | .064 | +201.724 |
| 96 | .064 | +1718.221 |
| 128 | .072 | +6310.930 |
| 256 | .08 | +48380.723 |

These failures do not imply that the code has low distance. They show
that tightening density persistence alone has not completed this method.
The next diagnostic should examine the remaining joint output/return
bounds and the coarse local-occupancy fallback, using the refined operators.
It should not infer their importance from the older ablations alone.

```sh
python -B research/workstreams/permutation_locality/independent_rows/row_verify.py --occupancies 64 80 96 128 256 --active-rows 4 --intervals 64:192 --tilts .052 .06 .064 .072 .08 --p 1/2 --cut 8 --full-feedback 6 --window-histogram 8 --joint-cancellation --class-tail --averaged-high --clamp --column-density --feedback-density 6 --screen-only --replay-below -52
```

The final local test suite passes 86 unit tests. Each production shape
build additionally checks 4356 atomic mixture identities and 726 direct
state inequalities. All five tested tilts passed those checks. The
earlier shared-shuffle proof files are unchanged.

## Unrestricted integration

The new density bounds need not remain restricted to Bernoulli row
classes. For the universal verifier, multiply each shape's coefficient
by its existing weight factor rho^(number of weight-four packets) a^W.
Then maximize over every shape of that local occupancy before tightening
the same unsplit density-column entries. This preserves the verifier's
full support coverage. Completeness of the shape census must be checked;
maximizing over only available shapes would not suffice. The feedback
coefficients already include alpha = 1/4, while the raw conditioned
column bound still needs that factor.

The adapter is now implemented in `universal_density.py` and enabled in
`verify.py` by `--column-density --feedback-density 6`. For each occupancy,
it checks the entire feedback shape-key set, intersects the applicable
per-shape bounds, applies the weight factors, and takes each coefficient's
maximum over all shapes. Above six packets, the conditioned-column bound
still applies through all 32 windows. Occupancy zero is unchanged, as are
all entries outside the six specified incoming density contributions.
The adapter rejects a split mature-density column.

Cache schema 2 includes both new options. Its source digest now includes
every Python module in this directory as well as the existing parent
helpers and authenticated maps. Thus previous memo files are misses,
not inputs to the revised run. The exact feedback census is shared across
the requested tilts and penalties. No floating-point snapshot is used as
certificate input.

The combined local suite passes 94 tests. New tests compare the adapter
with an independent higher-precision per-shape calculation, check both
weight factors and the single application of the lazy probability, verify
that intersections precede shape maxima, reject incomplete censuses, and
check unchanged envelope entries. Cache tests check source and option
invalidation. These implementation checks alone do not add any certified
occupancies; those require the separate outward cover replay.

### Universal support-class replay

At q = 64, the revised bounds were replayed with 192-bit outward
arithmetic for five homogeneous union-support vectors. Unlike the
four-row middle-weight class above, each vector here includes every
compatible row-weight pattern and every active-group location. It still
covers only the stated support vector, not all of occupancy 64.

| Union support in each active group | Outward log2 upper, rounded upward |
|---|---:|
| 176 | -532.069950 |
| 192 | -45.450765 |
| 208 | +72.092006 |
| 224 | -308.911856 |
| 256 | -4655.634041 |

All five selected witnesses use tilt .048 and rho = .9. The support-192
class now has more than 45.45 bits of margin, but misses the proposed
52-bit per-occupancy allocation even before including the other supports.
The positive support-208 result is vacuous as a probability bound. It
does not exhibit a low-weight codeword.

The CDF-only support-192 screen improved from the earlier +277.892
to -25.47655. Using the independently valid shell cap improves that
class further to the outward result above. Neither comparison transfers
the much larger margin of the differently restricted row-conditioned
class reported earlier.

Increasing the output-tilt grid to .052, .056, and .06 at rho = .9
subsequently gives the following 192-bit outward results. The table
rounds log2 uppers upward. All rows concern the same support classes as
above, without restrictions on the compatible row weights.

| Union support per group | Selected tilt | Outward log2 upper |
|---|---:|---:|
| 176 | .052 | -524.436189 |
| 192 | .052 | -162.333102 |
| 200 | .056 | -124.705334 |
| 208 | .056 | -238.447785 |
| 216 | .056 | -473.771563 |
| 224 | .056 | -856.033024 |
| 256 | .06 | -5880.767783 |

The earlier .048 witness remains better for support 176. These are
alternative bounds, so retain the smaller upper at each support rather
than replacing or adding overlapping results. In particular, support 208
is no longer an unresolved sampled point. Its improvement comes from a
proof witness, not an encoder change. The full q = 64 support cover is
still needed; the sampled points do not include heterogeneous vectors or
sum over neighboring support values.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .052 .056 .06 --penalties .9 --precision 192 --target-bits 52 --full-feedback 6 --window-histogram 8 --joint-cancellation --column-density --feedback-density 6 --probe-supports 176 192 200 208 216 224 256 --shell-points
```

A binary64 diagnostic also propagated the same region operators with
`cone_moment.log_moment`, enforcing the existing mature-subset inequalities
after each region. At tilt .048 and rho = .9, the support-192 and
support-208 scores changed by less than 3e-10 bits after reoptimizing the
Bernoulli witness. Thus region-boundary clamping is not a useful next
refinement at these points. This diagnostic was not installed in the
universal verifier and is not an additional certificate.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .032 .04 .048 --penalties .75 .9 1 --precision 192 --target-bits 52 --full-feedback 6 --window-histogram 8 --joint-cancellation --column-density --feedback-density 6 --probe-supports 176 192 208 224 256 --shell-points
```

A preliminary q = 35 full-cover attempt with only tilts .032 and .04
failed: its largest remaining boxes have union support 38–41. This
grid omits the smaller tilts used by the earlier sparse certificate;
its failure is not a regression of that certificate. Combining the
smaller sparse-side tilts with the revised dense-side operators now
passes complete outward covers at q = 35 and q = 40, with margins above
431.906849 and 326.709317 bits, respectively. The complete q = 41,...,48
batch also passes; its outward sum is below 2.722085e-24.
The larger tilt grid also passes a complete q = 49 cover with margin
above 142.369776 bits. The subsequent complete q = 50, 51, and 52
covers have margins above 100.745432, 70.748291, and 142.976724 bits.
Further q = 53, 54, 55, and 56 covers pass with margins above 119.291446,
98.749570, 77.846455, and 137.798272 bits. The complete q = 57 and 58
covers also pass, with margins above 106.081729 and 64.956003 bits.
Coverage therefore extends through occupancy 58,
with combined upper below 6.786362e-14.
[LOW_OCCUPANCIES.md](LOW_OCCUPANCIES.md) gives the
reproduction driver and aggregate calculation. The q = 59 cover stopped
at 175 splits with binary64 log2 upper -48.463562. It missed the -54
proposal gate and did not receive outward replay. The q = 60 cover
also stopped without replay, at 175 splits and 6677 leaves, with
binary64 log2 upper -15.829523. Neither result adds certified coverage.

### Full occupancy-64 cover: remaining gap

The sixteen-witness shell/CDF cover with prefix-rank folding reached its
200-split limit with 5968 retained leaves. Its binary64 log2 upper was
+159.302959, so it did not proceed to outward replay. The largest leaves
primarily contain supports 176–202, mixed with a few smaller and larger
supports. This is not a contradiction of the successful homogeneous
points: the complete calculation also sums neighboring support values
and includes heterogeneous vectors.

The current complete coverage therefore remains q = 1,...,58. A
selected-point run tested input-weight tilt a = 1.02 at q = 64 and 80,
with output tilts .052, .056, .064, and .072 and rho = .9. Both the inner
factor a^W and the reciprocal outer factor are retained. This changes
the proof witness, not the construction.

That direction worsens the bound. At q = 64, the outward log2 uppers
for supports 192 and 200 are +55.919898 and +96.814692, respectively.
At q = 80 they are +1046.678350 and +1210.148053. Positive values are
vacuous as probability bounds. The earlier a = 1 witnesses remain valid
and stronger; the new run does not replace them. The opposite direction,
a = .98, also loses. Its outward log2 uppers at supports 192 and 200
are +58.659876 and +62.963443 for q = 64, and +902.581774 and
+985.214162 for q = 80. This run uses the matching reciprocal outer
counts and the below-one safeguards in [WEIGHT_TILT.md](WEIGHT_TILT.md).
Thus neither tested input-weight adjustment improves the difficult points.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --occupancies 64 80 --tilts .052 .056 .064 .072 --penalties .9 --weight-tilt 1.02 --precision 192 --target-bits 52 --full-feedback 6 --window-histogram 8 --joint-cancellation --column-density --feedback-density 6 --probe-supports 176 184 192 200 208 216 224 --shell-points
```

Repeat with `--weight-tilt .98` for the second direction. The experimental
`candidates/weight_cover.py` controller can combine several weight witnesses
in one complete cover. It labels every operator with both weight factors
and uses only the matching outer counts. Four tests check this coupling,
duplicate rejection, and a small outward cover with opaque measure labels.
No complete real-parameter cover has used this controller yet.

[MASS_DENSITY.md](MASS_DENSITY.md) records a separate experimental
concentration bound. Its first screen gives a modest selective improvement,
but a blanket replacement loses substantially. It is not in the verifier
or in the aggregate certificate.

### Finer output-tilt and penalty grid

A 192-bit selected-point replay tested output tilts .052, .054, .056,
and .058 against penalties .85, .875, .9, .925, and .95. It tested
supports 176, 184, 192, 196, 200, 204, 208, 216, and 224 at q = 64.
Every selected winner still used rho = .9. The useful new witnesses are:

| Union support per group | Output tilt | Outward log2 upper, rounded upward |
|---|---:|---:|
| 196 | .054 | -138.480742 |
| 200 | .054 | -140.541282 |
| 204 | .054 | -168.541142 |
| 224 | .058 | -878.301833 |

At support 200 this improves the prior .056 result by about 15.84 bits.
It does not certify the full occupancy. The experiment gives no reason
to expand this nearby penalty grid further before strengthening the
local bounds or completing the support cover.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .052 .054 .056 .058 --penalties .85 .875 .9 .925 .95 --precision 192 --target-bits 52 --full-feedback 6 --window-histogram 8 --joint-cancellation --column-density --feedback-density 6 --probe-supports 176 184 192 196 200 204 208 216 224 --shell-points
```

The combined mass-based refinement now passes selected outward replays
at supports 192, 200, and 208. Its gain at support 200 is about 52.30
bits relative to the best finer-grid baseline above. A complete q = 64
cover is running with eighteen baseline witnesses and this alternative;
no successful occupancy-64 result is claimed yet. See
[MASS_DENSITY.md](MASS_DENSITY.md) for the operator argument and commands.

### Next local refinement: beyond the six-window cutoff

At q = 80, the best tested mass-based selected-point uppers remain
positive. Computing exact zero and expansion-class counts through ten
windows substantially improves them, even before obtaining exact peaks.
For supports 192, 200, and 208, the 192-bit outward log2 uppers are
+223.496836, +256.803777, and +158.901375, rounded upward. These are
still vacuous probability bounds. Fixed convex mixture optimization
improves the support-200 result by only another 0.61 bits.

This points to the local feedback census as the next useful refinement,
rather than more mixture tuning. An exact full-distribution extension
uses checked integer limbs to avoid overflow in Walsh inversion. Its
full-size replay passes all local checks and improves the selected
log2 uppers to +181.416331, +210.043130, and +60.370886, respectively.
[SPECTRAL_FEEDBACK.md](SPECTRAL_FEEDBACK.md)
gives the derivation, completed results, and commands. These experiments
do not change the complete occupancy coverage or its aggregate budget.

The original mass-based q = 64 cover ended at 400 splits and 12137
leaves, with binary64 log2 upper +97.719149. It did not receive outward
replay. A second complete cover uses exact feedback through ten windows
and joint cancellation/output counts at four windows, with four refined
tilts in addition to the eighteen baseline witnesses. Neither currently
has an outward occupancy certificate.

The second run subsequently ended at 300 splits and 8784 leaves, with
binary64 log2 upper -43.1317304404339. It missed the requested -54
replay gate and did not receive outward replay.

The joint cancellation/output enumeration now passes through four
windows. Together with exact ten-window feedback, it gives outward
selected q = 80 log2 uppers +141.953270, +164.918511, and +15.492140
at supports 192, 200, and 208. These improve the preceding bounds by
about 40--45 bits but remain vacuous. See [JOINT_FOUR.md](JOINT_FOUR.md).
The stronger complete-cover batch for q = 59,...,63 is now running
with these refinements. The separate density-convolution extension
targets the remaining persistence loss without changing the encoder.
Its completed 791-shape census tightens 24 scalar entries per tilt, with
mass inflation below 1.00000024. The selected q = 80 outward log2 uppers
improve to +136.909130, +158.598967, and +9.147664, respectively, rounded
upward. These are still vacuous probability bounds. Mixture optimization
improves the support-200 upper only to +158.187720;
[DENSITY_EXTENSION.md](DENSITY_EXTENSION.md)
explains which gains the current fixed replacement discards.

[OVERLAP_MEAN.md](OVERLAP_MEAN.md) gives a further local refinement.
An exact character calculation constrains the mean overlap between the
input and the expanded feedback. A convex bound then improves the
shape-maximized triangle estimates for zero returns by about 33--37%
at the tested tilts and local occupancies five through ten. Its complete
local census and exact small-model checks pass. Selected-point outward
replay at q = 80, support 200 gives log2 upper +133.142758 with the
zero-return refinement alone and +121.548507 with its translated-density
extension. Adding the quadratic zero-return bound from
[OVERLAP_SECOND.md](OVERLAP_SECOND.md) improves the latter to +85.881088.
All these uppers remain vacuous; no additional occupancy certificate
currently uses the overlap refinements.

The translated quadratic-density replay has also completed. Together
with quadratic zero returns it lowers the same selected-event upper to
+63.907200, rounded upward. It still does not close that event.

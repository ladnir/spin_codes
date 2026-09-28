# Joint cancellation and emitted weight

The existing bound sometimes combines a cancellation probability with a
triangle-inequality bound on the output weight. `cancellation_joint.py`
instead enumerates the emitted weight on the cancellation event, through
three active four-bit windows. The construction and setup distribution are
unchanged. This is an optional refinement of the existing proof envelope.

## What is counted

Fix the nonzero input-weight multiset b=(b_1,...,b_j), with j<=3.
Choose j distinct windows uniformly among 32, assign the weights uniformly,
and choose uniform lane masks of those weights. Write X for the resulting
128-bit input, S=BX for its feedback, and V=wt(ES).

The two-update state transition is lazy with probability alpha=1/4; otherwise
it refreshes uniformly over the nonzero states before adding S. Output is
emitted first. On a lazy return to zero, the entering state must equal S,
and the output is ES+X. Define W=wt(ES+X).

For this input multiset, the script enumerates the integer counts

    J_b(v,w) = number of inputs with V=v and W=w.

The denominator is

    D_b = binomial(32,j) j!/product_a n_a! * product_i binomial(4,b_i),

where n_a is the multiplicity of a in b. The enumeration includes all
window subsets, all assignments of the multiset, and all lane masks.
It checks that the counts sum to D_b and that |v-sum_i b_i|<=w<=v+sum_i b_i.

For a preceding one-window weight a, let A_a(s) count the one-window inputs
with feedback s. It has total D_a=32 binomial(4,a). The current production
feedback map makes each A_a(s) either zero or one; this is checked, not
assumed by a generic lookup routine. The script additionally counts

    H_{b,a}(w) = sum_{X: W=w} A_a(BX).

These joint counts retain the correlation between cancellation and output
weight that separate marginals discard.

## Coefficients in the existing envelope

Let z=exp(-lambda), with lambda>=0, and let N_v be the number of nonzero
states whose expansion has weight v. For a uniform source on this class,
the lazy return coefficient is exactly

    alpha/(D_b N_v) * sum_w J_b(v,w) z^w.

For the fresh source distribution A_a(s)/D_a, it is exactly

    alpha/(D_b D_a) * sum_w H_{b,a}(w) z^w.

The envelope takes the maximum over the four possible fresh-source weights.
For a general nonzero source with pointwise density at most c, the lazy
return contribution is at most c times

    alpha/D_b * sum_{v>0,w} J_b(v,w) z^w.

The last expression remains an upper bound, not an exact return moment for
an unspecified entering distribution. Excluding v=0 is valid because the
expansion map is injective and this source excludes state zero.

The fresh and uniform coordinates also have a refresh contribution to
zero. The implementation adds a valid upper obtained from their existing
uniform-class output coefficients, divided by the class sizes. It does not
discard refresh returns or subtract a loose lazy bound. All shape penalties
are applied before maximizing over input multisets.

## Verification and use

There are 4, 10, and 20 multisets at one, two, and three active windows.
All counts are exact integers. The total enumeration size bounds every
signed-64-bit accumulator; no weighted floating histogram is used.
Independent checks include:

- Direct scalar enumeration at selected window subsets for each occupancy.
- Every output-weight marginal against the previous cancellation census.
- 136 fresh-cancellation marginals against the overlap/kernel identities.
- 34 feedback marginals against the Walsh-distribution census when that
  census is enabled in the driver.

Arb evaluates the resulting positive moment sums outward. The standalone
calculation has been run at 192 and 384 bits of precision.

    python -B research/workstreams/permutation_locality/cancellation_joint.py
    python -B research/workstreams/permutation_locality/cancellation_joint.py --precision 384

The optional driver flag is `--joint-cancellation`. It leaves higher local
occupancies unchanged and leaves the default proof path disabled for this
new refinement. It changes neither production code nor the paper.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --window-histogram 8 --joint-cancellation --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties .75 --tilts .032 --probe-supports 128 144 160

This command evaluates selected points in binary64, not a complete outward
support cover. No full-code certificate follows from the local enumeration.

At tilt .032 and all-one penalty .75, 16 zero-return coefficients improve.
The resulting selected-point log2 upper bounds are:

| Equal support in all 64 groups | Before joint returns | After joint returns |
|---|---:|---:|
| 128 | +2824.859997 | +2821.181478 |
| 144 | +2391.926914 | +2388.230135 |
| 160 | +2240.114367 | +2236.411804 |

The gain is only about 3.7 bits. The earlier general-density cancellation
bound already used the exact output moment at two and three windows.
The new information mainly tightens fresh and uniform-class returns.
These results do not support extending brute-force enumeration to four
windows solely to close the current gap.

## Retuning after the count improvements

Keeping the same refinements, a separate selected-point search used tilt
.048 and penalties .75, .875, and .95. All four points selected penalty .95:

| Equal support in all 64 groups | Log2 upper at tilt .048, penalty .95 |
|---|---:|
| 128 | +3210.405477 |
| 144 | +1765.699104 |
| 160 | +745.253448 |
| 176 | +258.937318 |

These are binary64 diagnostics, not outward certificates. Retuning improves
the larger supports substantially; support 128 still prefers the earlier
witness among these trials. In particular, the 3.7-bit local refinement
gain above does not measure the remaining parameter slack. Because every
new point selected the largest tested penalty, the next bounded check used
penalty 1 and tilts .048, .056, and .064.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --window-histogram 8 --joint-cancellation --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties .75 .875 .95 --tilts .048 --probe-supports 128 144 160 176

The penalty-1 check completed. Every point selected tilt .048 among the
three tested tilts:

| Equal support in all 64 groups | Log2 upper at penalty 1 |
|---|---:|
| 144 | +2083.582650 |
| 160 | +1042.818779 |
| 176 | +609.371772 |
| 192 | +950.794469 |

At the three common supports, these bounds are worse than the penalty-.95
bounds. Thus removing the all-one penalty is not an improvement, and
larger tilts at penalty 1 do not fix these points. This does not rule out
other joint parameter choices. The next focused search should bracket
the useful interior penalty and tilt, retaining separate witnesses for
different supports rather than imposing one global choice.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --window-histogram 8 --joint-cancellation --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties 1 --tilts .048 .056 .064 --probe-supports 144 160 176 192

An interior grid then tested all nine pairs from tilts {.04,.048,.056}
and penalties {.9,.95,.975}. Joint conditioning was optimized for all nine
pairs (`--joint-top 9`), rather than only the two best initial proposals.

| Support | Best log2 upper in this grid | Tilt | Penalty |
|---|---:|---:|---:|
| 112 | +2747.450728 | .04 | .9 |
| 128 | +2461.339977 | .04 | .95 |
| 144 | +1293.226734 | .04 | .975 |
| 160 | +555.460137 | .04 | .975 |
| 176 | +258.937318 | .048 | .95 |
| 192 | +549.506575 | .048 | .95 |

The smaller tilt helps supports 128, 144, and 160. None of the six points
closes, and support 176 reproduces its previous best witness. Holding each
new witness fixed, an outer-count reduction of (score+40)/64 bits per group
would give a 40-bit bound for that selected support vector. For supports
128, 160, and 176, this requires about 39.08, 9.30, and 4.67 bits per group.
These are refinement targets, not available count improvements or full-cover
claims. All scores in this table remain binary64 diagnostics.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --window-histogram 8 --joint-cancellation --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --joint-top 9 --penalties .9 .95 .975 --tilts .04 .048 .056 --probe-supports 112 128 144 160 176 192

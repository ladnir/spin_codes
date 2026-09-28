# Remember low expansion weights within the mature component

The mature component previously paid for expansion weight 48 at every step.
Most states do not have that weight, and a new random window often leaves
the low-weight classes. Two extra upper bounds retain enough information
to exploit this without claiming that the mature state is uniform.

## Coordinates and output moment

Let mu be the nonnegative, output-weight-tilted measure assigned to the
mature component on nonzero states. The existing coordinates bound its
total mass M and its largest point mass C. Add L48 and L56 satisfying

    sum_q mu(q) <= M,                 mu(q) <= C,
    sum_{wt(Eq)<=48} mu(q) <= L48,
    sum_{wt(Eq)<=56} mu(q) <= L56.

The tails overlap: L48 is a subset of L56, and both are subsets of M.
They are auxiliary upper bounds, not additional state components. Their
terminal coefficients are zero, just like C. The terminal mass is still
Z+F+M+sum_v U_v.

Fix a current column-weight shape and let g(v) upper-bound its output
moment for every entering state whose expansion weight is v. The only
nonzero expansion weights are 48, 56, 64, 72, and 80. The multi-window
formula in WINDOW_AVERAGES.md supplies g. Define

    a = max(g(64),g(72),g(80)),
    b = max(0,g(56)-a),
    c = max(0,g(48)-a-b).

Pointwise, g(wt(Eq)) <= a+b[wt(Eq)<=56]+c[wt(Eq)<=48]. Consequently,

    sum_q mu(q) E_x[z^wt(Eq+x)] <= a M+b L56+c L48.

The lazy branch contributes at most half this expression to outgoing
mature mass. The refresh branch contributes at most this expression
divided by 2(2^19-1) to zero mass, and N_v times that amount to U_v.
For an empty input the zero contribution is exactly zero. Each shape's
all-one penalty is applied before maximizing its three coefficients.

The lower-complexity cutoff-56 variant uses a=max_{v>=56}g(v), b=0, and
c=max(0,g(48)-a). It still retains both tails, but only L48 affects the
total output moment.

## Propagate the tails

An empty lazy step leaves the state unchanged. Thus its tail contributions
are bounded by

    L48' <= z^48 L48/2,
    L56' <= (z^56 L56+(z^48-z^56)L48)/2.

Other old components do not enter mature mass in an empty step. Refresh
mass remains in the existing U_v components and is not also put in M.

For one occupied window of weight b, enumerate every entering state and
every uniform input of that weight. Record the probability that
q+Bx is nonzero and has expansion weight at most 48 or 56. Use its maximum
over nonzero q for an arbitrary mature input. For a fresh component,
average over each p_a and maximize over a. For U_v, average over its
expansion-weight class. Multiplying these probabilities by the corresponding
pointwise output factor and 1/2 gives outgoing tail bounds. Intersect each
with the old bound on that source's entire lazy mature mass.

For two or more occupied windows, the current implementation allows all
outgoing mature mass to enter either tail. This is conservative and leaves
a possible later refinement. The old density bound C is preserved; it
does not depend on whether the new mass estimates are tight.

## Exact census and checks

`mature_tail.py` enumerates all 2^19 entering states and the 480 nonzero
single-window inputs. Integer counts determine every tail probability.
It independently replays the maximizing state and selected other states.
The maximum probabilities, conditioned on column weight b, are:

| b | Nonzero next expansion weight <=48 | Nonzero next expansion weight <=56 |
|---|---:|---:|
| 1 | 9/128 | 53/128 |
| 2 | 5/96 | 35/96 |
| 3 | 9/128 | 51/128 |
| 4 | 5/32 | 19/32 |

The transfer check directly computes all component contributions for
representative states of every expansion weight, empty inputs, and each
single-window weight. It checks 330 coordinate inequalities per operator
choice, with higher-precision reference exponentials. Independent exact
small-population tests for the multi-window moment remain in force.

## Initial comparison and scope

At thirteen groups, each support 152, lambda=.008, and all-one penalty
.6875, the previous nine-coordinate estimate has log2 upper -39.44.
The cutoff-56 and cutoff-64 tail models give -100.00 and -138.91. All use
the same prefix-count, fresh-cancellation, and cancellation-output
refinements. These initial scores use binary64 matrix powers and do not
establish full support coverage.

    python -B research/workstreams/permutation_locality/mature_tail.py

`occupancy_cdf_cover.py --mature-tail 64` uses the eleven-coordinate model
and excludes both new coordinates from outward terminal evaluation.
The permutation distribution and encoder implementation are unchanged.

## Complete thirteen-group result

The first full cover used only tilts .0064, .008, and .01. Its remaining
largest boxes had small supports, including [38,41] in every group. Adding
lower tilts resolved this witness-selection gap. The resulting 366-leaf
cover, after 50 splits, includes all supports, ranks, and group locations
at exactly thirteen active groups. Its failure upper is below
2.266948e-23, or **75.2235952311 bits**, at both 192- and 384-bit precision.
The cumulative upper through thirteen groups remains below 3.188222e-15,
about 48.1562 bits. This is not a full-code certificate; fourteen or more
active groups remain outside the combined result.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --max-splits 150 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --retain-parents --joint-witness --penalties .5 .625 .75 1 --tilts .0032 .004 .005 .0064 .008 .01
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --precision 384 --max-splits 75 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --retain-parents --joint-witness --penalties .5 .625 .75 1 --tilts .0032 .004 .005 .0064 .008 .01

## Pair-conditioned tail probabilities

The next refinement no longer assigns every multi-window mature output to
both low-weight tails. For two current column weights a,b, let D_ab(s)
count ordered inputs in distinct windows whose feedback XOR is s. For
cutoff c in {48,56}, let f_c(s) indicate that s is nonzero and wt(Es)<=c.
The XOR convolution D_ab*f_c, evaluated at q, counts inputs taking q into
that tail. Its maximum over all q, divided by the number of inputs, gives
a uniform probability bound P_ab,c.

`pair_tail.py` computes these counts by integer Walsh--Hadamard transforms.
First convolve the unrestricted single-window counts, then subtract
same-window pairs using the mask-completion identity in THIRTEEN_GROUPS.md.
The inverse transform is exactly divisible by 2^19. Every intermediate is
bounded in magnitude by 2^19 times the input count times the tail size,
which is below 2^63. The implementation checks divisibility, nonnegativity,
totals, and equality with a direct enumeration of every distinct-window
input pair. It also directly sums the convolution at selected states and
its maximizing state. No floating-point transform is used.

For j>=2 active windows, fix all but two inputs and their positions. There
are 34-j windows left for the chosen pair. Its numerator is at most the
corresponding count across all 32 windows, even after the other feedback
has shifted q. Its tail probability is therefore at most

    min(1, P_ab,c * 32*31 / ((34-j)*(33-j))).

Choose the best pair of weights present in the fixed shape. This bound
holds conditionally for every other input, so averaging preserves it.
Multiply by the pointwise output factor and maximize over shapes with
their all-one penalties. Fresh, mature, and uniform sources use the lazy
factor 1/2; a zero source has no such factor. For a zero source at j=2,
use the exact convolution at zero instead of the maximum over q.

At fourteen groups, support 160, lambda=.008, and penalty .75, the original
nine-coordinate score is +95.38. Single-window tails reduce it to -22.79;
the pair-conditioned refinement reduces it to -58.31. This is a selected-
point binary64 screen, not a fourteen-group certificate. The pair refinement
is not used by the thirteen-group certificate above.

    python -B research/workstreams/permutation_locality/mature_tail.py --groups 14 --support 160 --tilt .008 --penalty .75 --pair-tail

The next step is a full fourteen-group cover with `--pair-tail`, including
lower tilts for small supports. Encoder performance remains unmodified.

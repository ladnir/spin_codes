# An alternative bound on mature-state concentration

This is an experimental local refinement, not an implemented certificate.
The construction and probability space are those in [BRIDGE.md](BRIDGE.md).
The current universal verifier propagates the mature component's maximum
pointwise density through its unsplit density column. Exact feedback
probabilities also give a different upper bound using the mature mass.
These estimates need not be tight in the same regime.

Fix a nonempty packet-weight shape h, with total input weight W. Let X
have its uniform distinct-window and lane-mask distribution. Write

    beta_h = max_r Pr[B X = r],
    f_v = exp(-lambda |v-W|).

The maximum defining beta_h includes r = 0. The existing full-feedback
census supplies this value as max(zero count, nonzero peak)/denominator
through six occupied windows.

Let mu be the entering mature measure, supported on nonzero states and
independent of X. For a nonzero target t, its lazy-branch density is

    alpha * sum_s mu(s) E[1_{B X=t+s} exp(-lambda wt(E(s)+X))],

where alpha = 1/4 for two transvections. The triangle inequality gives
wt(E(s)+X) >= |wt(E(s))-W|. Thus this density is at most

    alpha * beta_h * sum_v f_v M_v,

where M_v is the mass of states whose expansion has weight v. This uses
the joint event and output weight inside the expectation. It does not
assume that feedback and output tilt are independent.

The existing envelope stores total mature mass M and the nested masses
L48 and L56 below expansion weights 48 and 56. Define

    a = max(f_64, f_72, f_80),
    b = max(0, f_56-a),
    c = max(0, f_48-a-b).

Then sum_v f_v M_v <= a M + b L56 + c L48. All three coefficients are
nonnegative, so this is compatible with the positive-matrix envelope.
For each shape, multiply the coefficients by alpha beta_h and by the
existing input-weight factors. Maximize over every shape of that local
occupancy. This produces an alternative mature-to-density column whose
sources are M, L56, and L48, instead of C.

## Safe integration boundary

This is an alternative bound on the **whole mature contribution** to the
density coordinate. It must not be combined with the existing C column
by independently minimizing matrix entries: that could zero both valid
alternatives and invalidate the bound.

A fixed convex combination of the two complete columns is valid. Its
mixing coefficient can depend on the proof witness and local occupancy,
but not on unknown entering coordinates. Keep the Z, F, and uniform-class
contributions unchanged. Keep occupancy zero unchanged. Apply this choice
after the refinements that require an unsplit density column.

The prototype is isolated in `candidates/mass_density.py`; the production
verifier does not import it. Nine tests cover exact local models,
mixture endpoints, shape completeness, zero feedback atoms, and the
conditioning bound below. The existing 94 tests also pass. An improvement
still requires outward replay and cache keys that identify these choices
before it can contribute to a certificate.

## First numerical comparison

The screen regenerated the complete current baseline at lambda = .056,
rho = .9, and input-weight tilt one. It includes the six-window feedback
and density censuses, eight-window output histograms, and joint cancellation.
The following binary64 scores concern homogeneous group union support 200.
They include all compatible row weights and all group locations, but not
the other support vectors at either occupancy.

| Mature-density choice | q = 64 log2 upper | q = 80 log2 upper |
|---|---:|---:|
| Existing density column | -124.705 | +982.349 |
| One-quarter mass mixture at all j = 1,...,6 | +808.899 | +1966.111 |
| Mass column only at j = 4 | -78.274 | +1042.760 |
| Mass column only at j = 5 | -128.721 | +989.458 |
| Mass column only at j = 6 | -140.000 | +974.810 |

Thus the six-window alternative gives about 15.30 bits at this q = 64
point, but only 7.54 bits at q = 80. Applying it to sparse epochs is
counterproductive. This is evidence for a selective refinement, not a
replacement of the baseline or a resolution of the dense gap. The
baseline support-200 score agrees with its prior outward replay within
4e-9 bits; the alternative scores have not been replayed outward.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_density_screen.py --tilt .056 --groups 64 80 --supports 192 200 208
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates -p 'test_*.py'
```

## Extending a short census by conditioning

Fix a shape h containing j labeled packets and a submultiset s containing
m packets. Expose the other j-m packets, including their window positions
and masks. The remaining m packets have their usual uniform injection
conditioned on avoiding those positions. This conditioning event has
probability binom(32-j+m,m)/binom(32,m). Exposed feedback only shifts the
target. Consequently,

    beta_h <= min(1, beta_s * binom(32,m)/binom(32-j+m,m)).

Minimizing over all available submultisets preserves the upper bound.
The prototype computes those minima by an exact-rational recurrence,
rather than enumerating every full-shape/subshape pair. This extends the
six-window census to all shapes through any chosen cutoff at most 32.
The same triangle-weight factors and coupled input penalties then apply.
The local tests compare this bound with complete toy feedback distributions
and independently enumerate every submultiset in a second check.

A selected-point screen through twelve windows gives the following scores
for union support 200 at the same lambda and rho. The mass column is used
at every local occupancy from the listed start through twelve.

| First local occupancy using mass bound | q = 64 log2 upper | q = 80 log2 upper |
|---|---:|---:|
| 5 | -149.253 | +976.550 |
| 6 | -153.854 | +967.665 |
| 7 | -143.622 | +973.183 |

The best tested choice gains about 29.15 bits at q = 64 and 14.68 bits
at q = 80, relative to the unchanged baseline. It does not close the
q = 80 point or the full q = 64 support sum. No bound from this extension
has entered the complete occupancy aggregate.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_density_screen.py --tilt .056 --groups 64 80 --supports 192 200 208 --condition-through 12
```

## The same estimate for a return to zero

The local mass estimate also holds at target t = 0: the entering mature
measure still has no mass at source zero. Thus the same coefficients can
bound its lazy return to zero. Keeping beta_h's zero feedback atom is
conservative here.

In the current baseline, C->Z bounds the lazy mature return, while
M/L48/L56->Z bound its refresh contribution. A convex replacement must
retain those refresh coefficients, add the new lazy coefficients there,
and scale only C->Z by one minus the chosen fraction. Contributions from
the zero, fresh, and uniform components remain unchanged. This argument
uses that component decomposition; it is not a transformation for arbitrary
nonnegative matrices.

The prototype exposes this choice separately through `--zero-returns`.
The exact toy check includes target zero, and a matrix test checks that
the refresh coefficients survive unchanged. Its numerical screen gives
the following binary64 results for homogeneous union support 200.

| Local occupancies using the alternative return bound | q = 64 log2 upper | q = 80 log2 upper |
|---|---:|---:|
| 6 only | -139.726 | +972.967 |
| 7 only | -149.790 | +964.118 |
| 5 through 12 | -169.547 | +958.734 |
| 6 through 12 | -177.624 | +945.017 |
| 7 through 12 | -164.540 | +953.884 |

Replacing returns at occupancies six through twelve improves this q = 64
point by about 52.92 bits and the q = 80 point by 37.33 bits. Sparse
replacements again lose substantially. These gains do not certify the
full q = 64 event or resolve the q = 80 point. No outward replay or
complete cover uses the prototype yet.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_density_screen.py --tilt .056 --groups 64 80 --supports 192 200 208 --condition-through 12 --zero-returns
```

## Combining the density and zero-return alternatives

The two alternatives bound different outgoing coordinates of the same
mature measure. They may therefore be applied together. Each coordinate
still satisfies its own inequality; no independence between those bounds
is needed. The implementation replaces the zero-return column first and
then the density column, retaining the original refresh contributions.
A regression test checks both resulting columns against separate applications.

At lambda = .056, rho = .9, and homogeneous support 200, applying both
alternatives at local occupancies six through twelve gives these binary64
scores:

| Choice | q = 64 | q = 80 |
|---|---:|---:|
| Baseline | -124.705335 | +982.349100 |
| Zero return only | -177.624311 | +945.016890 |
| Zero return and density | -192.844929 | +936.015098 |

The combined improvement is about 68.14 bits at the selected q = 64 point
and 46.33 bits at q = 80. It does not resolve the latter point.
`candidates/mass_verify.py` regenerates the alternative with outward
arithmetic. Its default mode replays selected support vectors; its
`--full-cover` mode retains eighteen baseline witnesses, including the
finer tilts .054 and .058, and adds the alternative to the existing
complete shell/CDF cover. Experimental
operators are not loaded from or written to the production operator cache.
No complete occupancy result has used this alternative yet.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 64 --supports 192 200 208
```

This selected-point run passes at 192-bit precision. Its log2 uppers,
rounded upward, are -189.838543 at support 192, -192.844929 at support
200, and -300.691502 at support 208. These include every compatible
row-weight pattern and group location, but only the stated support vectors.
At support 192 the best prior tilt .052 already gave -162.333102;
the improvement over that stronger reference is about 27.51 bits, not
the 75.47-bit difference from the matched .056 baseline. The local
candidate tests and the unchanged production suite pass 14 and 94 tests,
respectively.

A subsequent finer baseline grid improves support 200 to -140.541282
using tilt .054. Against that best tested baseline, the combined mass
bound improves this point by about 52.30 bits. A full q = 64 cover is
now running with all eighteen baseline witnesses and the alternative:

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 64 --full-cover --max-splits 400
```

## Remaining transition sensitivity

The companion `candidates/transition_screen.py` also deletes selected
entries to measure their influence. These deletions are invalid as
probability bounds. At the same support-200 point, they give:

| Deleted active-epoch terms | q = 64 score | q = 80 score |
|---|---:|---:|
| Fresh-state return to zero | -187.483 | +926.414 |
| Fresh-state density contribution | -166.117 | +950.612 |
| Both fresh-state terms | -227.955 | +894.894 |
| Mature-density contribution to zero | -368.169 | +768.941 |
| Mature-density self-transition | -240.729 | +889.310 |
| Uniform-class returns to zero | -329.104 | +762.405 |

Thus improving only the fresh return and density coefficients cannot
close this q = 80 witness: even deleting both terms leaves a positive
score. This does not rule out a richer state decomposition that improves
other transitions too. The separate deletions are not additive estimates
of possible gains. The unchanged baseline still bounds every packet shape;
the screen does not identify an actual low-weight codeword.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/transition_screen.py
```

## Retuning and local sensitivity

The benefit grows at tilt .064. At q = 80 and support 200, the
baseline binary64 score is +924.579633. Replacing both concentration
columns at local occupancies six through twelve gives +589.648409,
an improvement of about 334.93 bits. The latter value passes a
192-bit outward replay, but remains vacuous as a probability bound.

A grid over tilts .06, .064, .068, .072, and .08 selects .064 for
supports 192 and 200, and .068 for support 208. The corresponding
outward log2 uppers, rounded upward, are +505.161937, +589.648409,
and +533.308303. None closes its selected class.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 80 --supports 192 200 208 --tilts .06 .064 .068 .072 .08
```

`candidates/local_sensitivity.py` differentiates the finite positive-matrix
calculation rather than replacing it by a limiting eigenvalue. It
propagates an adjoint through all 256 regions and the normalized
without-replacement recurrence for 64 epochs per region. For a local
coefficient T_j(s,t), it reports

    T_j(s,t) * partial log(moment) / partial T_j(s,t).

These elasticities describe the envelope calculation, not transition
frequencies in the actual encoder. Their sum must equal 256*64 = 16384,
because the moment is homogeneous of that degree in all local matrices.
Both production diagnostic runs satisfy this identity. Small
noncommuting examples also match finite-difference derivatives.

At the .064 baseline point, the C-to-C and C-to-Z elasticities sum to
about 659.88 and 458.86, respectively. The largest nonempty-epoch terms
occur at local occupancies roughly four through eight. For example,
C-to-Z contributes 83.27 at occupancy four and 86.83 at occupancy seven,
but only 0.16 at occupancy thirteen. This supports improving the middle
of the local census before extending its very high-occupancy tail.

The experimental `candidates/mass_optimize.py` can use these derivatives
to choose fixed convex fractions for the two alternative columns at each
local occupancy. It rounds the proposed fractions to rationals and
reconstructs the coupled columns outward before selected-point replay.
At q = 80, support 200, tilt .064, and rho = .9, two optimization
passes select fraction one for both columns at occupancies five through
sixteen, and zero elsewhere. The 192-bit outward log2 upper is
+584.204583, rounded upward. This gains about 5.44 bits over the
six-through-twelve choice, but still does not close the selected event.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_optimize.py --groups 80 --support 200 --tilt .064 --maximum 16 --iterations 30 --passes 2
```

The [compressed feedback extension](SPECTRAL_FEEDBACK.md) gives a much
larger improvement at these points. Combining that extension with this
optimizer gives only another 0.61 bits at support 200. Neither selected-point optimization nor
its rational outward replay is a complete support cover.

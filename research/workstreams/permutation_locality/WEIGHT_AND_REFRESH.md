# Joint weight counts and refresh sensitivity

The two-update, four-row candidate remains the proof target. These tests
ask where its intermediate-occupancy bound loses strength. They do not
change the encoder or extend the complete occupancy certificates.

## Coupling two outer weights

For a nonzero ordered tuple of four outer words, let u be its union support,
W its total binary weight, and J its number of all-one columns. The inner
operators can carry the factor a^W rho^J, where a >= 1 and 0 < rho <= 1.
The matching outer sum weights each tuple by a^-W rho^-J.

Previously, the implementation supported a > 1 only when rho = 1.
`joint_weight_caps.py` bounds their joint weighted sum using two marginal
count bounds. The fourth power of the ordinary outer enumerator gives a
CDF upper bound for W. The existing all-one-column flags bound the number
of tuples at each J. Both bounds are restricted to union support at most u.

First intersect their total masses with the existing support-count bound.
Then pack the W mass toward small values and the J mass toward large
values. Pair these two extremal marginals in decreasing order of their
weights. The rearrangement inequality makes this pairing an upper bound
on every feasible joint sum. All counts and marginal allocations are exact
integers; the final powers use outward Arb arithmetic and integer ceilings.

This is a relaxation, not a joint enumerator. In particular, it does not
enforce W >= 4J or the weights of all fifteen nonzero row combinations.
The local operators continue to retain their actual shape-specific weights
before taking coordinatewise maxima.

The tests check 4,374 exact rational inequalities for small count tables,
including inflated marginal caps. Another 51 inequalities check the full
outer-count pipeline against exhaustive enumeration of three small codes.

    python -B research/workstreams/permutation_locality/joint_weight_caps.py

The first trial used a = 1.03, two updates, containment moments, and full
feedback distributions through six active windows. The following matched
points all use 64 active groups, lambda = 0.032, and rho = 0.75. Each group
has the displayed support. Values are binary64 log2 contribution bounds,
including group locations, rather than outward certificates.

| Support | a = 1 | a = 1.03 |
|---|---:|---:|
| 80 | +801.6342 | +1008.3468 |
| 96 | +1847.0632 | +2052.3222 |
| 128 | +3380.7392 | +3629.9809 |

This trial does not improve the difficult points. It does not rule out
other tilts or a stronger joint count. At support 160 the trial selected
lambda = 0.048 instead, so that result is not a matched comparison with
the earlier lambda = 0.032 run.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --full-feedback 6 --weight-tilt 1.03 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --joint-witness --penalties .5 .75 --tilts .032 .048 --probe-supports 80 96 128 160

## Would more inner refresh resolve this bound?

For a fixed nonzero state, r independent transvection updates have a lazy
coefficient 2^-r and a uniform-refresh coefficient 1 - 2^-r. The diagnostic
rescales the two-update operators to several values of r. It also tests
the limiting distribution with no lazy contribution. This limit still
leaves the zero state fixed before feedback; it is not unconditional
replacement of every state by an independent random state.

The rescaling separates lazy and refresh contributions. Empty-epoch
uniform-class self-loops require subtracting their known lazy term before
rescaling. For active-epoch returns to zero, fresh and uniform source
coordinates use the existing refresh upper bound separately. The full
feedback refinements are then applied at the new value of r, not rescaled
after taking their minima. The r = 2 identity is checked exactly on the
cached binary64 arrays.

The cache is a diagnostic input, not an outward-arithmetic proof input.
The following points use 64 groups, support 128 in each group, lambda =
0.032, rho = 0.75, and optimized Bernoulli witness p. Neither the outer
count nor the tilts are reoptimized for the alternative ensembles.

| Updates | Baseline local bounds | Full feedback through six windows |
|---|---:|---:|
| 2 | +3809.0603 | +3380.7392 |
| 3 | +3354.3584 | +3210.1475 |
| 4 | +3219.3017 | +3146.4157 |
| 6 | +3134.2353 | +3106.6099 |
| Uniform-refresh limit | +3108.0285 | +3094.7639 |

Thus this bound still fails even in the limit. This observation concerns
the tested envelope and parameters; it does not prove that the limiting
ensemble lacks distance, or that another analysis cannot certify it.
Timings of the two-update implementation do not apply to more updates.

The ablation tests remove selected positive matrix entries. These modified
matrices are **not valid upper bounds**. They identify which terms dominate
this diagnostic, not which events the proof may discard.

| Deleted active-epoch returns to zero | Two updates | Uniform-refresh limit |
|---|---:|---:|
| None | +3380.7392 | +3094.7639 |
| From the zero source | +3379.5793 | +3093.5715 |
| From all nonzero sources | -4665.0768 | -4825.8599 |
| From uniform-class sources only | +3205.5898 | -4784.7351 |

Direct zero-feedback events from the zero source no longer explain the
gap. Returns from a nonzero state remain important. Uniform refresh does
not remove these events: the refreshed state can equal the incoming
feedback, producing zero after addition. Their probabilities must remain
in a valid proof.

Holding the refined two-update inner bound fixed, the support-128 point
would reach -40 only after reducing its weighted outer CDF cap by about
53.45 bits per group. This is a sensitivity threshold, not a proved count
improvement. Reaching -40 at this one point would still not cover all
support vectors or sum all occupancy contributions.

Recreate the baseline cache locally, then run the diagnostics:

    python -B research/workstreams/permutation_locality/occupancy_sensitivity.py --groups 64 --support 128 --tilt .032 --penalty .75 --snapshot tmp/r2-q64-sensitivity.npz
    python -B research/workstreams/permutation_locality/round_sensitivity.py
    python -B research/workstreams/permutation_locality/round_sensitivity.py --full-feedback 6 --ablate

The next priority is to couple the outer column patterns more tightly to
the inner cancellation analysis, or substantially improve their count.
Increasing the update count alone is not supported as the next step by
these tests. The measured two-update design and its partial certificates
are unchanged.

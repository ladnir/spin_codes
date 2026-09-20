# First local audit

The first audit asks whether a 64-bit state can produce a strong 256-bit live
output without a dense binary matrix.  It does not include the outer code or
the complete two-state recurrence.

## Expansion results

The table reports the Chernoff suppression obtained when all 8,192 epochs
start with a live randomized state.  The global weight threshold is
`floor(0.09*2^21)`.

| Expansion ensemble | XOR work per output | Suppression (bits) |
|:---|---:|---:|
| Four interleaved 64-bit accumulators | about 1 | 443,251 |
| Repeat-4 then one accumulator | 1 | 433,073 |
| Repeat-4 then two accumulators | 2 | 491,249 |
| Repeat-4 then three accumulators | 3 | 518,325 |
| Random sparse rows of weight 5 | 4 | 436,388 |
| Random sparse rows of weight 15 | 14 | 687,980 |
| Random sparse rows of weight 25 | 24 | 1,166,172 |
| Dense random 64-dimensional subspace | dense reference | 1,176,755 |

The dense outer envelope at full occupation is about 1,056,768 bits.  The
dense random reference clears that value.  The cheap accumulator and LDGM
surrogates do not.  Thus the 64-bit state is not the obstruction.  The
obstruction is obtaining a random-code-like 64-to-256 expansion with a cheap
encoder.

Even row weights are unsuitable for the sampled sparse generator ensemble.
The all-one state lies in the kernel of every even-weight output row.

## Compression results

The audit sampled 64-by-256 matrices with distinct fixed-weight columns.  It
counted kernel supports through weight four exactly.

| Column weight | Best rank | Weight-3 words | Weight-4 words |
|---:|---:|---:|---:|
| 3 | 64 | 0 | 3 |
| 5 | 64 | 0 | 0 |

Sixteen matrices were sampled at each degree.  The degree-5 result does not
prove a minimum distance above four.  It does show that sparse compression
need not create immediate weight-two or weight-four cancellations.

Irregular degree-2/3 matrices reduce the syndrome cost but introduce a small
explicit exceptional set:

| Average column weight | Weight-3 words | Weight-4 words |
|---:|---:|---:|
| 2.125 | 45 | 274 |
| 2.250 | 28 | 148 |
| 2.500 | 8 | 50 |

These multiplicities are not automatically fatal.  A proof can isolate the
exceptional kernel words and charge the outer permutations for aligning an
epoch input with one of them.  Average degree 2.5 is the current performance
candidate; uniform degree 3 is the cleaner proof candidate.

## Revised constituent target

A nested-code realization appears more promising than a naive LDGM map.
Choose binary codes

\[
C_1\subseteq C_0\subseteq\mathbb F_2^{256}
\]

with dimensions 64 and 192.  Let `B` be a sparse parity-check map for `C_0`,
and let `A` encode the 64-bit state into `C_1`.  Then `BA=0` follows from the
nesting.

The target combines:

1. an LDPC-like sparse syndrome for `C_0`;
2. a random-code-like spectrum for `C_1`;
3. a structured linear-time encoder for `C_1`.

## Fixed nested-pair result

The first fixed realization uses coordinates `(q,a,p)` of widths
`(64,128,64)` and the block-triangular encoder in `CONSTRUCTION.md`.  The
table compares the useful nearby choices.  The minimum is exhaustive over
all state inputs of weight at most four; it is not the global code distance.

| Auxiliary transform | `P` column weight | Minimum found | Logical `A` XOR/output |
|:---|---:|---:|---:|
| one accumulator | 6 | 42 | 3.496 |
| two accumulators with an interleaver | 6 | 62 | 3.992 |
| three accumulators with interleavers | 6 | 66 | 4.488 |
| two accumulators with an interleaver | 8 | 64 | 4.492 |

The depth-two, weight-six point was the first useful performance knee.  For
its selected fixed instance:

- `BA=0` verifies exactly;
- `rank(A)=64` and `rank(B)=64`;
- the exhaustive input-weight-at-most-four audit contains no output word of
  weight at most 55;
- the minimum in that audit is 62, attained by an input of weight four;
- 200,000 uniformly sampled 64-bit inputs have mean output weight 127.986
  and minimum output weight 90;
- 20,000 unrestricted greedy-search restarts found local minima no lower than
  88, so the 62-bit sparse-input witness remains the best word found;
- `ker(B)` has no words of weight two or three and has 293 words of weight
  four.

The logical count for this constituent is 3.992 XOR/output to form `A(q)`
and 2.250 XOR/input to compute `B(X)`, before applying `A(q)` to the input,
the field state randomizer, circuit sharing, fusion, or SIMD effects.

This result changes the accumulator conclusion above.  Repeated accumulation
alone is too weak as a 64-to-256 expansion, but a second accumulator inside
the systematic nested encoder removes the observed low-support defect at a
cost of only 0.496 XOR/output.  A third layer gives only four more bits in
the audited floor for the same incremental cost.

The generic MILP audit did not certify the global distance: after 120 seconds
its dual lower bound was only 3.  This is evidence against generic
branch-and-bound as a proof tool, not evidence against the code.

## Exact ensemble spectrum

The layered encoder admits an exact expected-spectrum calculation.  For an
information word of weight `u`, XORing the sampled fixed-weight columns of
`P` is a Markov chain on the current weight in 128 coordinates.  The random
interleaver makes each following accumulator depend only on that weight.
The final `S/R` parity is a second fixed-weight-column chain in 64
coordinates.  Combining these chains gives the expected weight enumerator of
the complete `[256,64]` image of `A`; no sampling is used in this calculation.

For the depth-two, `P`-degree-6, `S/R`-degree-3 ensemble, the expected numbers
of nonzero codewords below selected thresholds are:

| Maximum weight | Expected words |
|---:|---:|
| 23 | 0.000196 |
| 31 | 0.000795 |
| 39 | 0.003187 |
| 47 | 0.033830 |
| 55 | 10.699 |

The unexpurgated expected moment gives only 566,427 bits of suppression in
the idealized all-live recurrence.  Rare sampled codes with very low-weight
words dominate that average.  Expurgation removes this artifact.

All 256 columns of `B=[S R I]` have odd weight.  Therefore `ker(B)` has no
words of weight one or three.  If the 192 degree-3 `S/R` columns are distinct,
it also has no words of weight two.  This distinctness event has probability
about 0.6435 in the with-replacement ensemble.  Combining it with Markov
expurgation gives the following exact existence bounds:

| Expurgate `im(A)` through | Joint-event probability lower bound | Live suppression (bits) |
|---:|---:|---:|
| 23 | 0.643346 | 668,984 |
| 31 | 0.642747 | 983,985 |
| 39 | 0.640355 | 1,294,555 |
| 47 | 0.609712 | 1,579,400 |

For example, the row at 39 proves that there exists a fixed nested pair for
which `im(A)` has minimum distance at least 40, `ker(B)` has minimum distance
at least four, and the live-state Chernoff moment gives 1,294,555 bits of
suppression at the global 9% threshold.  This exceeds the earlier dense-outer
full-occupation envelope of about 1,056,768 bits by about 238,000 bits.

This remains a constituent result.  The complete recurrence must still
account for zero-state epochs, activation by `B(X)`, the rare termination
transition, outer multiplicities, and placement.  The receipt is
`receipts/nested_depth2_ensemble_spectrum.json`.

## Zero-state activation

Write an epoch input as `(x_m,x_p)`, where `x_m` contains the first three
packet lanes and `x_p` contains the parity lane.  Put
`m=wt(x_m)` and `v=wt(x_p)`.  For independently sampled degree-3 columns,
the probability that `B(X)=0` is

\[
\Pr[B(X)=0]
=\frac{\Pr[\operatorname{wt}(Hx_m)=v]}{\binom{64}{v}}.
\]

The numerator is an exact 65-state weight chain.  Conditioning the selected
columns to be distinct gives a valid upper bound by dividing by their
distinctness probability.

Every column of `B` has odd weight.  Odd-weight epoch inputs therefore
activate deterministically.  With distinct columns, weights one, two, and
three also activate deterministically.  The first possible failure has
weight four.  Its worst split is `(m,v)=(1,3)`, with probability
`1/C(64,3)`, or about `2^-15.35`.  Under a uniform 256-position support, the
average weight-four failure probability is about `2^-19.16`.  The uniform
support reference falls below `5e-10` at weight eight and approaches about
`2^-63` at high weight.

The table is saved in `receipts/zero_state_activation_table.json`.  It does
not multiply probabilities across epochs that share one fixed `B`.

## Packet-kernel wiring

The parity-lane permutation strengthens the fixed compressor without adding
arithmetic.  For the selected `B`:

- every four-lane packet map has rank four;
- no input supported on one, two, or three packets belongs to `ker(B)`;
- the exact minimum nonactivation packet support is four;
- all 293 coordinate-weight-four kernel words use four distinct packets.

For four labeled singleton packets, the worst lane composition hits a kernel
word with probability at most `2^-30.4108` under the uniform permutation of
2048 packet groups.  The exact audit is saved in
`receipts/fixed_nested_pair_depth2_packet_kernel.json`.

At the 210-region late-start prefix, the modeled weight-38 shell gives the
following margins before the post-activation output tail:

| Active outer blocks | Prefix nonactivation bits | Modeled multiplicity bits | Margin |
|---:|---:|---:|---:|
| 1 | 123.318 | 37.278 | 86.041 |
| 2 | 246.637 | 73.555 | 173.082 |
| 3 | 369.955 | 109.247 | 260.708 |
| 4 | 493.274 | 144.524 | 348.750 |

For one through three active blocks, the prefix must contain no active
region.  For four blocks, a nonempty prefix region must be a four-way
intersection whose four singleton packets hit a kernel pattern.  The latter
events are negligible here, so the empty-prefix term remains dominant.  See
`receipts/sparse_prefix_w38.json`.

## One-active recurrence

Using the conditional-expectation spectrum directly loses 0.635 bits in
every live epoch because it conditions on distinct `B` columns.  That proof
artifact makes the one-active recurrence fail by 865 bits.  The actual
unscaled spectrum gives 81.563 bits of modeled margin.

The deterministic secant bound in `PROOF_PLAN.md` removes the conditioning
artifact.  It gives 79.671 bits of modeled one-active margin.  The optimum
has `log_surprisal=-7.8333`, and outer weight 38 dominates.  The recorded
fixed pair has no zero output coordinates, so it satisfies the full-support
hypothesis.  Its distance-40 hypothesis remains an ensemble existence claim,
not a certificate for those recorded matrices.

## Preferred degree-7 pair

The preferred proof instance raises the column weight of `P` from six to
seven and reuses the degree-6 compressor `B` and its parity-lane wiring.  The
extra edge costs 0.250 logical XOR per output bit.  It removes the duplicate
output-coordinate form forced by even column parity in the degree-6 map.

For the recorded fixed pair:

- `BA=0` verifies exactly;
- no dependency is supported on one, two, or three output coordinates of
  `A`, and exactly one dependency is supported on four coordinates;
- no nonzero input supported on at most three packets belongs to `ker(B)`;
- the exact packet-support distance of `ker(B)` is four;
- the exhaustive input-weight-at-most-four floor for `A` is 58;
- sampled uniform inputs had minimum output weight 92, while unrestricted
  hill climbing found 86; and
- the logical costs are 4.242 XOR/output for `A` and 2.250 XOR/input for `B`.

The fixed-pair and packet audits are in
`receipts/fixed_nested_pair_depth2_p7_fixedB_search.json` and
`receipts/fixed_nested_pair_depth2_p7_fixedB_packet_kernel.json`.

The exact degree-7 ensemble has expected nonzero codeword count
`0.000304870` through weight 39.  After requiring distinct columns in `B`
and expurgating through weight 39, the all-live suppression is
1,318,370 bits.  This is an existence calculation; it does not certify
distance 40 for the recorded fixed pair.

## Occupation ladder

The complete two-state recurrence was compressed to a two-by-two transfer.
For `a` active outer blocks, an exact coefficient dynamic program places
their active packets among the 8,192 epochs.  A pointwise regular-code
spectrum envelope covers every modeled even outer weight from 38 through
218.  The calculation assumes that the active blocks occupy distinct
four-block packet groups and excludes the all-one outer word.

Using only the mean of a live output gives these margins at the 9% target:

| Active outer blocks | Margin (bits) |
|---:|---:|
| 1 | 72.699 |
| 2 | 149.992 |
| 4 | 297.855 |
| 8 | 569.580 |
| 16 | 1,025.865 |
| 32 | 1,654.804 |
| 64 | 1,641.202 |
| 80 | 1,246.070 |
| 96 | 807.065 |
| 112 | 57.071 |
| 128 | -438.635 |

The apparent failure at occupation 128 is an artifact of retaining only the
mean.  The preferred degree-7 map has no dependency among at most three
output coordinates.  Consequently a uniform nonzero live output has the
same first three factorial moments as independent balanced coordinates, up
to conditioning away the zero state.  Maximizing its Chernoff moment over
all distributions on weights 40 through 256 with those moments is a small
linear program.  This three-moment bound gives 9,929.636 bits of margin at
occupation 128.  It is pointwise no weaker than the mean-only bound, so the
same recurrence is positive at every tested regular occupation from 1
through 128.

The occupation receipts are
`receipts/occupation_ladder_regular_1_16.json`,
`receipts/occupation_ladder_regular_32_64.json`,
`receipts/occupation_ladder_regular_80_112.json`, and
`receipts/occupation_ladder_regular_128_moment3.json`.

This is not yet an end-to-end certificate.  It is conditional on a
distance-40 constituent with the recorded three-coordinate independence
property, and it leaves shared packet groups, the all-one outer word, and
occupations above 128 open.

## Shared packet groups

Write the active-block profile as `(n1,n2,n3,n4)`.  Here `nr` counts fixed
packet groups that contain exactly `r` active outer blocks.  For a fixed
profile, a four-variable coefficient transfer averages the packet
permutation exactly.  The transfer also averages the Bernoulli outer
envelope before applying the two-state epoch matrix.

Every profile was evaluated exactly at occupations 2, 4, 8, and 16.  The
distinct-group profile was worst after charging its block-set multiplicity.
Its three-moment margins were 153.512, 301.040, 575.738, and 1,031.748 bits,
respectively.  Thus shared groups do not create a low-occupation failure.

The fully packed occupation-128 profile has 32 groups of width four.  At the
distinct-profile Chernoff point it has negative margin.  After optimizing its
own point, it has 7,105.481 bits of profile-specific margin.  Replacing its
profile multiplicity by all `C(8192,128)` active block sets would retain
6,393.067 bits if a packing domination lemma were available.  Representative
mixed profiles are also positive, but these calculations do not cover all
16,335 occupation-128 profiles.

A direct local packing proof does not work.  Suppose two separated packets
are merged.  The split execution can terminate the live state at the first
packet and reactivate it at the second packet.  The packed execution has no
corresponding state path.  An arbitrary intervening transfer therefore makes
the local split-to-packed ratio infinite.  This support mismatch is caused
by the `U`-to-`Z` transition, whose probability contains the `2^-64`
termination factor.  It does not provide refutation evidence against the
complete construction, but it blocks the proposed composable dominance
lemma.

The exact low-occupation receipt is
`receipts/shared_group_low_occupation_profiles.json`; occupation 16 is in
`receipts/shared_group_occupation16_profiles.json`.  The optimized packed
stress test is in `receipts/shared_group_occupation128_quad_optimized.json`.
The failed local compression audit is
`receipts/packing_domination_ladder_audit.json`.

## Termination phase

Multiplying each `U`-to-`Z` transfer by `exp(y)` marks termination
transitions.  The derivative at `y=0` gives the expected termination count
under the Chernoff-tilted path measure.  This is not the expected count under
the encoder's unweighted setup distribution.

At the optimized occupation-128 all-quad point, the no-termination inner
moment is 1,007.389 bits smaller than the full moment.  The tilted path
measure has expected termination count 123.234 and variance 2.130.
Termination is therefore the dominant renewal mechanism in that profile.

The transition is sharp along the singleton-to-quad curve:

| Active groups | Quad groups | Expected terminations | Moment inflation (bits) |
|---:|---:|---:|---:|
| 80 | 16 | 0.000071 | 0.000102 |
| 68 | 20 | 0.011 | 0.016 |
| 62 | 22 | 0.146 | 0.211 |
| 56 | 24 | 1.855 | 2.705 |
| 50 | 26 | 19.294 | 31.422 |
| 44 | 28 | 75.868 | 202.748 |
| 38 | 30 | 112.068 | 564.191 |
| 32 | 32 | 123.234 | 1,007.389 |

Thus a proof cannot treat termination as a rare perturbation uniformly over
profiles.  The receipt is
`receipts/termination_fugacity_singleton_quad_curve.json`.

## Region-level packing target

The pathwise packing comparison fails, but the complete one-region transfer
has the opposite empirical behavior.  The maximally packed region matrix
entrywise dominates:

- every profile at occupations 2, 4, 8, and 16;
- every occupation-128 profile with at least 30 quad groups; and
- the tested singleton, pair, triple, and mixed curves between those sets.

The occupation-128 near-packed scan contains all 15 profiles with at least
30 quad groups.  The all-quad profile is worst and retains 7,105.481 bits of
profile-specific margin.  The receipt is
`receipts/shared_group_occupation128_near_packed.json`.

This suggests applying packing only after all 32 epoch transfers in one
region have been averaged.  A region-level lemma would retain termination–
reactivation paths and would cover every shared-group profile without a case
enumeration.  The current evidence does not yet prove that lemma.

The same data also pass the stronger local test induced by the six elementary
integer-partition moves

`1+1->2`, `1+2->3`, `1+3->4`, `2+2->4`,
`2+3->1+4`, and `3+3->2+4`.

For each recorded source and target profile related by one move, the target
region matrix entrywise dominates the source matrix.  There are zero
violations among 240 edges in the complete profile lattices at occupations
2, 4, 8, and 16, and zero violations among 30 recorded edges in the
occupation-128 packed corner.  This is floating-point evidence at the
recorded Chernoff points, not an outward-rounded proof.  The consolidated
receipt is `receipts/region_packing_recorded_profiles_audit.json`.

The order is not parameter-free.  A grid over all 30 occupation-8 edges
finds counterexamples at low Bernoulli density and at sufficiently aggressive
Chernoff tilt.  In contrast, all 30 edges pass at `p=1/2` for every tested
surprisal from 0.001 through 1, and all 204 occupation-16 edges pass at
`p=1/2` for every tested surprisal from 0.01 through 1.  The optimized
occupation-128 proof point uses `p=1/2` and surprisal 0.0288.  Thus the useful
conjecture is a symmetric-envelope, bounded-tilt lemma; it must not be stated
for arbitrary `p` and `z`.  The grid receipts are
`receipts/region_packing_parameter_grid.json`,
`receipts/region_packing_p_half_grid.json`, and
`receipts/region_packing_p_half_occupation16_grid.json`.

Conditioning on the placement of every background packet is too strong.  At
occupation 128 and the optimized symmetric proof point, a diagnostic over 228
canonical or pseudorandom background/move pairs finds four violations.  The
four pure merge moves each fail on one highly structured canonical
background; their largest entrywise ratio is 1.00605.  The two overflow moves
have no sampled violations.  The largest merge discrepancy occurs in a tiny
live-to-live matrix entry.  Thus the desired proof must preserve the uniform
averaging of the background packets rather than induct after conditioning on
their positions.  See `receipts/conditional_region_packing_audit.json`.

Joint averaging does not appear to require all 32 epochs at once.  Complete
profile audits at occupations 8 and 16 find no tolerance-level packing-edge
violation after 1, 2, 4, 8, 16, or 32 averaged epochs.  One epoch is
analytically equality because every packet contributes to the same epoch
mass.  At depths 2 and 4 the closest live-to-live entries are at float64
resolution, so these runs are evidence for an epoch-count induction rather
than numerical certificates.  The receipts are
`receipts/region_packing_epoch_depth.json` and
`receipts/region_packing_epoch_depth_occupation16.json`.

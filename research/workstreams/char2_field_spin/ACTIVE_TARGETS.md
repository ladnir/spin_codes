# Acceptance criteria for field-SPIN

The existing binary SPIN implementation is the baseline, not an open proof
problem. Extending its binary generator to GF(2^m) preserves its minimum
symbol distance. A new half-rate construction proving only 10% distance
therefore does not establish a new rate--distance guarantee.

The public RS8 baseline certifies 10% distance at rate 1/2 and measures
approximately 0.100 / 0.461 / 2.862 ms at K16 / K18 / K20. These are the
recorded precomputed transposed-encoding timings on 128-bit elements, not
new matched measurements. See [the baseline index](../../CERTIFIED_BASELINES.md).

## Required improvement

The new family must offer a useful tradeoff beyond that baseline:

1. At 10% symbol distance, substantially higher rate and competitive total
   encoding cost. The first screening target is rate 3/4; rate 2/3 is an
   intermediate candidate, not an automatic endpoint.
2. At rate 1/2, substantially higher symbol distance at competitive cost.
   Screen 30% first, with the original near-50% ambition retained for large
   fields. A small-subfield implementation does not acquire a larger distance
   merely because its ambient field grows.

Both tracks require exactly one global routing stage and exclude general
large-field multiplications. The requested fields remain GF(2^m) for
m in {8,32,64,128}. A coefficient set of low-degree native-field polynomials
is allowed; it must not be confused with a fixed small subfield.

Forty bits of whole-code setup-failure margin is the initial proof target.
An increase from forty to eighty bits, with no useful rate, distance, or
performance benefit, is not a substitute for either improvement above.

## Work gate

The [completed high-distance optimization round](HD_OPTIMIZED_CHECKPOINT.md)
has a validated approximately 5.5-ms K20 forward encoder at rate 1/2,
40% GF65536-symbol distance, and more than 56.65 bits of whole-code margin.
The sub-5-ms aspiration remains unmet. Its concrete proof and implementation
alternatives have been evaluated and recorded; the next recommendation is
checkpoint review and promotion, not another parameter sweep. The older
search instructions below are preserved as history, not outstanding work.

The high-rate round is complete for now. The user accepted the
[frozen winner](GF8_HIGH_RATE_WINNER.md): 83.3333% nominal rate, 10% distance,
57-bit K20 margin, and about 3.2 ms precomputed transposed encoding time.
Keep its implementation and proof fixed; do not continue tuning it merely
because another small optimization might exist.

The active objective is high distance, ideally close to 50%. Start at rate
1/2 and requested K20, with a 40-bit whole-code setup-failure target. Retain
one regional routing stage and 128-bit application records. Reconsider the
coefficient alphabet and local maps jointly: the near-50% aspiration does
not follow by scalar extension of the frozen GF256 code. New scripts and
implementations must remain separate from that checkpoint.
The [restart note](HIGH_DISTANCE_RESTART.md) records the first 36% GF256
class screens and the next larger-alphabet comparison.
The [initial larger-field round](HIGH_DISTANCE_INITIAL_RESULTS.md) is now
complete: positive sampled classes at 42% over GF(2^16) and 45% over
GF(2^32), with validated faster multiplication kernels. These are not
whole-code certificates or encoder timings. The user then prioritized
reducing internal expansion before a closure campaign. The
[lower-expansion round](HIGH_DISTANCE_LOWER_EXPANSION.md) favors the
GRS[256,64]/keep2, GF(2^16) candidate: rate 1/2, 4K raw symbols, and
positive sampled bounds at 40%. Its first pruned-FFT implementation measures
about 15.6 ms at K20; the outer takes about 7.5 ms. First fuse the outer
stages to reduce scratch traffic, then test the route-to-outer layout
separately, before investing in directed certification.
The 8K version remains a proof control, not the optimization target.

The following paragraphs preserve the preceding high-rate search decisions;
their next-step recommendations do not override this pivot.

The subsequent [rate climb](GF8_RATE_FRONTIER.md) closes nominal 37/44
(84.0909%) at K20, 10% symbol distance, and 183-bit setup-failure margin.
Its requested rate is 83.8657% after padding. This is close to the GF256
random-code benchmark, not a proof of a maximum achievable rate.
The [serial native comparison](GF8_RATE_NATIVE.md) selects the 84% profile
at 4.807--4.849 ms, with 115-bit certified margin. Fixed-six takes
4.903--4.918 ms and 84.09% takes 5.318--5.350 ms. The latter has more raw
inner work despite slightly fewer retained outputs. Next compare the same
winning design at smaller sizes and replay their certificates separately.
Further proof-only rate climbing should first tighten the output-ball bound.

The explicit [three-quarter-rate goal](GF8_THREE_QUARTER_GOAL.md) now has a
complete [K20 proof](GF8_PUNCTURED_WHOLE_CERTIFICATE.md): nominal rate3/4,
10% distance, and116-bit whole-code setup-failure margin. Its [native
implementation](GF8_PUNCTURED_NATIVE.md) takes5.57--5.82ms at K20, roughly
the cost of the preceding lower-rate field control with13.38% fewer outputs.
This meets the first rate target, not binary speed parity or a maximum-rate
claim. Keep this historical proof/native checkpoint separate from production;
the later rate climb and matched comparison above supersede its next-step plan.

The [repaired cubic K20 candidate](GF8_PACKET_CERTIFICATE.md) now passes
the complete proof gate: requested rate 64/99, 10% distance, and 65 bits.
Its width8/state4 paired map has a substantially smaller arithmetic
schedule, but more shuffles per output and a 4.762% expansion increase
over the old [128,87] point. The [native gate](GF8_CUBIC_NATIVE.md) now
passes correctness. The first implementation regressed; a paired schedule
recovers 14--16%, reaching about .237 / .94 / 5.71--5.86 ms. This is still
slower than binary SPIN and does not improve the rate frontier. Do not
promote this candidate on arithmetic counts or extra failure margin alone.
The newer punctured checkpoint above reaches3/4. Enlarging the outer alone
did not repair the unpunctured cubic candidate's bound.

Before implementing another candidate, record its actual rate and distance
target, proof obstruction or bound, and estimated full encoding cost.
Count the outer, inner, coefficient access, routing, conversions, and output
stores. A positive singleton screen is insufficient: inspect intermediate
active-row classes before investing in a kernel.

Only a candidate with a credible advantage proceeds to implementation,
complete directed certification, independent audit, and serial measurements
at K16/K18/K20. Compare equal useful message sizes and disclose padding,
output lengths, representation conversions, and setup exclusions.
No two benchmarks run concurrently.

## Preserved but parked work

The GF16 half-rate Toeplitz[64,32] certificates and implementation are retained
as mathematical and implementation controls. They are not a successful
answer to the field-generalization objective. Their stronger failure margins
do not compensate for reproducing an existing distance guarantee more slowly.

Half-rate Toeplitz[48,24] proof work and the partial native refactoring are
also parked. The latter has not received final validation. Do not resume its
kernel or benchmark work solely to complete that side branch. The binary
production profiles and the proved two-permutation reference remain intact.

## Current decisions after the useful-target screens

The latest [outer-dimension sweep](GF8_OUTER_RATE_SWEEP.md) holds the
eight-state inner fixed. Dimension88 closes at K18 with 72 bits and 10%
distance, increasing requested rate to32/47. Its output reduction is only
2.083%, and no native timing has been measured. K20 dimension88 remains
unclosed; dimension89 has sizable middle-bound gaps. Preserve this small
proved improvement, but seek a stronger design or proof change before a
larger rate/performance campaign. The implementation timings below remain
for dimension87.

The active [GF256 exploration goal](GF8_EXPLORATION.md) follows the
medium-alphabet comparison. Rank-history and source refinements target
a cheaper inner at10% distance. A parallel coefficient-basis outer removes
systematic interpolation while preserving the certified full code image.
Both native schedules pass correctness and sanitizer checks. Fused source
loads save about 15% at K16/K18; both schedules regress at K20, where the
systematic reference remains preferred.
The cheaper six-state inner remains unclosed. Neither result establishes
a speed advantage over binary SPIN or a new rate--distance certificate.

The user authorized [medium-alphabet exploration](MEDIUM_ALPHABETS.md)
for the 10% track. Keep 128-bit application records, comparing sixteen
GF256 lanes with four GF32 lanes. The first native arithmetic controls
and paired proof screens favor GF256 as the baseline. Full GF32 diagonals
do not rescue the tested smaller states or 3/4-rate construction under
the current cover. The next gate is a middle-class or geometry improvement,
not more singleton margin. Existing certificates below remain preserved.

The smaller GRS[128,87] candidate has complete 10% certificates at
K16/K18/K20 with 63/79/67 bits. See the [K16/K18 theorem and replays](FIELD128_SIZE_CERTIFICATES.md)
and preserved [K20 certificate](FIELD128_K20_CERTIFICATE.md). Requested
rates are 2/3, 2/3, and 128/189. Its explicit four-row-interleaved source
order preserves the full code and avoids the final record transpose.
[Latest matched results](FIELD128_NATIVE.md) are about 0.288/1.128 ms at
K16/K18, with K20 varying from 5.73 to 7.46 ms. Binary SPIN remains
0.099/0.46/3.17 ms. This passes the higher-rate proof gate, not the performance gate.

Keep the dense inner: both exact-map FFT alternatives were validated but
regressed, including the new lane-local layout. Compact-route indices give
modest direct gains, not a resolution of the full performance gap. The
previous round is complete; the medium-alphabet comparison above is the
subsequently authorized direction.
Do not spend the available margin solely on increasing confidence.
Reduced-state [128,86] has a middle-class gap, and several cheap postmixes
have actual obstructions. [The design status](CHEAPER_DESIGN_STATUS.md)
separates those cases. The higher-distance track remains open. Prior proofs
and production stay intact.

### Larger-outer reference

The K20 GRS[256,176], width-sixteen/state-eight construction now has a
[whole-code 10% certificate](FIELD256_K20_CERTIFICATE.md) with 119 bits
of setup-failure margin. Its padded rate is 11/16; its requested-input
rate is 32/47. This passes the useful rate/proof gate, not the performance
gate. The [first native encoder](FIELD256_NATIVE.md) passes independent
correctness checks, but its measured cost does not pass the performance gate.

The exact-map FFT, packed-inner and fused-store pass is complete. Its
matched times remain 2.3--2.8 times slower than binary SPIN. The next task
is to seek cheaper compatible maps while preserving the higher-rate
objective and controlling regional window obstructions. Measure full encoding, including routing,
packing and stores; do not infer a speedup from the shorter output or
field-product count. Preserve binary production. K16/K18 and the cheaper
width-thirty-two candidate require separate certificates.

### Preserved earlier decisions

The GRS[32,22], four-state profile is now **parked after two verified
realized-route obstructions**. [Unequal windows with different row
supports](WINDOW_DIMENSION_SEARCH.md) force output weight below 10% at
K20, irrespective of the inner coefficients. Both witnesses survive
removing the row-major padding. Their probability under ideal regional
permutations is not proved here; the sampled evidence nevertheless
changes the engineering decision. Do not continue coefficient-only
refinements as a repair for these routes.

The next gate at that checkpoint compared changed geometry at a useful rate. GRS[64,44]
retains rate 11/16 and the four-state circuits; GRS[33,22] retains the
small outer with rate 2/3. Increased state is another option, but its
circuit cost must be established. Screen arbitrary row supports and
unequal windows before a new certification or native implementation
campaign. Absence of a found window witness is not a distance proof.

The remainder records preceding decisions and preserved class results;
it does not override the current certified candidate or performance gate.

The rectangular Toeplitz and GRS[16,12] screens did not close at higher
rate. Exact joint local moments remove substantial numerical slack, but
not the failing intermediate classes. See [RATE_TARGETS.md](RATE_TARGETS.md).
The positive integer local enumerator is retained as a proof tool, not a
new rate--distance result.

A full uniform global permutation followed by a causal inner is now ruled
out for the proposed small GRS outers. A suffix retains exactly the outer
subspace shortened into it, independently of the inner coefficients or
width. See [the finite whole-code obstruction](HIGH_DISTANCE_GLOBAL.md).
One global routing stage does not mean that its distribution must be the
full uniform distribution. Keep regional coordinate spread.

GRS[32,24] with state mixing only at region boundaries passes the one-row
gate but still fails forty bits if injection labels repeat every sixteen
or thirty-two positions. Two shortened outer rows on seventeen matching
coordinates have a cancelling subspace for every diagonal and regional
matrix assignment. This is another actual construction obstruction.

Expanding the label pattern to sixty-four positions fixes that tested
pair mechanism, but does not make the proposed code viable. A complete
37-row window event gives failure probability above 2^-38 at K20, for
every coefficient assignment and every injection period. The outer
shortened space has dimension 37, while nine zero-state constraints can
remove at most 36 dimensions. This is a construction failure, not an
unresolved upper bound. See [RATE_DESIGN_REVIEW.md](RATE_DESIGN_REVIEW.md).

The new comparison favors more outer parity over a larger state. Equal
window lengths gave false encouragement: the broader unequal-length
window upper bounds have only 27.17 and 19.84 bits at K20 for [32,23] and
[33,24]. Those are insufficient upper bounds, not failure lower bounds.

The preceding proof candidates were GRS[32,22] and GRS[34,24], both with four
state symbols. They have rates 11/16 and 12/17. At K20 their tested
unequal-window event families have upper bounds below 2^-90 and 2^-74;
their one-row contributions have more than 70 and 46 bits, respectively.
See [the capacity comparison](REGIONAL_CAPACITY_COMPARISON.md). These
families do not cover all sparse multirow words.

A new [dense-class theorem](REGIONAL_GLOBAL_ROUTE.md) covers every bad
output whose inverse-inner input has weight at least .24N for [32,22],
or .25N for [34,24], with exact rational finite bounds. It applies to any
invertible inner independent of the route and coordinate diagonals.
The remaining task is a complete sparse multirow covering argument,
clipping conditional counts before averaging rare routes. No new native
kernel is justified by the currently separated classes alone.

An earlier refinement added the [dual-MDS emission repair](DUAL_MDS_PHASE.md)
to the [32,22] candidate. Reversing four state coordinates before emission
adds no field products and makes both period-64 maps MDS. This is a change
to the construction, not an exact-map implementation optimization. Its
full certificate remains open; the earlier construction is preserved.

The [sparse proof checkpoint](REGIONAL_SPARSE_COVER.md) now gives an
output-first packet factorization, explicit state-rank charges, and exact
zero/one-output rank transitions without reusing random diagonals.
The remaining task is to compose these local statements with outer-row
freshness into a finite operator covering every sparse class.

A separate [fresh-Toeplitz packet study](RANDOM_PACKET_GATE.md) has an exact
local law but does not close the tested rate-2/3 middle classes. Its rate-5/8
backup has positive sampled bounds, not a certificate. The
[source-envelope and cost check](NONEMPTY_SOURCE_ENVELOPE.md) does not
justify further implementation or a certificate campaign for that backup.

The 30% and near-50% distance targets remain open. Rejecting the uniform
route does not complete that track or justify another half-rate/10%
implementation campaign.

## Preserved construction and class proofs

The parked higher-rate candidate uses [Cauchy-form emission and a
fresh four-coordinate state diagonal](DIAGONAL_REFRESH_PROOF.md) in every
packet. It keeps GRS[32,22], the single regional route, and the 10% target.
A local rank-moment argument gives nearly the root-equation contraction
of a full GL4 refresh, without claiming the same subspace distribution.

Its direct fixed circuits need fewer field products: the complete count
is 169/88 rather than 181/88 GFNI512 equivalents per source record. Native
packing and runtime remain unmeasured. Directed replay now covers all
inverse input weights at least .11N at K16/K18/K20; nonzero inputs below
.11N with at least three active outer rows remain unresolved. Separate
one-row and [two-row](CAUCHY_Q2_TRANSFER.md) replays give 107/105/103 and
65/72/77 bits, respectively. Existing one-row/window results for the old
constant-state construction must not be silently transferred to this
fresh-diagonal construction.

The subsequent proof operator retains actual state-coordinate support. The new
joint-return bound improves the sparse transfer, but the two-state
zero/nonzero relaxation forgets the cost of reaching rare sparse states.
Exact regional counting and the positive zero-set rank expansion are
implemented. The [17-state composition](MASK_REGIONAL_SCREEN.md) was
tested and remains vacuous in larger active-row classes. Exact total-weight
truncation also gives negligible improvement there. Larger outers at the
same rate and width-eight inner packets were screened, not promoted.

Acceptance still requires the entire sparse class and matched native
performance, now for a geometry that addresses the new obstruction.
Neither a class certificate nor a multiplication count satisfies the
goal by itself. Preserve the original regional candidates and their
proof receipts.

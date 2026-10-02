# Packed GL32: Toward 10% Distance and 40 Bits

Implementation follow-on: [wide BCH stores and compact GL32 coefficients](../packed_bch_tune_REPORT.md)
now give about **5.85--5.90 ms** without changing the closed R4 construction
or its >10% / >52.05-bit certificate. The 6.212096 ms timings below record
the earlier proof-closing checkpoint.

2026-09-30. The retained encoder in [the first closure](FIRST_CLOSURE.md)
is unchanged: K=2^20, rate one-half, canonical GL32 mixing, four-bit routing
packets, and IMT(128,19) with two updates. Its retained precomputed transpose
time is 6.203224 ms. A separate R4 candidate and its new matched timings
are reported below. No production defaults have changed.

The retained R2 whole-code result is **>9.5% distance / >36.68-bit margin**.
Its entire sparse prefix q=1..32 closes at 10% with **47.58879650 bits**,
but its dense bound remains open. The separate R4 variant has now passed
fresh 384-bit replay of all 32 sparse occupancies and all 337 dense intervals.
The completed whole-code replay establishes **>10% distance with >52.05 bits
of setup-failure margin**. Its exact-dyadic aggregate has displayed margin
**52.0528843554777 bits** at cutoff 209715. The coordinator completed its
exact aggregation and source-provenance checks successfully.
The [fused R4 implementation](../packed_fused_r4_REPORT.md) measures
**6.212096 ms**, preserving exactly the same sampled R4 encoder.

## An Exact Outer-Code Improvement

Let H be the number of nonempty canonical eight-column blocks in a tuple
of four BCH rows. The old cumulative bound at H=5 allowed 3,020,640
nonzero tuples. At distance 9.6%, this term dominated the one-group bound
by approximately 25 bits over H=6.

For the actual production BCH generator, **no nonzero codeword is supported
on five or fewer canonical octets**. The new checker establishes this
without estimating the unknown BCH spectrum:

1. Read all 128 generator rows in the production coordinate order. Check
   binary rank 128, even parity, and divisibility of each punctured row by
   the primitive BCH polynomial with 36 consecutive roots. These checks
   establish minimum distance at least 38.
2. Build a linear quotient map whose kernel is exactly the generator span.
3. A word contained in five octets uses at most 40 coordinates. Even parity
   and distance 38 force weight 38 or 40. Thus it is the full 40-position
   union, possibly with exactly two positions deleted.
4. Enumerate all 201,376 five-octet unions and compare their syndromes with
   zero and the 32,640 pair-column syndromes. No codeword is found.

Four or fewer octets contain fewer than 38 coordinates. Consequently every
nonzero row, and every nonzero four-row tuple, occupies at least six octets.
The GL32 map is invertible within each octet block, so its output has at
least six active packets. The previous proof safely used the weaker floor
five; it is not invalidated by this improvement.

This statement is specific to the **canonical octets**. It does not show
that shortening on every arbitrary 40-coordinate set has dimension zero.
The uniform shortening-dimension table is deliberately left unchanged.

Implementation: `outer_hill_octets.py`. An independent review checked the
coordinate order, BCH-distance premise, quotient map, and enumeration.
Tests include exhaustive toy codes and positive examples, not just the
production zero result.

## Counting Larger Supports by Incidence

Let F_r(h) count nonzero ordered four-row tuples of row rank r whose union
occupies at most h canonical blocks. Let d(u) bound the shortened dimension
on every u-coordinate subset, and let R_r(d) count rank-r ordered four-row
tuples in a binary vector space of dimension d.

For each t>=h, count pairs consisting of such a tuple and a containing
t-block set. A tuple with s<=h active blocks belongs to
binom(32-s,t-s)>=binom(32-h,t-h) such sets. Therefore

    F_r(h) <= floor(binom(32,t) R_r(d(8t)) / binom(32-h,t-h)).

The checker minimizes over t, intersects the previous rank caps and the
exact H<=5 result, and enforces monotonicity. Division and flooring occur
on integer tuple counts, before the rational GL32 expectation transport.
For example, the H=6 rank-one codeword cap improves from 906,192 to 97,092
before the factor 15 for nonzero ordered four-row labels.

Implementation: `outer_hill_incidence.py`. Its fresh authenticated interface
returns both the transported CDF and a new premise schema. It does not
change the old replay interface or old receipts. The direct-shell comparison
is also rebuilt from the improved canonical counts; CDF differences are
never treated as shell bounds.

## R2 Sparse Progress

The first fresh 256-bit run covers occupancy q=1 at all three cutoffs below.
Here q counts nonzero four-row groups among the 2,048 groups. The run uses
the exact H5 exclusion, but not the additional incidence improvement.

| Relative-distance threshold | Integer bad-weight cutoff | q=1 margin (bits) |
|---|---:|---:|
| 9.5% | 199229 | 75.18694388 |
| 9.6% | 201326 | 74.56770350 |
| 10% | 209715 | 72.10727366 |

Thirteen output tilts were checked with outward arithmetic. The original
three tilts already gave 58.28946 bits at 10%; finer tilts improved this to
72.10727. Averaging before or after support transport alone gives no gain
for the old monotone support-event bounds. The improvement comes from
excluding H=5 and optimizing the event bound.

The local receipt is `tmp/packed-hill/q1-h5exact-tilts-d095-d096-d10-p256.json`;
SHA256 `b630665eaabb03a091086c4b501fba26c871e1334c6a1a0097e545288068f9da`.
Its expected-CDF hash is
`8ab2d6560ef6e59b685c12713626d7b71c1675306effaca2dd4e420e7f9bd1cc`.

A separate fresh run uses H5 exclusion plus incidence for q=2..32 at 10%.
Its requested per-occupancy target is 46 bits, leaving room to sum the
contributions. **All 31 requested occupancies pass.** The weakest is q=26
at 48.5046 bits; q=32 closes at 57.5021 bits after 41 support-box splits.
The exact sum of these endpoints and the q=1 endpoint has margin
47.588796500773 bits and is strictly below 2^-40. The two count refinements
differ but bound the same encoder; they may therefore be combined across
their disjoint message classes.

The q=2..32 receipt is
`tmp/packed-hill/sparse-q2-32-h5-incidence-d10-p256.json`, SHA256
`298f963586eb39916d725b051f54a23c06e35d30d650b8c3754e58ffb33cd40c`.
The exact aggregate is `tmp/packed-hill/sparse-q1-32-d10-aggregate.json`,
SHA256 `4b2b6aed6f7ec630087e12a04e9b7dc11d59d48d672f53cb59e2327cbbcb952f`.
The scope check verifies all 32 occupancies without overlap, cutoff 209715,
R2, support floor six, and the common production-generator SHA256
`c13131ca9d02861ac5b4688849a4550196120abdce5c5aa0bec7df71eac02424`.
This is still **sparse-only**. Replaying the sparse witnesses at 384-bit
precision remains part of assembling a new final whole-code certificate.
An independent audit rechecked the source and generator hashes, recomputed
the exact rational count transports, validated all 31 support-partition
trees, and reproduced the endpoint sum. Its upper endpoint is approximately
4.72437971065e-15. This audit did not replace the planned 384-bit replay.

## R2 Dense Progress and the Remaining Obstruction

A fresh direct-shell comparison at 10% still fails at several individual
comparison means. Positive log2 upper bounds below are **vacuous bounds**,
not probabilities or evidence that the code has low distance.

| Comparison mean | Initial fresh log2 first-moment upper bound |
|---|---:|
| .032 | +892.286 |
| .104 | +3027.532 |
| .114 | +4257.291 |
| .120 | +4366.486 |
| .208 | +3992.722 |

These are q=33..2048 point bounds, not exhaustive mean coverage. Their
receipt is `tmp/packed-hill/r2-d10-direct-points-p256.json`. Unlike the old
failed 10% probe, this run uses the tighter direct expected-shell majorant.

An isolated `hill_probe.py` retargets only the validated construction
metadata, freshly authenticates counts, and drops old accepted bounds and
point hints. It forces regional refinements, retains the ordinary and
variance-only alternatives, and searches both classified and uniform-source
variants over a wider output-tilt grid. Each selected result is independently
evaluated outward. Neither a floating search score nor a successful point
is an exhaustive certificate.

The broad grid from 1/8 through 5/4 of the IID starting tilt did not resolve
the first two difficult means. This coarse grid does not include every
smaller refinement in the original driver's grid, so its minimum must not
replace a better existing bound. The next probe combines the new outer
counts with 64 rather than 16 variance intervals. At mean .032 its first
checked bound is still +868.203 bits, only a modest improvement over the
initial +892.286. All five selected means still fail. The narrower tilt
grid used with 64 bins differs from the baseline grid, so worsening at a
point is not evidence that improved counts or subdivision weakened the
underlying best bound. None of these tested candidates closes 10%.

A final nearby-tilt test checks factors .85, .90, and .95 against both
regional variants. At mean .032, the classified variant at .90 improves
the fresh outward result to **+669.283890 bits**, still vacuous. Its receipt
is `tmp/packed-hill/d10-refined-counts-v64-near-tilts-p256.json`. This removes
the coarse-grid ambiguity at that point but does not locate a global optimum
or prove that the construction fails at 10%.

`hill_variance_diagnostics.py` exposes each variance interval's floating
contribution, separating outer count and inner terms. Its interval-deletion
counterfactuals deliberately omit messages and are never certificate
evidence. This attribution is intended to select the next refinement from
the actual dominating term rather than an indiscriminate parameter sweep.

At mean .032, the diagnostic exactly reproduces the regional proposal
+868.2030329874578. Its dominant interval is variance [.0095,.0100], with
outer variance multiplier gamma=0. The adjacent intervals have only about
14.34 and 2.20 bits of saved-dual endpoint sensitivity. Thus large outer
affine endpoint slack is not the evident obstruction. This does not rule
out improvements to the count-MGF or inner bounds within those intervals.

For that dominant interval, the pointwise count-probability caps sum to
C=1.0650208383. Imposing only normalization on those caps cannot repair the
gap: the caps divided by C are a feasible distribution for that relaxation,
so the optimal nonnegative weighted-sum bound improves by at most C per
region, or about 23.27 bits over 256 regions. Dividing the caps by C is
**not** itself a valid upper bound on the actual distribution. This is a
diagnostic of one proposed relaxation, not a certificate refinement.

Additional [outer refinements](OUTER_REFINEMENTS.md) are recorded separately.
An intersecting-family argument targets canonical supports near H=19/20;
it was not enabled in the preceding runs. It is now an optional, separately
authenticated refinement. Extra H6 packing bounds are lower priority
because the sparse prefix already exceeds the requested margin.

## R2 Joint Count Bound: Implemented, but Insufficient Alone

The regional count is one normalized random variable J, conditional on
the outer comparison-component choices. Bounding its atoms independently
can overestimate a weighted regional sum. `joint_count_dual.py` instead
optimizes that sum subject to normalization, atom caps, raw MGF bounds,
and optional first and second moments. A floating LP proposes a dual;
fresh outward arithmetic checks every integer count and charges any
rounding repair to its atom cap. Saved duals can be verified without the
solver. The normalized law is the law of J, before the factor tau^-J.

`joint_count_region.py` applies the bound separately to every matrix entry
and retains all variance intervals. It takes the entrywise minimum with
the old bound, then composes the resulting nonnegative matrix over all
256 regions. An independent review checked conditional normalization,
the raw MGF interface, and this composition. Tests also compare with
direct repeated-region calculations for small noncommuting matrices.

At mean .032, with 64 variance intervals and the retained .90 output-tilt
factor, the first pilot refines only intervals 18, 19, and 20:

| Variance interval index | Old log2 contribution | Joint-count log2 contribution |
|---|---:|---:|
| 18 | +658.0214 | +629.2363 |
| 19 | +667.8156 | +637.0473 |
| 20 | +668.6338 | +637.4091 |

All other intervals remain in the sum with their old valid bounds.
The full point bound improves from +669.283890 to **+659.232459**;
the unrefined interval 21 then dominates. This is a valid improvement,
but remains vacuous. The 29--31-bit improvements in the refined terms do
not justify expecting this relaxation alone to close a roughly 669-bit gap.
The receipt is `tmp/packed-hill/joint-count-d10-m032-parts18-20-p256.json`.

## Construction Candidates: Additional Inner Updates

The optional `hill_probe.py --updates 3` and `--updates 4` modes change the
encoder, not just its proof. They retain the same BCH outer, GL32 maps,
packet route, and fixed IMT expansion/feedback maps. Each step applies
three or four independent sampled transvections instead of two. The
local proof rebuilds the actual update law: its lazy mass is 2^-r for r
updates. It does not assign independent fresh bits to the whole state.

These candidates have distinct ensemble metadata. Neither the retained
R2 sparse certificate nor its 6.203224 ms timing applies to them. No
production default has changed. At 10%, the fresh R3/EKR probe clears
mean .032 with log2 upper bound -126.393982, but still gives +505.029740
at mean .104 and +1015.722207 at .114. Extra mixing materially improves
these bounds but three updates alone have not closed all tested points.
R4 clears the four selected higher-density points:

| Comparison mean | R3 log2 upper bound | R4 log2 upper bound |
|---|---:|---:|
| .104 | +505.030 | -1113.933 |
| .114 | +1015.722 | -1333.692 |
| .120 | +995.281 | -1442.383 |
| .208 | -348.485 | -3475.958 |

Both columns use fresh EKR counts and adaptive point searches. These are
not bounds on the intervening intervals. The R4 receipt is
`tmp/packed-hill/r4-ekr-d10-adaptive-p256.json`. Its favorable results
motivate an exhaustive dense search and a new q=1..32 search at 10%.
The dense per-cell target is 52 bits; the sparse per-occupancy target is
48 bits. Only their fresh aggregate can establish a whole-code margin.

The R4 sparse search completed every q=1..32 at 10%, and fresh 384-bit
replay now confirms an aggregate margin of **58.0241331220138 bits**.
The selected q=6 bound, rather than q=1, limits this sparse result.
The final whole-code sum must also include all dense occupancies.

The earlier [matched implementation comparison](../packed_updates_REPORT.md) measured
**6.784301 ms for R4 versus 6.170696 ms for R2** in one executable: +9.94%.
All R2/R3/R4 reference and adjoint checks pass at K=2^14 and K=2^20 for two
seeds. The extra measured cost is in the inner/route phase; the outer and
route distribution remain unchanged. All benchmarks ran serially. This
does not alter the retained earlier R2 measurement or certificate.

The subsequent [exact-map fusion](../packed_fused_r4_REPORT.md) reduces R4
to **6.212096 ms**, versus **6.753111 ms** for sequential R4 in that same
executable, an 8.01% reduction. Setup precomputes the product of the original
four transvections; it adds no randomness and changes no sampled map.
The measurements cover full precomputed transposed encoding at K=2^20
with 128-bit XOR elements, excluding setup and allocation. Reference,
adjoint, full-buffer, and preserved-suffix checks pass for both tested sizes
and seeds. These timings do not replace the earlier R2/R4 matched campaign.

`hill_cover.py` and `sparse_hill_variant.py` keep update counts and count
premises explicit. The dense helper starts coverage from the full comparison
domain and ignores point endpoints. Sparse operators are regenerated with
the requested update count. A new replay must cover all dense cells and
all 32 sparse occupancies before the final exact sum.

The first dense search exposed expensive repeated optimization in narrow
low-density cells. The search now partitions a frozen frontier among
independent workers and can reuse same-scope rational witness parameters.
Every accepted bound is evaluated afresh; the regional operator cache is
local to one search process. A checked merge retains every unresolved
subtree and discards saved numerical bounds. The final replay is unwrapped
and does not use these search caches.

Both higher-mean shards closed all 190 assigned intervals. A frozen
checkpoint then retained 248 accepted leaves and 89 unresolved frontier
cells. Four continuations closed all 89, giving **337 accepted leaves and
zero unresolved cells**. An independent geometry audit verifies an exact
partition of [2079/7829504,1], with depths six through eleven. This is
coverage of an auxiliary comparison interval, not a fraction of all
messages or failure probability. Reused witnesses supply only parameter
proposals; no interval inherits an old numerical bound.

The complete search is `tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json`,
SHA256 `a44142155d23553b589983cfb87d3ed3f8a0514217e2e46a0efb4477065bf01c`.
The first 384-bit replay independently confirmed all sparse occupancies
with 58.02413312 bits. Its serial dense phase was stopped to add parallel
verification; the partial directory `tmp/packed-hill/r4-d10-whole-p384/`
contains no whole-code receipt and must not be presented as a certificate.

The renewed run is `tmp/packed-hill/r4-d10-whole-parallel-p384/`.
Its sparse stage and four independent, uncached dense workers have finished.
Each worker freshly authenticated the count premises and rebuilt the actual
R4 inner. Every cell is accounted for, with no unresolved interval.
An independent audit summed their positive dyadic endpoints by exact integer
alignment, avoiding expensive rational-denominator reductions:

| Fresh replay scope | Audited margin (bits) |
|---|---:|
| Sparse q=1..32 | 58.0241331220138 |
| Dense q=33..2048, all 337 intervals | 52.0760654522247 |
| Sum of both scopes | 52.0528843554777 |

These margins are logarithmic displays of exact endpoint sums; exact
arithmetic decides the probability threshold. The coordinator completed its
separate exact sum and source-freeze checks and wrote a successful
`whole.json`, SHA256
`637a1ac049a7aa81441705e5a4fbeebf79e78fa80e63fe54f93cfe6dad05a8f0`.
Its outward endpoint U also satisfies the exact integer check
U^20 < 2^-1041, proving margin greater than 52.05 bits.
The resulting minimum distance is at least 209716 out of N=2097152,
strictly greater than 10%, except with probability at most U over independent
ideal setup. Setup is sampled once and fixed for every message. This claim
does not certify the implementation's pseudorandom setup generator.
See [the R4 closure record](R4_CLOSURE.md) for the full claim and provenance.

## Reproduction and Next Decision

Run from the repository root; use fresh output paths. Numerical receipts
belong under ignored `tmp/`, not in version control.

The complete R4 search witnesses support the following fresh whole replay:

```text
python -B research/workstreams/permutation_locality/packed_mixing/whole_hill_replay.py --dense tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json --sparse tmp/packed-hill/sparse-r4-ekr-d10-p256.json --output-dir tmp/packed-hill/new-r4-whole-p384 --precision 384 --target-bits 40 --dense-workers 4
```

The following commands reproduce the earlier R2 outer-count and sparse work:

```text
python -B research/workstreams/permutation_locality/packed_mixing/outer_hill_octets.py
python -B research/workstreams/permutation_locality/packed_mixing/sparse_hill_q1.py --h5-exact --output tmp/new-q1.json
python -B research/workstreams/permutation_locality/packed_mixing/sparse_hill_prefix.py --incidence --output tmp/new-prefix.json
```

For dense probes, pass the completed 9.5% search snapshot as the source to
`hill_probe.py`, with `--distance .1`, selected `--means`, and a new output.
Use `--outer-refinement --variance-bins 64` to rebuild the comparison from
the new counts. The new context schema is intentionally rejected by the
legacy dense replay, whose count premises differ.

Keep the retained R2 encoder and its old whole-code certificate fixed when
revisiting the R2 dense obstacle. The implemented jointly constrained regional count
bound provides a reusable proof refinement. The original method bounds each atom Pr[J=j] and
then sums the weighted bounds. A linear-programming dual can instead impose
normalization, atom caps, and several moment-generating-function (MGF)
inequalities together. For a nonnegative regional objective f(j), check

    f(j) <= alpha + sum_t beta_t exp(t*j) + sum_k gamma_k 1[j=k]

at every integer j from 0 through 2048, with beta_t,gamma_k nonnegative
and alpha unrestricted. This yields the upper bound
alpha + sum_t beta_t M_t + sum_k gamma_k c_k, where M_t bounds
E[exp(tJ)] and c_k bounds Pr[J=k]. Floating optimization may propose the
dual, but both its pointwise inequalities and final value need outward
verification over the full mean/variance scope.

Use **raw** bounds from `mgf_upper`, not values multiplied by the tilted-atom
correction: that correction bounds individual probabilities, not the total
MGF. A matrix or continuation-vector bound must remain valid for every
entering inner state; a fitted eigenvector alone does not prove a
composable regional inequality. Merely normalizing the atom caps is
insufficient, as quantified above.

The R4 milestone is complete: **>10% distance / >52.05-bit margin /
6.212096 ms** for precomputed transposed encoding. Preserve the whole receipt,
component endpoints, and the retained R2 archive before promotion or a new
proof target. The larger packed GL32/BCH phase is the next implementation
target, using the fused R4 baseline without changing the sampled distribution.

Validation with scoped witness reuse, shard merging, and fresh parallel
replay: **260 packed-mixing tests pass**. The local probability kernels, old certificate receipts,
retained R2 construction, and production defaults remain unchanged. The
research driver additionally instantiates R3/R4 for matched measurements.

# Field-SPIN checkpoints

The user rejected further half-rate/10% optimization: that rate--distance
guarantee already follows from binary SPIN. The older half-rate/10% results
below are retained controls, not evidence of a useful field-code advance.
The high-rate and high-distance checkpoints address the useful targets.
[ACTIVE_TARGETS.md](ACTIVE_TARGETS.md) defines the new-work gate
against the existing implementation and certificate.

## Current useful-target status

The user accepted the [3.2 ms high-rate winner](GF8_HIGH_RATE_WINNER.md) on
2026-10-05. A verified source-only snapshot preserves the implementation
and proof dependencies. The subsequent
[high-distance checkpoint](HD_OPTIMIZED_CHECKPOINT.md) now reaches 40%
GF65536-symbol distance at rate 1/2, with more than 56.65 bits of whole-code
margin and about 5.5 ms full precomputed forward encoding at K20.
Its certificate passes independent 384-bit replay. The near-50% ambition
and sub-5-ms aspiration remain unmet. Production is unchanged.

The following entries preserve earlier search stages; their open tasks
are not the current status of the optimized checkpoint.

The [initial high-distance exploration](HIGH_DISTANCE_INITIAL_RESULTS.md)
has positive sampled sparse, middle, and dense bounds at 42% over GF(2^16)
and 45% over GF(2^32), using GRS[64,8], width 4/state 4, keep 1. Coverage and
directed certification remain open. New GF(2^16) multiplication is close
to byte cost in the measured pass; precomputed GF(2^32) matrices improve
the old tower kernel by about 2.7x at K16/K18 and 25% in time at K20.
Those are arithmetic-only measurements, not complete encoder timings.
The 49% larger-field middle bound remains negative at tested points.

The [lower-expansion round](HIGH_DISTANCE_LOWER_EXPANSION.md) reduces the
promising 40% geometry to 4K raw symbols using GRS[256,64] and GF(2^16).
Its point screens remain uncertified. Full length-64 cost controls confirm
that the old 8K design is expensive, around 29/37 ms for 16/32-bit fields
in the first prototype launch. A separately proved window obstruction
rules out 40%-distance/40-bit margin for direct-2K GRS[64,32]/state4 under
the specified routing distribution, independently of the diagonal field.
The actual length-256 candidate now has a pruned-FFT implementation and
measures 15.62–15.63 ms in place, with similar out-of-place time. The outer
takes about 7.5 ms. It is a proof/cost checkpoint, not a promoted winner;
next investigate the outer and its memory layout.

The [speed-first round](GF8_SPEED_NATIVE.md) supersedes the old84% speed
recommendation with a fully certified [64,40], width4/state4 profile.
It gives 83.3333% nominal rate, 10% distance, and 57-bit K20 margin.
The [exact-map follow-up](GF8_SPEED_ROUTED.md) reaches3.184--3.201ms with
relocated diagonals, sixteen-packet scheduling, a fused outer, and padded
scratch columns. Restoring input before every call gives3.332--3.361ms,
about20% less time than the matched earlier schedule. Production is unchanged.

The [phase audit](GF8_SPEED_PROFILE.md) subsequently tested sixteen
exact-map controls. The inner and route each consume over one third of the
remaining time; no tested schedule materially beats the selected3.2ms path.

The [rate frontier](GF8_RATE_FRONTIER.md) reaches 84.0909% nominal rate at
10% distance and K20. The [native comparison](GF8_RATE_NATIVE.md) recommends
the 84% profile instead: 4.807--4.849 ms for precomputed transposed encoding
of 128-bit records, with 115-bit setup-failure margin. It reduces retained
output length by 10.714% and encoding time by about 9% versus the matched
75% frozen control. The 84.09% profile takes 5.318--5.350 ms for only
another 0.3265% output reduction. Both remain isolated research checkpoints.

The [cubic8 native round](GF8_CUBIC_NATIVE.md) is complete: correctness and
sanitizers pass. Pairing consecutive packets removes a register dependency
in the first slow implementation, reducing time by 14--16% to about
.237 / .94 / 5.71--5.86 ms at K16/K18/K20. The code and proof are unchanged
by this schedule. Its lower rate and slower encoding than binary SPIN
still prevent automatic production promotion.
The [higher-rate outer screen](GF8_HIGH_RATE_OUTER.md) also rejects outer
enlargement alone as a remedy for the present 3/4-rate proof gaps. These
are bound failures, not counterexamples to the large-outer codes.

The [repaired width8/state4 inner](GF8_PACKET_CERTIFICATE.md) has a complete
K20 certificate at requested rate 64/99, 10% symbol distance, and 65-bit
setup-failure margin. GRS[128,83] expands 4.762% more than the old [128,87]
profile; its packed inner uses 18 instead of 70 GFNI per sixteen
outputs. This is a proof and arithmetic-count result, not a native speedup.
The exact algebra, complete sparse/middle/dense joins, and 46 tests pass.
The native follow-up above preserves it as a research control.

The [outer-rate sweep](GF8_OUTER_RATE_SWEEP.md) adds a complete K18
GRS[128,88] certificate: requested rate32/47, 10% distance, and 72-bit
whole-code margin. N falls from393216 to385024 (2.083%). No encoder time
for this profile is claimed. The same dimension at K20 remains open;
the previously closed dimension87 checkpoints and timings below are intact.

The [GF256 exploration](GF8_EXPLORATION.md) preserves the certified
GRS[128,87] code image while replacing systematic source coordinates with
additive-polynomial coefficients. Its final fused-load path takes about
.244/.956 ms at K16/K18, about 15% faster within this higher-rate family.
Both new schedules regress at K20, where the systematic reference stays
preferred. Native correctness, sanitizer and large-size checks pass.
The new rank-history/source proof refinements do not close the proposed
six-state inner. These results do not change production or certificates.

The [medium-alphabet comparison](MEDIUM_ALPHABETS.md) now has validated
GF256 and two genuine GF32 arithmetic controls on 128-bit records.
At K18 their isolated coefficient passes take about85us versus338--341us.
These are not full encoders. The paired GF32 proof screens barely improve
the difficult middle classes and close no new rate/state configuration.
Seven proof tests and native Release/ASAN/UBSAN checks pass. Existing
whole-code certificates and production remain unchanged.

The [GRS128/87 K20 checkpoint](FIELD128_K20_CERTIFICATE.md) is now complete:
10% distance, 67-bit whole-code margin, effective rate 128/189. All active
row counts 1..1105 and all inverse weights 46449..1548288 are covered, with
the MDS minimum 42 closing the overlap. The component margins 67/408/13690
are summed as exact probability bounds, not combined by their minimum.

The naturally padded K16/K18 sizes now have separate complete certificates
with 63/79 bits at requested rate 2/3. The [size replay](FIELD128_SIZE_CERTIFICATES.md)
uses a sharper exact modal-conditioning factor and rechecks the dense
exponent. Independent full replays and all ten new tests pass.

The [latest native round](FIELD128_NATIVE.md) measures about 0.288/1.128 ms
at K16/K18. Compact regional indices and batched diagonals improve direct
K20 comparisons modestly, but final-build matched launches range from
5.73 to 7.46 ms. The slow launches are retained and unexplained; this is
not a dependable speed improvement over the earlier 6.01--6.10 ms checkpoint.
Matched binary SPIN takes 0.099/0.46/3.17 ms. Both FFT inner layouts lose
despite fewer products. All six variants pass Release, ASAN/UBSAN,
independent packet/interpolation, adjoint, boundary, and alias checks.
Production is unchanged. The round ends for a status discussion before
choosing further work.

The following GRS256/176 results remain the preserved reference, not the
latest best K20 research timing.

The [GRS[256,176] K20 checkpoint](FIELD256_K20_CERTIFICATE.md) closes the
complete 10% distance proof at rate 11/16. Its sparse, intermediate and
dense bounds have integer margins 119, 344 and 12072 bits. Their exact
sum remains below 2^-119. Every nonzero word is covered; this is not a
pointwise screen. Independent sparse and intermediate replays match.

Requested K=2^20 pads to 1,058,816 source symbols and 1,540,096 outputs.
The effective requested rate is 32/47, and the output is 26.5625% shorter
than half-rate SPIN on the same request. The guarantee extends from GF256
to the requested compatible fields without another failure union.

The [native implementation](FIELD256_NATIVE.md) passes independent
scalar/native and ASAN/UBSAN checks. Packed inner broadcasts, an additive-FFT
outer and fused systematic stores improve matched times to approximately
0.279/1.23/7.35 ms at K16/K18/K20. Matched binary controls take about
0.099/0.46/3.16 ms: the higher-rate candidate remains slower. Boolean work,
packing and scratch traffic matter beyond the [product counts](LARGE_OUTER_COST.md).
Next, seek cheaper compatible maps without losing the higher-rate proof
target. K16/K18 and width thirty-two are not certified by the K20 replay.

## Earlier higher-rate checkpoints

The four-state GRS[32,22] candidate is parked after a [verified window
obstruction](WINDOW_DIMENSION_SEARCH.md). At K20, two frozen routes have
shortened source dimensions 136 and 132 even after excluding padded
rows. Thirty-two window exits supply at most 128 homogeneous equations.
Thus every assignment of the inner coefficients has a nonzero word of
weight at most 152064, below the required 152576. These are exact
realized-code statements, not an ideal-setup failure-probability bound.

The search tested three route seeds; unequal-window optimization found
witnesses on two. Equal-window searches found none. The earlier
common-support window bounds did not cover these different row supports.
The compact fixtures replay without optimization and pass an independent
inverse-permutation membership audit. No production profile uses this
research construction.

The [17-state mask screen](MASK_REGIONAL_SCREEN.md) was also implemented
and audited. Its optimized K20 margins at 4/8/16 active rows are
70.870/-49.025/-249.904 bits, worse than the earlier two-state bound in
the latter classes. Separately retaining the feasible total input-weight
range improves that earlier bound by only .00363 bits at eight rows and
.000011 bits at sixteen rows. These diagnostics do not repair the actual
window witnesses.

The next candidates must change outer/state/routing geometry. The
existing [larger-outer cost study](OUTER_LENGTH_TRADEOFF.md) supplies a
rate-preserving GRS[64,44] option, not a certified performance result.
The following class proofs remain valid and are preserved as references.

The latest [Cauchy/diagonal-refresh candidate](DIAGONAL_REFRESH_PROOF.md)
has rate 11/16 and one regional permutation stage. At a 10% output-distance
threshold, directed replay covers every inverse-input weight in [.11N,.24N)
with class margins 202/1833/8535 bits at K16/K18/K20. The existing independent
dense theorem covers [.24N,N]. The remaining interval (0,.11N) is open;
these are not whole-code margins.

Inside that lower interval, a separate [exact one-row
replay](CAUCHY_FIRSTMOMENT_GATE.md) gives 107/105/103 bits at K16/K18/K20.
It uses the new diagonal-refresh distribution and exact GRS support counts.
The new [two-row certificate](CAUCHY_Q2_TRANSFER.md) gives 65/72/77 bits.
The remaining uncovered words have at least three active outer rows and
input weight below .11N. These class bounds do not supply the missing union.

The [joint return analysis](CAUCHY_JOINT_RETURNS.md) now treats two and
three packet inputs exactly. It finds genuine 255^-3 exceptions; each
emits sixteen outputs. Its separate rank-potential interface sums possible
incoming state supports and does not count field-valued witnesses.
Together with a positive selected-zero expansion and exact regional
occupancy, it improves the K20 exploratory margins for three/four active
rows to 102.21/71.98 bits. Eight and more active rows still fail the tested
upper bound. Those floating screens are not additional certificates.

The verified [five-product packet-map circuits](CAUCHY_DIAGONAL_PACKET_COST.md)
give 169/88 GFNI512 processed-byte equivalents per source record including
the outer, input diagonals, and state refresh. The previous eta-scaled
regional estimate was 181/88. No native runtime or speedup is claimed.

For the preserved once-per-region construction, the
[return-pivot proof](INTERMEDIATE_RETURN_CERTIFICATE.md) separately covers
[.19N,.24N), with class margins 646/3616/15647 bits. Exact freshness and
rank-Bellman operators were implemented and tested, but their separated
regional-profile envelopes do not close the remaining sparse class.
The [freshness study](OUTER_FRESHNESS_TRANSFER.md) identifies a concrete
loss: separate regional maxima can select empty profiles and omit the
mandatory active-region span of an outer row.

At that earlier checkpoint, no complete useful one-route certificate had
been established at higher rate or higher distance. The full uniform-route proposal is ruled out by
an exact [suffix obstruction](HIGH_DISTANCE_GLOBAL.md). Regional candidates
must address multi-row cancellations as well as individual-word spread.
No new encoder timing is claimed from the latest proof-only exploration.

The new [dense-class proof](REGIONAL_GLOBAL_ROUTE.md) covers all inverse
inputs of density at least .24 for GRS[32,22], or .25 for GRS[34,24], with
exact finite bounds. The [capacity comparison](REGIONAL_CAPACITY_COMPARISON.md)
selects these as higher-rate proof candidates: one-row and common-support
window events are controlled, but the remaining sparse multirow class is
not. The fixed-circuit [cost comparison](PARITY_STATE_COST.md) is arithmetic
accounting, not a benchmark or a full certificate.

An earlier proof candidate uses a [dual-MDS emission map](DUAL_MDS_PHASE.md).
It changes one four-coordinate reversal, preserving the multiplication
count while making both period maps MDS. Exact local/adjoint checks and
the [joint sparse rank lemmas](REGIONAL_SPARSE_COVER.md) are recorded;
there is still no new complete certificate or native timing.
The [random-packet alternative](RANDOM_PACKET_GATE.md) and
[nonempty-source envelope](NONEMPTY_SOURCE_ENVELOPE.md) remain proof
diagnostics, not promoted designs.

## Preserved checkpoint: a complete half-rate one-route K20 certificate

Independent systematic Toeplitz[64,32] blocks over GF16, one 64-region
permutation, and a three-symbol state give a complete 10% distance
certificate at K20. The failure bound is 4.471426459321588e-25 < 2^-80,
rounded upward here. The proof uses the complete average outer-value
counting measure, so it requires no input-scaling layer. See
[K20_RANK_CONTINUATION.md](K20_RANK_CONTINUATION.md). Separate directed
replays, independently audited, certify the identical construction at K16
and K18 with minimum distances 13108 and 52429 and failure below 2^-84 and
2^-82, respectively. Thus every size in the timing table has a complete proof.

Initial serial phase timings on Ryzen 7950X/GCC 15.2/CPU 2:

| Size | Inner plus route | Outer | Full transpose | Schedule |
|---|---:|---:|---:|---|
| K16 | 108.271 us | 169.424 us | 277.875 us | scatter |
| K18 | 458.530 us | 672.768 us | 1132.771 us | scatter |
| K20 | 3441.513 us | 3359.331 us | 6808.799 us | region-nt-all |

Each row uses seed 113, 51 calls after five warmups, precomputed setup,
128-bit records, and huge-page-preferred buffers. Separate phase medians
need not sum to the median total. All routing, arithmetic, and output
stores are included. Setup is approximately 1.94/7.76/34.08 ms and excluded.
The longer outer accounts for most small-size time and about half at K20.
This initial implementation is not a speed winner; the proof closes a
previously uncovered point.

The outer uses a fixed recursive Toeplitz circuit with 243 GFNI512 products
per batch of four source rows. The dense reference takes 647.331 us for
the K16 outer versus 169.424 us recursively, with identical output. Adding
the now-unneeded input scalars increases K16 full time to 291.500 us in the
same campaign. These two diagonal settings sample different codes; only
the dense/recursive and routing-schedule comparisons preserve the exact map.

Independent mathematical and native-source audits pass. All eight Release
and ASAN/UBSAN targets pass, as do full-size dense/recursive comparisons,
all four accepted alignments, natural-size tails, field adjoints, and
in-place suffix preservation. Production library profiles are unchanged.

## Earlier checkpoint: one routing stage, complete K16 and K18 certificates

The new K18 profile uses six GF16 state symbols. At K=262144, N=524288,
it certifies minimum symbol distance at least 52429, with setup-failure
probability at most 8.272540769181098e-14 < 2^-43. The independent audit
and full directed replay passed. See [GF16_SIZES.md](GF16_SIZES.md).
It retains the same outer and one regional permutation stage as K16.

The first six-state K18 native run takes 1.040 ms (0.877 ms inner/routing,
0.163 ms outer), versus 0.725 ms for four states in the matched phase run.
The faster four-state K18 row does not have a certificate. Six-state setup
is 14.3 ms, excluded from encoding. A new exact six-state FFT variant passes
correctness and sanitizer tests. Its subsequent full-mode holdouts measure
1.0134/1.0148 ms at K18 versus 1.0696/1.0732 ms for direct maps with the same
two seeds. A separate matched phase campaign gave 1.011 versus 1.040 ms.

The four-state kernel also gains a packed-coefficient variant, preserving
the code and K16 proof. Two-seed full-mode holdouts give 0.1763--0.1772 ms
at K16, 0.7390--0.7393 ms at K18, and 5.2818--5.2845 ms at K20. Compared
with `direct` in the same campaign, this saves about 4% at small sizes and
2% at K20. The earlier 5.05 ms result below came from a different campaign;
the new measurement is not presented as an improvement over it.

At that checkpoint, K20 had no complete one-route certificate. Exact small-return spectra,
a shorter inner step, and fresh random feedback did not close the tested
upper bounds. The [progress record](ONE_ROUTE_PROGRESS.md) records these
negative results so they are not mistaken for untried options.

Validation at this checkpoint: 120 Python tests (14 model, 34 lane-mixing,
72 theory) and all seven native targets in both Release and ASAN/UBSAN pass.
The full K18 receipt was replayed independently after the endpoint-integrity
fix. Binary production code and the preserved two-permutation reference
remain unchanged.

## Earlier checkpoint: complete K16 certificate and direct-map baseline

GF16 GRS[16,8], one regional permutation, and a width-16/four-symbol-state
inner now join the fast implementation with a complete finite proof.
At K=65536 and N=131072, the ideal ensemble has minimum symbol distance
at least 13108 except with probability at most 1.763459790707448e-13,
strictly below 2^-42. Scalar extension preserves this guarantee over
GF(2^8), GF(2^32), GF(2^64), and GF(2^128).

The direct-map implementation uses GFNI with an explicit GF16 embedding,
not general large-field products. On the same Ryzen 7950X, GCC 15.2,
CPU 2, serial precomputed 128-bit-record transpose measurements give:

| Candidate | K16 | K18 | K20 | Full finite proof |
|---|---:|---:|---:|---|
| GF16, four-state, direct maps | 0.1857--0.1861 ms | 0.7643--0.7645 ms | 5.0506--5.0638 ms | 10%, 42 bits at K16 |
| GF256, two-state, optimized exact schedule | 0.132 ms | 0.585 ms | 4.08--4.10 ms | Still open |

Small sizes use scatter; K20 uses a region-sized temporary and non-temporal
stores. Setup and page preparation are excluded. The GF16 setup itself is
about 2.3/9.4/39 ms at K16/K18/K20. The final GF16 holdout uses two setup
seeds and 101 measured calls per seed. The earlier baseline campaign measured
0.183/0.760--0.773/5.36--5.38 ms; it is preserved in the progress record.
All timings include routing, arithmetic, output stores, and required fences.
They use the explicit field basis documented by the implementation, with
no conversion from a different consumer basis included.

An exact additive-FFT realization of the inner maps passed correctness checks
but regressed the full encode by approximately 10--14%. It remains an optional
research A/B variant; the direct map is the retained implementation.

The [progress record](ONE_ROUTE_PROGRESS.md) and
[independent audit](SUBFIELD_AUDIT.md) describe the proof refinements,
negative experiments, and remaining work. The older results below remain
useful references; the binary production library has not changed.

## Proved control

The control uses a GRS outer, two uniform symbol permutations, and two
input-weighted accumulators. Each multiplier is sampled from polynomial bit
patterns 1,...,128, then fixed during encoding. These multipliers mix the
native field bits rather than applying the same XOR code independently to bits.
For m>8 they are not restricted to a GF256 subfield.

The [proof and verifier](SMALL_COEFFICIENTS.md) adapt the BAA companion's
shortened-space and quotient-rank argument. No odd-characteristic refinement
or integer Newton outer is imported. Each row bounds the probability that
the sampled code has a bad nonzero message; it is not a fixed-message bound.

| Outer | K | Relative distance | Setup-failure margin |
|---|---:|---:|---:|
| GRS[8,4] | 2^16 | at least 10% | 41.290 bits |
| GRS[8,4] | 2^18 | at least 10% | 45.298 bits |
| GRS[8,4] | 2^20 | at least 10% | 49.300 bits |
| GRS[16,8] | 2^20 | at least 30.5% | 58.865 bits |

Rate is 1/2. Distance counts nonzero symbols, and K counts message symbols.
The bounds apply separately to the specified ideal ensemble over each of
GF(2^8), GF(2^32), GF(2^64), and GF(2^128).

## Full precomputed transpose performance

Measurements use one GF(2^128) symbol per 128-bit record, on Ryzen 7950X
with GCC 15.2, pinned to CPU 2. Both scalar coefficient lookup and the small
CLMUL reduction are enabled. Setup and allocation are excluded; routing,
field arithmetic, and output stores are included. Benchmarks ran serially.

The table gives the range of two per-seed medians, with seeds 113 and 7919.
Each seed uses five warmups followed by 51 calls. Payloads at least 16 MiB
request transparent huge pages. This is a research measurement, not a
production-library performance promise.

| Construction | K16 | K18 | K20 | Proof status |
|---|---:|---:|---:|---|
| Two scans, GRS[8,4] | 0.3684–0.3685 ms | 1.907–1.921 ms | 21.025–21.085 ms | Complete, 10% |
| Two scans, GRS[16,8] | not in this holdout | not in this holdout | 23.751–23.841 ms | Complete, 30.5% at K20 |
| Batched t16, post-refresh, no diagonal | 0.2032–0.2034 ms | 0.8348–0.8386 ms | 7.325–7.369 ms | No whole-code certificate |
| Batched t16, pre-refresh and diagonal | 0.2794–0.2809 ms | 1.138–1.143 ms | 10.396–10.761 ms | No whole-code certificate |
| GF256 GRS[16,8], t16, exact SL2 pre-refresh and diagonal | 0.1326–0.1337 ms | 0.5722–0.5785 ms | 6.766–6.798 ms | Partial K16 bounds only |

The two native-field batched rows use the same GRS[8,4] outer but different inner ensembles.
They are not certified replacements for the first two rows. The no-diagonal
row uses destination prefetch; the diagonal row retains the baseline inner
loop because the explicit prefetch variant regressed in the small-size screen.

The final row is a different, implemented GF256 construction. It interprets
128-bit records in a declared GF256 basis of GF(2^128). It does not use the
native polynomial basis of the preceding rows. No basis conversion is included;
a consumer with a different field representation must account for that cost.

Raw holdout output is ignored at `tmp/char2-field-spin/final-holdout.txt`
and `tmp/char2-field-spin/subfield-holdout.txt`.
The [implementation guide](../../../spin/experiments/char2_field_spin/README.md)
gives the build and invocation parameters. A fresh campaign can reproduce
the experiment without that raw log.

The control's two global routes dominate K20 latency. Earlier phase measurements
put roughly 9–10 ms in each scan and about 2 ms in the GRS8 outer. A gather
implementation of the same map was slower. Expanded coefficient tables helped
K16 but worsened large-size memory pressure. Four-group SIMD processing improved
the batched outer; additional unrolling alone made little difference.

This separates two issues: small-field-polynomial multiplication is feasible,
but a complete field code still needs economical routing and a provable inner.
The recorded ~0.54 ns/record coefficient-only microbenchmark is not an
end-to-end encoder estimate or a dependency-chain latency.

## The implemented GF256 candidate

The outer evaluates a degree-below-eight polynomial at 16 distinct GF256
points. Each coordinate has its own region, with an independently sampled
uniform permutation of the outer-row labels. Independent nonzero GF256
scalars multiply the routed symbols. The inner consumes t=16 symbols X
with a two-symbol state Q:

    R = M Q,
    Y = X + A R,
    Q_next = R + C X.

Here row j of A is (1,j), C=A^T, and every step has an independent uniform
M from SL(2,256), sampled during setup. The state starts at zero, continues
across regions, and is discarded at the end. This is a causal invertible
inner, so the complete encoder has rate 1/2 for every setup.

GFNI performs the GF256 products bytewise. In the declared extension bases,
this is multiplication by a subfield scalar, not an arbitrary binary map.
The encoder is linear over each full field GF(2^m), m=8,32,64,128.
Its symbol distance equals the base-field code's distance; a larger m does
not by itself improve the guarantee. The optimized transpose uses 128-bit
records; scalar forward and transpose references cover all four widths.
No general 32-, 64-, or 128-bit field product appears in encoding.

The exact SL2 sampler parameterizes a nonzero first column and one free
GF256 element. This avoids rejection over complete matrices. Encoding uses
four independent small products for the matrix update, without the three
serial shear dependencies of the earlier native-field candidate.

In the initial phase run, inner arithmetic plus routing cost approximately
0.090/0.395/5.61 ms at K16/K18/K20; the outer cost 0.045/0.185/1.19 ms.
Setup was about 1.4/5.8/22 ms and remains outside precomputed timings.
These measurements identify inner plus routing as the next large-size target.

A matched prefetch variant preserved every checksum but regressed K16/K18
by roughly 5–6%. Its K20 improvement was below 1%. The baseline remains
the default; no size-dependent dispatch was added for that small gain.

## What the batched proof established

The implemented native-field candidate routes each outer coordinate to its own
region. Within each region, it independently permutes individual outer-group
labels. It does not move intact t-row cohorts. The inner holds two field symbols
and processes t=8 or t=16 input symbols together.

Three small-coefficient shears give an invertible state update. For a fixed
nonzero state, a prescribed target has probability at most Q^-2, where Q is
the coefficient-set size. A rank-two matrix of state forms has target atoms
at most Q^-3. These are local bounds, not a complete distance theorem.

A second, now implemented construction uses GF256 coefficients throughout, independent
nonzero input scales, and a uniform SL(2,256) state refresh before emission.
It has an exact joint transition for output weight and return to zero.
The transition initially depends on the active lane positions; uniform regional
routing permits averaging those positions before composing the regional operators.
The independent [audit](lane_mixing/GF256_JOINT_AUDIT.md) checks this distinction.

Scalar extension transfers a base-field distance bound to the four declared
fields. The separate final timing row above measures this implementation.

The first one-active-outer-row screen at K16 and 10% distance gives:

| GF256 model, t16 | One-row bound |
|---|---:|
| GRS[8,4] outer | approximately 19.88 bits |
| GRS[16,8] outer | approximately 49.79 bits |

These original screens motivated the larger outer. The new
[outward calculation](SUBFIELD_PROOF.md) proves an upper failure contribution
of 1.031016e-15 for one active row and 7.830859e-14 for two active rows.
Their sum is below 2^-43. The [independent dense bound](SUBFIELD_AUDIT.md)
covers all active-row counts 2879 through 8192 jointly below 2^-64.
These are partial bounds, not a complete 10% certificate.

There is now a precise obstruction to completing the same word-counting
argument. Some rare routes confine selected rows to short regional windows.
Only a few state-zero equations then leave a large subspace of low-output
messages. Counting every projective line in this subspace makes the expected
bad-word count exceed 2^80. This does not show a comparable code-failure
probability: many witnesses occur together in the same rare setup.

Increasing the state to four or eight symbols does not remove this counting
obstruction. Removing the approximate outer-support envelope also recovers
less than half a bit in the tested multirow classes. Neither is the right
next implementation change. The BAA-style quotient-rank proof avoids counting
all combinations of witnesses. Fresh-diagonal and SL2 return lemmas are
recorded in the proof note, together with conservative active-input bounds.
Their full regional composition remains unfinished.

For the current A/C maps, t16 also has a structural distance ceiling near 25%.
It cannot be the route to a near-50% half-rate design merely by enlarging the
field. Width t8 raises this ceiling near 50%; that does not prove it is attained.
See [BATCHED_PROOF.md](BATCHED_PROOF.md) for the exact constructions and scope.

## Validation and next step

Native scalar/SIMD differential checks, full-field linearity and adjoints,
all implementation variants, and ASAN/UBSAN checks pass. Separate exact tests
cover field irreducibility, GRS shortening, cancellation bounds, finite-certificate
replay and tampering, and the small-field local transition formulas.
At this checkpoint, all 78 Python tests pass across the three suites;
all four C++ targets pass in both Release and ASAN/UBSAN builds.
The binary production profiles are unchanged.

Next, implement and evaluate the conservative quotient-rank transition while
retaining the rank of the two state forms. Compose it with the regional support
count before changing the kernel again. The ordinary fixed-message transfer cannot simply
be reused after conditioning on earlier zero equations. Keep the two-scan
construction as a certified control and leave the binary production profiles
unchanged until the faster candidate has a complete proof and a measured use case.

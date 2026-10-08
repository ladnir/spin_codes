# SPIN over characteristic-two fields

This workstream starts a new code family over F = GF(2^m), for
m in {8, 32, 64, 128}. It does not replace the binary library profiles.
The completed RS8 library promotion is commit `3ce07a64`.

We want either less encoding work at 10% distance or substantially greater
distance at a comparable rate. Here K and N count field symbols,
R = K/N is the rate, and weight counts nonzero field symbols, not bits.
For the first track, less expansion means a higher R. The active comparison
therefore seeks higher rate at 10% distance and competitive full encoding
cost, rather than another half-rate scalar-extension guarantee.

The initial rate-1/2 high-distance target is close to 50%.

## Current direction: high distance

The user has redirected the search to the
[blockwise SPIN restart](HD_BLOCKWISE_RESTART.md): fast local outer,
one routing stage, and a batched stateful inner, without a pre-accumulator.
The target remains K20, rate 1/2, at least 40% symbol distance with a
40-bit setup-failure margin, and full forward encoding below 5 ms.

The proof target is now closed for a
[sixteen-state refinement](HD_SUMRANK_FEEDBACK_CERTIFICATE.md): 40% distance
at rate 1/2 and more than 59.87 bits of setup-failure margin. It uses
sum-rank-separated fixed maps and independent feedback randomization.
The directed replay covers all 8192 active-row counts. This construction
is preserved as the proof reference. The
[encoding optimization](HD_ENCODING_OPTIMIZATION.md) now has a native
[shared-diagonal RS inner](HD_SHARED_DIAGONAL.md) at about 6.1 ms, with a complete
40%-distance certificate above 49.33 bits. The independent-map
[6.3-ms checkpoint](HD_SYSTEMATIC_MDS.md), with 50.06 bits, is preserved.
The prior 7.4-ms layered LRS
version, its [FFT state maps](HD_FFT_STATE_MAPS.md), and its 48.32-bit
certificate remain preserved. The
sub-5-ms target remains open; this is research code, not a library promotion.

A [block-feedback refinement](HD_BLOCK_FEEDBACK.md) now reaches about
5.5 ms with streaming output, retaining a complete 256-bit directed
certificate above 51.04 bits
at 40% distance. It alternates four/eight-check local feedback, retaining
global RS emission. The 104 rational witnesses are pinned; a separate
384-bit replay confirms all 8192 row counts and 51.04261653 bits of margin.
The streaming implementation computes the same map; the prior
5.7--5.8-ms cached-output checkpoint remains preserved. All 45 native
variants pass Release and ASAN/UBSAN checks.
Pure four-check feedback is about 5.6 ms
but remains uncertified. [Cheaper block emission](HD_BLOCK_EMISSION.md)
now has a certified rank-two alternative, described below. A
[common-map refresh](HD_COMMON_REFRESH.md) reaches about 5.3 ms, but its
local proof has not yet yielded a complete distance certificate.
The [return-bound investigation](HD_RETURN_DESIGN_SCREEN.md) records the
corrected grid-rank bound and failed basis/scalar/state-refresh screens.
A rank-two refresh with eight-check feedback measures about5.6 ms with
cached output and 5.4--5.5 ms with streaming output.
An audit corrected its grouping of complementary transition probabilities.
The corrected 256-bit search now certifies all 8192 classes with more than
50.9966 bits at 40% distance. Its 133 rational witnesses are pinned, and a
separate 384-bit replay confirms the complete certificate.
Earlier, pre-correction partial receipts
remain withdrawn. It is a modest alternative, not a sub-5-ms result.
Exact-map route-layout, address-precomputation, and shared outer-load
controls have not yielded another repeatable improvement. Compacting the
per-step setup also gave no repeatable win. The period-four candidate uses
eight-check feedback only once per four packets and takes about 5.45 ms.
Its complete certificate gives more than 56.6591 bits at 40% distance.
All 119 rational witnesses are pinned; a separate 384-bit replay confirms
all 8192 classes. The simpler schedule has not established a substantial
timing advantage over the other two independently replayed alternatives.
The [optimization checkpoint](HD_OPTIMIZED_CHECKPOINT.md) consolidates the
three certified designs, final K20 validation, reproduction commands, and
the remaining limits. No sub-5-ms or production result is claimed.

The current baseline adds [independent GF65536 outer-output scalars](HD_RANDOMIZED_OUTER.md),
restoring the expected-measure comparison used in the binary packet proof.
The checked prototype costs 7.02ms versus 6.53--6.55ms for its identity
control. Directed [class bounds](HD_RANDOMIZED_OUTER_PROOF.md) close
all low-weight words with at least7548 active rows, and sparse/dense quiet
classes. The [exact sparse ladder](HD_RANDOMIZED_OUTER_SPARSE_LADDER.md)
now covers all words with one through five active rows, including nonquiet
histories, with a union below 2^-56.51. The
[quiet extension](HD_RANDOMIZED_OUTER_QUIET_EXACT.md) covers active-row
counts 1 through 962 and 5865 through 8192 using the combined receipts.
The [occupied-tile state moment](HD_RANDOMIZED_OUTER_FLAT_MOMENT.md) and
[composed local refinements](HD_RANDOMIZED_OUTER_LOCAL_REFINEMENTS.md)
improve the proof bounds without changing the encoder. The directed
[joint-input and coefficient-cap receipts](HD_RANDOMIZED_OUTER_JOINT_INPUT.md)
extend quiet coverage to 1508 through 8192. Quiet counts 963 through 1507
and nonquiet counts 6 through 7547 remain open.
The complete certificate for that four-state native construction remains open.
The latest [tile-conditioned proof](HD_RANDOMIZED_OUTER_QUIET_TILE_CONDITIONING.md)
recovers substantial loss from averaging packet layouts before fractional
powers, but its middle witnesses remain unclosed. A
[support-floor experiment](HD_RANDOMIZED_OUTER_QUIET_SUPPORT.md) identifies
why a coarse input-weight penalty does not resolve them. The algebraically
valid [FFT16 alternative](HD_RANDOMIZED_OUTER_FFT16_PROBE.md) also fails the
tested quiet bounds, so it is not selected or implemented. The next proof
step is to retain state histories longer in the composition. A distinct
[sum-rank-separated inner](HD_SUMRANK_INNER.md) has exact local mask ranks
and input/output moments. Its four-state quiet class is certified with
93.99 bits of margin, but its nonquiet classes remain open. The larger
feedback-randomized refinement above is the fully certified design. The
[closure goal log](PROOF_CLOSURE_GOAL.md) separates accepted coverage from
failed-bound diagnostics. All native timings above are
unchanged; this iteration performed proof computations only.
The [randomizer audit](HD_RANDOMIZED_OUTER_OPTIONS.md) identifies which
weaker distributions can preserve these proof statements. No extra local
shuffle is required at the present routed-symbol granularity.

The following fixed-outer work is preserved as a checkpoint. The
[shared eight-point circuit](HD_BLOCKWISE_SHARED_MAPS.md) adds MDS input
mixing while reusing its transforms for state emission and feedback. Its
exact algebra and structural checks pass. The new
[cube circuit](HD_BLOCKWISE_SHARED_FACTOR.md) shares intermediate XOR sums
across feedthrough, emission, and feedback:59 packed GFNI arithmetic calls
per64 outputs, versus82 for the shared FFT version. The
[native comparison](HD_BLOCKWISE_SHARED_NATIVE.md) now rejects this schedule
as a speed improvement:7.65ms versus6.54ms for the same-map FFT, with a
whole-vector rewrite at7.00ms. The working shift64 map removes a specific
[paired cancellation](HD_BLOCKWISE_SHARED_SUPPORT.md) present at shift8.
The [state-budget bound](HD_BLOCKWISE_SHARED_BUDGET.md) rules out excessively
infrequent state updates, and a [sparse quiet-class proof](HD_BLOCKWISE_SHARED_CUTS.md)
bounds every quiet word using at most256 outer rows below2^-2236.70.
An [all-density type bound](HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md) closes the
high-information quiet class below2^-80.87, and the
[joint-image transfer](HD_BLOCKWISE_SHARED_TRANSFER.md) caps whole-subspace
returns under the actual rank-one mixer. Low-information quiet words and
nonquiet histories remain open for that fixed-outer ensemble. The selected
next step is the randomized-outer proof above, rather than further tuning
of the current cube circuit.

The latest [direct-SPIN native control](HD_BLOCKWISE_NATIVE.md) measures
about6.35ms forward at K20: GRS[256,128], one regional route, and a
width16/state4 inner. It has no full distance certificate and is not promoted.
The first width64/state16 control remains frozen at about9.0ms. Its
[fixed maps and local state bound](HD_BLOCKWISE_MAPS.md) are checked; its
[capped global routing bound](HD_BLOCKWISE_GLOBAL.md) handles one complete-row
trap family. Those global numerical receipts do not automatically transfer
to the smaller-packet construction.

The [quiet-subcode bound](HD_BLOCKWISE_QUIET_SUBCODE.md) now rejects the
width64/state4 and width64/state8 controls for every setup, including with
the larger256-symbol outer. Their timings are cost controls, not viable
40%-distance results. Shrinking packet width together with state avoids
this particular obstruction. The [capped support composition](HD_BLOCKWISE_COMPOSITION.md)
closes all one-active-row messages, while many-row certification remains open.
For the selected t16/s4 geometry, the [interval-span bound](HD_BLOCKWISE_INTERVAL_SPANS.md)
also puts one dimension-forced complete-row trap family below 2^-528.39.
It counts routes once, avoiding a demonstrated uncapped-first-moment failure.
Partial-row shortening and general output histories are not covered.

The [accumulator detour is frozen](ACCUMULATOR_DETOUR_CHECKPOINT.md), with
a verified source archive, measured 8.3-ms control, and all proof results
indexed. It has no complete distance certificate and is not the active
implementation direction. The narrative below records that history.

The following entries describe the preserved pre-restart exploration.
The [high-distance speed goal](HIGH_DISTANCE_SPEED_GOAL.md) targets
40--45% distance below5ms at K20, then about3ms and higher distance, against
the16.6ms forward-RS baseline. Its first exact-map forward control is16.8ms;
no new complete distance certificate is claimed. Construction/proof changes
are needed for a substantial advantage.

The next [two-scalar-scan control](HD_STRIPED_NATIVE.md) reached about 12 ms
but is withdrawn: the [grouped-check bound](HD_BAA_PROOF.md) rules out its
40% target for every setup. Stronger [row-order cuts](HD_REGIONAL_THREE_PROOF.md)
also reject the next three-scan proposal, and a [run cut](HD_RATIONAL_POST_OBSTRUCTION.md)
rejects its rational scalar post-inner. The later
[two-state interval cut](HD_TWO_STATE_CUTS.md) rejects the proposed
[two-additive-state construction](HD_ADDITIVE_TWO_STATE.md) with setup-failure
probability greater than99.8%. Its local laws and
[regional proof machinery](HD_TWO_STATE_REGIONAL.md) remain indexed for reuse.
That search next prioritized fixed affine source-row placement. A directed
spectral bound disables the known quarter-row cut, without proving distance.
Its [forward cost checkpoint](HD_AFFINE_NATIVE.md) reaches 8.327--8.360 ms
with four-row affine groups and streaming stores. Grouping changes the
construction; the [bounded cut audit](HD_AFFINE_CUTS.md) finds no retained
counterexample but supplies no distance certificate. The
[exact source/rank interface](HD_AFFINE_PROOF.md) identifies the remaining
all-support proof requirement.
The [component contraction](HD_COMPONENT_RANK.md) now bounds the entire
permitted-singleton term below 2^-84.93, and the
[permitted-pair bound](HD_COMPONENT_PAIRS.md) below 2^-41.11. These are
positive expansion terms, not a full certificate. The
[dense tail](HD_COMPONENT_DENSE.md) now bounds all classes a>=29943 below
2^-322.58, but also proves that rare routes defeat the uncapped actual-code
first moment in the middle range. Exact MDS filtering cannot fix that issue.
The new [routing-condition gate](HD_ROUTING_CONDITION.md) bounds a family
of these rare structural events below 2^-155.23, uniformly over large row
sets and low-weight outputs with at most 64 runs. For that parked candidate,
further proof work would need to control rank-defect probability, not count
every line in rare large kernels.
The [MDS-aware interface](HD_MDS_COMPONENT.md),
[joint component geometry](HD_COMPONENT_GEOMETRY.md), and
[sparse/offset refinements](HD_SPARSE_REFINEMENTS.md) are preserved as
partial tools; none supplies a full certificate.
The shared upper operator also covers a [three-product post-inner](HD_THREE_PRODUCT.md),
which saves one field product but has not been benchmarked.
The [larger vector-state alternative](HD_VECTOR_OBSERVATION_PROOF.md) has
checked local laws but failed tested middle-class bounds. Its
[rank-refinement review](HD_VECTOR_SOURCE_REFINEMENT.md) rules out one
proposed repair and identifies useful transitions for cheaper inners.
The [composable Cauchy rank refinement](HD_CAUCHY_RANK.md) improves its
saved bounds only modestly; it does not close their large deficits.
Further [grouped cut diagnostics](HD_GROUPED_CUT_FOLLOWUP.md) found no new
counterexample, without bounding finite failure probabilities.
Complete high-distance certification and the sub-5-ms target remain open.

The user accepted the [3.2 ms high-rate winner](GF8_HIGH_RATE_WINNER.md).
Its implementation, proof dependencies, and research notes are preserved in
a verified source-only archive. The 83.33% rate, 10% distance, and 57-bit
K20 certificate are unchanged. Further high-rate tuning is parked.

The active search now targets high symbol distance, ideally close to 50%,
starting at rate 1/2 and K20. Keep 128-bit application records, one regional
routing stage, and a 40-bit whole-code setup-failure target. These are search
criteria, not a new certificate. The coefficient alphabet is open again:
wider records alone do not strengthen a GF256 scalar code's distance.
The [first high-distance screens](HIGH_DISTANCE_RESTART.md) identify a
36% GF256 candidate at rate 1/2 and distinguish it from the larger-alphabet
path toward 50%. Complete certification and native timing remain open.

The [larger-field initial results](HIGH_DISTANCE_INITIAL_RESULTS.md) now
identify promising 42% GF(2^16) and 45% GF(2^32) controls, still without
complete certificates. New validated multiplication kernels make the
larger fields cheaper; the proposed encoder's eightfold raw expansion
remains a major cost concern. The subsequent
[lower-expansion exploration](HIGH_DISTANCE_LOWER_EXPANSION.md) prioritizes
a 4K/[256,64] candidate with positive 40% screens, before completing another
certificate. It also proves a genuine window obstruction for the smaller
outer's direct-2K alternative. The 49% middle bound remains open.

## Preserved GF256 proof and encoding cost

The new [speed-first profile](GF8_SPEED_DESIGN.md) uses GRS[64,40], a
four-input/four-state inner, and three retained outputs. It closes K20 at
83.3333% nominal rate, 10% distance, and 57-bit setup-failure margin.
The [latest native comparison](GF8_SPEED_ROUTED.md) selects relocated
diagonals, sixteen-packet scheduling, a fused outer, and a padded scratch
stride: 3.184--3.201 ms at K20. Restoring input before every timed call
gives 3.332--3.361 ms. This is about20% less time than the matched
[eight-packet checkpoint](GF8_SPEED_NATIVE.md), without changing the code
or its certificate. Production is unchanged.

The [next profiling round](GF8_SPEED_PROFILE.md) splits the remaining work
into about1.25ms inner computation,1.17ms routing, and0.84ms outer work.
Sixteen further exact-map schedules did not produce a material improvement;
the3.2ms checkpoint remains selected.

The [rate frontier](GF8_RATE_FRONTIER.md) now reaches nominal rate 37/44
(84.0909%) at 10% symbol distance with a 183-bit whole-code setup-failure
margin at K20. This is 0.0537 percentage points below the GF256 random-code
benchmark; requested K20 rate is 83.8657% after padding. Four complete new
checkpoints cover fixed and periodic output retention. The strongest proof
uses the [exact MDS spectrum](GF8_EXACT_SPECTRUM_DENSE.md) in the dense class.
There is still one regional route and the same eight-input/eight-state maps.
The [native comparison](GF8_RATE_NATIVE.md) selects the 84% profile:
4.807--4.849 ms for full precomputed K20 transposed encoding of 128-bit
records, with 115-bit certified margin. Fixed-six takes 4.903--4.918 ms;
84.09% takes 5.318--5.350 ms for only 0.3265% fewer retained records.
Folded routing, page-aligned column strides, and fused outer loads/stores
preserve the sampled maps. Production and frozen proofs are unchanged.

The [three-quarter-rate goal](GF8_THREE_QUARTER_GOAL.md) now has a complete
[K20 certificate](GF8_PUNCTURED_WHOLE_CERTIFICATE.md): nominal rate3/4,
10% symbol distance, and116-bit whole-code setup-failure margin. A stronger
GRS[128,84] outer feeds a full-Cauchy8/8 inner; keeping seven of every eight
outputs supplies the rate increase. There is still only one regional route.
The [native research encoder](GF8_PUNCTURED_NATIVE.md) measures about
.240/.938/5.57--5.82ms at K16/K18/K20, with13.38% fewer K20 outputs than
the previous cubic field control at roughly comparable encoding cost.
Only K20 has the new whole-code certificate. No production profile changes.

The earlier unpunctured [full-Cauchy maps](GF8_MDS_PACKET.md),
[partial-live return lemma](PARTIAL_LIVE_RANK.md),
[excursion diagnostic](GF8_EXCURSION_DIAGNOSTIC.md), and [joint source
screens](GF8_JOINT_SOURCE.md) remain preserved. Larger states and wider packets alone did not
close the intermediate-support gap; local puncturing changes that tradeoff.

The [native eight-record implementation](GF8_CUBIC_NATIVE.md) now passes
scalar-reference, large-size, alias, and sanitizer checks. Its first packed
schedule regressed; pairing consecutive packets removes an artificial
state dependency and reaches .237 / .94 / 5.71--5.86 ms at K16/K18/K20.
This is a 14--16% exact-map improvement, but the construction has lower
rate than the previous field profile and remains slower than binary SPIN.
Keep both as research controls; no production change.

The user explicitly wants rates well above 2/3. The [larger-outer
screens](GF8_HIGH_RATE_OUTER.md) show that doubling the outer alone does
not reach 3/4 under the present bound. The GF256 random-code reference at
10% distance is about84.1% rate, not a guarantee for this structured family.
The earlier [inner/return investigation](GF8_HIGH_RATE_INNER.md) tests sharper
rank history, pair mixing, and width8/state6 maps. None of those unpunctured
candidates closes3/4; the newly certified punctured code above is different.

The [new eight-record inner](GF8_PACKET_CERTIFICATE.md) closes a cheaper
arithmetic candidate at K20: GRS[128,83], requested rate 64/99, 10% symbol
distance, and 65-bit whole-code setup-failure margin. A repaired pair of
cubic maps gives a proposed 18 rather than 70 GFNI instructions per sixteen
outputs, at the cost of 4.762% more output records. All 46 new tests and the
complete directed replay pass. Its native result above does not improve the
rate frontier; the certificate and existing profiles remain intact.

The latest [fixed-inner rate sweep](GF8_OUTER_RATE_SWEEP.md) closes a new
GRS[128,88] K18 point: 10% distance, 72-bit whole-code margin, requested
rate 32/47, and 2.083% fewer outputs. K20 dimension88 remains unclosed in
the sparse classes, and dimension89 has middle-class gaps. No new native
timing or production change accompanies this certificate. The timings
below remain for dimension87.

The [GF256 follow-up](GF8_EXPLORATION.md) finds a cheaper outer by changing
the source basis, preserving the entire certified full codeword set.
Fusing the source loads brings full-encoder gains to about 15% at both
K16 and K18 (.244/.956 ms), but both new schedules regress at K20.
The 10% distance, requested rates about 2/3, and
63/79/67-bit certificates are unchanged. This is not yet competitive with
binary half-rate encoding, and does not silently replace its interface.

The parallel [rank-history proof work](GF8_RANK_HISTORY_SCREEN.md) tightens
several bounds but leaves the six-state candidate unclosed. It is a
numerical screen with an independently audited transfer, not a new
whole-code certificate. Existing proofs and production remain intact.

### Preceding medium-alphabet comparison

The latest [GF256/GF32 comparison](MEDIUM_ALPHABETS.md) keeps 128-bit
application records and the useful higher-rate/10% target. Two genuine
GF32 arithmetic implementations pass native and sanitizer checks, but
their isolated coefficient pass is about four times slower at K16/K18.
The [paired proof screens](MEDIUM_FIELD_PROOF.md) do not move the sampled
rate/state frontier: larger coefficients barely help the middle classes.
Keep GF256 as the performance baseline and target the middle-class
analysis or source/packet geometry next. No new full encoder or certificate
is claimed, and no production profile changes.

## Preserved higher-rate checkpoint

**Higher-rate family closed at all three sizes:** GRS[128,87] has complete
10% certificates with 63/79/67 bits of setup-failure margin at K16/K18/K20.
Requested rates are 2/3, 2/3, and 128/189. The new
[K16/K18 replays](FIELD128_SIZE_CERTIFICATES.md) sharpen finite conditioning;
the [K20 certificate](FIELD128_K20_CERTIFICATE.md) remains unchanged.

The [validated native prototype](FIELD128_NATIVE.md) takes approximately
0.288/1.128 ms at K16/K18. Latest K20 matched launches range from 5.73 to
7.46 ms, versus about 3.17 ms for binary SPIN. A compact-route option helps
in direct A/B tests, but the process-level variation is unresolved. The
earlier 6.01--6.10 ms checkpoint is preserved. This is a proved rate benefit,
not a demonstrated encoding-speed win. No production profile changes.

Both exact-map FFT inners were measured and rejected, including a new
lane-local layout: fewer field products do not compensate for their other
costs. The dense inner stays selected. All six current variants pass native
and sanitizer checks. The round is complete for the requested status discussion.
[Design screening](CHEAPER_DESIGN_STATUS.md) records reduced-state and
high-distance gaps, plus actual obstructions for specific cheap mixers.

### Preserved larger-outer checkpoint

**Higher-rate certificate closed:** GRS[256,176], one regional routing
stage, and the symmetric width-sixteen/eight-state inner certify 10%
symbol distance at K20 with at least 119 bits of whole-code setup-failure
margin. The padded rate is 11/16; the requested-input rate is 32/47.
The complete sparse/intermediate/dense replay and scalar-extension claim
are indexed in [the certificate](FIELD256_K20_CERTIFICATE.md).

This is a rate improvement, not a new half-rate/10% result. The
[native implementation](FIELD256_NATIVE.md) passes independent correctness
and sanitizer checks. Exact-map optimizations reach about 0.279/1.23/7.35 ms
at K16/K18/K20, but matched binary SPIN remains faster. Its K20 output is
26.5625% shorter; outer and inner costs consume that benefit. The next step
is a cheaper compatible design, not margin inflation. K16/K18 and width
thirty-two remain open.

### Earlier higher-rate attempts and preserved evidence

The four-state GRS[32,22] profile remains parked. Two
replayed K20 routes force outputs below 10% for every assignment of the
inner coefficients. Unequal regional windows leave shortened source
dimensions 136 and 132 after removing padded rows, but their zero-state
exits impose at most 128 equations. See the [proof and compact
witnesses](WINDOW_DIMENSION_SEARCH.md). This is a realized-code failure,
not merely a loose upper bound. The sampled routes do not by themselves
give a theorem about the ideal ensemble's failure probability.

The replacement had to change outer redundancy, state capacity, or
routing geometry. More randomness in the same four-state inner cannot
repair these routes. GRS[64,44] preserves rate 11/16; GRS[33,22] is a
rate-2/3 alternative. Neither was certified or promoted. The [progress
record](ONE_ROUTE_PROGRESS.md) records this search and preserves all
earlier proof results.

**Objective correction:** the half-rate/10% results below reproduce a
guarantee already available by lifting binary SPIN. They are preserved
controls, not the intended improvement. Further work is gated on higher
rate at 10% distance or substantially higher distance at rate 1/2, with
competitive total cost. See [active targets](ACTIVE_TARGETS.md).

An earlier structural result is a [whole-code obstruction](HIGH_DISTANCE_GLOBAL.md):
a full uniform permutation cannot rescue small outer blocks with a causal
inner. The proposed half-rate/30% instance fails with overwhelming probability
at every requested size. Regional routing avoids that particular obstruction.
Higher-rate regional designs are under proof review. The GRS[32,24]
candidate with state mixing only at region boundaries also fails forty
bits when injection labels repeat every sixteen or thirty-two positions.
Sixty-four labels fix the tested pair mechanism, but a larger shortened
subspace still forces failure probability above 2^-38 at K20. This rules
out the four-state GRS[32,24] profile independently of the injection period;
see [the design review](RATE_DESIGN_REVIEW.md). The next design must change
outer redundancy, state capacity, or routing. The larger GRS[256,176]
checkpoint above supplies the first complete higher-rate certificate.

One reusable positive result is a [proved dense class](REGIONAL_GLOBAL_ROUTE.md)
that holds for any independently chosen invertible inner. The
[capacity comparison](REGIONAL_CAPACITY_COMPARISON.md) now favors [32,22]
and [34,24] four-state outers for the remaining sparse proof. Their rates
are 11/16 and 12/17, with exact one-row and broad window-event bounds, but
no complete multirow coverage. The [cost study](PARITY_STATE_COST.md) compares
their fixed circuits without claiming native timings.

The current [dual-MDS repair](DUAL_MDS_PHASE.md) reverses four state
coordinates before emission, with no extra field products. It strengthens
the period-64 emission geometry and retains the existing conservative
one-row bound. The [sparse proof checkpoint](REGIONAL_SPARSE_COVER.md)
records new joint rank lemmas; complete sparse coverage is still open.
The alternative [fresh random packet maps](RANDOM_PACKET_GATE.md) have
an exact local law but no competitive complete certificate. Their
[cheap source envelope](NONEMPTY_SOURCE_ENVELOPE.md) was tested and parked.

The now-parked candidate uses [Cauchy-form packet maps and fresh
state diagonals](DIAGONAL_REFRESH_PROOF.md), still at rate 11/16 with one
regional routing stage. The rank-moment proof now certifies every inverse
input of weight at least .11N. Separate directed one-row and
[two-row bounds](CAUCHY_Q2_TRANSFER.md) give 107/105/103 and 65/72/77 bits at
K16/K18/K20. The sparse class with at least three rows remains open. Its verified
[four-point circuits](CAUCHY_DIAGONAL_PACKET_COST.md) reduce the estimated
full product volume from 181/88 to 169/88 per source record. This is a
changed construction and an arithmetic count, not a measured speedup or
a complete distance certificate. The prior once-per-region version and
its [.19N-and-above proof](INTERMEDIATE_RETURN_CERTIFICATE.md) are preserved.

The latest [sparse rank transfer](SPARSE_MULTIROW_REVIEW.md) retains exact
regional occupancy and uses the [joint return equations](CAUCHY_JOINT_RETURNS.md)
without counting field-valued messages. Its three/four-row screens are
promising, but larger classes remained unresolved. The evaluated
[actual-state-support transfer](MASK_REGIONAL_SCREEN.md) does not improve
those larger bounds. The window obstruction now takes priority over
further proof refinement for this unchanged profile. No replacement
encoder is promoted.

### Preserved half-rate controls

The archived one-permutation work has a single construction certified at all
three sizes: independent
systematic Toeplitz[64,32] blocks over GF16, one 64-region routing stage, and
a width-16/three-symbol-state inner. It certifies at least 10% symbol distance
with setup failure below 2^-84 / 2^-82 / 2^-80 at K16/K18/K20, respectively.
The proof also removes the input-scaling layer.
See the [construction, complete proof, and replay](K20_RANK_CONTINUATION.md).
The first native timings are about 0.278/1.13/6.81 ms at K16/K18/K20 on
128-bit records. The longer outer is expensive. These certificates close
that ensemble's analysis, but reproduce an existing binary rate--distance
guarantee and are not a useful result for the revised objective.

The earlier complete K16 and K18 certificates are preserved. A GRS[16,8]
outer over GF16, one regional routing stage, and a width-16/four-symbol-state
inner achieve rate 1/2 and minimum symbol distance at least 13108 out of
131072, except with probability below 2^-42 over ideal setup. The inner uses
independent uniform GL4(GF16) refreshes and independent nonzero input scalars.
Scalar extension gives the same distance over all four requested fields.
At K18, six state symbols give minimum distance at least 52429 out of
524288, except with probability below 2^-43. The outer and single regional
route are unchanged; see the [K18 proof and replay](GF16_SIZES.md).
See the [complete proof and certificate](GF16_PROOF.md),
the [progress record](ONE_ROUTE_PROGRESS.md), and
[independent audit](SUBFIELD_AUDIT.md).

The packed four-state precomputed transpose takes about 0.177/0.739/5.28 ms
at K16/K18/K20 on 128-bit records. Only its K16 point has a certificate.
The certified six-state K18 profile measures about 1.014 ms using its exact
FFT alternative, a 3--5.5% gain across matched campaigns. K20 remains
unclosed for this smaller outer. See the [matched campaigns](RESULTS.md) for the
distinction between these timings and the earlier 5.05 ms result. Small sizes use
scatter; K20 uses region-sized reordering and non-temporal stores. A separate
GF256/two-state candidate reaches 4.08--4.10 ms at K20 but remains unproved.
The [implementation](../../../spin/experiments/char2_field_spin/README.md)
keeps both candidates and exact-map A/B schedules separate from production.

## Preserved reference constructions

The first complete construction uses small polynomial coefficients, a GRS
outer, and two input-weighted accumulators. It needs no full-width field
multiplication. The [finite proof](SMALL_COEFFICIENTS.md) certifies 10%
distance at K16/K18/K20 with 41.29/45.30/49.30 bits of setup-failure margin.
A larger [16,8] outer certifies 30.5% distance at K20 with 58.87 bits.
These bounds apply separately over GF(2^8), GF(2^32), GF(2^64), and GF(2^128).

The [implementation](../../../spin/experiments/char2_field_spin/README.md)
contains a scalar reference and optimized precomputed GF(2^128) transpose.
The two global permutations remain expensive. An earlier implementation uses
a GRS[16,8] outer, GF256 coefficients, regional routing, and a two-symbol
batched inner. It measures about 0.133/0.575/6.78 ms at K16/K18/K20 with
128-bit records. Its exact field bases and setup distribution are implemented
and tested; it still has no whole-code certificate.

At K16 and 10% distance, outward bounds cover one and two active outer rows
jointly below 2^-43. A separate bound covers the dense tail below 2^-64.
The middle remains open. An exact window argument shows that counting all
bad words cannot close this construction: rare setups contribute large
subspaces of bad words. The next proof must count independent constraints,
as in the BAA companion. See [proof status](SUBFIELD_PROOF.md), the
[independent audit](SUBFIELD_AUDIT.md), and [measurements](RESULTS.md).

The sections below preserve the initial design analysis. The original proposal
to use full-field GL matrices has been superseded by the user's requirement
to avoid large multiplications. Its algebra remains useful as a reference.

## What the companion contributes

The adjacent [field-universality project](../../../../gen-BAA/field-universality/README.md)
contains three useful ingredients:

1. [Scalar-extension invariance](../../../../gen-BAA/field-universality/paper/main.tex):
   a generator over a subfield has the same minimum symbol distance over
   every extension field, for each fixed setup.
2. [Shortened-space and equation-rank accounting](../../../../gen-BAA/field-universality/paper/distance-proof.tex):
   bound the existence of a bad message by counting support branches and
   their independent zero equations, rather than all field-valued messages.
3. [A convolution rank-mask transfer](../../../../gen-BAA/field-universality/expander-lift/CONVOLUTION_RANK_MASK_TRANSFER.md):
   fresh independent coefficients from a set of Q distinct nonzero elements
   cancel r independent forms with probability at most Q^(-r).

These references describe the local checkout inspected on 2026-10-03.
Its HEAD was `c6613d1fd62f568bfbe94f033fcef339a7196d7b`, with additional
uncommitted manuscript and proof work. That commit alone is not a snapshot
of the inspected files.

The rank and pivot arguments are characteristic-independent. The companion's
integer Newton outer and its opposite-pair refinements are not. For example,
the Newton [16,8] map reduced modulo two sends the eighth basis input to
a weight-two output. Use a genuine characteristic-two GRS constituent
with distinct field evaluation points instead.

Its global permutation formulas also cannot be imported unchanged: a SPIN
region route has different support-placement probabilities. Similarly,
the independent-coefficient pivot bound does not apply directly to the
dependent entries of a uniformly sampled invertible matrix.

## Coefficient field versus storage width

Using the current XOR-only encoder on wider records does not improve distance.
For a matrix G over a subfield E of F, expand each F-valued message in an
E-basis. The support of its image is the union of the supports of its
E-valued components. Consequently d_F(G) = d_E(G).

This also gives a useful implementation option. Every requested field
contains GF(256). One generator with coefficients in GF(256) therefore
has exactly the same symbol distance over all four fields. Multiplication
by its coefficients can act bytewise in an explicitly chosen GF(256)-basis.
An arbitrary consumer's polynomial basis need not have that layout; any
required conversion must be included in timings.

To benefit from the larger alphabet rather than only wider records,
some coefficients must leave that small subfield. One possible split uses
fixed GF(256) circuits for the outer and inner maps, with larger-field
coefficients only in the small persistent state. The current native-field
prototypes instead use low-degree polynomial coefficients throughout.
Those coefficients leave GF256 without requiring large products.

For orientation, the q-ary random-linear-code GV reference gives:

| m | Distance at R = 1/2 | Rate at distance 10% |
|---:|---:|---:|
| 8 | 38.0472% | 0.841446 |
| 32 | 46.8838% | 0.885344 |
| 64 | 48.4386% | 0.892672 |
| 128 | 49.2189% | 0.896336 |

These are asymptotic existence benchmarks, not upper bounds or SPIN
certificates. They use R = 1 - H_q(delta); see the
[q-ary GV derivation](https://www.cs.cmu.edu/~venkatg/teaching/codingtheory/notes/notes2.pdf).
The model recomputes the table with a stable entropy formula.

The first less-expansion screens should use R = 3/4 and 4/5 over GF(256).
At larger m, R = 7/8 is another useful target. These are research targets,
not a change to the public parameter API. Lowering the rate also opens a
small-field high-distance option: at 49% distance the GF256 GV reference
rate is about 0.385382. More output may be worthwhile if it removes enough
expensive full-field arithmetic; that tradeoff needs an end-to-end measurement.

## A cheap baseline and its obstruction

Consider a scalar-state inner. Sample a_i,b_i independently from F minus
zero during setup, then fix them for every message. For input X in F^N,
start Q_0 = 0 and set

    Y_i = X_i + Q_i,
    Q_(i+1) = a_i Q_i + b_i X_i.

The nonzero restriction matters: when X_i = 0, a nonzero state cannot die.
Allowing a_i = 0 would introduce a state-erasure probability of 1/q at
every zero-input step.

There is an equivalent distribution requiring no multiplication on the
loop-carried state. Sample independent nonzero scales delta_i and gamma_i:

    Y_i = X_i + delta_i * sum_(j<i) gamma_j X_j.

The unused delta_0 and gamma_(N-1) may be arbitrary. For i >= 1,
delta_i is the prefix product of a_0 through a_(i-1), and
gamma_j = b_j / delta_(j+1). This change of variables is bijective on the
used nonzero scales, so they can be sampled directly and independently.
The transpose uses a scaled suffix sum:

    (I^T Z)_i = Z_i + gamma_i * sum_(j>i) delta_j Z_j.

For each fixed input, its exact two-category transfer is small. Let z mark
output weight, p = 1/(q-1), and f(z) = p + (1-p)z. Rows and columns index
zero state and uniform nonzero state. Zero and nonzero input steps have
transfer matrices

    T_0(z) = [[1, 0], [0, z]],
    T_1(z) = [[0, z], [p*f(z), (1-p)*f(z)]].

Fresh coefficients make the next-state category independent of the output
category in the nonzero/nonzero case. This is a fixed-input enumerator;
it is not by itself a minimum-distance certificate.

More importantly, a small state allows an adversarial message to silence
long zero-input gaps. This limitation does not vanish as q grows.

**Quiet-gap bound.** Fix a k-dimensional outer constituent whose b outputs
route to positions p_1 < ... < p_b in {0,...,N-1}. Set every other message
block to zero. Suppose the inner is causal and homogeneous linear, with
zero initialization, an s-dimensional F-state for s >= 1,
nonzero scalar direct-feedthrough at every position, and
no output flush. Fix all its coefficients and the route.

Let ell_j = p_(j+1)-p_j-1 for j < b and ell_b = N-p_b-1. Order these b gap
lengths from largest to smallest. With r = floor((k-1)/s),

    d_min <= N - p_1 - sum_(j=1)^r ell_(j).

To prove the bound, require the state to be zero at the start of the r
selected gaps. These are at most rs < k homogeneous equations on the local
message. A nonzero message satisfies them all. Its outputs vanish in the
selected gaps and before p_1. Injectivity ensures that its remaining
output is nonzero. Counting those positions proves the inequality.

For a batched inner, apply this argument to wholly inactive steps and
charge the full output width of every remaining step. Do not apply the
symbolwise formula inside a step that mixes several input symbols.

This bound is a rejection test for actual routes. It is one reason to
retain SPIN's region structure and explore several field-state dimensions,
rather than commit to the cheapest scalar-state inner.

## First batched candidate

Take a genuine GRS[b,k] outer over GF(256), interpreted over F, and route
field symbols with a SPIN-style region permutation. Choose t = 8 or 16
symbols per inner step and begin with s = 2 state symbols. Larger s remains
available if the gap bound or certificate requires it.

Fix A in F^(t by s) and C in F^(s by t), initially over GF(256).
Sample independent M_i uniformly from GL(s,F) during setup. Given
X_i in F^t and Q_i in F^s, use

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The input/output map is unit block triangular. Hence the inner is
invertible even though its state starts at zero and is not flushed.
When Q_i is nonzero, M_i Q_i is exactly uniform among nonzero state
vectors. Zero input cannot extinguish that state.

Choose A as an injective [t,s] GRS map and C as a parity-check matrix for
an [t,t-s] MDS code. Then nonzero states spread across at least t-s+1
output coordinates, and nonzero inputs of weight at most s cannot lie in ker C.
The joint behavior of A, C, and routed inputs still needs analysis.

At s = 2, the state update uses four full-field multiplications per step,
plus additions. The A/C maps and outer have separate costs; four divided
by t is not the total arithmetic cost per output symbol. Test GFNI-oriented
GF(256) circuits and CLMUL-based larger-field state updates only after
fixing their bases and matching their scalar reference maps.

The formal transpose of one step, given adjoints V_i for Y_i and R_(i+1)
for the next state, is

    U_i = V_i + C^T R_(i+1),
    R_i = A^T V_i + M_i^T R_(i+1).

Start R_final = 0 and scan backwards. These are transposes over F, not
binary bit-matrix adjoints under an undeclared representation.

## Proof work to do

For a fixed source support S, let V_S be the shortened message space and
D_S its dimension. A branch imposes output-zero equations with rank r.
Every branch realized by a nonzero message has r <= D_S-1. Thus for tau >= 1,

    1{a bad message from S exists}
      <= tau^(D_S-1) * sum_(retained bad branches) tau^(-r).

For an MDS[b,k] constituent, shortening to j designated coordinates has
dimension max(0,j-(b-k)). This supplies the local support polynomial

    1 + sum_(j=b-k+1)^b binom(b,j) z^j tau^(j-(b-k)).

Retain its exact support information while deriving the region-placement
factor and inner continuation bound. A fixed-message uniform-state law
does not establish a uniform bound after earlier zero equations.

One useful local identity is available for the GL state maps. Let H be a
fixed s-by-D state-form matrix of rank h >= 1 in the current quotient.
Let a be a fixed nonzero row and let M be fresh uniform GL(s,F). Then aM
is uniform over the q^s-1 nonzero rows. Therefore aMH is zero with
probability (q^(s-h)-1)/(q^s-1), and equals each specified nonzero form
in the row space of H with probability q^(s-h)/(q^s-1).
This follows by counting the fibers of the rank-h map v -> vH.
Joint output rows require a joint-rank calculation, not independent use
of these marginal probabilities.

The immediate next experiment is a rank-aware continuation bound for this
batched inner, together with the quiet-gap screen. Only surviving designs
should receive optimized kernels. Evaluate both outer rate and state
dimension before tuning individual instructions.

## Validation and performance contract

The initial model has 14 tests, including exact GF4/GF256 transition
counts, transpose and inverse identities, and the GL(2,GF4) rank law.
The other directories add arithmetic, proof replay, and independent audits.
Run all three suites explicitly; discovery from the root does not recurse
into these non-package directories:

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin -p "test_*.py" -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/lane_mixing -p "test_*.py" -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p "test_*.py" -v
```

Keep precomputed transposed encoding as the first performance target.
Report setup and encoding separately, and compare both symbol count and
physical bytes. A 128-bit record contains sixteen GF256 lanes but only
one GF(2^128) symbol. Packing, basis conversion, routing, and output stores
must remain inside end-to-end encode timings when required by the API.
Run benchmarks serially. The initial model's 14 checks have been supplemented
by field arithmetic, finite-certificate, local-transfer, and native-implementation
tests. Current measurements and exact commands are in [RESULTS.md](RESULTS.md)
and the [implementation guide](../../../spin/experiments/char2_field_spin/README.md).

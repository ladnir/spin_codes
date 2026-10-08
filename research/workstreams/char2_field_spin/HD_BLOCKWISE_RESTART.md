# High-distance SPIN: return to a blockwise inner

Restarted on 2026-10-06 at the user's request. The
[accumulator detour](ACCUMULATOR_DETOUR_CHECKPOINT.md) is frozen and indexed;
the certified high-rate winner and production binary codes remain unchanged.

Current direction: [independent randomized outer outputs](HD_RANDOMIZED_OUTER.md)
restore the binary proof's expected-measure comparison. The checked prototype
adds about0.48ms to the shared FFT encoder, with new directed class bounds
but no complete40% certificate. Keep the measured6.35-ms identity-feedthrough
code as a reference. The [shared-transform proposal](HD_BLOCKWISE_SHARED_MAPS.md)
uses one eight-point input mixer to perform feedthrough, state emission,
and feedback together. Its [native comparison](HD_BLOCKWISE_SHARED_NATIVE.md)
now gives6.54ms for FFT,7.65ms for the59-product cube, and7.00ms for its
whole-vector rewrite. These fixed-outer controls are preserved unchanged.
The next proof work tracks weighted state densities and joint emission/syndrome
moments for the randomized outer. Lower multiplication counts alone have
not produced a speedup.

The latest [quiet-subcode obstruction](HD_BLOCKWISE_QUIET_SUBCODE.md) is
deterministic: identity feedthrough requires enough state per packet, not
merely enough state per regional window. At rate 1/2 and width64, states4
and8 have relative distance at most12.5% and25%, for every setup. A larger
outer cannot repair this obstruction. The t64/s16 baseline survives this
necessary check; smaller-packet controls t32/s8 and t16/s4 also survive it.
None has a complete high-distance certificate.

## Design and acceptance contract

Use a fast local outer, one sampled routing stage, and one batched stateful
inner. There is no separate pre-accumulator. Process 128-bit records with
SIMD over a characteristic-two field. Setup is sampled once and precomputed.
The inner state persists across packets and regions, starting at zero.

At K=2^20, first retain rate 1/2 and seek a complete finite certificate for
at least 40% symbol distance with a 40-bit setup-failure target. The same
construction must encode forward below 5 ms on the RS measurement host.
Then seek approximately 3 ms and 45--50% distance. Point bounds, operation
counts, and transposed timings cannot satisfy those acceptance conditions.

The primary geometry is direct K -> 2K, with every inner output retained.
Extra expansion followed by puncturing remains a separately labeled control,
not a hidden source of extra work. Do not return to a multi-scan accumulator
pipeline without explicitly revisiting this design choice with the user.

## First candidate family

Start with GRS[128,64] over embedded GF256, interpreted over GF65536.
The outer has 16384 independent rows at K20. Collect each coordinate into
one region and independently permute rows within each region. This is one
structured routing stage, not 128 sequential full-codeword permutations.
Eight GF65536 lanes share the scalar code in each 128-bit record.

Compare packet widths t=32 and 64 with s=16 or 32 state symbols. These
are search parameters, not selected certified operating points. In a
direct IMT-style family, fix maps A:F^s->F^t and C:F^t->F^s and compute

    Y = X + A Q;
    Q_next = M_i Q + C X.

Here X is the current routed packet and Q the incoming state. The maps
are fixed at setup. M_i is a sampled invertible state mixer, independent
across packets and from the route. Specify its distribution before using
any state-refresh probability. This packet map has identity input/output
feedthrough, so the full causal encoder is invertible even without flushing
the terminal state.

Two map options merit comparison: fast rectangular GRS-derived maps with
known rank/distance properties, and cheaper sparse maps with explicitly
computed support properties. Neither inherits the binary IMT certificate
merely by replacing bits with field elements. The exact map choices and
the relation between A and C must be recorded before certification.

The [field-valued IMT study](HD_BLOCKWISE_INNER.md) examines a rank-one
invertible update M_i. It is intended to amortize randomness over a packet,
not to assign independent field coefficients to several full-record scans.
The existing diagonal/full-Cauchy family is a useful comparison, but its
extra input-mixing map must be included in the arithmetic cost and proof.

The first checked mixer samples u uniformly from F^s excluding zero,
then samples v uniformly subject to v^T*u != 1, and sets M=I+u*v^T.
It costs one state dot product and one scaled-vector update: 2s field
products per packet. At (t,s)=(64,16), that is 0.5 field products per
output record, before the fixed A/C maps and routing. The parked pipeline
used five such products per output record. This tenfold reduction concerns
only variable-coefficient arithmetic, not total encoder time.

The new local law is exact but not an independent fresh state. For a
nonzero state, one round is a positive mixture of its original point,
its nonzero scalar line, and the uniform nonzero state distribution.
The total point/line mass is 1/|F|; after r consecutive independent
refresh factors before one feedback addition it is |F|^-r. This is not
a claim about r packets with intervening inputs. The associated
cancellation terms must remain in the proof.
One round is the initial cost control; a second round is a proof-strength
comparison, not a silently added requirement.

## Why begin with a 128-symbol outer?

The old direct GRS[64,32] construction with only four state symbols has
an actual 40%-distance obstruction. It is not a general impossibility
result for blockwise SPIN. For K20 and width16, its equal-window expected
shortened source dimension is about 2835, versus only 256 exit constraints.

For GRS[128,64], that expected dimension falls to about 331. A state of
dimension16 supplies many more exit constraints. This removes that specific
lower-bound rejection, but does not establish good distance. Fixed maps,
state refresh, and the full set of supports still require analysis.

The [geometry screen](HD_BLOCKWISE_GEOMETRY.md) tests packet-aligned,
balanced regional windows and accounts for the terminal boundary. Width
must not appear beneficial merely because a coarse window was rounded
down. GRS[256,128] is the next outer-size control, especially for 45%.
Increasing the outer is a measured tradeoff, not a free proof improvement.

The first prefix-window screen gave necessary thresholds of46 state symbols
for the64-symbol outer and three for the128-symbol outer at40% distance.
These historical numbers do not account for adjacent windows joining across
regional boundaries. The [partial-row analysis](HD_BLOCKWISE_PARTIAL.md)
includes those joins. More importantly, the quiet-subcode bound applies
to every route and supersedes these window screens for the rejected
t64/s4 and t64/s8 geometries.

## Proof and implementation sequence

1. Run finite dimension obstructions before expensive proof searches.
   A zero obstruction lower bound is only permission to continue.
2. Specify the full packet map and mixer distribution. Check field
   arithmetic, linearity, causal inverse, and exact small-field state laws.
3. Build local output/cancellation bounds retaining actual state dependence.
   Separate existence or rank tails from uncapped codeword counts where
   rare large kernels would dominate the latter.
4. Implement surviving maps with the existing outer and GF16 primitives.
   Fuse routed loads into packet arithmetic; avoid a separate shuffle pass.
   Compare against a scalar map before timing.
5. Benchmark full one-core forward calls serially, then close all finite
   source classes for the exact measured distribution.

The [cost review](HD_BLOCKWISE_COST.md) identifies reusable kernels. Fixed
embedded GF256 maps can remain bytewise GFNI even when random coefficients
use GF65536. Larger packets need bounded register tiles; a generic loop
or large live packet batch must not replace a tuned schedule by convenience.

The [fixed-map study](HD_BLOCKWISE_MAPS.md) selects the t64/s16 maps
and proves a local weighted-state bound. The [global analysis](HD_BLOCKWISE_GLOBAL.md)
caps one complete-row routing-trap family. The
[support composition](HD_BLOCKWISE_COMPOSITION.md) closes every one-active-row
message at40% distance, but the many-row classes remain open. A single
common fractional exponent is insufficient, even with exact local bounds.
No full certificate is claimed.
The subsequent [structural redesign](HD_BLOCKWISE_SHARED_MAPS.md) replaces
identity feedthrough by eight-point MDS mixing, shared with the state maps.
Its [cube implementation model](HD_BLOCKWISE_SHARED_FACTOR.md) needs59 packed
GFNI arithmetic calls per64 outputs, without a native timing claim. The
[joint-rank interface](HD_BLOCKWISE_SHARED_SUPPORT.md) and
[sparse quiet-class bound](HD_BLOCKWISE_SHARED_CUTS.md) now supply exact proof
progress, while [the tiled budget](HD_BLOCKWISE_SHARED_BUDGET.md) rules out
overly small or infrequently updated states. Shift64 is the current fixed-map
choice; the shift8 checkpoint and every earlier receipt remain preserved.
The [native forward checkpoints](HD_BLOCKWISE_NATIVE.md) retain the first
t64/s16 candidate at about9.0ms. A new GRS[256,128], t16/s4 control reaches
about6.35ms with non-temporal outer stores and a timed sequential region
copy before fused routed loads. Its one-active-row bound exceeds223 bits;
the many-row certificate remains open. The old8.3-ms number belongs to the
archived accumulator design, not either direct candidate.

## Restart validation

The restart adds ten exact geometry tests and eight algebraic inner tests.
All 176 high-distance helper tests pass. An additional independent
exhaustive comparison over the 16-element field and a two-coordinate
state confirms the rank-one image law. This small test field is not the
65536-element coefficient field proposed for the native candidate.

The initial restart made no native changes. The subsequent checkpoint
adds fixed-map and global-bound tests, bringing the high-distance helper
suite to190 passing tests. Native scalar, alias, route-width, guard, and
ASAN/UBSAN checks pass. No production edit, commit, or push was made.
The next checkpoint adds a deterministic quiet-subcode theorem, paired-window
checks, and capped support composition. The
[interval-span bound](HD_BLOCKWISE_INTERVAL_SPANS.md) then controls one
dimension-forced complete-row trap family for t16/s4 below 2^-528.39.
All219 high-distance helper tests
pass. The new t16/s4 native control passes Release and ASAN/UBSAN scalar,
alias, offset, guard, and route-width checks. Next hold that measured
construction fixed and address the many-row proof. No production promotion,
commit, or push was made.
